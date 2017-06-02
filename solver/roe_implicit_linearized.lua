local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'solver.roe'
local ffi = require 'ffi'
local math = require 'ext.math'
local CLGMRES = require 'solver.cl.gmres'

local ThisGMRES = class(CLGMRES)

function ThisGMRES:newBuffer(name)
	if not self.cache then self.cache = {} end
	local cached = self.cache[name]
	if cached then return cached end
	cached = ThisGMRES.super.newBuffer(self, name)
	cached:fill()
	self.cache[name] = cached
	return cached
end

-- technically this backwards Euler integrator works on any solver with calcDeriv implemented
-- which means maybe I could implement it as an integrator instead of a solver ...
local RoeImplicitLinearized = class(Roe)
RoeImplicitLinearized.name = 'RoeImplicitLinearized'

function RoeImplicitLinearized:refreshGridSize(...)
	RoeImplicitLinearized.super.refreshGridSize(self, ...)

	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	-- function that returns deriv when provided a state vector
	-- dUdtBuf and UBuf are cl.obj.buffer
	local function calc_dU_dt(dUdtBuf, UBuf)
		-- setting self.UBuf won't make a difference because the kernels are already bound to the original self.UBuf
		-- copying over self.UBuf will overwrite the contents of 'x'
		-- so I need to (a1) make a new buffer for 'x'
		-- (a2) here, copy UBuf into self.UBuf before running the kernel
		-- or (b) have Solver:calcDeriv accept a param for the 'UBuf' 
		-- ... but that'll be hard since UBuf is exchanged with getULRBuf with usePLM
--print'\nUBuf:' self:printBuf(UBuf) print(debug.traceback(),'\n\n')
		self.UBufObj:copyFrom(UBuf)
		dUdtBuf:fill()
		self:calcDeriv(dUdtBuf.obj, self.linearSolverDT)
--print'\ndUdtBuf:' self:printBuf(dUdtBuf) print(debug.traceback(),'\n\n')
--self:checkFinite(dUdtBuf)
	end

	local mulWithoutBorder = self.domain:kernel{
		name = 'RoeImplicitLinearized_mulWithoutBorder',
		header = self.codePrefix,
		argsOut = {
			{name='y', type=self.eqn.cons_t, obj=true},
		},
		argsIn = {
			{name='a', type=self.eqn.cons_t, obj=true},
			{name='b', type=self.eqn.cons_t, obj=true},
		},
		body = [[	
	if (OOB(numGhost, numGhost)) {
		for (int j = 0; j < numStates; ++j) {
			y[index].ptr[j] = 0;
		}
		return;
	}
	for (int j = 0; j < numStates; ++j) {
		y[index].ptr[j] = a[index].ptr[j] * b[index].ptr[j];
	}
]],
	}
	
	local numreals = self.volume * self.eqn.numStates
	local volumeWithoutBorder = tonumber(self.sizeWithoutBorder:volume())
	local numRealsWithoutBorder = volumeWithoutBorder * self.eqn.numStates
	
	local sum = self.app.env:reduce{
		size = numreals,
		op = function(x,y) return x..' + '..y end,
	}
	local dotWithoutBorder = function(a,b)
		mulWithoutBorder(sum.buffer, a, b)
		return sum()
	end

	local linearSolverArgs = {
		env = self.app.env,
		x = self.RoeImplicitLinear_xObj,
		size = numreals,
		epsilon = 1e-10,
		--maxiter = 1000,
		maxiter = 10 * numreals,
		restart = 10,
		-- logging:
		errorCallback = function(err, iter, x, rLenSq)
--print('gmres t', self.t, 'iter', iter, 'err', err, 'rLenSq', rLenSq)
			if not math.isfinite(err) then
				error("got non-finite err: "..err)
			end
		end,
		--[[ only scale down norm, keep dot the same
		dot = dotWithoutBorder,
		norm = function(v)
			local v2 = dotWithoutBorder(v,v)
			return math.sqrt(v2 
				-- need to divide by volume or else 2D will enter non-physical states
				/ numRealsWithoutBorder
				-- but if you divide 1D by volume then it runs really really slow
				* tonumber(self.sizeWithoutBorder.x) * self.eqn.numStates
			)
		end,
		--]]
		-- [[ scale down dot and let norm be sqrt(dot(v,v))
		dot = function(a,b)
			return dotWithoutBorder(a,b) / numRealsWithoutBorder
		end,
		--]]
	}

	-- [=[ backward Euler
	linearSolverArgs.b = self.RoeImplicitLinear_bObj
	linearSolverArgs.A = function(UNext, U)
		local dUdt = self.RoeImplicitLinear_dUdtObj
		calc_dU_dt(dUdt, self.RoeImplicitLinear_bObj)
		
		--UNext = U - dt * calc_dU_dt(lastU)
--print'\nU:' self:printBuf(U) print(debug.traceback(),'\n\n')
--print('self.linearSolverDT', self.linearSolverDT)		
		
--[[
2D GMRES fails.
1D works fine,
but 2D encounters differences with 1D in that,
if we use a L2 norm
and if our 2D data is an extrusion of our 1D data
then |r_2d| = sqrt(n) * |r_1d|	for extruded length n.
From there, v[1] = r/|r| is scaled down an extra factor of sqrt(n).
From there, A(w, v[1]), which builds the orthonormal basis,
accepts a *much smaller* v[1] = r/|r|, yet computes the *same* derivative as the 1D case,
which tends to step 'x' into non-physical values (negative density and energy).
In fact,
dot(w, v[1]) for 1D is positive
and dot(w, v[1]) for 2D is negative.

Maybe my norm function should be norm(r) = sqrt(sum_i r_i / n)?
Then an extruded norm is norm(r_2d) = sqrt(n sum_i r_i / n^2) = norm(r_1d).
Then the norms of extruded data would match that of the original data.
But does GMRES require a particular kind of norm?
The original paper doesn't specify a norm, but does use it hand-in-hand with an inner product.
I'm betting if the norm is changed then the inner product will probably also have to be changed.
For norm(r) = sqrt(dot(r,r)/n), we can define the inner product as dot(r,r)/n.
--]]
		self.linearSolver.args.mulAdd(UNext, U, dUdt.obj, -self.linearSolverDT)
--print'\nUNext:' self:printBuf(UNext) print(debug.traceback(),'\n\n')
	
		--[[ do I need to apply the boundary?
		-- doesn't seem to make a difference
		-- TODO don't even include the boundary cells in the GMRES
		-- in fact, should I even be storing ghost cells in memory, or just using conditions to provide their values?
		self.boundaryKernel:setArg(0, UNext)
		self:boundary()
		self.boundaryKernel:setArg(0, self.UBuf)
		--]]
	end
	--]=]
	--[=[ crank-nicolson - converges faster
	linearSolverArgs.b = (function(UBuf)
		UBuf = UBuf + .5 * dt * calc_dU_dt(UBuf)
		self.boundaryMethod(UBuf)
		return UBuf
	end)(UBuf)
	linearSolverArgs.A = function(UBuf)
		UBuf = UBuf - .5 * dt * calc_dU_dt(UBuf)
		self.boundaryMethod(UBuf)
		return UBuf
	end
	--]=]

	-- set up gmres solver here
	self.linearSolver = ThisGMRES(linearSolverArgs)
end

function RoeImplicitLinearized:createBuffers()
	RoeImplicitLinearized.super.createBuffers(self)
	self:clalloc('RoeImplicitLinear_b', self.volume * ffi.sizeof(self.eqn.cons_t))
	self:clalloc('RoeImplicitLinear_x', self.volume * ffi.sizeof(self.eqn.cons_t))
	self:clalloc('RoeImplicitLinear_dUdt', self.volume * ffi.sizeof(self.eqn.cons_t))
end

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
	self.linearSolverDT = dt
	-- UBuf needs to be overwritten to pass on to the calcFluxDeriv
	-- (TODO make calcFluxDeriv accept a parameter)
	self.RoeImplicitLinear_bObj:copyFrom(self.UBufObj)
--print'\nself.RoeImplicitLinear_bObj:' self:printBuf(self.RoeImplicitLinear_bObj) print(debug.traceback(),'\n\n')
	self.RoeImplicitLinear_xObj:copyFrom(self.UBufObj)
	self.linearSolver()
	self.UBufObj:copyFrom(self.RoeImplicitLinear_xObj)
end

return RoeImplicitLinearized
