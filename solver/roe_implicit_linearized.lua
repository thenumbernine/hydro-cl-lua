local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'solver.roe'
local ffi = require 'ffi'
local math = require 'ext.math'
local CLGMRES = require 'solver.cl.gmres'

--[[
local tolua = require 'ext.tolua'
local gputostr_n
local gputostr_cmds
function gputostr(x)
	local ptr = ffi.new('real['..gputostr_n..']')
	gputostr_cmds:enqueueReadBuffer{buffer=x, block=true, size=gputostr_n*ffi.sizeof'real', ptr=ptr}
	local s = ''
	for i=0,gputostr_n-1 do
		s = s .. '\t' ..ptr[i]
		if i % 6 == 5 then s = s .. '\n' end
	end
	return s 
end
--]]

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
	-- dUdtBuf is a cl.buffer
	-- UBuf is a cl.obj.buffer
	local function calc_dU_dt(dUdtBuf, UBuf)
		-- setting self.UBuf won't make a difference because the kernels are already bound to the original self.UBuf
		-- copying over self.UBuf will overwrite the contents of 'x'
		-- so I need to (a1) make a new buffer for 'x'
		-- (a2) here, copy UBuf into self.UBuf before running the kernel
		-- or (b) have Solver:calcDeriv accept a param for the 'UBuf' 
		-- ... but that'll be hard since UBuf is exchanged with getULRBuf with usePLM
		self.UBufObj:copyFrom(UBuf)
		self.app.cmds:enqueueFillBuffer{buffer=dUdtBuf, size=self.volume * ffi.sizeof(self.eqn.cons_t)}
		self:calcDeriv(dUdtBuf, self.linearSolverDT)
	end

	local size = require 'ext.range'(self.dim):map(function(i)
		return self.gridSize:ptr()[i-1]
	end)
	local env = self.app.env
	
	local domain = env:domain{size = size}

	local mulWithoutBorder = domain:kernel{
		header = self.codePrefix,
		size = size,
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
	local sum = env:reduce{
		size = numreals,
		op = function(x,y) return x..' + '..y end,
	}
	local dotWithoutBorder = function(a,b)
		mulWithoutBorder(sum.buffer, a, b)
		return sum()
	end

	local function wrapBufferObj(field)
		return setmetatable({
			env = env,
			name = field,
			type = 'real',
			size = numreals,
			obj = self[field],
		}, require 'cl.obj.buffer')
	end

	self.UBufObj = wrapBufferObj'UBuf'
	self.RoeImplicitLinear_bObj = wrapBufferObj'RoeImplicitLinear_b'
	self.RoeImplicitLinear_xBufObj = wrapBufferObj'RoeImplicitLinear_x'

	local linearSolverArgs = {
		env = env,
		x = self.RoeImplicitLinear_xBufObj,
		size = numreals,
		epsilon = 1e-10,
		--maxiter = 1000,
		maxiter = 10 * numreals,
		restart = 10,
		-- logging:
		errorCallback = function(err, iter, x, rLenSq, bLenSq)
			--print('t', self.t, 'iter', iter, 'err', err, 'rLenSq', rLenSq, 'bLenSq', bLenSq)
			if not math.isfinite(err) then
				error("got non-finite err: "..err)
			end
		end,
		dot = dotWithoutBorder,
	}

	-- [=[ backward Euler
	linearSolverArgs.b = self.RoeImplicitLinear_bObj
	-- usage of A:
	-- A(r, x)
	-- A(w, v[i])
	linearSolverArgs.A = function(UNext, U)
		-- if this is a linearized implicit solver
		-- then the matrix should be computed before invoking the iterative solver
		-- which means the matrix coeffiicents shouldn't be changing per-iteration
		-- which means calc_dU_dt() should be based on the initial state and not the iterative state
		-- use the integrator's deriv buffer and don't use this?
		local dUdt = self.integrator.derivBuf
		calc_dU_dt(dUdt, self.RoeImplicitLinear_bObj)
		
		--UNext = U - dt * calc_dU_dt(lastU)
		self.linearSolver.args.mulAdd(UNext, U, dUdt, -self.linearSolverDT)
	
		--[[ do I need to apply the boundary?
		self.UBufObj:copyFrom(UNextBuf)
		self:boundary()
		UNextBuf:copyFrom(self.UBufObj)
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
end

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
	self.linearSolverDT = dt
	-- UBuf needs to be overwritten to pass on to the calcFluxDeriv
	-- (TODO make calcFluxDeriv accept a parameter)
	self.RoeImplicitLinear_bObj:copyFrom(self.UBufObj)
	self.RoeImplicitLinear_xBufObj:copyFrom(self.UBufObj)
	self.linearSolver()
	self.UBufObj:copyFrom(self.RoeImplicitLinear_xBufObj)
end

return RoeImplicitLinearized
