--[[
backwards-Euler integrator
based on GMRES, easily swappable for any other OpenCL krylov solver of your choice
(found in solver.cl.*)
--]]

local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local math = require 'ext.math'
local tooltip = require 'tooltip'
local Integrator = require 'int.int'
local CLBuffer = require 'cl.obj.buffer'
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

local BackwardEuler = class(Integrator)

BackwardEuler.name = 'backward Euler'

function BackwardEuler:init(solver, args)
	self.solver = solver

	-- gui vars:
	self.last_err = 0
	self.last_iter = 0

-- formerly createBuffers

	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.cons_t)
	local realSize = ffi.sizeof(solver.app.env.real)
	for _,name in ipairs{
		'krylov_b',
		'krylov_x',
		'krylov_dUdt',
	} do
		self[name..'Obj'] = CLBuffer{
			env = solver.app.env,
			name = name,
			type = 'real',
			size = bufferSize / realSize,
		}
	end

-- formerly refreshGridSize
	
	-- the previous call in Solver:iterate is to solver:calcDT
	-- which calls solver.equation:calcInterfaceEigenBasis, which fills the solver.eigenvalues table

	-- function that returns deriv when provided a state vector
	-- dUdtBuf and UBuf are cl.obj.buffer
	local function calc_dU_dt(dUdtBuf, UBuf)
		-- setting solver.UBuf won't make a difference because the kernels are already bound to the original solver.UBuf
		-- copying over solver.UBuf will overwrite the contents of 'x'
		-- so I need to (a1) make a new buffer for 'x'
		-- (a2) here, copy UBuf into solver.UBuf before running the kernel
		-- or (b) have Solver:calcDeriv accept a param for the 'UBuf' 
		-- ... but that'll be hard since UBuf is exchanged with getULRBuf with usePLM
--print'\nUBuf:' solver:printBuf(UBuf) print(debug.traceback(),'\n\n')
		solver.UBufObj:copyFrom(UBuf)
		dUdtBuf:fill()

		-- TODO should 'dt' be passed along as well?
		-- typically only dU/dt is calculated which, for first order, doesn't depend on dt
		-- however for some 2nd order flux limiters, 'dt' is needed 
		-- ... is that the 'dt' of the overall step, 
		-- or of the individual step within the integrator overall step?
		self.integrateCallback(dUdtBuf)
		--solver:calcDeriv(dUdtBuf, self.linearSolverDT)

--print'\ndUdtBuf:' solver:printBuf(dUdtBuf) print(debug.traceback(),'\n\n')
--solver:checkFinite(dUdtBuf)
	end

	local mulWithoutBorder = solver.domain:kernel{
		name = 'RoeImplicitLinearized_mulWithoutBorder',
		header = solver.codePrefix,
		argsOut = {
			{name='y', type=solver.eqn.cons_t, obj=true},
		},
		argsIn = {
			{name='solver', type=solver.solver_t, obj=true},
			{name='a', type=solver.eqn.cons_t, obj=true},
			{name='b', type=solver.eqn.cons_t, obj=true},
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
	
	local numreals = solver.numCells * solver.eqn.numStates
	local volumeWithoutBorder = tonumber(solver.sizeWithoutBorder:volume())
	local numRealsWithoutBorder = volumeWithoutBorder * solver.eqn.numStates
	
	local sum = solver.app.env:reduce{
		size = numreals,
		op = function(x,y) return x..' + '..y end,
	}
	local dotWithoutBorder = function(a,b)
		mulWithoutBorder(sum.buffer, solver.solverBuf, a, b)
		return sum()
	end

	local restart = args and args.restart or 10

	local linearSolverArgs = {
		env = solver.app.env,
		x = self.krylov_xObj,
		size = numreals,
		epsilon = args and args.epsilon or 1e-10,
		--maxiter = 1000,
		restart = restart,
		maxiter = restart * numreals,
		-- logging:
		errorCallback = function(err, iter, x, rLenSq)
			self.last_err = err
			self.last_iter = iter
			
			if not math.isfinite(err) then
				print("got non-finite err: "..err)	-- error?
				return true	-- fail
			end
		end,
		dot = function(a,b)
			return dotWithoutBorder(a,b) / numRealsWithoutBorder
		end,
	}

	-- [=[ backward Euler
	linearSolverArgs.b = self.krylov_bObj
	linearSolverArgs.A = function(UNext, U)
		local dUdt = self.krylov_dUdtObj
		calc_dU_dt(dUdt, self.krylov_bObj)
		
		--UNext = U - dt * calc_dU_dt(lastU)
--print'\nU:' solver:printBuf(U) print(debug.traceback(),'\n\n')
--print('self.linearSolverDT', self.linearSolverDT)		
		
		self.linearSolver.args.mulAdd(UNext, U, dUdt.obj, -self.linearSolverDT)
--print'\nUNext:' solver:printBuf(UNext) print(debug.traceback(),'\n\n')
	
		--[[ do I need to apply the boundary?
		-- doesn't seem to make a difference
		-- TODO don't even include the boundary cells in the GMRES
		-- in fact, should I even be storing ghost cells in memory, or just using conditions to provide their values?
		for _,obj in ipairs(solver.boundaryKernelObjs) do
			obj.obj:setArg(1, UNext)
		end
		solver:boundary()
		for _,obj in ipairs(solver.boundaryKernelObjs) do
			obj.obj:setArg(1, solver.UBuf)
		end
		--]]
	end
	--]=]
	--[=[ crank-nicolson - converges faster
	linearSolverArgs.b = (function(UBuf)
		UBuf = UBuf + .5 * dt * calc_dU_dt(UBuf)
		solver.boundaryMethod(UBuf)
		return UBuf
	end)(UBuf)
	linearSolverArgs.A = function(UBuf)
		UBuf = UBuf - .5 * dt * calc_dU_dt(UBuf)
		solver.boundaryMethod(UBuf)
		return UBuf
	end
	--]=]

	-- set up gmres solver here
	self.linearSolver = ThisGMRES(linearSolverArgs)
end

-- step contains integrating flux and source terms
-- but not post iterate
function BackwardEuler:integrate(dt, callback)
	local solver = self.solver
	-- this is a call to 'calcDeriv'
	self.integrateCallback = callback
	
	self.linearSolverDT = dt
	-- UBuf needs to be overwritten to pass on to the calcFluxDeriv
	-- (TODO make calcFluxDeriv accept a parameter)
	self.krylov_bObj:copyFrom(solver.UBufObj)
--print'\nself.krylov_bObj:' self:printBuf(self.krylov_bObj) print(debug.traceback(),'\n\n')
	self.krylov_xObj:copyFrom(solver.UBufObj)
	self.linearSolver()
	solver.UBufObj:copyFrom(self.krylov_xObj)
end

function BackwardEuler:updateGUI()
	tooltip.numberTable('Krylov epsilon', self.linearSolver.args, 'epsilon')
	tooltip.intTable('GMRES restart', self.linearSolver.args, 'restart')
	tooltip.intTable('Krylov maxiter', self.linearSolver.args, 'maxiter')	-- typically restart * number of reals = restart * numCells * number of states
	-- read-only:
	ig.igText('err = '..self.last_err)	-- this is |r|
	ig.igText('iter = '..self.last_iter)
end

return BackwardEuler
