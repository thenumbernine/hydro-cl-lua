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
local template = require 'template'
local CLBuffer = require 'cl.obj.buffer'


--local CLKrylov = require 'solver.cl.conjgrad'
--local CLKrylov = require 'solver.cl.conjres'
local CLKrylov = require 'solver.cl.gmres'

local ThisKrylov = class(CLKrylov)

function ThisKrylov:newBuffer(name)
	if not self.cache then self.cache = {} end
	local cached = self.cache[name]
	if cached then return cached end
	cached = ThisKrylov.super.newBuffer(self, name)
	cached:fill()
	self.cache[name] = cached
	return cached
end

local BackwardEuler = class(Integrator)

BackwardEuler.name = 'backward Euler'

function BackwardEuler:init(solver, args)
	self.solver = solver
	self.verbose = cmdline.intVerbose or (args and args.verbose) or nil

	-- gui vars:
	self.lastResidual = 0
	self.lastIter = 0

-- formerly createBuffers

	local bufferSize = solver.volumeWithoutBorder * solver.eqn.numIntStates
	for _,name in ipairs{
		'krylov_b',
		'krylov_x',
		'krylov_dUdt',
	} do
		self[name..'Obj'] = CLBuffer{
			env = solver.app.env,
			name = name,
			type = solver.app.real,
			count = bufferSize,
		}
	end
	-- full buffer, with ghost cells and non-integratable state variables
	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.app.real,
		count = solver.numCells * solver.eqn.numStates,
	}

	local copyBufferWithOrWithoutGhostProgram = solver.Program{
		env = solver.app.env,
		code = template(
solver.codePrefix..[[
<? local range = require 'ext.range' ?>
kernel void copyBufferWithoutGhostToBufferWithGhost(
	constant <?=solver.solver_t?>* solver,
	global real* dst,
	global const real* src
) {
	SETBOUNDS_NOGHOST();
	for (int j = 0; j < numIntStates; ++j) {
		dst[j + numStates * index] = src[j 
			+ numIntStates * (
				(i.x - numGhost) 
<? if solver.dim > 1 then ?>				
				+ (solver->gridSize.x - 2 * numGhost) * (
					(i.y - numGhost) 
<? if solver.dim > 2 then ?>					
					+ (solver->gridSize.x - 2 * numGhost) * (i.z - numGhost)
<? end ?>				
				)
<? end ?>			
			)];
	}
}

kernel void copyBufferWithGhostToBufferWithoutGhost(
	constant <?=solver.solver_t?>* solver,
	global real* dst,
	global const real* src
) {
	SETBOUNDS_NOGHOST();
	for (int j = 0; j < numIntStates; ++j) {
		dst[j 
			+ numIntStates * (
				(i.x - numGhost) 
<? if solver.dim > 1 then ?>				
				+ (solver->gridSize.x - 2 * numGhost) * (
					(i.y - numGhost) 
<? if solver.dim > 2 then ?>					
					+ (solver->gridSize.x - 2 * numGhost) * (i.z - numGhost)
<? end ?>				
				)
<? end ?>			
			)] = src[j + numStates * index];
	}
}
]],		{
			solver = solver,
		})
	}
	copyBufferWithOrWithoutGhostProgram:compile()
	self.copyWithoutToWithGhostKernel = copyBufferWithOrWithoutGhostProgram:kernel{
		name='copyBufferWithoutGhostToBufferWithGhost',
		domain=solver.domainWithoutBorder,
	}
	self.copyWithToWithoutGhostKernel = copyBufferWithOrWithoutGhostProgram:kernel{
		name='copyBufferWithGhostToBufferWithoutGhost',
		domain=solver.domainWithoutBorder,
	}

-- formerly refreshGridSize
	
	-- the previous call in Solver:iterate is to solver:calcDT
	-- which calls solver.equation:calcInterfaceEigenBasis, which fills the solver.eigenvalues table

	-- function that returns deriv when provided a state vector
	-- dUdtBuf and UBuf are cl.obj.buffer
	local function calc_dU_dt(dUdtBuf, UBuf)
--print'\nUBuf:' solver:printBuf(UBuf, nil, solver.eqn.numIntStates)

		self.copyWithoutToWithGhostKernel(solver.solverBuf, solver.UBufObj, UBuf)
		
		solver:boundary()	
if solver.checkNaNs then assert(solver:checkFinite(derivBufObj)) end
		if solver.eqn.useConstrainU then
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
			solver.constrainUKernelObj(solver.solverBuf, solver.UBuf)
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
		end
		
		self.derivBufObj:fill()	-- fill with zero
		self.integrateCallback(self.derivBufObj)
		self.copyWithToWithoutGhostKernel(solver.solverBuf, dUdtBuf, self.derivBufObj)

--print'\ndUdtBuf:' solver:printBuf(dUdtBuf, nil, solver.eqn.numIntStates)
		if solver.checkNaNs then assert(solver:checkFinite(dUdtBuf)) end
	end

	local volumeWithoutBorder = solver.volumeWithoutBorder
	local numreals = volumeWithoutBorder * solver.eqn.numIntStates
	
	local restart = cmdline.intBERestart or (args and args.restart) or 20

	local linearSolverArgs = {
		env = solver.app.env,
		x = self.krylov_xObj,
		count = numreals,
		epsilon = cmdline.intBEEpsilon or (args and args.epsilon) or 1e-10,
		--maxiter = 1000,
		restart = restart,
		maxiter = cmdline.intBEMaxIter or (args and args.maxiter) or restart * numreals,
		-- logging:
		errorCallback = function(residual, iter, x, rLenSq)
			self.lastResidual = residual
			self.lastIter = iter
			if self.verbose then
				print('t', solver.t, 'iter', iter, 'residual', residual)
			end
--if iter < numreals then return false end
			if not math.isfinite(residual) then
				print("got non-finite residual: "..residual)	-- error?
				return true	-- fail
			end
--			if residual < self.linearSolver.args.epsilon then return true end
		end,
	}

	-- [=[ backward Euler
	linearSolverArgs.b = self.krylov_bObj
	linearSolverArgs.A = function(UNext, U)
		local dUdt = self.krylov_dUdtObj
		calc_dU_dt(dUdt, self.krylov_bObj)
		
		--UNext = U - dt * calc_dU_dt(lastU)
--print'\nU:' solver:printBuf(U, nil, solver.numIntStates) 
--print('self.linearSolverDT', self.linearSolverDT)		
		
		self.linearSolver.args.mulAdd(UNext, U, dUdt.obj, -self.linearSolverDT)
--print'\nUNext:' solver:printBuf(UNext, nil, solver.numIntStates)
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
	self.linearSolver = ThisKrylov(linearSolverArgs)

	local oldDot = self.linearSolver.args.dot
	self.linearSolver.args.dot = function(a,b)
		return oldDot(a,b) / math.sqrt(numreals)
	end
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
	self.copyWithToWithoutGhostKernel(solver.solverBuf, self.krylov_bObj, solver.UBufObj)
--print'\nself.krylov_bObj:' solver:printBuf(self.krylov_bObj, nil, solver.eqn.numIntStates)
	self.copyWithToWithoutGhostKernel(solver.solverBuf, self.krylov_xObj, solver.UBufObj)
	self.linearSolver()
	self.copyWithoutToWithGhostKernel(solver.solverBuf, solver.UBufObj, self.krylov_xObj)

	solver:boundary()	
if solver.checkNaNs then assert(solver:checkFinite(derivBufObj)) end
	if solver.eqn.useConstrainU then
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
		solver.constrainUKernelObj(solver.solverBuf, solver.UBuf)
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
	end

end

function BackwardEuler:updateGUI()
	tooltip.numberTable('Krylov epsilon', self.linearSolver.args, 'epsilon')
	tooltip.intTable('GMRES restart', self.linearSolver.args, 'restart')
	tooltip.intTable('Krylov maxiter', self.linearSolver.args, 'maxiter')	-- typically restart * number of reals = restart * numCells * number of states
	-- read-only:
	ig.igText('residual = '..self.lastResidual)	-- this is |r|
	ig.igText('iter = '..self.lastIter)
end

return BackwardEuler
