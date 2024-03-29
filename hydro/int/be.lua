--[[
backwards-Euler integrator
based on GMRES, easily swappable for any other OpenCL krylov solver of your choice
(found in solver.cl.*)
--]]

local ffi = require 'ffi'
local ig = require 'imgui'
local math = require 'ext.math'
local template = require 'template'
local CLBuffer = require 'cl.obj.buffer'
local Integrator = require 'hydro.int.int'

--local CLKrylov = require 'solver.cl.conjgrad'
local CLKrylov = require 'solver.cl.conjres'
--local CLKrylov = require 'solver.cl.gmres'

local ThisKrylov = CLKrylov:subclass()

function ThisKrylov:newBuffer(name)
	if not self.cache then self.cache = {} end
	local cached = self.cache[name]
	if cached then return cached end
	cached = ThisKrylov.super.newBuffer(self, name)
	cached:fill()
	self.cache[name] = cached
	return cached
end

local BackwardEuler = Integrator:subclass()

BackwardEuler.name = 'backward Euler'

function BackwardEuler:init(solver, args)
	self.solver = solver
	self.verbose = cmdline.intVerbose or (args and args.verbose) or nil

	-- gui vars:
	self.lastResidual = 0
	self.lastIter = 0

-- formerly createBuffers

	local bufferSizeWithoutBorder = solver.volumeWithoutBorder * solver.eqn.numIntStates
	for _,name in ipairs{
		'krylov_b',
		'krylov_x',
		'krylov_dUdt',
	} do
		self[name..'Obj'] = CLBuffer{
			env = solver.app.env,
			name = name,
			type = solver.app.real,
			count = bufferSizeWithoutBorder,
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
		name = 'int-be',
		code = solver.eqn:template(
			solver.modules:getCodeAndHeader(solver.sharedModulesEnabled:keys():append{
				solver.symbols.SETBOUNDS_NOGHOST,
				solver.symbols.solver_macros,
				assert(solver.solver_t),
			}:unpack())
			..[[
<? local range = require 'ext.range' ?>
kernel void copyBufferWithoutGhostToBufferWithGhost(
	constant <?=solver_t?> const * const solver,
	global real * const dst,
	global real const * const src
) {
	<?=SETBOUNDS_NOGHOST?>();
	for (int j = 0; j < numIntStates; ++j) {
		dst[j + numStates * index] = src[j
			+ numIntStates * (
				(i.x - solver->numGhost)
<? if solver.dim > 1 then ?>
				+ (solver->gridSize.x - 2 * solver->numGhost) * (
					(i.y - solver->numGhost) 
<? if solver.dim > 2 then ?>
					+ (solver->gridSize.y - 2 * solver->numGhost) * (i.z - solver->numGhost)
<? end ?>
				)
<? end ?>
			)];
	}
}

kernel void copyBufferWithGhostToBufferWithoutGhost(
	constant <?=solver_t?> const * const solver,
	global real * const dst,
	global real const * const src
) {
	<?=SETBOUNDS_NOGHOST?>();
	for (int j = 0; j < numIntStates; ++j) {
		dst[j 
			+ numIntStates * (
				(i.x - solver->numGhost) 
<? if solver.dim > 1 then ?>
				+ (solver->gridSize.x - 2 * solver->numGhost) * (
					(i.y - solver->numGhost) 
<? if solver.dim > 2 then ?>
					+ (solver->gridSize.y - 2 * solver->numGhost) * (i.z - solver->numGhost)
<? end ?>
				)
<? end ?>
			)] = src[j + numStates * index];
	}
}
]])
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
		if self.verbose then
			print('calc_dU_dt begin')
		end
		if cmdline.printBufsInt then
			print'\nUBuf:' solver:printBuf(UBuf, nil, solver.eqn.numIntStates)
		end

		self.copyWithoutToWithGhostKernel(solver.solverBuf, solver.UBufObj, UBuf)
		
		if solver.checkNaNs then assert(solver:checkFinite(derivBufObj)) end
		solver:constrainU()
		solver:boundary()
		
		self.derivBufObj:fill()	-- fill with zero
		self.integrateCallback(self.derivBufObj)
		self.copyWithToWithoutGhostKernel(solver.solverBuf, dUdtBuf, self.derivBufObj)

		if cmdline.printBufsInt then
			print'\ndUdtBuf:' solver:printBuf(dUdtBuf, nil, solver.eqn.numIntStates)
		end		
		if solver.checkNaNs then assert(solver:checkFinite(dUdtBuf)) end
		if self.verbose then
			print('calc_dU_dt end')
		end
	end

	local volumeWithoutBorder = solver.volumeWithoutBorder
	local numRealsWithoutBorder = volumeWithoutBorder * solver.eqn.numIntStates
	
	local restart = cmdline.intBERestart or (args and args.restart) or 20

	local linearSolverArgs = {
		env = solver.app.env,
		x = self.krylov_xObj,
		count = numRealsWithoutBorder,
		epsilon = cmdline.intBEEpsilon or (args and args.epsilon) or 1e-10,
		--maxiter = 1000,
		restart = restart,
		maxiter = cmdline.intBEMaxIter or (args and args.maxiter) or restart * numRealsWithoutBorder,
		-- logging:
		errorCallback = function(residual, iter, x, rLenSq)
			self.lastResidual = residual
			self.lastIter = iter
			if self.verbose then
				print('t', solver.t, 'iter', iter, 'residual', residual)
			end
--if iter < numRealsWithoutBorder then return false end
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
		if self.verbose then
			print('linearSolverArgs.A begin')
		end
		local dUdt = self.krylov_dUdtObj
		--[[ evolve based on U(t) ... so this is solving our Krylov solver against an unchanging matrix, i.e. static linear solver.
		-- this is stable, but not backward-Euler
		calc_dU_dt(dUdt, self.krylov_bObj)
		--]]
		-- [[ evolve based on U(t+dt) ... this is a correct backward-Euler solver
-- TODO ... by the 2nd time this is called, UBuf is all zeroes ... hmm ...		
-- the first time this is good
-- the second time this is zeroes
-- because it's the residual duh
-- so conjgrad on the nearly-zero-residual produces nans
-- conjres and gmres just don't do anything it seems
		calc_dU_dt(dUdt, U)
		--]]
		
		--UNext = U - dt * calc_dU_dt(lastU)
		if cmdline.printBufsInt then
			print'\nU:' solver:printBuf(U, nil, solver.eqn.numIntStates)
			print('self.linearSolverDT', self.linearSolverDT)
		end
		self.linearSolver.args.mulAdd(UNext, U, dUdt.obj, -self.linearSolverDT)
		if cmdline.printBufsInt then
			print'\nUNext:' solver:printBuf(UNext, nil, solver.eqn.numIntStates)
		end
		if self.verbose then
			print('linearSolverArgs.A end')
		end
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
		-- TODO getting nil values, esp from smaller grid sizes, 
		-- i think due to reduce() requiring a minimum amount of data, and for small grids it is overflowing?  maybe
		assert(a)
		assert(b)
		assert(numRealsWithoutBorder)
		local d = oldDot(a,b)
		assert(d)
		-- sqrt?  cbrt for 3D?
		return d / math.sqrt(numRealsWithoutBorder)
	end
end

-- step contains integrating flux and source terms
-- but not post iterate
function BackwardEuler:integrate(dt, callback)
	if self.verbose then
		print('BackwardEuler:integrate begin')
		print('dt', dt)
	end
	local solver = self.solver
	-- this is a call to 'calcDeriv'
	self.integrateCallback = callback
	
	self.linearSolverDT = dt
	-- UBuf needs to be overwritten to pass on to the calcFluxDeriv
	-- (TODO make calcFluxDeriv accept a parameter)
	self.copyWithToWithoutGhostKernel(solver.solverBuf, self.krylov_bObj, solver.UBufObj)
	if cmdline.printBufsInt then
		print'\nself.krylov_bObj:' solver:printBuf(self.krylov_bObj, nil, solver.eqn.numIntStates)
	end
	self.copyWithToWithoutGhostKernel(solver.solverBuf, self.krylov_xObj, solver.UBufObj)
	self.linearSolver()
	self.copyWithoutToWithGhostKernel(solver.solverBuf, solver.UBufObj, self.krylov_xObj)

	solver:boundary()	
	if solver.checkNaNs then assert(solver:checkFinite(derivBufObj)) end
	solver:constrainU()
	if self.verbose then
		print('BackwardEuler:integrate end')
	end
end

function BackwardEuler:updateGUI()
	ig.luatableTooltipInputFloatAsText('Krylov epsilon', self.linearSolver.args, 'epsilon')
	if self.linearSolver.args.restart then
		ig.luatableTooltipInputInt('GMRES restart', self.linearSolver.args, 'restart')
	end
	ig.luatableTooltipInputInt('Krylov maxiter', self.linearSolver.args, 'maxiter')	-- typically restart * number of reals = restart * numCells * number of states
	-- read-only:
	ig.igText('residual = '..self.lastResidual)	-- this is |r|
	ig.igText('iter = '..self.lastIter)
end

return BackwardEuler
