local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'


local PoissonSolver = class()

function PoissonSolver:init(solver)
	self.solver = solver
end

function PoissonSolver:createBuffers()
	local solver = self.solver
	solver:clalloc('ePotBuf', solver.volume * ffi.sizeof(solver.app.real))
end

function PoissonSolver:addConvertToTexs()
	self.solver:addConvertToTex{
		name = 'ePot',
		vars = {{['0'] = 'value = buf[index];'}},
	}
end

function PoissonSolver:getSolverCode()
	return require 'processcl'(
		file['solver/selfgrav.cl'], 
		table(self:getCodeParams(), {
			solver = self.solver,
		})
	)
end

function PoissonSolver:refreshSolverProgram()
	local solver = self.solver
	solver.initPotentialKernel = solver.solverProgram:kernel('initPotential', solver.ePotBuf, solver.UBuf)

	solver.solvePoissonKernel = solver.solverProgram:kernel('solvePoisson', solver.ePotBuf, solver.UBuf)
	
	solver.calcGravityDerivKernel = solver.solverProgram:kernel'calcGravityDeriv'
	solver.calcGravityDerivKernel:setArg(1, solver.UBuf)
	solver.calcGravityDerivKernel:setArg(2, solver.ePotBuf)	
end

function PoissonSolver:refreshBoundaryProgram()
	local solver = self.solver
	solver.potentialBoundaryProgram, solver.potentialBoundaryKernel =
		solver:createBoundaryProgramAndKernel{
			type = 'real',
			methods = table.map(solver.boundaryMethods, function(v,k)
				return solver.app.boundaryMethods[1+v[0]], k
			end),
		}
	solver.potentialBoundaryKernel:setArg(0, solver.ePotBuf)
end

function PoissonSolver:resetState()
	local solver = self.solver
	if solver.useGravity then 
		solver.app.cmds:enqueueNDRangeKernel{kernel=solver.initPotentialKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
		solver:potentialBoundary()
		for i=1,20 do
			solver.app.cmds:enqueueNDRangeKernel{kernel=solver.solvePoissonKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
			solver:potentialBoundary()
		end
	end
	
	-- TODO
	-- add potential energy into total energy
	-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
end

function PoissonSolver:step(dt)
	local solver = self.solver
	solver.integrator:integrate(dt, function(derivBuf)
		if solver.useGravity then
			for i=1,20 do
				solver:potentialBoundary()
				solver.app.cmds:enqueueNDRangeKernel{kernel=solver.solvePoissonKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
			end
		end
		
		solver.calcGravityDerivKernel:setArg(0, derivBuf)
		solver.app.cmds:enqueueNDRangeKernel{kernel=solver.calcGravityDerivKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
	end)
end

return PoissonSolver
