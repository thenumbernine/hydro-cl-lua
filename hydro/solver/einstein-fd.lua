local table = require 'ext.table'
local GridSolver = require 'hydro.solver.gridsolver'

local EinsteinFiniteDifferenceSolver = GridSolver:subclass()

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
EinsteinFiniteDifferenceSolver.numGhost = 3

EinsteinFiniteDifferenceSolver.name = 'EinsteinFiniteDifference'

function EinsteinFiniteDifferenceSolver:init(...)
	EinsteinFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function EinsteinFiniteDifferenceSolver:initCodeModules()
	EinsteinFiniteDifferenceSolver.super.initCodeModules(self)

	self.solverModulesEnabled[self.eqn.symbols.calcDeriv] = true
end

function EinsteinFiniteDifferenceSolver:refreshSolverProgram()
	EinsteinFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel(self.eqn.symbols.calcDeriv)
	self.calcDerivKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivKernelObj.obj:setArg(2, self.UBuf)
	self.calcDerivKernelObj.obj:setArg(3, self.cellBuf)
end

function EinsteinFiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivKernelObj()
end

return EinsteinFiniteDifferenceSolver
