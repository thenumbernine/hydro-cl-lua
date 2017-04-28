local class = require 'ext.class'
local BSSNOKFiniteDifferenceEquation = require 'eqn.bssnok-fd'
local Solver = require 'solver.solver'

local BSSNOKFiniteDifferenceForwardEulerSolver = class(Solver)
BSSNOKFiniteDifferenceForwardEulerSolver.name = 'BSSNOKFiniteDifferenceForwardEulerSolver'

function BSSNOKFiniteDifferenceForwardEulerSolver:refreshSolverProgram()
	self.calcDerivKernel = self.solverProgram:kernel'calcDeriv'
	self.calcDerivKernel:setArg(1, self.UBuf)
end

function BSSNOKFiniteDifferenceForwardEulerSolver:createEqn(eqn)
	self.eqn = BSSNOKFiniteDifferenceEquation(self)
end

function BSSNOKFiniteDifferenceForwardEulerSolver:getCalcDTCode() end
function BSSNOKFiniteDifferenceForwardEulerSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceForwardEulerSolver:calcDT()
	return self.fixedDT
end

function BSSNOKFiniteDifferenceForwardEulerSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

return BSSNOKFiniteDifferenceForwardEulerSolver
