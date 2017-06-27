local class = require 'ext.class'
local BSSNOKFiniteDifferenceEquation = require 'eqn.bssnok-fd'
local Solver = require 'solver.solver'

local BSSNOKFiniteDifferenceSolver = class(Solver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOKFiniteDifferenceSolver'

function BSSNOKFiniteDifferenceSolver:createEqn(eqn)
	self.eqn = BSSNOKFiniteDifferenceEquation(self)
end

function BSSNOKFiniteDifferenceSolver:refreshSolverProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernel = self.solverProgram:kernel'calcDeriv'
	self.calcDerivKernel:setArg(1, self.UBuf)

	self.constrainUKernel = self.solverProgram:kernel('constrainU', self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:refreshInitStateProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshInitStateProgram(self)

	self.initConnUBarKernel = self.initStateProgram:kernel('init_connBarU', self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:resetState()
	BSSNOKFiniteDifferenceSolver.super.resetState(self)
	
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initConnUBarKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self:boundary()
	self.app.cmds:enqueueNDRangeKernel{kernel=self.constrainUKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:finish()
end

function BSSNOKFiniteDifferenceSolver:getCalcDTCode() end
function BSSNOKFiniteDifferenceSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceSolver:calcDT()
	return self.fixedDT
end

function BSSNOKFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

function BSSNOKFiniteDifferenceSolver:step(dt)
	BSSNOKFiniteDifferenceSolver.super.step(self, dt)

	self.app.cmds:enqueueNDRangeKernel{kernel=self.constrainUKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

return BSSNOKFiniteDifferenceSolver
