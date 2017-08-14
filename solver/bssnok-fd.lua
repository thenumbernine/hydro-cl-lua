local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local BSSNOKFiniteDifferenceEquation = require 'eqn.bssnok-fd'
local Solver = require 'solver.solver'

local xNames = table{'x', 'y', 'z'}

local BSSNOKFiniteDifferenceSolver = class(Solver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOKFiniteDifferenceSolver'

function BSSNOKFiniteDifferenceSolver:init(...)
	BSSNOKFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function BSSNOKFiniteDifferenceSolver:createEqn(eqn)
	self.eqn = BSSNOKFiniteDifferenceEquation(self)
end

function BSSNOKFiniteDifferenceSolver:refreshSolverProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernel = self.solverProgram:kernel'calcDeriv'
	self.calcDerivKernel:setArg(1, self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:refreshInitStateProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshInitStateProgram(self)
	self.initConnUBarKernel = self.initStateProgram:kernel('init_connBarU', self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:resetState()
	BSSNOKFiniteDifferenceSolver.super.resetState(self)
	
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initConnUBarKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	self:boundary()
	self.app.cmds:finish()
end

function BSSNOKFiniteDifferenceSolver:getCalcDTCode() end
function BSSNOKFiniteDifferenceSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceSolver:calcDT() return self.fixedDT end

function BSSNOKFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
end

return BSSNOKFiniteDifferenceSolver
