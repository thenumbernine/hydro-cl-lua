--[[
from 2010 Colliander et al "Numerical Simulations ..."
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local Solver = require 'solver.solver'

local NLSSolver = class(Solver)
NLSSolver.name = 'NonLinearSchrodinger'
NLSSolver.fixedDT = 1e-6

function NLSSolver:createEqn(eqn)
	self.eqn = require 'eqn.nls'(self)
end

function NLSSolver:refreshSolverProgram()
	NLSSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernel = self.solverProgram:kernel'calcDeriv'
	self.calcDerivKernel:setArg(1, self.UBuf)
end

function NLSSolver:getCalcDTCode() end
function NLSSolver:refreshCalcDTKernel() end
function NLSSolver:calcDT() return self.fixedDT end

function NLSSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
end

return NLSSolver
