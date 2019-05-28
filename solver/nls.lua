--[[
from 2010 Colliander et al "Numerical Simulations ..."
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local GridSolver = require 'solver.gridsolver'

local NLSSolver = class(GridSolver)
NLSSolver.name = 'NonLinearSchrodinger'
NLSSolver.fixedDT = 1e-6
NLSSolver.eqnName = 'nls'

function NLSSolver:refreshSolverProgram()
	NLSSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivKernelObj.obj:setArg(2, self.UBuf)
end

function NLSSolver:refreshCalcDTKernel() end
function NLSSolver:calcDT() return self.fixedDT end

function NLSSolver:calcDeriv(derivBufObj, dt)
	self.calcDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivKernelObj()
end

return NLSSolver
