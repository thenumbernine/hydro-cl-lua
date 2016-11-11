local ffi = require 'ffi'
local class = require 'ext.class'
local Roe = require 'solver.roe'
local SelfGravitationBehavior = require 'solver.selfgrav'

local EulerRoe = class(SelfGravitationBehavior(Roe))

EulerRoe.name = 'Euler Roe'

local ConvertToTex_EulerRoe_U = class(EulerRoe.ConvertToTex)

function ConvertToTex_EulerRoe_U:setArgs(kernel, var)
	kernel:setArg(1, self.solver.UBuf)
	kernel:setArg(2, self.solver.ePotBuf)
end

function EulerRoe:addConvertToTexUBuf()
	self:addConvertToTex({
		name = 'U',
		type = 'cons_t',
		extraArgs = {'const __global real* ePotBuf'},	-- the code is in eqn/euler.lua, so maybe the 'addConvertToTexUBuf' function should be too?
		varCodePrefix = self.eqn.displayVarCodePrefix,
		vars = self.eqn.displayVars,
	}, ConvertToTex_EulerRoe_U)
end

function EulerRoe:createEqn()
	self.eqn = require 'eqn.euler'(self)
end
	
function EulerRoe:refreshInitStateProgram()
	EulerRoe.super.refreshInitStateProgram(self)
	self.initStateKernel:setArg(1, self.ePotBuf)
end

function EulerRoe:refreshSolverProgram()
	EulerRoe.super.refreshSolverProgram(self)

	self.calcDTKernel:setArg(2, self.ePotBuf)
	self.calcEigenBasisKernel:setArg(3, self.ePotBuf)
end

return EulerRoe
