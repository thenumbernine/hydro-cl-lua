local ffi = require 'ffi'
local class = require 'ext.class'
local SelfGrav = require 'solver.selfgrav'

-- TODO make this a behavior so Roe can be swapped with RoeImplicitLinear
local EulerRoe = class(SelfGrav(require 'solver.roe'))
--local EulerRoe = class(SelfGrav(require 'solver.roe_implicit_linearized'))

local ConvertToTex_EulerRoe_U = class(EulerRoe.ConvertToTex)

function ConvertToTex_EulerRoe_U:setArgs(kernel, var)
	kernel:setArg(1, self.solver.UBuf)
	kernel:setArg(2, self.solver.ePotBuf)
end

function EulerRoe:addConvertToTexUBuf()
	self:addConvertToTex({
		name = 'U',
		type = self.eqn.cons_t,
		extraArgs = {'const global real* ePotBuf'},	-- the code is in eqn/euler.lua, so maybe the 'addConvertToTexUBuf' function should be too?
		varCodePrefix = self.eqn:getDisplayVarCodePrefix(),
		vars = self.eqn:getDisplayVars(),
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
