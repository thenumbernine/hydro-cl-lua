local ffi = require 'ffi'
local class = require 'ext.class'
local Roe = require 'solver.roe'
local SelfGravitationBehavior = require 'solver.selfgrav'

local EulerRoe = class(SelfGravitationBehavior(Roe))

local ConvertToTex_EulerRoe_U = class(EulerRoe.ConvertToTex)

function ConvertToTex_EulerRoe_U:setArgs(kernel, var)
	kernel:setArg(1, ffi.new('int[1]', var.globalIndex))
	kernel:setArg(2, self.solver.UBuf)
	kernel:setArg(3, self.solver.ePotBuf)
end

function EulerRoe:addConvertToTexUBuf()
	self:addConvertToTex({
		name = 'U',
		type = 'cons_t',
		vars = assert(self.eqn.displayVars),
		extraArgs = {'const __global real* ePotBuf'},
		displayBody = self.eqn:getCalcDisplayVarCode(),
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
	self.addSourceKernel:setArg(2, self.ePotBuf)
end

return EulerRoe
