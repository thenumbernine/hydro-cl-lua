local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'solver.roe'

local SRHDRoe = class(Roe)

function SRHDRoe:createEqn()
	self.eqn = require 'eqn.srhd'(self)
end

function SRHDRoe:createBuffers()
	SRHDRoe.super.createBuffers(self)
	self:clalloc('primBuf', self.volume * self.dim * ffi.sizeof'prim_t')
end

local ConvertToTex_SRHD_U = class(SRHDRoe.ConvertToTex)

function ConvertToTex_SRHD_U:setArgs(kernel, var)
	kernel:setArg(1, self.solver.UBuf)
	kernel:setArg(2, self.solver.primBuf)
end

-- replace the U convertToTex with some custom code 
function SRHDRoe:addConvertToTexUBuf()
	self:addConvertToTex({
		name = 'U',
		type = 'cons_t',
		extraArgs = {'const __global prim_t* primBuf'},
-- the index vs dstindex stuff is shared in common with the main display code
		varCodePrefix = self.eqn.displayVarCodePrefix,
		vars = self.eqn.displayVars,
	}, ConvertToTex_SRHD_U)
end

function SRHDRoe:addConvertToTexs()
	SRHDRoe.super.addConvertToTexs(self)

	self:addConvertToTex{
		name = 'prim', 
		type = 'prim_t',
		varCodePrefix = self.eqn.primDisplayVarCodePrefix,
		vars = self.eqn.primDisplayVars,
	}
end

function SRHDRoe:refreshInitStateProgram()
	SRHDRoe.super.refreshInitStateProgram(self)
	self.initStateKernel:setArg(1, self.primBuf)
end

function SRHDRoe:refreshSolverProgram()
	-- createKernels in particular ...
	SRHDRoe.super.refreshSolverProgram(self)

	self.calcDTKernel:setArg(2, self.primBuf)
	self.calcEigenBasisKernel:setArg(2, self.primBuf)

	self.constrainUKernel = self.solverProgram:kernel('constrainU', self.UBuf)
	self.updatePrimsKernel = self.solverProgram:kernel('updatePrims', self.primBuf, self.UBuf)
end

function SRHDRoe:step(dt)
	SRHDRoe.super.step(self, dt)

	self.app.cmds:enqueueNDRangeKernel{kernel=self.constrainUKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueNDRangeKernel{kernel=self.updatePrimsKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

function SRHDRoe:boundary()
	-- U boundary
	SRHDRoe.super.boundary(self)
	-- prim boundary
	self.boundaryKernel:setArg(0, self.primBuf)
	SRHDRoe.super.boundary(self)
	self.boundaryKernel:setArg(0, self.UBuf)
end

return SRHDRoe
