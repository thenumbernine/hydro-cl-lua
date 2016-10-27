local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'solver.roe'
local SRHDEqn = require 'eqn.srhd'

local SRHDRoe = class(Roe)

SRHDRoe.eqn = SRHDEqn()

function SRHDRoe:createBuffers()
	SRHDRoe.super.createBuffers(self)
	
	ffi.cdef(self.eqn:getTypeCode())
	local primSize = ffi.sizeof'prim_t'
	
	local ctx = self.app.ctx
	self.primBuf = ctx:buffer{rw=true, size=self.volume * self.dim * primSize}
end

function SRHDRoe:createDisplayVars()
	SRHDRoe.super.createDisplayVars(self)
	self:addDisplayVarSet('prim', table.unpack(self.eqn.primDisplayVars))
end

function SRHDRoe:refreshInitStateProgram()
	SRHDRoe.super.refreshInitStateProgram(self)
	self.initStateKernel:setArg(1, self.primBuf)
end

function SRHDRoe:getCalcDTCode() end

function SRHDRoe:refreshSolverProgram()
	-- createKernels in particular ...
	SRHDRoe.super.refreshSolverProgram(self)

	self.calcDTKernel:setArg(2, self.primBuf)
	self.calcEigenBasisKernel:setArg(3, self.primBuf)
	
	self.constrainUKernel = self.solverProgram:kernel('constrainU', self.UBuf)
	self.updatePrimsKernel = self.solverProgram:kernel('updatePrims', self.primBuf, self.UBuf)
end

function SRHDRoe:update()
	SRHDRoe.super.update(self)

	self.app.cmds:enqueueNDRangeKernel{kernel=self.constrainUKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueNDRangeKernel{kernel=self.updatePrimsKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

return SRHDRoe
