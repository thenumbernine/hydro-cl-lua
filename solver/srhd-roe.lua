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

function SRHDRoe:addDisplayVarSets()
	SRHDRoe.super.addDisplayVarSets(self)
	
	self:addDisplayVarSet{
		name = 'prim', 
		vars = self.eqn.primDisplayVars,
		displayCode = [[
	const __global prim_t* prim = (const __global prim_t*)buf + index;
	switch (displayVar) {
	case display_prim_rho: value = prim->rho; break;
	case display_prim_vx: value = prim->vx; break;
	case display_prim_vy: value = prim->vy; break;
	case display_prim_vz: value = prim->vz; break;
	case display_prim_eInt: value = prim->eInt; break;
	case display_prim_P: value = calc_P(prim->rho, prim->eInt); break;
	case display_prim_h: value = calc_h(prim->rho, calc_P(prim->rho, prim->eInt), prim->eInt); break;
	}
]],
	}	
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
