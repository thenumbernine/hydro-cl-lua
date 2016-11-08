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
	kernel:setArg(1, ffi.new('int[1]', var.globalIndex))
	kernel:setArg(2, self.solver.UBuf)
	kernel:setArg(3, self.solver.primBuf)
end

-- replace the U convertToTex with some custom code 
function SRHDRoe:addConvertToTexUBuf()
	self:addConvertToTex({
		name = 'U',
		type = 'cons_t',
		vars = assert(self.eqn.displayVars),
		extraArgs = {'const __global prim_t* primBuf'},
-- the index vs dstindex stuff is shared in common with the main display code
		displayBodyCode = [[
	cons_t U = buf[index];
	prim_t prim = primBuf[index];
	switch (displayVar) {
	case display_U_D: value = U.D; break;
	case display_U_Sx: value = U.S.x; break;
	case display_U_Sy: value = U.S.y; break;
	case display_U_Sz: value = U.S.z; break;
	case display_U_S: value = coordLen(U.S); break;
	case display_U_tau: value = U.tau; break;
	case display_U_W: value = U.D / prim.rho; break;
	case display_U_primitive_reconstruction_error: 
		//prim have just been reconstructed from cons
		//so reconstruct cons from prims again and calculate the difference
		{
			cons_t U2 = consFromPrim(prim);
			value = 0;
			value += fabs(U.D - U2.D);
			value += fabs(U.S.x - U2.S.x);
			value += fabs(U.S.y - U2.S.y);
			value += fabs(U.S.z - U2.S.z);
			value += fabs(U.tau - U2.tau);
		}
		break;
	}
]]
	}, ConvertToTex_SRHD_U)
end

function SRHDRoe:addConvertToTexs()
	SRHDRoe.super.addConvertToTexs(self)

	self:addConvertToTex{
		name = 'prim', 
		type = 'prim_t',
		vars = self.eqn.primDisplayVars,
		displayBodyCode = [[
	prim_t prim = buf[index];
	switch (displayVar) {
	case display_prim_rho: value = prim.rho; break;
	case display_prim_vx: value = prim.v.x; break;
	case display_prim_vy: value = prim.v.y; break;
	case display_prim_vz: value = prim.v.z; break;
	case display_prim_v: value = coordLen(prim.v); break;
	case display_prim_eInt: value = prim.eInt; break;
	case display_prim_P: value = calc_P(prim.rho, prim.eInt); break;
	case display_prim_h: value = calc_h(prim.rho, calc_P(prim.rho, prim.eInt), prim.eInt); break;
	}
]],
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
