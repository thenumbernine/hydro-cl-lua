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
-- the index vs dstindex stuff is shared in common with the main display code
		displayCode = [[
__kernel void {name}(
	{input},
	int displayVar,
	const __global cons_t* UBuf,
	const __global prim_t* primBuf
) {
	SETBOUNDS(0,0);
	int dstindex = index;
	int4 dsti = i;

	//now constrain
	if (i.x < 2) i.x = 2;
	if (i.x > gridSize_x - 2) i.x = gridSize_x - 2;
#if dim >= 2
	if (i.y < 2) i.y = 2;
	if (i.y > gridSize_y - 2) i.y = gridSize_y - 2;
#endif
#if dim >= 3
	if (i.z < 2) i.z = 2;
	if (i.z > gridSize_z - 2) i.z = gridSize_z - 2;
#endif
	//and recalculate read index
	index = INDEXV(i);
	
	int side = 0;
	real value = 0;

	const __global cons_t* U = UBuf + index;
	const __global prim_t* prim = primBuf + index;
	
	switch (displayVar) {
	case display_U_D: value = U->D; break;
	case display_U_Sx: value = U->Sx; break;
	case display_U_Sy: value = U->Sy; break;
	case display_U_Sz: value = U->Sz; break;
	case display_U_S: value = sqrt(U->Sx*U->Sx + U->Sy*U->Sy + U->Sz*U->Sz); break;
	case display_U_tau: value = U->tau; break;
	case display_U_W: value = U->D / prim->rho; break;
	case display_U_primitive_reconstruction_error: 
		//prim have just been reconstructed from cons
		//so reconstruct cons from prims again and calculate the difference
		{
			cons_t U2 = consFromPrim(*prim);
			value = 0;
			value += fabs(U->D - U2.D);
			value += fabs(U->Sx - U2.Sx);
			value += fabs(U->Sy - U2.Sy);
			value += fabs(U->Sz - U2.Sz);
			value += fabs(U->tau - U2.tau);
		}
		break;
	}

{output}
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
	const __global prim_t* prim = buf + index;
	switch (displayVar) {
	case display_prim_rho: value = prim->rho; break;
	case display_prim_vx: value = prim->vx; break;
	case display_prim_vy: value = prim->vy; break;
	case display_prim_vz: value = prim->vz; break;
	case display_prim_v: value = sqrt(prim->vx*prim->vx + prim->vy*prim->vy + prim->vz*prim->vz); break;
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
	self.calcEigenBasisKernel:setArg(2, self.primBuf)
	if self.checkFluxError then
		self.calcEigenBasisKernel:setArg(3, self.fluxXformBuf)
	end

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
