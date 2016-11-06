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
		displayCode = [[
__kernel void <?=name?>(
	<?=input?>,
	int displayVar,
	const __global cons_t* UBuf,
	const __global real* ePotBuf 
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

	cons_t U = UBuf[index];
	real ePot = ePotBuf[index];
	
]]..self.eqn:getCalcDisplayVarCode()..[[	

<?=output?>
}
]],
	}, ConvertToTex_EulerRoe_U)
end

function EulerRoe:createEqn()
	self.eqn = require 'eqn.euler'(self)
end
	
function EulerRoe:refreshInitStateProgram()
	EulerRoe.super.refreshInitStateProgram(self)
	self.initStateKernel:setArg(1, self.ePotBuf)
end

function EulerRoe:getCalcDTCode() end

function EulerRoe:refreshSolverProgram()
	EulerRoe.super.refreshSolverProgram(self)

	self.calcDTKernel:setArg(2, self.ePotBuf)
	self.calcEigenBasisKernel:setArg(3, self.ePotBuf)
	if self.checkFluxError then
		self.calcEigenBasisKernel:setArg(4, self.fluxXformBuf)
	end
end

return EulerRoe
