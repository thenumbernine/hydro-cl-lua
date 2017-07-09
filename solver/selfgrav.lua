local class = require 'ext.class'
local template = require 'template'
local Poisson = require 'solver.poisson'

-- TODO guarantee all potential values are initially positive (and from then on?)

local SelfGrav = class(Poisson)

SelfGrav.gravitationConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

-- params for solver/poisson.cl 
function SelfGrav:getCodeParams()
	return {
		calcRho = template([[
#define gravitationalConstant <?=clnumber(self.gravitationConstant)?>
	global <?=eqn.cons_t?>* U = UBuf + index;
	//maybe a 4pi?  or is that only in the continuous case?
	rho = gravitationalConstant * U->rho;
]], 
		{
			self = self,
			solver = self.solver,
			eqn = self.solver.eqn,
			clnumber = require 'clnumber',
		}),
	}
end

function SelfGrav:getPoissonCode()
	return template(
		[[

kernel void calcGravityDeriv(
	global <?=eqn.cons_t?>* derivBuffer,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
	
		real gravity = (UBuf[indexR].<?=self.potentialField?> - UBuf[indexL].<?=self.potentialField?>) / (2. * dx<?=side?>_at(i));

		deriv->m.s[side] -= U->rho * gravity;
		deriv->ETotal -= U->rho * gravity * U->m.s[side];
	}<? end ?>
}

//TODO just use the display var kernels
kernel void reduce_ePot(
	global real* reduceBuf,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	reduceBuf[index] = UBuf[index].ePot;
}

//hmm, if ePot is negative then we get the cool turbulence effect
//but if it is positive, esp >1, then we get no extra behavior ... a static sphere
//and if it is too big then it explodes
kernel void offsetPotentialAndAddToTotal(
	global <?=eqn.cons_t?>* UBuf,
	real ePotMin
) {
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?>* U = UBuf + index;
	U->ePot += 0. - ePotMin;
	U->ETotal += U->rho * U->ePot;
}

]],
	{
		self = self,
		solver = self.solver,
		eqn = self.solver.eqn,
	})
end

function SelfGrav:refreshSolverProgram()
	SelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	self.calcGravityDerivKernel = solver.solverProgram:kernel'calcGravityDeriv'
	self.calcGravityDerivKernel:setArg(1, solver.UBuf)

	--TODO just use the display var kernels
	self.reduce_ePotKernel = solver.solverProgram:kernel('reduce_ePot', solver.reduceBuf, solver.UBuf)
	
	self.offsetPotentialAndAddToTotalKernel = solver.solverProgram:kernel('offsetPotentialAndAddToTotal', solver.UBuf)
end

local ffi = require 'ffi'
function SelfGrav:resetState()
	SelfGrav.super.resetState(self)

	local solver = self.solver
	
	-- TODO easier way to reduce min/max on display var ePot
	--		but that means compiling the ePot dispaly var code even if the ePot display var is turned off ...
	-- TODO stop calling them 'display var' since they're not just used for displaying
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.reduce_ePotKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	local ePotMin = solver.reduceMin()
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.reduce_ePotKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	local ePotMax = solver.reduceMax()
	
	self.offsetPotentialAndAddToTotalKernel:setArg(1, ffi.new('real[1]', ePotMin))
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.offsetPotentialAndAddToTotalKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.reduce_ePotKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	local new_ePotMin = solver.reduceMin()
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.reduce_ePotKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	local new_ePotMax = solver.reduceMax()

print('offsetting potential energy from '..ePotMin..','..ePotMax..' to '..new_ePotMin..','..new_ePotMax)
end

local field = 'gravityPoisson'
local enableField = 'useGravity'
return setmetatable({
	class = SelfGrav,
}, {
	__call = function(reqWrapper, parent)
		local apply = reqWrapper.class:createBehavior(field, enableField)
		local template = apply(parent)

		function template:step(dt)
			template.super.step(self, dt)
			
			if not self[enableField] then return end
			self.integrator:integrate(dt, function(derivBuf)
				self[field]:relax()
				self[field].calcGravityDerivKernel:setArg(0, derivBuf)
				self.app.cmds:enqueueNDRangeKernel{kernel=self[field].calcGravityDerivKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
			end)
		end

		return template
	end,
})
