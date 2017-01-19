local class = require 'ext.class'
local Poisson = require 'solver.poisson'

local SelfGrav = class(Poisson)

SelfGrav.gravityConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

-- params for solver/poisson.cl 
function SelfGrav:getCodeParams()
	return {
		args = 'const global '..self.solver.eqn.cons_t..'* UBuf',
		calcRho = require 'template'([[
#define gravitationalConstant <?=clnumber(self.gravityConstant)?>
	const global <?=eqn.cons_t?>* U = UBuf + index;
	rho = gravitationalConstant * U->rho;	//maybe a 4pi?  or is that only in the continuous case?
]], {
	self = self,
	eqn = self.solver.eqn,
	solver = self.solver,
	clnumber = require 'clnumber',
}),
	}
end

SelfGrav.extraCode = [[

kernel void calcGravityDeriv(
	global <?=eqn.cons_t?>* derivBuffer,
	const global <?=eqn.cons_t?>* UBuf,
	const global real* ePotBuf
) {
	SETBOUNDS(2,2);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
	
		real gradient = (ePotBuf[indexR] - ePotBuf[indexL]) / (2. * dx<?=side?>_at(i));
		real gravity = -gradient;

		deriv->m.s[side] -= U->rho * gravity;
		deriv->ETotal -= U->rho * gravity * U->m.s[side];
	}<? end ?>
}
]]

function SelfGrav:refreshSolverProgram()
	SelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	solver.calcGravityDerivKernel = solver.solverProgram:kernel'calcGravityDeriv'
	solver.calcGravityDerivKernel:setArg(1, solver.UBuf)
	solver.calcGravityDerivKernel:setArg(2, solver.ePotBuf)	
end

local field = 'gravityPoisson'
local enableField = 'useGravity'
local apply = SelfGrav:createBehavior(field, enableField)
return function(parent)
	local template = apply(parent)

	function template:step(dt)
		template.super.step(self, dt)
		
		if not self[enableField] then return end
		self.integrator:integrate(dt, function(derivBuf)
			self[field]:relax()
			self.calcGravityDerivKernel:setArg(0, derivBuf)
			self.app.cmds:enqueueNDRangeKernel{kernel=self.calcGravityDerivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		end)
	end

	return template
end
