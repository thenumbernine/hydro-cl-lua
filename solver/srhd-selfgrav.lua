local class = require 'ext.class'
local template = require 'template'
local Poisson = require 'solver.poisson'

local SRHDSelfGrav = class(Poisson)

SRHDSelfGrav.gravitationConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

function SRHDSelfGrav:getPotBufType()
	return self.solver.eqn.prim_t
end

function SRHDSelfGrav:getPotBuf()
	return self.solver.primBuf
end

function SRHDSelfGrav:refreshBoundaryProgram()
	SRHDSelfGrav.super.refreshBoundaryProgram(self)
	local solver = self.solver
	solver.potentialBoundaryKernel:setArg(0, solver.primBuf)
end

-- params for solver/poisson.cl 
function SRHDSelfGrav:getCodeParams()
	return {
		-- because solver/poisson.cl assumes it's UBuf, 
		-- we gotta keep the name 'UBuf'
		-- even though it's the primBuf ...
		-- TODO if args overrides the same thing that getPotBufType does
		--  then just use that.
		args = 'global '..self.solver.eqn.prim_t..'* UBuf',
		calcRho = template([[
#define gravitationalConstant <?=clnumber(self.gravitationConstant)?>
	global <?=eqn.prim_t?>* prim = UBuf + index;
	//maybe a 4pi?  or is that only in the continuous case?
	rho = -gravitationalConstant * prim->rho;
]], 
		{
			self = self,
			solver = self.solver,
			eqn = self.solver.eqn,
			clnumber = require 'clnumber',
		}),
	}
end

function SRHDSelfGrav:getPoissonCode()
	return template(
		[[

kernel void calcGravityDeriv(
	global <?=eqn.cons_t?>* derivBuffer,
	global const <?=eqn.cons_t?>* UBuf,
	global const <?=eqn.prim_t?>* primBuf
) {
	SETBOUNDS(2,2);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	const global <?=eqn.prim_t?>* prim = primBuf + index;

	real3 dv_dt = _real3(0,0,0);
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];

		real gradient = (primBuf[indexR].<?=self.potentialField?> - primBuf[indexL].<?=self.potentialField?>) / (2. * dx<?=side?>_at(i));
		dv_dt.s<?=side?> = -gradient;
	}<? end ?>

	//W = 1/sqrt(1 - v^2)
	//W,t = -1/2 1/(1 - v^2)^(3/2) * (-2 v dot v,t)
	//W,t = W^3 (v dot v,t)
	real W = U->D / prim->rho;
	real dW_dt = W * W * W * real3_dot(prim->v, dv_dt);

	real h = 1. + heatCapacityRatio * prim->eInt;

	//u = W v
	//u,t = = (W v),t = W,t v + W v,t
	//S = rho h W^2 v = rho h W u
	//assuming rho and h are constant ... 
	//S,t = rho h (2 W W,t v + W^2 v,t)
	
	deriv->S = real3_sub(deriv->S,
		real3_add(		
			real3_scale(dv_dt, prim->rho * h * W * W),
			real3_scale(prim->v, prim->rho * h * 2. * W * dW_dt)
		)
	);
	
	//tau = rho h W^2 - p - rho W
	//where does the velocity even factor into this? 
	//W is all I can see
	//tau,t = rho h (2 W W,t) - rho W,t
	//tau,t = rho W,t (2 h W - 1)
	
	deriv->tau -= prim->rho * dW_dt * (2. * h * W - 1.);
}

]],
	{
		self = self,
		solver = self.solver,
		eqn = self.solver.eqn,
	})
end

function SRHDSelfGrav:refreshSolverProgram()
	SRHDSelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	solver.calcGravityDerivKernel = solver.solverProgram:kernel'calcGravityDeriv'
	solver.calcGravityDerivKernel:setArg(1, solver.UBuf)
	solver.calcGravityDerivKernel:setArg(2, solver.primBuf)
end

local field = 'gravityPoisson'
local enableField = 'useGravity'
return setmetatable({
	class = SRHDSelfGrav,
}, {
	__call = function(reqWrapper, parent)
		local apply = reqWrapper.class:createBehavior(field, enableField)
		local template = apply(parent)

		function template:step(dt)
			template.super.step(self, dt)
		
			self.app.cmds:enqueueNDRangeKernel{kernel=self.constrainUKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
			self.app.cmds:enqueueNDRangeKernel{kernel=self.updatePrimsKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
			
			if not self[enableField] then return end
			self.integrator:integrate(dt, function(derivBuf)
				self[field]:relax()
				self.calcGravityDerivKernel:setArg(0, derivBuf)
				self.app.cmds:enqueueNDRangeKernel{kernel=self.calcGravityDerivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
			end)
		end

		return template
	end,
})
