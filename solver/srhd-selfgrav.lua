local class = require 'ext.class'
local template = require 'template'
local Poisson = require 'solver.poisson'

local SRHDSelfGrav = class(Poisson)

SRHDSelfGrav.gravitationConstant = 2	---- 6.67384e-11 m^3 / (kg s^2)

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

	real3 du_dt = _real3(0,0,0);
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];

		real gradient = (primBuf[indexR].<?=self.potentialField?> - primBuf[indexL].<?=self.potentialField?>) / (2. * dx<?=side?>_at(i));
		du_dt.s<?=side?> = -gradient;
	}<? end ?>

	//u = W v
	real W = U->D / prim->rho;
	real3 u = real3_scale(prim->v, W);
	
	//u,t = W,t v + W v,t
	//W = 1/sqrt(1 - v^2)
	//W,t = -1/2 1/(1 - v^2)^(3/2) * (-2 v dot v,t)
	//W,t = W^3 (v dot v,t)
	//u,t = W (W^2 v (v dot v,t) + v,t)
	//u,t = W (u outer u + I) * v,t

	//(u outer u + I) (I - alpha u outer u)
	//= u outer u - alpha u^2 u outer u + I - alpha u outer u
	//= I + u outer u (1 - alpha (1 + u^2))
	//is the inverse for 1 - alpha (1 + u^2) = 0
	// i.e. alpha = 1 / (1 + u^2)
	//so the inverse of (I + u outer u) is (I - u outer u / (1 + u^2))

	//v,t = 1/W (I - u outer u / (1 + u^2)) u,t
	//v,t = u,t / W - u (u dot u,t) / (1 + u^2)
	real3 dv_dt = real3_add(
		real3_scale(du_dt, 1. / W),
		real3_scale(u, real3_dot(u, du_dt) / (1. + real3_dot(u, u)))
	);

	real dW_dt = W * W * W * real3_dot(prim->v, dv_dt);
	real h = 1. + heatCapacityRatio * prim->eInt;

	//why am I integrating negative again?

	//D = W rho
	//D,t = W,t rho
	deriv->D -= dW_dt * prim->rho;

	//S = rho h W^2 v = rho h W u
	//assuming rho and h are constant ... 
	//S,t = rho h (W,t u + W u,t)
	deriv->S = real3_sub(deriv->S,
		real3_add(		
			real3_scale(u, prim->rho * h * dW_dt),
			real3_scale(du_dt, prim->rho * h * W)
		)
	);
	
	//tau = rho h W^2 - p - rho W
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
