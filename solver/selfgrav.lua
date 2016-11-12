local class = require 'ext.class'
local PoissonSolver = require 'solver.poisson'

local GravityPotential = class(PoissonSolver)

GravityPotential.gravityConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

-- params for solver/selfgrav.cl 
function GravityPotential:getCodeParams()
	return {
		args = 'const __global cons_t* UBuf',
		calcRho = '#define gravitationalConstant '..require 'clnumber'(self.gravityConstant)..'\n'..[[
	//TODO make this modular
	//4 pi G rho for gravity
	//div(E) for electromagnetism
	const __global cons_t* U = UBuf + index;
	rho = 4. * M_PI * gravitationalConstant * U->rho;
]],
	}
end

GravityPotential.extraCode = [[

__kernel void calcGravityDeriv(
	__global cons_t* derivBuffer,
	const __global cons_t* UBuf,
	const __global real* potentialBuf)
{
	SETBOUNDS(2,2);
	
	__global cons_t* deriv = derivBuffer + index;
	const __global cons_t* U = UBuf + index;

	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
	
		real gradient = (potentialBuf[indexR] - potentialBuf[indexL]) / (2. * dx<?=side?>_at(i));
		real gravity = -gradient;

		deriv->m.s[side] -= U->rho * gravity;
		deriv->ETotal -= U->rho * gravity * U->m.s[side];
	}<? end ?>
}
]]

return GravityPotential:createBehavior('gravityPoisson', 'useGravity')
