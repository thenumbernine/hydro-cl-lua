local class = require 'ext.class'
local PoissonSolver = require 'solver.poisson'

local GravityPotential = class(PoissonSolver)

GravityPotential.gravityConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

-- params for solver/poisson.cl 
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

function GravityPotential:refreshSolverProgram()
	GravityPotential.super.refreshSolverProgram(self)
	
	local solver = self.solver
	solver.calcGravityDerivKernel = solver.solverProgram:kernel'calcGravityDeriv'
	solver.calcGravityDerivKernel:setArg(1, solver.UBuf)
	solver.calcGravityDerivKernel:setArg(2, solver.ePotBuf)	
end

GravityPotential.extraCode = [[

__kernel void calcGravityDeriv(
	__global cons_t* derivBuffer,
	const __global cons_t* UBuf,
	const __global real* ePotBuf
) {
	SETBOUNDS(2,2);
	
	__global cons_t* deriv = derivBuffer + index;
	const __global cons_t* U = UBuf + index;

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

function GravityPotential:step(dt)
	local solver = self.solver
	if not solver.useGravity then return end
	solver.integrator:integrate(dt, function(derivBuf)
		for i=1,20 do
			solver:potentialBoundary()
			solver.app.cmds:enqueueNDRangeKernel{kernel=solver.solvePoissonKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
		end
		
		solver.calcGravityDerivKernel:setArg(0, derivBuf)
		solver.app.cmds:enqueueNDRangeKernel{kernel=solver.calcGravityDerivKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
	end)
end

local field = 'gravityPoisson'
local apply = GravityPotential:createBehavior(field, 'useGravity')
return function(parent)
	local template = apply(parent)

	function template:step(dt)
		template.super.step(self, dt)
		self[field]:step(dt)	
	end

	return template
end
