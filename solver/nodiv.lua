local class = require 'ext.class'
local PoissonSolver = require 'solver.poisson'

local RemoveDivergence = class(PoissonSolver)

function RemoveDivergence:getCodeParams()
	return {
		args = 'const __global cons_t* UBuf',
		calcRho = [[
	//TODO make this modular
	//4 pi G rho for gravity
	//div(B) for electromagnetism
	const __global cons_t* U = UBuf + index;
	real divB = .5 * ((U[stepsize.x].B.x - U[-stepsize.x].B.x) / grid_dx0
					+ (U[stepsize.y].B.y - U[-stepsize.y].B.y) / grid_dx1,
					+ (U[stepsize.z].B.z - U[-stepsize.z].B.z) / grid_dx2);
	rho = divB;
]],
	}
end

RemoveDivergence.extraCode = [[

__kernel void removeDiv(
	__global cons_t* UBuf,
	const __global real* ePotBuf
) {
	SETBOUNDS(2,2);
	__global cons_t* U = UBuf + index;
	const __global real* ePot = ePotBuf + index;
	U->ePot.x -= (ePot[stepsize.x] - ePot[-stepsize.x]) / (2. * grid_dx0);
	U->ePot.y -= (ePot[stepsize.y] - ePot[-stepsize.y]) / (2. * grid_dx1);
	U->ePot.z -= (ePot[stepsize.z] - ePot[-stepsize.z]) / (2. * grid_dx2);
}

]]

return RemoveDivergence:createBehavior'noDivPoisson' 
