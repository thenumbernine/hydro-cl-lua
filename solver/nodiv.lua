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
	rho = divB;	//times 4 pi?
]],
	}
end

RemoveDivergence.extraCode = [[


]]

return RemoveDivergence:createBehavior'noDivPoisson' 
