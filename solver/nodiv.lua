local class = require 'ext.class'
local Poisson = require 'solver.poisson'

local NoDiv = class(Poisson)

function NoDiv:getCodeParams()
	return {
		args = 'const __global cons_t* UBuf',
		calcRho = [[
	//TODO make this modular
	//4 pi G rho for gravity
	//div(B) for electromagnetism
	const __global cons_t* U = UBuf + index;
	real divB = .5 * ((U[stepsize.x].B.x - U[-stepsize.x].B.x) / grid_dx0
					+ (U[stepsize.y].B.y - U[-stepsize.y].B.y) / grid_dx1
					+ (U[stepsize.z].B.z - U[-stepsize.z].B.z) / grid_dx2);
	rho = -4. * M_PI * divB;
]],
	}
end

NoDiv.extraCode = [[

__kernel void noDiv(
	__global cons_t* UBuf,
	const __global real* ePotBuf
) {
	SETBOUNDS(2,2);
	__global cons_t* U = UBuf + index;
	const __global real* ePot = ePotBuf + index;
	U->B.x = -= (ePot[stepsize.x] - ePot[-stepsize.x]) / (2. * grid_dx0);
	U->B.y = -= (ePot[stepsize.y] - ePot[-stepsize.y]) / (2. * grid_dx1);
	U->B.z = -= (ePot[stepsize.z] - ePot[-stepsize.z]) / (2. * grid_dx2);
}

]]

function NoDiv:refreshSolverProgram()
	NoDiv.super.refreshSolverProgram(self)

	local solver = self.solver
	solver.noDivKernel = solver.solverProgram:kernel('noDiv', solver.UBuf, solver.ePotBuf)
end

local field = 'noDivPoisson' 
local apply = NoDiv:createBehavior(field)
return function(parent)
	local template = apply(parent)
	
	function template:step(dt)
		template.super.step(self, dt)

		self[field]:relax()
		self.app.cmds:enqueueNDRangeKernel{kernel=self.noDivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	end

	return template
end
