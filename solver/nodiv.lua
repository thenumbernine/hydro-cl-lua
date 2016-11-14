local class = require 'ext.class'
local Poisson = require 'solver.poisson'

local NoDiv = class(Poisson)

function NoDiv:getCodeParams()
	return {
		args = 'const __global cons_t* UBuf',
		calcRho = require 'processcl'([[
	const __global cons_t* U = UBuf + index;
	real divB = .5 * (
		(U[stepsize.x].B.x - U[-stepsize.x].B.x) / grid_dx0
<? if solver.dim > 1 then ?>	
		+ (U[stepsize.y].B.y - U[-stepsize.y].B.y) / grid_dx1
<? end
if solver.dim > 2 then ?>
		+ (U[stepsize.z].B.z - U[-stepsize.z].B.z) / grid_dx2
<? end ?>
	);
	//because this is the discrete case, no 4pi
	rho = divB;
]], {solver = self.solver}),
	}
end

NoDiv.extraCode = [[

__kernel void noDiv(
	__global cons_t* UBuf,
	const __global real* ePotBuf
) {
return;	
	SETBOUNDS(2,2);
	__global cons_t* U = UBuf + index;
	const __global real* ePot = ePotBuf + index;
<? for i=0,solver.dim-1 do ?> 
	U->B.s<?=i?> -= (ePot[stepsize.s<?=i?>] - ePot[-stepsize.s<?=i?>]) / (2. * grid_dx0);
<? end ?>
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
