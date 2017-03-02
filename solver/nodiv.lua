local class = require 'ext.class'
local Poisson = require 'solver.poisson'

local NoDiv = class(Poisson)

NoDiv.potentialField = 'BPot'

NoDiv.gaussSeidelMaxIters = 20

function NoDiv:getCodeParams()
	return {
		args = 'global '..self.solver.eqn.cons_t..'* UBuf',
		calcRho = require 'template'([[
	global <?=eqn.cons_t?>* U = UBuf + index;
	real divB = .5 * (0
<? for j=0,solver.dim-1 do ?>
		+ (U[stepsize.s<?=j?>].B.s<?=j?> - U[-stepsize.s<?=j?>].B.s<?=j?>) / grid_dx<?=j?>
<? end ?>
	);
	//because this is the discrete case, no 4pi
	rho = divB;
]], {
	eqn = self.solver.eqn,
	solver = self.solver,
}),
	}
end

NoDiv.extraCode = [[

kernel void noDiv(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	global <?=eqn.cons_t?>* U = UBuf + index;
<? for j=0,solver.dim-1 do ?> 
	U->B.s<?=j?> -= (UBuf[stepsize.s<?=j?>].BPot - UBuf[-stepsize.s<?=j?>].BPot) / (2. * grid_dx<?=j?>);
<? end ?>
}

]]

function NoDiv:refreshSolverProgram()
	NoDiv.super.refreshSolverProgram(self)

	local solver = self.solver
	solver.noDivKernel = solver.solverProgram:kernel('noDiv', solver.UBuf)
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
