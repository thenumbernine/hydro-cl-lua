local class = require 'ext.class'
local Poisson = require 'solver.poisson'
local template = require 'template'

local NoDiv = class(Poisson)

NoDiv.potentialField = 'BPot'

function NoDiv:getCodeParams()
	return {
		calcRho = template([[
	real divB = .5 * (0
<? 
for j=0,solver.dim-1 do
?>		+ (U[stepsize.s<?=j?>].B.s<?=j?> 
			- U[-stepsize.s<?=j?>].B.s<?=j?>
		) / grid_dx<?=j?>
<? 
end 
?>	);
	//because this is the discrete case, no 4pi
	rho = divB;
]], 
		{
			eqn = self.solver.eqn,
			solver = self.solver,
		}),
	}
end

function NoDiv:getPoissonCode()
	return [[
kernel void noDiv(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
<? for j=0,solver.dim-1 do ?> 
	U->B.s<?=j?> -= (U[stepsize.s<?=j?>].BPot - U[-stepsize.s<?=j?>].BPot) / (2. * grid_dx<?=j?>);
<? end ?>
}

]]
end

function NoDiv:refreshSolverProgram()
	NoDiv.super.refreshSolverProgram(self)

	local solver = self.solver
	self.noDivKernel = solver.solverProgram:kernel('noDiv', solver.UBuf)
end

local field = 'noDivPoisson' 
local apply = NoDiv:createBehavior(field)
return function(parent)
	local template = apply(parent)
	
	function template:step(dt)
		template.super.step(self, dt)

		self[field]:relax()
		self.app.cmds:enqueueNDRangeKernel{kernel=self[field].noDivKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	end

	return template
end
