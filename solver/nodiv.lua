local class = require 'ext.class'
local Poisson = require 'solver.poisson'
local template = require 'template'

local NoDiv = class(Poisson)

-- which cons_t field to store the solved potential value in
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
	self.noDivKernelObj = solver.solverProgramObj:kernel('noDiv', solver.UBuf)
end

function NoDiv:step(dt)
	local solver = self.solver
	self:relax()
	self.noDivKernelObj()
end

return NoDiv
