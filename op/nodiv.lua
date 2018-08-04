local class = require 'ext.class'
local Poisson = require 'op.poisson'
local template = require 'template'

local NoDiv = class(Poisson)

-- which cons_t field to store the solved potential value in
NoDiv.vectorField = 'B'
NoDiv.potentialField = 'BPot'
NoDiv.chargeField = nil	-- nil means zero

function NoDiv:init(args)
	NoDiv.super.init(self, args)
	self.chargeField = args.chargeField
end

-- template parameters forwarded back to getSolverCode
function NoDiv:getCalcRhoCode()
	if not self.chargeField then return end
	return template([[
	rho = <?=U?>-><?=poisson.chargeField?>; 
]], 
	{
		poisson = self,
		solver = self.solver,
	})
end

--[[
template parameters forwarded back to getSolverCode
solve del^2 BPot = delta . B for BPot
--]]
function NoDiv:getCalcRhoCode()
	return template([[
<?
local scalar = poisson.scalar
local zero = scalar..'_zero'
local add = scalar..'_add'
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'
?>
	<?=scalar?> divergence = <?=zero?>;
<? 
for j=0,solver.dim-1 do
?>	divergence = <?=add?>(
		divergence,
		<?=real_mul?>(
			<?=sub?>(
				U[stepsize.s<?=j?>].<?=poisson.vectorField?>.s<?=j?>,
				U[-stepsize.s<?=j?>].<?=poisson.vectorField?>.s<?=j?>
			),
			1. / grid_dx<?=j?>
		)
	);
<? 
end 
?>	
	
	divergence = <?=real_mul?>(divergence, .5);
	
	//because this is the discrete case, no 4pi
	rho = divergence;
<? if poisson.chargeField then ?>
	rho = <?=add?>(rho, U-><?=poisson.chargeField?>);
<? end ?>
]], 
	{
		poisson = self,
		solver = self.solver,
	})
end

--[[
subtract the gradient of the divergence potential from the vector field
so B' = B - grad BPot
so delta . B' = delta . B - delta . del^-2 delta . B = ...should be 0
--]]
function NoDiv:getPoissonCode()
	return template([[
<?
local scalar = poisson.scalar
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'
?>
kernel void noDiv<?=poisson.suffix?>(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
<? for j=0,solver.dim-1 do ?> 
	U-><?=poisson.vectorField?>.s<?=j?> = 
		<?=sub?>(
			U-><?=poisson.vectorField?>.s<?=j?>,
			<?=real_mul?>(
				<?=sub?>(
					U[stepsize.s<?=j?>].<?=poisson.potentialField?>,
					U[-stepsize.s<?=j?>].<?=poisson.potentialField?>
				), 1. / (2. * grid_dx<?=j?>)
			)
		);
<? end ?>
}

]], {
		poisson = self,
		solver = self.solver,
		eqn = self.solver.eqn,
	})
end

function NoDiv:refreshSolverProgram()
	NoDiv.super.refreshSolverProgram(self)
	local solver = self.solver
	self.noDivKernelObj = solver.solverProgramObj:kernel('noDiv'..self.suffix, solver.UBuf)
end

function NoDiv:step(dt)
	local solver = self.solver
	self:relax()
	self.noDivKernelObj()
end

return NoDiv
