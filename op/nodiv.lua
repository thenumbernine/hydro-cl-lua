local class = require 'ext.class'
local template = require 'template'

-- TODO make this a ctor parameter 
local Poisson = require(
	cmdline.noDivPoissonSolver 
	and 'op.poisson_'..cmdline.noDivPoissonSolver
	--or 'op.poisson_krylov'		-- Krylov
	or 'op.poisson_jacobi'		-- Jacobi
)

local NoDiv = class(Poisson)

-- which cons_t field to store the solved potential value in
NoDiv.vectorField = 'B'
NoDiv.potentialField = 'BPot'
NoDiv.chargeField = nil	-- nil means zero

function NoDiv:init(args)
	NoDiv.super.init(self, args)
	self.chargeField = args.chargeField
end

--[[
template parameters forwarded back to getCode
solve del^2 BPot = delta . B for BPot
--]]
function NoDiv:getPoissonDivCode()
	return template([[
<?
local solver = self.solver

local scalar = op.scalar
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
				U[solver->stepsize.s<?=j?>].<?=op.vectorField?>.s<?=j?>,
				U[-solver->stepsize.s<?=j?>].<?=op.vectorField?>.s<?=j?>
			),
			1. / solver->grid_dx.s<?=j?>
		)
	);
<? 
end 
?>	
	
	divergence = <?=real_mul?>(divergence, .5);
	
	//because this is the discrete case, no 4pi
	source = divergence;
<? if op.chargeField then ?>
	source = <?=add?>(source, U-><?=op.chargeField?>);
<? end ?>
]], 
	{
		op = self,
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
local scalar = op.scalar
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'
?>
kernel void noDiv<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
<? for j=0,solver.dim-1 do ?> 
	U-><?=op.vectorField?>.s<?=j?> = 
		<?=sub?>(
			U-><?=op.vectorField?>.s<?=j?>,
			<?=real_mul?>(
				<?=sub?>(
					U[solver->stepsize.s<?=j?>].<?=op.potentialField?>,
					U[-solver->stepsize.s<?=j?>].<?=op.potentialField?>
				), 1. / (2. * solver->grid_dx.s<?=j?>)
			)
		);
<? end ?>
}

]], {
		op = self,
		solver = self.solver,
		eqn = self.solver.eqn,
	})
end

function NoDiv:refreshSolverProgram()
	NoDiv.super.refreshSolverProgram(self)
	local solver = self.solver
	self.noDivKernelObj = solver.solverProgramObj:kernel('noDiv'..self.name, solver.solverBuf, solver.UBuf)
end

function NoDiv:step(dt)
	local solver = self.solver
	self:relax()
	self.noDivKernelObj()
end

return NoDiv
