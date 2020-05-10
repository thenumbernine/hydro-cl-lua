return function(args)
	local parentClass
	if cmdline.noDivPoissonSolver then
		parentClass = require('op.poisson_'..cmdline.noDivPoissonSolver)
	else
		parentClass = args and args.poissonSolver or require 'op.poisson_krylov'
	end

	local class = require 'ext.class'
	local template = require 'template'

	local NoDiv = class(parentClass)
	NoDiv.name = 'NoDiv'

	-- which cons_t field to store the solved potential value in
	NoDiv.vectorField = 'B'
	NoDiv.potentialField = 'psi'
	NoDiv.chargeField = nil	-- nil means zero

	function NoDiv:init(args)
		self.scalar = args.scalar
		NoDiv.super.init(self, args)
		self.vectorField = args.vectorField
		self.potentialField = args.potentialField
		self.chargeField = args.chargeField
	end

	--[[
	template parameters forwarded back to getCode
	solve del^2 psi = delta . B for psi
	--]]
	function NoDiv:getPoissonDivCode()
		return template([[
<?
local solver = op.solver

local scalar = op.scalar
local zero = scalar..'_zero'
local add = scalar..'_add'
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'

for j=0,solver.dim-1 do
?>	source = <?=add?>(
		source,
		<?=real_mul?>(
			<?=sub?>(
				U[solver->stepsize.s<?=j?>].<?=op.vectorField?>.s<?=j?>,
				U[-solver->stepsize.s<?=j?>].<?=op.vectorField?>.s<?=j?>
			),
			.5 / solver->grid_dx.s<?=j?>
		)
	);
<? 
end 
?>	
	
<? if op.chargeField then 
?>	source = <?=add?>(source, U-><?=op.chargeField?>);
<? end 
?>
]], 
		{
			op = self,
		})
	end

	--[[
	subtract the gradient of the divergence potential from the vector field
	so B' = B - grad psi
	so delta . B' = delta . B - delta . del^-2 delta . B = ...should be 0
	--]]
	function NoDiv:getPoissonCode()
		return template([[
<?
local solver = op.solver
local eqn = solver.eqn

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

//TODO just use the display var kernels
kernel void copyPotentialToReduce<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global real* reduceBuf,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}
]], 	{
			op = self,
		})
	end

	function NoDiv:refreshSolverProgram()
		NoDiv.super.refreshSolverProgram(self)
		local solver = self.solver
		self.noDivKernelObj = solver.solverProgramObj:kernel('noDiv'..self.name, solver.solverBuf, solver.UBuf)
		self.copyPotentialToReduceKernelObj = solver.solverProgramObj:kernel('copyPotentialToReduce'..self.name, solver.solverBuf, solver.reduceBuf, solver.UBuf)
	end

	function NoDiv:step(dt)
		local solver = self.solver
		self:relax()
		self.noDivKernelObj()
	end

	return NoDiv
end
