return function(args)
	local parentClass
	if cmdline.noDivPoissonSolver then
		parentClass = require('op.poisson_'..cmdline.noDivPoissonSolver)
	else
		parentClass = args.poissonSolver or require 'op.poisson_krylov'
	end

	local class = require 'ext.class'
	local template = require 'template'

	local NoDiv = class(parentClass)

	-- which cons_t field to store the solved potential value in
	NoDiv.vectorField = 'B'
	NoDiv.potentialField = 'divBPot'
	NoDiv.chargeField = nil	-- nil means zero

	function NoDiv:init(args)
		self.scalar = args.scalar
		NoDiv.super.init(self, args)
		self.chargeField = args.chargeField
	end

	--[[
	template parameters forwarded back to getCode
	solve del^2 divBPot = delta . B for divBPot
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
				1. / solver->grid_dx.s<?=j?>
			)
		);
	<? 
	end 
	?>	
		
		source = <?=real_mul?>(source, .5);
		
		//because this is the discrete case, no 4pi
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
	so B' = B - grad divBPot
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
	]], {
			op = self,
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
end
