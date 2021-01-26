return function(args)
	local parentClass
	if cmdline.noDivPoissonSolver then
		parentClass = require('hydro.op.poisson_'..cmdline.noDivPoissonSolver)
	else
		parentClass = args and args.poissonSolver or require 'hydro.op.poisson_krylov'
	end

	local class = require 'ext.class'

	local NoDiv = class(parentClass)
	NoDiv.name = 'NoDiv'

	-- which cons_t field to store the solved potential value in
	NoDiv.vectorField = 'B'
	NoDiv.potentialField = 'psi'
	NoDiv.chargeField = nil	-- nil means zero
	NoDiv.chargeCode = nil

	--[[
	args:
		scalar = scalar type, default 'real'
		
		vectorField = vector field name within UBuf 
		optionally you can skip vectorField's use by setting 
			readVectorField = function(offset) for returning code to read at offset from U pointer
			writeVectorField = function(dv) for returning code to subtract from ptr 'U' the value 'dv'

		potentialField = potential field name within UBuf
	
		chargeField = field of charge.  nil means zero.
		chargeCode = optionally, code for generating charge field value.  TODO just use type(chargeField)=='function'
	--]]
	function NoDiv:init(args)
		self.scalar = args.scalar
		NoDiv.super.init(self, args)
		self.vectorField = args.vectorField
		self.readVectorField = args.readVectorField	 -- or default below
		self.writeVectorField = args.writeVectorField	 -- or default below
		self.potentialField = args.potentialField
		self.chargeField = args.chargeField
		self.chargeCode = args.chargeCode
	end

	function NoDiv:getSymbolFields()
		return NoDiv.super.getSymbolFields(self):append{
			'noDiv',
			'copyPotentialToReduce',
		}
	end

	function NoDiv:readVectorField(offset)
		return 'U['..offset..'].'..self.vectorField
	end

	-- returns code that performs U->$v = U->$v - $dv
	function NoDiv:writeVectorField(dv)
		local scalar = self.scalar
		local UField = 'U->'..self.vectorField
		local sub = scalar..'3_sub'
		return UField..' = '..sub..'('..UField..', '..dv..');'
	end

	--[[
	template parameters forwarded back to getCode
	solve del^2 psi = delta . B for psi
	--]]
	function NoDiv:getPoissonDivCode()
		return self.solver.eqn:template([[
	if (<?=OOB?>(1,1)) {
		source = 0.;
	} else {
<?
local scalar = op.scalar
local zero = scalar..'_zero'
local add = scalar..'_add'
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'

for j=0,solver.dim-1 do
?>		source = <?=add?>(
			source,
			<?=real_mul?>(
				<?=sub?>(
					(<?=op:readVectorField(' solver->stepsize.s'..j)?>).s<?=j?>,
					(<?=op:readVectorField('-solver->stepsize.s'..j)?>).s<?=j?>
				),
				.5 / solver->grid_dx.s<?=j?>
			)
		);
<? 
end 
?>	
	
<? if op.chargeField then 
?>	source = <?=add?>(source, U-><?=op.chargeField?>);
<? elseif op.chargeCode then
?><?=op.chargeCode?><? 
end 
?>
	}
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
		-- to-be-templated by parent class
		return [[
//// MODULE_NAME: <?=noDiv?>
<?
local scalar = op.scalar
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'
?>
kernel void <?=noDiv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;

	real3 dv = real3_zero;
<? for j=0,solver.dim-1 do ?> 
	dv.s<?=j?> = <?=real_mul?>(
		<?=sub?>(
			U[solver->stepsize.s<?=j?>].<?=op.potentialField?>,
			U[-solver->stepsize.s<?=j?>].<?=op.potentialField?>
		), 1. / (2. * solver->grid_dx.s<?=j?>)
	);
<? end ?>
	
	<?=op:writeVectorField'dv'?>
}

//// MODULE_NAME: <?=copyPotentialToReduce?>

//TODO just use the display var kernels
kernel void <?=copyPotentialToReduce?>(
	constant <?=solver_t?> const * const solver,
	global real * const reduceBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}
]]
	end

	function NoDiv:initCodeModules()
		NoDiv.super.initCodeModules(self)
		local solver = self.solver
		solver.solverModulesEnabled[self.symbols.noDiv] = true
		solver.solverModulesEnabled[self.symbols.copyPotentialToReduce] = true
	end

	function NoDiv:refreshSolverProgram()
		NoDiv.super.refreshSolverProgram(self)
		local solver = self.solver
		self.noDivKernelObj = solver.solverProgramObj:kernel(self.symbols.noDiv, solver.solverBuf, solver.UBuf)
		self.copyPotentialToReduceKernelObj = solver.solverProgramObj:kernel(self.symbols.copyPotentialToReduce, solver.solverBuf, solver.reduceBuf, solver.UBuf)
	end

	function NoDiv:step(dt)
		local solver = self.solver
		self:relax()
		self.noDivKernelObj()
	end

	return NoDiv
end
