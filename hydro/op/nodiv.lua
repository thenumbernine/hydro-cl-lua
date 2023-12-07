return function(args)
	local parentClass
	if cmdline.noDivPoissonSolver then
		parentClass = require('hydro.op.poisson_'..cmdline.noDivPoissonSolver)
	else
		parentClass = args and args.poissonSolver or require 'hydro.op.poisson_krylov'
	end

	local NoDiv = parentClass:subclass()
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
		chargeCode = optionally, code for generating charge field value.  
			TODO just use type(chargeField)=='function'
			TODO TODO just override the default code generation function
			This is most convenient to implement, but prevents subclasses from overloading this function in parent and from reading this arg 
	--]]
	function NoDiv:init(args)
		self.scalar = args.scalar
		NoDiv.super.init(self, args)
		self.vectorField = args.vectorField
		self.readVectorField = args.readVectorField	 -- or default below
		self.writeVectorField = args.writeVectorField	 -- or default below
		self.potentialField = args.potentialField
		
		self.chargeField = args.chargeField			-- specify the field to use (nil means no charge)
		self.chargeCode = args.chargeCode		-- or alternatively just override the function and provide your own code
	end

	function NoDiv:getSymbolFields()
		return NoDiv.super.getSymbolFields(self):append{
			'noDiv',
			'copyPotentialToReduce',
		}
	end

	function NoDiv:readVectorField(offset, j)
		return 'U['..offset..'].'..self.vectorField..'.s'..j
	end

	-- returns code that performs U->$v = U->$v - $dv
	function NoDiv:writeVectorField(dv)
		local scalar = self.scalar
		local UField = 'U->'..self.vectorField
		local sub = scalar..'3_sub'
		return UField..' = '..sub..'('..UField..', '..dv..');'
	end

	function NoDiv:addChargeField()
		if self.chargeField then
			local add = self.scalar..'_add'
			return 'source = '..add..'(source, U->'..self.chargeField..');'
		elseif self.chargeCode then
			return self.chargeCode
		else
			return ''
		end
	end

	--[[
	template parameters forwarded back to getCode
	solve del^2 psi = delta . B for psi
	--]]
	function NoDiv:getPoissonDivCode()
		return self.solver.eqn:template([[
<?
local scalar = op.scalar
local zero = scalar..'_zero'
local add = scalar..'_add'
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'
?>
	if (<?=OOB?>(1,1)) {
		source = <?=zero?>;
	} else {
<?
for j=0,solver.dim-1 do
?>		source = <?=add?>(
			source,
			<?=real_mul?>(
				<?=sub?>(
					(<?=op:readVectorField(' solver->stepsize.s'..j, j)?>),
					(<?=op:readVectorField('-solver->stepsize.s'..j, j)?>)
				),
				.5 / solver->grid_dx.s<?=j?>
			)
		);
<? 
end 
?>	

<?=op:addChargeField()?>
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
<?
local scalar = op.scalar
local sub = scalar..'_sub'
local real_mul = scalar..'_real_mul'
?>

//// MODULE_NAME: <?=noDiv?>
kernel void <?=noDiv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	real3 const pt = cellBuf[index].pos;

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
		self.noDivKernelObj = solver.solverProgramObj:kernel(self.symbols.noDiv, solver.solverBuf, solver.UBuf, solver.cellBuf)
		self.copyPotentialToReduceKernelObj = solver.solverProgramObj:kernel(self.symbols.copyPotentialToReduce, solver.solverBuf, solver.reduceBuf, solver.UBuf)
	end

	function NoDiv:step(dt)
		local solver = self.solver
		self:relax()
		-- remove the divergence
		self.noDivKernelObj()
	end

	return NoDiv
end
