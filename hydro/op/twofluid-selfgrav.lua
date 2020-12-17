local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local SelfGrav = require 'hydro.op.selfgrav'

local TwoFluidSelfGrav = class(SelfGrav)

TwoFluidSelfGrav.name = 'twofluid_selfgrav'

function TwoFluidSelfGrav:getModuleDepends_Poisson()
	return table(TwoFluidSelfGrav.super.getModuleDepends_Poisson(self)):append{
		'cell_x',
	}
end

function TwoFluidSelfGrav:getPoissonDivCode()
	return [[
	source = 4. * M_PI * calc_rho_from_U(U)
		* solver->gravitationalConstant / unit_m3_per_kg_s2;	//'G'
]]
end

function TwoFluidSelfGrav:getPoissonCode()
	return self.solver.eqn:template([[
real3 <?=op.symbolPrefix?>_calcGravityAccel(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const x
) {
	real3 accel_g = real3_zero;
	
	<? for side=0,solver.dim-1 do ?>{
		// m/s^2
		accel_g.s<?=side?> = (
			U[solver->stepsize.s<?=side?>].<?=op.potentialField?>
			- U[-solver->stepsize.s<?=side?>].<?=op.potentialField?>
		) / (2. * cell_dx<?=side?>(x));
	}<? end ?>

	return accel_g;
}

kernel void <?=op.symbolPrefix?>_calcGravityDeriv(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuffer,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 const x = cell_x(i);
	
	global <?=cons_t?> * const deriv = derivBuffer + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
	real3 accel_g = <?=op.symbolPrefix?>_calcGravityAccel(solver, U, x);

	// kg/(m^2 s) = kg/m^3 * m/s^2
<? for _,fluid in ipairs(eqn.fluids) do	
?>	deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, real3_real_mul(accel_g, U-><?=fluid?>_rho));
<? end	
?>
	// kg/(m s^2) = (kg m^2 / s) * m/s^2
<? for _,fluid in ipairs(eqn.fluids) do	
?>	deriv-><?=fluid?>_ETotal -= real3_dot(U-><?=fluid?>_m, accel_g);
<? end	
?>
}

// this matches hydro/op/selfgrav.lua
kernel void <?=op.symbolPrefix?>_copyPotentialToReduce(
	constant <?=solver_t?> const * const solver,
	global real * const reduceBuf,
	global <?=cons_t?> const * const UBuf
) {
	SETBOUNDS(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}

// this matches hydro/op/selfgrav.lua
kernel void <?=op.symbolPrefix?>_offsetPotential(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	realparam const ePotMax
) {
	SETBOUNDS(0,0);
	global <?=cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> -= ePotMax;
}

]], {op=self})
end

return TwoFluidSelfGrav 
