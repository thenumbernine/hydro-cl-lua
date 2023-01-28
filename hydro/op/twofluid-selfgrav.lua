local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local SelfGrav = require 'hydro.op.selfgrav'

local TwoFluidSelfGrav = class(SelfGrav)

TwoFluidSelfGrav.name = 'twofluid_selfgrav'

function TwoFluidSelfGrav:getPoissonDivCode()
	return [[
	source = 4. * M_PI * calc_rho_from_U(U)
		* solver->gravitationalConstant / unit_m3_per_kg_s2;	//'G'
]]
end

function TwoFluidSelfGrav:getPoissonCode()
	return [[

//// MODULE_NAME: <?=calcGravityAccel?>
//// MODULE_DEPENDS: <?=cell_dx_i?>

#define <?=calcGravityAccel?>(\
	/*real3 * const */accel_g,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> * const */U,\
	/*real3 const */pt\
) {\
	*(accel_g) = {};\
\
	<? for side=0,solver.dim-1 do ?>{\
		/* m/s^2 */\
		(accel_g)->s<?=side?> = (\
			U[solver->stepsize.s<?=side?>].<?=op.potentialField?>\
			- U[-solver->stepsize.s<?=side?>].<?=op.potentialField?>\
		) / (2. * cell_dx<?=side?>(pt));\
	}<? end ?>\
}

//// MODULE_NAME: <?=calcGravityDeriv?>
//// MODULE_DEPENDS: <?=calcGravityAccel?>

kernel void <?=calcGravityDeriv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuffer,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const pt = cellBuf[index].pos;
	
	global <?=cons_t?> * const deriv = derivBuffer + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
	real3 accel_g;
	<?=calcGravityAccel?>(&accel_g, solver, U, pt);

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

//// MODULE_NAME: <?=copyPotentialToReduce?>

// this matches hydro/op/selfgrav.lua
kernel void <?=copyPotentialToReduce?>(
	constant <?=solver_t?> const * const solver,
	global real * const reduceBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}

//// MODULE_NAME: <?=offsetPotential?>
//// MODULE_DEPENDS: realparam

// this matches hydro/op/selfgrav.lua
kernel void <?=offsetPotential?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	realparam const ePotMax
) {
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> -= ePotMax;
}

]]
end

return TwoFluidSelfGrav 
