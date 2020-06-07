local class = require 'ext.class'
local template = require 'template'
local SelfGrav = require 'hydro.op.selfgrav'
local TwoFluidSelfGrav = class(SelfGrav)

TwoFluidSelfGrav.name = 'TwoFluidSelfGrav'

function TwoFluidSelfGrav:getPoissonDivCode()
	return [[
	source = 4. * M_PI * calc_rho_from_U(*U)
		* solver->gravitationalConstant / unit_m3_per_kg_s2;	//'G'
]]
end

function TwoFluidSelfGrav:getPoissonCode()
	return template([[
<?
local solver = op.solver
local eqn = solver.eqn
local fluids = eqn.fluids
?>

real3 calcGravityAccel<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global const <?=eqn.cons_t?>* U,
	real3 x
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

kernel void calcGravityDeriv<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuffer,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	real3 accel_g = calcGravityAccel<?=op.name?>(solver, U, x);

	// kg/(m^2 s) = kg/m^3 * m/s^2
<? for _,fluid in ipairs(fluids) do	
?>	deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, real3_real_mul(accel_g, U-><?=fluid?>_rho));
<? end	
?>
	// kg/(m s^2) = (kg m^2 / s) * m/s^2
<? for _,fluid in ipairs(fluids) do	
?>	deriv-><?=fluid?>_ETotal -= real3_dot(U-><?=fluid?>_m, accel_g);
<? end	
?>
}

// this matches hydro/op/selfgrav.lua
kernel void copyPotentialToReduce<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global real* reduceBuf,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}

// this matches hydro/op/selfgrav.lua
kernel void offsetPotential<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	realparam ePotMax
) {
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> -= ePotMax;
}

]], {op=self})
end

return TwoFluidSelfGrav 
