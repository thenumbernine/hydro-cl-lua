local class = require 'ext.class'
local template = require 'template'
local SelfGrav = require 'op.selfgrav'
local TwoFluidSelfGrav = class(SelfGrav)

TwoFluidSelfGrav.name = 'TwoFluidSelfGrav'

function TwoFluidSelfGrav:getPoissonDivCode()
	return [[
	source = solver->gravitationalConstant
		* unit_m3_per_kg_s2
		* (U->ion_rho + U->elec_rho);
]]
end

function TwoFluidSelfGrav:getPoissonCode()
	return template([[
<?
local solver = op.solver
local eqn = solver.eqn
local fluids = eqn.fluids
?>
kernel void calcGravityDeriv<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuffer,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
		
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - solver->stepsize.s<?=side?>;
		int indexR = index + solver->stepsize.s<?=side?>;

		// m/s^2
		real accel_g = (
			UBuf[indexR].<?=op.potentialField?> 
			- UBuf[indexL].<?=op.potentialField?>
		) / (2. * cell_dx<?=side?>(x));

		// kg/(m^2 s) = kg/m^3 * m/s^2
<? for _,fluid in ipairs(fluids) do	
?>		deriv-><?=fluid?>_m.s[side] -= U-><?=fluid?>_rho * accel_g;
<? end	
?>
		// kg/(m s^2) = (kg m^2 / s) * m/s^2
<? for _,fluid in ipairs(fluids) do	
?>		deriv-><?=fluid?>_ETotal -= U-><?=fluid?>_m.s[side] * accel_g;
<? end	
?>
	}<? end ?>
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

//keep energy positive
kernel void offsetPotentialAndAddToTotal<?=op.name?>(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	realparam ePotMin
) {
	const real basePotential = 0.;
	//const real basePotential = 1.;
	
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> += basePotential - ePotMin;
<? for _,fluid in ipairs(fluids) do
?>	U-><?=fluid?>_ETotal += U-><?=fluid?>_rho * U-><?=op.potentialField?>;
<? end
?>
}

]], {op=self})
end

return TwoFluidSelfGrav 
