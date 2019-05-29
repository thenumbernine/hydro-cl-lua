/*
using the interface flux
update the LR states by integrating dt/2 times the interface flux
*/
kernel void updateCTU(
	constant <?=solver.solver_t?>* solver,
	global <?=solver.getULRArg?>,
	const global <?=eqn.cons_t?>* fluxBuf,
	real dt
) {
	SETBOUNDS(0,1);
	real3 x = cell_x(i);
<? if eqn.weightFluxByGridVolume then ?>
	real volume = cell_volume(solver, x);
<? else ?>
	const real volume = 1.<? for i=0,solver.dim-1 do ?> * solver->grid_dx.s<?=i?><? end ?>;
<? end ?>

	<? 
for side=0,solver.dim-1 do 
	?>{
		const int side = <?=side?>;
		
		int indexIntL = side + dim * index;
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * solver->stepsize.s<?=side?>; 
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + indexIntR;
			
<? if eqn.weightFluxByGridVolume then ?>	
		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real volume_intL = cell_volume(solver, xIntL);
		real areaL = volume_intL / solver->grid_dx.s<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
		real volume_intR = cell_volume(solver, xIntR);
		real areaR = volume_intR / solver->grid_dx.s<?=side?>;
<? else ?>
		const real areaL = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
		const real areaR = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
<? end ?>

		<?
	for side2=0,solver.dim-1 do
		if side2 ~= side then
		?>{
<? if solver.getULRArg == eqn.consLR_t..'* ULRBuf' then
?>			int indexForSide = side + dim * index;
			global <?=eqn.consLR_t?>* ULR = ULRBuf + indexForSide;

			for (int j = 0; j < numIntStates; ++j) {
				real dF_dx = (
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL
				) / volume;

				ULR->L.ptr[j] -= .5 * dt * dF_dx;
				ULR->R.ptr[j] -= .5 * dt * dF_dx;
			}
<? elseif solver.getULRArg == eqn.cons_t..'* UBuf' then
?>			global <?=eqn.cons_t?>* U = UBuf + index;

			for (int j = 0; j < numIntStates; ++j) {
				real dF_dx = (
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL
				) / volume;
				U->ptr[j] -= .5 * dt * dF_dx;
			}
<? else error("can't handle getULRArg "..tostring(solver.getULRArg)) end
?>		}<? 
		end
	end ?>
	}<? 
end ?>
}
