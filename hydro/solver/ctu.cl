//// MODULE_NAME: updateCTU
//// MODULE_DEPENDS: <?=solver_t?> <?=cell_t?> <?=cons_t?> SETBOUNDS cell_sqrt_det_g
<? if solver.usePLM then ?>
//// MODULE_DEPENDS: consLR_t
<? end ?>

/*
using the interface flux
update the LR states by integrating dt/2 times the interface flux
*/
kernel void updateCTU(
	constant <?=solver_t?> const * const solver,
	global <?=cell_t?> const * const cellBuf,
	global <?=solver.getULRArg?>,
	global <?=cons_t?> const * const fluxBuf,
	real const dt
) {
	SETBOUNDS(0,1);
	real3 const x = cellBuf[index].pos;
<? if eqn.weightFluxByGridVolume then ?>
	real const volume = cell_sqrt_det_g(solver, x);
<? else ?>
	real const volume = 1.<? for i=0,solver.dim-1 do ?> * solver->grid_dx.s<?=i?><? end ?>;
<? end ?>

	<?
for side=0,solver.dim-1 do
	?>{
		int const side = <?=side?>;
		
		int const indexIntL = side + dim * index;
		global <?=cons_t?> const * const fluxL = fluxBuf + indexIntL;
		
		int const indexIntR = indexIntL + dim * solver->stepsize.s<?=side?>;
		global <?=cons_t?> const * const fluxR = fluxBuf + indexIntR;
		
<? if eqn.weightFluxByGridVolume then ?>
		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real const volume_intL = cell_sqrt_det_g(solver, xIntL);
		real const areaL = volume_intL / solver->grid_dx.s<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
		real const volume_intR = cell_sqrt_det_g(solver, xIntR);
		real const areaR = volume_intR / solver->grid_dx.s<?=side?>;
<? else ?>
		real const areaL = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
		real const areaR = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
<? end ?>

		<?
	for side2=0,solver.dim-1 do
		if side2 ~= side then
		?>{
<? if solver.getULRArg == consLR_t.."* ULRBuf" then
?>			int const indexForSide = side + dim * index;
			global <?=consLR_t?> * const ULR = ULRBuf + indexForSide;

			for (int j = 0; j < numIntStates; ++j) {
				real const dF_dx = (
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL
				) / volume;

				ULR->L.ptr[j] -= .5 * dt * dF_dx;
				ULR->R.ptr[j] -= .5 * dt * dF_dx;
			}
<? elseif solver.getULRArg == cons_t.."* UBuf" then
?>			global <?=cons_t?> * const U = UBuf + index;

			for (int j = 0; j < numIntStates; ++j) {
				real const dF_dx = (
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
