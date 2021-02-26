//// MODULE_NAME: calcFluxAtCell
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?> <?=solver_macros?> <?=fluxFromCons?>

kernel void calcFluxAtCell(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const fluxBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
<? for side=0,solver.dim-1 do 
?>	<?=fluxFromCons?>(fluxBuf + <?=side?> + dim * index, solver, U, cell, normal_forSide<?=side?>(cell->pos));
<? end
?>
}

//// MODULE_NAME: calcDerivFiniteDifference
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?> <?=solver_macros?>

kernel void calcDerivFiniteDifference(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const fluxBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const F = fluxBuf + dim * index;

	//2nd order.  TODO use n'th order.
	<? for j=0,solver.dim-1 do ?>{
		global <?=cons_t?> const * const FL = F - dim * solver->stepsize.s<?=j?> + <?=j?>;
		global <?=cons_t?> const * const FR = F + dim * solver->stepsize.s<?=j?> + <?=j?>;
		for (int k = 0; k < numIntStates; ++k) {
			deriv->ptr[k] -= (FR->ptr[k] - FL->ptr[k]) / (2. * solver->grid_dx.s<?=j?>);
		}
	}<? end ?>
}
