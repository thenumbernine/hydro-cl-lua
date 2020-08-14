
kernel void calcFluxAtCell(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = UBuf + index;
<? for side=0,solver.dim-1 do 
?>	fluxBuf[<?=side?> + dim * index] = fluxFromCons(solver, *U, x, normalInfo_forSide<?=side?>(x));
<? end
?>
}

kernel void calcDerivFiniteDifference(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* F = fluxBuf + dim * index;

	//2nd order.  TODO use n'th order.
	<? for j=0,solver.dim-1 do ?>{
		const global <?=eqn.cons_t?>* FL = F - dim * solver->stepsize.s<?=j?> + <?=j?>;
		const global <?=eqn.cons_t?>* FR = F + dim * solver->stepsize.s<?=j?> + <?=j?>;
		for (int k = 0; k < numIntStates; ++k) {
			deriv->ptr[k] -= (FR->ptr[k] - FL->ptr[k]) / (2. * solver->grid_dx.s<?=j?>);
		}
	}<? end ?>
}
