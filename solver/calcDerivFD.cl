kernel void calcFluxAtCell(
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = UBuf + index;
<? for j=0,solver.dim-1 do 
?>	fluxBuf[<?=j?> + dim * index] = fluxFromCons_<?=j?>(*U, x);
<? end
?>
}

/*
euler equations:
rho_,t = -(rho v^j)_,j
(rho v^k)_,t = -(rho v^j v^k + g^jk P)_,j
E_,t = -(H u^j)_,j
*/
kernel void calcDerivFiniteDifference(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS(numGhost,numGhost);
	
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* F = fluxBuf + dim * index;

	//2nd order.  TODO put flux in a grid and use n'th order.
	<? for j=0,solver.dim-1 do ?>{
		const global <?=eqn.cons_t?>* FL = F - dim * stepsize.s<?=j?> + <?=j?>;
		const global <?=eqn.cons_t?>* FR = F + dim * stepsize.s<?=j?> + <?=j?>;
		for (int k = 0; k < numIntStates; ++k) {
			deriv->ptr[k] -= (FR->ptr[k] - FL->ptr[k]) / grid_dx<?=j?>;
		}
	}<? end ?>
}
