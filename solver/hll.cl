kernel void calcFlux(
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - stepsize.s<?=side?>;
		real3 xL = xR;
		xL.s<?=side?> -= grid_dx<?=side?>;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		int indexInt = side + dim * index;
	
		<?=solver:getULRCode()?>
		
		// get min/max lambdas of UL, UR, and interface U (based on Roe averaging)
		// TODO this in a more computationally efficient way
		<?=eqn.eigen_t?> eigInt = eigen_forInterface(*UL, *UR, xInt, normalForSide<?=side?>());
		
		<?=eqn:eigenWaveCodePrefix(side, '&eigInt', 'xInt')?>
		real lambdaIntMin = <?=eqn:eigenMinWaveCode(side, '&eigInt', 'xInt')?>;
		real lambdaIntMax = <?=eqn:eigenMaxWaveCode(side, '&eigInt', 'xInt')?>;
	
#if 0	//Davis direct
		real sL = lambdaIntMin;
		real sR = lambdaIntMax;
#endif
#if 1	//David direct bounded
		range_t lambdaL = calcCellMinMaxEigenvalues_<?=side?>(UL, xL);
		range_t lambdaR = calcCellMinMaxEigenvalues_<?=side?>(UR, xR);
		real sL = min(lambdaL.min, lambdaIntMin);
		real sR = max(lambdaR.max, lambdaIntMax);
#endif

		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		if (0 <= sL) {
			<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
			*flux = FL;
		} else if (sR <= 0) {
			<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
			*flux = FR;
		} else if (sL <= 0 && 0 <= sR) {
			<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
			<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
			for (int j = 0; j < numIntStates; ++j) {
				flux->ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * (UR->ptr[j] - UL->ptr[j])) / (sR - sL);
			}
		}
	}<? end ?>
}
