// Roe solver:

//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>,
	const global <?=eqn.eigen_t?>* eigenBuf, 
	realparam dt
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / grid_dx<?=side?>;//dx<?=side?>_at(i);
		int indexL = index - stepsize.s<?=side?>;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		real3 xIntL = xInt;
		xIntL.s<?=side?> -= grid_dx<?=side?>;
		
		real3 xIntR = xInt;
		xIntR.s<?=side?> += grid_dx<?=side?>;
	
		<?=solver:getULRCode()?>
	
		int indexInt = side + dim * index;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		
		<?=eqn:eigenWaveCodePrefix(side, '*eig', 'xInt')?>

		<?=eqn.waves_t?> fluxEig;
<? if not eqn.roeUseFluxFromCons then ?>
		<?=eqn.cons_t?> UAvg;
		for (int j = 0; j < numIntStates; ++j) {
			UAvg.ptr[j] = .5 * (UL->ptr[j] + UR->ptr[j]);
		}
		fluxEig = eigen_leftTransform_<?=side?>(*eig, UAvg, xInt);
<? end ?>

		<?=eqn.cons_t?> deltaU;	
<? if solver.fluxLimiter > 1 then ?>
		int indexR2 = indexR + stepsize.s<?=side?>;
		int indexL2 = indexL - stepsize.s<?=side?>;
		<?=solver:getULRCode{indexL = 'indexL2', indexR = 'indexL', suffix='_L'}?>
		<?=solver:getULRCode{indexL = 'indexR', indexR = 'indexR2', suffix='_R'}?>	
			
		<?=eqn.cons_t?> deltaUL, deltaUR;
<? end ?>
		
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = UR->ptr[j] - UL->ptr[j];
<? if solver.fluxLimiter > 1 then ?>
			deltaUL.ptr[j] = UR_L->ptr[j] - UL_L->ptr[j];
			deltaUR.ptr[j] = UR_R->ptr[j] - UL_R->ptr[j];
<? end ?>	
		}
		
		<?=eqn.waves_t?> deltaUEig = eigen_leftTransform_<?=side?>(*eig, deltaU, xInt);
<? if solver.fluxLimiter > 1 then ?>
		const global <?=eqn.eigen_t?>* eigL = eig - dim * stepsize.s<?=side?>;
		const global <?=eqn.eigen_t?>* eigR = eig + dim * stepsize.s<?=side?>;
		<?=eqn.waves_t?> deltaUEigL = eigen_leftTransform_<?=side?>(*eigL, deltaUL, xIntL);
		<?=eqn.waves_t?> deltaUEigR = eigen_leftTransform_<?=side?>(*eigR, deltaUR, xIntR);
<? end ?>

		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real lambda = <?=eqn:eigenWaveCode(side, '*eig', 'xInt', j)?>;

<? if not eqn.roeUseFluxFromCons then ?>
			fluxEig.ptr[j] *= lambda;
<? else ?>
			fluxEig.ptr[j] = 0.;
<? end ?>
			real sgnLambda = lambda >= 0 ? 1 : -1;
		
<? if solver.fluxLimiter > 1 then ?>
			real rEig;
			if (deltaUEig.ptr[j] == 0) {
				rEig = 0;
			} else {
				if (lambda >= 0) {
					rEig = deltaUEigL.ptr[j] / deltaUEig.ptr[j];
				} else {
					rEig = deltaUEigR.ptr[j] / deltaUEig.ptr[j];
				}
			}
			real phi = fluxLimiter(rEig);
<? end ?>

			fluxEig.ptr[j] -= .5 * lambda * deltaUEig.ptr[j] * (sgnLambda
<? if solver.fluxLimiter > 1 then ?>
				+ phi * (lambda * dt_dx - sgnLambda)
<? end ?>			
			);
		}<? end ?>
		
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		*flux = eigen_rightTransform_<?=side?>(*eig, fluxEig, xInt);

<? if eqn.roeUseFluxFromCons then ?>
		
		real3 xL = xR;
		xL.s<?=side?> -= grid_dx<?=side?>;
		
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
		}
<? end ?>
	}<? end ?>
}
