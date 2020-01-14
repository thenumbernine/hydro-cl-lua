// Roe solver:

//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.getULRArg?>,
	realparam dt
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		int indexL = index - solver->stepsize.s<?=side?>;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		
		//this is used for the flux limiter
		//should it be using the coordinate dx or the grid dx?
		//real dt_dx = dt / cell_dx<?=side?>(xInt);
		real dt_dx = dt / solver->grid_dx.s<?=side?>;

		real3 xIntL = xInt;
		xIntL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		
		real3 xIntR = xInt;
		xIntR.s<?=side?> += solver->grid_dx.s<?=side?>;
	
<?=solver:getULRCode():gsub('\t', '\t\t')?>
	
		int indexInt = side + dim * index;
		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, normalForSide<?=side?>());
		
<?=eqn:eigenWaveCodePrefix(side, 'eig', 'xInt'):gsub('\t', '\t\t')?>

		<?=eqn.waves_t?> fluxEig;
<? if not eqn.roeUseFluxFromCons then 
?>		<?=eqn.cons_t?> UAvg;
		for (int j = 0; j < numIntStates; ++j) {
			UAvg.ptr[j] = .5 * (UL->ptr[j] + UR->ptr[j]);
		}
		fluxEig = eigen_leftTransform_<?=side?>(solver, eig, UAvg, xInt);
<? end 
?>
		<?=eqn.cons_t?> deltaU;	
<? if solver.fluxLimiter > 1 then 
?>		int indexR2 = indexR + solver->stepsize.s<?=side?>;
		int indexL2 = indexL - solver->stepsize.s<?=side?>;
		<?=solver:getULRCode{indexL = 'indexL2', indexR = 'indexL', suffix='_L'}?>
		<?=solver:getULRCode{indexL = 'indexR', indexR = 'indexR2', suffix='_R'}?>	
			
		<?=eqn.cons_t?> deltaUL, deltaUR;
<? end 
?>		
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = UR->ptr[j] - UL->ptr[j];
<? if solver.fluxLimiter > 1 then 
?>			deltaUL.ptr[j] = UR_L->ptr[j] - UL_L->ptr[j];
			deltaUR.ptr[j] = UR_R->ptr[j] - UL_R->ptr[j];
<? end 
?>		}
		
		<?=eqn.waves_t?> deltaUEig = eigen_leftTransform_<?=side?>(solver, eig, deltaU, xInt);
<? if solver.fluxLimiter > 1 then 
?>		<?=eqn.eigen_t?> eigL = eigen_forInterface(solver, *UL_L, *UR_L, xInt, normalForSide<?=side?>());
		<?=eqn.eigen_t?> eigR = eigen_forInterface(solver, *UL_R, *UR_R, xInt, normalForSide<?=side?>());
		<?=eqn.waves_t?> deltaUEigL = eigen_leftTransform_<?=side?>(solver, eigL, deltaUL, xIntL);
		<?=eqn.waves_t?> deltaUEigR = eigen_leftTransform_<?=side?>(solver, eigR, deltaUR, xIntR);
<? end 
?>
		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real lambda = <?=eqn:eigenWaveCode(side, 'eig', 'xInt', j)?>;

<? if not eqn.roeUseFluxFromCons then 
?>			fluxEig.ptr[j] *= lambda;
<? else 
?>			fluxEig.ptr[j] = 0.;
<? end 
?>			real sgnLambda = lambda >= 0 ? 1 : -1;
		
<? if solver.fluxLimiter > 1 then 
?>			real rEig;
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
<? end 
?>
			fluxEig.ptr[j] -= .5 * lambda * deltaUEig.ptr[j] * (sgnLambda
<? if solver.fluxLimiter > 1 then 
?>				+ phi * (lambda * dt_dx - sgnLambda)
<? end 
?>			);
		}<? end ?>
		
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		*flux = eigen_rightTransform_<?=side?>(solver, eig, fluxEig, xInt);

<? if eqn.roeUseFluxFromCons then 
?>		
		real3 xL = xR;
		xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(solver, *UL, xL);
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(solver, *UR, xR);
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
		}
<? end 
?>	}<? end ?>
}
