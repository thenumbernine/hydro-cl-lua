// Roe solver:

//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.getULRArg?>,
	realparam dt
) {
	typedef <?=eqn.cons_t?> cons_t;
	typedef <?=eqn.eigen_t?> eigen_t;
	typedef <?=eqn.waves_t?> waves_t;
	
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		int indexL = index - solver->stepsize.s<?=side?>;
	
		real dx = solver->grid_dx.s<?=side?>;
		
		real3 xL = xR;
		xL.s<?=side?> -= dx;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * dx;
		
		//this is used for the flux limiter
		//should it be using the coordinate dx or the grid dx?
		//real dt_dx = dt / cell_dx<?=side?>(xInt);
		real dt_dx = dt / dx;

		real3 xIntL = xInt;
		xIntL.s<?=side?> -= dx;
		
		real3 xIntR = xInt;
		xIntR.s<?=side?> += dx;
	
<?=solver:getULRCode():gsub('\t', '\t\t')?>

		// here is where parallel propagator comes into play
		cons_t pUL = cons_parallelPropagate<?=side?>(*UL, .5 * dx, xL);
		cons_t pUR = cons_parallelPropagate<?=side?>(*UR, -.5 * dx, xR);

		int indexInt = side + dim * index;
		eigen_t eig = eigen_forInterface(solver, pUL, pUR, xInt, normalForSide<?=side?>());
<?=eqn:eigenWaveCodePrefix(side, 'eig', 'xInt'):gsub('\t', '\t\t')?>

		waves_t fluxEig;
<? if not eqn.roeUseFluxFromCons then 
?>		cons_t UAvg;
		for (int j = 0; j < numIntStates; ++j) {
			UAvg.ptr[j] = .5 * (pUL.ptr[j] + pUR.ptr[j]);
		}
		fluxEig = eigen_leftTransform_<?=side?>(solver, eig, UAvg, xInt);
<? end 
?>
		cons_t deltaU;	
<? if solver.fluxLimiter > 1 then 
?>		int indexR2 = indexR + solver->stepsize.s<?=side?>;
		int indexL2 = indexL - solver->stepsize.s<?=side?>;
		<?=solver:getULRCode{indexL = 'indexL2', indexR = 'indexL', suffix='_L'}?>
		<?=solver:getULRCode{indexL = 'indexR', indexR = 'indexR2', suffix='_R'}?>	
		
		// here is where parallel propagator comes into play
		cons_t pUL_L = cons_parallelPropagate<?=side?>(*UL_L, 1.5 * dx, xIntL);		//xIntL2?
		cons_t pUL_R = cons_parallelPropagate<?=side?>(*UL_R, .5 * dx, xIntL);
		cons_t pUR_L = cons_parallelPropagate<?=side?>(*UR_L, -.5 * dx, xIntR);
		cons_t pUR_R = cons_parallelPropagate<?=side?>(*UR_R, -1.5 * dx, xIntR);	// xIntR2?
			
		cons_t deltaUL, deltaUR;
<? end 
?>		
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = pUR.ptr[j] - pUL.ptr[j];
<? if solver.fluxLimiter > 1 then 
?>			deltaUL.ptr[j] = pUR_L.ptr[j] - pUL_L.ptr[j];
			deltaUR.ptr[j] = pUR_R.ptr[j] - pUL_R.ptr[j];
<? end 
?>		}
		
		waves_t deltaUEig = eigen_leftTransform_<?=side?>(solver, eig, deltaU, xInt);
<? if solver.fluxLimiter > 1 then 
?>		eigen_t eigL = eigen_forInterface(solver, pUL_L, pUR_L, xInt, normalForSide<?=side?>());
		eigen_t eigR = eigen_forInterface(solver, pUL_R, pUR_R, xInt, normalForSide<?=side?>());
		waves_t deltaUEigL = eigen_leftTransform_<?=side?>(solver, eigL, deltaUL, xIntL);
		waves_t deltaUEigR = eigen_leftTransform_<?=side?>(solver, eigR, deltaUR, xIntR);
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
		
		global cons_t* flux = fluxBuf + indexInt;
		*flux = eigen_rightTransform_<?=side?>(solver, eig, fluxEig, xInt);

<? if eqn.roeUseFluxFromCons then 
?>		cons_t FL = fluxFromCons_<?=side?>(solver, *UL, xL);
		cons_t FR = fluxFromCons_<?=side?>(solver, *UR, xR);
		
		cons_t pFL = cons_parallelPropagate<?=side?>(FL, .5 * dx, xL);
		cons_t pFR = cons_parallelPropagate<?=side?>(FR, -.5 * dx, xR);
		
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] += .5 * (pFL.ptr[j] + pFR.ptr[j]);
		}
<? end 
?>	}<? end ?>
}
