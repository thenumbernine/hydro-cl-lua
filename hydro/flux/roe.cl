//// MODULE_NAME: calcFluxForInterface
//// MODULE_DEPENDS: solver.macros math eqn.waveCode fluxLimiter eigen_forInterface eigen_left/rightTransform
// Roe solver:

<? 
local useFlux = solver.fluxLimiter > 1 
	and flux.usesFluxLimiter -- just flux/roe.lua right now
?>

//TODO entropy fix ... for the Euler equations at least
<?=eqn.cons_t?> calcFluxForInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> pUL,
	<?=eqn.cons_t?> pUR,
	real3 xInt,
	normal_t n<? if useFlux then ?>,
	realparam dt_dx,
	<?=eqn.cons_t?> pUL_L,
	<?=eqn.cons_t?> pUL_R,
	<?=eqn.cons_t?> pUR_L,
	<?=eqn.cons_t?> pUR_R,
	real3 xIntL,
	real3 xIntR<? end ?>) {
	eigen_t eig = eigen_forInterface(solver, pUL, pUR, xInt, n);

	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'xInt'):gsub('\n', '\n\t\t')?>

	waves_t fluxEig;
<? if not eqn.roeUseFluxFromCons then 
?>	cons_t UAvg;
	for (int j = 0; j < numIntStates; ++j) {
		UAvg.ptr[j] = .5 * (pUL.ptr[j] + pUR.ptr[j]);
	}
	
	fluxEig = eigen_leftTransform(solver, eig, UAvg, xInt, n);
<? end
?>
	cons_t deltaU;	
<? if useFlux then 
?>	cons_t deltaUL, deltaUR;
<? end 
?>		
	for (int j = 0; j < numStates; ++j) {
		deltaU.ptr[j] = pUR.ptr[j] - pUL.ptr[j];
<? if useFlux then 
?>		deltaUL.ptr[j] = pUR_L.ptr[j] - pUL_L.ptr[j];
		deltaUR.ptr[j] = pUR_R.ptr[j] - pUL_R.ptr[j];
<? end 
?>	}

	waves_t deltaUEig = eigen_leftTransform(solver, eig, deltaU, xInt, n);
<? 	if useFlux then ?>
	eigen_t eigL = eigen_forInterface(solver, pUL_L, pUR_L, xInt, n);
	eigen_t eigR = eigen_forInterface(solver, pUL_R, pUR_R, xInt, n);
	waves_t deltaUEigL = eigen_leftTransform(solver, eigL, deltaUL, xIntL, n);
	waves_t deltaUEigR = eigen_leftTransform(solver, eigR, deltaUR, xIntR, n);
<? 	end ?>

	<? for j=0,eqn.numWaves-1 do ?>{
		const int j = <?=j?>;
		real lambda = <?=eqn:eigenWaveCode('n', 'eig', 'xInt', j)?>;

<? if not eqn.roeUseFluxFromCons then
?>		fluxEig.ptr[j] *= lambda;
<? else
?>		fluxEig.ptr[j] = 0.;
<? end
?>		real sgnLambda = lambda >= 0 ? 1 : -1;
	
<? if useFlux then
?>		real rEig;
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
<? if useFlux then 
?>			+ phi * (lambda * dt_dx - sgnLambda)
<? end 
?>		);
	}<? end ?>
	
	cons_t flux = eigen_rightTransform(solver, eig, fluxEig, xInt, n);

<? if eqn.roeUseFluxFromCons then 
-- TODO hmm, fluxFromCons vs eigen_fluxTransform using the 'eig' structure
-- fluxFromCons is using the left and right states to create their flux jacobian transform - applied to the left and right states to make the left and right flux vector
-- while eigen_fluxTransform would use the intermediate state to create the flux vector
?>
//// MODULE_DEPENDS: fluxFromCons
	cons_t FL = fluxFromCons(solver, pUL, xInt, n);
	cons_t FR = fluxFromCons(solver, pUR, xInt, n);

	for (int j = 0; j < numIntStates; ++j) {
		flux.ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
	}
<? end 
?>
	return flux;
}
