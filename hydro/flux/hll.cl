//// MODULE_NAME: calcFluxForInterface
//// MODULE_DEPENDS: solver.macros math eigen_forInterface eqn.waveCode fluxFromCons
//HLL solver:

<?=eqn.cons_t?> calcFluxForInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> pUL,
	<?=eqn.cons_t?> pUR,
	real3 xInt,
	normal_t n
) {
	// get min/max lambdas of UL, UR, and interface U (based on Roe averaging)
	// TODO this in a more computationally efficient way
	<?=eqn.eigen_t?> eigInt = eigen_forInterface(solver, pUL, pUR, xInt, n);
	
	<?=eqn:eigenWaveCodePrefix('n', 'eigInt', 'xInt')?>
	real lambdaIntMin = <?=eqn:eigenMinWaveCode('n', 'eigInt', 'xInt')?>;
	real lambdaIntMax = <?=eqn:eigenMaxWaveCode('n', 'eigInt', 'xInt')?>;

<? if solver.flux.hllCalcWaveMethod == 'Davis direct' then ?>
	real sL = lambdaIntMin;
	real sR = lambdaIntMax;
<? end ?>
<? if solver.flux.hllCalcWaveMethod == 'Davis direct bounded' then ?>
	real lambdaLMin;
	{
		<?=eqn:consWaveCodePrefix('n', 'pUL', 'xInt')?>
		lambdaLMin = <?=eqn:consMinWaveCode('n', 'pUL', 'xInt')?>;
	}

	real lambdaRMax;
	{
		<?=eqn:consWaveCodePrefix('n', 'pUR', 'xInt')?>
		lambdaRMax = <?=eqn:consMaxWaveCode('n', 'pUR', 'xInt')?>;
	}
	
	real sL = min(lambdaLMin, lambdaIntMin);
	real sR = max(lambdaRMax, lambdaIntMax);
<? end ?>

	if (0 <= sL) {
		<?=eqn.cons_t?> FL = fluxFromCons(solver, pUL, xInt, n);
		return FL;
	} else if (sR <= 0) {
		<?=eqn.cons_t?> FR = fluxFromCons(solver, pUR, xInt, n);
		return FR;
	} else if (sL <= 0 && 0 <= sR) {
		<?=eqn.cons_t?> flux;
		<?=eqn.cons_t?> FL = fluxFromCons(solver, pUL, xInt, n);
		<?=eqn.cons_t?> FR = fluxFromCons(solver, pUR, xInt, n);
		for (int j = 0; j < numIntStates; ++j) {
			flux.ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * (pUR.ptr[j] - pUL.ptr[j])) / (sR - sL);
		}
		return flux;
	}
}
