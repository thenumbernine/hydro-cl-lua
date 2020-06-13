//HLL solver:

<?=eqn.cons_t?> calcFluxForInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> pUL,
	<?=eqn.cons_t?> pUR,
	real3 xInt,
	normalInfo_t n
) {
	// get min/max lambdas of UL, UR, and interface U (based on Roe averaging)
	// TODO this in a more computationally efficient way
	<?=eqn.eigen_t?> eigInt = eigen_forInterface(solver, pUL, pUR, xInt, n);
	
	<?=eqn:eigenWaveCodePrefix('n', 'eigInt', 'xInt')?>
	real lambdaIntMin = <?=eqn:eigenMinWaveCode('n', 'eigInt', 'xInt')?>;
	real lambdaIntMax = <?=eqn:eigenMaxWaveCode('n', 'eigInt', 'xInt')?>;

<? if solver.calcWaveMethod == 'Davis direct' then ?>
	real sL = lambdaIntMin;
	real sR = lambdaIntMax;
<? end ?>
<? if solver.calcWaveMethod == 'Davis direct bounded' then ?>
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


<? if not require 'hydro.solver.meshsolver'.is(solver) then ?>

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.getULRArg?>,
	realparam dt	//not used by HLL, just making this match Roe / other FV solvers
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - solver->stepsize.s<?=side?>;
		real3 xL = xR;
		xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		int indexInt = side + dim * index;
	
		<?=solver:getULRCode():gsub('\n', '\n\t\t')?>
		
		cons_t pUL = cons_parallelPropagate<?=side?>(*UL, xL, .5 * dx);
		cons_t pUR = cons_parallelPropagate<?=side?>(*UR, xR, -.5 * dx);

		normalInfo_t n = normalInfo_forSide<?=side?>(x);

		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		*flux = calcFluxForInterface(solver, pUL, pUR, xInt, n);
	}<? end ?>
}

<? end -- mesh vs grid solver ?>
