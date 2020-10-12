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


<? if not require 'hydro.solver.meshsolver'.is(solver) then ?>

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.getULRArg?>,
	realparam dt,	//not used by HLL, just making this match Roe / other FV solvers
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	typedef <?=eqn.cons_t?> cons_t;

	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cellBuf[index].pos;
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;

		real dx = solver->grid_dx.s<?=side?>;
		
		int indexL = index - solver->stepsize.s<?=side?>;
		real3 xL = xR;
		xL.s<?=side?> -= dx;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * dx;
		
		int indexInt = side + dim * index;
		global cons_t* flux = fluxBuf + indexInt;


<? if solver.coord.vectorComponent == 'holonomic'
	or require 'hydro.coord.cartesian'.is(solver.coord)
	then ?>
		real area = 1.<?
	for i=0,solver.dim-1 do
		if i ~= side then
			?> * solver->grid_dx.s<?=i?><?
		end
	end
?>;
<? else ?>
		real area = cell_area<?=side?>(xInt);
<? end ?>
		if (area <= 1e-7) {
			for (int j = 0; j < numStates; ++j) {
				flux->ptr[j] = 0;
			}
			return;
		}


		<?=solver:getULRCode():gsub('\n', '\n\t\t')?>
		
		cons_t pUL = cons_parallelPropagate<?=side?>(*UL, xL, .5 * dx);
		cons_t pUR = cons_parallelPropagate<?=side?>(*UR, xR, -.5 * dx);

		normal_t n = normal_forSide<?=side?>(xInt);

		*flux = calcFluxForInterface(solver, pUL, pUR, xInt, n);
	}<? end ?>
}

<? end -- mesh vs grid solver ?>
