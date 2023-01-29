//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: <?=solver_macros?> math <?=eigen_forInterface?> <?=fluxFromCons?>
//HLL solver:

static inline <?=cons_t?> <?=calcFluxForInterface?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & UL,
	<?=cons_t?> const & UR,
	<?=cell_t?> const & cellL,
	<?=cell_t?> const & cellR,
	real3 const xInt,
	<?=normal_t?> const n
) {
	<?=cons_t?> resultFlux;
	/* 
	get min/max lambdas of UL, UR, and interface U (based on Roe averaging)
	TODO this in a more computationally efficient way
	*/
	<?=eigen_t?> eigInt = <?=eigen_forInterface?>(solver, UL, UR, cellL, cellR, xInt, n);
	<?=eqn:eigenWaveCodeMinMax{
		n = "n",
		eig = "eigInt",
		pt = "xInt",
		resultMin = "lambdaIntMin",
		resultMax = "lambdaIntMax",
		declare = true,
	}:gsub("\n", "\n\t")?>

<? if solver.flux.hllCalcWaveMethod == "Davis direct" then ?>
	real const sL = lambdaIntMin;
	real const sR = lambdaIntMax;
<? end ?>

<? if solver.flux.hllCalcWaveMethod == "Davis direct bounded" then ?>

	real lambdaLMin;
	{
		<?=eqn:consWaveCodeMinMax{
			n = "n",
			U = "UL",
			pt = "cellL.pos",
			resultMin = "lambdaLMin",
		}:gsub("\n", "\n\t\t")?>;
	}

	real lambdaRMax;
	{
		<?=eqn:consWaveCodeMinMax{
			n = "n",
			U = "UR",
			pt = "cellR.pos",
			resultMax = "lambdaRMax",
		}:gsub("\n", "\n\t\t")?>;
	}

	real const sL = min(lambdaLMin, lambdaIntMin);
	real const sR = max(lambdaRMax, lambdaIntMax);
<? end ?>

	if (0 <= sL) {
		resultFlux = <?=fluxFromCons?>(solver, UL, cellL, n);
	} else if (sR <= 0) {
		resultFlux = <?=fluxFromCons?>(solver, UR, cellR, n);
	} else if (sL <= 0 && 0 <= sR) {
		<?=cons_t?> FL = <?=fluxFromCons?>(solver, UL, cellL, n);
		<?=cons_t?> FR = <?=fluxFromCons?>(solver, UR, cellR, n);
		for (int j = 0; j < numIntStates; ++j) {
			resultFlux.ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * (UR.ptr[j] - UL.ptr[j])) / (sR - sL);
		}
	}
	return resultFlux;
}
