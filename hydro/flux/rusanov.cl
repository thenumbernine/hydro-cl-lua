//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: <?=solver_macros?> math <?=fluxFromCons?>

static inline  <?=cons_t?> <?=calcFluxForInterface?>(
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const * const UL,
	<?=cons_t?> const * const UR,
	<?=cell_t?> const * const cellL,
	<?=cell_t?> const * const cellR,
	real3 const xInt,
	<?=normal_t?> const n
) {
	<?=cons_t?> resultFlux;
	real lambdaMax;
	<? for _,lr in ipairs{"L", "R"} do ?>{
		<?=eqn:consWaveCodeMinMax{
			n = "n",
			U = "U"..lr,
			pt = "xInt",
			resultMin = "lambdaMin"..lr,
			resultMax = "lambdaMax"..lr,
			declare = true,
		}:gsub("\\*\n", "\\\n\t\t")?>
		lambdaMax = max(fabs(lambdaMin<?=lr?>), fabs(lambdaMax<?=lr?>));
	}<? end ?>
	<?=cons_t?> FL;
	<?=fluxFromCons?>(&FL, solver, UL, cellL, n);
	<?=cons_t?> FR;
	<?=fluxFromCons?>(&FR, solver, UR, cellR, n);
	for (int j = 0; j < numIntStates; ++j) {
		resultFlux.ptr[j] = .5 * (FL.ptr[j] + FR.ptr[j] - lambdaMax * ((UR)->ptr[j] - (UL)->ptr[j]));
	}
	return resultFlux;
}
