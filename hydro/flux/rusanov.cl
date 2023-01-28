//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: <?=solver_macros?> math <?=fluxFromCons?>

static inline  <?=cons_t?> <?=calcFluxForInterface?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & UL,
	<?=cons_t?> const & UR,
	<?=cell_t?> const & cellL,
	<?=cell_t?> const & cellR,
	real3 const xInt,
	<?=normal_t?> const n
) {
	<?=cons_t?> result;
	real lambdaMax;
	<? for _,lr in ipairs{"L", "R"} do ?>{
		<?=eqn:consWaveCodeMinMax{
			n = "n",
			U = "&U"..lr,
			pt = "xInt",
			resultMin = "lambdaMin"..lr,
			resultMax = "lambdaMax"..lr,
			declare = true,
		}:gsub("\\*\n", "\\\n\t\t")?>
		lambdaMax = max(fabs(lambdaMin<?=lr?>), fabs(lambdaMax<?=lr?>));
	}<? end ?>
	<?=cons_t?> FL = <?=fluxFromCons?>(solver, UL, cellL, n);
	<?=cons_t?> FR = <?=fluxFromCons?>(solver, UR, cellR, n);
	for (int j = 0; j < numIntStates; ++j) {
		result.ptr[j] = .5 * (FL.ptr[j] + FR.ptr[j] - lambdaMax * ((UR)->ptr[j] - (UL)->ptr[j]));
	}
	return result;
}
