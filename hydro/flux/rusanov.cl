//// MODULE_NAME: calcFluxForInterface
//// MODULE_DEPENDS: solver.macros math <?=eigen_forInterface?> <?=eqn_waveCode_depends?> <?=fluxFromCons?>

#define calcFluxForInterface(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */xInt,\
	/*normal_t const */n\
) {\
	<?=eigen_t?> eigInt;\
	<?=eigen_forInterface?>(&eigInt, solver, UL, UR, xInt, n);\
<?=eqn:eigenWaveCodePrefix("n", "&eigInt", "xInt"):gsub("\n", "\\\n")?>\
\
	real lambdaMax;\
	<? for _,U in ipairs{"UL", "UR"} do --\
		for _,minmax in ipairs{"Min", "Max"} do ?>{\
<?=eqn:consWaveCodePrefix("n", U, "xInt"):gsub("\t", "\t\t"):gsub("\n", "\\\n")?>\
		real const absLambda = fabs(<?=eqn["cons"..minmax.."WaveCode"](eqn, "n", U, "xInt")?>);\
		lambdaMax = max(lambdaMax, absLambda);\
	}<? end --\
	end --\
?>\
	<?=cons_t?> FL;\
	<?=fluxFromCons?>(&FL, solver, UL, xInt, n);\
	<?=cons_t?> FR;\
	<?=fluxFromCons?>(&FR, solver, UR, xInt, n);\
	for (int j = 0; j < numIntStates; ++j) {\
		(resultFlux)->ptr[j] = .5 * (FL.ptr[j] + FR.ptr[j] - lambdaMax * ((UR)->ptr[j] - (UL)->ptr[j]));\
	}\
}
