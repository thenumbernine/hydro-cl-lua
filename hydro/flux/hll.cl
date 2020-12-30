//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: <?=solver_macros?> math <?=eigen_forInterface?> <?=eqn_waveCode_depends?> <?=fluxFromCons?>
//HLL solver:

#define <?=calcFluxForInterface?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */xInt,\
	/*<?=normal_t?> const */n\
) {\
	/* \
	get min/max lambdas of UL, UR, and interface U (based on Roe averaging)\
	TODO this in a more computationally efficient way\
	*/\
	<?=eigen_t?> eigInt;\
	<?=eigen_forInterface?>(&eigInt, solver, UL, UR, xInt, n);\
\
	<?=eqn:eigenWaveCodePrefix("n", "&eigInt", "xInt"):gsub("\n", "\\\n")?>\
	real lambdaIntMin = <?=eqn:eigenMinWaveCode("n", "&eigInt", "xInt")?>;\
	real lambdaIntMax = <?=eqn:eigenMaxWaveCode("n", "&eigInt", "xInt")?>;\
\
<? if solver.flux.hllCalcWaveMethod == "Davis direct" then ?>\
	real sL = lambdaIntMin;\
	real sR = lambdaIntMax;\
<? end ?>\
<? if solver.flux.hllCalcWaveMethod == "Davis direct bounded" then ?>\
	real lambdaLMin;\
	{\
		<?=eqn:consWaveCodePrefix("n", "UL", "xInt"):gsub("\n", "\\\n")?>\
		lambdaLMin = <?=eqn:consMinWaveCode("n", "UL", "xInt")?>;\
	}\
\
	real lambdaRMax;\
	{\
		<?=eqn:consWaveCodePrefix("n", "UR", "xInt"):gsub("\n", "\\\n")?>\
		lambdaRMax = <?=eqn:consMaxWaveCode("n", "UR", "xInt")?>;\
	}\
\
	real sL = min(lambdaLMin, lambdaIntMin);\
	real sR = max(lambdaRMax, lambdaIntMax);\
<? end ?>\
\
	if (0 <= sL) {\
		<?=fluxFromCons?>(resultFlux, solver, UL, xInt, n);\
	} else if (sR <= 0) {\
		<?=fluxFromCons?>(resultFlux, solver, UR, xInt, n);\
	} else if (sL <= 0 && 0 <= sR) {\
		<?=cons_t?> FL;\
		<?=fluxFromCons?>(&FL, solver, UL, xInt, n);\
		<?=cons_t?> FR;\
		<?=fluxFromCons?>(&FR, solver, UR, xInt, n);\
		for (int j = 0; j < numIntStates; ++j) {\
			(resultFlux)->ptr[j] = (sR * FL.ptr[j] - sL * FR.ptr[j] + sL * sR * ((UR)->ptr[j] - (UL)->ptr[j])) / (sR - sL);\
		}\
	}\
}
