//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: math <?=solver_macros?>
// Marquina flux:

#define <?=calcFluxForInterface?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */xInt,\
	/*<?=normal_t?> const */n\
) {\
	/* get flux at left and right cells */\
	<?=cons_t?> FL, FR;\
	<?=fluxFromCons?>(&FL, solver, UL, cellL, n);\
	<?=fluxFromCons?>(&FR, solver, UR, cellR, n);\
\
	/* get flux eigensystem for the left and right cells */\
	<?=eigen_t?> eigL, eigR;\
	<?=eigen_forCell?>(&eigL, solver, UL, cellL, n);\
	<?=eigen_forCell?>(&eigR, solver, UR, cellL, n);\
\
	/* transform */\
	<?=waves_t?> omegaL, omegaR, phiL, phiR;\
	<?=eigen_leftTransform?>(&omegaL, solver, &eigL, UL, xInt, n);\
	<?=eigen_leftTransform?>(&omegaR, solver, &eigR, UR, xInt, n);\
	<?=eigen_leftTransform?>(&phiL, solver, &eigL, &FL, xInt, n);\
	<?=eigen_leftTransform?>(&phiR, solver, &eigR, &FR, xInt, n);\
\
	/* I thought eqn:eigenWaveCodePrefix was more flexible, allowing optional suffixes to variables to-be-paired-with multiple eqn:eigenWaveCode calls ... oh well I guess I hvae to store them ... */\
	<?=waves_t?> lambdasL;\
	{\
		<?=eqn:eigenWaveCodePrefix{ --\
			n = "n", --\
			eig = "&eigL", --\
			pt = "xInt", --\
		}:gsub("\\*\n", "\\\n\t")?>\
		<? for k=0,eqn.numWaves-1 do ?>{\
			int const k = <?=k?>;\
			lambdasL.ptr[k] = <?=eqn:eigenWaveCode{ --\
				n = "n", --\
				eig = "&eigL", --\
				pt = "xInt", --\
				waveIndex = k, --\
			}?>;\
		}<? end ?>\
	}\
\
	<?=waves_t?> lambdasR;\
	{\
		<?=eqn:eigenWaveCodePrefix{ --\
			n = "n", --\
			eig = "&eigR", --\
			pt = "xInt", --\
		}:gsub("\\*\n", "\\\n\t")?>\
		<? for k=0,eqn.numWaves-1 do ?>{\
			int const k = <?=k?>;\
			lambdasR.ptr[k] = <?=eqn:eigenWaveCode{ --\
				n = "n", --\
				eig = "&eigR", --\
				pt = "xInt", --\
				waveIndex = k, --\
			}?>;\
		}<? end ?>\
	}\
\
	<?=waves_t?> phiPlus, phiMinus;\
	for (int k = 0; k < numWaves; ++k) {\
		if (lambdasL.ptr[k] * lambdasR.ptr[k] > 0) {\
			if (lambdasL.ptr[k] > 0) {\
				phiPlus.ptr[k] = phiL.ptr[k];\
				phiMinus.ptr[k] = 0;\
			} else {\
				phiPlus.ptr[k] = 0;\
				phiMinus.ptr[k] = phiR.ptr[k];\
			}\
		} else {\
			real const alpha = max(fabs(lambdasL.ptr[k]), fabs(lambdasR.ptr[k]));\
			phiPlus.ptr[k] = (phiL.ptr[k] + alpha * omegaL.ptr[k]) * .5;\
			phiMinus.ptr[k] = (phiR.ptr[k] - alpha * omegaR.ptr[k]) * .5;\
		}\
	}\
\
	<?=cons_t?> resultFluxL, resultFluxR;\
	<?=eigen_rightTransform?>(&resultFluxL, solver, &eigL, &phiPlus, xInt, n);\
	<?=eigen_rightTransform?>(&resultFluxR, solver, &eigR, &phiMinus, xInt, n);\
\
	for (int k = 0; k < numStates; ++k) {\
		(resultFlux)->ptr[k] = resultFluxL.ptr[k] + resultFluxR.ptr[k];\
	}\
}
