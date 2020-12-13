//// MODULE_NAME: calcFluxForInterface
//// MODULE_DEPENDS: solver.macros math fluxLimiter <?=eqn_waveCode_depends?> <?=eigen_forInterface?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> <?=waves_t?>

<? if eqn.roeUseFluxFromCons then 
-- this was inline'd before I made the function into a giant macro, then I can't use the //// comments to inline anymore so TODO change the MODULE_ markup to handle /* */ instead/aswell?
?>
//// MODULE_DEPENDS: <?=fluxFromCons?>
<? end ?>

// Roe solver:

<? 
local useFlux = solver.fluxLimiter > 1 
	and flux.usesFluxLimiter -- just flux/roe.lua right now
?>

//TODO entropy fix ... for the Euler equations at least
#define calcFluxForInterface(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant solver_t const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */xInt,\
	/*normal_t const */n<? if useFlux then ?>,\
	/*realparam const */dt_dx,\
	/*<?=cons_t?> const * const */UL_L,\
	/*<?=cons_t?> const * const */UL_R,\
	/*<?=cons_t?> const * const */UR_L,\
	/*<?=cons_t?> const * const */UR_R,\
	/*real3 const */xIntL,\
	/*real3 const */xIntR<? end ?>\
) {\
	<?=eigen_t?> eig;\
	<?=eigen_forInterface?>(&eig, solver, UL, UR, xInt, n);\
\
<?=eqn:eigenWaveCodePrefix("n", "&eig", "xInt"):gsub("\n", "\\\n")?>\
\
	<?=waves_t?> fluxEig;\
<? if not eqn.roeUseFluxFromCons then --\
?>	<?=cons_t?> UAvg;\
	for (int j = 0; j < numIntStates; ++j) {\
		UAvg.ptr[j] = .5 * ((UL)->ptr[j] + (UR)->ptr[j]);\
	}\
\
	<?=eigen_leftTransform?>(&fluxEig, solver, &eig, &UAvg, xInt, n);\
<? end --\
?>\
	<?=cons_t?> deltaU;\
<? if useFlux then  --\
?>	<?=cons_t?> deltaUL, deltaUR;\
<? end  --\
?>\
	for (int j = 0; j < numStates; ++j) {\
		deltaU.ptr[j] = (UR)->ptr[j] - (UL)->ptr[j];\
<? if useFlux then  --\
?>		deltaUL.ptr[j] = (UR_L)->ptr[j] - (UL_L)->ptr[j];\
		deltaUR.ptr[j] = (UR_R)->ptr[j] - (UL_R)->ptr[j];\
<? end  --\
?>	}\
\
	<?=waves_t?> deltaUEig;\
	<?=eigen_leftTransform?>(&deltaUEig, solver, &eig, &deltaU, xInt, n);\
<? 	if useFlux then ?>\
	<?=eigen_t?> eigL;\
	<?=eigen_forInterface?>(&eigL, solver, UL_L, UR_L, xInt, n);\
	<?=eigen_t?> eigR;\
	<?=eigen_forInterface?>(&eigR, solver, UL_R, UR_R, xInt, n);\
	<?=waves_t?> deltaUEigL;\
	<?=eigen_leftTransform?>(&deltaUEigL, solver, &eigL, &deltaUL, xIntL, n);\
	<?=waves_t?> deltaUEigR;\
	<?=eigen_leftTransform?>(&deltaUEigR, solver, &eigR, &deltaUR, xIntR, n);\
<? 	end ?>\
\
	<? for j=0,eqn.numWaves-1 do ?>{\
		const int j = <?=j?>;\
		real lambda = <?=eqn:eigenWaveCode("n", "&eig", "xInt", j)?>;\
\
<? if not eqn.roeUseFluxFromCons then --\
?>		fluxEig.ptr[j] *= lambda;\
<? else --\
?>		fluxEig.ptr[j] = 0.;\
<? end --\
?>		real sgnLambda = lambda >= 0 ? 1 : -1;\
\
<? if useFlux then --\
?>		real rEig;\
		if (deltaUEig.ptr[j] == 0) {\
			rEig = 0;\
		} else {\
			if (lambda >= 0) {\
				rEig = deltaUEigL.ptr[j] / deltaUEig.ptr[j];\
			} else {\
				rEig = deltaUEigR.ptr[j] / deltaUEig.ptr[j];\
			}\
		}\
		real phi = fluxLimiter(rEig);\
<? end --\
?>\
		fluxEig.ptr[j] -= .5 * lambda * deltaUEig.ptr[j] * (sgnLambda\
<? if useFlux then  --\
?>			+ phi * (lambda * (dt_dx) - sgnLambda)\
<? end  --\
?>		);\
	}<? end ?>\
\
	<?=eigen_rightTransform?>(resultFlux, solver, &eig, &fluxEig, xInt, n);\
\
<? if eqn.roeUseFluxFromCons then  --\
-- TODO hmm, fluxFromCons vs eigen_fluxTransform using the 'eig' structure --\
-- fluxFromCons is using the left and right states to create their flux jacobian transform - applied to the left and right states to make the left and right flux vector --\
-- while eigen_fluxTransform would use the intermediate state to create the flux vector --\
?>\
	<?=cons_t?> FL;\
	<?=fluxFromCons?>(&FL, solver, UL, xInt, n);\
	<?=cons_t?> FR;\
	<?=fluxFromCons?>(&FR, solver, UR, xInt, n);\
\
	for (int j = 0; j < numIntStates; ++j) {\
		(resultFlux)->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);\
	}\
<? end --\
?>\
}
