//// MODULE_NAME: <?=calcFluxForInterface?>
//// MODULE_DEPENDS: math <?=solver_macros?> <?=fluxLimiter?> <?=eqn_waveCode_depends?> <?=eigen_forInterface?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> <?=waves_t?>

<? if eqn.roeUseFluxFromCons then 
-- this was inline'd before I made the function into a giant macro, then I can't use the //// comments to inline anymore so TODO change the MODULE_ markup to handle /* */ instead/aswell?
?>
//// MODULE_DEPENDS: <?=fluxFromCons?>
<? end ?>

// Roe solver:

<? 
local useFluxLimiter = solver.fluxLimiter > 1 
	and flux.usesFluxLimiter -- just flux/roe.lua right now
?>

//TODO entropy fix ... for the Euler equations at least
#define <?=calcFluxForInterface?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */xInt,\
	/*<?=normal_t?> const */n<? if useFluxLimiter then ?>,\
	/*realparam const */dt_dx,\
	/*<?=cons_t?> const * const */UL_L,\
	/*<?=cons_t?> const * const */UL_R,\
	/*<?=cell_t?> const * const */cellL_L,\
	/*<?=cell_t?> const * const */cellL_R,\
	/*real3 const */xIntL,\
	/*<?=cons_t?> const * const */UR_L,\
	/*<?=cons_t?> const * const */UR_R,\
	/*<?=cell_t?> const * const */cellR_L,\
	/*<?=cell_t?> const * const */cellR_R,\
	/*real3 const */xIntR<? end ?>\
) {\
	<?=eigen_t?> eig;\
	<?=eigen_forInterface?>(&eig, solver, UL, UR, cellL, cellR, xInt, n);\
\
	<?=eqn:eigenWaveCodePrefix{ --\
		n = "n", --\
		eig = "&eig", --\
		pt = "xInt", --\
	}:gsub("\\*\n", "\\\n\t")?>\
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
<? if useFluxLimiter then  --\
?>	<?=cons_t?> deltaUL, deltaUR;\
<? end  --\
?>\
	for (int j = 0; j < numStates; ++j) {\
		deltaU.ptr[j] = (UR)->ptr[j] - (UL)->ptr[j];\
<? if useFluxLimiter then  --\
?>		deltaUL.ptr[j] = (UR_L)->ptr[j] - (UL_L)->ptr[j];\
		deltaUR.ptr[j] = (UR_R)->ptr[j] - (UL_R)->ptr[j];\
<? end  --\
?>	}\
\
	<?=waves_t?> deltaUEig;\
	<?=eigen_leftTransform?>(&deltaUEig, solver, &eig, &deltaU, xInt, n);\
<? 	if useFluxLimiter then ?>\
	<?=eigen_t?> eigL;\
	<?=eigen_forInterface?>(&eigL, solver, UL_L, UR_L, cellL_L, cellR_L, xIntL, n);\
	<?=eigen_t?> eigR;\
	<?=eigen_forInterface?>(&eigR, solver, UL_R, UR_R, cellL_R, cellR_R, xIntR, n);\
	<?=waves_t?> deltaUEigL;\
	<?=eigen_leftTransform?>(&deltaUEigL, solver, &eigL, &deltaUL, xIntL, n);\
	<?=waves_t?> deltaUEigR;\
	<?=eigen_leftTransform?>(&deltaUEigR, solver, &eigR, &deltaUR, xIntR, n);\
<? 	end ?>\
\
	<? for j=0,eqn.numWaves-1 do ?>{\
		int const j = <?=j?>;\
		real const lambda = <?=eqn:eigenWaveCode{ --\
			n = "n", --\
			eig = "&eig", --\
			pt = "xInt", --\
			waveIndex = j, --\
		}:gsub("\\*\n", "\\\n\t\t")?>;\
\
<? if not eqn.roeUseFluxFromCons then --\
?>		fluxEig.ptr[j] *= lambda;\
<? else --\
?>		fluxEig.ptr[j] = 0.;\
<? end --\
?>		real sgnLambda = lambda >= 0 ? 1 : -1;\
\
<? if useFluxLimiter then --\
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
		real phi = <?=fluxLimiter?>(rEig);\
<? end --\
?>\
		fluxEig.ptr[j] -= .5 * lambda * deltaUEig.ptr[j] * (sgnLambda\
<? if useFluxLimiter then  --\
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
-- which is F(U) --\
-- while eigen_fluxTransform would use the intermediate state to create the flux vector --\
-- which is dF/dU * dU/dx --\
-- It turns out dF/dU * U = F *ONLY FOR* the Euler fluid equations, and this is why lots of literature assumes it is true for all equations (when it's not).--\
-- But that means roeUseFluxFromCons==true's option is the correct one, but not ==false? --\
-- But Dullemond eqn 6.56 says that the FL and FR come from R Lambda L U = dF/dU * U calculations ... --\
-- ... which are *ONLY* equal to F in the case of the Euler fluid equations, and systems with linear flux jacobians, (NOT shallow-water). --\
-- So that would mean this part is wrong. --\
-- But maybe the description in Dullemond and other sources is wrong. --\
-- Because this value itself is diff'd across the cell, so it is a 'dF', so you wouldn't want to compute a second 'd' of it via R Lambda L, because that would give you 'd dF'. --\
-- so roeUseFluxFromCons==true is good. --\
?>\
	<?=cons_t?> FL;\
	<?=fluxFromCons?>(&FL, solver, UL, cellL, n);\
	<?=cons_t?> FR;\
	<?=fluxFromCons?>(&FR, solver, UR, cellR, n);\
\
	for (int j = 0; j < numIntStates; ++j) {\
		(resultFlux)->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);\
	}\
<? end --\
?>\
}
