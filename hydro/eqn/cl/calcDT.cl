//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=normal_t?> <?=eqn_waveCode_depends?> <?=SETBOUNDS?>

<? if not require "hydro.solver.meshsolver".is(solver) then ?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell\
) {\
	real3 const x = (cell)->pos;\
	<? for side=0,solver.dim-1 do ?>{\
<? --\
if solver.coord.vectorComponent == "holonomic" --\
or require "hydro.coord.cartesian".is(solver.coord) --\
then --\
?>		real const dx = solver->grid_dx.s<?=side?>;\
<? else --\
?>		real const dx = cell_dx<?=side?>(x);\
<? end --\
?>		if (dx > 1e-7) {\
			<?=normal_t?> const n = normal_forSide<?=side?>(x);\
			/* use cell-centered eigenvalues */\
			<?=eqn:consWaveCodePrefix("n", "U", "x"):gsub("\n", "\\\n\t\t\t")?>\
			real const lambdaMin = <?=eqn:consMinWaveCode("n", "U", "x")?>;\
			real const lambdaMax = <?=eqn:consMaxWaveCode("n", "U", "x")?>;\
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
			absLambdaMax = max((real)1e-9, absLambdaMax);\
			*(dt) = (real)min(*(dt), dx / absLambdaMax);\
		}\
	}<? end ?>\
}

<? else -- meshsolver ?>
//// MODULE_DEPENDS: <?=face_t?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell,\
	/*global <?=face_t?> const * const */faces,		/* [numFaces] */\
	/*global int const * const */cellFaceIndexes	/* [numCellFaceIndexes] */\
) {\
	real3 const x = (cell)->pos;\
	for (int i = 0; i < (cell)->faceCount; ++i) {\
		global <?=face_t?> const * const face = faces + cellFaceIndexes[i + (cell)->faceOffset];\
		real const dx = face->area;	/* face->cellDist? */\
		if (dx > 1e-7 && face->cells.x != -1 && face->cells.y != -1) {\
			<?=normal_t?> const n = normal_forFace(face);\
			/* all sides? or only the most prominent side? */\
			/* which should we pick eigenvalues from? */\
			/* use cell-centered eigenvalues */\
			<?=eqn:consWaveCodePrefix("n", "U", "x"):gsub("\n", "\\\n\t\t\t")?>\
			real const lambdaMin = <?=eqn:consMinWaveCode("n", "U", "x")?>;\
			real const lambdaMax = <?=eqn:consMaxWaveCode("n", "U", "x")?>;\
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
			absLambdaMax = max((real)1e-9, absLambdaMax);\
			*(dt) = (real)min(*(dt), dx / absLambdaMax);\
		}\
	}\
}

<? end -- meshsolver ?>
