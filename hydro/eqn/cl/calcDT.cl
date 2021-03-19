//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=normal_t?> <?=eqn_waveCode_depends?> <?=SETBOUNDS?> <?=cell_dx_i?>

<? if not require "hydro.solver.meshsolver":isa(solver) then ?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell\
) {\
	<? for side=0,solver.dim-1 do ?>{\
<? --\
if solver.coord.vectorComponent == "holonomic" --\
or require "hydro.coord.cartesian":isa(solver.coord) --\
then --\
?>		real const dx = solver->grid_dx.s<?=side?>;\
<? else --\
?>		real const dx = cell_dx<?=side?>((cell)->pos);\
<? end --\
?>		if (dx > 1e-7) {\
			<?=normal_t?> const n = normal_forSide<?=side?>((cell)->pos);\
			/* use cell-centered eigenvalues */\
			<?=eqn:consWaveCodePrefix("n", "U", "(cell)->pos"):gsub("\n", "\\\n\t\t\t")?>\
			real const lambdaMin = <?=eqn:consMinWaveCode("n", "U", "(cell)->pos")?>;\
			real const lambdaMax = <?=eqn:consMaxWaveCode("n", "U", "(cell)->pos")?>;\
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
	for (int i = 0; i < (cell)->faceCount; ++i) {\
		global <?=face_t?> const * const face = faces + cellFaceIndexes[i + (cell)->faceOffset];\
		real const dx = face->area;	/* face->cellDist? */\
		if (dx > 1e-7 && face->cells.x != -1 && face->cells.y != -1) {\
			<?=normal_t?> const n = normal_forFace(face);\
			/* all sides? or only the most prominent side? */\
			/* which should we pick eigenvalues from? */\
			/* use cell-centered eigenvalues */\
			<?=eqn:consWaveCodePrefix("n", "U", "(cell)->pos"):gsub("\n", "\\\n\t\t\t")?>\
			real const lambdaMin = <?=eqn:consMinWaveCode("n", "U", "(cell)->pos")?>;\
			real const lambdaMax = <?=eqn:consMaxWaveCode("n", "U", "(cell)->pos")?>;\
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
			absLambdaMax = max((real)1e-9, absLambdaMax);\
			*(dt) = (real)min(*(dt), dx / absLambdaMax);\
		}\
	}\
}

<? end -- meshsolver ?>
