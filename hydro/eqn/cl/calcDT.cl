<? if not require 'hydro.solver.meshsolver'.is(solver) then ?>
/*
run across each cell
*/

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cellBuf[index].pos;

	const global <?=eqn.cons_t?>* U = UBuf + index;

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
<? 
if solver.coord.vectorComponent == 'holonomic'
or require 'hydro.coord.cartesian'.is(solver.coord)
then 
?>		real dx = solver->grid_dx.s<?=side?>;
<? else 
?>		real dx = cell_dx<?=side?>(x); 
<? end 
?>
		if (dx > 1e-7) {
			normal_t n = normal_forSide<?=side?>(x);
			//use cell-centered eigenvalues
			<?=eqn:consWaveCodePrefix('n', '*U', 'x')?>
			real lambdaMin = <?=eqn:consMinWaveCode('n', '*U', 'x')?>;
			real lambdaMax = <?=eqn:consMaxWaveCode('n', '*U', 'x')?>;
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
			absLambdaMax = max((real)1e-9, absLambdaMax);
			dt = (real)min(dt, dx / absLambdaMax);
		}
	}<? end ?>
	dtBuf[index] = dt;
}

<? else -- meshsolver ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,					//[numCells]
	const global <?=eqn.cons_t?>* UBuf,	//[numCells]
	const global <?=solver.coord.cell_t?>* cells,			//[numCells]
	const global <?=solver.coord.face_t?>* faces,			//[numFaces]
	const global int* cellFaceIndexes	//[numCellFaceIndexes]
) {
	SETBOUNDS(0,0);
	const global <?=solver.coord.cell_t?>* cell = cells + index;
	real3 x = cell->pos;
	
	const global <?=eqn.cons_t?>* U = UBuf + index;

	real dt = INFINITY;
	for (int i = 0; i < cell->faceCount; ++i) {
		const global <?=solver.coord.face_t?>* face = faces + cellFaceIndexes[i + cell->faceOffset];
		real dx = face->area;	//face->cellDist?
		if (dx > 1e-7 && face->cells.x != -1 && face->cells.y != -1) {
			normal_t n = normal_forFace(face);
			//all sides? or only the most prominent side?
			//which should we pick eigenvalues from?
			//use cell-centered eigenvalues
			<?=eqn:consWaveCodePrefix('n', '*U', 'x')?>
			real lambdaMin = <?=eqn:consMinWaveCode('n', '*U', 'x')?>;
			real lambdaMax = <?=eqn:consMaxWaveCode('n', '*U', 'x')?>;
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
			absLambdaMax = max((real)1e-9, absLambdaMax);
			dt = (real)min(dt, dx / absLambdaMax);
		}
	}
	dtBuf[index] = dt;
}

<? end -- mesh vs grid solver ?>
