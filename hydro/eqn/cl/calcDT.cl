<? if not require 'hydro.solver.meshsolver'.is(solver) then ?>
/*
run across each cell
*/

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);

	const global <?=eqn.cons_t?>* U = UBuf + index;

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		normalInfo_t n = normalInfo_forSide<?=side?>(x);
		//use cell-centered eigenvalues
		<?=eqn:consWaveCodePrefix('n', '*U', 'x')?>
		real lambdaMin = <?=eqn:consMinWaveCode('n', '*U', 'x')?>;
		real lambdaMax = <?=eqn:consMaxWaveCode('n', '*U', 'x')?>;
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);

<? 
if solver.coord.vectorComponent == 'cartesian' 
and not require 'hydro.coord.cartesian'.is(solver.coord)
then 
?>		real dx = cell_dx<?=side?>(x); 
<? else 
?>		real dx = solver->grid_dx.s<?=side?>;
<? end 
?>

		dt = (real)min(dt, dx / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt;
}

<? else -- meshsolver ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,					//[numCells]
	const global <?=eqn.cons_t?>* UBuf,	//[numCells]
	const global cell_t* cells,			//[numCells]
	const global face_t* faces,			//[numFaces]
	const global int* cellFaceIndexes	//[numCellFaceIndexes]
) {
	int cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	const global <?=eqn.cons_t?>* U = UBuf + cellIndex;
	const global cell_t* cell = cells + cellIndex;
	real3 x = cell->pos;

	real dt = INFINITY;
	for (int i = 0; i < cell->faceCount; ++i) {
		const global face_t* face = faces + cellFaceIndexes[i + cell->faceOffset];
		real dx = face->area;
		if (dx > 1e-7 && f->cells.x != -1 && f->cells.y != -1) {
			normalInfo_t n = normalInfo_forFace(face);
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
	dtBuf[cellIndex] = dt;
}

<? end -- mesh vs grid solver ?>
