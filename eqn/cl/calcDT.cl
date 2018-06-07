<? if require 'solver.gridsolver'.is(solver) then ?>

/*
run across each cell
*/

kernel void calcDT(
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
		//use cell-centered eigenvalues
		range_t lambda = calcCellMinMaxEigenvalues_<?=side?>(U, x); 
		lambda.min = (real)min((real)0., lambda.min);
		lambda.max = (real)max((real)0., lambda.max);
		
		dt = (real)min((real)dt, (real)(grid_dx<?=side?> / (fabs(lambda.max - lambda.min) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt;
}

<? else -- mesh solver ?>

kernel void calcDT(
	global real* dtBuf,					//[numCells]
	const global <?=eqn.cons_t?>* UBuf,	//[numCells]
	const global cell_t* cells,			//[numCells]
	const global iface_t* ifaces		//[numInterfaces]
) {
	int cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	const global <?=eqn.cons_t?>* U = UBuf + cellIndex;
	const global cell_t* cell = cells + cellIndex;

	real dt = INFINITY;
	for (int i = 0; i < cell->numSides; ++i) {
		const global iface_t* iface = ifaces + cell->ifaces[i];
		//all sides? or only the most prominent side?
		//which should we pick eigenvalues from?
		<? for side=0,solver.dim-1 do ?>{
			//use cell-centered eigenvalues
			range_t lambda = calcCellMinMaxEigenvalues_<?=side?>(U, cell->x); 
			lambda.min = min(0., lambda.min);
			lambda.max = max(0., lambda.max);
			
			dt = min(dt, (iface->dist / (fabs(lambda.max - lambda.min) + (real)1e-9)));
		}<? end ?>
	}
	dtBuf[cellIndex] = dt;
}

<? end -- mesh vs grid solver ?>
