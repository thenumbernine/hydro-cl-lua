kernel void calcDT(
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(2,2)) {
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
		// anholonomic normalized
		//dt = (real)min(dt, (real)dx<?=side?>_at(i) / (fabs(lambda.max - lambda.min) + (real)1e-9));
		// holonomic
		dt = (real)min((real)dt, (real)(grid_dx<?=side?> / (fabs(lambda.max - lambda.min) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt;
}
