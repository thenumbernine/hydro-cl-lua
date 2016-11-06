__kernel void calcDT(
	__global real* dtBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(2,2)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const __global cons_t* U = UBuf + index;

	real dt = INFINITY;
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		range_t lambda = calcCellMinMaxEigenvalues(U, <?=side?>); 
		lambda.min = min((real)0., lambda.min);
		lambda.max = max((real)0., lambda.max);
		dt = min(dt, dx<?=side?>_at(i) / (fabs(lambda.max - lambda.min) + (real)1e-9));
	}<? end ?>
	dtBuf[index] = dt; 
}
