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
	for (int side = 0; side < dim; ++side) {
		range_t lambda = calcCellMinMaxEigenvalues(U, side); 
		lambda.min = min((real)0., lambda.min);
		lambda.max = max((real)0., lambda.max);
		dt = min(dt, dxs[side] / (fabs(lambda.max - lambda.min) + (real)1e-9));
	}
	dtBuf[index] = dt; 
}
