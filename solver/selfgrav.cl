/*
used with Euler
but works with anything that has a cons_t::rho, m[], ETotal
*/

__kernel void solvePoisson(
	__global real* potentialBuf,
	const __global cons_t* UBuf)
{
	SETBOUNDS(2,2);

	real dx = dx_at0(i);
#if dim > 1
	real dy = dx_at1(i);
#endif
#if dim > 2
	real dz = dx_at2(i);
#endif

	real skewSum = 0;
	if (i.x < gridSize_x-3) skewSum += potentialBuf[index + stepsize.x] / (dx * dx);
	if (i.x > 2) skewSum += potentialBuf[index - stepsize.x] / (dx * dx);
#if dim > 1
	if (i.y < gridSize_y-3) skewSum += potentialBuf[index + stepsize.y] / (dy * dy);
	if (i.y > 2) skewSum += potentialBuf[index - stepsize.y] / (dy * dy);
#endif
#if dim > 2
	if (i.z < gridSize_z-3) skewSum += potentialBuf[index + stepsize.z] / (dz * dz);
	if (i.z > 2) skewSum += potentialBuf[index - stepsize.z] / (dz * dz);
#endif

	const real diag = -2. * (1. / (dx * dx)
#if dim > 1
		+ 1. / (dy * dy)
#endif
#if dim > 2
		+ 1. / (dz * dz)
#endif
	);

	const __global cons_t* U = UBuf + index;
	potentialBuf[index] = (4. * M_PI * gravitationalConstant * U->rho - skewSum) / diag;
}

__kernel void calcGravityDeriv(
	__global cons_t* derivBuffer,
	const __global cons_t* UBuf,
	const __global real* potentialBuf)
{
	SETBOUNDS(2,2);
	
	__global cons_t* deriv = derivBuffer + index;
	const __global cons_t* U = UBuf + index;

	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
	
		real gradient = (potentialBuf[indexR] - potentialBuf[indexL]) / (2. * dx_at<?=side?>(i));
		real gravity = -gradient;

		deriv->m.s[side] -= U->rho * gravity;
		deriv->ETotal -= U->rho * gravity * U->m.s[side];
	}<? end ?>
}
