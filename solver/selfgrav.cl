/*
used with Euler
but works with anything that has a cons_t::rho, m[], ETotal
*/

__kernel void initPotential(
	__global real* ePotBuf,
	const __global cons_t* UBuf)
{
	SETBOUNDS(0,0);
	ePotBuf[index] = -UBuf[index].rho;
}

__kernel void solvePoisson(
	__global real* potentialBuf,
	<?=args?>
) {
	SETBOUNDS(2,2);

	real dx = dx0_at(i);
<? if solver.dim > 1 then ?>
	real dy = dx1_at(i);
<? end
if solver.dim > 2 then ?>
	real dz = dx2_at(i);
<? end ?>

	real skewSum = 0;
	if (i.x < gridSize_x-3) skewSum += potentialBuf[index + stepsize.x] / (dx * dx);
	if (i.x > 2) skewSum += potentialBuf[index - stepsize.x] / (dx * dx);
<? if solver.dim > 1 then ?>
	if (i.y < gridSize_y-3) skewSum += potentialBuf[index + stepsize.y] / (dy * dy);
	if (i.y > 2) skewSum += potentialBuf[index - stepsize.y] / (dy * dy);
<? end
if solver.dim > 2 then ?>
	if (i.z < gridSize_z-3) skewSum += potentialBuf[index + stepsize.z] / (dz * dz);
	if (i.z > 2) skewSum += potentialBuf[index - stepsize.z] / (dz * dz);
<? end ?>

	const real diag = -2. * (1. / (dx * dx)
<? if solver.dim > 1 then ?>
		+ 1. / (dy * dy)
<? end
if solver.dim > 2 then ?>
		+ 1. / (dz * dz)
<? end ?>
	);

	real rho = 0;
	<?=calcRho?>

	potentialBuf[index] = (rho - skewSum) / diag;
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
	
		real gradient = (potentialBuf[indexR] - potentialBuf[indexL]) / (2. * dx<?=side?>_at(i));
		real gravity = -gradient;

		deriv->m.s[side] -= U->rho * gravity;
		deriv->ETotal -= U->rho * gravity * U->m.s[side];
	}<? end ?>
}
