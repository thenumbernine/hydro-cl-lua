real calcMaxEigenvalue(real alpha, real gamma_xx) {
	real f = calc_f(alpha);
	real lambda = alpha * sqrt(f / gamma_xx);
	return lambda;
}

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const __global cons_t* U
) {
	real lambda = calcMaxEigenvalue(U->alpha, U->gamma_xx);
	return (range_t){.min=-lambda, .max=lambda};
}
<? end ?>

eigen_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	real gamma_xx = .5 * (UL.gamma_xx + UR.gamma_xx);
	real f = calc_f(alpha);
	return (eigen_t){
		.alpha = alpha,
		.sqrt_f_over_gamma_xx = sqrt(f / gamma_xx),
	};
}

__kernel void calcEigenBasis(
	__global real* waveBuf,
	__global eigen_t* eigenBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];

		eigen_t eig = calcEigenBasisSide(UBuf[indexL], UBuf[indexR]);

		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		real lambda = eig.alpha * eig.sqrt_f_over_gamma_xx;
		wave[0] = -lambda;
		wave[1] = 0;
		wave[2] = lambda;
	
		eigenBuf[intindex] = eig;
	}
}

<? for side=0,2 do ?>

void eigen_leftTransform_<?=side?>(
	real* y,
	const __global eigen_t* eigen,
	const real* x
) {
	real gamma_xx_over_f = 1. / (eigen->sqrt_f_over_gamma_xx * eigen->sqrt_f_over_gamma_xx);
	y[0] = .5 * (x[2] / eigen->sqrt_f_over_gamma_xx - x[4]);
	y[1] = x[3] - x[2] * gamma_xx_over_f;
	y[2] = .5 * (x[2] / eigen->sqrt_f_over_gamma_xx + x[4]);
}

void eigen_rightTransform_<?=side?>(
	real* y,
	const __global eigen_t* eigen,
	const real* x
) {
	y[0] = 0;
	y[1] = 0;
	y[2] = (x[0] + x[2]) * eigen->sqrt_f_over_gamma_xx;
	y[3] = (x[0] + x[2]) / eigen->sqrt_f_over_gamma_xx + x[1]; 
	y[4] = x[2] - x[0];
}

<?	if solver.checkFluxError then ?>
void fluxTransform_<?=side?>(
	real* y,
	const __global eigen_t* eigen,
	const real* x
) {
	real f_over_gamma_xx = eigen->sqrt_f_over_gamma_xx * eigen->sqrt_f_over_gamma_xx;
	
	y[0] = 0;
	y[1] = 0;
	y[2] = x[4] * eigen->alpha * f_over_gamma_xx;
	y[3] = x[4] * eigen->alpha; 
	y[4] = x[2] * eigen->alpha;
}
<?
	end
end
?>

kernel void addSource(
	__global cons_t* derivBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(2,2);
	__global cons_t* deriv = derivBuf + index;
	const __global cons_t* U = UBuf + index;
	
	real alpha = U->alpha;
	real gamma_xx = U->gamma_xx;
	real a_x = U->a_x;
	real d_xxx = U->d_xxx;
	real K_xx = U->K_xx;
	
	real f = calc_f(alpha);
	real dalpha_f = calc_dalpha_f(alpha);
	
	deriv->alpha -= alpha * alpha * f * K_xx / gamma_xx;
	deriv->gamma_xx -= 2. * alpha * K_xx;
	deriv->K_xx += alpha / gamma_xx * (a_x * d_xxx - K_xx * K_xx);
// terms that mysteriously disappear when you compare the linearized flux matrix terms moved to source, vs the source that Alcubierre uses in his 1997 paper
// adding these neglected terms back in make things blow up
#if 0 
	deriv->a_x += alpha * K_xx / gamma_xx * (f * (2. * d_xxx / gamma_xx - a_x) - a_x * alpha * dalpha_f);
	deriv->d_xxx -= alpha * a_x * K_xx; 
	deriv->K_xx -= alpha * a_x * a_x; 
#endif
}
