real calcMaxEigenvalue(real alpha, real gamma_xx) {
	real f = calc_f(alpha);
	real lambda = alpha * sqrt(f / gamma_xx);
	return lambda;
}

range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	int side
) {
	real lambda = calcMaxEigenvalue(U->alpha, U->gamma_xx);
	return (range_t){.min=-lambda, .max=lambda};
}

//alpha, gamma_xx, f
typedef fluxXform_t Roe_t;

Roe_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	real gamma_xx = .5 * (UL.gamma_xx + UR.gamma_xx);
	real f = calc_f(alpha);
	return (Roe_t){.alpha=alpha, .gamma_xx=gamma_xx, .f=f};
}

__kernel void calcEigenBasis(
	__global real* waveBuf,
	__global eigen_t* eigenBuf,
	const __global cons_t* UBuf
#if defined(checkFluxError)
	, __global fluxXform_t* fluxXformBuf
#endif
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];

		Roe_t roe = calcEigenBasisSide(UBuf[indexL], UBuf[indexR]);

		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		real lambda = calcMaxEigenvalue(roe.alpha, roe.gamma_xx); 
		wave[0] = -lambda;
		wave[1] = 0;
		wave[2] = lambda;
	
		eigenBuf[intindex].sqrt_g_xx_over_f = sqrt(roe.gamma_xx / roe.f);

#if defined(checkFluxError)
		//only used for eigen basis reconstruction validity testing
		fluxXformBuf[intindex] = roe;
#endif
	}
}

void fluxTransform(
	real* y,
	const __global fluxXform_t* flux,
	const real* x,
	int side
) {
	y[0] = 0;
	y[1] = 0;
	y[2] = x[4] * flux->alpha * flux->f / flux->gamma_xx;
	y[3] = x[4] * flux->alpha; 
	y[4] = x[2] * flux->alpha;
}

void eigen_leftTransform(
	real* y,
	const __global eigen_t* eigen,
	const real* x,
	int side
) {
	y[0] = .5 * (eigen->sqrt_g_xx_over_f  * x[2] - x[4]);
	y[1] = x[3] - eigen->sqrt_g_xx_over_f * eigen->sqrt_g_xx_over_f * x[2];
	y[2] = .5 * (eigen->sqrt_g_xx_over_f * x[2] + x[4]);
}

void eigen_rightTransform(
	real* y,
	const __global eigen_t* eigen,
	const real* x,
	int side
) {
	y[0] = 0;
	y[1] = 0;
	y[2] = (x[0] + x[2]) / eigen->sqrt_g_xx_over_f;
	y[3] = (x[0] + x[2]) * eigen->sqrt_g_xx_over_f + x[1]; 
	y[4] = x[2] - x[0];
}

kernel void addSourceTerm(
	__global cons_t* derivBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	const __global cons_t* U = UBuf + index;
	__global cons_t* deriv = derivBuf + index;
	
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
