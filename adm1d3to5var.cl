//TODO make parameters out of this somehow:
//real calc_f(alpha) { return 1.; }
//real calc_f(alpha) { return .49; }
//real calc_f(alpha) { return 1.69; }
//real calc_dalpha_f(alpha) { return 0.; }
constant real kappa = 1.;
real calc_f(alpha) { return 1. + kappa / (alpha * alpha); }
real calc_dalpha_f(alpha) { return -kappa / (alpha * alpha * alpha); }

__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real4 xs = CELL_X(i);
	real x = xs[0];
	__global cons_t* U = UBuf + index;
	U->alpha = init_calc_alpha(x);
	U->gamma_xx = init_calc_gamma_xx(x);
	U->a_x = init_calc_dx_alpha(x) / init_calc_alpha(x);
	U->d_xxx = .5 * init_calc_dx_gamma_xx(x);
	real K_xx = init_calc_K_xx(x);
	U->KTilde_xx = K_xx / sqrt(U->gamma_xx);
}

real calcDisplayVar_UBuf(int displayVar, const __global real* U_) {
	const __global cons_t* U = (const __global cons_t*)U_;
	switch (displayVar) {
	case display_U_alpha: return U->alpha;
	case display_U_gamma_xx: return U->gamma_xx;
	case display_U_a_x: return U->a_x;
	case display_U_d_xxx: return U->d_xxx;
	case display_U_KTilde_xx: return U->KTilde_xx;
	case display_U_K_xx: return U->KTilde_xx * sqrt(U->gamma_xx);
	case display_U_volume: return U->alpha * sqrt(U->gamma_xx);
	}
	return 0;
}

real calcMaxEigenvalue(real alpha, real gamma_xx) {
	real f = calc_f(alpha);
	real lambda = alpha * sqrt(f / gamma_xx);
	return lambda;
}

real2 calcCellMinMaxEigenvalues(const __global real* U_) {
	const __global cons_t* U = (const __global cons_t*)U_;
	real lambda = calcMaxEigenvalue(U->alpha, U->gamma_xx);
	return (real2)(-lambda, lambda);
}

typedef struct {
	real alpha, gamma_xx, f;
} Roe_t;

Roe_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	real gamma_xx = .5 * (UL.gamma_xx + UR.gamma_xx);
	real f = calc_f(alpha);
	return (Roe_t){.alpha=alpha, .gamma_xx=gamma_xx, .f=f};
}

__kernel void calcEigenBasis(
	__global real* waveBuf,
	__global eigen_t* eigenBuf,
	__global real* fluxMatrixBuf,
	const __global cons_t* UBuf
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
		wave[2] = 0;
		wave[3] = 0;
		wave[4] = lambda;
	
		eigenBuf[intindex].f = roe.f;	
	
		// flux matrix? fill with alpha, gamma_xx, f
	}
}

void eigen_leftTransform(
	real* y,
	const __global eigen_t* eigen,
	real* x
) {
	real f = eigen->f;
	real sqrt_f = sqrt(f);
	y[0] = x[2] / (2. * f) - x[4] / (2. * sqrt_f);
	y[1] = x[0];
	y[2] = x[1];
	y[3] = -2. * x[2] / f + x[3];
	y[4] = x[2] / (2. * f) + x[4] / (2. * sqrt_f);
}

void eigen_rightTransform(
	real* y,
	const __global eigen_t* eigen,
	real* x
) {
	real f = eigen->f;
	y[0] = x[1];
	y[1] = x[2];
	y[2] = (x[0] + x[4]) * f;
	y[3] = 2. * x[0] + x[3] + 2. * x[4];
	y[4] = sqrt(f) * (x[4] - x[0]);
}

real eigen_calcDisplayVar(
	int displayVar,
	const __global eigen_t* eigen
) {
	return eigen->f;
}
