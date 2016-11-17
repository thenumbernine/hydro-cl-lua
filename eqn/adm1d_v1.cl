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
		.f = f,
		.alpha = alpha,
		.gamma_xx = gamma_xx,
	};
}

__kernel void calcEigenBasis(
	__global real* waveBuf,
	__global eigen_t* eigenBuf,
	const __global consLR_t* ULRBuf
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];

		cons_t UL = ULRBuf[side + dim * indexL].R;
		cons_t UR = ULRBuf[side + dim * indexR].L;
		eigen_t eig = calcEigenBasisSide(UL, UR);

		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		real lambda = calcMaxEigenvalue(eig.alpha, eig.gamma_xx); 
		wave[0] = -lambda;
		wave[1] = 0;
		wave[2] = lambda;
	
		eigenBuf[intindex].f = eig.f;	
	}
}

<? for side=0,2 do ?>

void eigen_leftTransform_<?=side?>(
	real* y,
	const __global eigen_t* eigen,
	const real* x
) {
	real sqrt_f = sqrt(eigen->f);
	y[0] = (x[2] / eigen->f - x[4] / sqrt_f) / 2.;
	y[1] = -2. * x[2] / eigen->f + x[3];
	y[2] = (x[2] / eigen->f + x[4] / sqrt_f) / 2.;
}

void eigen_rightTransform_<?=side?>(
	real* y,
	const __global eigen_t* eigen,
	const real* x
) {
	y[0] = 0;
	y[1] = 0;
	y[2] = (x[0] + x[2]) * eigen->f;
	y[3] = 2. * x[0] + x[1] + 2. * x[2];
	y[4] = sqrt(eigen->f) * (x[2] - x[0]);
}

<?	if solver.checkFluxError then ?>
void fluxTransform_<?=side?>(
	real* y,
	const __global eigen_t* eigen,
	const real* x
) {
	real alpha_sqrt_gamma_xx = eigen->alpha / sqrt(eigen->gamma_xx);
	y[0] = 0;
	y[1] = 0;
	y[2] = x[4] * eigen->f / alpha_sqrt_gamma_xx;
	y[3] = x[4] * 2. * alpha_sqrt_gamma_xx;
	y[4] = x[2] * alpha_sqrt_gamma_xx;
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
	real D_g = U->D_g;
	real KTilde_xx = U->KTilde_xx;
	
	real sqrt_gamma_xx = sqrt(gamma_xx);
	real K_xx = KTilde_xx / sqrt_gamma_xx;
	
	real f = calc_f(alpha);
	real dalpha_f = calc_dalpha_f(alpha);
	
	deriv->alpha -= alpha * alpha * f * K_xx / gamma_xx;
	deriv->gamma_xx -= 2. * alpha * K_xx;
	deriv->a_x -= alpha * K_xx * (f * (.5 * D_g - a_x) - a_x * alpha * dalpha_f);
	deriv->D_g -= 2. * alpha * K_xx * (.5 * D_g - a_x);
	deriv->KTilde_xx -= alpha * a_x / sqrt_gamma_xx * (.5 * D_g - a_x);
}
