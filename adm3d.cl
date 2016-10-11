real calcMaxEigenvalue(real alpha, real gammaUii) {
	return alpha * sqrt(calc_f(alpha) * gammaUii);
}

range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	int side
) {
	real gamma = symMatDet_prefix(U->gamma_);
	real lambda;
	switch (side) {
	case 0: 
		{
			real gammaUxx = (U->gamma_yy * U->gamma_zz - U->gamma_yz * U->gamma_yz) / gamma;
			lambda = calcMaxEigenvalue(U->alpha, gammaUxx);
		}
		break;
	case 1:
		{
			real gammaUyy = (U->gamma_xx * U->gamma_zz - U->gamma_xz * U->gamma_xz) / gamma;
			lambda = calcMaxEigenvalue(U->alpha, gammaUyy);
		}
		break;
	case 2:
		{
			real gammaUzz = (U->gamma_xx * U->gamma_yy - U->gamma_xy * U->gamma_xy) / gamma;
			lambda = calcMaxEigenvalue(U->alpha, gammaUzz);
		}
		break;
	}
	return (range_t){.min=-lambda, .max=lambda};
}

//alpha, gammaUxx, f
typedef eigen_t Roe_t;

Roe_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	Roe_t results;
	for (int i = 0; i < numStates; ++i) {
		((real*)&results.U)[i] = .5 * (((real*)&UL)[i] + ((real*)&UR)[i]);
	}

	results.gamma = symMatDet_prefix(results.U.gamma_);
	symMatInv_prefix((real*)&results.gammaUxx, results.gamma, results.gamma_);
	results.f = calc_f(results.alpha);

	return results;
}

__kernel void calcEigenBasis(
	__global real* waveBuf,
	__global eigen_t* eigenBuf,
	__global fluxXform_t* fluxXformBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];

		Roe_t roe = calcEigenBasisSide(UBuf[indexL], UBuf[indexR]);

		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;

		real lambdaLight;
		switch (side) {
		case 0:
			lambdaLight = roe.alpha * sqrt(roe.gammaUxx); 
			break;
		case 1:
			lambdaLight = roe.alpha * sqrt(roe.gammaUyy); 
			break;
		case 2:
			lambdaLight = roe.alpha * sqrt(roe.gammaUzz); 
			break;
		}
		real lambdaGauge = lambdaLight * sqrt(roe.f);
		eigenvalues[0] = -lambdaGauge;
		eigenvalues[1] = -lambdaLight;
		eigenvalues[2] = -lambdaLight;
		eigenvalues[3] = -lambdaLight;
		eigenvalues[4] = -lambdaLight;
		eigenvalues[5] = -lambdaLight;
		eigenvalues[6] = 0.f;
		eigenvalues[7] = 0.f;
		eigenvalues[8] = 0.f;
		eigenvalues[9] = 0.f;
		eigenvalues[10] = 0.f;
		eigenvalues[11] = 0.f;
		eigenvalues[12] = 0.f;
		eigenvalues[13] = 0.f;
		eigenvalues[14] = 0.f;
		eigenvalues[15] = 0.f;
		eigenvalues[16] = 0.f;
		eigenvalues[17] = 0.f;
		eigenvalues[18] = 0.f;
		eigenvalues[19] = 0.f;
		eigenvalues[20] = 0.f;
		eigenvalues[21] = 0.f;
		eigenvalues[22] = 0.f;
		eigenvalues[23] = 0.f;
		eigenvalues[24] = lambdaLight;
		eigenvalues[25] = lambdaLight;
		eigenvalues[26] = lambdaLight;
		eigenvalues[27] = lambdaLight;
		eigenvalues[28] = lambdaLight;
		eigenvalues[29] = lambdaGauge;

		eigenBuf[intindex] = roe;

		//only used for eigen basis reconstruction validity testing
		fluxXformBuf[intindex] = roe;
	}
}


