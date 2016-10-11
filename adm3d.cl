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
typedef struct {
	eigen_t eig;
	real alpha;
} Roe_t;

Roe_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	Roe_t roe;
	roe.alpha = .5 * (UL.alpha + UR.alpha);
#define AVERAGE(var) real avg_##var = .5 * (UL.var + UR.var)
	AVERAGE(gamma_xx);
	AVERAGE(gamma_xy);
	AVERAGE(gamma_xz);
	AVERAGE(gamma_yy);
	AVERAGE(gamma_yz);
	AVERAGE(gamma_zz);
#undef AVERAGE

	real avg_gamma = symMatDet_prefix(avg_gamma_);
	symMatInv_prefix(roe.eig.gammaU, avg_gamma, avg_gamma_);
	roe.eig.f = calc_f(roe.alpha);

	return roe;
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
			lambdaLight = roe.alpha * sqrt(roe.eig.gammaUxx); 
			break;
		case 1:
			lambdaLight = roe.alpha * sqrt(roe.eig.gammaUyy); 
			break;
		case 2:
			lambdaLight = roe.alpha * sqrt(roe.eig.gammaUzz); 
			break;
		}
		real lambdaGauge = lambdaLight * sqrt(roe.eig.f);
		wave[0] = -lambdaGauge;
		wave[1] = -lambdaLight;
		wave[2] = -lambdaLight;
		wave[3] = -lambdaLight;
		wave[4] = -lambdaLight;
		wave[5] = -lambdaLight;
		wave[6] = 0.;
		wave[7] = 0.;
		wave[8] = 0.;
		wave[9] = 0.;
		wave[10] = 0.;
		wave[11] = 0.;
		wave[12] = 0.;
		wave[13] = 0.;
		wave[14] = 0.;
		wave[15] = 0.;
		wave[16] = 0.;
		wave[17] = 0.;
		wave[18] = 0.;
		wave[19] = 0.;
		wave[20] = 0.;
		wave[21] = 0.;
		wave[22] = 0.;
		wave[23] = 0.;
		wave[24] = lambdaLight;
		wave[25] = lambdaLight;
		wave[26] = lambdaLight;
		wave[27] = lambdaLight;
		wave[28] = lambdaLight;
		wave[29] = lambdaGauge;

		eigenBuf[intindex] = roe.eig;

		//only used for eigen basis reconstruction validity testing
		fluxXformBuf[intindex] = roe.eig;
	}
}

void fluxTransform(
	real* y,
	const __global fluxXform_t* flux,
	const real* x,
	int side
) {
	for (int i = 0; i < numStates; ++i) {
		*y = 0;
		++y;
	}
}

void eigen_leftTransform(
	real* results,
	const __global eigen_t* eigen,
	const real* input,
	int side
) {
	real gammaUxx = eigen->gammaUxx;
	real gammaUxy = eigen->gammaUxy;
	real gammaUxz = eigen->gammaUxz;
	real gammaUyy = eigen->gammaUyy;
	real gammaUyz = eigen->gammaUyz;
	real gammaUzz = eigen->gammaUzz;
	real f = eigen->f;
	real sqrt_f = sqrt(f);
	
	//input of left eigenvectors is the state
	//so skip the first 7 of input
	input += 7;
	
	// left eigenvectors in x:
	if (side == 0) {
		real sqrt_gammaUxx = sqrt(gammaUxx);
		real gammaUxx_toThe_3_2 = sqrt_gammaUxx * gammaUxx;
		
		results[0] = ((sqrt_f * gammaUxx_toThe_3_2 * input[21]) + (2. * sqrt_f * sqrt_gammaUxx * gammaUxy * input[22]) + (2. * sqrt_f * sqrt_gammaUxx * gammaUxz * input[23]) + (sqrt_f * sqrt_gammaUxx * gammaUyy * input[24]) + (2. * sqrt_f * sqrt_gammaUxx * gammaUyz * input[25]) + (((((((sqrt_f * sqrt_gammaUxx * gammaUzz * input[26]) - (gammaUxx * input[0])) - (2. * gammaUxx * input[27])) - (gammaUxy * input[1])) - (2. * gammaUxy * input[28])) - (gammaUxz * input[2])) - (2. * gammaUxz * input[29])));
		results[1] = ((-(input[1] + (2. * input[28]) + ((((2. * gammaUxx * input[4]) - (gammaUxx * input[9])) - (2. * sqrt_gammaUxx * input[22])) - (2. * gammaUxz * input[11])) + ((((2. * gammaUxz * input[16]) - (gammaUyy * input[12])) - (2. * gammaUyz * input[13])) - (gammaUzz * input[14])))) / 2.);
		results[2] = ((-(input[2] + (2. * input[29]) + (((2. * gammaUxx * input[5]) - (gammaUxx * input[15])) - (2. * sqrt_gammaUxx * input[23])) + (((((2. * gammaUxy * input[11]) - (2. * gammaUxy * input[16])) - (gammaUyy * input[18])) - (2. * gammaUyz * input[19])) - (gammaUzz * input[20])))) / 2.);
		results[3] = (-(((gammaUxx * input[6]) - (sqrt_gammaUxx * input[24])) + (gammaUxy * input[12]) + (gammaUxz * input[18])));
		results[4] = (-(((gammaUxx * input[7]) - (sqrt_gammaUxx * input[25])) + (gammaUxy * input[13]) + (gammaUxz * input[19])));
		results[5] = (-(((gammaUxx * input[8]) - (sqrt_gammaUxx * input[26])) + (gammaUxy * input[14]) + (gammaUxz * input[20])));
		results[6] = input[1];
		results[7] = input[2];
		results[8] = input[9];
		results[9] = input[10];
		results[10] = input[11];
		results[11] = input[12];
		results[12] = input[13];
		results[13] = input[14];
		results[14] = input[15];
		results[15] = input[16];
		results[16] = input[17];
		results[17] = input[18];
		results[18] = input[19];
		results[19] = input[20];
		results[20] = input[27];
		results[21] = input[28];
		results[22] = input[29];
		results[23] = ((((((input[0] - (f * gammaUxx * input[3])) - (2. * f * gammaUxy * input[4])) - (2. * f * gammaUxz * input[5])) - (f * gammaUyy * input[6])) - (2. * f * gammaUyz * input[7])) - (f * gammaUzz * input[8]));
		results[24] = ((input[1] + (2. * input[28]) + ((2. * gammaUxx * input[4]) - (gammaUxx * input[9])) + ((2. * sqrt_gammaUxx * input[22]) - (2. * gammaUxz * input[11])) + ((((2. * gammaUxz * input[16]) - (gammaUyy * input[12])) - (2. * gammaUyz * input[13])) - (gammaUzz * input[14]))) / 2.);
		results[25] = ((input[2] + (2. * input[29]) + ((2. * gammaUxx * input[5]) - (gammaUxx * input[15])) + (2. * sqrt_gammaUxx * input[23]) + (((((2. * gammaUxy * input[11]) - (2. * gammaUxy * input[16])) - (gammaUyy * input[18])) - (2. * gammaUyz * input[19])) - (gammaUzz * input[20]))) / 2.);
		results[26] = ((gammaUxx * input[6]) + (sqrt_gammaUxx * input[24]) + (gammaUxy * input[12]) + (gammaUxz * input[18]));
		results[27] = ((gammaUxx * input[7]) + (sqrt_gammaUxx * input[25]) + (gammaUxy * input[13]) + (gammaUxz * input[19]));
		results[28] = ((gammaUxx * input[8]) + (sqrt_gammaUxx * input[26]) + (gammaUxy * input[14]) + (gammaUxz * input[20]));
		results[29] = ((sqrt_f * gammaUxx_toThe_3_2 * input[21]) + (2. * sqrt_f * sqrt_gammaUxx * gammaUxy * input[22]) + (2. * sqrt_f * sqrt_gammaUxx * gammaUxz * input[23]) + (sqrt_f * sqrt_gammaUxx * gammaUyy * input[24]) + (2. * sqrt_f * sqrt_gammaUxx * gammaUyz * input[25]) + (sqrt_f * sqrt_gammaUxx * gammaUzz * input[26]) + (gammaUxx * input[0]) + (2. * gammaUxx * input[27]) + (gammaUxy * input[1]) + (2. * gammaUxy * input[28]) + (gammaUxz * input[2]) + (2. * gammaUxz * input[29]));
	// left eigenvectors in y:
	} else if (side == 1) {
		real sqrt_gammaUyy = sqrt(gammaUyy);
		real gammaUyy_toThe_3_2 = sqrt_gammaUyy * gammaUyy;
		
		results[0] = ((sqrt_f * gammaUyy_toThe_3_2 * input[24]) + (sqrt_f * sqrt_gammaUyy * gammaUxx * input[21]) + (2. * sqrt_f * sqrt_gammaUyy * gammaUxy * input[22]) + (2. * sqrt_f * sqrt_gammaUyy * gammaUxz * input[23]) + (2. * sqrt_f * sqrt_gammaUyy * gammaUyz * input[25]) + (((((((sqrt_f * sqrt_gammaUyy * gammaUzz * input[26]) - (gammaUxy * input[0])) - (2. * gammaUxy * input[27])) - (gammaUyy * input[1])) - (2. * gammaUyy * input[28])) - (gammaUyz * input[2])) - (2. * gammaUyz * input[29])));
		results[1] = (-((gammaUxy * input[3]) + ((gammaUyy * input[9]) - (sqrt_gammaUyy * input[21])) + (gammaUyz * input[15])));
		results[2] = ((-(input[0] + ((((2. * input[27]) - (gammaUxx * input[3])) - (2. * gammaUxz * input[5])) - (gammaUyy * input[6])) + (((2. * gammaUyy * input[10]) - (2. * sqrt_gammaUyy * input[22])) - (2. * gammaUyz * input[7])) + ((2. * gammaUyz * input[16]) - (gammaUzz * input[8])))) / 2.);
		results[3] = (-((gammaUxy * input[5]) + ((gammaUyy * input[11]) - (sqrt_gammaUyy * input[23])) + (gammaUyz * input[17])));
		results[4] = ((-(input[2] + ((2. * input[29]) - (gammaUxx * input[15])) + (((2. * gammaUxy * input[7]) - (2. * gammaUxy * input[16])) - (2. * gammaUxz * input[17])) + ((((2. * gammaUyy * input[13]) - (gammaUyy * input[18])) - (2. * sqrt_gammaUyy * input[25])) - (gammaUzz * input[20])))) / 2.);
		results[5] = (-((gammaUxy * input[8]) + ((gammaUyy * input[14]) - (sqrt_gammaUyy * input[26])) + (gammaUyz * input[20])));
		results[6] = input[0];
		results[7] = input[2];
		results[8] = input[3];
		results[9] = input[4];
		results[10] = input[5];
		results[11] = input[6];
		results[12] = input[7];
		results[13] = input[8];
		results[14] = input[15];
		results[15] = input[16];
		results[16] = input[17];
		results[17] = input[18];
		results[18] = input[19];
		results[19] = input[20];
		results[20] = input[27];
		results[21] = input[28];
		results[22] = input[29];
		results[23] = ((((((input[1] - (f * gammaUxx * input[9])) - (2. * f * gammaUxy * input[10])) - (2. * f * gammaUxz * input[11])) - (f * gammaUyy * input[12])) - (2. * f * gammaUyz * input[13])) - (f * gammaUzz * input[14]));
		results[24] = ((gammaUxy * input[3]) + (gammaUyy * input[9]) + (sqrt_gammaUyy * input[21]) + (gammaUyz * input[15]));
		results[25] = ((input[0] + ((((2. * input[27]) - (gammaUxx * input[3])) - (2. * gammaUxz * input[5])) - (gammaUyy * input[6])) + (2. * gammaUyy * input[10]) + ((2. * sqrt_gammaUyy * input[22]) - (2. * gammaUyz * input[7])) + ((2. * gammaUyz * input[16]) - (gammaUzz * input[8]))) / 2.);
		results[26] = ((gammaUxy * input[5]) + (gammaUyy * input[11]) + (sqrt_gammaUyy * input[23]) + (gammaUyz * input[17]));
		results[27] = ((input[2] + ((2. * input[29]) - (gammaUxx * input[15])) + (((2. * gammaUxy * input[7]) - (2. * gammaUxy * input[16])) - (2. * gammaUxz * input[17])) + ((2. * gammaUyy * input[13]) - (gammaUyy * input[18])) + ((2. * sqrt_gammaUyy * input[25]) - (gammaUzz * input[20]))) / 2.);
		results[28] = ((gammaUxy * input[8]) + (gammaUyy * input[14]) + (sqrt_gammaUyy * input[26]) + (gammaUyz * input[20]));
		results[29] = ((sqrt_f * gammaUyy_toThe_3_2 * input[24]) + (sqrt_f * sqrt_gammaUyy * gammaUxx * input[21]) + (2. * sqrt_f * sqrt_gammaUyy * gammaUxy * input[22]) + (2. * sqrt_f * sqrt_gammaUyy * gammaUxz * input[23]) + (2. * sqrt_f * sqrt_gammaUyy * gammaUyz * input[25]) + (sqrt_f * sqrt_gammaUyy * gammaUzz * input[26]) + (gammaUxy * input[0]) + (2. * gammaUxy * input[27]) + (gammaUyy * input[1]) + (2. * gammaUyy * input[28]) + (gammaUyz * input[2]) + (2. * gammaUyz * input[29]));
	// left eigenvectors in z:
	} else if (side == 2) {
		real sqrt_gammaUzz = sqrt(gammaUzz);
		real gammaUzz_toThe_3_2 = sqrt_gammaUzz * gammaUzz;
		
		results[0] = ((sqrt_f * gammaUzz_toThe_3_2 * input[26]) + (sqrt_f * sqrt_gammaUzz * gammaUxx * input[21]) + (2. * sqrt_f * sqrt_gammaUzz * gammaUxy * input[22]) + (2. * sqrt_f * sqrt_gammaUzz * gammaUxz * input[23]) + (sqrt_f * sqrt_gammaUzz * gammaUyy * input[24]) + (((((((2. * sqrt_f * sqrt_gammaUzz * gammaUyz * input[25]) - (gammaUxz * input[0])) - (2. * gammaUxz * input[27])) - (gammaUyz * input[1])) - (2. * gammaUyz * input[28])) - (gammaUzz * input[2])) - (2. * gammaUzz * input[29])));
		results[1] = (-((gammaUxz * input[3]) + (gammaUyz * input[9]) + ((gammaUzz * input[15]) - (sqrt_gammaUzz * input[21]))));
		results[2] = (-((gammaUxz * input[4]) + (gammaUyz * input[10]) + ((gammaUzz * input[16]) - (sqrt_gammaUzz * input[22]))));
		results[3] = ((-(input[0] + (((((2. * input[27]) - (gammaUxx * input[3])) - (2. * gammaUxy * input[4])) - (gammaUyy * input[6])) - (2. * gammaUyz * input[7])) + ((2. * gammaUyz * input[11]) - (gammaUzz * input[8])) + ((2. * gammaUzz * input[17]) - (2. * sqrt_gammaUzz * input[23])))) / 2.);
		results[4] = (-((gammaUxz * input[6]) + (gammaUyz * input[12]) + ((gammaUzz * input[18]) - (sqrt_gammaUzz * input[24]))));
		results[5] = ((-(input[1] + (((2. * input[28]) - (gammaUxx * input[9])) - (2. * gammaUxy * input[10])) + ((((2. * gammaUxz * input[7]) - (2. * gammaUxz * input[11])) - (gammaUyy * input[12])) - (gammaUzz * input[14])) + ((2. * gammaUzz * input[19]) - (2. * sqrt_gammaUzz * input[25])))) / 2.);
		results[6] = input[0];
		results[7] = input[1];
		results[8] = input[3];
		results[9] = input[4];
		results[10] = input[5];
		results[11] = input[6];
		results[12] = input[7];
		results[13] = input[8];
		results[14] = input[9];
		results[15] = input[10];
		results[16] = input[11];
		results[17] = input[12];
		results[18] = input[13];
		results[19] = input[14];
		results[20] = input[27];
		results[21] = input[28];
		results[22] = input[29];
		results[23] = ((((((input[2] - (f * gammaUxx * input[15])) - (2. * f * gammaUxy * input[16])) - (2. * f * gammaUxz * input[17])) - (f * gammaUyy * input[18])) - (2. * f * gammaUyz * input[19])) - (f * gammaUzz * input[20]));
		results[24] = ((gammaUxz * input[3]) + (gammaUyz * input[9]) + (gammaUzz * input[15]) + (sqrt_gammaUzz * input[21]));
		results[25] = ((gammaUxz * input[4]) + (gammaUyz * input[10]) + (gammaUzz * input[16]) + (sqrt_gammaUzz * input[22]));
		results[26] = ((input[0] + (((((2. * input[27]) - (gammaUxx * input[3])) - (2. * gammaUxy * input[4])) - (gammaUyy * input[6])) - (2. * gammaUyz * input[7])) + ((2. * gammaUyz * input[11]) - (gammaUzz * input[8])) + (2. * gammaUzz * input[17]) + (2. * sqrt_gammaUzz * input[23])) / 2.);
		results[27] = ((gammaUxz * input[6]) + (gammaUyz * input[12]) + (gammaUzz * input[18]) + (sqrt_gammaUzz * input[24]));
		results[28] = ((input[1] + (((2. * input[28]) - (gammaUxx * input[9])) - (2. * gammaUxy * input[10])) + ((((2. * gammaUxz * input[7]) - (2. * gammaUxz * input[11])) - (gammaUyy * input[12])) - (gammaUzz * input[14])) + (2. * gammaUzz * input[19]) + (2. * sqrt_gammaUzz * input[25])) / 2.);
		results[29] = ((sqrt_f * gammaUzz_toThe_3_2 * input[26]) + (sqrt_f * sqrt_gammaUzz * gammaUxx * input[21]) + (2. * sqrt_f * sqrt_gammaUzz * gammaUxy * input[22]) + (2. * sqrt_f * sqrt_gammaUzz * gammaUxz * input[23]) + (sqrt_f * sqrt_gammaUzz * gammaUyy * input[24]) + (2. * sqrt_f * sqrt_gammaUzz * gammaUyz * input[25]) + (gammaUxz * input[0]) + (2. * gammaUxz * input[27]) + (gammaUyz * input[1]) + (2. * gammaUyz * input[28]) + (gammaUzz * input[2]) + (2. * gammaUzz * input[29]));
	}
}

void eigen_rightTransform(
	real* results,
	const __global eigen_t* eigen,
	const real* input,
	int side)
{
	real gammaUxx = eigen->gammaUxx;
	real gammaUxy = eigen->gammaUxy;
	real gammaUxz = eigen->gammaUxz;
	real gammaUyy = eigen->gammaUyy;
	real gammaUyz = eigen->gammaUyz;
	real gammaUzz = eigen->gammaUzz;
	real f = eigen->f;
	real sqrt_f = sqrt(f);

	//write zeros to the alpha and gammaLL terms
	for (int i = 0; i < 7; ++i) {
		*results = 0;
		++results;
	}

	// right eigenvectors in x:
	if (side == 0) {
		real sqrt_gammaUxx = sqrt(gammaUxx);
		real gammaUxx_toThe_3_2 = sqrt_gammaUxx * gammaUxx;
		real gammaUxxSq = gammaUxx * gammaUxx;
		
		results[0] = ((input[0] + ((4. * input[20] * gammaUxx) - input[29]) + (2. * gammaUxy * input[6]) + (4. * gammaUxy * input[21]) + (2. * gammaUxz * input[7]) + (4. * gammaUxz * input[22])) / (-(2. * gammaUxx)));
		results[1] = input[6];
		results[2] = input[7];
		results[3] = ((-(input[0] + (4. * input[20] * gammaUxx) + (((2. * input[23] * gammaUxx) - input[29]) - (2. * gammaUxy * input[1] * f)) + (2. * gammaUxy * input[6]) + (2. * gammaUxy * input[8] * gammaUxx * f) + (4. * gammaUxy * input[21]) + ((((2. * gammaUxy * input[24] * f) - (2. * gammaUxy * f * input[6])) - (4. * gammaUxy * f * input[21])) - (2. * gammaUxz * input[2] * f)) + (2. * gammaUxz * input[7]) + (2. * gammaUxz * input[14] * gammaUxx * f) + (4. * gammaUxz * input[22]) + ((((2. * gammaUxz * input[25] * f) - (2. * gammaUxz * f * input[7])) - (4. * gammaUxz * f * input[22])) - (gammaUyy * input[3] * f)) + ((gammaUyy * input[26] * f) - (2. * gammaUyz * input[4] * f)) + ((2. * gammaUyz * input[27] * f) - (gammaUzz * input[5] * f)) + (gammaUzz * input[28] * f))) / (2. * gammaUxxSq * f));
		results[4] = ((input[1] + (input[6] - (input[8] * gammaUxx)) + (((2. * input[21]) - input[24]) - (2. * gammaUxz * input[10])) + ((((2. * gammaUxz * input[15]) - (gammaUyy * input[11])) - (2. * gammaUyz * input[12])) - (gammaUzz * input[13]))) / (-(2. * gammaUxx)));
		results[5] = ((input[2] + (input[7] - (input[14] * gammaUxx)) + ((2. * input[22]) - input[25]) + (((((2. * gammaUxy * input[10]) - (2. * gammaUxy * input[15])) - (gammaUyy * input[17])) - (2. * gammaUyz * input[18])) - (gammaUzz * input[19]))) / (-(2. * gammaUxx)));
		results[6] = (((input[3] - input[26]) + (2. * gammaUxy * input[11]) + (2. * gammaUxz * input[17])) / (-(2. * gammaUxx)));
		results[7] = (((input[4] - input[27]) + (2. * gammaUxy * input[12]) + (2. * gammaUxz * input[18])) / (-(2. * gammaUxx)));
		results[8] = (((input[5] - input[28]) + (2. * gammaUxy * input[13]) + (2. * gammaUxz * input[19])) / (-(2. * gammaUxx)));
		results[9] = input[8];
		results[10] = input[9];
		results[11] = input[10];
		results[12] = input[11];
		results[13] = input[12];
		results[14] = input[13];
		results[15] = input[14];
		results[16] = input[15];
		results[17] = input[16];
		results[18] = input[17];
		results[19] = input[18];
		results[20] = input[19];
		results[21] = ((input[0] + ((((((((((input[29] - (2. * gammaUxy * input[1] * sqrt_f)) - (2. * gammaUxy * input[24] * sqrt_f)) - (2. * gammaUxz * input[2] * sqrt_f)) - (2. * gammaUxz * input[25] * sqrt_f)) - (gammaUyy * input[3] * sqrt_f)) - (gammaUyy * input[26] * sqrt_f)) - (2. * gammaUyz * input[4] * sqrt_f)) - (2. * gammaUyz * input[27] * sqrt_f)) - (gammaUzz * input[5] * sqrt_f)) - (gammaUzz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUxx_toThe_3_2));
		results[22] = ((input[1] + input[24]) / (2. * sqrt_gammaUxx));
		results[23] = ((input[2] + input[25]) / (2. * sqrt_gammaUxx));
		results[24] = ((input[3] + input[26]) / (2. * sqrt_gammaUxx));
		results[25] = ((input[4] + input[27]) / (2. * sqrt_gammaUxx));
		results[26] = ((input[5] + input[28]) / (2. * sqrt_gammaUxx));
		results[27] = input[20];
		results[28] = input[21];
		results[29] = input[22];
	// right eigenvectors in y:
	} else if (side == 1) {
		real sqrt_gammaUyy = sqrt(gammaUyy);
		real gammaUyy_toThe_3_2 = sqrt_gammaUyy * gammaUyy;
		real gammaUyySq = gammaUyy * gammaUyy;
		
		results[0] = input[6];
		results[1] = ((-(input[0] + ((4. * input[21] * gammaUyy) - input[29]) + (2. * gammaUxy * input[6]) + (4. * gammaUxy * input[20]) + (2. * gammaUyz * input[7]) + (4. * gammaUyz * input[22]))) / (2. * gammaUyy));
		results[2] = input[7];
		results[3] = input[8];
		results[4] = input[9];
		results[5] = input[10];
		results[6] = input[11];
		results[7] = input[12];
		results[8] = input[13];
		results[9] = (((input[1] - input[24]) + (2. * gammaUxy * input[8]) + (2. * gammaUyz * input[14])) / (-(2. * gammaUyy)));
		results[10] = ((-(input[2] + (input[6] - (input[11] * gammaUyy)) + (((((2. * input[20]) - input[25]) - (gammaUxx * input[8])) - (2. * gammaUxz * input[10])) - (2. * gammaUyz * input[12])) + ((2. * gammaUyz * input[15]) - (gammaUzz * input[13])))) / (2. * gammaUyy));
		results[11] = (((input[3] - input[26]) + (2. * gammaUxy * input[10]) + (2. * gammaUyz * input[16])) / (-(2. * gammaUyy)));
		results[12] = ((-(input[0] + (4. * input[21] * gammaUyy) + (((2. * input[23] * gammaUyy) - input[29]) - (gammaUxx * input[1] * f)) + ((gammaUxx * input[24] * f) - (2. * gammaUxy * input[2] * f)) + (2. * gammaUxy * input[6]) + (2. * gammaUxy * input[11] * gammaUyy * f) + (4. * gammaUxy * input[20]) + ((((2. * gammaUxy * input[25] * f) - (2. * gammaUxy * f * input[6])) - (4. * gammaUxy * f * input[20])) - (2. * gammaUxz * input[3] * f)) + ((2. * gammaUxz * input[26] * f) - (2. * gammaUyz * input[4] * f)) + (2. * gammaUyz * input[7]) + (2. * gammaUyz * input[17] * gammaUyy * f) + (4. * gammaUyz * input[22]) + ((((2. * gammaUyz * input[27] * f) - (2. * gammaUyz * f * input[7])) - (4. * gammaUyz * f * input[22])) - (gammaUzz * input[5] * f)) + (gammaUzz * input[28] * f))) / (2. * gammaUyySq * f));
		results[13] = ((input[4] + (input[7] - (input[17] * gammaUyy)) + (((2. * input[22]) - input[27]) - (gammaUxx * input[14])) + ((((2. * gammaUxy * input[12]) - (2. * gammaUxy * input[15])) - (2. * gammaUxz * input[16])) - (gammaUzz * input[19]))) / (-(2. * gammaUyy)));
		results[14] = (((input[5] - input[28]) + (2. * gammaUxy * input[13]) + (2. * gammaUyz * input[19])) / (-(2. * gammaUyy)));
		results[15] = input[14];
		results[16] = input[15];
		results[17] = input[16];
		results[18] = input[17];
		results[19] = input[18];
		results[20] = input[19];
		results[21] = ((input[1] + input[24]) / (2. * sqrt_gammaUyy));
		results[22] = ((input[2] + input[25]) / (2. * sqrt_gammaUyy));
		results[23] = ((input[3] + input[26]) / (2. * sqrt_gammaUyy));
		results[24] = ((input[0] + ((((((((((input[29] - (gammaUxx * input[1] * sqrt_f)) - (gammaUxx * input[24] * sqrt_f)) - (2. * gammaUxy * input[2] * sqrt_f)) - (2. * gammaUxy * input[25] * sqrt_f)) - (2. * gammaUxz * input[3] * sqrt_f)) - (2. * gammaUxz * input[26] * sqrt_f)) - (2. * gammaUyz * input[4] * sqrt_f)) - (2. * gammaUyz * input[27] * sqrt_f)) - (gammaUzz * input[5] * sqrt_f)) - (gammaUzz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUyy_toThe_3_2));
		results[25] = ((input[4] + input[27]) / (2. * sqrt_gammaUyy));
		results[26] = ((input[5] + input[28]) / (2. * sqrt_gammaUyy));
		results[27] = input[20];
		results[28] = input[21];
		results[29] = input[22];
	// right eigenvectors in z:
	} else if (side == 2) {
		real sqrt_gammaUzz = sqrt(gammaUzz);
		real gammaUzz_toThe_3_2 = sqrt_gammaUzz * gammaUzz;
		real gammaUzzSq = gammaUzz * gammaUzz;
		
		results[0] = input[6];
		results[1] = input[7];
		results[2] = ((-(input[0] + ((4. * input[22] * gammaUzz) - input[29]) + (2. * gammaUxz * input[6]) + (4. * gammaUxz * input[20]) + (2. * gammaUyz * input[7]) + (4. * gammaUyz * input[21]))) / (2. * gammaUzz));
		results[3] = input[8];
		results[4] = input[9];
		results[5] = input[10];
		results[6] = input[11];
		results[7] = input[12];
		results[8] = input[13];
		results[9] = input[14];
		results[10] = input[15];
		results[11] = input[16];
		results[12] = input[17];
		results[13] = input[18];
		results[14] = input[19];
		results[15] = (((input[1] - input[24]) + (2. * gammaUxz * input[8]) + (2. * gammaUyz * input[14])) / (-(2. * gammaUzz)));
		results[16] = (((input[2] - input[25]) + (2. * gammaUxz * input[9]) + (2. * gammaUyz * input[15])) / (-(2. * gammaUzz)));
		results[17] = ((-(input[3] + (input[6] - (input[13] * gammaUzz)) + ((((((2. * input[20]) - input[26]) - (gammaUxx * input[8])) - (2. * gammaUxy * input[9])) - (gammaUyy * input[11])) - (2. * gammaUyz * input[12])) + (2. * gammaUyz * input[16]))) / (2. * gammaUzz));
		results[18] = (((input[4] - input[27]) + (2. * gammaUxz * input[11]) + (2. * gammaUyz * input[17])) / (-(2. * gammaUzz)));
		results[19] = ((-(input[5] + (input[7] - (input[19] * gammaUzz)) + ((((2. * input[21]) - input[28]) - (gammaUxx * input[14])) - (2. * gammaUxy * input[15])) + (((2. * gammaUxz * input[12]) - (2. * gammaUxz * input[16])) - (gammaUyy * input[17])))) / (2. * gammaUzz));
		results[20] = ((-(input[0] + (4. * input[22] * gammaUzz) + (((2. * input[23] * gammaUzz) - input[29]) - (gammaUxx * input[1] * f)) + ((gammaUxx * input[24] * f) - (2. * gammaUxy * input[2] * f)) + ((2. * gammaUxy * input[25] * f) - (2. * gammaUxz * input[3] * f)) + (2. * gammaUxz * input[6]) + (2. * gammaUxz * input[13] * gammaUzz * f) + (4. * gammaUxz * input[20]) + ((((2. * gammaUxz * input[26] * f) - (2. * gammaUxz * f * input[6])) - (4. * gammaUxz * f * input[20])) - (gammaUyy * input[4] * f)) + ((gammaUyy * input[27] * f) - (2. * gammaUyz * input[5] * f)) + (2. * gammaUyz * input[7]) + (2. * gammaUyz * input[19] * gammaUzz * f) + (4. * gammaUyz * input[21]) + (((2. * gammaUyz * input[28] * f) - (2. * gammaUyz * f * input[7])) - (4. * gammaUyz * f * input[21])))) / (2. * gammaUzzSq * f));
		results[21] = ((input[1] + input[24]) / (2. * sqrt_gammaUzz));
		results[22] = ((input[2] + input[25]) / (2. * sqrt_gammaUzz));
		results[23] = ((input[3] + input[26]) / (2. * sqrt_gammaUzz));
		results[24] = ((input[4] + input[27]) / (2. * sqrt_gammaUzz));
		results[25] = ((input[5] + input[28]) / (2. * sqrt_gammaUzz));
		results[26] = ((input[0] + ((((((((((input[29] - (gammaUxx * input[1] * sqrt_f)) - (gammaUxx * input[24] * sqrt_f)) - (2. * gammaUxy * input[2] * sqrt_f)) - (2. * gammaUxy * input[25] * sqrt_f)) - (2. * gammaUxz * input[3] * sqrt_f)) - (2. * gammaUxz * input[26] * sqrt_f)) - (gammaUyy * input[4] * sqrt_f)) - (gammaUyy * input[27] * sqrt_f)) - (2. * gammaUyz * input[5] * sqrt_f)) - (2. * gammaUyz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUzz_toThe_3_2));
		results[27] = input[20];
		results[28] = input[21];
		results[29] = input[22];
	}
}

__kernel void addSourceTerm(
	__global cons_t* derivBuf,
	const __global cons_t* UBuf)
{
	SETBOUNDS(2,2);
	
	__global cons_t* deriv = derivBuf + index;
	const __global cons_t* U = UBuf + index;

	real gamma = symMatDet_prefix(U->gamma_);
	real gammaInv[6];
	symMatInv_prefix(gammaInv, gamma, U->gamma_);
	real gammaUxx = gammaInv[0], gammaUxy = gammaInv[1], gammaUxz = gammaInv[2], gammaUyy = gammaInv[3], gammaUyz = gammaInv[4], gammaUzz = gammaInv[5];
	real f = calc_f(U->alpha);	//could be based on alpha...

real density = 0;//state[STATE_DENSITY];
real pressure = 0;//state[STATE_PRESSURE];
real4 vel3 = (real4)(0,0,0,0);//(real4)(state[STATE_VELOCITY_X], state[STATE_VELOCITY_Y], state[STATE_VELOCITY_Z], 0.);
real vel3Sq = dot(vel3, vel3);
real LorentzFactor = 1. / sqrt(1. - vel3Sq);
real4 vel4_ = vel3 * LorentzFactor;
vel4_.w = LorentzFactor;

real4 beta_ = (real4)(0., 0., 0., 0.);

// source terms
real KUL[3][3] = {
{gammaUxx * U->K_xx + gammaUxy * U->K_xy + gammaUxz * U->K_xz,
gammaUxx * U->K_xy + gammaUxy * U->K_yy + gammaUxz * U->K_yz,
gammaUxx * U->K_xz + gammaUxy * U->K_yz + gammaUxz * U->K_zz,
},{gammaUxy * U->K_xx + gammaUyy * U->K_xy + gammaUyz * U->K_xz,
gammaUxy * U->K_xy + gammaUyy * U->K_yy + gammaUyz * U->K_yz,
gammaUxy * U->K_xz + gammaUyy * U->K_yz + gammaUyz * U->K_zz,
},{gammaUxz * U->K_xx + gammaUyz * U->K_xy + gammaUzz * U->K_xz,
gammaUxz * U->K_xy + gammaUyz * U->K_yy + gammaUzz * U->K_yz,
gammaUxz * U->K_xz + gammaUyz * U->K_yz + gammaUzz * U->K_zz,
},};
real trK = KUL[0][0] + KUL[1][1] + KUL[2][2];
real KSqSymLL[6] = {
U->K_xx * KUL[0][0] + U->K_xy * KUL[1][0] + U->K_xz * KUL[2][0],
U->K_xx * KUL[0][1] + U->K_xy * KUL[1][1] + U->K_xz * KUL[2][1],
U->K_xx * KUL[0][2] + U->K_xy * KUL[1][2] + U->K_xz * KUL[2][2],
U->K_xy * KUL[0][1] + U->K_yy * KUL[1][1] + U->K_yz * KUL[2][1],
U->K_xy * KUL[0][2] + U->K_yy * KUL[1][2] + U->K_yz * KUL[2][2],
U->K_xz * KUL[0][2] + U->K_yz * KUL[1][2] + U->K_zz * KUL[2][2],
};
real DLUL[3][3][3] = {
{{U->d_xxx * gammaUxx + U->d_xxy * gammaUxy + U->d_xxz * gammaUxz,
U->d_xxy * gammaUxx + U->d_xyy * gammaUxy + U->d_xyz * gammaUxz,
U->d_xxz * gammaUxx + U->d_xyz * gammaUxy + U->d_xzz * gammaUxz,
},{U->d_xxx * gammaUxy + U->d_xxy * gammaUyy + U->d_xxz * gammaUyz,
U->d_xxy * gammaUxy + U->d_xyy * gammaUyy + U->d_xyz * gammaUyz,
U->d_xxz * gammaUxy + U->d_xyz * gammaUyy + U->d_xzz * gammaUyz,
},{U->d_xxx * gammaUxz + U->d_xxy * gammaUyz + U->d_xxz * gammaUzz,
U->d_xxy * gammaUxz + U->d_xyy * gammaUyz + U->d_xyz * gammaUzz,
U->d_xxz * gammaUxz + U->d_xyz * gammaUyz + U->d_xzz * gammaUzz,
},},{{U->d_yxx * gammaUxx + U->d_yxy * gammaUxy + U->d_yxz * gammaUxz,
U->d_yxy * gammaUxx + U->d_yyy * gammaUxy + U->d_yyz * gammaUxz,
U->d_yxz * gammaUxx + U->d_yyz * gammaUxy + U->d_yzz * gammaUxz,
},{U->d_yxx * gammaUxy + U->d_yxy * gammaUyy + U->d_yxz * gammaUyz,
U->d_yxy * gammaUxy + U->d_yyy * gammaUyy + U->d_yyz * gammaUyz,
U->d_yxz * gammaUxy + U->d_yyz * gammaUyy + U->d_yzz * gammaUyz,
},{U->d_yxx * gammaUxz + U->d_yxy * gammaUyz + U->d_yxz * gammaUzz,
U->d_yxy * gammaUxz + U->d_yyy * gammaUyz + U->d_yyz * gammaUzz,
U->d_yxz * gammaUxz + U->d_yyz * gammaUyz + U->d_yzz * gammaUzz,
},},{{U->d_zxx * gammaUxx + U->d_zxy * gammaUxy + U->d_zxz * gammaUxz,
U->d_zxy * gammaUxx + U->d_zyy * gammaUxy + U->d_zyz * gammaUxz,
U->d_zxz * gammaUxx + U->d_zyz * gammaUxy + U->d_zzz * gammaUxz,
},{U->d_zxx * gammaUxy + U->d_zxy * gammaUyy + U->d_zxz * gammaUyz,
U->d_zxy * gammaUxy + U->d_zyy * gammaUyy + U->d_zyz * gammaUyz,
U->d_zxz * gammaUxy + U->d_zyz * gammaUyy + U->d_zzz * gammaUyz,
},{U->d_zxx * gammaUxz + U->d_zxy * gammaUyz + U->d_zxz * gammaUzz,
U->d_zxy * gammaUxz + U->d_zyy * gammaUyz + U->d_zyz * gammaUzz,
U->d_zxz * gammaUxz + U->d_zyz * gammaUyz + U->d_zzz * gammaUzz,
},},};
real D1L[3] = {
DLUL[0][0][0] + DLUL[0][1][1] + DLUL[0][2][2],
DLUL[1][0][0] + DLUL[1][1][1] + DLUL[1][2][2],
DLUL[2][0][0] + DLUL[2][1][1] + DLUL[2][2][2],
};
real D3L[3] = {
DLUL[0][0][0] + DLUL[1][1][0] + DLUL[2][2][0],
DLUL[0][0][1] + DLUL[1][1][1] + DLUL[2][2][1],
DLUL[0][0][2] + DLUL[1][1][2] + DLUL[2][2][2],
};
real DUUL[3][3][3] = {
{{DLUL[0][0][0] * gammaUxx + DLUL[1][0][0] * gammaUxy + DLUL[2][0][0] * gammaUxz,
DLUL[0][0][1] * gammaUxx + DLUL[1][0][1] * gammaUxy + DLUL[2][0][1] * gammaUxz,
DLUL[0][0][2] * gammaUxx + DLUL[1][0][2] * gammaUxy + DLUL[2][0][2] * gammaUxz,
},{DLUL[0][1][0] * gammaUxx + DLUL[1][1][0] * gammaUxy + DLUL[2][1][0] * gammaUxz,
DLUL[0][1][1] * gammaUxx + DLUL[1][1][1] * gammaUxy + DLUL[2][1][1] * gammaUxz,
DLUL[0][1][2] * gammaUxx + DLUL[1][1][2] * gammaUxy + DLUL[2][1][2] * gammaUxz,
},{DLUL[0][2][0] * gammaUxx + DLUL[1][2][0] * gammaUxy + DLUL[2][2][0] * gammaUxz,
DLUL[0][2][1] * gammaUxx + DLUL[1][2][1] * gammaUxy + DLUL[2][2][1] * gammaUxz,
DLUL[0][2][2] * gammaUxx + DLUL[1][2][2] * gammaUxy + DLUL[2][2][2] * gammaUxz,
},},{{DLUL[0][0][0] * gammaUxy + DLUL[1][0][0] * gammaUyy + DLUL[2][0][0] * gammaUyz,
DLUL[0][0][1] * gammaUxy + DLUL[1][0][1] * gammaUyy + DLUL[2][0][1] * gammaUyz,
DLUL[0][0][2] * gammaUxy + DLUL[1][0][2] * gammaUyy + DLUL[2][0][2] * gammaUyz,
},{DLUL[0][1][0] * gammaUxy + DLUL[1][1][0] * gammaUyy + DLUL[2][1][0] * gammaUyz,
DLUL[0][1][1] * gammaUxy + DLUL[1][1][1] * gammaUyy + DLUL[2][1][1] * gammaUyz,
DLUL[0][1][2] * gammaUxy + DLUL[1][1][2] * gammaUyy + DLUL[2][1][2] * gammaUyz,
},{DLUL[0][2][0] * gammaUxy + DLUL[1][2][0] * gammaUyy + DLUL[2][2][0] * gammaUyz,
DLUL[0][2][1] * gammaUxy + DLUL[1][2][1] * gammaUyy + DLUL[2][2][1] * gammaUyz,
DLUL[0][2][2] * gammaUxy + DLUL[1][2][2] * gammaUyy + DLUL[2][2][2] * gammaUyz,
},},{{DLUL[0][0][0] * gammaUxz + DLUL[1][0][0] * gammaUyz + DLUL[2][0][0] * gammaUzz,
DLUL[0][0][1] * gammaUxz + DLUL[1][0][1] * gammaUyz + DLUL[2][0][1] * gammaUzz,
DLUL[0][0][2] * gammaUxz + DLUL[1][0][2] * gammaUyz + DLUL[2][0][2] * gammaUzz,
},{DLUL[0][1][0] * gammaUxz + DLUL[1][1][0] * gammaUyz + DLUL[2][1][0] * gammaUzz,
DLUL[0][1][1] * gammaUxz + DLUL[1][1][1] * gammaUyz + DLUL[2][1][1] * gammaUzz,
DLUL[0][1][2] * gammaUxz + DLUL[1][1][2] * gammaUyz + DLUL[2][1][2] * gammaUzz,
},{DLUL[0][2][0] * gammaUxz + DLUL[1][2][0] * gammaUyz + DLUL[2][2][0] * gammaUzz,
DLUL[0][2][1] * gammaUxz + DLUL[1][2][1] * gammaUyz + DLUL[2][2][1] * gammaUzz,
DLUL[0][2][2] * gammaUxz + DLUL[1][2][2] * gammaUyz + DLUL[2][2][2] * gammaUzz,
},},};
real D12SymLL[6] = {
U->d_xxx * DUUL[0][0][0] + U->d_xxy * DUUL[0][1][0] + U->d_xxz * DUUL[0][2][0] + U->d_yxx * DUUL[1][0][0] + U->d_yxy * DUUL[1][1][0] + U->d_yxz * DUUL[1][2][0] + U->d_zxx * DUUL[2][0][0] + U->d_zxy * DUUL[2][1][0] + U->d_zxz * DUUL[2][2][0],
U->d_xxy * DUUL[0][0][0] + U->d_xyy * DUUL[0][1][0] + U->d_xyz * DUUL[0][2][0] + U->d_yxy * DUUL[1][0][0] + U->d_yyy * DUUL[1][1][0] + U->d_yyz * DUUL[1][2][0] + U->d_zxy * DUUL[2][0][0] + U->d_zyy * DUUL[2][1][0] + U->d_zyz * DUUL[2][2][0],
U->d_xxz * DUUL[0][0][0] + U->d_xyz * DUUL[0][1][0] + U->d_xzz * DUUL[0][2][0] + U->d_yxz * DUUL[1][0][0] + U->d_yyz * DUUL[1][1][0] + U->d_yzz * DUUL[1][2][0] + U->d_zxz * DUUL[2][0][0] + U->d_zyz * DUUL[2][1][0] + U->d_zzz * DUUL[2][2][0],
U->d_xxy * DUUL[0][0][1] + U->d_xyy * DUUL[0][1][1] + U->d_xyz * DUUL[0][2][1] + U->d_yxy * DUUL[1][0][1] + U->d_yyy * DUUL[1][1][1] + U->d_yyz * DUUL[1][2][1] + U->d_zxy * DUUL[2][0][1] + U->d_zyy * DUUL[2][1][1] + U->d_zyz * DUUL[2][2][1],
U->d_xxz * DUUL[0][0][1] + U->d_xyz * DUUL[0][1][1] + U->d_xzz * DUUL[0][2][1] + U->d_yxz * DUUL[1][0][1] + U->d_yyz * DUUL[1][1][1] + U->d_yzz * DUUL[1][2][1] + U->d_zxz * DUUL[2][0][1] + U->d_zyz * DUUL[2][1][1] + U->d_zzz * DUUL[2][2][1],
U->d_xxz * DUUL[0][0][2] + U->d_xyz * DUUL[0][1][2] + U->d_xzz * DUUL[0][2][2] + U->d_yxz * DUUL[1][0][2] + U->d_yyz * DUUL[1][1][2] + U->d_yzz * DUUL[1][2][2] + U->d_zxz * DUUL[2][0][2] + U->d_zyz * DUUL[2][1][2] + U->d_zzz * DUUL[2][2][2],
};
real GammaLSymLL[3][6] = {
{U->d_xxx,
U->d_yxx,
U->d_zxx,
((2 * U->d_yxy) - U->d_xyy),
(U->d_zxy + (U->d_yxz - U->d_xyz)),
((2 * U->d_zxz) - U->d_xzz),
},{((2 * U->d_xxy) - U->d_yxx),
U->d_xyy,
(U->d_zxy + (U->d_xyz - U->d_yxz)),
U->d_yyy,
U->d_zyy,
((2 * U->d_zyz) - U->d_yzz),
},{((2 * U->d_xxz) - U->d_zxx),
(U->d_yxz + (U->d_xyz - U->d_zxy)),
U->d_xzz,
((2 * U->d_yyz) - U->d_zyy),
U->d_yzz,
U->d_zzz,
},};
real GammaUSymLL[3][6] = {
{gammaUxx * GammaLSymLL[0][0] + gammaUxy * GammaLSymLL[1][0] + gammaUxz * GammaLSymLL[2][0],
gammaUxx * GammaLSymLL[0][1] + gammaUxy * GammaLSymLL[1][1] + gammaUxz * GammaLSymLL[2][1],
gammaUxx * GammaLSymLL[0][2] + gammaUxy * GammaLSymLL[1][2] + gammaUxz * GammaLSymLL[2][2],
gammaUxx * GammaLSymLL[0][3] + gammaUxy * GammaLSymLL[1][3] + gammaUxz * GammaLSymLL[2][3],
gammaUxx * GammaLSymLL[0][4] + gammaUxy * GammaLSymLL[1][4] + gammaUxz * GammaLSymLL[2][4],
gammaUxx * GammaLSymLL[0][5] + gammaUxy * GammaLSymLL[1][5] + gammaUxz * GammaLSymLL[2][5],
},{gammaUxy * GammaLSymLL[0][0] + gammaUyy * GammaLSymLL[1][0] + gammaUyz * GammaLSymLL[2][0],
gammaUxy * GammaLSymLL[0][1] + gammaUyy * GammaLSymLL[1][1] + gammaUyz * GammaLSymLL[2][1],
gammaUxy * GammaLSymLL[0][2] + gammaUyy * GammaLSymLL[1][2] + gammaUyz * GammaLSymLL[2][2],
gammaUxy * GammaLSymLL[0][3] + gammaUyy * GammaLSymLL[1][3] + gammaUyz * GammaLSymLL[2][3],
gammaUxy * GammaLSymLL[0][4] + gammaUyy * GammaLSymLL[1][4] + gammaUyz * GammaLSymLL[2][4],
gammaUxy * GammaLSymLL[0][5] + gammaUyy * GammaLSymLL[1][5] + gammaUyz * GammaLSymLL[2][5],
},{gammaUxz * GammaLSymLL[0][0] + gammaUyz * GammaLSymLL[1][0] + gammaUzz * GammaLSymLL[2][0],
gammaUxz * GammaLSymLL[0][1] + gammaUyz * GammaLSymLL[1][1] + gammaUzz * GammaLSymLL[2][1],
gammaUxz * GammaLSymLL[0][2] + gammaUyz * GammaLSymLL[1][2] + gammaUzz * GammaLSymLL[2][2],
gammaUxz * GammaLSymLL[0][3] + gammaUyz * GammaLSymLL[1][3] + gammaUzz * GammaLSymLL[2][3],
gammaUxz * GammaLSymLL[0][4] + gammaUyz * GammaLSymLL[1][4] + gammaUzz * GammaLSymLL[2][4],
gammaUxz * GammaLSymLL[0][5] + gammaUyz * GammaLSymLL[1][5] + gammaUzz * GammaLSymLL[2][5],
},};
real Gamma3L[3] = {
GammaUSymLL[0][0] + GammaUSymLL[1][1] + GammaUSymLL[2][2],
GammaUSymLL[0][1] + GammaUSymLL[1][3] + GammaUSymLL[2][4],
GammaUSymLL[0][2] + GammaUSymLL[1][4] + GammaUSymLL[2][5],
};
real Gamma31SymLL[6] = {
Gamma3L[0] * GammaUSymLL[0][0] + Gamma3L[1] * GammaUSymLL[1][0] + Gamma3L[2] * GammaUSymLL[2][0],
Gamma3L[0] * GammaUSymLL[0][1] + Gamma3L[1] * GammaUSymLL[1][1] + Gamma3L[2] * GammaUSymLL[2][1],
Gamma3L[0] * GammaUSymLL[0][2] + Gamma3L[1] * GammaUSymLL[1][2] + Gamma3L[2] * GammaUSymLL[2][2],
Gamma3L[0] * GammaUSymLL[0][3] + Gamma3L[1] * GammaUSymLL[1][3] + Gamma3L[2] * GammaUSymLL[2][3],
Gamma3L[0] * GammaUSymLL[0][4] + Gamma3L[1] * GammaUSymLL[1][4] + Gamma3L[2] * GammaUSymLL[2][4],
Gamma3L[0] * GammaUSymLL[0][5] + Gamma3L[1] * GammaUSymLL[1][5] + Gamma3L[2] * GammaUSymLL[2][5],
};
real GammaLUL[3][3][3] = {
{{gammaUxx * GammaLSymLL[0][0] + gammaUxy * GammaLSymLL[0][1] + gammaUxz * GammaLSymLL[0][2],
gammaUxx * GammaLSymLL[0][1] + gammaUxy * GammaLSymLL[0][3] + gammaUxz * GammaLSymLL[0][4],
gammaUxx * GammaLSymLL[0][2] + gammaUxy * GammaLSymLL[0][4] + gammaUxz * GammaLSymLL[0][5],
},{gammaUxy * GammaLSymLL[0][0] + gammaUyy * GammaLSymLL[0][1] + gammaUyz * GammaLSymLL[0][2],
gammaUxy * GammaLSymLL[0][1] + gammaUyy * GammaLSymLL[0][3] + gammaUyz * GammaLSymLL[0][4],
gammaUxy * GammaLSymLL[0][2] + gammaUyy * GammaLSymLL[0][4] + gammaUyz * GammaLSymLL[0][5],
},{gammaUxz * GammaLSymLL[0][0] + gammaUyz * GammaLSymLL[0][1] + gammaUzz * GammaLSymLL[0][2],
gammaUxz * GammaLSymLL[0][1] + gammaUyz * GammaLSymLL[0][3] + gammaUzz * GammaLSymLL[0][4],
gammaUxz * GammaLSymLL[0][2] + gammaUyz * GammaLSymLL[0][4] + gammaUzz * GammaLSymLL[0][5],
},},{{gammaUxx * GammaLSymLL[1][0] + gammaUxy * GammaLSymLL[1][1] + gammaUxz * GammaLSymLL[1][2],
gammaUxx * GammaLSymLL[1][1] + gammaUxy * GammaLSymLL[1][3] + gammaUxz * GammaLSymLL[1][4],
gammaUxx * GammaLSymLL[1][2] + gammaUxy * GammaLSymLL[1][4] + gammaUxz * GammaLSymLL[1][5],
},{gammaUxy * GammaLSymLL[1][0] + gammaUyy * GammaLSymLL[1][1] + gammaUyz * GammaLSymLL[1][2],
gammaUxy * GammaLSymLL[1][1] + gammaUyy * GammaLSymLL[1][3] + gammaUyz * GammaLSymLL[1][4],
gammaUxy * GammaLSymLL[1][2] + gammaUyy * GammaLSymLL[1][4] + gammaUyz * GammaLSymLL[1][5],
},{gammaUxz * GammaLSymLL[1][0] + gammaUyz * GammaLSymLL[1][1] + gammaUzz * GammaLSymLL[1][2],
gammaUxz * GammaLSymLL[1][1] + gammaUyz * GammaLSymLL[1][3] + gammaUzz * GammaLSymLL[1][4],
gammaUxz * GammaLSymLL[1][2] + gammaUyz * GammaLSymLL[1][4] + gammaUzz * GammaLSymLL[1][5],
},},{{gammaUxx * GammaLSymLL[2][0] + gammaUxy * GammaLSymLL[2][1] + gammaUxz * GammaLSymLL[2][2],
gammaUxx * GammaLSymLL[2][1] + gammaUxy * GammaLSymLL[2][3] + gammaUxz * GammaLSymLL[2][4],
gammaUxx * GammaLSymLL[2][2] + gammaUxy * GammaLSymLL[2][4] + gammaUxz * GammaLSymLL[2][5],
},{gammaUxy * GammaLSymLL[2][0] + gammaUyy * GammaLSymLL[2][1] + gammaUyz * GammaLSymLL[2][2],
gammaUxy * GammaLSymLL[2][1] + gammaUyy * GammaLSymLL[2][3] + gammaUyz * GammaLSymLL[2][4],
gammaUxy * GammaLSymLL[2][2] + gammaUyy * GammaLSymLL[2][4] + gammaUyz * GammaLSymLL[2][5],
},{gammaUxz * GammaLSymLL[2][0] + gammaUyz * GammaLSymLL[2][1] + gammaUzz * GammaLSymLL[2][2],
gammaUxz * GammaLSymLL[2][1] + gammaUyz * GammaLSymLL[2][3] + gammaUzz * GammaLSymLL[2][4],
gammaUxz * GammaLSymLL[2][2] + gammaUyz * GammaLSymLL[2][4] + gammaUzz * GammaLSymLL[2][5],
},},};
real GammaLSymUU[3][6] = {
{gammaUxx * GammaLUL[0][0][0] + gammaUxy * GammaLUL[0][0][1] + gammaUxz * GammaLUL[0][0][2],
gammaUxy * GammaLUL[0][0][0] + gammaUyy * GammaLUL[0][0][1] + gammaUyz * GammaLUL[0][0][2],
gammaUxz * GammaLUL[0][0][0] + gammaUyz * GammaLUL[0][0][1] + gammaUzz * GammaLUL[0][0][2],
gammaUxy * GammaLUL[0][1][0] + gammaUyy * GammaLUL[0][1][1] + gammaUyz * GammaLUL[0][1][2],
gammaUxz * GammaLUL[0][1][0] + gammaUyz * GammaLUL[0][1][1] + gammaUzz * GammaLUL[0][1][2],
gammaUxz * GammaLUL[0][2][0] + gammaUyz * GammaLUL[0][2][1] + gammaUzz * GammaLUL[0][2][2],
},{gammaUxx * GammaLUL[1][0][0] + gammaUxy * GammaLUL[1][0][1] + gammaUxz * GammaLUL[1][0][2],
gammaUxy * GammaLUL[1][0][0] + gammaUyy * GammaLUL[1][0][1] + gammaUyz * GammaLUL[1][0][2],
gammaUxz * GammaLUL[1][0][0] + gammaUyz * GammaLUL[1][0][1] + gammaUzz * GammaLUL[1][0][2],
gammaUxy * GammaLUL[1][1][0] + gammaUyy * GammaLUL[1][1][1] + gammaUyz * GammaLUL[1][1][2],
gammaUxz * GammaLUL[1][1][0] + gammaUyz * GammaLUL[1][1][1] + gammaUzz * GammaLUL[1][1][2],
gammaUxz * GammaLUL[1][2][0] + gammaUyz * GammaLUL[1][2][1] + gammaUzz * GammaLUL[1][2][2],
},{gammaUxx * GammaLUL[2][0][0] + gammaUxy * GammaLUL[2][0][1] + gammaUxz * GammaLUL[2][0][2],
gammaUxy * GammaLUL[2][0][0] + gammaUyy * GammaLUL[2][0][1] + gammaUyz * GammaLUL[2][0][2],
gammaUxz * GammaLUL[2][0][0] + gammaUyz * GammaLUL[2][0][1] + gammaUzz * GammaLUL[2][0][2],
gammaUxy * GammaLUL[2][1][0] + gammaUyy * GammaLUL[2][1][1] + gammaUyz * GammaLUL[2][1][2],
gammaUxz * GammaLUL[2][1][0] + gammaUyz * GammaLUL[2][1][1] + gammaUzz * GammaLUL[2][1][2],
gammaUxz * GammaLUL[2][2][0] + gammaUyz * GammaLUL[2][2][1] + gammaUzz * GammaLUL[2][2][2],
},};
real Gamma11SymLL[6] = {
GammaLSymLL[0][0] * GammaLSymUU[0][0] + GammaLSymLL[0][1] * GammaLSymUU[0][1] + GammaLSymLL[0][2] * GammaLSymUU[0][2] + GammaLSymLL[0][1] * GammaLSymUU[0][1] + GammaLSymLL[0][3] * GammaLSymUU[0][3] + GammaLSymLL[0][4] * GammaLSymUU[0][4] + GammaLSymLL[0][2] * GammaLSymUU[0][2] + GammaLSymLL[0][4] * GammaLSymUU[0][4] + GammaLSymLL[0][5] * GammaLSymUU[0][5],
GammaLSymLL[0][0] * GammaLSymUU[1][0] + GammaLSymLL[0][1] * GammaLSymUU[1][1] + GammaLSymLL[0][2] * GammaLSymUU[1][2] + GammaLSymLL[0][1] * GammaLSymUU[1][1] + GammaLSymLL[0][3] * GammaLSymUU[1][3] + GammaLSymLL[0][4] * GammaLSymUU[1][4] + GammaLSymLL[0][2] * GammaLSymUU[1][2] + GammaLSymLL[0][4] * GammaLSymUU[1][4] + GammaLSymLL[0][5] * GammaLSymUU[1][5],
GammaLSymLL[0][0] * GammaLSymUU[2][0] + GammaLSymLL[0][1] * GammaLSymUU[2][1] + GammaLSymLL[0][2] * GammaLSymUU[2][2] + GammaLSymLL[0][1] * GammaLSymUU[2][1] + GammaLSymLL[0][3] * GammaLSymUU[2][3] + GammaLSymLL[0][4] * GammaLSymUU[2][4] + GammaLSymLL[0][2] * GammaLSymUU[2][2] + GammaLSymLL[0][4] * GammaLSymUU[2][4] + GammaLSymLL[0][5] * GammaLSymUU[2][5],
GammaLSymLL[1][0] * GammaLSymUU[1][0] + GammaLSymLL[1][1] * GammaLSymUU[1][1] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][1] * GammaLSymUU[1][1] + GammaLSymLL[1][3] * GammaLSymUU[1][3] + GammaLSymLL[1][4] * GammaLSymUU[1][4] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][4] * GammaLSymUU[1][4] + GammaLSymLL[1][5] * GammaLSymUU[1][5],
GammaLSymLL[1][0] * GammaLSymUU[2][0] + GammaLSymLL[1][1] * GammaLSymUU[2][1] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][1] * GammaLSymUU[2][1] + GammaLSymLL[1][3] * GammaLSymUU[2][3] + GammaLSymLL[1][4] * GammaLSymUU[2][4] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][4] * GammaLSymUU[2][4] + GammaLSymLL[1][5] * GammaLSymUU[2][5],
GammaLSymLL[2][0] * GammaLSymUU[2][0] + GammaLSymLL[2][1] * GammaLSymUU[2][1] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][1] * GammaLSymUU[2][1] + GammaLSymLL[2][3] * GammaLSymUU[2][3] + GammaLSymLL[2][4] * GammaLSymUU[2][4] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][4] * GammaLSymUU[2][4] + GammaLSymLL[2][5] * GammaLSymUU[2][5],
};
real ADL[3] = {
U->a_x - 2 * D3L[0],
U->a_y - 2 * D3L[1],
U->a_z - 2 * D3L[2],
};
real ADU[3] = {
gammaUxx * ADL[0] + gammaUxy * ADL[1] + gammaUxz * ADL[2],
gammaUxy * ADL[0] + gammaUyy * ADL[1] + gammaUyz * ADL[2],
gammaUxz * ADL[0] + gammaUyz * ADL[1] + gammaUzz * ADL[2],
};
real ADDSymLL[6] = {
ADU[0] * (2 * U->d_xxx) + ADU[1] * (2 * U->d_xxy) + ADU[2] * (2 * U->d_xxz),
ADU[0] * (U->d_xxy + U->d_yxx) + ADU[1] * (U->d_xyy + U->d_yxy) + ADU[2] * (U->d_xyz + U->d_yxz),
ADU[0] * (U->d_xxz + U->d_zxx) + ADU[1] * (U->d_xyz + U->d_zxy) + ADU[2] * (U->d_xzz + U->d_zxz),
ADU[0] * (2 * U->d_yxy) + ADU[1] * (2 * U->d_yyy) + ADU[2] * (2 * U->d_yyz),
ADU[0] * (U->d_yxz + U->d_zxy) + ADU[1] * (U->d_yyz + U->d_zyy) + ADU[2] * (U->d_yzz + U->d_zyz),
ADU[0] * (2 * U->d_zxz) + ADU[1] * (2 * U->d_zyz) + ADU[2] * (2 * U->d_zzz),
};
real R4SymLL[6] = {
8. * M_PI * ((density + pressure) * vel4_.x * vel4_.x + .5 * (density - pressure) * U->gamma_xx),
8. * M_PI * ((density + pressure) * vel4_.x * vel4_.y + .5 * (density - pressure) * U->gamma_xy),
8. * M_PI * ((density + pressure) * vel4_.x * vel4_.z + .5 * (density - pressure) * U->gamma_xz),
8. * M_PI * ((density + pressure) * vel4_.y * vel4_.y + .5 * (density - pressure) * U->gamma_yy),
8. * M_PI * ((density + pressure) * vel4_.y * vel4_.z + .5 * (density - pressure) * U->gamma_yz),
8. * M_PI * ((density + pressure) * vel4_.z * vel4_.z + .5 * (density - pressure) * U->gamma_zz),
};
real SSymLL[6] = {
-R4SymLL[0] + trK * U->K_xx - 2 * KSqSymLL[0] + 4 * D12SymLL[0] + Gamma31SymLL[0] - Gamma11SymLL[0] + ADDSymLL[0] + (U->a_x * ((2 * U->V_x) - D1L[0])),
-R4SymLL[1] + trK * U->K_xy - 2 * KSqSymLL[1] + 4 * D12SymLL[1] + Gamma31SymLL[1] - Gamma11SymLL[1] + ADDSymLL[1] + ((((2 * U->a_y * U->V_x) - (U->a_y * D1L[0])) + ((2 * U->a_x * U->V_y) - (U->a_x * D1L[1]))) / 2),
-R4SymLL[2] + trK * U->K_xz - 2 * KSqSymLL[2] + 4 * D12SymLL[2] + Gamma31SymLL[2] - Gamma11SymLL[2] + ADDSymLL[2] + ((((2 * U->a_z * U->V_x) - (U->a_z * D1L[0])) + ((2 * U->a_x * U->V_z) - (U->a_x * D1L[2]))) / 2),
-R4SymLL[3] + trK * U->K_yy - 2 * KSqSymLL[3] + 4 * D12SymLL[3] + Gamma31SymLL[3] - Gamma11SymLL[3] + ADDSymLL[3] + (U->a_y * ((2 * U->V_y) - D1L[1])),
-R4SymLL[4] + trK * U->K_yz - 2 * KSqSymLL[4] + 4 * D12SymLL[4] + Gamma31SymLL[4] - Gamma11SymLL[4] + ADDSymLL[4] + ((((2 * U->a_z * U->V_y) - (U->a_z * D1L[1])) + ((2 * U->a_y * U->V_z) - (U->a_y * D1L[2]))) / 2),
-R4SymLL[5] + trK * U->K_zz - 2 * KSqSymLL[5] + 4 * D12SymLL[5] + Gamma31SymLL[5] - Gamma11SymLL[5] + ADDSymLL[5] + (U->a_z * ((2 * U->V_z) - D1L[2])),
};
real GU0L[3] = {
8. * M_PI * ((density + pressure) * vel4_.w * vel4_.x + pressure * beta_.x),
8. * M_PI * ((density + pressure) * vel4_.w * vel4_.y + pressure * beta_.y),
8. * M_PI * ((density + pressure) * vel4_.w * vel4_.z + pressure * beta_.z),
};
real AKL[3] = {
U->a_x * KUL[0][0] + U->a_y * KUL[1][0] + U->a_z * KUL[2][0],
U->a_x * KUL[0][1] + U->a_y * KUL[1][1] + U->a_z * KUL[2][1],
U->a_x * KUL[0][2] + U->a_y * KUL[1][2] + U->a_z * KUL[2][2],
};
real K12D23L[3] = {
KUL[0][0] * DLUL[0][0][0] +KUL[0][1] * DLUL[0][1][0] +KUL[0][2] * DLUL[0][2][0] + KUL[1][0] * DLUL[0][0][1] +KUL[1][1] * DLUL[0][1][1] +KUL[1][2] * DLUL[0][2][1] + KUL[2][0] * DLUL[0][0][2] +KUL[2][1] * DLUL[0][1][2] +KUL[2][2] * DLUL[0][2][2],
KUL[0][0] * DLUL[1][0][0] +KUL[0][1] * DLUL[1][1][0] +KUL[0][2] * DLUL[1][2][0] + KUL[1][0] * DLUL[1][0][1] +KUL[1][1] * DLUL[1][1][1] +KUL[1][2] * DLUL[1][2][1] + KUL[2][0] * DLUL[1][0][2] +KUL[2][1] * DLUL[1][1][2] +KUL[2][2] * DLUL[1][2][2],
KUL[0][0] * DLUL[2][0][0] +KUL[0][1] * DLUL[2][1][0] +KUL[0][2] * DLUL[2][2][0] + KUL[1][0] * DLUL[2][0][1] +KUL[1][1] * DLUL[2][1][1] +KUL[1][2] * DLUL[2][2][1] + KUL[2][0] * DLUL[2][0][2] +KUL[2][1] * DLUL[2][1][2] +KUL[2][2] * DLUL[2][2][2],
};
real KD23L[3] = {
KUL[0][0] * D1L[0] + KUL[1][0] * D1L[1] + KUL[2][0] * D1L[2],
KUL[0][1] * D1L[0] + KUL[1][1] * D1L[1] + KUL[2][1] * D1L[2],
KUL[0][2] * D1L[0] + KUL[1][2] * D1L[1] + KUL[2][2] * D1L[2],
};
real K12D12L[3] = {
KUL[0][0] * DLUL[0][0][0] + KUL[0][1] * DLUL[0][1][0] + KUL[0][2] * DLUL[0][2][0] + KUL[1][0] * DLUL[1][0][0] + KUL[1][1] * DLUL[1][1][0] + KUL[1][2] * DLUL[1][2][0] + KUL[2][0] * DLUL[2][0][0] + KUL[2][1] * DLUL[2][1][0] + KUL[2][2] * DLUL[2][2][0],
KUL[0][0] * DLUL[0][0][1] + KUL[0][1] * DLUL[0][1][1] + KUL[0][2] * DLUL[0][2][1] + KUL[1][0] * DLUL[1][0][1] + KUL[1][1] * DLUL[1][1][1] + KUL[1][2] * DLUL[1][2][1] + KUL[2][0] * DLUL[2][0][1] + KUL[2][1] * DLUL[2][1][1] + KUL[2][2] * DLUL[2][2][1],
KUL[0][0] * DLUL[0][0][2] + KUL[0][1] * DLUL[0][1][2] + KUL[0][2] * DLUL[0][2][2] + KUL[1][0] * DLUL[1][0][2] + KUL[1][1] * DLUL[1][1][2] + KUL[1][2] * DLUL[1][2][2] + KUL[2][0] * DLUL[2][0][2] + KUL[2][1] * DLUL[2][1][2] + KUL[2][2] * DLUL[2][2][2],
};
real KD12L[3] = {
KUL[0][0] * D3L[0] + KUL[1][0] * D3L[1] + KUL[2][0] * D3L[2],
KUL[0][1] * D3L[0] + KUL[1][1] * D3L[1] + KUL[2][1] * D3L[2],
KUL[0][2] * D3L[0] + KUL[1][2] * D3L[1] + KUL[2][2] * D3L[2],
};
real PL[3] = {
GU0L[0] + AKL[0] - U->a_x * trK + K12D23L[0] + KD23L[0] - 2 * K12D12L[0] + 2 * KD12L[0],
GU0L[1] + AKL[1] - U->a_y * trK + K12D23L[1] + KD23L[1] - 2 * K12D12L[1] + 2 * KD12L[1],
GU0L[2] + AKL[2] - U->a_z * trK + K12D23L[2] + KD23L[2] - 2 * K12D12L[2] + 2 * KD12L[2],
};

	/*
	TODO now that alpha and gamma are moved from the flux eigensystem
	how would we still incorporate the shift terms with them?
	the easy solution is to re-add them back in (even though shift computation is supposed to go on apart from the eigensystem)
	typically, when they are included, they get a flux between cells equal to the shift at that cell
	maybe that has to be done, even if they are not a part of the eigensystem?
	then again, maybe I should be keeping alpha and gamma in the eigensystem,
	and maybe aux vars like the density should be in the eigensystem as well,
	all for no reason more than to be influenced by the shift
	
	Then again, maybe this is an argument for the solver to specify the flux vector size 
	-- especially if it is allowed a custom RoeFluxDeriv function.
	*/

	deriv->alpha += -U->alpha * U->alpha * f * trK;
	deriv->gamma_xx += -2. * U->alpha * U->K_xx;
	deriv->gamma_xy += -2. * U->alpha * U->K_xy;
	deriv->gamma_xz += -2. * U->alpha * U->K_xz;
	deriv->gamma_yy += -2. * U->alpha * U->K_yy;
	deriv->gamma_yz += -2. * U->alpha * U->K_yz;
	deriv->gamma_zz += -2. * U->alpha * U->K_zz;
	
	deriv->K_xx += U->alpha * SSymLL[0];
	deriv->K_xy += U->alpha * SSymLL[1];
	deriv->K_xz += U->alpha * SSymLL[2];
	deriv->K_yy += U->alpha * SSymLL[3];
	deriv->K_yz += U->alpha * SSymLL[4];
	deriv->K_zz += U->alpha * SSymLL[5];
	deriv->V_x += U->alpha * PL[0];
	deriv->V_y += U->alpha * PL[1];
	deriv->V_z += U->alpha * PL[2];
}

// the 1D version has no problems, but at 2D we get instabilities ... 
__kernel void constrain(
	__global cons_t* UBuf
) {
#if 0	//use constraints at all?
	
	int4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	if (i.x < 2 || i.x >= SIZE_X - 2 
#if DIM > 1
		|| i.y < 2 || i.y >= SIZE_Y - 2 
#endif
#if DIM > 2
		|| i.z < 2 || i.z >= SIZE_Z - 2
#endif
	) {
		return;
	}
	
	int index = INDEXV(i);
	__global cons_t* U = UBuf + index;

	real gamma = symMatDet_prefix(U->gamma_);
	real gammaInv[6];
	symMatInv_prefix(gammaInv, gamma, U->gamma_);
	real gammaUxx = gammaInv[0], gammaUxy = gammaInv[1], gammaUxz = gammaInv[2], gammaUyy = gammaInv[3], gammaUyz = gammaInv[4], gammaUzz = gammaInv[5];

	real D3_D1_x = 
		(gammaUxy * U->d_xxy)
		+ (gammaUxz * U->d_xxz)
		+ (gammaUyy * U->d_xyy)
		+ (2. * gammaUyz * U->d_xyz)
		+ (gammaUzz * U->d_xzz)
		- (gammaUxy * U->d_yxx)
		- (gammaUxz * U->d_zxx)
		- (gammaUyy * U->d_yxy)
		- (gammaUyz * U->d_zxy)
		- (gammaUyz * U->d_yxz)
		- (gammaUzz * U->d_zxz);
	real D3_D1_y = 
		(gammaUxx * U->d_yxx)
		+ (gammaUxy * U->d_yxy)
		+ (2. * gammaUxz * U->d_yxz)
		+ (gammaUyz * U->d_yyz)
		+ (gammaUzz * U->d_yzz)
		- (gammaUxx * U->d_xxy)
		- (gammaUxz * U->d_zxy)
		- (gammaUxy * U->d_xyy)
		- (gammaUyz * U->d_zyy)
		- (gammaUxz * U->d_xyz)
		- (gammaUzz * U->d_zyz);
	real D3_D1_z = 
		(gammaUxx * U->d_zxx)
		+ (2. * gammaUxy * U->d_zxy)
		+ (gammaUxz * U->d_zxz)
		+ (gammaUyy * U->d_zyy)
		+ (gammaUyz * U->d_zyz)
		- (gammaUxx * U->d_xxz)
		- (gammaUxy * U->d_yxz)
		- (gammaUxy * U->d_xyz)
		- (gammaUyy * U->d_yyz)
		- (gammaUxz * U->d_xzz)
		- (gammaUyz * U->d_yzz);

#if 0	//directly assign V_i's
	U->V_x = D3_D1_x;
	U->V_y = D3_D1_y;
	U->V_z = D3_D1_z;
#endif
#if 0	//linearly project out the [V_i, U->d_ijk] vector
#endif
#if 0	//do a single gradient descent step
/*
V_i = d_ik^k - d^k_ki
f_i = V_i - (d_ijk - d_jki) gamma^jk
*/
	U->V_x -= epsilon;
	U->V_y -= epsilon;
#endif

#endif
}
