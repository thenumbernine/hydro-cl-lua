#define sqrt_1_2 0.70710678118654757273731092936941422522068023681641

real ESq(cons_t U) { 
	return (U.epsEx * U.epsEx + U.epsEy * U.epsEy + U.epsEz * U.epsEz) / (eps0 * eps0);
}

real BSq(cons_t U) {
	return U.Bx * U.Bx + U.By * U.By + U.Bz * U.Bz;
}

real calcDisplayVar_UBuf(
	int displayVar, 
	const __global real* U_
) {
	const __global cons_t* U = (const __global cons_t*)U_;
	switch (displayVar) {
	case display_U_Ex: return U->epsEx / eps0;
	case display_U_Ey: return U->epsEy / eps0;
	case display_U_Ez: return U->epsEz / eps0;
	case display_U_E: return sqrt(ESq(*U));
	case display_U_Bx: return U->Bx;
	case display_U_By: return U->By;
	case display_U_Bz: return U->Bz;
	case display_U_B: return sqrt(BSq(*U));
	case display_U_energy: return .5 * (ESq(*U) * eps0 + BSq(*U) / mu0);
	}
	return 0;
}

real calcEigenvalue() { 
	return 1./(sqrt_eps0 * sqrt_mu0);
}

range_t calcCellMinMaxEigenvalues(const __global cons_t* U, int side) {
	real lambda = calcEigenvalue();
	return (range_t){-lambda, lambda};
}

void fill(__global real* ptr, int step, real a, real b, real c, real d, real e, real f) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
	ptr[3*step] = d;
	ptr[4*step] = e;
	ptr[5*step] = f;
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
		real lambda = calcEigenvalue();
		
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		fill(wave, 1, -lambda, -lambda, 0, 0, lambda, lambda);
	
		//no eigenbuf info since waves are unrelated to state
	}
}

void fluxTransform(
	real* y,
	const __global fluxXform_t* flux,
	const real* x,
	int side
) {
	//swap input dim x<->side
	real4 epsE = (real4)(x[0], x[1], x[2], 0);
	real4 B = (real4)(x[3], x[4], x[5], 0);
	
	switch (side) {
	case 0:
		y[0] = 0;
		y[1] = B.z / mu0;
		y[2] = -B.y / mu0;
		y[3] = 0;
		y[4] = -epsE.z / eps0;
		y[5] = epsE.y / eps0;
		break;
	case 1:
		y[0] = -B.z / mu0;
		y[1] = 0;
		y[2] = B.x / mu0;
		y[3] = epsE.z / eps0;
		y[4] = 0;
		y[5] = -epsE.x / eps0;
		break;
	case 2:
		y[0] = B.y / mu0;
		y[1] = -B.x / mu0;
		y[2] = 0;
		y[3] = -epsE.y / eps0;
		y[4] = epsE.x / eps0;
		y[5] = 0;
		break;
	} 
}

void eigen_leftTransform(
	real* y,
	const __global eigen_t* eigen,
	const real* x,
	int side
) {
	const real ise = sqrt_1_2/sqrt_eps0;
	const real isu = sqrt_1_2/sqrt_mu0;
	
	//swap input dim x<->side
	real4 epsE = (real4)(x[0], x[1], x[2], 0);
	real4 B = (real4)(x[3], x[4], x[5], 0);
	
	switch (side) {
	default:
	case 0:
		y[0] = epsE.z * ise + B.y * isu;
		y[1] = epsE.y * -ise + B.z * isu;
		y[2] = epsE.x * -ise + B.x * isu;
		y[3] = epsE.x * ise + B.x * isu;
		y[4] = epsE.y * ise + B.z * isu;
		y[5] = epsE.z * -ise + B.y * isu;
		break;
	}
}

void eigen_rightTransform(
	real* y,
	const __global eigen_t* eigen,
	const real* x,
	int side
) {
	const real se = sqrt_1_2 * sqrt_eps0;
	const real su = sqrt_1_2 * sqrt_mu0;

	switch (side) {
	default:
	case 0:
		y[0] = -se * x[2] + se * x[3];
		y[1] = -se * x[1] + se * x[4];
		y[2] = se * x[0] - se * x[5];
		y[3] = su * x[2] + su * x[3];
		y[4] = su * x[0] + su * x[5];
		y[5] = su * x[1] + su * x[4];
		break;
	}
}

__kernel void addSourceTerm(
	__global cons_t* derivBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	const __global cons_t* U = UBuf + index;
	__global cons_t* deriv = derivBuf + index;
	deriv->epsEx -= U->epsEx / eps0 * sigma;
	deriv->epsEy -= U->epsEy / eps0 * sigma;
	deriv->epsEz -= U->epsEz / eps0 * sigma;
}
