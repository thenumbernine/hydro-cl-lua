#define sqrt_1_2 0.70710678118654757273731092936941422522068023681641

real calcEigenvalue() { 
	return 1./(sqrt_eps0 * sqrt_mu0);
}

range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	int side
) {
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
	const real ise = sqrt_1_2 / sqrt_eps0;
	const real isu = sqrt_1_2 / sqrt_mu0;
	
	switch (side) {
	case 0:
		y[0] = x[2] *  ise + x[4] * isu;
		y[1] = x[1] * -ise + x[5] * isu;
		y[2] = x[0] * -ise + x[3] * isu;
		y[3] = x[0] *  ise + x[3] * isu;
		y[4] = x[1] *  ise + x[5] * isu;
		y[5] = x[2] * -ise + x[4] * isu;
		break;
	case 1:
		y[0] = x[0] *  ise + x[5] * isu;
		y[1] = x[2] * -ise + x[3] * isu;
		y[2] = x[1] * -ise + x[4] * isu;
		y[3] = x[1] *  ise + x[4] * isu;
		y[4] = x[2] *  ise + x[3] * isu;
		y[5] = x[0] * -ise + x[5] * isu;
		break;
	case 2:
		y[0] = x[1] *  ise + x[3] * isu;
		y[1] = x[0] * -ise + x[4] * isu;
		y[2] = x[2] * -ise + x[5] * isu;
		y[3] = x[2] *  ise + x[5] * isu;
		y[4] = x[0] *  ise + x[4] * isu;
		y[5] = x[1] * -ise + x[3] * isu;
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
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
		y[0] = se * (-x[2] + x[3]);
		y[1] = se * (-x[1] + x[4]);
		y[2] = se * (x[0] + -x[5]);
		y[3] = su * (x[2] + x[3]);
		y[4] = su * (x[0] + x[5]);
		y[5] = su * (x[1] + x[4]);
		break;
	case 1:
/*
x, -z, -y, y, z, -x
z,  x,  y, y, x,  z

1  0  0 0 0 -1
0  0 -1 1 0  0
0 -1  0 0 1  0
0  1  0 0 1  0
0  0  1 1 0  0
1  0  0 0 0  1
*/
		y[0] = se * (x[0] - x[5]);
		y[1] = se * (-x[2] + x[3]);
		y[2] = se * (-x[1] + x[4]);
		y[3] = su * (x[1] + x[4]);
		y[4] = su * (x[2] + x[3]);
		y[5] = su * (x[0] + x[5]);
		break;
	case 2:
/*
y, -x, -z, z, x, -y
x,  y,  z, z,  y,  x

0 -1  0 0 1  0
1  0  0 0 0 -1
0  0 -1 1 0  0
1  0  0 0 0  1
0  1  0 0 1  0
0  0  1 1 0  0
*/
		y[0] = se * (-x[1] + x[4]);
		y[1] = se * (x[0] - x[5]);
		y[2] = se * (-x[2] + x[3]);
		y[3] = su * (x[0] + x[5]);
		y[4] = su * (x[1] + x[4]);
		y[5] = su * (x[2] + x[3]);
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
