__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real4 x = CELL_X(i);
	real4 mids = (real).5 * (mins + maxs);
	bool lhs = x[0] < mids[0]
#if dim > 1
		&& x[1] < mids[1]
#endif
#if dim > 2
		&& x[2] < mids[2]
#endif
	;
	__global cons_t* U = UBuf + index;
	U->epsEx = 0;
	U->epsEy = 0;
	U->epsEz = 1 * eps0;
	U->Bx = 1;
	U->By = lhs ? 1 : -1;
	U->Bz = 0;
}

real ESq(cons_t U) { 
	return sqrt(U.epsEx*U.epsEx + U.epsEy*U.epsEy + U.epsEz*U.epsEz) / eps0;
}

real BSq(cons_t U) {
	return sqrt(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz);
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
	case display_U_E: return ESq(*U);
	case display_U_Bx: return U->Bx;
	case display_U_By: return U->By;
	case display_U_Bz: return U->Bz;
	case display_U_B: return BSq(*U);
	case display_U_energy: return .5 * (ESq(*U) * eps0 + BSq(*U) / mu0);
	}
	return 0;
}

real calcEigenvalue() { 
	return 1./sqrt(eps0 * mu0);
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
	__global real* fluxMatrixBuf,
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
		__global real* dF_dU = fluxMatrixBuf + numStates * numStates * intindex;;
		fill(dF_dU+0,6,	0, 0, 0, 0, 0, 0);
		fill(dF_dU+1,6,	0, 0, 0, 0, 0, 1/mu0);
		fill(dF_dU+2,6,	0, 0, 0, 0, -1/mu0, 0);
		fill(dF_dU+3,6,	0, 0, 0, 0, 0, 0);
		fill(dF_dU+4,6,	0, 0, -1/eps0, 0, 0, 0);
		fill(dF_dU+5,6,	0, 1/eps0, 0, 0, 0, 0);
	
		__global eigen_t* eigen = eigenBuf + intindex;
#define _2 0.70710678118654757273731092936941422522068023681641
	}
}

void eigen_leftTransform(
	real* y,
	const __global eigen_t* eigen,
	real* x
) {
	real se = _2 * sqrt(eps0);
	real ise = 1./se;
	real su = _2 * sqrt(mu0);
	real isu = 1./su;
	y[0] = ise * x[2] + isu * x[4];
	y[1] = -ise * x[1] + isu * x[5];
	y[2] = -ise * x[0] + isu * x[3];
	y[3] = ise * x[0] + isu * x[3];
	y[4] = ise * x[1] + isu * x[5];
	y[5] = -ise * x[2] + isu * x[4];
}

void eigen_rightTransform(
	real* y,
	const __global eigen_t* eigen,
	real* x
) {
	real se = _2 * sqrt(eps0);
	real su = _2 * sqrt(mu0);
	y[0] = -se * x[2] + se * x[3];
	y[1] = -se * x[1] + se * x[4];
	y[2] = se * x[0] - se * x[5];
	y[3] = su * x[2] + su * x[3];
	y[4] = su * x[0] + su * x[5];
	y[5] = su * x[1] + su * x[4];
}
