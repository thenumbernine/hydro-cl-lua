#define gamma_1 (gamma-1.)
#define gamma_3 (gamma-3.)

real calc_hTotal(real rho, real P, real ETotal) {
	return (P + ETotal) / rho;
}

cons_t consFromPrim(prim_t W) {
	return (cons_t){
		.rho = W.rho,
		.mx = W.rho * W.vx,
		.ETotal = .5 * W.rho * W.vx * W.vx + W.P / (gamma - 1),
	};
}

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
	real rho = 0;
	real vx = 0;
	real P = 0;

	INIT_STATE_CODE

	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .P=P});
}

prim_t primFromCons(cons_t U) {
	real EInt = U.ETotal - .5 * U.mx * U.mx / U.rho;
	return (prim_t){
		.rho = U.rho,
		.vx = U.mx / U.rho,
		.P = EInt / (gamma - 1),
	};
}

//called from calcDisplayVar
real calcDisplayVar_UBuf(int displayVar, const __global real* U_) {
	const __global cons_t* U = (const __global cons_t*)U_;
	prim_t W = primFromCons(*U);
	switch (displayVar) {
	case display_U_rho: return W.rho;
	case display_U_vx: return W.vx;
	case display_U_P: return W.P;
	case display_U_mx: return U->mx;
	case display_U_eInt: return W.P / (W.rho * gamma_1);
	case display_U_eKin: return .5 * W.vx * W.vx;
	case display_U_eTotal: return U->ETotal / W.rho;
	case display_U_EInt: return W.P / gamma_1;
	case display_U_EKin: return .5 * W.rho * W.vx * W.vx;
	case display_U_ETotal: return U->ETotal;
	case display_U_S: return W.P / pow(W.rho, (real)gamma);
	case display_U_H: return W.P * gamma / gamma_1;
	case display_U_h: return W.P * gamma / gamma_1 / W.rho;
	case display_U_HTotal: return W.P * gamma / gamma_1 + .5 * W.rho * W.vx * W.vx;
	case display_U_hTotal: return W.P * gamma / gamma_1 / W.rho + .5 * W.vx * W.vx;
	}
	return 0;
}

//called from calcDT
range_t calcCellMinMaxEigenvalues(const __global cons_t* U, int side) {
	prim_t W = primFromCons(*U);
	real Cs = sqrt(gamma * W.P / W.rho);
	return (range_t){.min=W.vx - Cs, .max=W.vx + Cs};
}

typedef struct {
	real rho, vx, hTotal;
} Roe_t;

Roe_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	prim_t WL = primFromCons(UL);
	real sqrtRhoL = sqrt(WL.rho);
	real vxL = WL.vx;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL.ETotal);

	prim_t WR = primFromCons(UR);
	real sqrtRhoR = sqrt(UR.rho);
	real vxR = WR.vx;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR.ETotal);

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);

	return (Roe_t){
		.rho = sqrtRhoL * sqrtRhoR,
		.vx = invDenom * (sqrtRhoL * vxL + sqrtRhoR * vxR),
		.hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR),
	};	
}

void fillRow(__global real* ptr, int step, real a, real b, real c) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
}
	
__kernel void calcEigenBasis(
	__global real* waveBuf,			//[volume][dim][numWaves]
	__global eigen_t* eigenBuf,		//[volume][dim]
	__global real* fluxMatrixBuf,	//[volume][dim][numStates][numStates]
	const __global cons_t *UBuf		//[volume]
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
		
		Roe_t roe = calcEigenBasisSide(UBuf[indexL], UBuf[indexR]);
		real vx = roe.vx;
		real hTotal = roe.hTotal;
		
		real vxSq = vx * vx;	
		real CsSq = gamma_1 * (hTotal - .5 * vx * vx);
		real Cs = sqrt(CsSq);
	
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		fillRow(wave, 1, vx - Cs, vx, vx + Cs);

		__global real* dF_dU = fluxMatrixBuf + numStates * numStates * intindex;
		fillRow(dF_dU+0,3,	0, 									1, 							0			);
		fillRow(dF_dU+1,3,	.5 * gamma_3 * vxSq, 				-gamma_3 * vx, 				gamma_1		);
		fillRow(dF_dU+2,3, 	vx * (.5 * gamma_1 * vxSq - hTotal), hTotal - gamma_1 * vxSq,	gamma*vx	);

		__global eigen_t* eigen = eigenBuf + intindex;
		__global real* evL = eigen->evL; 
		fillRow(evL+0, 3, (.5 * gamma_1 * vxSq + Cs * vx) / (2. * CsSq),	-(Cs + gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);
		fillRow(evL+1, 3, 1. - gamma_1 * vxSq / (2. * CsSq),				gamma_1 * vx / CsSq,				-gamma_1 / CsSq			);
		fillRow(evL+2, 3, (.5 * gamma_1 * vxSq - Cs * vx) / (2. * CsSq),	(Cs - gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);
		__global real* evR = eigen->evR;
		fillRow(evR+0, 3, 1., 				1., 		1.				);
		fillRow(evR+1, 3, vx - Cs, 			vx, 		vx + Cs			);
		fillRow(evR+2, 3, hTotal - Cs * vx, .5 * vxSq, 	hTotal + Cs * vx);
	}
}
