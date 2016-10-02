constant real gamma_1 = gamma - 1;
constant real gamma_3 = gamma - 3;

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

prim_t primFromCons(cons_t U) {
	real EInt = U.ETotal - .5 * U.mx * U.mx / U.rho;
	return (prim_t){
		.rho = U.rho,
		.vx = U.mx / U.rho,
		.P = EInt / (gamma - 1),
	};
}

//called from convertToTex
real convertToTex_UBuf(int displayVar, const __global real* U_) {
	cons_t U = *(const __global cons_t*)U_;
	prim_t W = primFromCons(U);
	real rho = W.rho;
	real vx = W.vx;
	real P = W.P;
	switch (displayVar) {
	case display_U_rho: return rho;
	case display_U_vx: return vx;
	case display_U_P: return P;
	case display_U_eInt: return P / (rho * gamma_1);
	case display_U_eKin: return .5 * vx * vx;
	case display_U_eTotal: return U.ETotal / rho;
	}
	return 0;
}

//called from calcDT
range_t calcCellMinMaxEigenvalues(cons_t U) {
	prim_t W = primFromCons(U);
	real Cs = sqrt(gamma * W.P / W.rho);
	return (range_t){
		.min = W.vx - Cs,
		.max = W.vx + Cs,
	};
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

__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real x = (real)(i.x + .5) * dx + xmin;
	UBuf[index] = consFromPrim((prim_t){
		.rho = x < 0 ? 1 : .125,
		.vx = 0,
		.P = x < 0 ? 1 : .1,
	}); 
}

void fillRow(__global real* ptr, int step, real a, real b, real c) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
}
	
__kernel void calcEigenBasis(
	__global real* waveBuf,			//[volume][dim][numWaves]
	__global real* eigenBuf,		//[volume][dim][numEigen]
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

		__global real* fluxMatrix = fluxMatrixBuf + numStates * numStates * intindex;
		fillRow(fluxMatrix+0, 3,	0, 									1, 							0			);
		fillRow(fluxMatrix+1, 3,	.5 * gamma_3 * vxSq, 				-gamma_3 * vx, 				gamma_1		);
		fillRow(fluxMatrix+2, 3, 	vx * (.5 * gamma_1 * vxSq - hTotal), hTotal - gamma_1 * vxSq,	gamma*vx	);

		__global real* evL = eigenBuf + intindex * numEigen;
		fillRow(evL+0, 3, (.5 * gamma_1 * vxSq + Cs * vx) / (2. * CsSq),	-(Cs + gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);
		fillRow(evL+1, 3, 1. - gamma_1 * vxSq / (2. * CsSq),				gamma_1 * vx / CsSq,				-gamma_1 / CsSq			);
		fillRow(evL+2, 3, (.5 * gamma_1 * vxSq - Cs * vx) / (2. * CsSq),	(Cs - gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);

		__global real* evR = evL + 3*3;
		fillRow(evR+0, 3, 1., 				1., 		1.				);
		fillRow(evR+1, 3, vx - Cs, 			vx, 		vx + Cs			);
		fillRow(evR+2, 3, hTotal - Cs * vx, .5 * vxSq, 	hTotal + Cs * vx);
	}
}


