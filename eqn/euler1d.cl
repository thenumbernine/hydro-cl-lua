//called from calcDT
range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	int side
) {
	prim_t W = primFromCons(*U);
	real Cs = sqrt(gamma * W.P / W.rho);
	return (range_t){.min=W.vx - Cs, .max=W.vx + Cs};
}

real calc_hTotal(real rho, real P, real ETotal) {
	return (P + ETotal) / rho;
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

void fill(__global real* ptr, int step, real a, real b, real c) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
}
	
__kernel void calcEigenBasis(
	__global real* waveBuf,			//[volume][dim][numWaves]
	__global eigen_t* eigenBuf,		//[volume][dim]
	__global fluxXform_t* fluxXformBuf,	//[volume][dim]
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
		fill(wave, 1, vx - Cs, vx, vx + Cs);

		__global real* dF_dU = fluxXformBuf[intindex].A;
		fill(dF_dU+0,3,	0, 									1, 							0			);
		fill(dF_dU+1,3,	.5 * gamma_3 * vxSq, 				-gamma_3 * vx, 				gamma_1		);
		fill(dF_dU+2,3, 	vx * (.5 * gamma_1 * vxSq - hTotal), hTotal - gamma_1 * vxSq,	gamma*vx	);

		__global eigen_t* eigen = eigenBuf + intindex;
		__global real* evL = eigen->evL; 
		fill(evL+0, 3, (.5 * gamma_1 * vxSq + Cs * vx) / (2. * CsSq),	-(Cs + gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);
		fill(evL+1, 3, 1. - gamma_1 * vxSq / (2. * CsSq),				gamma_1 * vx / CsSq,				-gamma_1 / CsSq			);
		fill(evL+2, 3, (.5 * gamma_1 * vxSq - Cs * vx) / (2. * CsSq),	(Cs - gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);
		__global real* evR = eigen->evR;
		fill(evR+0, 3, 1., 				1., 		1.				);
		fill(evR+1, 3, vx - Cs, 			vx, 		vx + Cs			);
		fill(evR+2, 3, hTotal - Cs * vx, .5 * vxSq, 	hTotal + Cs * vx);
	}
}
