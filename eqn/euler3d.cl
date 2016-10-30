//called from calcDT
range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	int side
) {
	prim_t W = primFromCons(*U);
	real Cs = sqrt(gamma * W.P / W.rho);
	real v = W.v[side];
	return (range_t){.min = v - Cs, .max = v + Cs};
}

typedef struct {
	real rho;
	real vx, vy, vz;
	real hTotal;
} Roe_t;

Roe_t calcEigenBasisSide(cons_t UL, cons_t UR) {
	prim_t WL = primFromCons(UL);
	real sqrtRhoL = sqrt(WL.rho);
	real vxL = WL.vx;
	real vyL = WL.vy;
	real vzL = WL.vz;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL.ETotal);

	prim_t WR = primFromCons(UR);
	real sqrtRhoR = sqrt(UR.rho);
	real vxR = WR.vx;
	real vyR = WR.vy;
	real vzR = WR.vz;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR.ETotal);

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);

	return (Roe_t){
		.rho = sqrtRhoL * sqrtRhoR,
		.vx = invDenom * (sqrtRhoL * vxL + sqrtRhoR * vxR),
		.vy = invDenom * (sqrtRhoL * vyL + sqrtRhoR * vyR),
		.vz = invDenom * (sqrtRhoL * vzL + sqrtRhoR * vzR),
		.hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR),
	};	
}

void fill(__global real* ptr, int step, real a, real b, real c, real d, real e) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
	ptr[3*step] = d;
	ptr[4*step] = e;
}
	
__kernel void calcEigenBasis(
	__global real* waveBuf,			//[volume][dim][numWaves]
	__global eigen_t* eigenBuf,		//[volume][dim]
	__global fluxXform_t* fluxXformBuf,	//[volume][dim]
	const __global cons_t* UBuf		//[volume]
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
	
		cons_t UL = UBuf[indexL];
		cons_t UR = UBuf[indexR];

		real tmp;
		tmp = UL.mx; UL.mx = UL.m[side]; UL.m[side] = tmp;
		tmp = UR.mx; UR.mx = UR.m[side]; UR.m[side] = tmp;

		Roe_t roe = calcEigenBasisSide(UL, UR);
		real vx = roe.vx;
		real vy = roe.vy;
		real vz = roe.vz;
		real hTotal = roe.hTotal;
	
		real vSq = vx * vx + vy * vy + vz * vz;
		real CsSq = gamma_1 * (hTotal - .5 * vSq);
		real Cs = sqrt(CsSq);
	
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		fill(wave, 1, vx - Cs, vx, vx, vx, vx + Cs);

		__global real* dF_dU = fluxXformBuf[intindex].A;
		fill(dF_dU+0,5,	0, 									1, 								0, 					0, 					0			);
		fill(dF_dU+1,5,	-vx * vx + .5 * gamma_1 * vSq,		-vx * gamma_3,					-vy * gamma_1,		-vz * gamma_1,		gamma - 1	);
		fill(dF_dU+2,5,	-vx * vy,							vy, 							vx, 				0, 					0			);
		fill(dF_dU+3,5,	-vx * vz, 							vz, 							0, 					vx, 				0			);
		fill(dF_dU+4,5,	vx * (.5 * vSq * gamma_1 - hTotal),	-gamma_1 * vx * vx + hTotal,	-gamma_1 * vx * vy,	-gamma_1 * vx * vz,	gamma * vx	);

		__global eigen_t* eigen = eigenBuf + intindex;
		
		__global real* evL = eigen->evL; 
		real invDenom = .5 / CsSq;
		fill(evL+0,5,(.5 * gamma_1 * vSq + Cs * vx) * invDenom,	-(Cs + gamma_1 * vx) * invDenom,	-gamma_1 * vy * invDenom,		-gamma_1 * vz * invDenom,		gamma_1 * invDenom		);
		fill(evL+1,5,1 - gamma_1 * vSq * invDenom,				gamma_1 * vx * 2 * invDenom,		gamma_1 * vy * 2 * invDenom,	gamma_1 * vz * 2 * invDenom,	-gamma_1 * 2 * invDenom	);
		fill(evL+2,5,-vy, 										0, 									1, 								0, 								0						);
		fill(evL+3,5,-vz, 										0, 									0, 								1, 								0						);
		fill(evL+4,5,(.5 * gamma_1 * vSq - Cs * vx) * invDenom,	(Cs - gamma_1 * vx) * invDenom,		-gamma_1 * vy * invDenom,		-gamma_1 * vz * invDenom,		gamma_1 * invDenom		);

		__global real* evR = eigen->evR;
		fill(evR+0,5,1, 				1, 			0, 	0, 	1				);
		fill(evR+1,5,vx - Cs, 			vx, 		0, 	0, 	vx + Cs			);
		fill(evR+2,5,vy, 				vy, 		1, 	0, 	vy				);
		fill(evR+3,5,vz, 				vz, 		0, 	1, 	vz				);
		fill(evR+4,5,hTotal - Cs * vx, .5 * vSq, 	vy, vz, hTotal + Cs * vx);

#if dim > 1
	if (side == 1) {
		for (int i = 0; i < numStates; ++i) {
			real tmp;
			//each row's xy <- yx
			tmp = evL[i + numStates * cons_mx];
			evL[i + numStates * cons_mx] = evL[i + numStates * cons_my];
			evL[i + numStates * cons_my] = tmp;
			//each column's xy <- yx
			tmp = evR[cons_mx + numStates * i];
			evR[cons_mx + numStates * i] = evR[cons_my + numStates * i];
			evR[cons_my + numStates * i] = tmp;
		}
	}
#endif
#if dim > 2
	else if (side == 2) {
		for (int i = 0; i < numStates; ++i) {
			real tmp;
			tmp = evL[i + numStates * cons_mx];
			evL[i + numStates * cons_mx] = evL[i + numStates * cons_mz];
			evL[i + numStates * cons_mz] = tmp;
			tmp = evR[cons_mx + numStates * i];
			evR[cons_mx + numStates * i] = evR[cons_mz + numStates * i];
			evR[cons_mz + numStates * i] = tmp;
		}
	}
#endif
	}
}
