#define gamma_1 (gamma-1.)
#define gamma_3 (gamma-3.)

real calc_hTotal(real rho, real P, real ETotal) { return (P + ETotal) / rho; }
real calc_eKin(prim_t W) { return .5 * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz); }
real calc_EKin(prim_t W) { return W.rho * calc_eKin(W); }
real calc_EInt(prim_t W) { return W.P / gamma_1; }
real calc_eInt(prim_t W) { return calc_EInt(W) / W.rho; }
real calc_ETotal(prim_t W) { return calc_EKin(W) + calc_EInt(W); }

cons_t consFromPrim(prim_t W) {
	return (cons_t){
		.rho = W.rho,
		.mx = W.rho * W.vx,
		.my = W.rho * W.vy,
		.mz = W.rho * W.vz,
		.ETotal = calc_ETotal(W),
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
	real vy = 0;
	real vz = 0;
	real P = 0;
	
	INIT_STATE_CODE
	
	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .vy=vy, .vz=vz, .P=P});
}

prim_t primFromCons(cons_t U) {
	real EInt = U.ETotal - .5 * (U.mx * U.mx + U.my * U.my + U.mz * U.mz) / U.rho;
	return (prim_t){
		.rho = U.rho,
		.vx = U.mx / U.rho,
		.vy = U.my / U.rho,
		.vz = U.mz / U.rho,
		.P = EInt / gamma_1,
	};
}

//called from calcDisplayVar
real calcDisplayVar_UBuf(int displayVar, const __global real* U_) {
	const __global cons_t* U = (const __global cons_t*)U_;
	prim_t W = primFromCons(*U);
	switch (displayVar) {
	case display_U_rho: return W.rho;
	case display_U_vx: return W.vx;
	case display_U_vy: return W.vy;
	case display_U_vz: return W.vz;
	case display_U_v: return sqrt(W.vx * W.vx + W.vy * W.vy + W.vz * W.vz);
	case display_U_mx: return U->mx;
	case display_U_my: return U->my;
	case display_U_mz: return U->mz;
	case display_U_m: return sqrt(U->mx * U->mx + U->my * U->my + U->mz * U->mz);
	case display_U_P: return W.P;
	case display_U_eInt: return calc_eInt(W);
	case display_U_eKin: return calc_eKin(W);
	case display_U_eTotal: return U->ETotal / W.rho;
	case display_U_EInt: return calc_EInt(W);
	case display_U_EKin: return calc_EKin(W);
	case display_U_ETotal: return U->ETotal;
	case display_U_S: return W.P / pow(W.rho, (real)gamma);
	case display_U_H: return W.P * gamma / gamma_1;
	case display_U_h: return W.P * gamma / gamma_1 / W.rho;
	case display_U_HTotal: return W.P * gamma / gamma_1 + .5 * W.rho * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz);
	case display_U_hTotal: return W.P * gamma / gamma_1 / W.rho + .5 * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz);
	}
	return 0;
}

//called from calcDT
range_t calcCellMinMaxEigenvalues(const __global cons_t* U, int side) {
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
	__global real* fluxMatrixBuf,	//[volume][dim][numStates][numStates]
	const __global cons_t *UBuf		//[volume]
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

		__global real* dF_dU = fluxMatrixBuf + numStates * numStates * intindex;
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
