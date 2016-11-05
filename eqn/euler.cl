//called from calcDT
range_t calcCellMinMaxEigenvalues(
	cons_t U,
	real ePot,
	int side
) {
	prim_t W = primFromCons(U, ePot);
	real Cs = calc_Cs(W);
	real v = W.v.s[side];
	return (range_t){.min = v - Cs, .max = v + Cs};
}

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
__kernel void calcDT(
	__global real* dtBuf,
	const __global cons_t* UBuf,
	const __global real* ePotBuf 
) {
	SETBOUNDS(0,0);
	if (OOB(2,2)) {
		dtBuf[index] = INFINITY;
		return;
	}
	
	cons_t U = UBuf[index];
	real ePot = ePotBuf[index];

	real dt = INFINITY;

	<? for side=0,solver.dim-1 do ?> {
		range_t lambda = calcCellMinMaxEigenvalues(U, ePot, <?=side?>);
		lambda.min = min((real)0., lambda.min);
		lambda.max = max((real)0., lambda.max);
		dt = min(dt, dx_at<?=side?>(i) / (fabs(lambda.max - lambda.min) + (real)1e-9));
	}<? end ?>
	dtBuf[index] = dt; 
}

typedef struct {
	real rho;
	real3 v;
	real hTotal;
} Roe_t;

Roe_t calcEigenBasisSide(
	cons_t UL,
	real ePotL, 
	cons_t UR,
	real ePotR
) {
	prim_t WL = primFromCons(UL, ePotL);
	real sqrtRhoL = sqrt(WL.rho);
	real3 vL = WL.v;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL.ETotal);

	prim_t WR = primFromCons(UR, ePotR);
	real sqrtRhoR = sqrt(UR.rho);
	real3 vR = WR.v;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR.ETotal);

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);

	return (Roe_t){
		.rho = sqrtRhoL * sqrtRhoR,
		.v = real3_add(
			real3_scale(vL, sqrtRhoL * invDenom),
			real3_scale(vR, sqrtRhoR * invDenom)),
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
	const __global cons_t* UBuf,	//[volume]
	const __global real* ePotBuf
#if defined(checkFluxError)
	, __global fluxXform_t* fluxXformBuf	//[volume][dim]
#endif
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
	
		cons_t UL = UBuf[indexL];
		cons_t UR = UBuf[indexR];

		real tmp;
		tmp = UL.m.s0; UL.m.s0 = UL.m.s[side]; UL.m.s[side] = tmp;
		tmp = UR.m.s0; UR.m.s0 = UR.m.s[side]; UR.m.s[side] = tmp;

		real ePotL = ePotBuf[indexL];
		real ePotR = ePotBuf[indexR];

		Roe_t roe = calcEigenBasisSide(UL, ePotL, UR, ePotR);
		real3 v = roe.v;
		real hTotal = roe.hTotal;
	
		real vSq = v.x*v.x + v.y*v.y + v.z*v.z;//coordLenSq(v);
		real CsSq = gamma_1 * (hTotal - .5 * vSq);
		real Cs = sqrt(CsSq);
	
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		fill(wave, 1, v.s0 - Cs, v.s0, v.s0, v.s0, v.s0 + Cs);

		__global eigen_t* eigen = eigenBuf + intindex;
		
		__global real* evL = eigen->evL; 
		real invDenom = .5 / CsSq;
		fill(evL+0,5,	(.5 * gamma_1 * vSq + Cs * v.s0) * invDenom,	-(Cs + gamma_1 * v.s0) * invDenom,	-gamma_1 * v.s1 * invDenom,		-gamma_1 * v.s2 * invDenom,		gamma_1 * invDenom		);
		fill(evL+1,5,	1 - gamma_1 * vSq * invDenom,					gamma_1 * v.s0 * 2 * invDenom,		gamma_1 * v.s1 * 2 * invDenom,	gamma_1 * v.s2 * 2 * invDenom,	-gamma_1 * 2 * invDenom	);
		fill(evL+2,5,	-v.s1, 											0, 									1, 								0, 								0						);
		fill(evL+3,5,	-v.s2, 											0, 									0, 								1, 								0						);
		fill(evL+4,5,	(.5 * gamma_1 * vSq - Cs * v.s0) * invDenom,	(Cs - gamma_1 * v.s0) * invDenom,	-gamma_1 * v.s1 * invDenom,		-gamma_1 * v.s2 * invDenom,		gamma_1 * invDenom		);

		__global real* evR = eigen->evR;
		fill(evR+0,5,	1, 					1, 			0, 		0, 		1					);
		fill(evR+1,5,	v.s0 - Cs, 			v.s0, 		0, 		0, 		v.s0 + Cs			);
		fill(evR+2,5,	v.s1, 				v.s1, 		1, 		0, 		v.s1				);
		fill(evR+3,5,	v.s2, 				v.s2, 		0, 		1, 		v.s2				);
		fill(evR+4,5,	hTotal - Cs * v.s0,	.5 * vSq, 	v.s1,	v.s2,	hTotal + Cs * v.s0	);

#if defined(checkFluxError)
		__global real* dF_dU = fluxXformBuf[intindex].A;
		fill(dF_dU+0,5,	0, 										1, 									0, 						0, 						0				);
		fill(dF_dU+1,5,	-v.s0 * v.s0 + .5 * gamma_1 * vSq,		-v.s0 * gamma_3,					-v.s1 * gamma_1,		-v.s2 * gamma_1,		gamma - 1		);
		fill(dF_dU+2,5,	-v.s0 * v.s1,							v.s1, 								v.s0, 					0, 						0				);
		fill(dF_dU+3,5,	-v.s0 * v.s2, 							v.s2, 								0, 						v.s0, 					0				);
		fill(dF_dU+4,5,	v.s0 * (.5 * vSq * gamma_1 - hTotal),	-gamma_1 * v.s0 * v.s0 + hTotal,	-gamma_1 * v.s0 * v.s1,	-gamma_1 * v.s0 * v.s2,	gamma * v.s0	);
#endif

#if dim > 1
	if (side == 1) {
		for (int i = 0; i < numStates; ++i) {
			real tmp;
			//each row's xy <- yx
			tmp = evL[i + numStates * cons_m0];
			evL[i + numStates * cons_m0] = evL[i + numStates * cons_m1];
			evL[i + numStates * cons_m1] = tmp;
			//each column's xy <- yx
			tmp = evR[cons_m0 + numStates * i];
			evR[cons_m0 + numStates * i] = evR[cons_m1 + numStates * i];
			evR[cons_m1 + numStates * i] = tmp;
#if defined(checkFluxError)
			tmp = dF_dU[i + numStates * cons_m0];
			dF_dU[i + numStates * cons_m0] = dF_dU[i + numStates * cons_m1];
			dF_dU[i + numStates * cons_m1] = tmp;
			tmp = dF_dU[cons_m0 + numStates * i];
			dF_dU[cons_m0 + numStates * i] = dF_dU[cons_m1 + numStates * i];
			dF_dU[cons_m1 + numStates * i] = tmp;
#endif
		}
	}
#endif	//dim > 1
#if dim > 2
	else if (side == 2) {
		for (int i = 0; i < numStates; ++i) {
			real tmp;
			tmp = evL[i + numStates * cons_m0];
			evL[i + numStates * cons_m0] = evL[i + numStates * cons_m2];
			evL[i + numStates * cons_m2] = tmp;
			tmp = evR[cons_m0 + numStates * i];
			evR[cons_m0 + numStates * i] = evR[cons_m2 + numStates * i];
			evR[cons_m2 + numStates * i] = tmp;
#if defined(checkFluxError)
			tmp = dF_dU[i + numStates * cons_m0];
			dF_dU[i + numStates * cons_m0] = dF_dU[i + numStates * cons_m2];
			dF_dU[i + numStates * cons_m2] = tmp;
			tmp = dF_dU[cons_m0 + numStates * i];
			dF_dU[cons_m0 + numStates * i] = dF_dU[cons_m2 + numStates * i];
			dF_dU[cons_m2 + numStates * i] = tmp;
#endif
		}
	}
#endif	//dim > 2
	}
}
