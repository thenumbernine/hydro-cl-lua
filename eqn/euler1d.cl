//called from calcDT
range_t calcCellMinMaxEigenvalues(
	const __global cons_t* U,
	int side
) {
	prim_t W = primFromCons(*U);
	real Cs = sqrt(heatCapacityRatio * W.P / W.rho);
	return (range_t){.min=W.vx - Cs, .max=W.vx + Cs};
}

inline real calc_hTotal(real rho, real P, real ETotal) {
	return (P + ETotal) / rho;
}

__kernel void calcEigenBasis(
	__global real* waveBuf,			//[volume][dim][numWaves]
	__global eigen_t* eigenBuf,		//[volume][dim]
	const __global cons_t *UBuf		//[volume]
) {
	SETBOUNDS(2,1);
	int indexR = index;
	cons_t UR = UBuf[indexR];
	prim_t WR = primFromCons(UR);
	real sqrtRhoR = sqrt(UR.rho);
	real vxR = WR.vx;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR.ETotal);

	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];

		cons_t UL = UBuf[indexL];
		prim_t WL = primFromCons(UL);
		real sqrtRhoL = sqrt(WL.rho);
		real vxL = WL.vx;
		real hTotalL = calc_hTotal(WL.rho, WL.P, UL.ETotal);

		real invDenom = 1./(sqrtRhoL + sqrtRhoR);

		real rho = sqrtRhoL * sqrtRhoR;
		real vx = invDenom * (sqrtRhoL * vxL + sqrtRhoR * vxR);
		real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);
		
		real vxSq = vx * vx;	
		real CsSq = (heatCapacityRatio - 1.) * (hTotal - .5 * vx * vx);
		real Cs = sqrt(CsSq);
	
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		wave[0] = vx - Cs;
		wave[1] = vx;
		wave[2] = vx + Cs;

		__global eigen_t* eigen = eigenBuf + intindex;
		__global real* evR = eigen->evR;
		evR[0+3*0] = 1.;
		evR[1+3*0] = vx - Cs;
		evR[2+3*0] = hTotal - Cs * vx;
		evR[0+3*1] = 1.;
		evR[1+3*1] = vx;
		evR[2+3*1] = .5 * vxSq;
		evR[0+3*2] = 1.;
		evR[1+3*2] = vx + Cs;
		evR[2+3*2] = hTotal + Cs * vx;
		__global real* evL = eigen->evL; 
		evL[0+3*0] = (.5 * (heatCapacityRatio - 1.) * vxSq + Cs * vx) / (2. * CsSq);
		evL[0+3*1] = -(Cs + (heatCapacityRatio - 1.) * vx) / (2. * CsSq);
		evL[0+3*2] = (heatCapacityRatio - 1.) / (2. * CsSq);
		evL[1+3*0] = 1. - (heatCapacityRatio - 1.) * vxSq / (2. * CsSq);
		evL[1+3*1] = (heatCapacityRatio - 1.) * vx / CsSq;
		evL[1+3*2] = -(heatCapacityRatio - 1.) / CsSq;
		evL[2+3*0] = (.5 * (heatCapacityRatio - 1.) * vxSq - Cs * vx) / (2. * CsSq);
		evL[2+3*1] = (Cs - (heatCapacityRatio - 1.) * vx) / (2. * CsSq);
		evL[2+3*2] = (heatCapacityRatio - 1.) / (2. * CsSq);
<? if solver.checkFluxError then ?>
		__global real* A = eigen->A;
		A[0+3*0] = 0;
		A[0+3*1] = 1;
		A[0+3*2] = 0;
		A[1+3*0] = .5 * (heatCapacityRatio-3.) * vxSq;
		A[1+3*1] = -(heatCapacityRatio-3.) * vx;
		A[1+3*2] = (heatCapacityRatio - 1.);
		A[2+3*0] = vx * (.5 * (heatCapacityRatio - 1.) * vxSq - hTotal);
		A[2+3*1] = hTotal - (heatCapacityRatio - 1.) * vxSq;
		A[2+3*2] = heatCapacityRatio*vx;
<? end ?>
	}
}
