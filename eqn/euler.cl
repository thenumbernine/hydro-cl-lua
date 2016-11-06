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

	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
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
	
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
	
		cons_t UL = UBuf[indexL];
		cons_t UR = UBuf[indexR];

		//TODO real3 n = normalForSide(side);
#if 1
		//normal
		real3 n = _real3(0,0,0);
		n.s[side] = 1;

		//tangent space
		real3 n1 = _real3(0,0,0);
		n1.s[(side+1)%3] = 1;

		real3 n2 = _real3(0,0,0);
		n2.s[(side+2)%3] = 1;
#else
		real3 n = e<?=side?>unit_at(i);
		real3 n1 = e<?=(side+1)%3?>unit_at(i);
		real3 n2 = e<?=(side+2)%3?>unit_at(i);
#endif

		real ePotL = ePotBuf[indexL];
		real ePotR = ePotBuf[indexR];

		Roe_t roe = calcEigenBasisSide(UL, ePotL, UR, ePotR);
		real3 v = roe.v;
		real hTotal = roe.hTotal;

		real v_n = real3_dot(v,n);
		real v_n1 = real3_dot(v,n1);
		real v_n2 = real3_dot(v,n2);
		real vSq = real3_dot(v,v);
		real eKin = .5 * vSq;
		real CsSq = gamma_1 * (hTotal - eKin);
		real Cs = sqrt(CsSq);
	
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		fill(wave, 1, v_n - Cs, v_n, v_n, v_n, v_n + Cs);

		__global eigen_t* eigen = eigenBuf + intindex;

		__global real* evR = eigen->evR;
		fill(evR+0,5,	1,					1,		0,		0,		1					);
		fill(evR+1,5,	v.x - n.x * Cs,		v.x,	n1.x,	n2.x,	v.x + n.x * Cs		);
		fill(evR+2,5,	v.y - n.y * Cs,		v.y,	n1.y,	n2.y,	v.y + n.y * Cs		);
		fill(evR+3,5,	v.z - n.z * Cs,		v.z,	n1.z,	n2.z,	v.z + n.z * Cs		);
		fill(evR+4,5,	hTotal - v_n * Cs,	eKin,	v_n1,	v_n2,	hTotal + v_n * Cs	);

		__global real* evL = eigen->evL; 
		real invDenom = .5 / CsSq;
		fill(evL+0,5,	(gamma_1 * eKin + Cs * v_n) * invDenom,	-(n.x * Cs + gamma_1 * v.x) * invDenom,	-(n.y * Cs + gamma_1 * v.y) * invDenom,	-(n.z * Cs + gamma_1 * v.z) * invDenom,	gamma_1 * invDenom			);
		fill(evL+1,5,	1 - gamma_1 * vSq * invDenom,			gamma_1 * v.x * 2 * invDenom,			gamma_1 * v.y * 2 * invDenom,			gamma_1 * v.z * 2 * invDenom,			-gamma_1 * 2 * invDenom		);
		fill(evL+2,5,	-v_n1,									n1.x,									n1.y,									n1.z,									0							);
		fill(evL+3,5,	-v_n2,									n2.x,									n2.y,									n2.z,									0							);
		fill(evL+4,5,	(gamma_1 * eKin - Cs * v_n) * invDenom,	(n.x * Cs - gamma_1 * v.x) * invDenom,	(n.y * Cs - gamma_1 * v.y) * invDenom,	(n.z * Cs - gamma_1 * v.z) * invDenom,	gamma_1 * invDenom			);

#if defined(checkFluxError)
		__global real* dF_dU = fluxXformBuf[intindex].A;
		fill(dF_dU+0,5,	0,									n.x,									n.y,									n.z,									0				);
		fill(dF_dU+1,5,	-v_n * v.x + gamma_1 * eKin * n.x,	v.x * n.x - gamma_1 * n.x * v.x + v_n,	v.x * n.y - gamma_1 * n.x * v.y,		v.x * n.z - gamma_1 * n.x * v.z,		gamma_1 * n.x	);
		fill(dF_dU+2,5,	-v_n * v.y + gamma_1 * eKin * n.y,	v.y * n.x - gamma_1 * n.y * v.x,		v.y * n.y - gamma_1 * n.y * v.y + v_n,	v.y * n.z - gamma_1 * n.y * v.z,		gamma_1 * n.y	);
		fill(dF_dU+3,5,	-v_n * v.z + gamma_1 * eKin * n.z,	v.z * n.x - gamma_1 * n.z * v.x,		v.z * n.y - gamma_1 * n.z * v.y,		v.z * n.z - gamma_1 * n.z * v.z + v_n,	gamma_1 * n.z	);
		fill(dF_dU+4,5,	v_n * (gamma_1 * eKin - hTotal),	-gamma_1 * v_n * v.x + n.x * hTotal,	-gamma_1 * v_n * v.y + n.y * hTotal,	-gamma_1 * v_n * v.z + n.z * hTotal,	gamma * v_n		);
#endif
	}<? end ?>
}
