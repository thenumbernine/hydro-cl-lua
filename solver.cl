// common types:

typedef struct {
	real min, max;
} range_t;

// common:

//http://developer.amd.com/resources/documentation-articles/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
//calculate min of all elements on buffer[0..length-1]
__kernel void reduceMin(
	const __global real* buffer,
	__local real* scratch,
	__const int length,
	__global real* result)
{
	int global_index = get_global_id(0);
	real accumulator = INFINITY;
	
	// Loop sequentially over chunks of input vector
	while (global_index < length) {
		real element = buffer[global_index];
		accumulator = (accumulator < element) ? accumulator : element;
		global_index += get_global_size(0);
	}

	// Perform parallel reduction
	int local_index = get_local_id(0);
	scratch[local_index] = accumulator;
	barrier(CLK_LOCAL_MEM_FENCE);
	for (int offset = get_local_size(0) / 2; offset > 0; offset = offset / 2) {
		if (local_index < offset) {
			real other = scratch[local_index + offset];
			real mine = scratch[local_index];
			scratch[local_index] = (mine < other) ? mine : other;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (local_index == 0) {
		result[get_group_id(0)] = scratch[0];
	}
}

// private:

constant real gamma_1 = gamma - 1;
constant real gamma_3 = gamma - 3;
real calc_hTotal(real rho, real P, real ETotal) {
	return (P + ETotal) / rho;
}

// equation/solver API: 

constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);
constant real4 dxs = (real4)(dx, 0, 0, 0);

real slopeLimiter(real r) { return 0; }
	
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

__kernel void convertToTex(
	__write_only dstimage_t tex,
	const __global cons_t* UBuf,
	int displayVar
) {
	SETBOUNDS(0,0);
	cons_t U = UBuf[index];
	prim_t W = primFromCons(U);
	real rho = W.rho;
	real vx = W.vx;
	real P = W.P;
	
	real value = 0;
	switch (displayVar) {
	case display_rho: value = rho; break;
	case display_vx: value = vx; break;
	case display_P: value = P; break;
	case display_eInt: value = P / (rho * gamma_1); break;
	case display_eKin: value = .5 * vx * vx; break;
	case display_eTotal: value = U.ETotal / rho; break;
	}
	write_imagef(tex, WRITEIMAGEARGS, (float4)(value, 0., 0., 0.));
}

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

// Roe solver:

__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real x = (real)(i.x + .5) / (real)gridSize_x * (xmax - xmin) + xmin;
	UBuf[index] = consFromPrim((prim_t){
		.rho = x < 0 ? 1 : .125,
		.vx = 0,
		.P = x < 0 ? 1 : .1,
	}); 
}

__kernel void calcDT(
	__global real* dtBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(2,2);
	range_t lambda = calcCellMinMaxEigenvalues(UBuf[index]); 
	lambda.min = min(0., lambda.min);
	lambda.max = max(0., lambda.max);
	dtBuf[index] = dx / (fabs(lambda.max - lambda.min) + 1e-9);
}

void fill3(__global real* ptr, real a, real b, real c) {
	ptr[0] = a;
	ptr[1] = b;
	ptr[2] = c;
}
	
__kernel void calcEigenBasis(
	__global real* waveBuf,	// wave buffer, size [volume][dim][numWaves]
	__global real* eigenBuf,	// eigen buffer, size [volume][dim][numEigen]
	const __global cons_t *UBuf
) {
	SETBOUNDS(2,1);
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int4 iL = i;
		iL[side] = (iL[side]-1+gridSize[side]) % gridSize[side];
		int indexL = INDEXV(iL);
		
		Roe_t roe = calcEigenBasisSide(UBuf[indexL], UBuf[indexR]);
		real vx = roe.vx;
		real hTotal = roe.hTotal;
		
		real vxSq = vx * vx;	
		real CsSq = gamma_1 * (hTotal - .5 * vx * vx);
		real Cs = sqrt(CsSq);
	
		int intindex = side + dim * index;	
		__global real* wave = waveBuf + numWaves * intindex;
		fill3(wave, vx - Cs, vx, vx + Cs);

		/*
		fill(dF_dU[1], 0, 									1, 							0			)
		fill(dF_dU[2], .5 * gamma_3 * vxSq, 				-gamma_3 * vx, 				gamma_1		)
		fill(dF_dU[3], vx * (.5 * gamma_1 * vxSq - hTotal), hTotal - gamma_1 * vxSq,	gamma*vx	)
		*/

		__global real* evL = eigenBuf + numEigen * intindex;	
		__global real* evR = evL + 3*3;

		fill3(evL+0, (.5 * gamma_1 * vxSq + Cs * vx) / (2. * CsSq),	-(Cs + gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);
		fill3(evL+3, 1. - gamma_1 * vxSq / (2. * CsSq),				gamma_1 * vx / CsSq,				-gamma_1 / CsSq			);
		fill3(evL+6, (.5 * gamma_1 * vxSq - Cs * vx) / (2. * CsSq),	(Cs - gamma_1 * vx) / (2. * CsSq),	gamma_1 / (2. * CsSq)	);

		fill3(evR+0, 1., 				1., 		1.				);
		fill3(evR+3, vx - Cs, 			vx, 		vx + Cs			);
		fill3(evR+6, hTotal - Cs * vx, .5 * vxSq, 	hTotal + Cs * vx);
	}
}

void eigenLeftTransform(
	real* y,
	const __global real* eigen,
	real* x)
{
	const __global real* A = eigen;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

void eigenRightTransform(
	real* y,
	const __global real* eigen,
	real* x)
{
	const __global real* A = eigen + numStates * numWaves;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}
__kernel void calcDeltaUTilde(
	__global real* deltaUTildeBuf,
	const __global real* UBuf,
	const __global real* eigenBuf
) {
	SETBOUNDS(2,1);	
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int4 iL = i;
		iL[side] = (iL[side]-1+gridSize[side]) % gridSize[side];
		int indexL = INDEXV(iL);
	
		const __global real* UL = (const __global real*)(UBuf + indexL);
		const __global real* UR = (const __global real*)(UBuf + indexR);
	
		real deltaU[numStates];
		for (int j = 0; j < numStates; ++j) {
			deltaU[j] = UR[j] - UL[j];
		}
	
		int intindex = side + dim * index;	
		real deltaUTilde[numWaves];
		eigenLeftTransform(
			deltaUTilde,
			eigenBuf + intindex * numEigen,
			deltaU);
	
		__global real* deltaUTilde_ = deltaUTildeBuf + intindex * numWaves;
		//TODO memcpy
		for (int j = 0; j < numWaves; ++j) {
			deltaUTilde_[j] = deltaUTilde[j];
		}
	}
}

constant int4 stepsize = (int4)(1, 
	gridSize_x, 
	gridSize_x * gridSize_y, 
	gridSize_x * gridSize_y * gridSize_z);

__kernel void calcRTilde(
	__global real* rTildeBuf,
	const __global real* deltaUTildeBuf,
	const __global real* waveBuf
) {
	SETBOUNDS(2,1);
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
		int intindex = side + dim * index;
		int intindexL = side + dim * indexL;
		int intindexR = side + dim * indexR;
		__global real* rTilde = rTildeBuf + intindex * numWaves;
		const __global real* deltaUTilde = deltaUTildeBuf + intindex * numWaves;
		const __global real* deltaUTildeL = deltaUTildeBuf + intindexL * numWaves;
		const __global real* deltaUTildeR = deltaUTildeBuf + intindexR * numWaves;
		const __global real* wave = waveBuf + intindex * numWaves;
		for (int j = 0; j < numWaves; ++j) {
			if (deltaUTilde[j] == 0) {
				rTilde[j] = 0;
			} else {
				if (wave[j] >= 0) {
					rTilde[j] = deltaUTildeR[j] / deltaUTilde[j];
				} else {
					rTilde[j] = deltaUTildeL[j] / deltaUTilde[j];
				}
			}
		}
	}
}

__kernel void calcFlux(
	__global cons_t* fluxBuf,
	const __global cons_t* UBuf,
	const __global real* waveBuf, 
	const __global real* eigenBuf, 
	const __global real* deltaUTildeBuf,
	const __global real* rTildeBuf,
	real dt
) {
	SETBOUNDS(2,1);
	for (int side = 0; side < dim; ++side) {
		real dt_dx = dt / dxs[side];
		
		int indexL = index - stepsize[side];
		int indexR = index;
		const __global real* UL = (const __global real*)(UBuf + indexL);
		const __global real* UR = (const __global real*)(UBuf + indexR);
		
		real UAvg[numStates];
		for (int j = 0; j < numStates; ++j) {
			UAvg[j] = .5 * (UL[j] + UR[j]);
		}
		
		int intindex = side + dim * index;
		const __global real* eigen = eigenBuf + intindex * numEigen;

		real fluxTilde[numWaves];
		eigenLeftTransform(fluxTilde, eigen, UAvg);

		const __global real* lambdas = waveBuf + intindex * numWaves;
		const __global real* deltaUTilde = deltaUTildeBuf + intindex * numWaves;

		for (int j = 0; j < numWaves; ++j) {
			real lambda = lambdas[j];
			fluxTilde[j] *= lambda;
			real rTilde = rTildeBuf[j + numWaves * intindex];
			real theta = lambda >= 0 ? 1 : -1;
			real phi = slopeLimiter(rTilde);
			real epsilon = lambda * dt_dx;
			real deltaFluxTilde = lambda * deltaUTilde[j];
			fluxTilde[j] -= .5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
		}

		real flux[numStates];
		eigenRightTransform(flux, eigen, fluxTilde);

		__global real* flux_ = (__global real*)(fluxBuf + intindex);
		for (int j = 0; j < numStates; ++j) {
			flux_[j] = flux[j];
		}
	}
}

__kernel void calcDerivFromFlux(
	__global real* derivBuf,
	const __global real* fluxBuf
) {
	SETBOUNDS(2,2);
	__global real* deriv = derivBuf + numStates * index;	
	for (int side = 0; side < dim; ++side) {
		int intindexL = side + dim * index;
		int intindexR = intindexL + dim * stepsize[side]; 
		const __global real* fluxL = fluxBuf + numStates * intindexL;
		const __global real* fluxR = fluxBuf + numStates * intindexR;
		for (int j = 0; j < numStates; ++j) {
			real deltaFlux = fluxR[j] - fluxL[j];
			deriv[j] -= deltaFlux / dxs[side];
		}
	}
}

__kernel void multAddTo(
	__global real* a,
	const __global real* b,
	real c
) {
	size_t i = get_global_id(0);
	if (i >= get_global_size(0)) return;
	a[i] += b[i] * c;
}
