// Roe solver:

__kernel void calcDT(
	__global real* dtBuf,
	const __global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	if (i.x < 2 || i.x >= gridSize_x - 2 
#if dim > 1
		|| i.y < 2 || i.y >= gridSize_y - 2
#endif
#if dim > 2
		|| i.z < 2 || i.z >= gridSize_z - 2
#endif
	) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const __global cons_t* U = UBuf + index;

	real dt = INFINITY;
	for (int side = 0; side < dim; ++side) {
		range_t lambda = calcCellMinMaxEigenvalues(U, side); 
		lambda.min = min((real)0., lambda.min);
		lambda.max = max((real)0., lambda.max);
		dt = min(dt, dxs[side] / (fabs(lambda.max - lambda.min) + (real)1e-9));
	}
	dtBuf[index] = dt; 
}

//1/log(10)
#define _1_LN_10 0.4342944819032517611567811854911269620060920715332

__kernel void calcErrors(
	__global error_t* errorBuf,
	const __global real* waveBuf,
	const __global eigen_t* eigenBuf,
	const __global real* fluxMatrixBuf
) {
	SETBOUNDS(2,1);
	for (int side = 0; side < dim; ++side) {
		int intindex = side + dim * index;
		const __global real* wave = waveBuf + numWaves * intindex;
		const __global eigen_t* eigen = eigenBuf + intindex;
		const __global real* fluxMatrix = fluxMatrixBuf + numStates * numStates * intindex;

		real orthoError = 0;
		real fluxError = 0;
		for (int k = 0; k < numStates; ++k) {
			
			real src[numStates];
			for (int j = 0; j < numStates; ++j) {
				src[j] = k == j ? 1 : 0;
			}
			
			real mid[numWaves];
			eigen_leftTransform(mid, eigen, src);
			
			real scaled[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				scaled[j] = mid[j] * wave[j];
			}
			
			real identCheck[numStates];
			eigen_rightTransform(identCheck, eigen, mid);
			
			real fluxCheck[numStates];
			eigen_rightTransform(fluxCheck, eigen, scaled);
			
			for (int j = 0; j < numStates; ++j) {
				orthoError += fabs(identCheck[j] - src[j]);
				fluxError += fabs(fluxCheck[j] - fluxMatrix[j + numStates * k]);
			}
		}
		errorBuf[intindex] = (error_t){
			.ortho = log(orthoError) * _1_LN_10,
			.flux = log(fluxError) * _1_LN_10,
		};
	}
}

__kernel void calcDeltaUTilde(
	__global real* deltaUTildeBuf,
	const __global real* UBuf,
	const __global eigen_t* eigenBuf
) {
	SETBOUNDS(2,1);	
	int indexR = index;
	for (int side = 0; side < dim; ++side) {
		int indexL = index - stepsize[side];
	
		const __global real* UL = UBuf + indexL * numStates;
		const __global real* UR = UBuf + indexR * numStates;
	
		real deltaU[numStates];
		for (int j = 0; j < numStates; ++j) {
			deltaU[j] = UR[j] - UL[j];
		}
	
		int intindex = side + dim * index;	
		real deltaUTilde[numWaves];
		eigen_leftTransform(
			deltaUTilde,
			eigenBuf + intindex,
			deltaU);
	
		//TODO memcpy
		__global real* deltaUTilde_ = deltaUTildeBuf + intindex * numWaves;
		for (int j = 0; j < numWaves; ++j) {
			deltaUTilde_[j] = deltaUTilde[j];
		}
	}
}

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
					rTilde[j] = deltaUTildeL[j] / deltaUTilde[j];
				} else {
					rTilde[j] = deltaUTildeR[j] / deltaUTilde[j];
				}
			}
		}
	}
}

__kernel void calcFlux(
	__global real* fluxBuf,
	const __global real* UBuf,
	const __global real* waveBuf, 
	const __global eigen_t* eigenBuf, 
	const __global real* deltaUTildeBuf,
	const __global real* rTildeBuf,
	real dt
) {
	SETBOUNDS(2,1);
	for (int side = 0; side < dim; ++side) {
		real dt_dx = dt / dxs[side];
		
		int indexL = index - stepsize[side];
		int indexR = index;
		const __global real* UL = UBuf + indexL * numStates;
		const __global real* UR = UBuf + indexR * numStates;
		
		real UAvg[numStates];
		for (int j = 0; j < numStates; ++j) {
			UAvg[j] = .5 * (UL[j] + UR[j]);
		}
		
		int intindex = side + dim * index;
		const __global eigen_t* eigen = eigenBuf + intindex;

		real fluxTilde[numWaves];
		eigen_leftTransform(fluxTilde, eigen, UAvg);

		const __global real* lambdas = waveBuf + numWaves * intindex;
		const __global real* deltaUTilde = deltaUTildeBuf + numWaves * intindex;
		const __global real* rTilde = rTildeBuf + numWaves * intindex;

		for (int j = 0; j < numWaves; ++j) {
			real lambda = lambdas[j];
			fluxTilde[j] *= lambda;
			real theta = lambda >= 0 ? 1 : -1;
			real phi = slopeLimiter(rTilde[j]);
			real epsilon = lambda * dt_dx;
			real deltaFluxTilde = lambda * deltaUTilde[j];
			fluxTilde[j] -= .5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
		}

		real flux[numStates];
		eigen_rightTransform(flux, eigen, fluxTilde);

		__global real* flux_ = fluxBuf + intindex * numStates;
		for (int j = 0; j < numStates; ++j) {
			flux_[j] = flux[j];
		}
	}
}

__kernel void calcDerivFromFlux(
	__global real* derivBuf,
	const __global real* fluxBuf
) {
	SETBOUNDS(0,0);
	
	__global real* deriv = derivBuf + numStates * index;
	for (int j = 0; j < numStates; ++j) {
		deriv[j] = 0;
	}
	
	if (i.x < 2 || i.x >= gridSize_x - 2 
#if dim > 1
		|| i.y < 2 || i.y >= gridSize_y - 2
#endif
#if dim > 2
		|| i.z < 2 || i.z >= gridSize_z - 2
#endif
	) {
		return;
	}
	
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
