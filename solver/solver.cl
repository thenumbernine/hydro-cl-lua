// Roe solver:

#if defined(checkFluxError) || defined(checkOrthoError)
__kernel void calcErrors(
	__global error_t* errorBuf,
	const __global real* waveBuf,
	const __global eigen_t* eigenBuf,
	const __global fluxXform_t* fluxXformBuf
) {
	SETBOUNDS(0,0);

	for (int side = 0; side < dim; ++side) {
		int intindex = side + dim * index;
		const __global real* wave = waveBuf + numWaves * intindex;
		const __global eigen_t* eigen = eigenBuf + intindex;
		const __global fluxXform_t* fluxXform = fluxXformBuf + intindex;

		real orthoError = 0;
		real fluxError = 0;
		
		//the flux transform is F v = R Lambda L v, I = R L
		//but if numWaves < numStates then certain v will map to the nullspace 
		//so to test orthogonality for only numWaves dimensions, I will verify that Qinv Q v = v 
		//I = L R
		for (int k = 0; k < numWaves; ++k) {
			real basis[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				basis[j] = k == j ? 1 : 0;
			}
			
			real eigenInvCoords[numStates];
			eigen_rightTransform(eigenInvCoords, eigen, basis, side);
		
			real newbasis[numWaves];
			eigen_leftTransform(newbasis, eigen, eigenInvCoords, side);
			
			for (int j = 0; j < numWaves; ++j) {
				orthoError += fabs(newbasis[j] - basis[j]);
			}
		}
		
		for (int k = 0; k < numStates; ++k) {
			real basis[numStates];
			for (int j = 0; j < numStates; ++j) {
				basis[j] = k == j ? 1 : 0;
			}

			real eigenCoords[numWaves];
			eigen_leftTransform(eigenCoords, eigen, basis, side);

			real eigenScaled[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				eigenScaled[j] = eigenCoords[j] * wave[j];
			}
			
			real newtransformed[numStates];
			eigen_rightTransform(newtransformed, eigen, eigenScaled, side);
			
			real transformed[numStates];
			fluxTransform(transformed, fluxXform, basis, side);
			
			for (int j = 0; j < numStates; ++j) {
				fluxError += fabs(newtransformed[j] - transformed[j]);
			}
		}
		
		errorBuf[intindex] = (error_t){
			.ortho = orthoError,
			.flux = fluxError,
		};
	}
}
#endif	//defined(checkFluxError) || defined(checkOrthoError)

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
			deltaU,
			side);
	
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
		real dt_dx = dt / dx_at(i, side);
		
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
		eigen_leftTransform(fluxTilde, eigen, UAvg, side);

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
		eigen_rightTransform(flux, eigen, fluxTilde, side);

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

	//would it be faster to fill the buffer with zeros beforehand?
	__global real* deriv = derivBuf + numStates * index;
	for (int j = 0; j < numStates; ++j) {
		deriv[j] = 0;
	}
	
	if (OOB(2,2)) return;
	
	for (int side = 0; side < dim; ++side) {
		int intindexL = side + dim * index;
		int intindexR = intindexL + dim * stepsize[side]; 
		const __global real* fluxL = fluxBuf + numStates * intindexL;
		const __global real* fluxR = fluxBuf + numStates * intindexR;
		for (int j = 0; j < numStates; ++j) {
			real deltaFlux = fluxR[j] - fluxL[j];
			deriv[j] -= deltaFlux / dx_at(i,side);
		}
	}
}
