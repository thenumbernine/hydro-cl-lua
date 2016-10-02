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

//donor cell
//real slopeLimiter(real r) { return 0.; }
//Lax-Wendroff:
//real slopeLimiter(real r) { return 1.; }
//Superbee
real slopeLimiter(real r) { return max(0., max(min(1., 2. * r), min(2., r))); }

// private:

__kernel void convertToTex(
	__write_only dstimage_t tex,
	int displayVar,
	const __global real* buf
) {
	SETBOUNDS(0,0);
	real value = 0;
	int intindex = dim * index;	//side 0
	if (displayVar >= display_wave_0 && displayVar < display_wave_0 + numWaves) {
		const __global real* wave = buf + intindex * numWaves;
		value = wave[displayVar - display_wave_0];
	} else if (displayVar >= display_eigen_0 && displayVar < display_eigen_0 + numEigen) {
		const __global real* eigen = buf + intindex * numEigen;
		value = eigen[displayVar - display_eigen_0];
	} else if (displayVar >= display_deltaUTilde_0 && displayVar < display_deltaUTilde_0 + numWaves) {
		const __global real* deltaUTilde = buf + intindex * numWaves;
		value = deltaUTilde[displayVar - display_deltaUTilde_0];
	} else if (displayVar >= display_rTilde_0 && displayVar < display_rTilde_0 + numWaves) {
		const __global real* rTilde = buf + intindex * numWaves;
		value = rTilde[displayVar - display_rTilde_0];
	} else if (displayVar >= display_flux_0 && displayVar < display_flux_0 + numStates) {
		const __global real* flux = buf + intindex * numStates;
		value = flux[displayVar - display_flux_0];
	} else if (displayVar >= display_deriv_0 && displayVar < display_deriv_0 + numStates) {
		const __global real* deriv = buf + index * numStates;
		value = deriv[displayVar - display_deriv_0];
	} else if (displayVar == display_dt_0) {
		value = buf[index];
	} else if (displayVar == display_orthoError_0) {
		value = buf[intindex];
	} else {
		value = convertToTex_UBuf(displayVar, buf + numStates * index);
	}
	write_imagef(tex, WRITEIMAGEARGS, (float4)(value, 0., 0., 0.));
}

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
	
	range_t lambda = calcCellMinMaxEigenvalues(UBuf[index]); 
	lambda.min = min(0., lambda.min);
	lambda.max = max(0., lambda.max);
	dtBuf[index] = dx / (fabs(lambda.max - lambda.min) + 1e-9);
}

// the default eigen transforms, using eigen struct as a dense matrix:

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
	real* x
) {
	const __global real* A = eigen + numStates * numWaves;
	for (int i = 0; i < numWaves; ++i) {
		real sum = 0;
		for (int j = 0; j < numStates; ++j) {
			sum += A[i+numWaves*j] * x[j];
		}
		y[i] = sum;
	}
}

__kernel void calcErrors(
	__global real* orthoErrorBuf,
	__global real* fluxErrorBuf,
	const __global real* waveBuf,
	const __global real* eigenBuf,
	const __global real* fluxMatrixBuf
) {
	SETBOUNDS(2,1);
	for (int side = 0; side < dim; ++side) {
		int intindex = side + dim * index;
		const __global real* wave = waveBuf + numWaves * intindex;
		const __global real* eigen = eigenBuf + numEigen * intindex;
		const __global real* fluxMatrix = fluxMatrixBuf + numStates * numStates * intindex;

		real orthoError = 0;
		real fluxError = 0;
		for (int k = 0; k < numStates; ++k) {
			
			real src[numStates];
			for (int j = 0; j < numStates; ++j) {
				src[j] = k == j ? 1 : 0;
			}
			
			real mid[numWaves];
			eigenLeftTransform(mid, eigen, src);
			
			real scaled[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				scaled[j] = mid[j] * wave[j];
			}
			
			real identCheck[numStates];
			eigenRightTransform(identCheck, eigen, mid);
			
			real fluxCheck[numStates];
			eigenRightTransform(fluxCheck, eigen, scaled);
			
			for (int j = 0; j < numStates; ++j) {
				orthoError += fabs(identCheck[j] - src[j]);
				fluxError += fabs(fluxCheck[j] - fluxMatrix[j + numStates * k]);
			}
		}
		orthoErrorBuf[intindex] = log(orthoError) / log(10.);
		fluxErrorBuf[intindex] = log(fluxError) / log(10.);
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
		int indexL = index - stepsize[side];
	
		const __global real* UL = UBuf + indexL * numStates;
		const __global real* UR = UBuf + indexR * numStates;
	
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
		const __global real* UL = UBuf + indexL * numStates;
		const __global real* UR = UBuf + indexR * numStates;
		
		real UAvg[numStates];
		for (int j = 0; j < numStates; ++j) {
			UAvg[j] = .5 * (UL[j] + UR[j]);
		}
		
		int intindex = side + dim * index;
		const __global real* eigen = eigenBuf + intindex * numEigen;

		real fluxTilde[numWaves];
		eigenLeftTransform(fluxTilde, eigen, UAvg);

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
		eigenRightTransform(flux, eigen, fluxTilde);

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

__kernel void multAdd(
	__global real* a,
	const __global real* b,
	const __global real* c,
	real d
) {
	size_t i = get_global_id(0);
	if (i >= get_global_size(0)) return;
	a[i] = b[i] + c[i] * d;
}
