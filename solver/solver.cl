// Roe solver:

kernel void calcLR(
	global consLR_t* ULRBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	
	//TODO skip this lr stuff if we're doing piecewise-constant
	//...and just use the original buffers
	<? for side=0,solver.dim-1 do ?>{
		int intindex = <?=side?> + dim * index;
		//constant
		ULRBuf[intindex].L = ULRBuf[intindex].R = UBuf[index];
	}<? end ?>
}

<? if solver.checkFluxError or solver.checkOrthoError then ?>
kernel void calcErrors(
	global error_t* errorBuf,
	const global real* waveBuf,
	const global eigen_t* eigenBuf
) {
	SETBOUNDS(0,0);

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int intindex = side + dim * index;
		const global real* wave = waveBuf + numWaves * intindex;
		const global eigen_t* eigen = eigenBuf + intindex;

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
			eigen_rightTransform_<?=side?>(eigenInvCoords, eigen, basis);
		
			real newbasis[numWaves];
			eigen_leftTransform_<?=side?>(newbasis, eigen, eigenInvCoords);
			
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
			eigen_leftTransform_<?=side?>(eigenCoords, eigen, basis);

			real eigenScaled[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				eigenScaled[j] = eigenCoords[j] * wave[j];
			}
			
			real newtransformed[numStates];
			eigen_rightTransform_<?=side?>(newtransformed, eigen, eigenScaled);
			
			real transformed[numStates];
			fluxTransform_<?=side?>(transformed, eigen, basis);
			
			for (int j = 0; j < numStates; ++j) {
				fluxError += fabs(newtransformed[j] - transformed[j]);
			}
		}
		
		errorBuf[intindex] = (error_t){
			.ortho = orthoError,
			.flux = fluxError,
		};
	}<? end ?>
}
<? end ?>

kernel void calcDeltaUTilde(
	global real* deltaUTildeBuf,
	const global consLR_t* ULRBuf,
	const global eigen_t* eigenBuf
) {
	SETBOUNDS(2,1);	
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		
		const global real* UL = (const global real*)&(ULRBuf + side + dim * indexL)->R;
		const global real* UR = (const global real*)&(ULRBuf + side + dim * indexR)->L;

		real deltaU[numStates];
		for (int j = 0; j < numStates; ++j) {
			deltaU[j] = UR[j] - UL[j];
		}
	
		int intindex = side + dim * index;	
		real deltaUTilde[numWaves];
		eigen_leftTransform_<?=side?>(
			deltaUTilde,
			eigenBuf + intindex,
			deltaU);
	
		//TODO memcpy
		global real* deltaUTilde_ = deltaUTildeBuf + intindex * numWaves;
		for (int j = 0; j < numWaves; ++j) {
			deltaUTilde_[j] = deltaUTilde[j];
		}
	}<? end ?>
}

kernel void calcRTilde(
	global real* rTildeBuf,
	const global real* deltaUTildeBuf,
	const global real* waveBuf
) {
	SETBOUNDS(2,1);
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
		int intindex = side + dim * index;
		int intindexL = side + dim * indexL;
		int intindexR = side + dim * indexR;
		global real* rTilde = rTildeBuf + intindex * numWaves;
		const global real* deltaUTilde = deltaUTildeBuf + intindex * numWaves;
		const global real* deltaUTildeL = deltaUTildeBuf + intindexL * numWaves;
		const global real* deltaUTildeR = deltaUTildeBuf + intindexR * numWaves;
		const global real* wave = waveBuf + intindex * numWaves;
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
	}<? end ?>
}

kernel void calcFlux(
	global real* fluxBuf,
	const global consLR_t* ULRBuf,
	const global real* waveBuf, 
	const global eigen_t* eigenBuf, 
	const global real* deltaUTildeBuf,
	const global real* rTildeBuf,	//not needed with slope limiters 
	real dt
) {
	SETBOUNDS(2,1);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / dx<?=side?>_at(i);
		int indexL = index - stepsize[side];
		const global real* UL = (const global real*)&(ULRBuf + side + dim * indexL)->R;
		const global real* UR = (const global real*)&(ULRBuf + side + dim * indexR)->L;
		
		real UAvg[numStates];
		for (int j = 0; j < numStates; ++j) {
			UAvg[j] = .5 * (UL[j] + UR[j]);
		}
		
		int intindex = side + dim * index;
		const global eigen_t* eigen = eigenBuf + intindex;

		real fluxTilde[numWaves];
		eigen_leftTransform_<?=side?>(fluxTilde, eigen, UAvg);

		const global real* lambdas = waveBuf + numWaves * intindex;
		const global real* deltaUTilde = deltaUTildeBuf + numWaves * intindex;
		const global real* rTilde = rTildeBuf + numWaves * intindex;

		for (int j = 0; j < numWaves; ++j) {
			real lambda = lambdas[j];
			fluxTilde[j] *= lambda;
			real theta = lambda >= 0 ? 1 : -1;
			real phi = fluxLimiter(rTilde[j]);
			real epsilon = lambda * dt_dx;
			real deltaFluxTilde = lambda * deltaUTilde[j];
			fluxTilde[j] -= .5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
		}

		real flux[numStates];
		eigen_rightTransform_<?=side?>(flux, eigen, fluxTilde);

		real3 interfaceI = _real3(i.x, i.y, i.z);
		interfaceI.s[side] -= .5;
		real3 interfaceX = cell_x(interfaceI);
		real volume = volume_at(interfaceX);

		global real* flux_ = fluxBuf + intindex * numStates;
		for (int j = 0; j < numStates; ++j) {
			flux_[j] = volume * flux[j];
		}
	}<? end ?>
}

kernel void calcDerivFromFlux(
	global real* derivBuf,
	const global real* fluxBuf
) {
	SETBOUNDS(2,2);
	global real* deriv = derivBuf + numStates * index;
		
	real volume = volume_at(cell_x(i));
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int intindexL = side + dim * index;
		int intindexR = intindexL + dim * stepsize[side]; 
		const global real* fluxL = fluxBuf + intindexL * numStates;
		const global real* fluxR = fluxBuf + intindexR * numStates;
		for (int j = 0; j < numStates; ++j) {
			real deltaFlux = fluxR[j] - fluxL[j];
			deriv[j] -= deltaFlux / (volume * grid_dx<?=side?>);
		}
	}<? end ?>
}
