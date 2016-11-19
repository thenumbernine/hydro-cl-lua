// Roe solver:

kernel void calcLR(
	global consLR_t* ULRBuf,
	const global cons_t* UBuf,
	real dt
) {
	SETBOUNDS(1,1);
	const global cons_t* U = UBuf + index;
	
	//TODO skip this lr stuff if we're doing piecewise-constant
	//...and just use the original buffers
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int intindex = side + dim * index;
		global consLR_t* ULR = ULRBuf + intindex;	
		
<? if not solver.usePLM then ?>
		
		//constant
		ULRBuf[intindex].L = ULRBuf[intindex].R = *U;

<? else ?>

		//piecewise-linear
		//based on https://arxiv.org/pdf/0804.0402v1.pdf
		
		//calc eigen values and vectors at cell center
		eigen_t eig;
		eigen_forCell_<?=side?>(&eig, U);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig);
		
		//1) calc delta q's ... l r c (eqn 36)
		const global cons_t* UL = U - stepsize[side];
		const global cons_t* UR = U + stepsize[side];
		cons_t dUL, dUR, dUC;
		for (int j = 0; j < numStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
			dUC.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}

		//2) calc eigenspace delta qs (eqn 37)
		real dULEig[numWaves], dUREig[numWaves], dUCEig[numWaves];
		eigen_leftTransform_<?=side?>___(dULEig, &eig, dUL.ptr);
		eigen_leftTransform_<?=side?>___(dUREig, &eig, dUR.ptr);
		eigen_leftTransform_<?=side?>___(dUCEig, &eig, dUC.ptr);

		//3) do the limiter (eqn 38)
		real dUMEig[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			dUMEig[j] = min(
				min(2. * fabs(dULEig[j]),
					2. * fabs(dUREig[j])),
				fabs(dUCEig[j]))
				* (dUCEig[j] >= 0. ? 1. : -1.)
				* (dULEig[j] * dUREig[j] > 0. ? 1. : 0.)
			;
		}
	
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		real pl[numWaves], pr[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			pl[j] = wave[j] >= 0 ? dUMEig[j] * .5 * (1. - dt_dx * wave[j]) : 0;
			pr[j] = wave[j] <= 0 ? dUMEig[j] * .5 * (1. + dt_dx * wave[j]) : 0;
		}

		//convert back
		cons_t ql, qr;
		eigen_rightTransform_<?=side?>___(ql.ptr, &eig, pl);
		eigen_rightTransform_<?=side?>___(qr.ptr, &eig, pr);
		
		for (int j = 0; j < numStates; ++j) {
			ULR->L.ptr[j] = U->ptr[j] - qr.ptr[j];
			ULR->R.ptr[j] = U->ptr[j] + ql.ptr[j];
		}

<? end ?>
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
		const global eigen_t* eig = eigenBuf + intindex;

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
			eigen_rightTransform_<?=side?>__global_(eigenInvCoords, eig, basis);
		
			real newbasis[numWaves];
			eigen_leftTransform_<?=side?>__global_(newbasis, eig, eigenInvCoords);
			
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
			eigen_leftTransform_<?=side?>__global_(eigenCoords, eig, basis);

			real eigenScaled[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				eigenScaled[j] = eigenCoords[j] * wave[j];
			}
			
			real newtransformed[numStates];
			eigen_rightTransform_<?=side?>__global_(newtransformed, eig, eigenScaled);
			
			real transformed[numStates];
			eigen_fluxTransform_<?=side?>__global_(transformed, eig, basis);
			
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

kernel void calcDeltaUEig(
	global real* deltaUEigBuf,
	<?= solver.getULRArg ?>,
	const global eigen_t* eigenBuf
) {
	SETBOUNDS(2,1);	
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		<?= solver.getULRCode ?>

		cons_t deltaU;
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}
	
		int intindex = side + dim * index;	
		global real* deltaUEig = deltaUEigBuf + intindex * numWaves;
		const global eigen_t* eig = eigenBuf + intindex;
		eigen_leftTransform_<?=side?>_global_global_(deltaUEig, eig, deltaU.ptr);
	}<? end ?>
}

<? if solver.fluxLimiter[0] > 0 then ?>
kernel void calcREig(
	global real* rEigBuf,
	const global real* deltaUEigBuf,
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
		global real* rEig = rEigBuf + intindex * numWaves;
		const global real* deltaUEig = deltaUEigBuf + intindex * numWaves;
		const global real* deltaUEigL = deltaUEigBuf + intindexL * numWaves;
		const global real* deltaUEigR = deltaUEigBuf + intindexR * numWaves;
		const global real* wave = waveBuf + intindex * numWaves;
		for (int j = 0; j < numWaves; ++j) {
			if (deltaUEig[j] == 0) {
				rEig[j] = 0;
			} else {
				if (wave[j] >= 0) {
					rEig[j] = deltaUEigL[j] / deltaUEig[j];
				} else {
					rEig[j] = deltaUEigR[j] / deltaUEig[j];
				}
			}
		}
	}<? end ?>
}
<? end ?>

kernel void calcFlux(
	global real* fluxBuf,
	<?= solver.getULRArg ?>,
	const global real* waveBuf, 
	const global eigen_t* eigenBuf, 
	const global real* deltaUEigBuf,
	real dt
<? if solver.fluxLimiter[0] > 0 then ?>
	,const global real* rEigBuf
<? end ?>
) {
	SETBOUNDS(2,1);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / dx<?=side?>_at(i);
		int indexL = index - stepsize[side];

		<?= solver.getULRCode ?>

		cons_t UAvg;
		for (int j = 0; j < numStates; ++j) {
			UAvg.ptr[j] = .5 * (UL->ptr[j] + UR->ptr[j]);
		}
		
		int intindex = side + dim * index;
		const global eigen_t* eig = eigenBuf + intindex;

		real fluxEig[numWaves];
		eigen_leftTransform_<?=side?>__global_(fluxEig, eig, UAvg.ptr);

		const global real* lambdas = waveBuf + numWaves * intindex;
		const global real* deltaUEig = deltaUEigBuf + numWaves * intindex;
<? if solver.fluxLimiter[0] > 0 then ?>
		const global real* rEig = rEigBuf + numWaves * intindex;
<? end ?>

		for (int j = 0; j < numWaves; ++j) {
			real lambda = lambdas[j];
			fluxEig[j] *= lambda;
			real theta = lambda >= 0 ? 1 : -1;
		
<? if solver.fluxLimiter[0] > 0 then ?>
			real phi = fluxLimiter(rEig[j]);
<? end ?>

			real epsilon = lambda * dt_dx;
			real deltaFluxEig = lambda * deltaUEig[j];
			fluxEig[j] -= .5 * deltaFluxEig * (theta
<? if solver.fluxLimiter[0] > 0 then ?>				
				+ phi * (epsilon - theta)
<? end ?>			
			);
		}

		global real* flux = fluxBuf + intindex * numStates;
		eigen_rightTransform_<?=side?>_global_global_(flux, eig, fluxEig);

		real3 interfaceI = _real3(i.x, i.y, i.z);
		interfaceI.s[side] -= .5;
		real3 interfaceX = cell_x(interfaceI);
		real volume = volume_at(interfaceX);

		for (int j = 0; j < numStates; ++j) {
			flux[j] *= volume;
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
