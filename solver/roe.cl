// Roe solver:

<? if solver.checkFluxError or solver.checkOrthoError then ?>
kernel void calcErrors(
	global error_t* errorBuf,
	const global real* waveBuf,
	const global <?=eqn.eigen_t?>* eigenBuf
) {
	SETBOUNDS(0,0);

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int intindex = side + dim * index;
		const global real* wave = waveBuf + numWaves * intindex;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + intindex;

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
	const global <?=eqn.eigen_t?>* eigenBuf
) {
	SETBOUNDS(2,1);	
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		<?= solver.getULRCode ?>

		<?=eqn.cons_t?> deltaU;
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}
	
		int intindex = side + dim * index;	
		global real* deltaUEig = deltaUEigBuf + intindex * numWaves;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + intindex;
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
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>,
	const global real* waveBuf, 
	const global <?=eqn.eigen_t?>* eigenBuf, 
	const global real* deltaUEigBuf,
	real dt
<? if solver.fluxLimiter[0] > 0 then ?>
	,const global real* rEigBuf
<? end ?>
) {
	SETBOUNDS(2,1);
	int indexR = index;
//	real3 x = cell_x(i);
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / dx<?=side?>_at(i);//grid_dx<?=side?>;
		int indexL = index - stepsize[side];

		<?= solver.getULRCode ?>
		
		int intindex = side + dim * index;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + intindex;

		real fluxEig[numWaves];
<? if not eqn.hasFluxFromCons then ?>
		<?=eqn.cons_t?> UAvg;
		for (int j = 0; j < numStates; ++j) {
			UAvg.ptr[j] = .5 * (UL->ptr[j] + UR->ptr[j]);
		}
		eigen_leftTransform_<?=side?>__global_(fluxEig, eig, UAvg.ptr);
<? end ?>

		const global real* lambdas = waveBuf + numWaves * intindex;
		const global real* deltaUEig = deltaUEigBuf + numWaves * intindex;
<? if solver.fluxLimiter[0] > 0 then ?>
		const global real* rEig = rEigBuf + numWaves * intindex;
<? end ?>

		for (int j = 0; j < numWaves; ++j) {
			real lambda = lambdas[j];
<? if not eqn.hasFluxFromCons then ?>
			fluxEig[j] *= lambda;
<? else ?>
			fluxEig[j] = 0.;
<? end ?>
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
		
		global <?=eqn.cons_t?>* flux = fluxBuf + intindex;
		eigen_rightTransform_<?=side?>_global_global_(flux->ptr, eig, fluxEig);

<? if eqn.hasFluxFromCons then ?>
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL);
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR);
		for (int j = 0; j < numStates; ++j) {
			flux->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
		}
<? end ?>

		real3 interfaceI = _real3(i.x, i.y, i.z);
		interfaceI.s[side] -= .5;
		real3 interfaceX = cell_x(interfaceI);
		real volume = volume_at(interfaceX);

		for (int j = 0; j < numStates; ++j) {
			flux->ptr[j] *= volume;
//NOTICE
//you only need to do this if you are going to represent your state coordinates in anholonomic normalized vectors
//if you're going to represent them as holonomic coordinate system vector components then don't worry, volume is enough.
//			flux->ptr[j] /= coordHolBasisLen<?=side?>(x);
		}
	}<? end ?>
}

kernel void calcDerivFromFlux(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS(2,2);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
		
	real volume = volume_at(cell_x(i));
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
	
		real3 interfaceIL = _real3(i.x, i.y, i.z);
		interfaceIL.s[side] -= .5;
		real3 interfaceXL = cell_x(interfaceIL);
		real volumeL = volume_at(interfaceXL);

		real3 interfaceIR = _real3(i.x, i.y, i.z);
		interfaceIR.s[side] -= .5;
		real3 interfaceXR = cell_x(interfaceIR);
		real volumeR = volume_at(interfaceXR);
	
		int intindexL = side + dim * index;
		int intindexR = intindexL + dim * stepsize[side]; 
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + intindexL;
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + intindexR;
		for (int j = 0; j < numStates; ++j) {
			real deltaFlux = fluxR->ptr[j] - fluxL->ptr[j];
			//going off of Trangenstein's examples
			//(p.466 for spherical, p.474 for cylindrical) 
			//it looks like, instead of dividing by volume, 
			//I should be dividing by (x^j rhs * volume at rhs - x^j lhs * volume at lhs)
			deriv->ptr[j] -= deltaFlux 
				/ (volumeR * cell_x<?=side?>((real)i.s<?=side?> + .5) 
				- volumeL * cell_x<?=side?>((real)i.s<?=side?> - .5));
		}
	}<? end ?>
}
