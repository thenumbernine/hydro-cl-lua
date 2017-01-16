// Roe solver:

inline real minmod(real a, real b) {
	if (a * b <= 0) return 0;
	return fabs(a) < fabs(b) ? a : b;
}

inline real maxmod(real a, real b) {
	if (a * b <= 0) return 0;
	return fabs(a) > fabs(b) ? a : b;
}

<? if solver.usePLM then ?>
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
		
		//piecewise-linear
		
#if 0	//Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)
		//and https://en.wikipedia.org/wiki/MUSCL_scheme
		//Works for Euler Sod 1D and 2D
		//Failing for adm1d_v1

		const global cons_t* UL = U - stepsize[side];
		const global cons_t* UR = U + stepsize[side];
		cons_t dUL, dUR, dUC;
		for (int j = 0; j < numStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		}

		cons_t UHalfL, UHalfR;
		for (int j = 0; j < numStates; ++j) {
			
			//Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)	
			//https://en.wikipedia.org/wiki/MUSCL_scheme
			
			real r = dUR.ptr[j] == 0 ? 0 : (dUL.ptr[j] / dUR.ptr[j]);
			real phi = slopeLimiter(r);	//works good with minmod, bad with superbee
			
			real sigma = phi * dUR.ptr[j];
			
			//q^n_i-1/2,R = q^n_i - 1/2 dx sigma	(Hydrodynamics II 6.58)
			UHalfL.ptr[j] = U->ptr[j] - .5 * sigma;
			
			//q^n_i+1/2,L = q^n_i + 1/2 dx sigma	(Hydrodynamics II 6.59)
			UHalfR.ptr[j] = U->ptr[j] + .5 * sigma;
		}
	
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		cons_t FHalfL = fluxForCons_<?=side?>(UHalfL);
		cons_t FHalfR = fluxForCons_<?=side?>(UHalfR);
		
		for (int j = 0; j < numStates; ++j) {
			real dF = FHalfR.ptr[j] - FHalfL.ptr[j];

			//U-cell-L = q^n+1/2_i-1/2,R (Hydrodynamics II 6.62)
			// = q^n_i-1/2,R + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
			ULR->L.ptr[j] = UHalfL.ptr[j] + .5 * dt_dx * dF;

			//U-cell R = q^n+1/2_i+1/2,L (Hydrodynamics II 6.63)
			// = q^n_i+1/2,L + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
			ULR->R.ptr[j] = UHalfR.ptr[j] + .5 * dt_dx * dF;
		}

#elif 1	//based on https://arxiv.org/pdf/0804.0402v1.pdf 
		//and Trangenstein "Numeric Simulation of Hyperbolic Conservation Laws" section 6.2.5
		//except I'm projecting the differences in conservative values instead of primitive values.
		//This also needs modular slope limiter support.
		//This works for adm1d_v1 and 1D Euler Sod 
		//For 2D Euler Sod this gets strange behavior and slowly diverges.

		//1) calc delta q's ... l r c (eqn 36)
		const global cons_t* UL = U - stepsize[side];
		const global cons_t* UR = U + stepsize[side];
		cons_t dUL, dUR, dUC;
		for (int j = 0; j < numStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
			dUC.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}

		//calc eigen values and vectors at cell center
		eigen_t eig;
		eigen_forCell_<?=side?>(&eig, U);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig);
		
		real dULEig[numWaves], dUREig[numWaves], dUCEig[numWaves];
		eigen_leftTransform_<?=side?>___(dULEig, &eig, dUL.ptr);
		eigen_leftTransform_<?=side?>___(dUREig, &eig, dUR.ptr);
		eigen_leftTransform_<?=side?>___(dUCEig, &eig, dUC.ptr);

		real dUMEig[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			//MUSCL slope of characteristic variables
			dUMEig[j] = dULEig[j] * dUREig[j] < 0 ? 0 : (
				(dUCEig[j] >= 0. ? 1. : -1.)
				* min(
					2. * min(
						fabs(dULEig[j]),
						fabs(dUREig[j])),
					fabs(dUCEig[j])
				)
			);	
		}
	
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		real pl[numWaves], pr[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			pl[j] = wave[j] < 0 ? 0 : 
				dUMEig[j] * .5 * (1. - wave[j] * dt_dx);
			pr[j] = wave[j] > 0 ? 0 : 
				dUMEig[j] * .5 * (1. + wave[j] * dt_dx);
		}

		//convert back
		cons_t ql, qr;
		eigen_rightTransform_<?=side?>___(ql.ptr, &eig, pl);
		eigen_rightTransform_<?=side?>___(qr.ptr, &eig, pr);
		
		for (int j = 0; j < numStates; ++j) {
			ULR->L.ptr[j] = U->ptr[j] - qr.ptr[j];
			ULR->R.ptr[j] = U->ptr[j] + ql.ptr[j];
		}
#elif 0	//Trangenstein, Athena, etc, except working on primitives like it says to
		
		const real ePot = 0, ePotL = 0, ePotR = 0;	//TODO fixme
		//this requires standardizing primFromCons (can't pass ePot so easily)
		//...which might mean incorporating ePot into cons_t
		//...which might mean some cons_t variables don't get integrated
		//...which might mean custom mul & add functions for the integrators (to skip the non-integrated fields)
		//...and will mean sizeof(cons_t) >= sizeof(real[numStates])

		//1) calc delta q's ... l r c (eqn 36)
		const global cons_t* UL = U - stepsize[side];
		const global cons_t* UR = U + stepsize[side];
		prim_t W = primFromCons(*U, ePot);
		prim_t WL = primFromCons(*UL, ePotL);
		prim_t WR = primFromCons(*UR, ePotR);
		prim_t dWL, dWR, dWC;
		for (int j = 0; j < numStates; ++j) {
			dWL.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWR.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWC.ptr[j] = WR.ptr[j] - WL.ptr[j];
		}

		//calc eigen values and vectors at cell center
		eigen_t eig;
		eigen_forCell_<?=side?>(&eig, U);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig);
	
		//apply dU/dW before applying left/right eigenvectors so the eigenvectors are of the flux wrt primitives 
		//RW = dW/dU RU, LW = LU dU/dW
		cons_t tmp;
		
		apply_dU_dW(&tmp, &W, &dWL);
		real dWLEig[numWaves];
		eigen_leftTransform_<?=side?>___(dWLEig, &eig, tmp.ptr);
		
		apply_dU_dW(&tmp, &W, &dWR);
		real dWREig[numWaves];
		eigen_leftTransform_<?=side?>___(dWREig, &eig, tmp.ptr);
		
		apply_dU_dW(&tmp, &W, &dWC);
		real dWCEig[numWaves];
		eigen_leftTransform_<?=side?>___(dWCEig, &eig, tmp.ptr);

		real dWMEig[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			//MUSCL slope of characteristic variables
			dWMEig[j] = dWLEig[j] * dWREig[j] < 0 ? 0 : (
				(dWCEig[j] >= 0. ? 1. : -1.)
				* min(
					2. * min(
						fabs(dWLEig[j]),
						fabs(dWREig[j])),
					fabs(dWCEig[j])
				)
			);
		}
	
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		real pl[numWaves], pr[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			pl[j] = wave[j] < 0 ? 0 : 
				dWMEig[j] * .5 * (1. - wave[j] * dt_dx);
			pr[j] = wave[j] > 0 ? 0 : 
				dWMEig[j] * .5 * (1. + wave[j] * dt_dx);
		}

		//convert back
		
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, pl);
		prim_t ql;
		apply_dW_dU(&ql, &W, &tmp);
		
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, pr);
		prim_t qr;
		apply_dW_dU(&qr, &W, &tmp);
	
		prim_t W2L, W2R;
		for (int j = 0; j < numStates; ++j) {
			W2L.ptr[j] = W.ptr[j] - qr.ptr[j];
			W2R.ptr[j] = W.ptr[j] + ql.ptr[j];
		}
		ULR->L = consFromPrim(W2L, ePotL);
		ULR->R = consFromPrim(W2R, ePotR);
#endif

	}<? end ?>
}
<? end ?>

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
	global cons_t* fluxBuf,
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
		
		int intindex = side + dim * index;
		const global eigen_t* eig = eigenBuf + intindex;

		real fluxEig[numWaves];
<? if not solver.eqn.hasFluxFromCons then ?>
		cons_t UAvg;
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
<? if not solver.eqn.hasFluxFromCons then ?>
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
		
		global cons_t* flux = fluxBuf + intindex;
		eigen_rightTransform_<?=side?>_global_global_(flux->ptr, eig, fluxEig);

<? if solver.eqn.hasFluxFromCons then ?>
		cons_t FL = fluxFromCons_<?=side?>(*UL);
		cons_t FR = fluxFromCons_<?=side?>(*UR);
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
		}
	}<? end ?>
}

kernel void calcDerivFromFlux(
	global cons_t* derivBuf,
	const global cons_t* fluxBuf
) {
	SETBOUNDS(2,2);
	global cons_t* deriv = derivBuf + index;
		
	real volume = volume_at(cell_x(i));
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int intindexL = side + dim * index;
		int intindexR = intindexL + dim * stepsize[side]; 
		const global cons_t* fluxL = fluxBuf + intindexL;
		const global cons_t* fluxR = fluxBuf + intindexR;
		for (int j = 0; j < numStates; ++j) {
			real deltaFlux = fluxR->ptr[j] - fluxL->ptr[j];
			deriv->ptr[j] -= deltaFlux / (volume * grid_dx<?=side?>);
		}
	}<? end ?>
}
