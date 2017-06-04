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
	global <?=eqn.consLR_t?>* ULRBuf,
	const global <?=eqn.cons_t?>* UBuf,
	real dt
) {
	SETBOUNDS(1,1);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	//TODO skip this lr stuff if we're doing piecewise-constant
	//...and just use the original buffers
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int intindex = side + dim * index;
		global <?=eqn.consLR_t?>* ULR = ULRBuf + intindex;	
		
		//piecewise-linear
		
#if 0	//Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)
		//and https://en.wikipedia.org/wiki/MUSCL_scheme
		//Works for Euler Sod 1D and 2D
		//Failing for adm1d_v1

		const global <?=eqn.cons_t?>* UL = U - stepsize[side];
		const global <?=eqn.cons_t?>* UR = U + stepsize[side];
		<?=eqn.cons_t?> dUL, dUR, dUC;
		for (int j = 0; j < numStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		}

		<?=eqn.cons_t?> UHalfL, UHalfR;
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

		<?=eqn.cons_t?> FHalfL = fluxForCons_<?=side?>(UHalfL);
		<?=eqn.cons_t?> FHalfR = fluxForCons_<?=side?>(UHalfR);
		
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
		const global <?=eqn.cons_t?>* UL = U - stepsize[side];
		const global <?=eqn.cons_t?>* UR = U + stepsize[side];
		<?=eqn.cons_t?> dUL, dUR, dUC;
		for (int j = 0; j < numStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
			dUC.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig;
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
		<?=eqn.cons_t?> ql, qr;
		eigen_rightTransform_<?=side?>___(ql.ptr, &eig, pl);
		eigen_rightTransform_<?=side?>___(qr.ptr, &eig, pr);
		
		for (int j = 0; j < numStates; ++j) {
			ULR->L.ptr[j] = U->ptr[j] - qr.ptr[j];
			ULR->R.ptr[j] = U->ptr[j] + ql.ptr[j];
		}
#elif 0	//Trangenstein, Athena, etc, except working on primitives like it says to
		
		//1) calc delta q's ... l r c (eqn 36)
		const global <?=eqn.cons_t?>* UL = U - stepsize[side];
		const global <?=eqn.cons_t?>* UR = U + stepsize[side];
		<?=eqn.prim_t?> W = primFromCons(*U);
		<?=eqn.prim_t?> WL = primFromCons(*UL);
		<?=eqn.prim_t?> WR = primFromCons(*UR);
		<?=eqn.prim_t?> dWL, dWR, dWC;
		for (int j = 0; j < numStates; ++j) {
			dWL.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWR.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWC.ptr[j] = WR.ptr[j] - WL.ptr[j];
		}

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig;
		eigen_forCell_<?=side?>(&eig, U);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig);
	
		//apply dU/dW before applying left/right eigenvectors so the eigenvectors are of the flux wrt primitives 
		//RW = dW/dU RU, LW = LU dU/dW
		<?=eqn.cons_t?> tmp;
		
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
		<?=eqn.prim_t?> ql;
		apply_dW_dU(&ql, &W, &tmp);
		
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, pr);
		<?=eqn.prim_t?> qr;
		apply_dW_dU(&qr, &W, &tmp);
	
		<?=eqn.prim_t?> W2L, W2R;
		for (int j = 0; j < numStates; ++j) {
			W2L.ptr[j] = W.ptr[j] - qr.ptr[j];
			W2R.ptr[j] = W.ptr[j] + ql.ptr[j];
		}
		ULR->L = consFromPrim(W2L);
		ULR->R = consFromPrim(W2R);
#endif

	}<? end ?>
}
<? end ?>
