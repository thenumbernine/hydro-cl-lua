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
	real3 x = cell_x(i);

	//TODO skip this lr stuff if we're doing piecewise-constant
	//...and just use the original buffers
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexInt = side + dim * index;
		global <?=eqn.consLR_t?>* ULR = ULRBuf + indexInt;	
		
		//piecewise-linear
		
#if 0	//Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)
		//and https://en.wikipedia.org/wiki/MUSCL_scheme
		//Works with oscillations for Euler Sod 1D
		//Works for Euler Sod 2D
		//Works with oscillations for MHD Brio-Wu 1D
		//Works with some oscillations for adm1d_v1 freeflow (fails for mirror)
		//Works for Maxwell

		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
		<?=eqn.cons_t?> dUL, dUR;
		for (int j = 0; j < numIntStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		}
		
		real3 xIntL = x;
		xIntL.s<?=side?> -= grid_dx<?=side?>;
		real3 xIntR = x;
		xIntR.s<?=side?> += grid_dx<?=side?>;

		<?=eqn.cons_t?> UHalfL, UHalfR;
		for (int j = 0; j < numIntStates; ++j) {
			
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

		<?=eqn.cons_t?> FHalfL = fluxFromCons_<?=side?>(UHalfL, xIntL);
		<?=eqn.cons_t?> FHalfR = fluxFromCons_<?=side?>(UHalfR, xIntR);
		
		for (int j = 0; j < numIntStates; ++j) {
			real dF = FHalfR.ptr[j] - FHalfL.ptr[j];

			//U-cell-L = q^n+1/2_i-1/2,R (Hydrodynamics II 6.62)
			// = q^n_i-1/2,R + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
			ULR->L.ptr[j] = UHalfL.ptr[j] + .5 * dt_dx * dF;

			//U-cell R = q^n+1/2_i+1/2,L (Hydrodynamics II 6.63)
			// = q^n_i+1/2,L + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
			ULR->R.ptr[j] = UHalfR.ptr[j] + .5 * dt_dx * dF;
		}

#elif 0	//based on https://arxiv.org/pdf/0804.0402v1.pdf 
		//and Trangenstein "Numeric Simulation of Hyperbolic Conservation Laws" section 6.2.5
		//except I'm projecting the differences in conservative values instead of primitive values.
		//This also needs modular slope limiter support.
		//This works for adm1d_v1 and 1D Euler Sod 
		//(fails for adm1d_v1 mirror)
		//For 2D Euler Sod this gets strange behavior and slowly diverges.
		//This fails for Maxwell 

		//1) calc delta q's ... l r c (eqn 36)
		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
		<?=eqn.cons_t?> dUL, dUR, dUC;
		for (int j = 0; j < numIntStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
			dUC.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}

		real3 xIntL = x;
		xIntL.s<?=side?> -= grid_dx<?=side?>;
		real3 xIntR = x;
		xIntR.s<?=side?> += grid_dx<?=side?>;

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig;
		eigen_forCell_<?=side?>(&eig, U, x);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig, x);
		
		real dULEig[numWaves], dUREig[numWaves], dUCEig[numWaves];
		eigen_leftTransform_<?=side?>___(dULEig, &eig, dUL.ptr, xIntL);
		eigen_leftTransform_<?=side?>___(dUREig, &eig, dUR.ptr, xIntR);
		eigen_leftTransform_<?=side?>___(dUCEig, &eig, dUC.ptr, x);

		//MUSCL slope of characteristic variables
		//TODO make this based on the choice of slope limiter
		real dUMEig[numWaves];
		for (int j = 0; j < numWaves; ++j) {
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

		// slopes in characteristic space
		real aL[numWaves], aR[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			aL[j] = wave[j] < 0 ? 0 : dUMEig[j] * .5 * (1. - wave[j] * dt_dx);
			aR[j] = wave[j] > 0 ? 0 : dUMEig[j] * .5 * (1. + wave[j] * dt_dx);
		}

		//convert back to conservation variable space
		<?=eqn.cons_t?> sL, sR;
		eigen_rightTransform_<?=side?>___(sL.ptr, &eig, aL, x);
		eigen_rightTransform_<?=side?>___(sR.ptr, &eig, aR, x);
		
		for (int j = 0; j < numIntStates; ++j) {
			ULR->L.ptr[j] = U->ptr[j] - sR.ptr[j];
			ULR->R.ptr[j] = U->ptr[j] + sL.ptr[j];
		}
#elif 1	//Trangenstein, Athena, etc, except working on primitives like it says to
		//fails for Maxwell

		real3 xL = x;
		xL.s<?=side?> -= grid_dx<?=side?>;
		real3 xR = x;
		xR.s<?=side?> += grid_dx<?=side?>;
		
		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		
		//1) calc delta q's ... l r c (eqn 36)
		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
		<?=eqn.prim_t?> W = primFromCons(*U, x);
		<?=eqn.prim_t?> WL = primFromCons(*UL, xL);
		<?=eqn.prim_t?> WR = primFromCons(*UR, xR);
		<?=eqn.prim_t?> dWL, dWR, dWC;
		for (int j = 0; j < numIntStates; ++j) {
			dWL.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWR.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWC.ptr[j] = WR.ptr[j] - WL.ptr[j];
		}

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig;
		eigen_forCell_<?=side?>(&eig, U, x);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig, x);
	
		//apply dU/dW before applying left/right eigenvectors so the eigenvectors are of the flux wrt primitives 
		//RW = dW/dU RU, LW = LU dU/dW
		<?=eqn.cons_t?> tmp;

		apply_dU_dW(&tmp, &W, &dWL, xIntL);
		real dWLEig[numWaves];
		eigen_leftTransform_<?=side?>___(dWLEig, &eig, tmp.ptr, xIntL);
		
		apply_dU_dW(&tmp, &W, &dWR, xIntR);
		real dWREig[numWaves];
		eigen_leftTransform_<?=side?>___(dWREig, &eig, tmp.ptr, xIntR);
		
		apply_dU_dW(&tmp, &W, &dWC, x);
		real dWCEig[numWaves];
		eigen_leftTransform_<?=side?>___(dWCEig, &eig, tmp.ptr, x);

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
		
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, pl, xIntL);
		<?=eqn.prim_t?> ql;
		apply_dW_dU(&ql, &W, &tmp, xIntL);
		
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, pr, xIntR);
		<?=eqn.prim_t?> qr;
		apply_dW_dU(&qr, &W, &tmp, xIntR);
	
		<?=eqn.prim_t?> W2L, W2R;
		for (int j = 0; j < numIntStates; ++j) {
			W2L.ptr[j] = W.ptr[j] - qr.ptr[j];
			W2R.ptr[j] = W.ptr[j] + ql.ptr[j];
		}
		ULR->L = consFromPrim(W2L, xIntL);
		ULR->R = consFromPrim(W2R, xIntR);

#elif 0	//based on Athena

		real3 xIntL = x;
		xIntL.s<?=side?> -= grid_dx<?=side?>;
		real3 xIntR = x;
		xIntR.s<?=side?> += grid_dx<?=side?>;
		
		real3 xL = x;
		xL.s<?=side?> -= grid_dx<?=side?>;
		real3 xR = x;
		xR.s<?=side?> += grid_dx<?=side?>;

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig;
		eigen_forCell_<?=side?>(&eig, U, x);
		real wave[numWaves];
		eigen_calcWaves_<?=side?>__(wave, &eig, x);

		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		//1) calc delta q's ... l r c (eqn 36)
		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;

		<?=eqn.prim_t?> W = primFromCons(*U, x);
		<?=eqn.prim_t?> WL = primFromCons(*UL, xL);
		<?=eqn.prim_t?> WR = primFromCons(*UR, xR);
		
		<?=eqn.prim_t?> dWl, dWr, dWc, dWg;
		for (int j = 0; j < numIntStates; ++j) {
			dWl.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWr.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWc.ptr[j] = WR.ptr[j] - WL.ptr[j];
			dWg.ptr[j] = (dWl.ptr[j] * dWr.ptr[j]) <= 0. ? 0. : (
				2. * dWl.ptr[j] * dWr.ptr[j] / (dWl.ptr[j] + dWr.ptr[j])
			);
		}
		
		real dal[numWaves], dar[numWaves], dac[numWaves], dag[numWaves];
		eigen_leftTransform_<?=side?>___(dal, &eig, dWl.ptr, xIntL);
		eigen_leftTransform_<?=side?>___(dar, &eig, dWr.ptr, xIntR);
		eigen_leftTransform_<?=side?>___(dac, &eig, dWc.ptr, x);
		eigen_leftTransform_<?=side?>___(dag, &eig, dWg.ptr, x);

		real da[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			da[j] = 0;
			if (dal[j] * dar[j] > 0) {
				real lim_slope1 = min(fabs(dal[j]), fabs(dar[j]));
				real lim_slope2 = min(.5 * fabs(dac[j]), fabs(dag[j]));
				da[j] = sign(dac[j]) * min(2. * lim_slope1, lim_slope2);
			}
		}
		
		<?=eqn.prim_t?> dWm;
		eigen_rightTransform_<?=side?>___(dWm.ptr, &eig, da, x);

		<?=eqn.prim_t?> Wlv, Wrv;
		for (int j = 0; j < numWaves; ++j) {
			Wlv.ptr[j] = W.ptr[j] - .5 * dWm.ptr[j];
			Wrv.ptr[j] = W.ptr[j] + .5 * dWm.ptr[j];
			real C = Wrv.ptr[j] + Wlv.ptr[j];
			
			Wlv.ptr[j] = max(min(W.ptr[j], WL.ptr[j]), Wlv.ptr[j]);
			Wlv.ptr[j] = min(max(W.ptr[j], WL.ptr[j]), Wlv.ptr[j]);
			Wrv.ptr[j] = C - Wlv.ptr[j];

			Wrv.ptr[j] = max(min(W.ptr[j],WR.ptr[j]),Wrv.ptr[j]);
			Wrv.ptr[j] = min(max(W.ptr[j],WR.ptr[j]),Wrv.ptr[j]);
			Wlv.ptr[j] = C - Wrv.ptr[j];
		}

		ULR->L = consFromPrim(Wlv, xIntL);
		ULR->R = consFromPrim(Wrv, xIntR);

#elif 0	//here's my attempt at Trangenstein section 5.12 PPM
	
		alphas = evL * right side of U[i-1]
		betas = evL * (left side of U[i+1] - right side of U[i-1])
		gamma = 6 * (U[i] - alpha - .5 * beta)
		...but then there are conditions on the definitions of left side and right side of U
			based on beta .... but beta is defined in terms of the left and right sides of U
			... so there's a circular definition.
		Unless (in the text after 5.1) beta is specified for a particular situation
			... if so, then what is the definition of beta?
		hmm, section 5.9.4 defines left and right U as:
			left side of U[i+1/2] = U[i] + .5 (1 - lambda[i] dt/dx) s[j] dx : lambda[i] > 0
									= U[i] : lambda <= 0
			right side of U[i+1/2] = U[i] - .5 (1 - lambda[i+1] dt/dx) s[j] dx : lambda[i+1] < 0
									= U[i] : lambda[i+1] >= 0
		...but this is for muscl, and s[j] is the minmod slope limiter ...
#endif

	}<? end ?>
}
<? end ?>
