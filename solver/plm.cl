<? if not solver.usePLM then return end ?>

real min3(real x, real y, real z) {
	return min(min(x, y), z);
}

real minmod(real a, real b) {
	if (a * b < 0) return 0;
	return fabs(a) < fabs(b) ? a : b;
}

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
		
		//cell-centered index for a particular side...
		int indexForSide = side + dim * index;
		global <?=eqn.consLR_t?>* ULR = ULRBuf + indexForSide;
		
		//piecewise-linear

<? if solver.usePLM == 'plm-cons' then ?>
/*
#1: slope based on conservative variables
----------------------------------------

phi(x) = minmod(x) = clamp(r, 0, 1)

dUL = U_i - U_i-1
dUR = U_i+1 - U_i
sigma = phi( dUL / dUR ) * dUR

dF = F(U_i+1/2) - F(U_i-1/2)

dx = (x_i+1/2 - x_i-1/2)
U_i-1/2,R = U_i - 1/2 (sigma - dt/dx dF)
U_i+1/2,L = U_i + 1/2 (sigma + dt/dx dF)

Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)
and https://en.wikipedia.org/wiki/MUSCL_scheme
works for euler Sod 1D with oscillations 
works for euler Sod 2D
works for mhd Brio-Wu 1D with oscillations 
works for maxwell with oscillations 
works for adm1d_v1 freeflow with oscillations (fails for mirror)
*/

		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
		<?=eqn.cons_t?> dUL, dUR;
		for (int j = 0; j < numIntStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		}
		for (int j = numIntStates; j < numStates; ++j) {
			dUL.ptr[j] = dUR.ptr[j] = 0;
		}
		
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		
		<?=eqn.cons_t?> UHalfL, UHalfR;
		for (int j = 0; j < numIntStates; ++j) {
			//Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)	
			//https://en.wikipedia.org/wiki/MUSCL_scheme
			
/*
minmod(r) = max(0, min(r, 1))
superbee(r) = max(0, max(
					min(1, 2*r),
					min(2, r)
				))

minmod2(a,b) = a * b < 0  
		? 0
		: (|a| < |b| ? a : b)

a = u[i+1] - u[i]
b = u[i] - u[i-1]
r = a / b = (u[i+1] - u[i]) / (u[i] - u[i-1])

minmod2(a,b) = minmod2(rb, b) = |b| minmod2(r sgn b, sgn b)
TODO oops I forgot r can be negative, and so this sgn b < 0 case is bad ...
for sgn b < 0: 
	minmod2(a,b) = |b| minmod2(-r, -1) = 0
for sgn b >= 0:
	minmod2(a,b) = |b| minmod2(r, 1) = 
		|b| * (r * 1 < 0 ? 0
						: (|r| < 1 ? r : 1))
	minmod2(a,b) = |b| sgn r * max(0, min(|r|, 1))
	minmod2(a,b) = |b| sgn r * minmod(|r|)
*/
#if 1
			real r = dUR.ptr[j] == 0 ? 0 : (dUL.ptr[j] / dUR.ptr[j]);
			real phi = slopeLimiter(r);	//works good with minmod, bad with superbee
			real sigma = phi * dUR.ptr[j];
#else
			real sigma = minmod(minmod(fabs(dUL.ptr[j]), 
                                    fabs(dUR.ptr[j])),
                             fabs(.5 * (dUL.ptr[j] - dUR.ptr[j]))
			) * sign(dUL.ptr[j] - dUR.ptr[j]);
#endif
			//q^n_i-1/2,R = q^n_i - 1/2 dx sigma	(Hydrodynamics II 6.58)
			UHalfL.ptr[j] = U->ptr[j] - .5 * sigma;
			
			//q^n_i+1/2,L = q^n_i + 1/2 dx sigma	(Hydrodynamics II 6.59)
			UHalfR.ptr[j] = U->ptr[j] + .5 * sigma;
		}
		for (int j = numIntStates; j < numStates; ++j) {
			UHalfL.ptr[j] = UHalfR.ptr[j] = U->ptr[j];
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
		for (int j = numIntStates; j < numStates; ++j) {
			ULR->L.ptr[j] = ULR->R.ptr[j] = U->ptr[j];
		}


<? elseif solver.usePLM == 'plm-eig' then ?>
/*
#2: next step, project into eigenspace
----------------------------------------

evL, evR, wave = eigensystem(dF/dU at i)

dULe = evL( U_i - U_i-1 )
dURe = evL( U_i+1 - U_i )
dUCe = evL( (U_i+1 - U_i-1)/2 )
dUMe_j = step(dULe_j * dURe_j) * sign(dUCe_j) * 2 * min(|dULe_j|, |dURe_j|, |dUCe_j|) 

aL = step(wave_j) dUMe_j (1 - wave_j dt/dx)
aR = step(-wave_j) dUMe_j (1 + wave_j dt/dx)

U_i-1/2,R = U_i - evR 1/2 aR
U_i+1/2,L = U_i + evR 1/2 aL

based on https://arxiv.org/pdf/0804.0402v1.pdf 
and Trangenstein "Numeric Simulation of Hyperbolic Conservation Laws" section 6.2.5
except I'm projecting the differences in conservative values instead of primitive values.
This also needs modular slope limiter support.
works for Euler 1D Sod
euler 2D Sod this gets strange behavior and slowly diverges.
works for MHD Brio-Wu
works for maxwell 
works for adm1d_v1 
*/

		//1) calc delta q's ... l r c (eqn 36)
		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
		<?=eqn.cons_t?> dUL, dUR, dUC;
		for (int j = 0; j < numIntStates; ++j) {
			dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
			dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
			dUC.ptr[j] = .5 * (UR->ptr[j] - UL->ptr[j]);
		}
		for (int j = numIntStates; j < numStates; ++j) {
			dUL.ptr[j] = dUR.ptr[j] = dUC.ptr[j] = 0;
		}

		real3 xIntL = x; xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * grid_dx<?=side?>;

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(U, x);
			
		real dULEig[numWaves], dUREig[numWaves], dUCEig[numWaves];
		eigen_leftTransform_<?=side?>___(dULEig, &eig, dUL.ptr, xIntL);
		eigen_leftTransform_<?=side?>___(dUREig, &eig, dUR.ptr, xIntR);
		eigen_leftTransform_<?=side?>___(dUCEig, &eig, dUC.ptr, x);

		//MUSCL slope of characteristic variables
		//TODO make this based on the choice of slope limiter
		real dUMEig[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			//dUMEig[j] = minmod(minmod(2. * dULEig[j], 2. * dUREig[j]), dUCEig[j]);
			dUMEig[j] = dULEig[j] * dUREig[j] < 0 ? 0 : (
				(dUCEig[j] >= 0. ? 1. : -1.)
				* 2. * min3(
					fabs(dULEig[j]),
					fabs(dUREig[j]),
					fabs(dUCEig[j]))
			);	
		}
	
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;
	
		<?=eqn:eigenWaveCodePrefix(side, '&eig', 'x')?>

		// slopes in characteristic space
		real aL[numWaves], aR[numWaves];
		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real wave_j = <?=eqn:eigenWaveCode(side, '&eig', 'x', j)?>;
			aL[j] = wave_j < 0 ? 0 : dUMEig[j] * .5 * (1. - wave_j * dt_dx);
			aR[j] = wave_j > 0 ? 0 : dUMEig[j] * .5 * (1. + wave_j * dt_dx);
		}<? end ?>
		
		//convert back to conservation variable space
		<?=eqn.cons_t?> sL, sR;
		eigen_rightTransform_<?=side?>___(sL.ptr, &eig, aL, x);
		eigen_rightTransform_<?=side?>___(sR.ptr, &eig, aR, x);
		
		for (int j = 0; j < numIntStates; ++j) {
			ULR->L.ptr[j] = U->ptr[j] - sR.ptr[j];
			ULR->R.ptr[j] = U->ptr[j] + sL.ptr[j];
		}
		for (int j = numIntStates; j < numStates; ++j) {
			ULR->L.ptr[j] = ULR->R.ptr[j] = U->ptr[j];
		}

<? elseif solver.usePLM == 'plm-eig-prim' or solver.usePLM == 'plm-eig-prim-ref' then ?>
/*
#3a: next step, convert to primitives
----------------------------------------

evL, evR, wave = eigensystem(dF/dU at i)

W_i = W(U_i)

dWLe = evL dU/dW ( W_i - W_i-1 )
dWRe = evL dU/dW ( W_i+1 - W_i )
dWCe = evL dU/dW ( (W_i+1 - W_i-1)/2 )
dWMe_j = step(dWL_j * dWRe_j) * sign(dWCe_j) * 2 * min(|dWLe_j|, |dWRe_j|, |dWCe_j|)

aL = step(wave_j) dWMe_j (1 - wave_j dt/dx)
aR = step(-wave_j) dWMe_j (1 + wave_j dt/dx)

U_i-1/2,R = U(W_i - 1/2 dW/dU evR aR)
U_i+1/2,L = U(W_i + 1/2 dW/dU evR aL)


#3b: next step, subtract out reference state
----------------------------------------

evL, evR, wave = eigensystem(dF/dU at i)

W_i = W(U_i)

dWLe = evL dU/dW ( W_i - W_i-1 )
dWRe = evL dU/dW ( W_i+1 - W_i )
dWCe = evL dU/dW ( (W_i+1 - W_i-1)/2 )
dWMe_j = step(dWL_j * dWRe_j) * sign(dWCe_j) * 2 * min(|dWLe_j|, |dWRe_j|, |dWCe_j|)

aL = step(wave_j) dWMe_j dt/dx (max(wave) - wave_j)
aR = step(-wave_j) dWMe_j dt/dx (wave_j - min(wave))

U_i-1/2,R = U( W_i - 1/2 dW/dU evR (aR + (1 + dt/dx min(wave)) dWMe) )
U_i+1/2,L = U( W_i + 1/2 dW/dU evR (aL + (1 - dt/dx max(wave)) dWMe) )
		

based on Trangenstein, Athena, etc, except working on primitives like it says to
*/

		real3 xL = x; xL.s<?=side?> -= grid_dx<?=side?>;
		real3 xR = x; xR.s<?=side?> += grid_dx<?=side?>;
		
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		
		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
		
		//1) calc delta q's ... l r c (eqn 36)
		<?=eqn.prim_t?> W = primFromCons(*U, x);
		<?=eqn.prim_t?> WL = primFromCons(*UL, xL);
		<?=eqn.prim_t?> WR = primFromCons(*UR, xR);
		<?=eqn.prim_t?> dWL, dWR, dWC;
		for (int j = 0; j < numIntStates; ++j) {
			dWL.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWR.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWC.ptr[j] = .5 * (WR.ptr[j] - WL.ptr[j]);
		}
		for (int j = numIntStates; j < numStates; ++j) {
			dWL.ptr[j] = dWR.ptr[j] = dWC.ptr[j] = 0.;
		}

		//calc eigen values and vectors at cell center
		//TODO calculate the eigenstate wrt W instead of U - to save some computations
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(U, x);
			
		//apply dU/dW before applying left/right eigenvectors so the eigenvectors are of the flux wrt primitives 
		//RW = dW/dU RU, LW = LU dU/dW
		<?=eqn.cons_t?> tmp;

		apply_dU_dW(&tmp, &W, &dWL, xIntL); real dWLEig[numWaves]; eigen_leftTransform_<?=side?>___(dWLEig, &eig, tmp.ptr, xIntL);
		apply_dU_dW(&tmp, &W, &dWR, xIntR); real dWREig[numWaves]; eigen_leftTransform_<?=side?>___(dWREig, &eig, tmp.ptr, xIntR);
		apply_dU_dW(&tmp, &W, &dWC, x); real dWCEig[numWaves]; eigen_leftTransform_<?=side?>___(dWCEig, &eig, tmp.ptr, x);

		//MUSCL slope of characteristic variables
		real dWMEig[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			dWMEig[j] = dWLEig[j] * dWREig[j] < 0 ? 0 : (
				(dWCEig[j] >= 0. ? 1. : -1.)
				* 2. * min3(
					fabs(dWLEig[j]),
					fabs(dWREig[j]),
					fabs(dWCEig[j])
				)
			);
		}
	
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		<?=eqn:eigenWaveCodePrefix(side, '&eig', 'x')?>

<? 	if solver.usePLM == 'plm-eig-prim' then ?>
		//without reference state

		// calculate left and right slopes in characteristic space
 		real aL[numWaves], aR[numWaves];
 		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real wave_j = <?=eqn:eigenWaveCode(side, '&eig', 'x', j)?>;
			aL[j] = wave_j < 0 ? 0 : dWMEig[j] * .5 * (1. - wave_j * dt_dx);
			aR[j] = wave_j > 0 ? 0 : dWMEig[j] * .5 * (1. + wave_j * dt_dx);
		}<? end ?>

		// transform slopes back to conserved variable space
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, aL, xIntL); <?=eqn.prim_t?> sL; apply_dW_dU(&sL, &W, &tmp, xIntL);
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, aR, xIntR); <?=eqn.prim_t?> sR; apply_dW_dU(&sR, &W, &tmp, xIntR);
	
		// linearly extrapolate the slopes forward and backward from the cell center
 		<?=eqn.prim_t?> W2L, W2R;
 		for (int j = 0; j < numIntStates; ++j) {
			W2L.ptr[j] = W.ptr[j] - sR.ptr[j];
			W2R.ptr[j] = W.ptr[j] + sL.ptr[j];
 		}
		for (int j = numIntStates; j < numStates; ++j) {
			W2L.ptr[j] = W2R.ptr[j] = W.ptr[j];
		}
		ULR->L = consFromPrim(W2L, xIntL);
		ULR->R = consFromPrim(W2R, xIntR);


<?	elseif solver.usePLM == 'plm-eig-prim-ref' then ?>
		//with reference state

		//min and max waves
		//TODO use calcCellMinMaxEigenvalues ... except based on eigen_t
		// so something like calcEigenMinMaxWaves ... 
		real waveMin = min(0., <?=eqn:eigenWaveCode(side, '&eig', 'x', 0)?>);
		real waveMax = max(0., <?=eqn:eigenWaveCode(side, '&eig', 'x', eqn.numWaves-1)?>);

		//limited slope in primitive variable space
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, dWMEig, x);
		<?=eqn.prim_t?> dWM;
		apply_dW_dU(&dWM, &W, &tmp, x);

		//left and right reference states
		<?=eqn.prim_t?> WLRef, WRRef;
		for (int j = 0; j < numIntStates; ++j) {
			WLRef.ptr[j] = W.ptr[j] + .5 * (1. - dt_dx * waveMax) * dWM.ptr[j];
			WRRef.ptr[j] = W.ptr[j] - .5 * (1. + dt_dx * waveMin) * dWM.ptr[j];
		}
		for (int j = numIntStates; j < numStates; ++j) {
			WLRef.ptr[j] = WRRef.ptr[j] = 0;
		}

		// calculate left and right slopes in characteristic space
		real aL[numWaves], aR[numWaves];
		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real wave_j = <?=eqn:eigenWaveCode(side, '&eig', 'x', j)?>;
			aL[j] = wave_j < 0 ? 0 : (dWMEig[j] * dt_dx * (waveMax - wave_j));
			aR[j] = wave_j > 0 ? 0 : (dWMEig[j] * dt_dx * (waveMin - wave_j));
		}<? end ?>

		// transform slopes back to conserved variable space
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, aL, xIntL); <?=eqn.prim_t?> sL; apply_dW_dU(&sL, &W, &tmp, xIntL);
		eigen_rightTransform_<?=side?>___(tmp.ptr, &eig, aR, xIntR); <?=eqn.prim_t?> sR; apply_dW_dU(&sR, &W, &tmp, xIntR);
	
		// linearly extrapolate the slopes forward and backward from the cell center
		<?=eqn.prim_t?> W2L, W2R;
		for (int j = 0; j < numIntStates; ++j) {
			W2R.ptr[j] = WRRef.ptr[j] + .5 * sR.ptr[j];
			W2L.ptr[j] = WLRef.ptr[j] + .5 * sL.ptr[j];
		}
		for (int j = numIntStates; j < numStates; ++j) {
			W2R.ptr[j] = W2L.ptr[j] = W.ptr[j];
		}
		//TODO fix the x's
		ULR->L = consFromPrim(W2R, xIntR);
		ULR->R = consFromPrim(W2L, xIntL);


<? 	end	-- solver.usePLM
elseif solver.usePLM == 'plm-athena' then 
?>
		//based on Athena

		real3 xIntL = x; xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		
		real3 xL = x; xL.s<?=side?> -= grid_dx<?=side?>;
		real3 xR = x; xR.s<?=side?> += grid_dx<?=side?>;

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(U, x);

		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		//1) calc delta q's ... l r c (eqn 36)
		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;

		<?=eqn.prim_t?> W = primFromCons(*U, x);
		<?=eqn.prim_t?> WL = primFromCons(*UL, xL);
		<?=eqn.prim_t?> WR = primFromCons(*UR, xR);
		
		<?=eqn.prim_t?> dWL, dWR, dWC, dWG;
		for (int j = 0; j < numIntStates; ++j) {
			dWL.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWR.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWC.ptr[j] = .5 * (WR.ptr[j] - WL.ptr[j]);
			dWG.ptr[j] = (dWL.ptr[j] * dWR.ptr[j]) <= 0. ? 0. : (
				2. * dWL.ptr[j] * dWR.ptr[j] / (dWL.ptr[j] + dWR.ptr[j])
			);
		}
		for (int j = numIntStates; j < numStates; ++j) {
			dWL.ptr[j] = dWR.ptr[j] = dWC.ptr[j] = dWG.ptr[j] = 0.;
		}
		
		real dal[numWaves], dar[numWaves], dac[numWaves], dag[numWaves];
		eigen_leftTransform_<?=side?>___(dal, &eig, dWL.ptr, xIntL);
		eigen_leftTransform_<?=side?>___(dar, &eig, dWR.ptr, xIntR);
		eigen_leftTransform_<?=side?>___(dac, &eig, dWC.ptr, x);
		eigen_leftTransform_<?=side?>___(dag, &eig, dWG.ptr, x);

		real da[numWaves];
		for (int j = 0; j < numWaves; ++j) {
			da[j] = 0;
			if (dal[j] * dar[j] > 0) {
				real lim_slope1 = min(fabs(dal[j]), fabs(dar[j]));
				real lim_slope2 = min(fabs(dac[j]), fabs(dag[j]));
				da[j] = sign(dac[j]) * min(2. * lim_slope1, lim_slope2);
			}
		}
		
		<?=eqn.prim_t?> dWm;
		eigen_rightTransform_<?=side?>___(dWm.ptr, &eig, da, x);

		<?=eqn.prim_t?> Wlv, Wrv;
		for (int j = 0; j < numIntStates; ++j) {
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
		for (int j = numIntStates; j < numStates; ++j) {
			Wlv.ptr[j] = W.ptr[j];
			Wrv.ptr[j] = W.ptr[j];
		}

		ULR->L = consFromPrim(Wlv, xIntL);
		ULR->R = consFromPrim(Wrv, xIntR);


<? elseif solver.usePLM == 'ppm-experimental' then ?>
//here's my attempt at Trangenstein section 5.12 PPM
//with help from Zingale's ppm code

		real3 xL = x; xL.s<?=side?> -= grid_dx<?=side?>;
		real3 xR = x; xR.s<?=side?> += grid_dx<?=side?>;

		const global <?=eqn.cons_t?>* UL = U - stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + stepsize.s<?=side?>;
	
		//1984 Collela & Woodward eqn 1.7
		<?=eqn.prim_t?> W = primFromCons(*U, x);
		<?=eqn.prim_t?> WL = primFromCons(*UL, xL);
		<?=eqn.prim_t?> WR = primFromCons(*UR, xR);
		<?=eqn.prim_t?> WR2 = primFromCons(*UR, xR);
		<?=eqn.prim_t?> dWL, dWR, dWR2, dWCL, dWCR;
		for (int j = 0; j < numIntStates; ++j) {
			dWL.ptr[j] = W.ptr[j] - WL.ptr[j];
			dWR.ptr[j] = WR.ptr[j] - W.ptr[j];
			dWR2.ptr[j] = WR2.ptr[j] - WR.ptr[j];
			dWCL.ptr[j] = .5 * (WR.ptr[j] - WL.ptr[j]);
			dWCR.ptr[j] = .5 * (WR2.ptr[j] - W.ptr[j]);
		}
		for (int j = numIntStates; j < numStates; ++j) {
			dWL.ptr[j] = dWR.ptr[j] = dWR2.ptr[j] = dWCL.ptr[j] = dWCR.ptr[j] = 0.;
		}

		//slope limiter
		//1984 Collela & Woodward eqn 1.8
		real dWML[numStates];
		for (int j = 0; j < numStates; ++j) {
			dWML[j] = dWL.ptr[j] * dWR.ptr[j] < 0 ? 0 : (
				(dWCL.ptr[j] >= 0. ? 1. : -1.)
				* 2. * min3(
					fabs(dWL.ptr[j]),
					fabs(dWR.ptr[j]),
					fabs(dWCL.ptr[j])
				)
			);
		}
		real dWMR[numStates];
		for (int j = 0; j < numStates; ++j) {
			dWMR[j] = dWR.ptr[j] * dWR2.ptr[j] < 0 ? 0 : (
				(dWCR.ptr[j] >= 0. ? 1. : -1.)
				* 2. * min3(
					fabs(dWR.ptr[j]),
					fabs(dWR2.ptr[j]),
					fabs(dWCR.ptr[j])
				)
			);
		}

		//1984 Collela & Woodward eqn 1.6
		<?=eqn.prim_t?> Wp, Wm;
		for (int j = 0; j < numStates; ++j) {
			Wp.ptr[j] = Wm.ptr[j] = .5 * dWR.ptr[j] - 1./6. * (dWMR[j] - dWML[j]);
		}

		//1984 Collela & Woodward eqn 1.10
		for (int j = 0; j < numStates; ++j) {
			real dWp = Wp.ptr[j] - W.ptr[j];
			real dWm = W.ptr[j] - Wm.ptr[j];
			real dWpm = Wp.ptr[j] - Wm.ptr[j];
			real avgWpm = .5 * (Wp.ptr[j] + Wm.ptr[j]);
			if (dWp * dWm <= 0) {
				Wp.ptr[j] = Wm.ptr[j] = W.ptr[j];
			} else if (dWpm * (W.ptr[j] - avgWpm) > dWpm * dWpm / 6.) {
				Wm.ptr[j] = 3. * W.ptr[j] - 2. * Wp.ptr[j];
			} else if (-dWpm * dWpm / 6. > dWpm * (W.ptr[j] - avgWpm)) {
				Wp.ptr[j] = 3. * W.ptr[j] - 2. * Wm.ptr[j];
			}
		}
		
		for (int j = 0; j < numStates; ++j) {
			Wm.ptr[j] = .5 * (W.ptr[j] +  Wm.ptr[j]);
			Wp.ptr[j] = .5 * (W.ptr[j] +  Wp.ptr[j]);
		}

		//1984 Collela & Woodward eqn 1.5
		real W6[numStates];
		for (int j = 0; j < numStates; ++j) {
			W6[j] = 6. * W.ptr[j] - 3. * (Wp.ptr[j] + Wm.ptr[j]);
		}
		
		real dx = dx<?=side?>_at(i);
		real dt_dx = dt / dx;

		//calc eigen values and vectors at cell center
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(U, x);
		
		<?=eqn:eigenWaveCodePrefix(side, '&eig', 'x')?>
		real eval[numWaves] = {
<? for j=0,eqn.numWaves-1 do
?>			<?=eqn:eigenWaveCode(side, '&eig', 'x', j)?>,
<? end
?>		};

		<?=eqn.prim_t?> Iplus[numWaves];
		<?=eqn.prim_t?> Iminus[numWaves];
		for (int m = 0; m < numWaves; ++m) {
			real sigma = fabs(eval[m]) * dt_dx;
			if (eval[m] >= 0) {
				for (int n = 0; n < numStates; ++n) {
					Iplus[m].ptr[n] = Wp.ptr[n] - .5 * sigma * (Wp.ptr[n] - Wm.ptr[n] - (1. - 2./3.*sigma) * W6[n]);
				}
			} else {
				Iplus[m] = W;
			}
			if (eval[m] <= 0) {
				for (int n = 0; n < numStates; ++n) {
					Iminus[m].ptr[n] = Wm.ptr[n] + .5 * sigma * (Wp.ptr[n] - Wm.ptr[n] + (1. - 2./3.*sigma) * W6[n]);
				}
			} else {
				Iminus[m] = W;
			}
		}

		<?=eqn.prim_t?> Wref_xp = eval[numWaves-1] < 0 ? W : Iplus[numWaves-1];
		<?=eqn.prim_t?> Wref_xm = eval[0] > 0 ? W : Iminus[0];
	
		real beta_xm[numWaves];
		real beta_xp[numWaves];
		for (int m = 0; m < numWaves; ++m) {
			beta_xm[m] = beta_xp[m] = 0;
		}
		for (int m = 0; m < numWaves; ++m) {
			<?=eqn.prim_t?> Wref_xm_minus_Im, Wref_xp_minus_Ip;
			for (int n = 0; n < numStates; ++n) {
				Wref_xm_minus_Im.ptr[n] = Wref_xm.ptr[n] - Iminus[m].ptr[n];
				Wref_xp_minus_Ip.ptr[n] = Wref_xp.ptr[n] - Iplus[m].ptr[n];
			}
			
			if (eval[m] >= 0) {
				real beta_xp_add[numWaves];
				eigen_leftTransform_<?=side?>___(beta_xp_add, &eig, Wref_xm_minus_Im.ptr, x);
				for (int n = 0; n < numWaves; ++n) {
					beta_xp[n] += beta_xp_add[n];
				}
			}
			
			if (eval[m] <= 0) {
				real beta_xm_add[numWaves];
				eigen_leftTransform_<?=side?>___(beta_xm_add, &eig, Wref_xp_minus_Ip.ptr, x);
				for (int n = 0; n < numWaves; ++n) {
					beta_xm[n] += beta_xm_add[n];
				}
			}
		}
		
		for (int j = 0; j < numWaves; ++j) {
			beta_xp[j] = max(0., beta_xp[j]);
			beta_xm[j] = min(0., beta_xm[j]);
		}

		<?=eqn.prim_t?> WresL, WresR;
		eigen_rightTransform_<?=side?>___(WresL.ptr, &eig, beta_xp, x);
		eigen_rightTransform_<?=side?>___(WresR.ptr, &eig, beta_xm, x);
		
		for (int j = 0; j < numStates; ++j) {
			WresL.ptr[j] = Wref_xp.ptr[j] - WresL.ptr[j];
			WresR.ptr[j] = Wref_xm.ptr[j] - WresR.ptr[j];
		}
		
		ULR->L = consFromPrim(WresL, x);
		ULR->R = consFromPrim(WresR, x);

<? end ?>

	}<? end ?>
}
