//// MODULE_NAME: calcLR
//// MODULE_DEPENDS: consLR_t solver_t cons_t normal_t cell_x cell_dx#

// TODO incorporate parallel propagators

<? if not solver.usePLM then return end ?>

real min3(real x, real y, real z) {
	return min(min(x, y), z);
}

real minmod(real a, real b) {
	if (a * b < 0) return 0;
	return fabs(a) < fabs(b) ? a : b;
}

<? 
for side=0,solver.dim-1 do

	-- piecewise-linear
	if solver.usePLM == 'plm-cons' then 
?>
//// MODULE_DEPENDS: fluxFromCons slopeLimiter

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
<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;
	
	real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
	<?=eqn.cons_t?> UHalfL, UHalfR;
	for (int j = 0; j < numIntStates; ++j) {
		real dUL = U->ptr[j] - UL->ptr[j];
		real dUR = UR->ptr[j] - U->ptr[j];
		
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
		real r = dUR == 0 ? 0 : (dUL / dUR);
		real phi = slopeLimiter(r);	//works good with minmod, bad with superbee
		real sigma = phi * dUR;
#else	//isn't as accurate anyways
		real sigma = minmod(minmod(fabs(dUL), 
			fabs(dUR.ptr[j])),
			fabs(.5 * (dUL - dUR))
		) * sign(dUL - dUR);
#endif
		//q^n_i-1/2,R = q^n_i - 1/2 dx sigma	(Hydrodynamics II 6.58)
		UHalfL.ptr[j] = U->ptr[j] - .5 * sigma;
		
		//q^n_i+1/2,L = q^n_i + 1/2 dx sigma	(Hydrodynamics II 6.59)
		UHalfR.ptr[j] = U->ptr[j] + .5 * sigma;
	}
	for (int j = numIntStates; j < numStates; ++j) {
		UHalfL.ptr[j] = UHalfR.ptr[j] = U->ptr[j];
	}

	real dx = cell_dx<?=side?>(x);
	real dt_dx = dt / dx;

	<?=eqn.cons_t?> FHalfL = fluxFromCons(solver, UHalfL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=eqn.cons_t?> FHalfR = fluxFromCons(solver, UHalfR, xIntR, normal_forSide<?=side?>(xIntR));

	<?=eqn.consLR_t?> ULR = {
		.L = *U,
		.R = *U,
	};
	for (int j = 0; j < numIntStates; ++j) {
		real dF = FHalfR.ptr[j] - FHalfL.ptr[j];

		//U-cell-L = q^n+1/2_i-1/2,R (Hydrodynamics II 6.62)
		// = q^n_i-1/2,R + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
		ULR.L.ptr[j] = UHalfL.ptr[j] + .5 * dt_dx * dF;

		//U-cell R = q^n+1/2_i+1/2,L (Hydrodynamics II 6.63)
		// = q^n_i+1/2,L + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
		ULR.R.ptr[j] = UHalfR.ptr[j] + .5 * dt_dx * dF;
	}
	return ULR;
}

<? 
	-- technically i should call the one above 'plm-cons-with-flux' 
	-- and this one 'plm-cons'
	elseif solver.usePLM == 'plm-cons-alone' then 
?>
//// MODULE_DEPENDS: slopeLimiter

<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	// extrapolate slopes in consered variable space
	real dx = solver->grid_dx.s<?=side?>;

	real3 xL = x; xL.s<?=side?> -= dx;
	real3 xR = x; xR.s<?=side?> += dx;

	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;

	<?=eqn.consLR_t?> ULR = {
		.L = *U,
		.R = *U,
	};
	// minmod
	for (int j = 0; j < numIntStates; ++j) {
#if 1	// taken from above / from Dullemond "Hydrodynamics II"
		real dUR = UR->ptr[j] - U->ptr[j];
		real dUL = U->ptr[j] - UL->ptr[j];
		real r = dUR == 0 ? 0 : (dUL / dUR);
		real phi = slopeLimiter(r);	//works good with minmod, bad with superbee
		real sigma = phi * dUR;
		ULR.L.ptr[j] -= .5 * sigma;
		ULR.R.ptr[j] += .5 * sigma;
#endif
#if 0 	// Based on Mara's PLM.  Looks to be asymmetric, favoring the x+ and y+ directions. Maybe I copied it wrong.
		real dUL = U->ptr[j] - UL->ptr[j];
		real dUC = .5 * (UR->ptr[j] - UL->ptr[j]);
		real dUR = UR->ptr[j] - U->ptr[j];
		real adwl = fabs(dUL);
		real adwr = fabs(dUR);
		real adwc = fabs(dUC);
		real mindw = min3(adwl, adwr, adwc);
		real dUM = .25 * fabs(sign(dUL) + sign(dUC)) * fabs(sign(dUL) + sign(dUR)) * mindw;
		ULR.L.ptr[j] -= .5 * dUM; 
		ULR.R.ptr[j] += .5 * dUM; 
#endif
	}
	return ULR;
}

<? 
	elseif solver.usePLM == 'plm-prim-alone' then 
?>
//// MODULE_DEPENDS: slopeLimiter consFromPrim

<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	// extrapolate slopes in primitive space
	real dx = solver->grid_dx.s<?=side?>;

	real3 xL = x; xL.s<?=side?> -= dx;
	real3 xR = x; xR.s<?=side?> += dx;

	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;

	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	<?=eqn.prim_t?> WL = primFromCons(solver, *UL, xL);
	<?=eqn.prim_t?> WR = primFromCons(solver, *UR, xR);
	
	// minmod
	<?=eqn.prim_t?> nWL = W;
	<?=eqn.prim_t?> nWR = W;
	for (int j = 0; j < numIntStates; ++j) {
#if 1	// taken from above / from Dullemond "Hydrodynamics II"
		real dWR = WR.ptr[j] - W.ptr[j];
		real dWL = W.ptr[j] - WL.ptr[j];
		real r = dWR == 0 ? 0 : (dWL / dWR);
		real phi = slopeLimiter(r);
		real sigma = phi * dWR;
		nWL.ptr[j] -= .5 * sigma;
		nWR.ptr[j] += .5 * sigma;
#endif
#if 0	// Based on Mara's PLM.  Looks to be asymmetric, favoring the x+ and y+ directions. Maybe I copied it wrong.
		real dWL = W.ptr[j] - WL.ptr[j];	//TODO SOR or whatever ... scalar here x2?
		real dWC = .5 * (WR.ptr[j] - WL.ptr[j]);
		real dWR = WR.ptr[j] - W.ptr[j];
		real adwl = fabs(dWL);
		real adwr = fabs(dWR);
		real adwc = fabs(dWC);
		real mindw = min(min(adwl, adwr), adwc);
		real dWM = .25 * fabs(sign(dWL) + sign(dWC)) * fabs(sign(dWL) + sign(dWR)) * mindw;
		nWL.ptr[j] -= .5 * dWM; 
		nWR.ptr[j] += .5 * dWM; 
#endif	
	}
	
	real3 xIntL = x; xIntL.s<?=side?> -= .5 * dx;
	real3 xIntR = x; xIntR.s<?=side?> += .5 * dx;

	//converting from cons->prim and then prim->cons might lose us accuracy? 
	return (<?=eqn.consLR_t?>){
		.L = consFromPrim(solver, nWL, xIntL),
		.R = consFromPrim(solver, nWR, xIntR),
	};
}

<? 
	elseif solver.usePLM == 'plm-eig' then 
?>
//// MODULE_DEPENDS: eigen_forCell eigen_left/rightTransform

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

<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	//1) calc delta q's ... l r c (eqn 36)
	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;
	<?=eqn.cons_t?> dUL, dUR, dUC;
	for (int j = 0; j < numIntStates; ++j) {
		dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
		dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		dUC.ptr[j] = .5 * (UR->ptr[j] - UL->ptr[j]);
	}
	for (int j = numIntStates; j < numStates; ++j) {
		dUL.ptr[j] = dUR.ptr[j] = dUC.ptr[j] = 0;
	}

	real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

	//calc eigen values and vectors at cell center
	<?=eqn.eigen_t?> eig = eigen_forCell(solver, *U, x, n);
		
	<?=eqn.waves_t?> dULEig = eigen_leftTransform(solver, eig, dUL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=eqn.waves_t?> dUREig = eigen_leftTransform(solver, eig, dUR, xIntR, normal_forSide<?=side?>(xIntR));
	<?=eqn.waves_t?> dUCEig = eigen_leftTransform(solver, eig, dUC, x, n);

	//MUSCL slope of characteristic variables
	//TODO make this based on the choice of slope limiter
	<?=eqn.waves_t?> dUMEig;
	for (int j = 0; j < numWaves; ++j) {

// This is adapted from above to use the modular slopeLimiter
#if 0
		real r = dUREig.ptr[j] == 0 ? 0 : (dULEig.ptr[j] / dUREig.ptr[j]);
		real phi = slopeLimiter(r);
		real sigma = phi * dUREig.ptr[j];
		dUMEig.ptr[j] = sigma;

#endif

//This one is symmetric, unlike the above center-slope-based methods.
// It also takes into account the center slope, unlike the slopeLimiter() options.
// However it is fixed to minmod and does not use the modular slopeLimiter() option.
#if 1	
		//dUMEig.ptr[j] = minmod(minmod(2. * dULEig.ptr[j], 2. * dUREig.ptr[j]), dUCEig.ptr[j]);
		dUMEig.ptr[j] = dULEig.ptr[j] * dUREig.ptr[j] < 0 ? 0 : (
			(dUCEig.ptr[j] >= 0. ? 1. : -1.)
			* 2. * min3(
				fabs(dULEig.ptr[j]),
				fabs(dUREig.ptr[j]),
				fabs(dUCEig.ptr[j]))
		);
#endif	
	}

	real dx = cell_dx<?=side?>(x);
	real dt_dx = dt / dx;

	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'x')?>

	// slopes in characteristic space
	<?=eqn.waves_t?> aL;
	<?=eqn.waves_t?> aR;
	<? for j=0,eqn.numWaves-1 do ?>{
		const int j = <?=j?>;
		real wave_j = <?=eqn:eigenWaveCode('n', 'eig', 'x', j)?>;
		aL.ptr[j] = wave_j < 0 ? 0 : dUMEig.ptr[j] * .5 * (1. - wave_j * dt_dx);
		aR.ptr[j] = wave_j > 0 ? 0 : dUMEig.ptr[j] * .5 * (1. + wave_j * dt_dx);
	}<? end ?>
	
	//convert back to conservation variable space
	<?=eqn.cons_t?> sL = eigen_rightTransform(solver, eig, aL, x, n);
	<?=eqn.cons_t?> sR = eigen_rightTransform(solver, eig, aR, x, n);

	<?=eqn.consLR_t?> ULR = {
		.L = *U,
		.R = *U,
	};
	for (int j = 0; j < numIntStates; ++j) {
		ULR.L.ptr[j] += sL.ptr[j];
		ULR.R.ptr[j] -= sR.ptr[j];
	}
	return ULR;
}

<? 
	elseif solver.usePLM == 'plm-eig-prim' 
	or solver.usePLM == 'plm-eig-prim-ref' 
	then 
?>

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

<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	real3 xL = x; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
	real3 xR = x; xR.s<?=side?> += solver->grid_dx.s<?=side?>;
	
	real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;
	
	//1) calc delta q's ... l r c (eqn 36)
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	<?=eqn.prim_t?> WL = primFromCons(solver, *UL, xL);
	<?=eqn.prim_t?> WR = primFromCons(solver, *UR, xR);
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
	<?=eqn.eigen_t?> eig = eigen_forCell(solver, *U, x, n);
		
	//apply dU/dW before applying left/right eigenvectors so the eigenvectors are of the flux wrt primitives 
	//RW = dW/dU RU, LW = LU dU/dW
	<?=eqn.cons_t?> tmp;

	tmp = apply_dU_dW(solver, W, dWL, xIntL); 
	<?=eqn.waves_t?> dWLEig = eigen_leftTransform(solver, eig, tmp, xIntL, normal_forSide<?=side?>(xIntL));
	
	tmp = apply_dU_dW(solver, W, dWR, xIntR); 
	<?=eqn.waves_t?> dWREig = eigen_leftTransform(solver, eig, tmp, xIntR, normal_forSide<?=side?>(xIntR));
	
	tmp = apply_dU_dW(solver, W, dWC, x); 
	<?=eqn.waves_t?> dWCEig = eigen_leftTransform(solver, eig, tmp, x, n);

	//MUSCL slope of characteristic variables
	<?=eqn.waves_t?> dWMEig;
	for (int j = 0; j < numWaves; ++j) {
		dWMEig.ptr[j] = dWLEig.ptr[j] * dWREig.ptr[j] < 0. ? 0. : (
			(dWCEig.ptr[j] >= 0. ? 1. : -1.)
			* 2. * min3(
				fabs(dWLEig.ptr[j]),
				fabs(dWREig.ptr[j]),
				fabs(dWCEig.ptr[j])
			)
		);
	}

	real dx = cell_dx<?=side?>(x);
	real dt_dx = dt / dx;

	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'x')?>

<? 	
		if solver.usePLM == 'plm-eig-prim' then 
?>
//// MODULE_DEPENDS: apply_dU_dW apply_dW_dU eigen_left/rightTransform

	//without reference state

	// calculate left and right slopes in characteristic space
	<?=eqn.waves_t?> aL, aR;
	<? for j=0,eqn.numWaves-1 do ?>{
		const int j = <?=j?>;
		real wave_j = <?=eqn:eigenWaveCode(side, 'eig', 'x', j)?>;
		aL.ptr[j] = wave_j < 0 ? 0 : dWMEig.ptr[j] * .5 * (1. - wave_j * dt_dx);
		aR.ptr[j] = wave_j > 0 ? 0 : dWMEig.ptr[j] * .5 * (1. + wave_j * dt_dx);
	}<? end ?>

	// transform slopes back to conserved variable space
	tmp = eigen_rightTransform(solver, eig, aL, xIntL, normal_forSide<?=side?>(xIntL)); 
	<?=eqn.prim_t?> sL = apply_dW_dU(solver, W, tmp, xIntL);
	
	tmp = eigen_rightTransform(solver, eig, aR, xIntR, normal_forSide<?=side?>(xIntR)); 
	<?=eqn.prim_t?> sR = apply_dW_dU(solver, W, tmp, xIntR);

	// linearly extrapolate the slopes forward and backward from the cell center
	<?=eqn.prim_t?> W2L, W2R;
	for (int j = 0; j < numIntStates; ++j) {
		W2L.ptr[j] = W.ptr[j] + sL.ptr[j];
		W2R.ptr[j] = W.ptr[j] - sR.ptr[j];
	}
	for (int j = numIntStates; j < numStates; ++j) {
		W2L.ptr[j] = W2R.ptr[j] = W.ptr[j];
	}
	return (<?=eqn.consLR_t?>){
		.L = consFromPrim(solver, W2L, xIntL),
		.R = consFromPrim(solver, W2R, xIntR),
	};
<?	
		elseif solver.usePLM == 'plm-eig-prim-ref' then 
?>
//// MODULE_DEPENDS: eigen_forCell apply_dU_dW apply_dW_dU eigen_left/rightTransform

	//with reference state

	//min and max waves
	//TODO use calcCellMinMaxEigenvalues ... except based on eigen_t
	// so something like calcEigenMinMaxWaves ... 
	real waveMin = min((real)0., <?=eqn:eigenWaveCode('n', 'eig', 'x', 0)?>);
	real waveMax = max((real)0., <?=eqn:eigenWaveCode('n', 'eig', 'x', eqn.numWaves-1)?>);

	//limited slope in primitive variable space
	tmp = eigen_rightTransform(solver, eig, dWMEig, x, n);
	<?=eqn.prim_t?> dWM = apply_dW_dU(solver, W, tmp, x);

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
	<?=eqn.waves_t?> aL, aR;
	<? for j=0,eqn.numWaves-1 do ?>{
		const int j = <?=j?>;
		real wave_j = <?=eqn:eigenWaveCode(side, 'eig', 'x', j)?>;
		aL.ptr[j] = wave_j < 0 ? 0 : (dWMEig.ptr[j] * dt_dx * (waveMax - wave_j));
		aR.ptr[j] = wave_j > 0 ? 0 : (dWMEig.ptr[j] * dt_dx * (waveMin - wave_j));
	}<? end ?>

	// transform slopes back to conserved variable space
	tmp = eigen_rightTransform(solver, eig, aL, xIntL, normal_forSide<?=side?>(xIntL)); 
	<?=eqn.prim_t?> sL = apply_dW_dU(solver, W, tmp, xIntL);
	
	tmp = eigen_rightTransform(solver, eig, aR, xIntR, normal_forSide<?=side?>(xIntR)); 
	<?=eqn.prim_t?> sR = apply_dW_dU(solver, W, tmp, xIntR);

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
	return (<?=eqn.consLR_t?>){
		.L = consFromPrim(solver, W2L, xIntL),
		.R = consFromPrim(solver, W2R, xIntR),
	};
<?
		end	-- solver.usePLM
?>
}

<?

	elseif solver.usePLM == 'plm-athena' then 

?>
//// MODULE_DEPENDS: eigen_forCell eigen_left/rightTransform

//based on Athena
<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
	real3 xL = x; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
	real3 xR = x; xR.s<?=side?> += solver->grid_dx.s<?=side?>;

	//calc eigen values and vectors at cell center
	<?=eqn.eigen_t?> eig = eigen_forCell(solver, *U, x, n);

	real dx = cell_dx<?=side?>(x);
	real dt_dx = dt / dx;

	//1) calc delta q's ... l r c (eqn 36)
	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;

	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	<?=eqn.prim_t?> WL = primFromCons(solver, *UL, xL);
	<?=eqn.prim_t?> WR = primFromCons(solver, *UR, xR);
	
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

	//TODO shouldn't this be the dA/dW eigen transformation?
	<?=eqn.waves_t?> dal = eigen_leftTransform(solver, eig, *(<?=eqn.cons_t?>*)&dWL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=eqn.waves_t?> dar = eigen_leftTransform(solver, eig, *(<?=eqn.cons_t?>*)&dWR, xIntR, normal_forSide<?=side?>(xIntR));
	<?=eqn.waves_t?> dac = eigen_leftTransform(solver, eig, *(<?=eqn.cons_t?>*)&dWC, x, n);
	<?=eqn.waves_t?> dag = eigen_leftTransform(solver, eig, *(<?=eqn.cons_t?>*)&dWG, x, n);

	<?=eqn.waves_t?> da;
	for (int j = 0; j < numWaves; ++j) {
		da.ptr[j] = 0;
		if (dal.ptr[j] * dar.ptr[j] > 0) {
			real lim_slope1 = min(fabs(dal.ptr[j]), fabs(dar.ptr[j]));
			real lim_slope2 = min(fabs(dac.ptr[j]), fabs(dag.ptr[j]));
			da.ptr[j] = sign(dac.ptr[j]) * min((real)2. * lim_slope1, lim_slope2);
		}
	}

	<?=eqn.cons_t?> dWm_tmp = eigen_rightTransform(solver, eig, da, x, n);
	<?=eqn.prim_t?> dWm = *(<?=eqn.prim_t?>*)&dWm_tmp;

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

	return (<?=eqn.consLR_t?>){
		.L = consFromPrim(solver, Wrv, xIntL),
		.R = consFromPrim(solver, Wlv, xIntR),
	};
}

<? 
	elseif solver.usePLM == 'ppm-experimental' then 
?>
//// MODULE_DEPENDS: eigen_forCell eigen_left/rightTransform

//here's my attempt at Trangenstein section 5.12 PPM
//with help from Zingale's ppm code
<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {
	real3 xL = x; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
	real3 xR = x; xR.s<?=side?> += solver->grid_dx.s<?=side?>;

	const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;
	const global <?=eqn.cons_t?>* UR2 = UR + solver->stepsize.s<?=side?>;

	<?=eqn.prim_t?> WL = primFromCons(solver, *UL, xL);
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	<?=eqn.prim_t?> WR = primFromCons(solver, *UR, xR);
	<?=eqn.prim_t?> WR2 = primFromCons(solver, *UR2, xR);

	<?=eqn.prim_t?> Qplus, Qminus;
	for (int j = 0; j < numStates; ++j) {
		real dq0 = .5 * (WR.ptr[j] - WL.ptr[j]);
		real dqp = .5 * (WR2.ptr[j] - W.ptr[j]);
	
		if ((WR.ptr[j] - W.ptr[j]) * (W.ptr[j] - WL.ptr[j]) > 0) {
			dq0 = (dq0 > 0 ? 1. : -1.) * min(
				fabs(dq0), min(
					2. * fabs(W.ptr[j] - WL.ptr[j]),
					2. * fabs(WR.ptr[j] - W.ptr[j])
				)
			);
		} else {
			dq0 = 0.;
		}
	
		if ((WR2.ptr[j] - WR.ptr[j]) * (WR.ptr[j] - W.ptr[j]) > 0) {
			dqp = (dqp > 0 ? 1. : -1.) * min(
				fabs(dqp), min(
					2. * fabs(WR.ptr[j] - W.ptr[j]),
					2. * fabs(WR2.ptr[j] - WR.ptr[j])
				)
			);
		} else {
			dqp = 0.;
		}

		//Qminus[i+1]
		Qplus.ptr[j] = Qminus.ptr[j] = .5 * (W.ptr[j] + WR.ptr[j]) - 1./6. * (dqp - dq0);
	
		Qplus.ptr[j] = max(Qplus.ptr[j], min(W.ptr[j], WR.ptr[j]));
		Qplus.ptr[j] = min(Qplus.ptr[j], max(W.ptr[j], WR.ptr[j]));
	
		Qminus.ptr[j] = max(Qminus.ptr[j], min(W.ptr[j], WR.ptr[j]));
		Qminus.ptr[j] = min(Qminus.ptr[j], max(W.ptr[j], WR.ptr[j]));
	}

	for (int j = 0; j < numStates; ++j) {
		real dWR = WR.ptr[j] - W.ptr[j];
		real dWL = W.ptr[j] - WL.ptr[j];
		real dWRm = WR.ptr[j] - WL.ptr[j];
		real avgWRm = .5 * (WR.ptr[j] + WL.ptr[j]);
		if (dWR * dWL <= 0) {
			WR.ptr[j] = WL.ptr[j] = W.ptr[j];
		} else if (dWRm * (W.ptr[j] - avgWRm) > dWRm * dWRm / 6.) {
			WL.ptr[j] = 3. * W.ptr[j] - 2. * WR.ptr[j];
		} else if (-dWRm * dWRm / 6. > dWRm * (W.ptr[j] - avgWRm)) {
			WR.ptr[j] = 3. * W.ptr[j] - 2. * WL.ptr[j];
		}
	}
	
	for (int j = 0; j < numStates; ++j) {
		WL.ptr[j] = .5 * (W.ptr[j] +  WL.ptr[j]);
		WR.ptr[j] = .5 * (W.ptr[j] +  WR.ptr[j]);
	}

	//1984 Collela & Woodward eqn 1.5
	real W6[numStates];
	for (int j = 0; j < numStates; ++j) {
		W6[j] = 6. * W.ptr[j] - 3. * (WR.ptr[j] + WL.ptr[j]);
	}
	
	real dx = cell_dx<?=side?>(x);
	real dt_dx = dt / dx;

	//calc eigen values and vectors at cell center
	<?=eqn.eigen_t?> eig = eigen_forCell(solver, *U, x, n);
	
	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'x')?>
	real eval[numWaves] = {
<? for j=0,eqn.numWaves-1 do
?>		<?=eqn:eigenWaveCode('n', 'eig', 'x', j)?>,
<? end
?>	};

	<?=eqn.prim_t?> Iplus[numWaves];
	<?=eqn.prim_t?> Iminus[numWaves];
	for (int p = 0; p < numWaves; ++p) {
		real sigma = fabs(eval[p]) * dt_dx;
		if (eval[p] >= 0) {
			for (int q = 0; q < numStates; ++q) {
				Iplus[p].ptr[q] = WR.ptr[q] - .5 * sigma * (WR.ptr[q] - WL.ptr[q] - (1. - 2./3.*sigma) * W6[q]);
			}
		} else {
			Iplus[p] = W;
		}
		if (eval[p] <= 0) {
			for (int q = 0; q < numStates; ++q) {
				Iminus[p].ptr[q] = WL.ptr[q] + .5 * sigma * (WR.ptr[q] - WL.ptr[q] + (1. - 2./3.*sigma) * W6[q]);
			}
		} else {
			Iminus[p] = W;
		}
	}

	<?=eqn.prim_t?> Wref_xp = eval[numWaves-1] < 0 ? W : Iplus[numWaves-1];
	<?=eqn.prim_t?> Wref_xm = eval[0] > 0 ? W : Iminus[0];

	//TODO eigenvectors wrt dA/dW 
	<?=eqn.waves_t?> beta_xm, beta_xp;
	for (int p = 0; p < numWaves; ++p) {
		beta_xm.ptr[p] = beta_xp.ptr[p] = 0;
	}
	for (int p = 0; p < numWaves; ++p) {
		<?=eqn.prim_t?> Wref_xm_minus_Im, Wref_xp_minus_Ip;
		for (int q = 0; q < numStates; ++q) {
			Wref_xm_minus_Im.ptr[q] = Wref_xm.ptr[q] - Iminus[p].ptr[q];
			Wref_xp_minus_Ip.ptr[q] = Wref_xp.ptr[q] - Iplus[p].ptr[q];
		}
		
		if (eval[p] >= 0) {
			<?=eqn.waves_t?> beta_xp_add = eigen_leftTransform(solver, eig, *(<?=eqn.cons_t?>*)&Wref_xm_minus_Im, x, n);
			for (int q = 0; q < numWaves; ++q) {
				beta_xp.ptr[q] += beta_xp_add.ptr[q];
			}
		}
		
		if (eval[p] <= 0) {
			<?=eqn.waves_t?> beta_xm_add = eigen_leftTransform(solver, eig, *(<?=eqn.cons_t?>*)&Wref_xp_minus_Ip, x, n);
			for (int q = 0; q < numWaves; ++q) {
				beta_xm.ptr[q] += beta_xm_add.ptr[q];
			}
		}
	}
	
	for (int j = 0; j < numWaves; ++j) {
		beta_xp.ptr[j] = max(0., beta_xp.ptr[j]);
		beta_xm.ptr[j] = min(0., beta_xm.ptr[j]);
	}

	<?=eqn.cons_t?> tmp;
	tmp = eigen_rightTransform(solver, eig, beta_xp, x, n);
	<?=eqn.prim_t?> WresL = *(<?=eqn.prim_t?>*)&tmp;
	tmp = eigen_rightTransform(solver, eig, beta_xm, x, n);
	<?=eqn.prim_t?> WresR = *(<?=eqn.prim_t?>*)&tmp;
	
	for (int j = 0; j < numIntStates; ++j) {
		WresL.ptr[j] = Wref_xp.ptr[j] - WresL.ptr[j];
		WresR.ptr[j] = Wref_xm.ptr[j] - WresR.ptr[j];
	}
	for (int j = numIntStates; j < numStates; ++j) {
		WresL.ptr[j] = WresR.ptr[j] = 0;
	}

	return (<?=eqn.consLR_t?>){
		.L = consFromPrim(solver, WresL, x),
		.R = consFromPrim(solver, WresR, x),
	};
}

<?

	elseif solver.usePLM == 'plm-athena' then 

-- TODO correlate between calcCellLR's U access and SETBOUNDS of calcLR()
-- to make sure there's no OOB reads
?>

<?=eqn.consLR_t?> calcCellLR_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real dt,
	real3 x,
	int4 i,
	normal_t n
) {


}

<?
	end
end
?>

kernel void calcLR(
	constant <?=solver.solver_t?>* solver,
	const global <?=solver.coord.cell_t?>* cellBuf,
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
		normal_t n = normal_forSide<?=side?>(x);
		
		//cell-centered index for a particular side...
		int indexForSide = <?=side?> + dim * index;
		ULRBuf[indexForSide] = calcCellLR_<?=side?>(solver, U, dt, x, i, n);
	}<? end ?>
}

