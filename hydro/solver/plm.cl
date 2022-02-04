//// MODULE_NAME: <?=calcLR?>
//// MODULE_DEPENDS: <?=consLR_t?> <?=solver_t?> <?=cons_t?> <?=normal_t?> <?=cell_dx_i?>

// TODO incorporate parallel propagators

<? 
if not solver.usePLM then return end

for side=0,solver.dim-1 do
	
	-- piecewise-constant
	if solver.usePLM == "piecewise-constant" then 
?>

#define calcCellLR_<?=side?>(\
	/*<?=consLR_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell,\
	/*real const */dt,\
	/*int4 const */i,\
	/*<?=normal_t?> const */n\
) {\
	(result)->L = *(U);\
	(result)->R = *(U);\
}

<?
	-- piecewise-linear
	elseif solver.usePLM == "plm-cons" then 
?>

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
void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
	
	global <?=cell_t?> const * const cellL = cell - solver->stepsize.s<?=side?>;
	global <?=cell_t?> const * const cellR = cell + solver->stepsize.s<?=side?>;

	real3 xIntL = cell->pos; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = cell->pos; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
	<?=cons_t?> UHalfL, UHalfR;
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
//// MODULE_DEPENDS: <?=slopeLimiter?>
		real phi = <?=slopeLimiter?>(r);	//works good with minmod, bad with superbee
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

	real dx = cell_dx<?=side?>(cell->pos);
	real dt_dx = dt / dx;

//// MODULE_DEPENDS: <?=fluxFromCons?> 
<? print("TODO cell_t averaging here?") ?>
	<?=cons_t?> FHalfL;
	<?=fluxFromCons?>(&FHalfL, solver, &UHalfL, cellL, normal_forSide<?=side?>(xIntL));
	<?=cons_t?> FHalfR;
	<?=fluxFromCons?>(&FHalfR, solver, &UHalfR, cellR, normal_forSide<?=side?>(xIntR));

	result->L = *U;
	result->R = *U;
	for (int j = 0; j < numIntStates; ++j) {
		real dF = FHalfR.ptr[j] - FHalfL.ptr[j];

		//U-cell-L = q^n+1/2_i-1/2,R (Hydrodynamics II 6.62)
		// = q^n_i-1/2,R + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
		result->L.ptr[j] = UHalfL.ptr[j] + .5 * dt_dx * dF;

		//U-cell R = q^n+1/2_i+1/2,L (Hydrodynamics II 6.63)
		// = q^n_i+1/2,L + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
		result->R.ptr[j] = UHalfR.ptr[j] + .5 * dt_dx * dF;
	}
}

<? 
	-- technically i should call the one above "plm-cons-with-flux" 
	-- and this one "plm-cons"
	elseif solver.usePLM == "plm-cons-alone" then 
?>

void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	// extrapolate slopes in consered variable space
	real dx = solver->grid_dx.s<?=side?>;

	real3 xL = cell->pos; xL.s<?=side?> -= dx;
	real3 xR = cell->pos; xR.s<?=side?> += dx;

	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;

	result->L = *U;
	result->R = *U;
	// minmod
	for (int j = 0; j < numIntStates; ++j) {
#if 1	// taken from above / from Dullemond "Hydrodynamics II"
		real dUR = UR->ptr[j] - U->ptr[j];
		real dUL = U->ptr[j] - UL->ptr[j];
		real r = dUR == 0 ? 0 : (dUL / dUR);
//// MODULE_DEPENDS: <?=slopeLimiter?>
		real phi = <?=slopeLimiter?>(r);	//works good with minmod, bad with superbee
		real sigma = phi * dUR;
		result->L.ptr[j] -= .5 * sigma;
		result->R.ptr[j] += .5 * sigma;
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
		result->L.ptr[j] -= .5 * dUM; 
		result->R.ptr[j] += .5 * dUM; 
#endif
	}
}

<? 
	elseif solver.usePLM == "plm-prim-alone" then 
?>
//// MODULE_DEPENDS: <?=consFromPrim?> <?=primFromCons?>

void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	// extrapolate slopes in primitive space
	real dx = solver->grid_dx.s<?=side?>;

	real3 xL = cell->pos; xL.s<?=side?> -= dx;
	real3 xR = cell->pos; xR.s<?=side?> += dx;

	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;

	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, cell->pos);
	<?=prim_t?> WL;
	<?=primFromCons?>(&WL, solver, UL, xL);
	<?=prim_t?> WR;
	<?=primFromCons?>(&WR, solver, UR, xR);
	
	// minmod
	<?=prim_t?> nWL = W;
	<?=prim_t?> nWR = W;
	for (int j = 0; j < numIntStates; ++j) {
#if 1	// taken from above / from Dullemond "Hydrodynamics II"
		real dWR = WR.ptr[j] - W.ptr[j];
		real dWL = W.ptr[j] - WL.ptr[j];
		real r = dWR == 0 ? 0 : (dWL / dWR);
//// MODULE_DEPENDS: <?=slopeLimiter?>
		real phi = <?=slopeLimiter?>(r);
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
	
	real3 xIntL = cell->pos; xIntL.s<?=side?> -= .5 * dx;
	real3 xIntR = cell->pos; xIntR.s<?=side?> += .5 * dx;

	//converting from cons->prim and then prim->cons might lose us accuracy? 
	<?=consFromPrim?>(&result->L, solver, &nWL, xIntL);
	<?=consFromPrim?>(&result->R, solver, &nWR, xIntR);
}

<? 
	elseif solver.usePLM == "plm-eig" then 
?>
//// MODULE_DEPENDS: <?=eigen_forCell?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> min3

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

void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	//1) calc delta q's ... l r c (eqn 36)
	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
	<?=cons_t?> dUL, dUR, dUC;
	for (int j = 0; j < numIntStates; ++j) {
		dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
		dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		dUC.ptr[j] = .5 * (UR->ptr[j] - UL->ptr[j]);
	}
	for (int j = numIntStates; j < numStates; ++j) {
		dUL.ptr[j] = dUR.ptr[j] = dUC.ptr[j] = 0;
	}

	real3 xIntL = cell->pos; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = cell->pos; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

	//calc eigen values and vectors at cell center
	<?=eigen_t?> eig;
	<?=eigen_forCell?>(&eig, solver, U, cell, n);
		
	<?=waves_t?> dULEig;
	<?=eigen_leftTransform?>(&dULEig, solver, &eig, &dUL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=waves_t?> dUREig;
	<?=eigen_leftTransform?>(&dUREig, solver, &eig, &dUR, xIntR, normal_forSide<?=side?>(xIntR));
	<?=waves_t?> dUCEig;
	<?=eigen_leftTransform?>(&dUCEig, solver, &eig, &dUC, cell->pos, n);

	//MUSCL slope of characteristic variables
	//TODO make this based on the choice of slope limiter
	<?=waves_t?> dUMEig;
	for (int j = 0; j < numWaves; ++j) {

// This is adapted from above to use the modular slopeLimiter
<? if false then ?>
#if 0
		real r = dUREig.ptr[j] == 0 ? 0 : (dULEig.ptr[j] / dUREig.ptr[j]);
//// MODULE_DEPENDS: <?=slopeLimiter?>
		real phi = <?=slopeLimiter?>(r);
		real sigma = phi * dUREig.ptr[j];
		dUMEig.ptr[j] = sigma;

#endif
<? end ?>

//This one is symmetric, unlike the above center-slope-based methods.
// It also takes into account the center slope, unlike the <?=slopeLimiter?>() options.
// However it is fixed to minmod and does not use the modular <?=slopeLimiter?>() option.
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

	real dx = cell_dx<?=side?>(cell->pos);
	real dt_dx = dt / dx;

	<?=eqn:eigenWaveCodePrefix{
		n = "n",
		eig = "&eig",
		pt = "cell->pos",
	}:gsub("\n", "\n\t")?>

	// slopes in characteristic space
	<?=waves_t?> aL;
	<?=waves_t?> aR;
	<? for j=0,eqn.numWaves-1 do ?>{
		int const j = <?=j?>;
		real const wave_j = <?=eqn:eigenWaveCode{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
			waveIndex = j,
		}:gsub("\n", "\n\t\t")?>;
		aL.ptr[j] = wave_j < 0 ? 0 : dUMEig.ptr[j] * .5 * (1. - wave_j * dt_dx);
		aR.ptr[j] = wave_j > 0 ? 0 : dUMEig.ptr[j] * .5 * (1. + wave_j * dt_dx);
	}<? end ?>
	
	//convert back to conservation variable space
	<?=cons_t?> sL;
	<?=eigen_rightTransform?>(&sL, solver, &eig, &aL, cell->pos, n);
	<?=cons_t?> sR;
	<?=eigen_rightTransform?>(&sR, solver, &eig, &aR, cell->pos, n);

	result->L = *U;
	result->R = *U;
	for (int j = 0; j < numIntStates; ++j) {
		result->L.ptr[j] += sL.ptr[j];
		result->R.ptr[j] -= sR.ptr[j];
	}
}

<? 
	elseif solver.usePLM == "plm-eig-prim" 
	or solver.usePLM == "plm-eig-prim-ref" 
	then 
?>
//// MODULE_DEPENDS: <?=eigen_forCell?> min3 <?=primFromCons?>

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

void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	real3 xL = cell->pos; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
	real3 xR = cell->pos; xR.s<?=side?> += solver->grid_dx.s<?=side?>;
	
	real3 xIntL = cell->pos; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = cell->pos; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
	
	//1) calc delta q's ... l r c (eqn 36)
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, cell->pos);
	<?=prim_t?> WL;
	<?=primFromCons?>(&WL, solver, UL, xL);
	<?=prim_t?> WR;
	<?=primFromCons?>(&WR, solver, UR, xR);
	<?=prim_t?> dWL, dWR, dWC;
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
	<?=eigen_t?> eig;
	<?=eigen_forCell?>(&eig, solver, U, cell, n);
		
	//apply dU/dW before applying left/right eigenvectors so the eigenvectors are of the flux wrt primitives 
	//RW = dW/dU RU, LW = LU dU/dW
	<?=cons_t?> tmp;

	<?=apply_dU_dW?>(&tmp, solver, &W, &dWL, xIntL); 
	<?=waves_t?> dWLEig;
	<?=eigen_leftTransform?>(&dWLEig, solver, &eig, &tmp, xIntL, normal_forSide<?=side?>(xIntL));
	
	<?=apply_dU_dW?>(&tmp, solver, &W, &dWR, xIntR); 
	<?=waves_t?> dWREig;
	<?=eigen_leftTransform?>(&dWREig, solver, &eig, &tmp, xIntR, normal_forSide<?=side?>(xIntR));
	
	<?=apply_dU_dW?>(&tmp, solver, &W, &dWC, cell->pos); 
	<?=waves_t?> dWCEig;
	<?=eigen_leftTransform?>(&dWCEig, solver, &eig, &tmp, cell->pos, n);

	//MUSCL slope of characteristic variables
	<?=waves_t?> dWMEig;
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

	real dx = cell_dx<?=side?>(cell->pos);
	real dt_dx = dt / dx;

	<?=eqn:eigenWaveCodePrefix{
		n = "n",
		eig = "&eig",
		pt = "cell->pos",
	}:gsub("\n", "\n\t")?>

<? 	
		if solver.usePLM == "plm-eig-prim" then 
?>
//// MODULE_DEPENDS: <?=apply_dU_dW?> <?=apply_dW_dU?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> <?=consFromPrim?>

	//without reference state

	// calculate left and right slopes in characteristic space
	<?=waves_t?> aL, aR;
	<? for j=0,eqn.numWaves-1 do ?>{
		int const j = <?=j?>;
		real const wave_j = <?=eqn:eigenWaveCode{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
			waveIndex = j,
		}:gsub("\n", "\n\t\t")?>;
		aL.ptr[j] = wave_j < 0 ? 0 : dWMEig.ptr[j] * .5 * (1. - wave_j * dt_dx);
		aR.ptr[j] = wave_j > 0 ? 0 : dWMEig.ptr[j] * .5 * (1. + wave_j * dt_dx);
	}<? end ?>

	// transform slopes back to conserved variable space
	<?=eigen_rightTransform?>(&tmp, solver, &eig, &aL, xIntL, normal_forSide<?=side?>(xIntL)); 
	<?=prim_t?> sL;
	<?=apply_dW_dU?>(&sL, solver, &W, &tmp, xIntL);
	
	<?=eigen_rightTransform?>(&tmp, solver, &eig, &aR, xIntR, normal_forSide<?=side?>(xIntR)); 
	<?=prim_t?> sR;
	<?=apply_dW_dU?>(&sR, solver, &W, &tmp, xIntR);

	// linearly extrapolate the slopes forward and backward from the cell center
	<?=prim_t?> W2L, W2R;
	for (int j = 0; j < numIntStates; ++j) {
		W2L.ptr[j] = W.ptr[j] + sL.ptr[j];
		W2R.ptr[j] = W.ptr[j] - sR.ptr[j];
	}
	for (int j = numIntStates; j < numStates; ++j) {
		W2L.ptr[j] = W2R.ptr[j] = W.ptr[j];
	}
	<?=consFromPrim?>(&result->L, solver, &W2L, xIntL);
	<?=consFromPrim?>(&result->R, solver, &W2R, xIntR);
<?	
		elseif solver.usePLM == "plm-eig-prim-ref" then 
?>
//// MODULE_DEPENDS: <?=apply_dU_dW?> <?=apply_dW_dU?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> <?=consFromPrim?>

	//with reference state

	//min and max waves
	//TODO use eigenWaveCodeMinMax ... 
	// but the same vars have been declared in eigenWaveCodePrefix above ...
	// so which is worse? 5 min() calcs or 1 sqrt() ?
	real waveMin, waveMax;
	{
		<?=eqn:eigenWaveCodeMinMax{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
			resultMin = "waveMin",
			resultMax = "waveMax",
		}:gsub("\n", "\n\t\t")?>
	}
	waveMin = min(0., waveMin);
	waveMax = max(0., waveMax);

	//limited slope in primitive variable space
	<?=eigen_rightTransform?>(&tmp, solver, &eig, &dWMEig, cell->pos, n);
	<?=prim_t?> dWM;
	<?=apply_dW_dU?>(&dWM, solver, &W, &tmp, cell->pos);

	//left and right reference states
	<?=prim_t?> WLRef, WRRef;
	for (int j = 0; j < numIntStates; ++j) {
		WLRef.ptr[j] = W.ptr[j] + .5 * (1. - dt_dx * waveMax) * dWM.ptr[j];
		WRRef.ptr[j] = W.ptr[j] - .5 * (1. + dt_dx * waveMin) * dWM.ptr[j];
	}
	for (int j = numIntStates; j < numStates; ++j) {
		WLRef.ptr[j] = WRRef.ptr[j] = 0;
	}

	// calculate left and right slopes in characteristic space
	<?=waves_t?> aL, aR;
	<? for j=0,eqn.numWaves-1 do ?>{
		int const j = <?=j?>;
		real const wave_j = <?=eqn:eigenWaveCode{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
			waveIndex = j,
		}:gsub("\n", "\n\t\t")?>;
		aL.ptr[j] = wave_j < 0 ? 0 : (dWMEig.ptr[j] * dt_dx * (waveMax - wave_j));
		aR.ptr[j] = wave_j > 0 ? 0 : (dWMEig.ptr[j] * dt_dx * (waveMin - wave_j));
	}<? end ?>

	// transform slopes back to conserved variable space
	<?=eigen_rightTransform?>(&tmp, solver, &eig, &aL, xIntL, normal_forSide<?=side?>(xIntL)); 
	<?=prim_t?> sL;
	<?=apply_dW_dU?>(&sL, solver, &W, &tmp, xIntL);
	
	<?=eigen_rightTransform?>(&tmp, solver, &eig, &aR, xIntR, normal_forSide<?=side?>(xIntR)); 
	<?=prim_t?> sR;
	<?=apply_dW_dU?>(&sR, solver, &W, &tmp, xIntR);

	// linearly extrapolate the slopes forward and backward from the cell center
	<?=prim_t?> W2L, W2R;
	for (int j = 0; j < numIntStates; ++j) {
		W2R.ptr[j] = WRRef.ptr[j] + .5 * sR.ptr[j];
		W2L.ptr[j] = WLRef.ptr[j] + .5 * sL.ptr[j];
	}
	for (int j = numIntStates; j < numStates; ++j) {
		W2R.ptr[j] = W2L.ptr[j] = W.ptr[j];
	}
	//TODO fix the x's
	<?=consFromPrim?>(&result->L, solver, &W2L, xIntL);
	<?=consFromPrim?>(&result->R, solver, &W2R, xIntR);
<?
		end	-- solver.usePLM
?>
}

<?

	elseif solver.usePLM == "plm-athena" then 

?>
//// MODULE_DEPENDS: <?=eigen_forCell?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> <?=consFromPrim?> <?=primFromCons?>

//based on Athena
void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	real3 xIntL = cell->pos; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = cell->pos; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
	real3 xL = cell->pos; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
	real3 xR = cell->pos; xR.s<?=side?> += solver->grid_dx.s<?=side?>;

	//calc eigen values and vectors at cell center
	<?=eigen_t?> eig;
	<?=eigen_forCell?>(&eig, solver, U, cell, n);

	real const dx = cell_dx<?=side?>(cell->pos);
	real const dt_dx = dt / dx;

	//1) calc delta q's ... l r c (eqn 36)
	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;

	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, cell->pos);
	<?=prim_t?> WL;
	<?=primFromCons?>(&WL, solver, UL, xL);
	<?=prim_t?> WR;
	<?=primFromCons?>(&WR, solver, UR, xR);
	
	<?=prim_t?> dWL, dWR, dWC, dWG;
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
	<?=waves_t?> dal;
	<?=eigen_leftTransform?>(&dal, solver, &eig, (<?=cons_t?>*)&dWL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=waves_t?> dar;
	<?=eigen_leftTransform?>(&dar, solver, &eig, (<?=cons_t?>*)&dWR, xIntR, normal_forSide<?=side?>(xIntR));
	<?=waves_t?> dac;
	<?=eigen_leftTransform?>(&dac, solver, &eig, (<?=cons_t?>*)&dWC, cell->pos, n);
	<?=waves_t?> dag;
	<?=eigen_leftTransform?>(&dag, solver, &eig, (<?=cons_t?>*)&dWG, cell->pos, n);

	<?=waves_t?> da;
	for (int j = 0; j < numWaves; ++j) {
		da.ptr[j] = 0;
		if (dal.ptr[j] * dar.ptr[j] > 0) {
			real lim_slope1 = min(fabs(dal.ptr[j]), fabs(dar.ptr[j]));
			real lim_slope2 = min(fabs(dac.ptr[j]), fabs(dag.ptr[j]));
			da.ptr[j] = sign(dac.ptr[j]) * min((real)2. * lim_slope1, lim_slope2);
		}
	}

	<?=cons_t?> dWm_tmp;
	<?=eigen_rightTransform?>(&dWm_tmp, solver, &eig, &da, cell->pos, n);
	<?=prim_t?> dWm = *(<?=prim_t?>*)&dWm_tmp;

	<?=prim_t?> Wlv, Wrv;
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

	<?=consFromPrim?>(&result->L, solver, &Wrv, xIntL);
	<?=consFromPrim?>(&result->R, solver, &Wlv, xIntR);
}

<? 
	elseif solver.usePLM == 'ppm-wip' then 
?>
//// MODULE_DEPENDS: <?=eigen_forCell?> <?=eigen_leftTransform?> <?=eigen_rightTransform?> <?=consFromPrim?> <?=primFromCons?> min3 sqr <?=apply_dU_dW?> <?=apply_dW_dU?>

// based on http://zingale.github.io/hydro1d/  ppm code

void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {
	real3 xL = cell->pos; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
	real3 xR = cell->pos; xR.s<?=side?> += solver->grid_dx.s<?=side?>;

	real const z0 = .75;
	real const z1 = .85;
	real const delta = .33;
	real const minP = 1e-10;
	
	<?=prim_t?> W_im3;
	<?=primFromCons?>(&W_im3, solver, U - 3 * solver->stepsize.s<?=side?>, cell->pos);
	<?=prim_t?> W_im2;
	<?=primFromCons?>(&W_im2, solver, U - 2 * solver->stepsize.s<?=side?>, cell->pos);
	<?=prim_t?> W_im1;
	<?=primFromCons?>(&W_im1, solver, U - 1 * solver->stepsize.s<?=side?>, cell->pos);
	<?=prim_t?> W_i;
	<?=primFromCons?>(&W_i, solver, U, cell->pos);
	<?=prim_t?> W_ip1;
	<?=primFromCons?>(&W_ip1, solver, U + 1 * solver->stepsize.s<?=side?>, cell->pos);
	<?=prim_t?> W_ip2;
	<?=primFromCons?>(&W_ip2, solver, U + 2 * solver->stepsize.s<?=side?>, cell->pos);
	<?=prim_t?> W_ip3;
	<?=primFromCons?>(&W_ip3, solver, U + 3 * solver->stepsize.s<?=side?>, cell->pos);

	// this is based on pressure ... how to generalize for all PDEs?
	real const dP_i = fabs(W_ip1.P - W_im1.P);
	real const dP2_i = fabs(W_ip2.P - W_im2.P);
	real const z_i = dP_i / max(dP2_i, minP);
	
	real xi_t_i;
	if ((W_im1.v.s<?=side?> - W_ip1.v.s<?=side?> > 0) && (dP_i / min(W_ip1.P, W_im1.P) > delta)) {
		xi_t_i = clamp(1. - (z_i - z0) / (z1 - z0), 0., 1.);
	} else {
		xi_t_i = 0;
	}


	real const dP_ip1 = fabs(W_ip2.P - W_i.P);
	real const dP2_ip1 = fabs(W_ip3.P - W_im1.P);
	real const z_ip1 = dP_ip1 / max(dP2_ip1, minP);
	
	real xi_t_ip1;
	if ((W_i.v.s<?=side?> - W_ip2.v.s<?=side?> > 0) && (dP_ip1 / min(W_ip2.P, W_i.P) > delta)) {
		xi_t_ip1 = clamp(1. - (z_ip1 - z0) / (z1 - z0), 0., 1.);
	} else {
		xi_t_ip1 = 0;
	}


	real const dP_im1 = fabs(W_i.P - W_im2.P);
	real const dP2_im1 = fabs(W_ip1.P - W_im3.P);
	real const z_im1 = dP_im1 / max(dP2_im1, minP);
	
	real xi_t_im1;
	if ((W_im2.v.s<?=side?> - W_i.v.s<?=side?> > 0) && (dP_im1 / min(W_i.P, W_im2.P) > delta)) {
		xi_t_im1 = clamp(1. - (z_im1 - z0) / (z1 - z0), 0., 1.);
	} else {
		xi_t_im1 = 0;
	}


	real xi_i;
	if (W_ip1.P - W_im1.P > 0) {
		xi_i = min(xi_t_i, xi_t_im1);
	} else if (W_ip1.P - W_im1.P < 0) {
		xi_i = min(xi_t_i, xi_t_ip1);
	} else {
		xi_i = xi_t_i;
	}

	<?=cons_t?> Wplus_i;
	<?=cons_t?> Wminus_i;
	for (int j = 0; j < numStates; ++j) {
		// 1984 collela & woodward eqn 1.7
		real dq0_i = .5 * (W_ip1.ptr[j] - W_im1.ptr[j]);
		
		// 1984 collela & woodward eqn 1.8
		if ((W_ip1.ptr[j] - W_i.ptr[j]) * (W_i.ptr[j] - W_im1.ptr[j]) > 0) {
			dq0_i = sign(dq0_i) * min3( 
				fabs(dq0_i),
				2. * fabs(W_i.ptr[j] - W_im1.ptr[j]),
				2. * fabs(W_ip1.ptr[j] - W_i.ptr[j]));
		} else {
			dq0_i = 0;
		}
	
		// 1984 collela & woodward eqn 1.7
		real dqp_i = .5 * (W_ip2.ptr[j] - W_i.ptr[j]);
		
		// 1984 collela & woodward eqn 1.8
		if ((W_ip2.ptr[j] - W_ip1.ptr[j]) * (W_ip1.ptr[j] - W_i.ptr[j]) > 0) {
			dqp_i = sign(dqp_i) * min3(
				fabs(dqp_i),
				2. * fabs(W_ip1.ptr[j] - W_i.ptr[j]),
				2. * fabs(W_ip2.ptr[j] - W_ip1.ptr[j]));
		} else {
			dqp_i = 0;
		}


		// 1984 collela & woodward eqn 1.7
		real dq0_im1 = .5 * (W_i.ptr[j] - W_im2.ptr[j]);
		
		// 1984 collela & woodward eqn 1.8
		if ((W_i.ptr[j] - W_im1.ptr[j]) * (W_im1.ptr[j] - W_im2.ptr[j]) > 0) {
			dq0_im1 = sign(dq0_im1) * min3( 
				fabs(dq0_im1),
				2. * fabs(W_im1.ptr[j] - W_im2.ptr[j]),
				2. * fabs(W_i.ptr[j] - W_im1.ptr[j]));
		} else {
			dq0_im1 = 0;
		}
	
		// 1984 collela & woodward eqn 1.7
		real dqp_im1 = .5 * (W_ip1.ptr[j] - W_im1.ptr[j]);
		
		// 1984 collela & woodward eqn 1.8
		if ((W_ip1.ptr[j] - W_i.ptr[j]) * (W_i.ptr[j] - W_im1.ptr[j]) > 0) {
			dqp_im1 = sign(dqp_im1) * min3(
				fabs(dqp_im1),
				2. * fabs(W_i.ptr[j] - W_im1.ptr[j]),
				2. * fabs(W_ip1.ptr[j] - W_i.ptr[j]));
		} else {
			dqp_im1 = 0;
		}


		// 1984 collela & woodward, eqn 1.6
		Wplus_i.ptr[j] = .5 * (W_i.ptr[j] + W_ip1.ptr[j]) - (1./6.) * (dqp_i - dq0_i);
		Wminus_i.ptr[j] = .5 * (W_im1.ptr[j] + W_i.ptr[j]) - (1./6.) * (dqp_im1 - dq0_im1);
	
		// 2008 Colella & Sekora
		Wplus_i.ptr[j] = max(Wplus_i.ptr[j], min(W_i.ptr[j], W_ip1.ptr[j]));
		Wplus_i.ptr[j] = min(Wplus_i.ptr[j], max(W_i.ptr[j], W_ip1.ptr[j]));
		
		Wminus_i.ptr[j] = max(Wminus_i.ptr[j], min(W_im1.ptr[j], W_i.ptr[j]));
		Wminus_i.ptr[j] = min(Wminus_i.ptr[j], max(W_im1.ptr[j], W_i.ptr[j]));
	}

	// 1984 Collela & Woodward eqn 1.10 is somewhere in here:
	for (int j = 0; j < numStates; ++j) {
		if ((Wplus_i.ptr[j] - W_i.ptr[j]) * (W_i.ptr[j] - Wminus_i.ptr[j]) <= 0) {
			Wminus_i.ptr[j] = W_i.ptr[j];
			Wplus_i.ptr[j] = W_i.ptr[j];
		} else if ((Wplus_i.ptr[j] - Wminus_i.ptr[j]) * (W_i.ptr[j] - .5 * (Wminus_i.ptr[j] + Wplus_i.ptr[j])) > sqr(Wplus_i.ptr[j] - Wminus_i.ptr[j]) / 6.) {
			Wminus_i.ptr[j] = 3. * W_i.ptr[j] - 2. * Wplus_i.ptr[j];
		} else if (-sqr(Wplus_i.ptr[j] - Wminus_i.ptr[j]) / 6. > (Wplus_i.ptr[j] - Wminus_i.ptr[j]) * (W_i.ptr[j] - .5 * (Wminus_i.ptr[j] + Wplus_i.ptr[j]))) {
			Wplus_i.ptr[j] = 3. * W_i.ptr[j] - 2. * Wminus_i.ptr[j];
		}
	}

	// flatten
	for (int j = 0; j < numStates; ++j) {
		Wminus_i.ptr[j] = (1. - xi_i) * W_i.ptr[j] + xi_i * Wminus_i.ptr[j];
		Wplus_i.ptr[j] = (1. - xi_i) * W_i.ptr[j] + xi_i * Wplus_i.ptr[j];
	}
	
	// W6
	<?=cons_t?> W6_i;
	for (int j = 0; j < numStates; ++j) {
		W6_i.ptr[j] = 6. * W_i.ptr[j] - 3. * (Wminus_i.ptr[j] + Wplus_i.ptr[j]);
	}

	<?=eigen_t?> eig;
	<?=eigen_forCell?>(&eig, solver, U, cell, n);

	<?=waves_t?> lambda;
	real lambdaMin = INFINITY;
	real lambdaMax = -INFINITY;
	{
		<?=eqn:eigenWaveCodePrefix{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
		}:gsub("\n", "\n\t\t")?>
		<? for j=0,eqn.numWaves-1 do ?>{
			real const lambda_j = <?=eqn:eigenWaveCode{
				n = "n",
				eig = "&eig",	-- TODO once again this was set to U when the prefix was based on eig ...
				"cell->pos",
				waveIndex = j,
			}:gsub("\n", "\n\t\t\t")?>;
			lambdaMin = min(lambdaMin, lambda_j);
			lambdaMax = min(lambdaMax, lambda_j);
			lambda.ptr[<?=j?>] = lambda_j;
		}<? end ?>	
	}

	real dx = cell_dx<?=side?>(cell->pos);
	real dt_dx = dt / dx;
	
	<?=prim_t?> Iplus[numWaves];
	<?=prim_t?> Iminus[numWaves];
	for (int j = 0; j < numWaves; ++j) {
		real sigma = fabs(lambda.ptr[j]) * dt_dx;
		for (int k = 0; k < numStates; ++k) {
			if (lambda.ptr[j] >= 0) {
				Iplus[j].ptr[k] = Wplus_i.ptr[k] - .5 * sigma * (Wplus_i.ptr[k] - Wminus_i.ptr[k] - (1 - (2./3.)*sigma) * W6_i.ptr[k]);
			} else {
				Iplus[j].ptr[k] = W_i.ptr[k];
			}
			if (lambda.ptr[j] <= 0) {
				Iminus[j].ptr[k] = Wminus_i.ptr[k] + .5 * sigma * (Wplus_i.ptr[k] - Wminus_i.ptr[k] + (1 - (2./3.)*sigma) * W6_i.ptr[k]);
			} else {
				Iminus[j].ptr[k] = W_i.ptr[k];
			}
		}
	}

	<?=prim_t?> Wref_xp;
	if (lambdaMax >= 0) {
		Wref_xp = Iplus[numWaves-1];	//this assumes max == numWaves-1
	} else {
		Wref_xp = W_i;
	}

	<?=prim_t?> Wref_xm;
	if (lambdaMin <= 0) {
		Wref_xm = Iminus[0];			//this assumes min == 0
	} else {
		Wref_xm = W_i;
	}

	// beta± = L*(Wref± - I±)

	// dF = dF/dU dU = R_F Lambda_F L_F dU
	// 

	//TODO ... convert prim to cons, then cons eig transform
	<?=cons_t?> Uref_xp;
	<?=apply_dU_dW?>(&Uref_xp, solver, &W_i, &Wref_xp, cell->pos);
	<?=waves_t?> beta_xp;
	<?=eigen_leftTransform?>(&beta_xp, solver, &eig, &Uref_xp, cell->pos, n);
	
	<?=cons_t?> Uref_xm;
	<?=apply_dU_dW?>(&Uref_xm, solver, &W_i, &Wref_xm, cell->pos);
	<?=waves_t?> beta_xm;
	<?=eigen_leftTransform?>(&beta_xm, solver, &eig, &Uref_xm, cell->pos, n);
	
	for (int m = 0; m < numWaves; ++m) {
		{
			<?=cons_t?> Uminus;
			<?=apply_dU_dW?>(&Uminus, solver, &W_i, Iminus + m, cell->pos);
			<?=waves_t?> waves;
			<?=eigen_leftTransform?>(&waves, solver, &eig, &Uminus, cell->pos, n);
			beta_xm.ptr[m] -= waves.ptr[m];
		}

		{
			<?=cons_t?> Uplus;
			<?=apply_dU_dW?>(&Uplus, solver, &W_i, Iplus + m, cell->pos);
			<?=waves_t?> waves;
			<?=eigen_leftTransform?>(&waves, solver, &eig, &Uplus, cell->pos, n);
			beta_xp.ptr[m] -= waves.ptr[m];
		}
	}

	for (int j = 0; j < numWaves; ++j) {
		if (lambda.ptr[j] < 0) beta_xp.ptr[j] = 0;
		if (lambda.ptr[j] > 0) beta_xm.ptr[j] = 0;
	}
	
	<?=cons_t?> U_l_ip1;
	<?=eigen_rightTransform?>(&U_l_ip1, solver, &eig, &beta_xp, cell->pos, n);
	<?=cons_t?> U_r_i;
	<?=eigen_rightTransform?>(&U_r_i, solver, &eig, &beta_xm, cell->pos, n);

	<?=cons_t?> W_l_ip1;
	<?=apply_dW_dU?>(&W_l_ip1, solver, &W_i, &U_l_ip1, cell->pos);
	<?=cons_t?> W_r_i;
	<?=apply_dW_dU?>(&W_r_i, solver, &W_i, &U_r_i, cell->pos);

	for (int j = 0; j < numStates; ++j) {
		W_l_ip1.ptr[j] = Wref_xp.ptr[j] - W_l_ip1.ptr[j];
		W_r_i.ptr[j] = Wref_xm.ptr[j] - W_r_i.ptr[j];
	}

	<?=consFromPrim?>(&result->L, solver, &W_r_i, cell->pos);
	<?=consFromPrim?>(&result->R, solver, &W_l_ip1, cell->pos);
}

<?

	elseif solver.usePLM == 'plm-athena' then 

-- TODO correlate between calcCellLR's U access and <?=SETBOUNDS?> of <?=calcLR?>()
-- to make sure there's no OOB reads
?>

void calcCellLR_<?=side?>(
	global <?=consLR_t?> * const result,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	real const dt,
	int4 const i,
	<?=normal_t?> const n
) {


}

<?
	end
end
?>

kernel void <?=calcLR?>(
	constant <?=solver_t?> const * const solver,
	global <?=cell_t?> const * const cellBuf,
	global <?=consLR_t?> * const ULRBuf,
	global <?=cons_t?> const * const UBuf,
	real const dt
) {
	<?=SETBOUNDS?>(1,1);
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;

	//TODO skip this lr stuff if we're doing piecewise-constant
	//...and just use the original buffers
	<? for side=0,solver.dim-1 do ?>{
		<?=normal_t?> const n = normal_forSide<?=side?>(cell->pos);
		
		//cell-centered index for a particular side...
		int const indexForSide = <?=side?> + dim * index;
		calcCellLR_<?=side?>(ULRBuf + indexForSide, solver, U, cell, dt, i, n);
	}<? end ?>
}

