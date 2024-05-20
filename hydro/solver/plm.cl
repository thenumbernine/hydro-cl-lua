//// MODULE_NAME: <?=calcLR?>

// TODO incorporate parallel propagators

<? 
if not solver.usePLM then return end

for side=0,solver.dim-1 do
	
	if solver.usePLM == "piecewise constant" then 
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
	elseif solver.usePLM == "plm cons" then 
?>

// https://en.wikipedia.org/wiki/MUSCL_scheme
// linear extrapolate slopes to the edges.  apply limiters.
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
	real const dx = solver->grid_dx.s<?=side?>;

	real3 xL = cell->pos; xL.s<?=side?> -= dx;
	real3 xR = cell->pos; xR.s<?=side?> += dx;

	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;

	result->L = *U;
	result->R = *U;
	for (int j = 0; j < numIntStates; ++j) {
#if 1	// taken from above / from Dullemond "Hydrodynamics II"
		// oscillates too much
		real const dUR = UR->ptr[j] - U->ptr[j];
		real const dUL = U->ptr[j] - UL->ptr[j];
		// TODO base on wavespeed at cell? 
		// if lambda >= 0 use L/R, (or L/C?) if lambda < 0 use R/L (or R/C) ?
		// but if you want to use wavespeeds then you want to use FVS ?
		// maybe I can do like HLL and just take the leading wavespeed for the eqn?
		real sigma;
		real const dUC = .5 * (dUR - dUL);
		if (dUC >= 0) {	// is it good to do? idk.
			real const r = dUR == 0 ? 0 : (dUL / dUR);
			real const phi = <?=slopeLimiter?>(r);	//works good with minmod, bad with superbee
			sigma = phi * dUR;
		} else {
			real const r = dUL == 0 ? 0 : (dUR / dUL);
			real const phi = <?=slopeLimiter?>(r);
			sigma = phi * dUL;
		}
		result->R.ptr[j] += .5 * sigma;	//uL_{i+1/2} ... right-side of the cell = left-side of the right edge
		result->L.ptr[j] -= .5 * sigma;	//uR_{i-1/2} ... left-side of the cell = right-side of the left edge
#endif
#if 0 	// minmod
		// Based on Mara's PLM.  Looks to be asymmetric, favoring the x+ and y+ directions. Maybe I copied it wrong.
		// also looks like other MUSCL schemes I've looked at that use minmod() limiter
		// even worse than above ...
		real const dUR = UR->ptr[j] - U->ptr[j];
		real const dUL = U->ptr[j] - UL->ptr[j];
		real const dUC = .5 * (dUR + dUL);
		real const adwl = fabs(dUL);
		real const adwr = fabs(dUR);
		real const adwc = fabs(dUC);
		real const mindw = min3(adwl, adwr, adwc);
		real const dUM = .25 * fabs(sign(dUL) + sign(dUC)) * fabs(sign(dUL) + sign(dUR)) * mindw;
		result->L.ptr[j] -= .5 * dUM; 
		result->R.ptr[j] += .5 * dUM; 
#endif
	}
}


<?
	elseif solver.usePLM == "plm cons with flux" then 
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
	
	<?=cons_t?> UHalfL = *U;
	<?=cons_t?> UHalfR = *U;
	for (int j = 0; j < numIntStates; ++j) {
		real const dUL = U->ptr[j] - UL->ptr[j];
		real const dUR = UR->ptr[j] - U->ptr[j];
		
		//Hydrodynamics II slope-limiters (4.4.2) and MUSCL-Hancock (6.6)	
		//https://en.wikipedia.org/wiki/MUSCL_scheme
#if 1
		real const r = dUR == 0 ? 0 : (dUL / dUR);
		// wait, how come I get the feeling that defining limiters in terms of the ratio is something they only and always do in textbooks but never in code ... ???
		real const phi = <?=slopeLimiter?>(r);	//works good with minmod, bad with superbee
		real const sigma = phi * dUR;
#else	//isn't as accurate, or is it?
		real const sigma = minmod(
			minmod(
				fabs(dUL), 
				fabs(dUR)
			),
			fabs(.5 * (dUL - dUR))
		) * sign(dUL - dUR);
#endif
		//q^n_i-1/2,R = q^n_i - 1/2 dx sigma	(Hydrodynamics II 6.58)
		// cell i's left edge's right state
		UHalfL.ptr[j] -= .5 * sigma;
	
		//q^n_i+1/2,L = q^n_i + 1/2 dx sigma	(Hydrodynamics II 6.59)
		// cell i's right edge's left state
		UHalfR.ptr[j] += .5 * sigma;
	}

//// MODULE_DEPENDS: <?=cell_dx_i?>
	real const dx = cell_dx<?=side?>(cell->pos);
	real const dt_dx = dt / dx;

<? print("TODO cell_t averaging here?") ?>
	<?=cons_t?> FHalfL;
	<?=fluxFromCons?>(&FHalfL, solver, &UHalfL, cellL, normal_forSide<?=side?>(xIntL));
	<?=cons_t?> FHalfR;
	<?=fluxFromCons?>(&FHalfR, solver, &UHalfR, cellR, normal_forSide<?=side?>(xIntR));

	result->L = UHalfL;
	result->R = UHalfR;
	for (int j = 0; j < numIntStates; ++j) {
		real const dF = FHalfR.ptr[j] - FHalfL.ptr[j];

		//U-cell-L = q^n+1/2_i-1/2,R (Hydrodynamics II 6.62)
		// = q^n_i-1/2,R + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
		result->L.ptr[j] += .5 * dt_dx * dF;

		//U-cell-R = q^n+1/2_i+1/2,L (Hydrodynamics II 6.63)
		// = q^n_i+1/2,L + 1/2 dt/dx (f_k(q^n_i+1/2,L) - f_k(q_i-1/2,R))
		result->R.ptr[j] += .5 * dt_dx * dF;
	}
}

<? 
	elseif solver.usePLM == "plm prim" then 
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
		real const dWR = WR.ptr[j] - W.ptr[j];
		real const dWL = W.ptr[j] - WL.ptr[j];
		real const r = dWR == 0 ? 0 : (dWL / dWR);
		real const phi = <?=slopeLimiter?>(r);
		real const sigma = phi * dWR;
		nWL.ptr[j] -= .5 * sigma;
		nWR.ptr[j] += .5 * sigma;
#endif
#if 0	// Based on Mara's PLM.  Looks to be asymmetric, favoring the x+ and y+ directions. Maybe I copied it wrong.
		// muuuch weaker than above, hmm
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
	elseif solver.usePLM == "plm eig" then 
?>

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
#if 1



	// this is me trying to transform "plm cons with flux" into "plm eig" and see where things are going wrong ...



//// MODULE_DEPENDS: <?=cell_dx_i?>
	global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
	
	global <?=cell_t?> const * const cellL = cell - solver->stepsize.s<?=side?>;
	global <?=cell_t?> const * const cellR = cell + solver->stepsize.s<?=side?>;

	real3 xIntL = cell->pos; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real3 xIntR = cell->pos; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

	<?=cons_t?> dUL;
	<?=cons_t?> dUR;
	<?=cons_t?> dUC;
	for (int j = 0; j < numIntStates; ++j) {
		dUL.ptr[j] = U->ptr[j] - UL->ptr[j];
		dUR.ptr[j] = UR->ptr[j] - U->ptr[j];
		dUC.ptr[j] = .5 * (UR->ptr[j] - UL->ptr[j]);
	}
	for (int j = numIntStates; j < numStates; ++j) {
		dUL.ptr[j] = dUR.ptr[j] = dUC.ptr[j] = 0;
	}

	//calc eigen values and vectors at cell center
	<?=eigen_t?> eig;
	<?=eigen_forCell?>(&eig, solver, U, cell, n);

	// transform each into wave-space of the cell-centered flux
	<?=waves_t?> dULEig;
	<?=eigen_leftTransform?>(&dULEig, solver, &eig, &dUL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=waves_t?> dUREig;
	<?=eigen_leftTransform?>(&dUREig, solver, &eig, &dUR, xIntR, normal_forSide<?=side?>(xIntR));
	<?=waves_t?> dUCEig;
	<?=eigen_leftTransform?>(&dUCEig, solver, &eig, &dUC, cell->pos, n);

	// prefix-code for the cell-centered waves
	<?=eqn:eigenWaveCodePrefix{
		n = "n",
		eig = "&eig",
		pt = "cell->pos",
	}:gsub("\n", "\n\t")?>

	// now per-wave, zero the components that are going left or right depeneding on the wavespeed
	//MUSCL slope of characteristic variables
	<? for j=0,eqn.numWaves-1 do ?>{
		int const j = <?=j?>;
		real const wave_j = <?=eqn:eigenWaveCode{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
			waveIndex = j,
		}:gsub("\n", "\n\t\t")?>;

#if 1	// how about slope-limiter on char space vars?
		real const dULEig_j = dULEig.ptr[j];
		real const dUREig_j = dUREig.ptr[j];
		real const rEig = dUREig_j == 0 ? 0 : (dULEig_j / dUREig_j);
		real const phi = <?=slopeLimiter?>(rEig);
		real const sigma = phi * dUREig_j;
		dULEig.ptr[j] = sigma;
		dUREig.ptr[j] = sigma;
#endif

#if 0
		real const dULEig_j = dULEig.ptr[j];
		real const dUREig_j = dUREig.ptr[j];
		real const dUCEig_j = dUCEig.ptr[j];
		real const rEig = dUCEig_j == 0
			? 0
			: (
				wave_j >= 0
				? dULEig_j / dUCEig_j
				: dUREig_j / dUCEig_j
			);
		real const phi = <?=slopeLimiter?>(rEig);
		real const sigma = phi * dUCEig_j;
		// limit it in both directions?  if so then why store separate variables?
		dULEig.ptr[j] = sigma;
		dUREig.ptr[j] = sigma;
#endif

#if 1	// FVS ... on PLM ... 
		if (wave_j >= 0) {
			dUREig.ptr[j] = 0;
		}
		if (wave_j <= 0) {
			dULEig.ptr[j] = 0;
		}
#endif
	}<? end ?>


	//convert back to conservation variable space
	<?=cons_t?> sL;
	<?=eigen_rightTransform?>(&sL, solver, &eig, &dULEig, cell->pos, n);
	<?=cons_t?> sR;
	<?=eigen_rightTransform?>(&sR, solver, &eig, &dUREig, cell->pos, n);

	<?=cons_t?> UHalfL = *U;
	<?=cons_t?> UHalfR = *U;
	for (int j = 0; j < numIntStates; ++j) {
#if 0	//slope-limiter on delta cons vars		
		real const dUL_j = sL.ptr[j];
		real const dUR_j = sR.ptr[j];
		real const r = dUR_j == 0 ? 0 : (dUL_j / dUR_j);
		real const phi = <?=slopeLimiter?>(r);
		real const sigma = phi * dUR_j;
		UHalfL.ptr[j] -= .5 * sigma;
		UHalfR.ptr[j] += .5 * sigma;
#else	//no slope limiter
		// UHalfL = U[i] - 1/2 EigR limit(EigL (U[i] - U[i-1]))
		// = EigR EigL U[i] - 1/2 EigR limit(EigL (U[i] - U[i-1]))
		// = EigR (EigL U[i] - 1/2 limit(EigL (U[i] - U[i-1])))
		// ~= EigR (EigL U[i] - 1/2 EigL U[i] + 1/2 EigL U[i-1])
		// = 1/2 EigR (EigL U[i] + EigL U[i-1])
		UHalfL.ptr[j] -= .5 * sL.ptr[j];
		UHalfR.ptr[j] += .5 * sR.ptr[j];
#endif
	}

	real const dx = cell_dx<?=side?>(cell->pos);
	real const dt_dx = dt / dx;

	<?=cons_t?> FHalfL;
	<?=fluxFromCons?>(&FHalfL, solver, &UHalfL, cellL, normal_forSide<?=side?>(xIntL));
	<?=cons_t?> FHalfR;
	<?=fluxFromCons?>(&FHalfR, solver, &UHalfR, cellR, normal_forSide<?=side?>(xIntR));

	//UHalfL = 
	result->L = UHalfL;
	result->R = UHalfR;
	for (int j = 0; j < numIntStates; ++j) {
		real const dF = FHalfR.ptr[j] - FHalfL.ptr[j];
		result->L.ptr[j] += .5 * dt_dx * dF;
		result->R.ptr[j] += .5 * dt_dx * dF;
	}



#else
//// MODULE_DEPENDS: <?=cell_dx_i?>
	real const dx = cell_dx<?=side?>(cell->pos);
	real const dt_dx = dt / dx;

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
	
	// transform each into wave-space of the cell-centered flux
	<?=waves_t?> dULEig;
	<?=eigen_leftTransform?>(&dULEig, solver, &eig, &dUL, xIntL, normal_forSide<?=side?>(xIntL));
	<?=waves_t?> dUREig;
	<?=eigen_leftTransform?>(&dUREig, solver, &eig, &dUR, xIntR, normal_forSide<?=side?>(xIntR));
	<?=waves_t?> dUCEig;
	<?=eigen_leftTransform?>(&dUCEig, solver, &eig, &dUC, cell->pos, n);

	// prefix-code for the cell-centered waves
	<?=eqn:eigenWaveCodePrefix{
		n = "n",
		eig = "&eig",
		pt = "cell->pos",
	}:gsub("\n", "\n\t")?>
	
	// slopes in characteristic space
	<?=waves_t?> aL, aR;

	//MUSCL slope of characteristic variables
	//TODO make this based on the choice of slope limiter
	<? for j=0,eqn.numWaves-1 do ?>{
		int const j = <?=j?>;
		real const wave_j = <?=eqn:eigenWaveCode{
			n = "n",
			eig = "&eig",
			pt = "cell->pos",
			waveIndex = j,
		}:gsub("\n", "\n\t\t")?>;

#if 1	// This is adapted from above to use the modular slopeLimiter
		real const r = dUCEig.ptr[j] == 0 
			? 0
			: (wave_j >= 0
				? (dULEig.ptr[j] / dUCEig.ptr[j])
				: (dUREig.ptr[j] / dUCEig.ptr[j])
			);
		real const phi = <?=slopeLimiter?>(r);
		real const sigma = phi * dUCEig.ptr[j];
		real const dUMEig_j = sigma;
#endif

//This one is symmetric, unlike the above center-slope-based methods.
// It also takes into account the center slope, unlike the <?=slopeLimiter?>() options.
// However it is fixed to minmod and does not use the modular <?=slopeLimiter?>() option.
#if 0	
		//real const dUMEig_j = minmod(minmod(2. * dULEig.ptr[j], 2. * dUREig.ptr[j]), dUCEig.ptr[j]);
		real const dUMEig_j = dULEig.ptr[j] * dUREig.ptr[j] < 0 ? 0 : (
			(dUCEig.ptr[j] >= 0. ? 1. : -1.)
			* 2. * min3(
				fabs(dULEig.ptr[j]),
				fabs(dUREig.ptr[j]),
				fabs(dUCEig.ptr[j]))
		);
#endif	

		aL.ptr[j] = wave_j < 0 ? 0 : dUMEig_j * .5 * (1. - wave_j * dt_dx);
		aR.ptr[j] = wave_j > 0 ? 0 : dUMEig_j * .5 * (1. + wave_j * dt_dx);
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
#endif
}

<? 
	elseif solver.usePLM == "plm eig prim" 
	or solver.usePLM == "plm eig prim ref" 
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

//// MODULE_DEPENDS: <?=cell_dx_i?>
	real dx = cell_dx<?=side?>(cell->pos);
	real dt_dx = dt / dx;

	<?=eqn:eigenWaveCodePrefix{
		n = "n",
		eig = "&eig",
		pt = "cell->pos",
	}:gsub("\n", "\n\t")?>

<? 	
		if solver.usePLM == "plm eig prim" then 
?>

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
		elseif solver.usePLM == "plm eig prim ref" then 
?>

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

	elseif solver.usePLM == "plm athena" then 

?>

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

//// MODULE_DEPENDS: <?=cell_dx_i?>
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
	elseif solver.usePLM == 'ppm' then 
		assert(solver.numGhost >= 2)
?>

// https://en.wikipedia.org/wiki/MUSCL_scheme#Piecewise_parabolic_reconstruction
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
	global <?=cons_t?> const * const UL2 = UL - solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const UR2 = UR + solver->stepsize.s<?=side?>;

	<?=cons_t?> dUL2, dUL, dUR, dUR2;
	for (int j = 0; j < numIntStates; ++j) {
		dUL2.ptr[j] = UL->ptr[j] - UL2->ptr[j];	// δu_{i-3/2}
		dUL.ptr[j] = U->ptr[j] - UL->ptr[j];	// δu_{i-1/2}
		dUR.ptr[j] = UR->ptr[j] - U->ptr[j];	// δu_{i+1/2}
		dUR2.ptr[j] = UR2->ptr[j] - UR->ptr[j];	// δu_{i+3/2}
	}
	for (int j = numIntStates; j < numStates; ++j) {
		dUL2.ptr[j] = 0;
		dUL.ptr[j] = 0;
		dUR.ptr[j] = 0;
		dUR2.ptr[j] = 0;
	}

	real const kappa = 1./3.;

	<?=cons_t?> U_edgeR_L = *U;		// uL_{i+1/2}
	//<?=cons_t?> U_edgeR_R = *UR;	// uR_{i+1/2}
	//<?=cons_t?> U_edgeL_L = *UL;	// uL_{i-1/2}
	<?=cons_t?> U_edgeL_R = *U;		// uR_{i-1/2}
	for (int j = 0; j < numIntStates; ++j) {
		real const r = dUR.ptr[j] == 0 ? 0 : dUL.ptr[j] / dUR.ptr[j];
		real const phi = <?=slopeLimiter?>(r);
		//real const rL = dUL.ptr[j] == 0 ? 0 : dUL2.ptr[j] / dUL.ptr[j];
		//real const phiL = <?=slopeLimiter?>(rL);
		//real const rR = dUR2.ptr[j] == 0 ? 0 : dUR.ptr[j] / dUR2.ptr[j];
		//real const phiR = <?=slopeLimiter?>(rR);
		U_edgeR_L.ptr[j] += .25 * phi * ((1. - kappa) * dUL.ptr[j] + (1. + kappa) * dUR.ptr[j]);
		//U_edgeR_R.ptr[j] -= .25 * phiR * ((1. - kappa) * dUR2.ptr[j] + (1. + kappa) * dUR.ptr[j]);
		//U_edgeL_L.ptr[j] += .25 * phiL * ((1. - kappa) * dUL2.ptr[j] + (1. + kappa) * dUL.ptr[j]);
		U_edgeL_R.ptr[j] -= .25 * phi * ((1. - kappa) * dUR.ptr[j] + (1. + kappa) * dUL.ptr[j]);
	}

// TODO I think I have to store the flux on the left and right per-cell ... and don't assume they will match the neighbor's edge flux
// or why not just make it a distinct gridsolver?
// and if I did, I wouldn't even be using 'calcFluxKernelObj' because calcFluxKernelObj assumes the flux at an interface is matching for the cells on either side (is this the case for PPM?)
#if 0
	result->L = interfaceState(U_edgeL_L, U_edgeL_R);
	result->R = interfaceState(U_edgeR_L, U_edgeR_R);
#else
	result->L = U_edgeL_R;
	result->R = U_edgeR_L;
#endif

// TODO too much?
}

<?

	elseif solver.usePLM == 'plm athena' then 

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
	else
		error("got unknown solver.usePLM = "..require 'ext.tolua'(solver.usePLM))
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

	//TODO skip this lr stuff if we're doing 'piecewise constant'
	//...and just use the original buffers
	<? for side=0,solver.dim-1 do ?>{
		<?=normal_t?> const n = normal_forSide<?=side?>(cell->pos);
		
		//cell-centered index for a particular side...
		int const indexForSide = <?=side?> + dim * index;
		calcCellLR_<?=side?>(ULRBuf + indexForSide, solver, U, cell, dt, i, n);
	}<? end ?>
}

