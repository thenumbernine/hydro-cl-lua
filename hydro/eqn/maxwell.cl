//// MODULE_NAME: <?=sqrt_2_and_1_2?>

#define sqrt_2 <?=("%.50f"):format(math.sqrt(2))?>
#define sqrt_1_2 <?=("%.50f"):format(math.sqrt(.5))?>

//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: <?=coordLenSq?> <?=cartesianToCoord?> <?=coord_lower?>

<? if scalar == "real" then ?>

#define eqn_coordLenSq coordLenSq
#define eqn_cartesianToCoord cartesianToCoord
#define eqn_coord_lower coord_lower

<? elseif scalar == "cplx" then ?>

real eqn_coordLenSq(cplx3 v, real3 x) {
	return coordLenSq(cplx3_re(v), x)
		+ coordLenSq(cplx3_im(v), x);
}

cplx3 eqn_cartesianToCoord(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		cartesianToCoord(cplx3_re(v), x),
		cartesianToCoord(cplx3_im(v), x));
}

cplx3 eqn_coord_lower(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		coord_lower(cplx3_re(v), x),
		coord_lower(cplx3_im(v), x));
}

<? end -- scalar ?>

#define /*<?=vec3?>*/ calc_E(\
	/*global <?=cons_t?> const * const */U\
)  (<?=vec3?>_<?=susc_t?>_mul((U)->D, <?=susc_t?>_mul((U)->sqrt_1_eps, (U)->sqrt_1_eps)))

#define /*<?=vec3?>*/ calc_H(\
	/*global <?=cons_t?> const * const */U\
) 	(<?=vec3?>_<?=susc_t?>_mul((U)->B, <?=susc_t?>_mul((U)->sqrt_1_mu, (U)->sqrt_1_mu)))

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=eqn_common?> 

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool const lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;

	/* used */
	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=scalar?>_from_real(1.);
	<?=susc_t?> permittivity = <?=susc_t?>_from_real(1.);
	<?=susc_t?> permeability = <?=susc_t?>_from_real(1.);
	<?=scalar?> rhoCharge = <?=zero?>;

	/* throw-away */
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;

	<?=initCode()?>

	U->D = eqn_cartesianToCoord(D, x);
	U->B = eqn_cartesianToCoord(B, x);
	U->phi = <?=zero?>;
	U->psi = <?=zero?>;
	U->J = <?=vec3?>_zero;
	U->rhoCharge = rhoCharge;
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permittivity));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permeability));
}


//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=normal_t?> <?=cons_t?> <?=prim_t?> <?=eqn_common?>
//eqn_common has calc_E, calc_H

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultF,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	<?=vec3?> const E = calc_E(U);\
	<?=vec3?> const H = calc_H(U);\
	if (n.side == 0) {\
		(resultF)->D = _<?=vec3?>(<?=zero?>, H.z, <?=neg?>(H.y));\
		(resultF)->B = _<?=vec3?>(<?=zero?>, <?=neg?>(E.z), E.y);\
	} else if (n.side == 1) {\
		(resultF)->D = _<?=vec3?>(<?=neg?>(H.z), <?=zero?>, H.x);\
		(resultF)->B = _<?=vec3?>(E.z, <?=zero?>, <?=neg?>(E.x));\
	} else if (n.side == 2) {\
		(resultF)->D = _<?=vec3?>(H.y, <?=neg?>(H.x), <?=zero?>);\
		(resultF)->B = _<?=vec3?>(<?=neg?>(E.y), E.x, <?=zero?>);\
	}\
	(resultF)->phi = <?=zero?>;\
	(resultF)->psi = <?=zero?>;\
	(resultF)->D = <?=vec3?>_zero;\
	(resultF)->rhoCharge = <?=zero?>;\
	(resultF)->sqrt_1_eps = <?=susc_t?>_zero;\
	(resultF)->sqrt_1_mu = <?=susc_t?>_zero;\
}

//// MODULE_NAME: <?=eigen_forInterface?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	/* this will fail with tensor susceptibility */\
	/* but it doesn't belong here -- this is only the scalar case */\
	/* for tensors, I should be eigen-decomposing the levi-civita times the tensor	 */\
	(eig)->sqrt_1_eps = <?=susc_t?>_real_mul(<?=susc_t?>_add((UL)->sqrt_1_eps, (UR)->sqrt_1_eps), .5);\
	(eig)->sqrt_1_mu = <?=susc_t?>_real_mul(<?=susc_t?>_add((UL)->sqrt_1_mu, (UR)->sqrt_1_mu), .5);\
}

//// MODULE_NAME: <?=eigen_forCell?>

//used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	(eig)->sqrt_1_eps = U.sqrt_1_eps;\
	(eig)->sqrt_1_mu = U.sqrt_1_mu:\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=sqrt_2_and_1_2?>

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/
#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */x,\
	/*<?=normal_t?> const */n\
) {\
<? if false then -- TODO this, which means changing ops to macros, and determining the global-ness of the input params ... ?>\
	<?=scalar?> * const Yp = (<?=scalar?>*)(Y)->ptr;\
	<?=scalar?> const * const Xp = (<?=scalar?>*)(X)->ptr;\
<? end ?>\
\
	<?=susc_t?> const sqrt_1_eps = (eig)->sqrt_1_eps;					/* (m^3 kg)^.5/(C s) */\
	<?=susc_t?> const sqrt_eps = <?=susc_t?>_inv(sqrt_1_eps);			/* (C s)/(m^3 kg)^.5 */\
	<?=susc_t?> const sqrt_1_mu = (eig)->sqrt_1_mu;						/* C/(kg m)^.5 */\
	<?=susc_t?> const sqrt_mu = <?=susc_t?>_inv(sqrt_1_mu);			/* (kg m)^.5/C */\
	<?=susc_t?> const v_p = <?=susc_t?>_mul(sqrt_1_eps, sqrt_1_mu);	/* m/s */\
\
	if (n.side == 0) {\
\
		(Y)->ptr[0] = ((((X)->ptr[2] * sqrt_mu) + ((X)->ptr[4] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(Y)->ptr[1] = ((((X)->ptr[1] * sqrt_mu) - ((X)->ptr[5] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(Y)->ptr[2] = (sqrt_2 * (X)->ptr[0]);\
		(Y)->ptr[3] = (sqrt_2 * (X)->ptr[3]);\
		(Y)->ptr[4] = ((-(((X)->ptr[2] * sqrt_mu) - ((X)->ptr[4] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(Y)->ptr[5] = ((((X)->ptr[1] * sqrt_mu) + ((X)->ptr[5] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
\
	} else if (n.side == 1) {\
\
		(Y)->ptr[0] = ((((X)->ptr[2] * sqrt_mu) - ((X)->ptr[3] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(Y)->ptr[1] = ((((X)->ptr[0] * sqrt_mu) + ((X)->ptr[5] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(Y)->ptr[2] = (sqrt_2 * (X)->ptr[1]);\
		(Y)->ptr[3] = (sqrt_2 * (X)->ptr[4]);\
		(Y)->ptr[4] = ((((X)->ptr[2] * sqrt_mu) + ((X)->ptr[3] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(Y)->ptr[5] = ((-(((X)->ptr[0] * sqrt_mu) - ((X)->ptr[5] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
\
	} else if (n.side == 2) {\
\
		(Y)->ptr[0] = ((((X)->ptr[1] * sqrt_mu) + ((X)->ptr[3] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(Y)->ptr[1] = ((((X)->ptr[0] * sqrt_mu) - ((X)->ptr[4] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(Y)->ptr[2] = (sqrt_2 * (X)->ptr[2]);\
		(Y)->ptr[3] = (sqrt_2 * (X)->ptr[5]);\
		(Y)->ptr[4] = ((-(((X)->ptr[1] * sqrt_mu) - ((X)->ptr[3] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(Y)->ptr[5] = ((((X)->ptr[0] * sqrt_mu) + ((X)->ptr[4] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
\
	}\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=sqrt_2_and_1_2?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */x,\
	/*<?=normal_t?> const */n\
) {\
<? if false then -- TODO this, which means changing ops to macros, and determining the global-ness of the input params ... ?>\
	<?=scalar?> * const Yp = (<?=scalar?>*)(Y)->ptr;\
	<?=scalar?> const * const Xp = (<?=scalar?>*)(X)->ptr;\
<? end ?>\
\
	<?=susc_t?> const sqrt_1_eps = (eig)->sqrt_1_eps;					/* (m^3 kg)^.5/(C s) */\
	<?=susc_t?> const sqrt_eps = <?=susc_t?>_inv(sqrt_1_eps);			/* (C s)/(m^3 kg)^.5 */\
	<?=susc_t?> const sqrt_1_mu = (eig)->sqrt_1_mu;						/* C/(kg m)^.5 */\
	<?=susc_t?> const sqrt_mu = <?=susc_t?>_inv(sqrt_1_mu);			/* (kg m)^.5/C */\
	<?=susc_t?> const v_p = <?=susc_t?>_mul(sqrt_1_eps, sqrt_1_mu);	/* m/s */\
\
	if (n.side == 0) {\
\
		(Y)->ptr[0] = ((X)->ptr[2] / sqrt_2);\
		(Y)->ptr[1] = ((-(sqrt_eps * ((X)->ptr[1] - (X)->ptr[5]))) / sqrt_2);\
		(Y)->ptr[2] = ((sqrt_eps * ((X)->ptr[0] - (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[3] = ((X)->ptr[3] / sqrt_2);\
		(Y)->ptr[4] = ((sqrt_mu * ((X)->ptr[0] + (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[5] = ((sqrt_mu * ((X)->ptr[1] + (X)->ptr[5])) / sqrt_2);\
\
	} else if (n.side == 1) {\
\
		(Y)->ptr[0] = ((sqrt_eps * ((X)->ptr[1] - (X)->ptr[5])) / sqrt_2);\
		(Y)->ptr[1] = ((X)->ptr[2] / sqrt_2);\
		(Y)->ptr[2] = ((-(sqrt_eps * ((X)->ptr[0] - (X)->ptr[4]))) / sqrt_2);\
		(Y)->ptr[3] = ((sqrt_mu * ((X)->ptr[0] + (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[4] = ((X)->ptr[3] / sqrt_2);\
		(Y)->ptr[5] = ((sqrt_mu * ((X)->ptr[1] + (X)->ptr[5])) / sqrt_2);\
\
	} else if (n.side == 2) {\
\
		(Y)->ptr[0] = ((-(sqrt_eps * ((X)->ptr[1] - (X)->ptr[5]))) / sqrt_2);\
		(Y)->ptr[1] = ((sqrt_eps * ((X)->ptr[0] - (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[2] = ((X)->ptr[2] / sqrt_2);\
		(Y)->ptr[3] = ((sqrt_mu * ((X)->ptr[0] + (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[4] = ((sqrt_mu * ((X)->ptr[1] + (X)->ptr[5])) / sqrt_2);\
		(Y)->ptr[5] = ((X)->ptr[3] / sqrt_2);\
\
	}\
\
	for (int i = numWaves; i < numStates; ++i) {\
		(Y)->ptr[i] = 0;\
	}\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X_,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
)\
	<?=fluxFromCons?>(result, solver, X_, cell, n)

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=coord_sqrt_det_g?> <?=eqn_common?> 

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;

	/* TODO J = J_f + J_b = J_f + J_P + J_M = J_f + dP/dt + curl M */
	deriv->D = <?=vec3?>_sub(deriv->D, U->J);

	/* for non-time-varying susceptibilities, here's the source term: */
	/* D_i,t ... = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_,j B_k - J_i */
	/* B_i,t ... = 1/sqrt(g) g_il epsBar^ljk (1/eps)_,j B_k */
	
	<?=vec3?> grad_1_mu = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[solver->stepsize.<?=xj?>].sqrt_1_mu, U[solver->stepsize.<?=xj?>].sqrt_1_mu),
			<?=mul?>(U[-solver->stepsize.<?=xj?>].sqrt_1_mu, U[-solver->stepsize.<?=xj?>].sqrt_1_mu)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>
	
	<?=vec3?> grad_1_eps = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[solver->stepsize.<?=xj?>].sqrt_1_eps, U[solver->stepsize.<?=xj?>].sqrt_1_eps),
			<?=mul?>(U[-solver->stepsize.<?=xj?>].sqrt_1_eps, U[-solver->stepsize.<?=xj?>].sqrt_1_eps)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>
	
	real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>{
		<?=cons_t?> flux;
//// MODULE_DEPENDS: <?=fluxFromCons?>
		<?=fluxFromCons?>(&flux, solver, U, cell, normal_forSide<?=j?>(x));
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = <?=sub?>(deriv->D.<?=xj?>, <?=vec3?>_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = <?=sub?>(deriv->B.<?=xj?>, <?=vec3?>_dot(flux.B, grad_1_eps));
	}<? end ?>
}
