//// MODULE_NAME: <?=sqrt_1_2?>

#define sqrt_1_2 <?=("%.50f"):format(math.sqrt(.5))?>

//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: coordLenSq cartesianToCoord coord_lower

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

<?=vec3?> calc_E(global <?=cons_t?> const * const U) { 
	return <?=vec3?>_<?=susc_t?>_mul((U)->D, <?=susc_t?>_mul((U)->sqrt_1_eps, (U)->sqrt_1_eps));
}
<?=vec3?> calc_H(global <?=cons_t?> const * const U) { 
	return <?=vec3?>_<?=susc_t?>_mul((U)->B, <?=susc_t?>_mul((U)->sqrt_1_mu, (U)->sqrt_1_mu));
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=eqn_common?>

kernel void <?=applyInitCondCell?>(
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

	//used
	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=scalar?>_from_real(1.);
	<?=susc_t?> permittivity = <?=susc_t?>_from_real(1.);
	<?=susc_t?> permeability = <?=susc_t?>_from_real(1.);
	
	//throw-away
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
	
	<?=initCode()?>
	
	U->D = eqn_cartesianToCoord(D, x);
	U->B = eqn_cartesianToCoord(B, x);
	U->phi = <?=zero?>;
	U->psi = <?=zero?>;
	U->sigma = conductivity;
	U->rhoCharge = <?=zero?>;
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permittivity));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permeability));
}


//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: normal_t units <?=solver_t?> <?=cons_t?> <?=eqn_common?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */F,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=vec3?> const E = calc_E(U);\
	<?=vec3?> const H = calc_H(U);\
	if (n.side == 0) {\
		(F)->D = _<?=vec3?>(<?=real_mul?>((U)->phi, solver->divPhiWavespeed / unit_m_per_s),  H.z, <?=neg?>(H.y));\
		(F)->B = _<?=vec3?>(<?=real_mul?>((U)->psi, solver->divPsiWavespeed / unit_m_per_s), <?=neg?>(E.z),  E.y);\
	} else if (n.side == 1) {\
		(F)->D = _<?=vec3?>(<?=neg?>(H.z), <?=real_mul?>((U)->phi, solver->divPhiWavespeed / unit_m_per_s),  H.x);\
		(F)->B = _<?=vec3?>( E.z, <?=real_mul?>((U)->psi, solver->divPsiWavespeed / unit_m_per_s), <?=neg?>(E.x));\
	} else if (n.side == 2) {\
		(F)->D = _<?=vec3?>( H.y, <?=neg?>(H.x), <?=real_mul?>((U)->phi, solver->divPhiWavespeed / unit_m_per_s));\
		(F)->B = _<?=vec3?>(<?=neg?>(E.y),  E.x, <?=real_mul?>((U)->psi, solver->divPsiWavespeed / unit_m_per_s));\
	}\
	real const D_n = normal_vecDotN1(n, (U)->D);\
	real const B_n = normal_vecDotN1(n, (U)->B);\
	(F)->phi = <?=real_mul?>(D_n, solver->divPhiWavespeed / unit_m_per_s);\
	(F)->psi = <?=real_mul?>(B_n, solver->divPsiWavespeed / unit_m_per_s);\
	(F)->sigma = <?=zero?>;\
	(F)->rhoCharge = <?=zero?>;\
	(F)->sqrt_1_eps = <?=susc_t?>_zero;\
	(F)->sqrt_1_mu = <?=susc_t?>_zero;\
}


//// MODULE_NAME: <?=eigen_forInterface?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */pt,\
	/*normal_t const */n\
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
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	(eig)->sqrt_1_eps = (U)->sqrt_1_eps;\
	(eig)->sqrt_1_mu = (U)->sqrt_1_mu;\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=waves_t?> <?=sqrt_1_2?>

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/
#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=susc_t?> const sqrt_1_eps = (eig)->sqrt_1_eps;					/* (m^3 kg)^.5/(C s) */\
	<?=susc_t?> const sqrt_eps = <?=susc_t?>_inv(sqrt_1_eps);			/* (C s)/(m^3 kg)^.5 */\
	<?=susc_t?> const sqrt_1_mu = (eig)->sqrt_1_mu;						/* C/(kg m)^.5 */\
	<?=susc_t?> const sqrt_mu = <?=susc_t?>_inv(sqrt_1_mu);			/* (kg m)^.5/C */\
	<?=susc_t?> const v_p = <?=susc_t?>_mul(sqrt_1_eps, sqrt_1_mu);	/* m/s */\
\
	if (n.side == 0) {\
\
		(Y)->ptr[0] = ((X)->ptr[6] - (X)->ptr[0]) * sqrt_eps * sqrt_1_2;					/* (C^2 s)/(m^3.5 kg^.5) */\
		(Y)->ptr[1] = ((X)->ptr[7] - (X)->ptr[3]) * sqrt_eps * sqrt_1_2;					/* kg^.5/m^1.5 */\
		(Y)->ptr[2] = ((X)->ptr[4] * sqrt_eps + (X)->ptr[2] * sqrt_mu) * v_p * sqrt_1_2;	/* kg^.5/(m^.5 s)  */\
		(Y)->ptr[3] = ((X)->ptr[5] * sqrt_eps - (X)->ptr[1] * sqrt_mu) * v_p * sqrt_1_2;	/* kg^.5/(m^.5 s)  */\
		(Y)->ptr[4] = ((X)->ptr[4] * sqrt_eps - (X)->ptr[2] * sqrt_mu) * v_p * sqrt_1_2;	/* kg^.5/(m^.5 s)  */\
		(Y)->ptr[5] = ((X)->ptr[5] * sqrt_eps + (X)->ptr[1] * sqrt_mu) * v_p * sqrt_1_2;	/* kg^.5/(m^.5 s)  */\
		(Y)->ptr[6] = ((X)->ptr[0] + (X)->ptr[6]) * sqrt_eps * sqrt_1_2;					/* (C^2 s)/(m^3.5 kg^.5) */\
		(Y)->ptr[7] = ((X)->ptr[3] + (X)->ptr[7]) * sqrt_eps * sqrt_1_2;					/* kg^.5/m^1.5 */\
\
	} else if (n.side == 1) {\
\
		/* same units as x dir */\
		(Y)->ptr[0] = ((X)->ptr[6] - (X)->ptr[1]) * sqrt_eps * sqrt_1_2;\
		(Y)->ptr[1] = ((X)->ptr[7] - (X)->ptr[4]) * sqrt_eps * sqrt_1_2;\
		(Y)->ptr[2] = ((X)->ptr[3] * sqrt_eps - (X)->ptr[2] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[3] = ((X)->ptr[5] * sqrt_eps + (X)->ptr[0] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[4] = ((X)->ptr[3] * sqrt_eps + (X)->ptr[2] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[5] = ((X)->ptr[5] * sqrt_eps - (X)->ptr[0] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[6] = ((X)->ptr[1] + (X)->ptr[6]) * sqrt_eps * sqrt_1_2;\
		(Y)->ptr[7] = ((X)->ptr[4] + (X)->ptr[7]) * sqrt_eps * sqrt_1_2;\
\
	} else if (n.side == 2) {\
	\
		/* same units as x dir */\
		(Y)->ptr[0] = ((X)->ptr[6] - (X)->ptr[2]) * sqrt_eps * sqrt_1_2;\
		(Y)->ptr[1] = ((X)->ptr[7] - (X)->ptr[5]) * sqrt_eps * sqrt_1_2;\
		(Y)->ptr[2] = ((X)->ptr[3] * sqrt_eps + (X)->ptr[1] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[3] = ((X)->ptr[4] * sqrt_eps - (X)->ptr[0] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[4] = ((X)->ptr[3] * sqrt_eps - (X)->ptr[1] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[5] = ((X)->ptr[4] * sqrt_eps + (X)->ptr[0] * sqrt_mu) * v_p * sqrt_1_2;\
		(Y)->ptr[6] = ((X)->ptr[2] + (X)->ptr[6]) * sqrt_eps * sqrt_1_2;\
		(Y)->ptr[7] = ((X)->ptr[5] + (X)->ptr[7]) * sqrt_eps * sqrt_1_2;\
\
	}\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=waves_t?> <?=sqrt_1_2?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=scalar?> const sqrt_1_eps = (eig)->sqrt_1_eps;\
	<?=scalar?> const sqrt_1_mu = (eig)->sqrt_1_mu;\
	<?=scalar?> const sqrt_eps = <?=inv?>((eig)->sqrt_1_eps);\
	<?=scalar?> const sqrt_mu = <?=inv?>((eig)->sqrt_1_mu);\
	<?=scalar?> const sqrt_2 = <?=inv?>(sqrt_1_2);\
\
	if (n.side == 0) {\
\
		(Y)->ptr[0] = ((-((X)->ptr[0] - (X)->ptr[6])) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[1] = ((-(sqrt_eps * ((X)->ptr[3] - (X)->ptr[5]))) / sqrt_2);\
		(Y)->ptr[2] = ((sqrt_eps * ((X)->ptr[2] - (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[3] = ((-((X)->ptr[1] - (X)->ptr[7])) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[4] = ((sqrt_mu * ((X)->ptr[2] + (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[5] = ((sqrt_mu * ((X)->ptr[3] + (X)->ptr[5])) / sqrt_2);\
		(Y)->ptr[6] = (((X)->ptr[0] + (X)->ptr[6]) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[7] = (((X)->ptr[1] + (X)->ptr[7]) / (sqrt_2 * sqrt_eps));\
\
	} else if (n.side == 1) {\
\
		(Y)->ptr[0] = ((sqrt_eps * ((X)->ptr[3] - (X)->ptr[5])) / sqrt_2);\
		(Y)->ptr[1] = ((-((X)->ptr[0] - (X)->ptr[6])) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[2] = ((-(sqrt_eps * ((X)->ptr[2] - (X)->ptr[4]))) / sqrt_2);\
		(Y)->ptr[3] = ((sqrt_mu * ((X)->ptr[2] + (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[4] = ((-((X)->ptr[1] - (X)->ptr[7])) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[5] = ((sqrt_mu * ((X)->ptr[3] + (X)->ptr[5])) / sqrt_2);\
		(Y)->ptr[6] = (((X)->ptr[0] + (X)->ptr[6]) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[7] = (((X)->ptr[1] + (X)->ptr[7]) / (sqrt_2 * sqrt_eps));\
\
	} else if (n.side == 2) {\
\
		(Y)->ptr[0] = ((-(sqrt_eps * ((X)->ptr[3] - (X)->ptr[5]))) / sqrt_2);\
		(Y)->ptr[1] = ((sqrt_eps * ((X)->ptr[2] - (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[2] = ((-((X)->ptr[0] - (X)->ptr[6])) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[3] = ((sqrt_mu * ((X)->ptr[2] + (X)->ptr[4])) / sqrt_2);\
		(Y)->ptr[4] = ((sqrt_mu * ((X)->ptr[3] + (X)->ptr[5])) / sqrt_2);\
		(Y)->ptr[5] = ((-((X)->ptr[1] - (X)->ptr[7])) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[6] = (((X)->ptr[0] + (X)->ptr[6]) / (sqrt_2 * sqrt_eps));\
		(Y)->ptr[7] = (((X)->ptr[1] + (X)->ptr[7]) / (sqrt_2 * sqrt_eps));\
\
	}\
\
	for (int i = <?=eqn.numWaves?>; i < <?=eqn.numStates?>; ++i) {\
		(Y)->ptr[i] = 0;\
	}\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

#define <?=eigen_fluxTransform?>(Y, solver, eig, X, x, n) <?=fluxFromCons?>(Y, solver, X, x, n)

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: coord_sqrt_det_g <?=fluxFromCons?> SETBOUNDS_NOGHOST

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const derivBuf,
	const global <?=cons_t?>* UBuf,
	const global <?=cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 const x = cellBuf[index].pos;
	
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

	/* TODO J = J_f + J_b = J_f + J_P + J_M = J_f + dP/dt + curl M */
	deriv->D = <?=vec3?>_sub(
		deriv->D, 
		<?=vec3?>_<?=scalar?>_mul(
			U->D, 
			<?=mul?>(<?=mul?>(U->sqrt_1_eps, U->sqrt_1_eps), U->sigma)
		)
	);
	
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
		<?=fluxFromCons?>(&flux, solver, U, x, normal_forSide<?=j?>(x));
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = <?=sub?>(deriv->D.<?=xj?>, <?=vec3?>_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = <?=sub?>(deriv->B.<?=xj?>, <?=vec3?>_dot(flux.B, grad_1_eps));
	}<? end ?>
	
	deriv->phi = <?=add?>(deriv->phi, <?=real_mul?>(U->rhoCharge, solver->divPhiWavespeed / unit_m_per_s));
}
