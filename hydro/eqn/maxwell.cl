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
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
\
	(resultF)->D.x = H.y * nz - H.z * ny;	/* F_D^i = -eps^ijk n_j H_k */\
	(resultF)->B.x = E.z * ny - E.y * nz;	/* F_B^i = +eps^ijk n_j B_k */\
\
	(resultF)->D.y = H.z * nx - H.x * nz;\
	(resultF)->B.y = E.x * nz - E.z * nx;\
\
	(resultF)->D.z = H.x * ny - H.y * nx;\
	(resultF)->B.z = E.y * nx - E.x * ny;\
\
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
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	/* TODO add extra macro params for globa/not-global and cast pointers to 'scalar' ... */\
	/* ... and replace all operations with macros based on 'scalar' */\
\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = 1. / 2.;\
	real const tmp2 = n3_l.y * (X)->ptr[5];\
	real const tmp4 = n2_l.y * (X)->ptr[4];\
	real const tmp6 = n_l.y * (X)->ptr[3];\
	real const tmp7 = 1. / (eig)->sqrt_1_mu;\
	real const tmp8 = n3_l.z * (X)->ptr[2];\
	real const tmp9 = tmp7 * tmp8;\
	real const tmp11 = (eig)->sqrt_1_eps * tmp9;\
	real const tmp13 = n2_l.z * (X)->ptr[1];\
	real const tmp14 = tmp13 * tmp7;\
	real const tmp16 = (eig)->sqrt_1_eps * tmp14;\
	real const tmp18 = n_l.z * (X)->ptr[0];\
	real const tmp19 = tmp18 * tmp7;\
	real const tmp21 = (eig)->sqrt_1_eps * tmp19;\
	real const tmp22 = tmp16 * tmp1;\
	real const tmp23 = tmp21 * tmp1;\
	real const tmp24 = tmp11 * tmp1;\
	real const tmp26 = tmp6 * tmp1;\
	real const tmp28 = tmp4 * tmp1;\
	real const tmp30 = tmp1 * tmp2;\
	real const tmp33 = n_l.y * (X)->ptr[0];\
	real const tmp34 = tmp33 * tmp7;\
	real const tmp36 = (eig)->sqrt_1_eps * tmp34;\
	real const tmp37 = tmp36 * tmp1;\
	real const tmp39 = n2_l.y * (X)->ptr[1];\
	real const tmp40 = tmp39 * tmp7;\
	real const tmp42 = (eig)->sqrt_1_eps * tmp40;\
	real const tmp44 = tmp42 * tmp1;\
	real const tmp46 = n3_l.y * (X)->ptr[2];\
	real const tmp47 = tmp46 * tmp7;\
	real const tmp49 = (eig)->sqrt_1_eps * tmp47;\
	real const tmp51 = tmp49 * tmp1;\
	real const tmp53 = n_l.z * (X)->ptr[3];\
	real const tmp55 = n2_l.z * (X)->ptr[4];\
	real const tmp57 = n3_l.z * (X)->ptr[5];\
	real const tmp58 = tmp55 * tmp1;\
	real const tmp59 = tmp57 * tmp1;\
	real const tmp60 = tmp53 * tmp1;\
	(Y)->ptr[0] = tmp30 + tmp28 + tmp26 + tmp24 + tmp22 + tmp23;\
	(Y)->ptr[1] = tmp60 + tmp58 + tmp59 + -tmp37 - tmp44 - tmp51;\
	(Y)->ptr[2] = n3_l.x * (X)->ptr[2] + n2_l.x * (X)->ptr[1] + n_l.x * (X)->ptr[0];\
	(Y)->ptr[3] = n3_l.x * (X)->ptr[5] + n2_l.x * (X)->ptr[4] + n_l.x * (X)->ptr[3];\
	(Y)->ptr[4] = tmp30 + tmp28 + tmp26 - tmp24 - tmp22 - tmp23;\
	(Y)->ptr[5] = tmp59 + tmp58 + tmp60 + tmp51 + tmp44 + tmp37;\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=sqrt_2_and_1_2?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = 1. / (eig)->sqrt_1_eps;\
	(Y)->ptr[0] = (eig)->sqrt_1_mu * tmp1 * n_l.z * (X)->ptr[0] + (eig)->sqrt_1_mu * tmp1 * n_l.y * (X)->ptr[5] - (eig)->sqrt_1_mu * tmp1 * n_l.z * (X)->ptr[4] - (eig)->sqrt_1_mu * tmp1 * n_l.y * (X)->ptr[1] + n_l.x * (X)->ptr[2];\
	(Y)->ptr[1] = (eig)->sqrt_1_mu * tmp1 * n2_l.z * (X)->ptr[0] + (eig)->sqrt_1_mu * tmp1 * n2_l.y * (X)->ptr[5] - (eig)->sqrt_1_mu * tmp1 * n2_l.z * (X)->ptr[4] - (eig)->sqrt_1_mu * tmp1 * n2_l.y * (X)->ptr[1] + n2_l.x * (X)->ptr[2];\
	(Y)->ptr[2] = (eig)->sqrt_1_mu * tmp1 * n3_l.z * (X)->ptr[0] + (eig)->sqrt_1_mu * tmp1 * n3_l.y * (X)->ptr[5] - (eig)->sqrt_1_mu * tmp1 * n3_l.z * (X)->ptr[4] - (eig)->sqrt_1_mu * tmp1 * n3_l.y * (X)->ptr[1] + n3_l.x * (X)->ptr[2];\
	(Y)->ptr[3] = n_l.z * (X)->ptr[5] + n_l.y * (X)->ptr[4] + n_l.x * (X)->ptr[3] + n_l.z * (X)->ptr[1] + n_l.y * (X)->ptr[0];\
	(Y)->ptr[4] = n2_l.z * (X)->ptr[5] + n2_l.y * (X)->ptr[4] + n2_l.x * (X)->ptr[3] + n2_l.z * (X)->ptr[1] + n2_l.y * (X)->ptr[0];\
	(Y)->ptr[5] = n3_l.z * (X)->ptr[5] + n3_l.y * (X)->ptr[4] + n3_l.x * (X)->ptr[3] + n3_l.z * (X)->ptr[1] + n3_l.y * (X)->ptr[0];\
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
