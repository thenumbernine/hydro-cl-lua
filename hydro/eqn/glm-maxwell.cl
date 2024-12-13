//// MODULE_NAME: <?=sqrt_1_2?>

#define sqrt_1_2 <?=("%.50f"):format(math.sqrt(.5))?>

//// MODULE_NAME: <?=cons_setEB?>

#define <?=cons_setEB?>(/*cons_t & */U, /*real3 */Eval, /*real3 */Bval) {\
	(U).B = Bval;\
	(U).D = <?=vec3?>_real_mul(Eval, <?=mul?>((U).sqrt_1_eps, (U).sqrt_1_eps));\
}

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

<?=vec3?> calc_E(global <?=cons_t?> const * const U) { 
	return <?=vec3?>_<?=susc_t?>_mul((U)->D, <?=susc_t?>_mul((U)->sqrt_1_eps, (U)->sqrt_1_eps));
}
<?=vec3?> calc_H(global <?=cons_t?> const * const U) { 
	return <?=vec3?>_<?=susc_t?>_mul((U)->B, <?=susc_t?>_mul((U)->sqrt_1_mu, (U)->sqrt_1_mu));
}

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

	//used
	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> rhoCharge = <?=zero?>;
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
	U->rhoCharge = rhoCharge;
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permittivity));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permeability));
}


//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=normal_t?> units <?=solver_t?> <?=cons_t?> <?=eqn_common?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
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
	(resultFlux)->D.x = H.y * nz - H.z * ny + nx * (U)->phi * solver->divPhiWavespeed / unit_m_per_s;	/* F_D^i = -eps^ijk n_j H_k */\
	(resultFlux)->B.x = E.z * ny - E.y * nz + nx * (U)->psi * solver->divPsiWavespeed / unit_m_per_s;	/* F_B^i = +eps^ijk n_j B_k */\
\
	(resultFlux)->D.y = H.z * nx - H.x * nz + ny * (U)->phi * solver->divPhiWavespeed / unit_m_per_s;\
	(resultFlux)->B.y = E.x * nz - E.z * nx + ny * (U)->psi * solver->divPsiWavespeed / unit_m_per_s;\
\
	(resultFlux)->D.z = H.x * ny - H.y * nx + nz * (U)->phi * solver->divPhiWavespeed / unit_m_per_s;\
	(resultFlux)->B.z = E.y * nx - E.x * ny + nz * (U)->psi * solver->divPsiWavespeed / unit_m_per_s;\
\
	real const D_n = normal_vecDotN1(n, (U)->D);\
	real const B_n = normal_vecDotN1(n, (U)->B);\
	(resultFlux)->phi = <?=real_mul?>(D_n, solver->divPhiWavespeed / unit_m_per_s);\
	(resultFlux)->psi = <?=real_mul?>(B_n, solver->divPsiWavespeed / unit_m_per_s);\
	(resultFlux)->sigma = <?=zero?>;\
	(resultFlux)->rhoCharge = <?=zero?>;\
	(resultFlux)->sqrt_1_eps = <?=susc_t?>_zero;\
	(resultFlux)->sqrt_1_mu = <?=susc_t?>_zero;\
}


//// MODULE_NAME: <?=eigen_forInterface?>

//this is here for GLMMaxwell:eigenMinMaxWhatever , same as hydro/eqn/maxwell.lua
//// MODULE_DEPENDS: <?=coord_sqrt_det_g?>

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
	/*<?=normal_t?> const */n\
) {\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = 1. / (eig)->sqrt_1_eps;\
	real const tmp2 = tmp1 * n_l.x;\
	real const tmp3 = 1. / 2.;\
	real const tmp4 = (X)->ptr[0] * tmp2;\
	real const tmp5 = tmp3 * tmp4;\
	real const tmp7 = n_l.y * tmp1;\
	real const tmp9 = (X)->ptr[1] * tmp7;\
	real const tmp11 = tmp9 * tmp3;\
	real const tmp13 = n_l.z * tmp1;\
	real const tmp15 = (X)->ptr[2] * tmp13;\
	real const tmp17 = tmp15 * tmp3;\
	real const tmp20 = (X)->ptr[6] * tmp1;\
	real const tmp22 = tmp20 * tmp3;\
	real const tmp26 = (X)->ptr[3] * tmp2;\
	real const tmp27 = tmp26 * tmp3;\
	real const tmp31 = (X)->ptr[4] * tmp7;\
	real const tmp33 = tmp31 * tmp3;\
	real const tmp37 = (X)->ptr[5] * tmp13;\
	real const tmp39 = tmp37 * tmp3;\
	real const tmp42 = (X)->ptr[7] * tmp1;\
	real const tmp44 = tmp42 * tmp3;\
	real const tmp45 = (eig)->sqrt_1_mu * n2_l.z;\
	real const tmp47 = (X)->ptr[5] * tmp45;\
	real const tmp48 = (eig)->sqrt_1_mu * n2_l.y;\
	real const tmp50 = (X)->ptr[4] * tmp48;\
	real const tmp51 = (eig)->sqrt_1_mu * n2_l.x;\
	real const tmp53 = (X)->ptr[3] * tmp51;\
	real const tmp54 = (eig)->sqrt_1_eps * n3_l.z;\
	real const tmp56 = (X)->ptr[2] * tmp54;\
	real const tmp57 = (eig)->sqrt_1_eps * n3_l.y;\
	real const tmp59 = (X)->ptr[1] * tmp57;\
	real const tmp60 = (eig)->sqrt_1_eps * n3_l.x;\
	real const tmp62 = (X)->ptr[0] * tmp60;\
	real const tmp63 = tmp59 * tmp3;\
	real const tmp64 = tmp62 * tmp3;\
	real const tmp65 = tmp56 * tmp3;\
	real const tmp67 = tmp53 * tmp3;\
	real const tmp69 = tmp50 * tmp3;\
	real const tmp71 = tmp47 * tmp3;\
	real const tmp73 = (eig)->sqrt_1_eps * n2_l.x;\
	real const tmp75 = (X)->ptr[0] * tmp73;\
	real const tmp76 = tmp75 * tmp3;\
	real const tmp77 = (eig)->sqrt_1_eps * n2_l.y;\
	real const tmp79 = (X)->ptr[1] * tmp77;\
	real const tmp81 = tmp79 * tmp3;\
	real const tmp82 = (eig)->sqrt_1_eps * n2_l.z;\
	real const tmp84 = (X)->ptr[2] * tmp82;\
	real const tmp86 = tmp84 * tmp3;\
	real const tmp87 = (eig)->sqrt_1_mu * n3_l.x;\
	real const tmp89 = (X)->ptr[3] * tmp87;\
	real const tmp90 = (eig)->sqrt_1_mu * n3_l.y;\
	real const tmp92 = (X)->ptr[4] * tmp90;\
	real const tmp93 = (eig)->sqrt_1_mu * n3_l.z;\
	real const tmp95 = (X)->ptr[5] * tmp93;\
	real const tmp96 = tmp92 * tmp3;\
	real const tmp97 = tmp95 * tmp3;\
	real const tmp98 = tmp89 * tmp3;\
	(Y)->ptr[0] = tmp22 + -tmp5 - tmp11 - tmp17;\
	(Y)->ptr[1] = tmp44 + -tmp27 - tmp33 - tmp39;\
	(Y)->ptr[2] = tmp71 + tmp69 + tmp67 + tmp65 + tmp63 + tmp64;\
	(Y)->ptr[3] = tmp98 + tmp96 + tmp97 + -tmp76 - tmp81 - tmp86;\
	(Y)->ptr[4] = tmp67 + tmp69 + tmp71 + -tmp64 - tmp63 - tmp65;\
	(Y)->ptr[5] = tmp97 + tmp96 + tmp98 + tmp86 + tmp81 + tmp76;\
	(Y)->ptr[6] = tmp22 + tmp17 + tmp11 + tmp5;\
	(Y)->ptr[7] = tmp44 + tmp39 + tmp33 + tmp27;\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=waves_t?> <?=sqrt_1_2?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = (eig)->sqrt_1_eps * n_l.x;\
	real const tmp3 = 1. / (eig)->sqrt_1_eps;\
	real const tmp4 = tmp3 * n3_l.x;\
	real const tmp6 = n2_l.x * tmp3;\
	real const tmp22 = (eig)->sqrt_1_eps * n_l.y;\
	real const tmp25 = n3_l.y * tmp3;\
	real const tmp27 = n2_l.y * tmp3;\
	real const tmp43 = (eig)->sqrt_1_eps * n_l.z;\
	real const tmp46 = n3_l.z * tmp3;\
	real const tmp48 = n2_l.z * tmp3;\
	real const tmp66 = 1. / (eig)->sqrt_1_mu;\
	real const tmp67 = tmp66 * n2_l.x;\
	real const tmp69 = n3_l.x * tmp66;\
	real const tmp88 = n2_l.y * tmp66;\
	real const tmp90 = n3_l.y * tmp66;\
	real const tmp109 = n2_l.z * tmp66;\
	real const tmp111 = n3_l.z * tmp66;\
	(Y)->ptr[0] = (X)->ptr[6] * tmp1 + (X)->ptr[5] * tmp6 + (X)->ptr[2] * tmp4 - (X)->ptr[3] * tmp6 - (X)->ptr[4] * tmp4 + -(X)->ptr[0] * tmp1;\
	(Y)->ptr[1] = (X)->ptr[6] * tmp22 + (X)->ptr[5] * tmp27 + (X)->ptr[2] * tmp25 - (X)->ptr[3] * tmp27 - (X)->ptr[4] * tmp25 + -(X)->ptr[0] * tmp22;\
	(Y)->ptr[2] = (X)->ptr[6] * tmp43 + (X)->ptr[5] * tmp48 + (X)->ptr[2] * tmp46 - (X)->ptr[3] * tmp48 - (X)->ptr[4] * tmp46 + -(X)->ptr[0] * tmp43;\
	(Y)->ptr[3] = (X)->ptr[7] * tmp1 + (X)->ptr[5] * tmp69 + (X)->ptr[4] * tmp67 + (X)->ptr[3] * tmp69 + (X)->ptr[2] * tmp67 + -(X)->ptr[1] * tmp1;\
	(Y)->ptr[4] = (X)->ptr[7] * tmp22 + (X)->ptr[5] * tmp90 + (X)->ptr[4] * tmp88 + (X)->ptr[3] * tmp90 + (X)->ptr[2] * tmp88 + -(X)->ptr[1] * tmp22;\
	(Y)->ptr[5] = (X)->ptr[7] * tmp43 + (X)->ptr[5] * tmp111 + (X)->ptr[4] * tmp109 + (X)->ptr[3] * tmp111 + (X)->ptr[2] * tmp109 + -(X)->ptr[1] * tmp43;\
	(Y)->ptr[6] = (X)->ptr[6] * (eig)->sqrt_1_eps + (X)->ptr[0] * (eig)->sqrt_1_eps;\
	(Y)->ptr[7] = (X)->ptr[7] * (eig)->sqrt_1_eps + (X)->ptr[1] * (eig)->sqrt_1_eps;\
\
	for (int i = <?=eqn.numWaves?>; i < <?=eqn.numStates?>; ++i) {\
		(Y)->ptr[i] = 0;\
	}\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: <?=fluxFromCons?> 

#error need to do this one now.  dF/dU * U != U ... it is just a coincidence of Euler fluid equations.
#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X_,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) \
	<?=fluxFromCons?>(Y, solver, X, cell, n)

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=SETBOUNDS_NOGHOST?>

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
	// TODO TODO should it be D -= D / eps * sigma?  or D -= D * sigma, since it is originally E -= E * sigma, right?
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
	
//// MODULE_DEPENDS: <?=coord_sqrt_det_g?>
	real const _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
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
	
	<?=scalar?> rhoCharge = U->rhoCharge;
<?=initCond.rhoChargeCode and initCond.rhoChargeCode(solver) or ""?>
	deriv->phi = <?=add?>(deriv->phi, <?=real_mul?>(rhoCharge, solver->divPhiWavespeed / unit_m_per_s));
}
