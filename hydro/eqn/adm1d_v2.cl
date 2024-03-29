//// MODULE_NAME: <?=calc_gamma_ll?>

#define /*real3s3*/ <?=calc_gamma_ll?>(\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) _real3s3((U)->gamma_xx, 0, 0, 1, 0, 1)

//// MODULE_NAME: <?=calc_gamma_uu?>

#define /*real3s3*/ <?=calc_gamma_uu?>(\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) _real3s3(1. / (U)->gamma_xx, 0, 0, 1, 0, 1)

//// MODULE_NAME: <?=setFlatSpace?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?>

void <?=setFlatSpace?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	(U)->alpha = 1; 
	(U)->gamma_xx = 1;
	(U)->a_x = 0;
	(U)->d_xxx = 0;
	(U)->K_xx = 0;
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: real3s3 <?=coordMap?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	real3s3 gamma_ll = real3s3_ident;
	real3s3 K_ll = real3s3_zero;

	<?=initCode()?>

	U->alpha = alpha;
	U->gamma_xx = gamma_ll.xx;
	U->K_xx = K_ll.xx;
}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	
	real const dx_alpha = (U[1].alpha - U[-1].alpha) / solver->grid_dx.x;
	real const dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / solver->grid_dx.x;
	
	U->a_x = dx_alpha / U->alpha;
	U->d_xxx = .5 * dx_gamma_xx;
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultF,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const f = calc_f((U)->alpha);\
	(resultF)->alpha = 0;\
	(resultF)->gamma_xx = 0;\
	(resultF)->a_x = (U)->alpha * (U)->K_xx * f / (U)->gamma_xx;\
	(resultF)->d_xxx = (U)->alpha * (U)->K_xx;\
	(resultF)->K_xx = (U)->alpha * (U)->a_x;\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=initCond_codeprefix?>

//used for interface eigen basis
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
	(eig)->alpha = .5 * ((UL)->alpha + (UR)->alpha);\
	real const gamma_xx = .5 * ((UL)->gamma_xx + (UR)->gamma_xx);\
	real const f = calc_f((eig)->alpha);\
	(eig)->sqrt_f_over_gamma_xx = sqrt(f / gamma_xx);\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=initCond_codeprefix?>

//used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const f = calc_f((U)->alpha);\
	(eig)->alpha = (U)->alpha;\
	(eig)->sqrt_f_over_gamma_xx = sqrt(f / (U)->gamma_xx);\
}


//// MODULE_NAME: <?=eigen_leftTransform?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */x,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const gamma_xx_over_f = 1. / ((eig)->sqrt_f_over_gamma_xx * (eig)->sqrt_f_over_gamma_xx);\
	(result)->ptr[0] = .5 * ((x)->ptr[2] / (eig)->sqrt_f_over_gamma_xx - (x)->ptr[4]);\
	(result)->ptr[1] = (x)->ptr[3] - (x)->ptr[2] * gamma_xx_over_f;\
	(result)->ptr[2] = .5 * ((x)->ptr[2] / (eig)->sqrt_f_over_gamma_xx + (x)->ptr[4]);\
}

//// MODULE_NAME: <?=eigen_rightTransform?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */x,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	(result)->ptr[0] = 0;\
	(result)->ptr[1] = 0;\
	(result)->ptr[2] = ((x)->ptr[0] + (x)->ptr[2]) * (eig)->sqrt_f_over_gamma_xx;\
	(result)->ptr[3] = ((x)->ptr[0] + (x)->ptr[2]) / (eig)->sqrt_f_over_gamma_xx + (x)->ptr[1];\
	(result)->ptr[4] = (x)->ptr[2] - (x)->ptr[0];\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */x,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const f_over_gamma_xx = (eig)->sqrt_f_over_gamma_xx * (eig)->sqrt_f_over_gamma_xx;\
	(result)->ptr[0] = 0;\
	(result)->ptr[1] = 0;\
	(result)->ptr[2] = (x)->ptr[4] * (eig)->alpha * f_over_gamma_xx;\
	(result)->ptr[3] = (x)->ptr[4] * (eig)->alpha;\
	(result)->ptr[4] = (x)->ptr[2] * (eig)->alpha;\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=initCond_codeprefix?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
	real const K = U->K_xx / U->gamma_xx;	
	real const alphaSq_f = calc_f_alphaSq(U->alpha);
	
	deriv->alpha -= alphaSq_f * K;
	deriv->gamma_xx -= 2. * U->alpha * U->K_xx;
	deriv->K_xx += U->alpha / U->gamma_xx * (U->a_x * U->d_xxx - U->K_xx * U->K_xx);
// terms that mysteriously disappear when you compare the linearized flux matrix terms moved to source, vs the source that Alcubierre uses in his 1997 paper
// adding these neglected terms back in make things blow up ... 
// ... or it used to I guess, maybe now it doesn't because I replaced the f (~ 1/alpha) times alpha with analytically simplified versions?
#if 1 
	real const alphaSq_dalpha_f = calc_alphaSq_dalpha_f(U->alpha);
	real const alpha_f = calc_f_alpha(U->alpha);
	deriv->a_x += (
		(2. * U->d_xxx / U->gamma_xx - U->a_x) * alpha_f
		- alphaSq_dalpha_f * U->a_x
	) * K;
	deriv->d_xxx -= U->alpha * U->a_x * U->K_xx;
	deriv->K_xx -= U->alpha * U->a_x * U->a_x; 
#endif

	// and now for the first-order constraints
	
	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	real const dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
	deriv->a_x += solver->a_x_convCoeff * (dx_alpha / U->alpha - U->a_x);
	
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	real const dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
	deriv->d_xxx += solver->d_xxx_convCoeff * (.5 * dx_gamma_xx - U->d_xxx);

	//Kreiss-Oligar diffusion, for stability's sake?
}
