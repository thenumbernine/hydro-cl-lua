//// MODULE_NAME: <?=setFlatSpace?>

void <?=setFlatSpace?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	(U)->alpha = 1; 
	(U)->gamma_xx = 1;
	(U)->a_x = 0;
	(U)->D_g = 0;
	(U)->KTilde = 0;
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: sym3 coordMap

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
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;

	<?=initCode()?>

	U->alpha = alpha;
	U->gamma_xx = gamma_ll.xx;
	U->KTilde = K_ll.xx / sqrt(gamma_ll.xx);

}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> SETBOUNDS numGhost

kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	
	real const dx_alpha = (U[1].alpha - U[-1].alpha) / solver->grid_dx.x;
	real const dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / solver->grid_dx.x;

	U->a_x = dx_alpha / U->alpha;
	U->D_g = dx_gamma_xx / U->gamma_xx;
}


//// MODULE_NAME: <?=eigen_forInterface?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> const * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	(eig)->alpha = .5 * ((UL)->alpha + (UR)->alpha);\
	real const gamma_xx = .5 * ((UL)->gamma_xx + (UR)->gamma_xx);\
	(eig)->sqrt_gamma_xx = sqrt(gamma_xx);\
	(eig)->sqrt_f = sqrt(calc_f((eig)->alpha));\
}

//// MODULE_NAME: <?=eigen_forCell?>

//used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> const * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	(eig)->alpha = (U)->alpha;\
	(eig)->sqrt_gamma_xx = sqrt((U)->gamma_xx);\
	(eig)->sqrt_f = sqrt(calc_f((U)->alpha));\
}

//// MODULE_NAME: <?=eigen_leftTransform?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */x,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const f = (eig)->sqrt_f * (eig)->sqrt_f;\
	(result)->ptr[0] = ((x)->ptr[2] / f - (x)->ptr[4] / (eig)->sqrt_f) / 2.;\
	(result)->ptr[1] = -2. * (x)->ptr[2] / f + (x)->ptr[3];\
	(result)->ptr[2] = ((x)->ptr[2] / f + (x)->ptr[4] / (eig)->sqrt_f) / 2.;\
}

//// MODULE_NAME: <?=eigen_rightTransform?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */x,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const f = (eig)->sqrt_f * (eig)->sqrt_f;\
	(result)->ptr[0] = 0;\
	(result)->ptr[1] = 0;\
	(result)->ptr[2] = ((x)->ptr[0] + (x)->ptr[2]) * f;\
	(result)->ptr[3] = 2. * ((x)->ptr[0] + (x)->ptr[2]) + (x)->ptr[1];\
	(result)->ptr[4] = (eig)->sqrt_f * ((x)->ptr[2] - (x)->ptr[0]);\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */x,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const f = (eig)->sqrt_f * (eig)->sqrt_f;\
	real const alpha_over_sqrt_gamma_xx = (eig)->alpha / (eig)->sqrt_gamma_xx;\
	(result)->ptr[0] = 0;\
	(result)->ptr[1] = 0;\
	(result)->ptr[2] = (x)->ptr[4] * f * alpha_over_sqrt_gamma_xx;\
	(result)->ptr[3] = (x)->ptr[4] * 2. * alpha_over_sqrt_gamma_xx;\
	(result)->ptr[4] = (x)->ptr[2] * alpha_over_sqrt_gamma_xx;\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: initCond.codeprefix

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
	real const sqrt_gamma_xx = sqrt(U->gamma_xx);
	real const K_xx = U->KTilde / sqrt_gamma_xx;
	real const K = U->KTilde / sqrt_gamma_xx;

	real const f = calc_f(U->alpha);
	real const alphaSq_f = calc_f_alphaSq(U->alpha);
	real const alpha_dalpha_f = calc_alpha_dalpha_f(U->alpha);
	
	deriv->alpha -= alphaSq_f * K;
	deriv->gamma_xx -= 2. * U->alpha * U->gamma_xx * K;
	deriv->a_x -= ((.5 * U->D_g - U->a_x) * f - alpha_dalpha_f * U->a_x) * U->alpha * K;
	deriv->D_g -= (.5 * U->D_g - U->a_x) * 2 * U->alpha * K;
	deriv->KTilde -= (.5 * U->D_g - U->a_x) * U->a_x * U->alpha / sqrt_gamma_xx;

	// and now for the first-order constraints
	
	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	real const dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
	deriv->a_x += solver->a_x_convCoeff * (dx_alpha / U->alpha - U->a_x);
	
	// D_g = gamma_xx,x / gamma_xx <=> D_g += eta (gamma_xx,x / gamma_xx - D_g)
	real const dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
	deriv->D_g += solver->D_g_convCoeff * (dx_gamma_xx / U->gamma_xx - U->D_g);

	//Kreiss-Oligar diffusion, for stability's sake?
}
