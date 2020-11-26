//// MODULE_NAME: setFlatSpace
//// MODULE_DEPENDS: solver_t cons_t

void setFlatSpace(
	constant solver_t const * const solver,
	global cons_t * const U,
	real3 const x
) {
	(U)->alpha = 1; 
	(U)->gamma_xx = 1;
	(U)->a_x = 0;
	(U)->d_xxx = 0;
	(U)->K_xx = 0;
}

//// MODULE_NAME: applyInitCond
//// MODULE_DEPENDS: sym3 coordMap

kernel void applyInitCond(
	constant solver_t const * const solver,
	constant initCond_t const * const initCond,
	global cons_t * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS(0,0);
	real3 const x = cellBuf[index].pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global cons_t * const U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;

	<?=initCode()?>

	U->alpha = alpha;
	U->gamma_xx = gamma_ll.xx;
	U->K_xx = K_ll.xx;
}

//// MODULE_NAME: initDerivs
//// MODULE_DEPENDS: solver_t cons_t cell_t SETBOUNDS numGhost

kernel void initDerivs(
	constant solver_t const * const solver,
	global cons_t * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global cons_t * const U = UBuf + index;
	
	real const dx_alpha = (U[1].alpha - U[-1].alpha) / solver->grid_dx.x;
	real const dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / solver->grid_dx.x;
	
	U->a_x = dx_alpha / U->alpha;
	U->d_xxx = .5 * dx_gamma_xx;
}

//// MODULE_NAME: fluxFromCons
//// MODULE_DEPENDS: solver_t cons_t normal_t initCond.codeprefix

#define fluxFromCons(\
	/*cons_t * const */F,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	real const f = calc_f((U)->alpha);\
	(F)->alpha = 0;\
	(F)->gamma_xx = 0;\
	(F)->a_x = (U)->alpha * (U)->K_xx * f / (U)->gamma_xx;\
	(F)->d_xxx = (U)->alpha * (U)->K_xx;\
	(F)->K_xx = (U)->alpha * (U)->a_x;\
}

//// MODULE_NAME: eigen_forInterface
//// MODULE_DEPENDS: initCond.codeprefix

//used for interface eigen basis
#define eigen_forInterface(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */UL,\
	/*cons_t const * const */UR,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	(eig)->alpha = .5 * ((UL)->alpha + (UR)->alpha);\
	real const gamma_xx = .5 * ((UL)->gamma_xx + (UR)->gamma_xx);\
	real const f = calc_f((eig)->alpha);\
	(eig)->sqrt_f_over_gamma_xx = sqrt(f / gamma_xx);\
}

//// MODULE_NAME: eigen_forCell
//// MODULE_DEPENDS: initCond.codeprefix

//used by PLM
#define eigen_forCell(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	real const f = calc_f((U)->alpha);\
	(eig)->alpha = (U)->alpha;\
	(eig)->sqrt_f_over_gamma_xx = sqrt(f / (U)->gamma_xx);\
}


//// MODULE_NAME: eigen_left/rightTransform

#define eigen_leftTransform(\
	/*waves_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*cons_t const * const */x,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const gamma_xx_over_f = 1. / ((eig)->sqrt_f_over_gamma_xx * (eig)->sqrt_f_over_gamma_xx);\
	(result)->ptr[0] = .5 * ((x)->ptr[2] / (eig)->sqrt_f_over_gamma_xx - (x)->ptr[4]);\
	(result)->ptr[1] = (x)->ptr[3] - (x)->ptr[2] * gamma_xx_over_f;\
	(result)->ptr[2] = .5 * ((x)->ptr[2] / (eig)->sqrt_f_over_gamma_xx + (x)->ptr[4]);\
}

#define eigen_rightTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*waves_t const * const */x,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	(result)->ptr[0] = 0;\
	(result)->ptr[1] = 0;\
	(result)->ptr[2] = ((x)->ptr[0] + (x)->ptr[2]) * (eig)->sqrt_f_over_gamma_xx;\
	(result)->ptr[3] = ((x)->ptr[0] + (x)->ptr[2]) / (eig)->sqrt_f_over_gamma_xx + (x)->ptr[1];\
	(result)->ptr[4] = (x)->ptr[2] - (x)->ptr[0];\
}

//// MODULE_NAME: eigen_fluxTransform

#define eigen_fluxTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*cons_t const * const */x,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const f_over_gamma_xx = (eig)->sqrt_f_over_gamma_xx * (eig)->sqrt_f_over_gamma_xx;\
	(result)->ptr[0] = 0;\
	(result)->ptr[1] = 0;\
	(result)->ptr[2] = (x)->ptr[4] * (eig)->alpha * f_over_gamma_xx;\
	(result)->ptr[3] = (x)->ptr[4] * (eig)->alpha;\
	(result)->ptr[4] = (x)->ptr[2] * (eig)->alpha;\
}

//// MODULE_NAME: addSource
//// MODULE_DEPENDS: initCond.codeprefix

kernel void addSource(
	constant solver_t const * const solver,
	global cons_t * const derivBuf,
	global cons_t const * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	global cons_t * const deriv = derivBuf + index;
	global cons_t const * const U = UBuf + index;
	
	real const alpha = U->alpha;
	real const gamma_xx = U->gamma_xx;
	real const a_x = U->a_x;
	real const d_xxx = U->d_xxx;
	real const K_xx = U->K_xx;
	real const K = K_xx / gamma_xx;	
	real const f = calc_f(alpha);
	real const dalpha_f = calc_dalpha_f(alpha);
	
	deriv->alpha -= alpha * alpha * f * K;
	deriv->gamma_xx -= 2. * alpha * K_xx;
	deriv->K_xx += alpha / gamma_xx * (a_x * d_xxx - K_xx * K_xx);
// terms that mysteriously disappear when you compare the linearized flux matrix terms moved to source, vs the source that Alcubierre uses in his 1997 paper
// adding these neglected terms back in make things blow up
#if 0 
	deriv->a_x += ((2. * d_xxx / gamma_xx - a_x) * f - alpha * dalpha_f * a_x) * alpha * K;
	deriv->d_xxx -= alpha * a_x * K_xx;
	deriv->K_xx -= alpha * a_x * a_x; 
#endif

	// and now for the first-order constraints
	
	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	real const dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
	deriv->a_x += solver->a_x_convCoeff * (dx_alpha / alpha - a_x);
	
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	real const dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
	deriv->d_xxx += solver->d_xxx_convCoeff * (.5 * dx_gamma_xx - d_xxx);

	//Kreiss-Oligar diffusion, for stability's sake?
}
