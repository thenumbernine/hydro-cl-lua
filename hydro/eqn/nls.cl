//// MODULE_NAME: <?=calcDT?>

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cellBuf[index].pos;

	real r = fabs(x.x);
	cplx q = cplx_zero;
	<?=initCode()?>
	U->q = q;
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS_NOGHOST?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	real3 const x = cellBuf[index].pos;

	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

	//only 1D for now ...
	real const r = fabs(x.x);
	real const dr = solver->grid_dx.x;
	
	cplx const q_2L = U[-2].q;
	cplx const q_1L = U[-1].q;
	cplx const q = U[0].q;
	cplx const q_1R = U[1].q;
	cplx const q_2R = U[2].q;
	cplx const d2q_dx2 = cplx_real_mul(
		cplx_add5(
			cplx_neg(q_2L),
			cplx_real_mul(q_1L, 16.),
			cplx_real_mul(q, -30.),
			cplx_real_mul(q_1R, 16.),
			cplx_neg(q_2R)),
		 1. / (12. * dr * dr));
	cplx const dq_dx = cplx_real_mul(
		cplx_add4(
			q_2L,
			cplx_real_mul(q_1L, -8.),
			cplx_real_mul(q_1R, 8.),
			cplx_neg(q_2R)),
		 1. / (12. * dr));
	real const q2 = cplx_lenSq(q);
	real const q4 = q2 * q2;
	deriv->q = cplx_add(
		deriv->q,
		cplx_mul( 
			_cplx(0., 1.),
			cplx_add3(
				d2q_dx2,
				cplx_real_mul(dq_dx, 4. / r),
				cplx_real_mul(q, -q4))
		)
	);
}
