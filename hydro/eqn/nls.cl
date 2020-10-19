typedef <?=solver.solver_t?> solver_t;
typedef <?=eqn.cons_t?> cons_t;

// this is defined ... in cplx3
#define cplx_add3(a,b,c)		(cplx_add(a, cplx_add(b,c)))
#define cplx_add4(a,b,c,d)		(cplx_add(a, cplx_add3(b,c,d)))
#define cplx_add5(a,b,c,d,e)	(cplx_add(a, cplx_add4(b,c,d,e)))

<? if moduleName == nil then ?>
<? elseif moduleName == "calcDT" then ?>
<? elseif moduleName == "applyInitCond" then ?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;

	real r = fabs(x.x);
	cplx q = cplx_zero;
	<?=initCode()?>
	UBuf[index].q = q;
}

<? elseif moduleName == "addSource" then ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);

	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	//only 1D for now ...
	real r = fabs(x.x);
	const real dr = solver->grid_dx.x;
	
	cplx q_2L = U[-2].q;
	cplx q_1L = U[-1].q;
	cplx q = U[0].q;
	cplx q_1R = U[1].q;
	cplx q_2R = U[2].q;
	cplx d2q_dx2 = cplx_real_mul(
		cplx_add5(
			cplx_neg(q_2L),
			cplx_real_mul(q_1L, 16.),
			cplx_real_mul(q, -30.),
			cplx_real_mul(q_1R, 16.),
			cplx_neg(q_2R)),
		 1. / (12. * dr * dr));
	cplx dq_dx = cplx_real_mul(
		cplx_add4(
			q_2L,
			cplx_real_mul(q_1L, -8.),
			cplx_real_mul(q_1R, 8.),
			cplx_neg(q_2R)),
		 1. / (12. * dr));
	real q2 = cplx_lenSq(q);
	real q4 = q2 * q2;
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

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
