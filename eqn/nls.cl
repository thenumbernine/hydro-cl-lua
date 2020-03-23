typedef <?=eqn.cons_t?> cons_t;

#define cplx_add5(a,b,c,d,e) 	cplx_add(cplx_add(cplx_add(a,b),cplx_add(c,d)),e)

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(2,2);
	real3 x = cell_x(i);
	
	//only 1D for now ...
	real r = fabs(x.x);
	const real dr = solver->grid_dx.x;
	
	cplx q_2L = UBuf[index-2].q;
	cplx q_1L = UBuf[index-1].q;
	cplx q = UBuf[index].q;
	cplx q_1R = UBuf[index+1].q;
	cplx q_2R = UBuf[index+2].q;
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
	derivBuf[index].q = cplx_mul( 
		_cplx(0., 1.),
		cplx_add3(
			d2q_dx2,
			cplx_real_mul(dq_dx, 4. / r),
			cplx_real_mul(q, -q4))
	);
}
