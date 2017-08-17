typedef <?=eqn.cons_t?> U_t;

real cplx_norm(real2 z) {
	return z.x*z.x + z.y*z.y;
}

real2 cplx_mul(real2 a, real2 b) {
	return (real2)(
		a.x * b.x - a.y * b.y,
		a.x * b.y + a.y * b.x);
}

kernel void calcDeriv(
	global real2* derivBuf,
	const global real2* UBuf
) {
	SETBOUNDS(2,2);
	real3 x = cell_x(i);
	
	//only 1D for now ...
	real r = x.x;
	const real dr = grid_dx0;
	
	real2 q_2L = UBuf[index-2];
	real2 q_1L = UBuf[index-1];
	real2 q = UBuf[index];
	real2 q_1R = UBuf[index+1];
	real2 q_2R = UBuf[index+2];
	real2 d2q_dx2 = (-q_2L + 16 * q_1L - 30 * q + 16 * q_1R - q_2R) / (12 * dr * dr);
	real2 dq_dx = (q_2L - 8 * q_1L + 8 * q_1R - q_2R) / (12 * dr);
	real q2 = cplx_norm(q);
	real q4 = q2 * q2;
	derivBuf[index] = cplx_mul( (real2)(0,1),
		d2q_dx2 + (4. / r) * dq_dx - q4 * q
	);
}
