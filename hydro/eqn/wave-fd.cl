//// MODULE_NAME: <?=calcDT?>

//// MODULE_NAME: <?=applyInitCond?>
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?> SETBOUNDS

kernel void <?=applyInitCond?>(
	constant <?=solver_t?>* solver,
	constant <?=initCond_t?>* initCond,
	global <?=cons_t?>* UBuf,
	const global <?=cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real r = fabs(x.x);
	cplx q = cplx_zero;
	<?=initCode()?>
	UBuf[index] = (<?=cons_t?>){
		.psi = q,
		.zeta = cplx_zero,
	};
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> SETBOUNDS_NOGHOST cell_x

kernel void <?=addSource?>(
	constant <?=solver_t?>* solver,
	global <?=cons_t?>* derivBuf,
	const global <?=cons_t?>* UBuf,
	const global <?=cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	//only 1D for now ...
	real r = fabs(x.x);
	real r2 = r * r;
	real r4 = r2 * r2;
	
	const real h = solver->grid_dx.x;
	const real ih2 = 1. / (h * h);

	const real C = solver->C;
	const real m = solver->m;
	const real Vm = (m*m - .5) / r2;

	global <?=cons_t?>* deriv = derivBuf + index;
	const global <?=cons_t?>* U = UBuf + index;
	
	deriv->psi = cplx_add(deriv->psi, U->zeta);
	deriv->zeta = cplx_add4(
		deriv->zeta,
		cplx_real_mul(cplx_mul(cplx_i, U[0].zeta), -2 * C * m / r2),
		cplx_real_mul(U[0].psi, C * C * m * m / r4 - Vm),
		cplx_real_mul(
			cplx_add3(
				U[1].psi,
				cplx_real_mul(U[0].psi, -2),
				U[-1].psi
			),
			ih2
		)
	);
}
