//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real const r = fabs(x.x);
	cplx q = cplx_zero;
<?=initCode()
?>	U->psi = q;
	U->zeta = cplx_zero;
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
	
	//only 1D for now ...
	real const r = fabs(x.x);
	real const r2 = r * r;
	real const r4 = r2 * r2;
	
	real const h = solver->grid_dx.x;
	real const ih2 = 1. / (h * h);

	real const C = solver->C;
	real const m = solver->m;
	real const Vm = (m*m - .5) / r2;

	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
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
