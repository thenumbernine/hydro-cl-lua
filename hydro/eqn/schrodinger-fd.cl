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
// hydro.init.nls sets cplx q
<?=initCode()
?>	U->psi = q;
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
	
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;

	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

#if 0
<?=eqn:makePartial1"psi"?>	// you don't need partial1 if the connection is zero
<?=eqn:makePartial2"psi"?>	// does this even work as complex-symmetric ?
#endif

// TODO maybe make this a grid operator ...
// begin copy from hydro/op/poisson_jacobi.cl

//// MODULE_DEPENDS: <?=cell_dx_i?>
<? for j=0,solver.dim-1 do
?>	real const dx<?=j?> = cell_dx<?=j?>(x);
<? end
?>

	real3 xInt = x;
	real3 volL, volR;
<? for j=0,solver.dim-1 do 
?>	xInt.s<?=j?> = x.s<?=j?> - .5 * solver->grid_dx.s<?=j?>;
	// TODO instead of volume_intL as the avg between two cell volumes, and then divide by dx to get the face, instead, just store the face.
	real const volume_intL<?=j?> = .5 * (cell->volume + cell[-solver->stepsize.s<?=j?>].volume);
	volL.s<?=j?> = volume_intL<?=j?>;
	xInt.s<?=j?> = x.s<?=j?> + .5 * solver->grid_dx.s<?=j?>;
	// TODO instead of volume_intL as the avg between two cell volumes, and then divide by dx to get the face, instead, just store the face.
	real const volume_intR<?=j?> = .5 * (cell->volume + cell[solver->stepsize.s<?=j?>].volume);
	volR.s<?=j?> = volume_intR<?=j?>;
	xInt.s<?=j?> = x.s<?=j?>;
<? end 
?>	real const volAtX = cell->volume;

	cplx lapPsi = cplx_zero;

<? for j=0,solver.dim-1 do 
?>	lapPsi = cplx_add3(lapPsi,
		cplx_real_mul(U[solver->stepsize.s<?=j?>].psi, volR.s<?=j?> / (dx<?=j?> * dx<?=j?>)),
		cplx_real_mul(U[-solver->stepsize.s<?=j?>].psi, volL.s<?=j?> / (dx<?=j?> * dx<?=j?>)));
<? end 
?>	lapPsi = cplx_real_mul(lapPsi, 1. / volAtX);
	real const diag = (0.
<? for j=0,solver.dim-1 do 
?>		- (volR.s<?=j?> + volL.s<?=j?>) / (dx<?=j?> * dx<?=j?>)
<? end 
?>	) / volAtX;
	
	real const hBar = solver->hBar;
	real const m = solver->m;
	real const V = 0.;	// TODO a field of this

	// i ℏ Ψ_,t = - ℏ^2/(2 m) Ψ_;j^j + V Ψ
	// Ψ_,t = i ℏ/(2 m) Ψ_;j^j - i V / ℏ Ψ
	// Ψ_,t = i ℏ/(2 m) (g^ij Ψ_,ij - Γ^k_ij Ψ_,k) - i V / ℏ Ψ
	deriv->psi = cplx_add3(
		deriv->psi,
		// TODO is particle-mass constant?  after all ... you're talking about the mass *of a field* ...
		cplx_real_mul(cplx_mul(cplx_i, lapPsi), hBar / (2 * m)),
		// TODO is potential energy V complex?
		cplx_real_mul(cplx_mul(cplx_i, U->psi), -V / hBar)
	);
}
