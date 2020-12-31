kernel void diffuse(
	global <?=cons_t?> * const UNextBuf,
	global <?=cons_t?> const * const UBuf,
	real dt
) {
	<?=SETBOUNDS_NOGHOST?>()

	real const viscosity = .001;
	real const alpha = dt * <?=clnumber(solver.numCells)?> * viscosity;
	real const diag = (1. + 2. * dim * alpha);
	UNextBuf = (UBuf + source * dt + alpha * sumSkew) / diag;
}

kernel void advect(
	global <?=cons_t?> * const UNextBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
#error finishme
}

kernel void calcDiv(
	global real * const divBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> const * const U = UBuf + index;

	real div = 0.;
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		div -= .5 * (U[solver->stepsize.s<?=side?>].v.s<?=side?>
				- U[-solver->stepsize.s<?=side?>].v.s<?=side?>) / grid_dx<?=side?>;
	}<? end ?>
	divBuf[index] = div;
}

kernel void diffusePressure(
	global real * const PBuf,
	global real const * const divBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global real* P = PBuf + index;
	global real const * const div = divBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		*P += .25 * (div[solver->stepsize.s<?=side?>] - div[-solver->stepsize.s<?=side?>]) / grid_dx<?=side?>;
	}<? end ?>
}

kernel void project(
	global <?=cons_t?> * const UBuf,
	global real const * const PBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global real const * const P = PBuf + index;
	<?=cons_t?> * const U = UBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		U->v.s<?=side?> -= (P[solver->stepsize.s<?=side?>] - P[-solver->stepsize.s<?=side?>]) / grid_dx<?=side?>;
	}<? end ?>
}
