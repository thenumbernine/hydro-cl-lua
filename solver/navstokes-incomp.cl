kernel void diffuse(
	global <?=eqn.cons_t?>* UNextBuf,
	const global <?=eqn.cons_t?>* UBuf,
	real dt
) {
	SETBOUNDS_NOGHOST()

	const real viscosity = .001;
	const real alpha = dt * <?=clnumber(solver.volume)?> * viscosity;
	const real diag = (1. + 2. * dim * alpha);
	UNextBuf = (UBuf + source * dt + alpha * sumSkew) / diag;
}

kernel void advect(
	global <?=eqn.cons_t?>* UNextBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
#error finishme
}

kernel void calcDiv(
	global real* divBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	const global <?=eqn.cons_t?>* U = UBuf + index;

	real div = 0.;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		div -= .5 * (U[stepsize[side]].v.s<?=side?>
				- U[-stepsize[side]].v.s<?=side?>) / grid_dx<?=side?>;
	}<? end ?>
	divBuf[index] = div;
}

kernel void diffusePressure(
	global real* PBuf,
	const global real* divBuf
) {
	SETBOUNDS_NOGHOST();
	global real* P = PBuf + index;
	const global real* div = divBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		*P += .25 * (div[stepsize[side]] - div[-stepsize[side]]) / grid_dx<?=side?>;
	}<? end ?>
}

kernel void project(
	global <?=eqn.cons_t?>* UBuf,
	const global real* PBuf
) {
	SETBOUNDS_NOGHOST();
	const global real* P = PBuf + index;
	const <?=eqn.cons_t?>* U = UBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		U->v.s<?=side?> -= (P[stepsize[side]] - P[-stepsize[side]]) / grid_dx<?=side?>;
	}<? end ?>
}
