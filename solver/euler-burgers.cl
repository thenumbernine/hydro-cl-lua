//I could save on more calculations by calculting the step after the interface velocity
// but the interface velocity is performed doing the calcDeriv, which ultimately depends on dt
//So I could save one pass by flaggin whether the int vel is up-to-date or not.
kernel void calcDT(
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);

	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real Cs = calc_Cs(&W);
	
	real dt = INFINITY;
<? for side=0,solver.dim-1 do 
?>	dt = min(dt, grid_dx<?=side?> / (Cs + fabs(W.v.s<?=side?>)));
<? end
?>	dtBuf[index] = dt;
}

kernel void calcIntVel(
	global real* intVelBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	
	int indexR = index;
	const global <?=eqn.cons_t?>* UR = UBuf + indexR;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;

		// eventually: < ?= solver.getULRCode ? >
	
		int indexInt = side + dim * index;
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		real rhoL = UL->rho;
		real vL = UL->m.s<?=side?> / rhoL;
		
		real rhoR = UR->rho;
		real vR = UR->m.s<?=side?> / rhoR;
	
		intVelBuf[indexInt] = .5 * (vL + vR);
	}<? end ?>
}

kernel void calcFlux(
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global real* intVelBuf,
	real dt
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	const global <?=eqn.cons_t?>* UR = UBuf + indexR;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		real dt_dx = dt / grid_dx<?=side?>;//dx<?=side?>_at(i);
		
		int indexL = index - stepsize[side];
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;
		
		int indexL2 = indexL - stepsize[side];
		const global <?=eqn.cons_t?>* UL2 = UBuf + indexL2;
		
		int indexR2 = index + stepsize[side];
		const global <?=eqn.cons_t?>* UR2 = UBuf + indexR2;

		// eventually: < ?= solver.getULRCode ? >
	
		int indexInt = side + dim * index;
		real intVel = intVelBuf[indexInt];
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		for (int j = 0; j < numStates; ++j) {
			real deltaL = UL->ptr[j] - UL2->ptr[j];
			real delta = UR->ptr[j] - UL->ptr[j];
			real deltaR = UR2->ptr[j] - UR->ptr[j];
		
			real theta = intVel >= 0. ? 1. : -1.;
			real r = (intVel >= 0. ? deltaL : deltaR) / delta;
		
<? if solver.fluxLimiter[0] > 0 then ?>
			real phi = fluxLimiter(r);
<? end ?>
			flux->ptr[j] = .5 * intVel * ((1. + theta) * UL->ptr[j] + (1. - theta) * UR->ptr[j])
<? if solver.fluxLimiter[0] > 0 then ?>
				+ .5 * delta * phi * fabs(intVel) * (1. - fabs(intVel * dt_dx))
<? end ?>
			;
		}
	}<? end ?>
}

kernel void computePressure(
	global real* PBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost-1,numGhost-2);
	real3 x = cell_x(i);
	
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real P = W.P;

<? if false then -- Von Newmannartificial viscosity
?>	real dvSqSum = 0.;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		const global <?=eqn.cons_t?>* UL = U - stepsize[side];
		const global <?=eqn.cons_t?>* UR = U + stepsize[side];
		
		real dv = UR->m.s<?=side?> / UR->rho
				- UL->m.s<?=side?> / UL->rho;
		dvSqSum += dv * dv;
	}<? end ?>
	const real zeta = 2.;
	P += .25 * zeta * zeta * U->rho * dvSqSum;
<? end
?>	PBuf[index] = P;
}

kernel void diffuseMomentum(
	global <?=eqn.cons_t?>* derivBuf,
	const global real* PBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index; 
	const global real* P = PBuf + index; 

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		real dP = P[stepsize[side]] - P[-stepsize[side]];
		deriv->m.s<?=side?> -= .5 * dP / grid_dx<?=side?>;
	}<? end ?>
}

kernel void diffuseWork(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global real* PBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index; 
	const global <?=eqn.cons_t?>* U = UBuf + index;
	const global real* P = PBuf + index;

	real dE = 0.;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		const global <?=eqn.cons_t?>* UL = U - stepsize[side];
		const global <?=eqn.cons_t?>* UR = U + stepsize[side];
		
		real vR = UR->m.s<?=side?> / UR->rho;
		real PR = P[stepsize[side]];
		real vL = UL->m.s<?=side?> / UL->rho;
		real PL = P[-stepsize[side]];
		dE -= .5 * (PR * vR - PL * vL) / grid_dx<?=side?>;
	}<? end ?>

	deriv->ETotal += dE;
}