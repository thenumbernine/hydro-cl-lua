kernel void calcIntVel(
	constant <?=solver.solver_t?>* solver,
	global real* intVelBuf,
	const global <?=solver.getULRArg?>,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - solver->stepsize.s<?=side?>;

		<?=solver:getULRCode()?>
	
		int indexInt = side + dim * index;
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		real rhoL = UL->rho;
		real vL = UL->m.s<?=side?> / rhoL;
		
		real rhoR = UR->rho;
		real vR = UR->m.s<?=side?> / rhoR;
	
		intVelBuf[indexInt] = .5 * (vL + vR);
	}<? end ?>
}

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.getULRArg?>,
	const global real* intVelBuf,
	const global <?=solver.coord.cell_t?>* cellBuf,
	realparam dt
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		real dt_dx = dt / cell_dx<?=side?>(i);
		
		int indexL = index - solver->stepsize.s<?=side?>;
		<?=solver:getULRCode()?>
		
		int indexL2 = indexL - solver->stepsize.s<?=side?>;
		int indexR2 = index + solver->stepsize.s<?=side?>;
		<?=solver:getULRCode{suffix='2'}?>
	
		int indexInt = side + dim * index;
		real intVel = intVelBuf[indexInt];
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		for (int j = 0; j < numIntStates; ++j) {
			real deltaL = UL->ptr[j] - UL2->ptr[j];
			real delta = UR->ptr[j] - UL->ptr[j];
			real deltaR = UR2->ptr[j] - UR->ptr[j];
		
			real theta = intVel >= 0. ? 1. : -1.;
			real r = (intVel >= 0. ? deltaL : deltaR) / delta;
		
<? if solver.fluxLimiter > 1 then ?>
			real phi = fluxLimiter(r);
<? end ?>
			flux->ptr[j] = .5 * intVel * ((1. + theta) * UL->ptr[j] + (1. - theta) * UR->ptr[j])
<? if solver.fluxLimiter > 1 then ?>
				+ .5 * delta * phi * fabs(intVel) * (1. - fabs(intVel * dt_dx))
<? end ?>
			;
		}
	}<? end ?>
}

kernel void computePressure(
	constant <?=solver.solver_t?>* solver,
	global real* PBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost-1,numGhost-2);
	real3 x = cell_x(i);
	
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real P = W.P;

<? if false then -- Von Newmannartificial viscosity
?>	real dvSqSum = 0.;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;
		
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
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global real* PBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index; 
	const global real* P = PBuf + index; 

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		real dP = P[solver->stepsize.s<?=side?>] - P[-solver->stepsize.s<?=side?>];
		deriv->m.s<?=side?> -= .5 * dP / solver->grid_dx.s<?=side?>;
	}<? end ?>
}

kernel void diffuseWork(
	constant <?=solver.solver_t?>* solver,
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
		const global <?=eqn.cons_t?>* UL = U - solver->stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UR = U + solver->stepsize.s<?=side?>;
		
		real vR = UR->m.s<?=side?> / UR->rho;
		real PR = P[solver->stepsize.s<?=side?>];
		real vL = UL->m.s<?=side?> / UL->rho;
		real PL = P[-solver->stepsize.s<?=side?>];
		dE -= .5 * (PR * vR - PL * vL) / solver->grid_dx.s<?=side?>;
	}<? end ?>

	deriv->ETotal += dE;
}
