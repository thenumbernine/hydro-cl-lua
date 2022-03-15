//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=primFromCons?> <?=eqn_guiVars_compileTime?>

<? if require "hydro.solver.gridsolver":isa(solver) then ?>

kernel void <?=calcDTCell?>(
	global real * const dt,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	real const Cs = <?=calc_Cs?>(solver, &W);

<? for side=0,solver.dim-1 do 
?>	*(dt) = min(*(dt), (real)solver->grid_dx.s<?=side?> / (Cs + fabs(W.v.s<?=side?>)));
<? end
?>
}

<? else -- mesh solver ?>

kernel void <?=calcDTCell?>(
	global real * const dt,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	global <?=face_t?> const * const faceBuf,	//[numFaces]
	global int const * const cellFaceIndexes	//[numCellFaceIndexes]
) {
	real3 const x = cell->pos;
	
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	real const Cs = <?=calc_Cs?>(solver, &W);

	for (int i = 0; i < cell->faceCount; ++i) {
		global <?=face_t?> const * const face = faceBuf + cellFaceIndexes[i + cell->faceOffset];
		<?=normal_t?> const n = normal_forFace(face);
		real const v_n = normal_vecDotN1(n, W.v);
		real const dx = face->area;
		*(dt) = (real)min(*(dt), dx / (Cs + fabs(v_n)));
	}
}

<? end -- mesh solver ?>

//// MODULE_NAME: EulerBurgers.solver
//// MODULE_DEPENDS: <?=SETBOUNDS?> <?=fluxLimiter?> <?=primFromCons?>

kernel void calcIntVel(
	constant <?=solver_t?> const * const solver,
	global real * const intVelBuf,
	global const <?=solver.getULRArg?>,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost-1);
	real3 const xR = cellBuf[index].pos;
	int const indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		int const indexL = index - solver->stepsize.s<?=side?>;

		<?=solver:getULRCode()?>
	
		int const indexInt = side + dim * index;
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		real const rhoL = UL->rho;
		real const vL = UL->m.s<?=side?> / rhoL;
		
		real const rhoR = UR->rho;
		real const vR = UR->m.s<?=side?> / rhoR;
	
		intVelBuf[indexInt] = .5 * (vL + vR);
	}<? end ?>
}

kernel void <?=calcFlux?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const fluxBuf,
	global const <?=solver.getULRArg?>,
	global real const * const intVelBuf,
	global <?=cell_t?> const * const cellBuf,
	realparam const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost-1);
	real3 const xR = cellBuf[index].pos;
	int const indexR = index;

	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
//// MODULE_DEPENDS: <?=cell_dx_i?>
		real const dt_dx = dt / cell_dx<?=side?>(i);
		
		int const indexL = index - solver->stepsize.s<?=side?>;
		<?=solver:getULRCode()?>
		
		int const indexL2 = indexL - solver->stepsize.s<?=side?>;
		int const indexR2 = index + solver->stepsize.s<?=side?>;
		<?=solver:getULRCode{suffix="2"}?>
	
		int const indexInt = side + dim * index;
		real const intVel = intVelBuf[indexInt];
		global <?=cons_t?> * const flux = fluxBuf + indexInt;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		for (int j = 0; j < numIntStates; ++j) {
			real const deltaL = UL->ptr[j] - UL2->ptr[j];
			real const delta = UR->ptr[j] - UL->ptr[j];
			real const deltaR = UR2->ptr[j] - UR->ptr[j];
		
			real const theta = intVel >= 0. ? 1. : -1.;
			real const r = (intVel >= 0. ? deltaL : deltaR) / delta;
		
<? if solver.fluxLimiter > 1 then ?>
			real const phi = <?=fluxLimiter?>(r);
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
	constant <?=solver_t?> const * const solver,
	global real * const PBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost-1, solver->numGhost-2);
	real3 const x = cellBuf[index].pos;
	
	global <?=cons_t?> const * const U = UBuf + index;
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	real P = W.P;

<? if false then -- Von Newmann artificial viscosity
?>	real dvSqSum = 0.;
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
		
		real const dv = UR->m.s<?=side?> / UR->rho
			- UL->m.s<?=side?> / UL->rho;
		dvSqSum += dv * dv;
	}<? end ?>
	real const zeta = 2.;
	P += .25 * zeta * zeta * U->rho * dvSqSum;
<? end
?>	PBuf[index] = P;
}

kernel void diffuseMomentum(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global real const * const PBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> * const deriv = derivBuf + index; 
	global real const * const P = PBuf + index; 

	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		real const dP = P[solver->stepsize.s<?=side?>] - P[-solver->stepsize.s<?=side?>];
		deriv->m.s<?=side?> -= .5 * dP / solver->grid_dx.s<?=side?>;
	}<? end ?>
}

kernel void diffuseWork(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global real const * const PBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> * const deriv = derivBuf + index; 
	global <?=cons_t?> const * const U = UBuf + index;
	global real const * const P = PBuf + index;

	real dE = 0.;
	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
		
		real const vR = UR->m.s<?=side?> / UR->rho;
		real const PR = P[solver->stepsize.s<?=side?>];
		real const vL = UL->m.s<?=side?> / UL->rho;
		real const PL = P[-solver->stepsize.s<?=side?>];
		dE -= .5 * (PR * vR - PL * vL) / solver->grid_dx.s<?=side?>;
	}<? end ?>

	deriv->ETotal += dE;
}
