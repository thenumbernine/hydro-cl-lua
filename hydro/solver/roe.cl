// Roe solver:

<?=eqn.cons_t?> calcFluxForInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> pUL,
	<?=eqn.cons_t?> pUR,
	real3 xInt,
	normalInfo_t n
<? if solver.fluxLimiter > 1 then 
?>	,realparam dt_dx,
	<?=eqn.cons_t?> pUL_L,
	<?=eqn.cons_t?> pUL_R,
	<?=eqn.cons_t?> pUR_L,
	<?=eqn.cons_t?> pUR_R,
	real3 xIntL,
	real3 xIntR
<? end
?>) {
	eigen_t eig = eigen_forInterface(solver, pUL, pUR, xInt, n);

	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'xInt'):gsub('\n', '\n\t\t')?>

	waves_t fluxEig;
<? if not eqn.roeUseFluxFromCons then 
?>	cons_t UAvg;
	for (int j = 0; j < numIntStates; ++j) {
		UAvg.ptr[j] = .5 * (pUL.ptr[j] + pUR.ptr[j]);
	}
	
	fluxEig = eigen_leftTransform(solver, eig, UAvg, xInt, n);
<? end
?>
	cons_t deltaU;	
<? if solver.fluxLimiter > 1 then 
?>	cons_t deltaUL, deltaUR;
<? end 
?>		
	for (int j = 0; j < numStates; ++j) {
		deltaU.ptr[j] = pUR.ptr[j] - pUL.ptr[j];
<? if solver.fluxLimiter > 1 then 
?>		deltaUL.ptr[j] = pUR_L.ptr[j] - pUL_L.ptr[j];
		deltaUR.ptr[j] = pUR_R.ptr[j] - pUL_R.ptr[j];
<? end 
?>	}

	waves_t deltaUEig = eigen_leftTransform(solver, eig, deltaU, xInt, n);
<? 	if solver.fluxLimiter > 1 then ?>
	eigen_t eigL = eigen_forInterface(solver, pUL_L, pUR_L, xInt, n);
	eigen_t eigR = eigen_forInterface(solver, pUL_R, pUR_R, xInt, n);
	waves_t deltaUEigL = eigen_leftTransform(solver, eigL, deltaUL, xIntL, n);
	waves_t deltaUEigR = eigen_leftTransform(solver, eigR, deltaUR, xIntR, n);
<? 	end ?>

	<? for j=0,eqn.numWaves-1 do ?>{
		const int j = <?=j?>;
		real lambda = <?=eqn:eigenWaveCode('n', 'eig', 'xInt', j)?>;

<? if not eqn.roeUseFluxFromCons then 
?>		fluxEig.ptr[j] *= lambda;
<? else 
?>		fluxEig.ptr[j] = 0.;
<? end 
?>		real sgnLambda = lambda >= 0 ? 1 : -1;
	
<? if solver.fluxLimiter > 1 then 
?>		real rEig;
		if (deltaUEig.ptr[j] == 0) {
			rEig = 0;
		} else {
			if (lambda >= 0) {
				rEig = deltaUEigL.ptr[j] / deltaUEig.ptr[j];
			} else {
				rEig = deltaUEigR.ptr[j] / deltaUEig.ptr[j];
			}
		}
		real phi = fluxLimiter(rEig);
<? end 
?>
		fluxEig.ptr[j] -= .5 * lambda * deltaUEig.ptr[j] * (sgnLambda
<? if solver.fluxLimiter > 1 then 
?>			+ phi * (lambda * dt_dx - sgnLambda)
<? end 
?>		);
	}<? end ?>
	
	cons_t flux = eigen_rightTransform(solver, eig, fluxEig, xInt, n);

<? if eqn.roeUseFluxFromCons then ?>
	cons_t FL = fluxFromCons(solver, pUL, xInt, n);
	cons_t FR = fluxFromCons(solver, pUR, xInt, n);

	for (int j = 0; j < numIntStates; ++j) {
		flux.ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
	}
<? end 
?>
	return flux;
}


<? if not require 'hydro.solver.meshsolver'.is(solver) then ?>


//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.getULRArg?>,
	realparam dt
) {
	typedef <?=eqn.cons_t?> cons_t;
	typedef <?=eqn.eigen_t?> eigen_t;
	typedef <?=eqn.waves_t?> waves_t;
	
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		int indexL = index - solver->stepsize.s<?=side?>;

		real dx = solver->grid_dx.s<?=side?>;

		real3 xL = xR;
		xL.s<?=side?> -= dx;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * dx;
		
		//this is used for the flux limiter
		//should it be using the coordinate dx or the grid dx?
		//real dt_dx = dt / cell_dx<?=side?>(xInt);
<? 
if solver.coord.vectorComponent == 'cartesian' 
and not require 'hydro.coord.cartesian'.is(solver.coord)
then 
?>		real dt_dx = dt / cell_dx<?=side?>(xInt);
<? else 
?>		real dt_dx = dt / dx;
<? end 
?>	
		real3 xIntL = xInt;
		xIntL.s<?=side?> -= dx;
		
		real3 xIntR = xInt;
		xIntR.s<?=side?> += dx;
	
		<?=solver:getULRCode():gsub('\n', '\n\t\t')?>
		
		// here is where parallel propagator comes into play
		cons_t pUL = cons_parallelPropagate<?=side?>(*UL, xL, .5 * dx);
		cons_t pUR = cons_parallelPropagate<?=side?>(*UR, xR, -.5 * dx);
	
		normalInfo_t n = normalInfo_forSide<?=side?>(xInt);

<? if solver.fluxLimiter > 1 then 
?>		int indexR2 = indexR + solver->stepsize.s<?=side?>;
		int indexL2 = indexL - solver->stepsize.s<?=side?>;
		<?=solver:getULRCode{indexL = 'indexL2', indexR = 'indexL', suffix='_L'}:gsub('\n', '\n\t\t')?>
		<?=solver:getULRCode{indexL = 'indexR', indexR = 'indexR2', suffix='_R'}:gsub('\n', '\n\t\t')?>

		// here is where parallel propagator comes into play
		cons_t pUL_L = cons_parallelPropagate<?=side?>(*UL_L, xIntL, 1.5 * dx);		//xIntL2?
		cons_t pUL_R = cons_parallelPropagate<?=side?>(*UL_R, xIntL, .5 * dx);
		cons_t pUR_L = cons_parallelPropagate<?=side?>(*UR_L, xIntR, -.5 * dx);
		cons_t pUR_R = cons_parallelPropagate<?=side?>(*UR_R, xIntR, -1.5 * dx);	// xIntR2?
<? end
?>
		int indexInt = side + dim * index;
		global cons_t* flux = fluxBuf + indexInt;
		*flux = calcFluxForInterface(
			solver, pUL, pUR, xInt, n
<? if solver.fluxLimiter > 1 then
?>			,dt_dx, pUL_L, pUL_R, pUR_L, pUR_R, xIntL, xIntR
<? end
?>		);
	}<? end ?>
}

<? else -- meshsolver 
	if solver.fluxLimiter > 1 then
		error("MeshSolver doesn't work with fluxLimiter")
	end
?>

void getEdgeStates(
	const global face_t* e,
	<?=eqn.cons_t?>* UL,
	<?=eqn.cons_t?>* UR,
	const global <?=eqn.cons_t?>* UBuf		//[numCells]
	//const global cell_t* cells			//[numCells]
) {
	//const real resitution = 0.;
	int iL = e->cells.s0;
	int iR = e->cells.s1;
	//cell_t* cL = iL == -1 ? NULL : cells + e->cells.s0;
	//cell_t* cR = iR == -1 ? NULL : cells + e->cells.s1;
	if (iL != -1 && iR != -1) {
		*UL = UBuf[iL];
		*UR = UBuf[iR];
	} else if (iL != -1) {
		*UL = UBuf[iL];
		//TODO  
		*UR = *UL;//eqn.reflect(cL->U, e->normal, restitution);
	} else if (iR != -1) {
		*UR = UBuf[iR];
		//TODO  
		*UL = *UR;//eqn.reflect(cR->U, e->normal, restitution);
	} else {	// both iL and iR are null ...
		//error
		for (int i = 0; i < numStates; ++i) {
			UL->ptr[i] = UR->ptr[i] = 0./0.;
		}
	}
}

//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf,
	realparam dt,
//mesh-specific parameters:	
	const global cell_t* cells,			//[numCells]
	const global face_t* faces,			//[numFaces]
	const global int* cellFaceIndexes	//[numCellFaceIndexes]
) {
	typedef <?=eqn.cons_t?> cons_t;
	typedef <?=eqn.eigen_t?> eigen_t;
	typedef <?=eqn.waves_t?> waves_t;

	int faceIndex = get_global_id(0);
	if (faceIndex >= get_global_size(0)) return;
	
	const global face_t* face = faces + faceIndex;
	if (face->area <= 1e-7) return;
	
	real3 x = face->pos;
	normalInfo_t n = normalInfo_forFace(face);
	
	cons_t UL, UR;	
	getEdgeStates(face, &UL, &UR, UBuf);

	//TODO option to rotate to align fluxes?
	// then you'd have to build a new normalInfo_t based on the aligned (x-axis) normal.

	fluxBuf[faceIndex] = calcFluxForInterface(solver, UL, UR, x, n);
}

<? end -- mesh vs grid solver ?>