// Roe solver:

<? 
-- for cartesian, use the normal in cartesian coordinates 
-- however that 'x' is still being used for who knows what 
--  ... and could be used for grid-metric-based norms
-- but you can't use 'cartesianFromCoord' because we're using the cartesian basis already 
-- for coordinate or coordinate-orthonormal, use the coordinate-aligned normal
local gridIsCartesian = require 'coord.cartesian'.is(solver.coord)
local nonAxisAlignedNormals = solver.coord.vectorComponent == 'cartesian' and not gridIsCartesian
?>


//call this, initializing nL to the normalBasisForSide<?=side?>
// this calculates the rest based on the x coordinate
// TODO put this somewhere above where display functions can use it
void calcNormalBasisForGrid(
	real3x3* nL,
	real3x3* nU,
	real *nLen,
	real3 x
) {
	// if we are non-cartesian grid with cartesian components then the normal of the non-cartesian grid needs to be in cartesian components as well ...
<? if nonAxisAlignedNormals then ?>
	nL.x = coord_cartesianFromCoord(nL.x, xInt);
	nL.y = coord_cartesianFromCoord(nL.y, xInt);
	nL.z = coord_cartesianFromCoord(nL.z, xInt);
<? end ?>
<? if gridIsCartesian then -- or if the metric is identity ... cylinder surface ?>
	*nU = *nL;
	*nLen = 1;
<? else ?>
	//if we have a non-cartesian grid + cartesian component normal ... or if we have a non cartesian component normal ... then raising is different than the lowered version
	*nU = (real3x3){
		.x = coord_raise(nL->x, xInt),
		.y = coord_raise(nL->y, xInt),
		.z = coord_raise(nL->z, xInt),
	};
	*nLen = sqrt(real3_dot(nL->x, nU->x));
<? end ?>
}


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
		real dt_dx = dt / dx;

		real3 xIntL = xInt;
		xIntL.s<?=side?> -= dx;
		
		real3 xIntR = xInt;
		xIntR.s<?=side?> += dx;
	
<?=solver:getULRCode():gsub('\t', '\t\t')?>

		// here is where parallel propagator comes into play
		cons_t pUL = cons_parallelPropagate<?=side?>(*UL, xL, .5 * dx);
		cons_t pUR = cons_parallelPropagate<?=side?>(*UR, xR, -.5 * dx);

		int indexInt = side + dim * index;
		
		real3x3 nL = normalBasisForSide<?=side?>;
		real3x3 nU;
		real nLen;
		calcNormalBasisForGrid(&nL, &nU, &nLen, xInt);

		eigen_t eig = eigen_forInterface(solver, pUL, pUR, xInt, nL, nU, nLen);

<? if nonAxisAlignedNormals then ?>
<?=eqn:eigenWaveCodePrefixForNormal('nL', 'nU', 'nLen', 'eig', 'xInt'):gsub('\t', '\t\t')?>
<? else ?>
<?=eqn:eigenWaveCodePrefixForSide(side, 'eig', 'xInt'):gsub('\t', '\t\t')?>
<? end ?>

		waves_t fluxEig;
<? if not eqn.roeUseFluxFromCons then 
?>		cons_t UAvg;
		for (int j = 0; j < numIntStates; ++j) {
			UAvg.ptr[j] = .5 * (pUL.ptr[j] + pUR.ptr[j]);
		}
		
<? 	if nonAxisAlignedNormals then ?>
#error here
		fluxEig = eigen_leftTransformForNormal(solver, eig, UAvg, xInt, nL, nU, nLen);
<?	else ?>
		fluxEig = eigen_leftTransformForSide<?=side?>(solver, eig, UAvg, xInt);
<?	end
end 
?>
		cons_t deltaU;	
<? if solver.fluxLimiter > 1 then 
?>		int indexR2 = indexR + solver->stepsize.s<?=side?>;
		int indexL2 = indexL - solver->stepsize.s<?=side?>;
		<?=solver:getULRCode{indexL = 'indexL2', indexR = 'indexL', suffix='_L'}?>
		<?=solver:getULRCode{indexL = 'indexR', indexR = 'indexR2', suffix='_R'}?>	
		
		// here is where parallel propagator comes into play
		cons_t pUL_L = cons_parallelPropagate<?=side?>(*UL_L, xIntL, 1.5 * dx);		//xIntL2?
		cons_t pUL_R = cons_parallelPropagate<?=side?>(*UL_R, xIntL, .5 * dx);
		cons_t pUR_L = cons_parallelPropagate<?=side?>(*UR_L, xIntR, -.5 * dx);
		cons_t pUR_R = cons_parallelPropagate<?=side?>(*UR_R, xIntR, -1.5 * dx);	// xIntR2?
			
		cons_t deltaUL, deltaUR;
<? end 
?>		
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = pUR.ptr[j] - pUL.ptr[j];
<? if solver.fluxLimiter > 1 then 
?>			deltaUL.ptr[j] = pUR_L.ptr[j] - pUL_L.ptr[j];
			deltaUR.ptr[j] = pUR_R.ptr[j] - pUL_R.ptr[j];
<? end 
?>		}

<? if nonAxisAlignedNormals then ?>
		waves_t deltaUEig = eigen_leftTransformForNormal(solver, eig, deltaU, xInt, nL, nU, nLen);
<? 	if solver.fluxLimiter > 1 then ?>		
		eigen_t eigL = eigen_forInterface(solver, pUL_L, pUR_L, xInt, nL, nU, nLen);
		eigen_t eigR = eigen_forInterface(solver, pUL_R, pUR_R, xInt, nL, nU, nLen);
		waves_t deltaUEigL = eigen_leftTransformForNormal(solver, eigL, deltaUL, xIntL, nL, nU, nLen);
		waves_t deltaUEigR = eigen_leftTransformForNormal(solver, eigR, deltaUR, xIntR, nL, nU, nLen);
<? 	end ?>
<? else ?>
		waves_t deltaUEig = eigen_leftTransformForSide<?=side?>(solver, eig, deltaU, xInt);
<? 	if solver.fluxLimiter > 1 then ?>		
		eigen_t eigL = eigen_forInterface(solver, pUL_L, pUR_L, xInt, nL, nU, nLen);
		eigen_t eigR = eigen_forInterface(solver, pUL_R, pUR_R, xInt, nL, nU, nLen);
		waves_t deltaUEigL = eigen_leftTransformForSide<?=side?>(solver, eigL, deltaUL, xIntL);
		waves_t deltaUEigR = eigen_leftTransformForSide<?=side?>(solver, eigR, deltaUR, xIntR);
<? 	end ?>
<? end ?>

		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real lambda = <?=eqn:eigenWaveCode(side, 'eig', 'xInt', j)?>;

<? if not eqn.roeUseFluxFromCons then 
?>			fluxEig.ptr[j] *= lambda;
<? else 
?>			fluxEig.ptr[j] = 0.;
<? end 
?>			real sgnLambda = lambda >= 0 ? 1 : -1;
		
<? if solver.fluxLimiter > 1 then 
?>			real rEig;
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
?>				+ phi * (lambda * dt_dx - sgnLambda)
<? end 
?>			);
		}<? end ?>
		
		global cons_t* flux = fluxBuf + indexInt;

<? if nonAxisAlignedNormals then ?>
		*flux = eigen_rightTransformForNormal(solver, eig, fluxEig, xInt, nL, nU, nLen);
<? else ?>
		*flux = eigen_rightTransformForSide<?=side?>(solver, eig, fluxEig, xInt);
<? end ?>

<? if eqn.roeUseFluxFromCons then ?>
<? 	if nonAxisAlignedNormals then ?>
		cons_t FL = fluxFromConsForNormal(solver, *UL, xL, nL, nU, nLen);
		cons_t FR = fluxFromConsForNormal(solver, *UR, xR, nL, nU, nLen);
<?	else ?>
		cons_t FL = fluxFromConsForSide<?=side?>(solver, *UL, xL);
		cons_t FR = fluxFromConsForSide<?=side?>(solver, *UR, xR);
<?	end ?>

		cons_t pFL = cons_parallelPropagate<?=side?>(FL, xL, .5 * dx);
		cons_t pFR = cons_parallelPropagate<?=side?>(FR, xR, -.5 * dx);
		
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] += .5 * (pFL.ptr[j] + pFR.ptr[j]);
		}
<? end 
?>	}<? end ?>
}
