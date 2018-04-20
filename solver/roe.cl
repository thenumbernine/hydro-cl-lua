// Roe solver:

<? if solver.checkFluxError or solver.checkOrthoError then ?>
kernel void calcErrors(
	global error_t* errorBuf,
	const global <?=eqn.eigen_t?>* eigenBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexInt = side + dim * index;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
			
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		real orthoError = 0;
		real fluxError = 0;

<?	if solver.checkOrthoError then
?>		//the flux transform is F v = R Lambda L v, I = R L
		//but if numWaves < numIntStates then certain v will map to the nullspace 
		//so to test orthogonality for only numWaves dimensions, I will verify that Qinv Q v = v 
		//I = L R
		for (int k = 0; k < numWaves; ++k) {
			<?=eqn.cons_t?> basis;
			for (int j = 0; j < numStates; ++j) {
				basis.ptr[j] = k == j ? 1 : 0;
			}
			
			<?=eqn.waves_t?> eigenCoords = eigen_leftTransform_<?=side?>(*eig, basis, xInt);
			<?=eqn.cons_t?> newbasis = eigen_rightTransform_<?=side?>(*eig, eigenCoords, xInt);
		
			for (int j = 0; j < numWaves; ++j) {
				orthoError += fabs(newbasis.ptr[j] - basis.ptr[j]);
			}
		}
<? 	end
	if solver.checkFluxError then	
?>
		<?=eqn:eigenWaveCodePrefix(side, 'eig', 'xInt')?>

		for (int k = 0; k < numIntStates; ++k) {
			
//TODO find out which left/right/fluxTransform functions are writing more than they should
//I see errors in mhd and in adm3d
			//this only needs to be numIntStates in size
			//but just in case the left/right transforms are reaching past that memory boundary ...
			<?=eqn.cons_t?> basis;
			for (int j = 0; j < numStates; ++j) {
				basis.ptr[j] = k == j ? 1 : 0;
			}

			<?=eqn.waves_t?> eigenCoords = eigen_leftTransform_<?=side?>(*eig, basis, xInt);

			<?=eqn.waves_t?> eigenScaled;
			<? for j=0,eqn.numWaves-1 do ?>{
				const int j = <?=j?>;
				real wave_j = <?=eqn:eigenWaveCode(side, 'eig', 'xInt', j)?>;
				eigenScaled.ptr[j] = eigenCoords.ptr[j] * wave_j;
			}<? end ?>
		
			//once again, only needs to be numIntStates
			<?=eqn.cons_t?> newtransformed = eigen_rightTransform_<?=side?>(*eig, eigenScaled, xInt);

//this shouldn't need to be reset here
// but it will if leftTransform does anything destructive
for (int j = 0; j < numStates; ++j) {
	basis.ptr[j] = k == j ? 1 : 0;
}

			//once again, only needs to be numIntStates
			<?=eqn.cons_t?> transformed = eigen_fluxTransform_<?=side?>(*eig, basis, xInt);
			
			for (int j = 0; j < numIntStates; ++j) {
				fluxError += fabs(newtransformed.ptr[j] - transformed.ptr[j]);
			}
		}
<?	end
?>		errorBuf[indexInt] = (error_t){
			.ortho = orthoError,
			.flux = fluxError,
		};
	}<? end ?>
}
<? end ?>

kernel void calcDeltaUEig(
	global real* deltaUEigBuf,
	<?= solver.getULRArg ?>,
	const global <?=eqn.eigen_t?>* eigenBuf
) {
	SETBOUNDS(numGhost,numGhost-1);	
	real3 x = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize.s<?=side?>;
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		<?= solver.getULRCode ?>

		<?=eqn.cons_t?> deltaU;
		for (int j = 0; j < numIntStates; ++j) {
			deltaU.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}
	
		int indexInt = side + dim * index;	
		global real* deltaUEig = deltaUEigBuf + indexInt * numWaves;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		*(global <?=eqn.waves_t?>*)deltaUEig = eigen_leftTransform_<?=side?>(*eig, deltaU, xInt);
	}<? end ?>
}

<? if solver.fluxLimiter > 1 then ?>
kernel void calcREig(
	global real* rEigBuf,
	const global real* deltaUEigBuf,
	const global <?=eqn.eigen_t?>* eigenBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize.s<?=side?>;
		int indexR = index + stepsize.s<?=side?>;
		int indexInt = side + dim * index;
		int indexIntL = side + dim * indexL;
		int indexIntR = side + dim * indexR;
		global real* rEig = rEigBuf + indexInt * numWaves;
		const global real* deltaUEig = deltaUEigBuf + indexInt * numWaves;
		const global real* deltaUEigL = deltaUEigBuf + indexIntL * numWaves;
		const global real* deltaUEigR = deltaUEigBuf + indexIntR * numWaves;
		
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		<?=eqn:eigenWaveCodePrefix(side, 'eig', 'xInt')?>
		
		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			if (deltaUEig[j] == 0) {
				rEig[j] = 0;
			} else {
				real wave_j = <?=eqn:eigenWaveCode(side, 'eig', 'xInt', j)?>;
				if (wave_j >= 0) {
					rEig[j] = deltaUEigL[j] / deltaUEig[j];
				} else {
					rEig[j] = deltaUEigR[j] / deltaUEig[j];
				}
			}
		}<? end ?>
	}<? end ?>
}
<? end ?>

//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>,
	const global <?=eqn.eigen_t?>* eigenBuf, 
	const global real* deltaUEigBuf,
	real dt
<? if solver.fluxLimiter > 1 then ?>
	,const global real* rEigBuf
<? end ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / grid_dx<?=side?>;//dx<?=side?>_at(i);
		int indexL = index - stepsize.s<?=side?>;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
	
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;

		<?=eqn:eigenWaveCodePrefix(side, 'eig', 'xInt')?>

		<?=eqn.waves_t?> fluxEig;
<? if not eqn.roeUseFluxFromCons then ?>
		<?=eqn.cons_t?> UAvg;
		for (int j = 0; j < numIntStates; ++j) {
			UAvg.ptr[j] = .5 * (UL->ptr[j] + UR->ptr[j]);
		}
		fluxEig = eigen_leftTransform_<?=side?>(*eig, UAvg, xInt);
<? end ?>

		const global real* deltaUEig = deltaUEigBuf + numWaves * indexInt;
<? if solver.fluxLimiter > 1 then ?>
		const global real* rEig = rEigBuf + numWaves * indexInt;
<? end ?>

		<? for j=0,eqn.numWaves-1 do ?>{
			const int j = <?=j?>;
			real lambda = <?=eqn:eigenWaveCode(side, 'eig', 'xInt', j)?>;
<? if not eqn.roeUseFluxFromCons then ?>
			fluxEig.ptr[j] *= lambda;
<? else ?>
			fluxEig.ptr[j] = 0.;
<? end ?>
			real sgnLambda = lambda >= 0 ? 1 : -1;
		
<? if solver.fluxLimiter > 1 then ?>
			real phi = fluxLimiter(rEig[j]);
<? end ?>

			fluxEig.ptr[j] -= .5 * lambda * deltaUEig[j] * (sgnLambda
<? if solver.fluxLimiter > 1 then ?>
				+ phi * (lambda * dt_dx - sgnLambda)
<? end ?>			
			);
		}<? end ?>
		
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		*flux = eigen_rightTransform_<?=side?>(*eig, fluxEig, xInt);

<? if eqn.roeUseFluxFromCons then ?>
		
		real3 xL = xR;
		xL.s<?=side?> -= grid_dx<?=side?>;
		
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
		for (int j = 0; j < numIntStates; ++j) {
			flux->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
		}
<? end ?>
	}<? end ?>
}
