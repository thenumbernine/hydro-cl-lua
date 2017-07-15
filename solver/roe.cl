// Roe solver:

<? if solver.checkFluxError or solver.checkOrthoError then ?>
kernel void calcErrors(
	global error_t* errorBuf,
	const global real* waveBuf,
	const global <?=eqn.eigen_t?>* eigenBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexInt = side + dim * index;
		const global real* wave = waveBuf + numWaves * indexInt;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		real orthoError = 0;
		real fluxError = 0;

<?	if solver.checkOrthoError then
?>		//the flux transform is F v = R Lambda L v, I = R L
		//but if numWaves < numStates then certain v will map to the nullspace 
		//so to test orthogonality for only numWaves dimensions, I will verify that Qinv Q v = v 
		//I = L R
		for (int k = 0; k < numWaves; ++k) {
			real basis[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				basis[j] = k == j ? 1 : 0;
			}
			
			real eigenInvCoords[numStates];
			eigen_rightTransform_<?=side?>__global_(eigenInvCoords, eig, basis, xInt);
		
			real newbasis[numWaves];
			eigen_leftTransform_<?=side?>__global_(newbasis, eig, eigenInvCoords, xInt);
			
			for (int j = 0; j < numWaves; ++j) {
				orthoError += fabs(newbasis[j] - basis[j]);
			}
		}
<? 	end
	if solver.checkFluxError then	
?>		for (int k = 0; k < numStates; ++k) {
			real basis[numStates];
			for (int j = 0; j < numStates; ++j) {
				basis[j] = k == j ? 1 : 0;
			}

			real eigenCoords[numWaves];
			eigen_leftTransform_<?=side?>__global_(eigenCoords, eig, basis, xInt);

			real eigenScaled[numWaves];
			for (int j = 0; j < numWaves; ++j) {
				eigenScaled[j] = eigenCoords[j] * wave[j];
			}
			
			real newtransformed[numStates];
			eigen_rightTransform_<?=side?>__global_(newtransformed, eig, eigenScaled, xInt);
			
			real transformed[numStates];
			eigen_fluxTransform_<?=side?>__global_(transformed, eig, basis, xInt);
			
			for (int j = 0; j < numStates; ++j) {
				fluxError += fabs(newtransformed[j] - transformed[j]);
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
	SETBOUNDS(2,1);	
	real3 x = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		<?= solver.getULRCode ?>

		<?=eqn.cons_t?> deltaU;
		for (int j = 0; j < numStates; ++j) {
			deltaU.ptr[j] = UR->ptr[j] - UL->ptr[j];
		}
	
		int indexInt = side + dim * index;	
		global real* deltaUEig = deltaUEigBuf + indexInt * numWaves;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		eigen_leftTransform_<?=side?>_global_global_(deltaUEig, eig, deltaU.ptr, xInt);
	}<? end ?>
}

<? if solver.fluxLimiter[0] > 0 then ?>
kernel void calcREig(
	global real* rEigBuf,
	const global real* deltaUEigBuf,
	const global real* waveBuf
) {
	SETBOUNDS(2,1);
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		int indexR = index + stepsize[side];
		int indexInt = side + dim * index;
		int indexIntL = side + dim * indexL;
		int indexIntR = side + dim * indexR;
		global real* rEig = rEigBuf + indexInt * numWaves;
		const global real* deltaUEig = deltaUEigBuf + indexInt * numWaves;
		const global real* deltaUEigL = deltaUEigBuf + indexIntL * numWaves;
		const global real* deltaUEigR = deltaUEigBuf + indexIntR * numWaves;
		const global real* wave = waveBuf + indexInt * numWaves;
		for (int j = 0; j < numWaves; ++j) {
			if (deltaUEig[j] == 0) {
				rEig[j] = 0;
			} else {
				if (wave[j] >= 0) {
					rEig[j] = deltaUEigL[j] / deltaUEig[j];
				} else {
					rEig[j] = deltaUEigR[j] / deltaUEig[j];
				}
			}
		}
	}<? end ?>
}
<? end ?>

//TODO entropy fix ... for the Euler equations at least
kernel void calcFlux(
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>,
	const global real* waveBuf, 
	const global <?=eqn.eigen_t?>* eigenBuf, 
	const global real* deltaUEigBuf,
	real dt
<? if solver.fluxLimiter[0] > 0 then ?>
	,const global real* rEigBuf
<? end ?>
) {
	SETBOUNDS(2,1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / grid_dx<?=side?>;//dx<?=side?>_at(i);
		int indexL = index - stepsize[side];

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
	
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;
		const global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;

		real fluxEig[numWaves];
<? if not eqn.hasFluxFromCons then ?>
		<?=eqn.cons_t?> UAvg;
		for (int j = 0; j < numStates; ++j) {
			UAvg.ptr[j] = .5 * (UL->ptr[j] + UR->ptr[j]);
		}
		eigen_leftTransform_<?=side?>__global_(fluxEig, eig, UAvg.ptr, xInt);
<? end ?>

		const global real* lambdas = waveBuf + numWaves * indexInt;
		const global real* deltaUEig = deltaUEigBuf + numWaves * indexInt;
<? if solver.fluxLimiter[0] > 0 then ?>
		const global real* rEig = rEigBuf + numWaves * indexInt;
<? end ?>

		for (int j = 0; j < numWaves; ++j) {
			real lambda = lambdas[j];
<? if not eqn.hasFluxFromCons then ?>
			fluxEig[j] *= lambda;
<? else ?>
			fluxEig[j] = 0.;
<? end ?>
			real sgnLambda = lambda >= 0 ? 1 : -1;
		
<? if solver.fluxLimiter[0] > 0 then ?>
			real phi = fluxLimiter(rEig[j]);
<? end ?>

			fluxEig[j] -= .5 * lambda * deltaUEig[j] * (sgnLambda
<? if solver.fluxLimiter[0] > 0 then ?>
				+ phi * (lambda * dt_dx - sgnLambda)
<? end ?>			
			);
		}
			
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		eigen_rightTransform_<?=side?>_global_global_(flux->ptr, eig, fluxEig, xInt);
		
<? if eqn.hasFluxFromCons then ?>
		
		real3 xL = xR;
		xL.s<?=side?> -= grid_dx<?=side?>;
		
		<?=eqn.cons_t?> FL = fluxFromCons_<?=side?>(*UL, xL);
		<?=eqn.cons_t?> FR = fluxFromCons_<?=side?>(*UR, xR);
		for (int j = 0; j < numStates; ++j) {
			flux->ptr[j] += .5 * (FL.ptr[j] + FR.ptr[j]);
		}
<? end ?>
	}<? end ?>
}

kernel void calcDerivFromFlux(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS(2,2);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	
	real3 x = cell_x(i);
	real volume = volume_at(x);

	<? for side=0,solver.dim-1 do ?>{
		int indexIntL = <?=side?> + dim * index;
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * stepsize[<?=side?>]; 
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + indexIntR;

<? -- trying to apply Trangenstein to my holonomic coordinate code:
if false
and require 'geom.cylinder'.is(solver.geometry) 
and solver.dim == 2
then
?>
		real rL = x.x - .5 * grid_dx0;
		real rR = x.x + .5 * grid_dx0;
		real thetaL = x.y - .5 * grid_dx1;
		real thetaR = x.y + .5 * grid_dx1;
		real dr = rR - rL;
		real dtheta = thetaR - thetaL;
		
		//the volume of a discrete cell
		// is the integral of the jacobian determinant across the coord range of the cell
		//so integral from r1 to r2 of r = r^2/2 at r2 - r1 = (r2^2-r1^2)/2
<?	if side==0 then	-- r+ and r-
?>		real areaL = rL * dtheta;
		real areaR = rR * dtheta;
<?	elseif side==1 then	-- θ+ and θ-	
?>		real areaL = dr;
		real areaR = dr;
<?	end
?>
		//so the next question is -- how to incorporate the d/dxi ln(sqrt(g)) ...
		// is it magically already considered?
		//or do I have to also scale by that here (in addition to scaling by the sides) ?
		real volume = .5 * (rR*rR - rL*rL) * (thetaR - thetaL);

<? else
?>		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real volumeIntL = volume_at(xIntL);
		real areaL = volumeIntL / grid_dx<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		real volumeIntR = volume_at(xIntR);
		real areaR = volumeIntR / grid_dx<?=side?>;
<? end
?>
		for (int j = 0; j < numStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
	}<? end ?>
}
