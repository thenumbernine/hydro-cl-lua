#define divPhiWavespeed 	1.
#define divPsiWavespeed		1.
#define sqrt_mu0 			1.
#define sqrt_eps0 			1.
#define mu0					(sqrt_mu0 * sqrt_mu0)
#define eps0				(sqrt_eps0 * sqrt_eps0)
#define speedOfLight		(1./(sqrt_mu0 * sqrt_eps0))
#define speedOfLightSq		(speedOfLight * speedOfLight)

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	real lambda = max(max(divPsiWavespeed, divPhiWavespeed), 1.) * speedOfLight;
	return (range_t){-lambda, lambda};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	return (<?=eqn.eigen_t?>){};
}
<? end ?>

//same as in eqn/euler.cl
kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
		
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		*eig = eigen_forSide_<?=side?>(UL, UR, xInt);
	}<? end ?>
}
		
<? for side=0,solver.dim-1 do ?>

<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=eqn.waves_t?> Y;

	<? if side==0 then ?>

	Y.ptr[0] = (X.ptr[0] - X.ptr[6] * speedOfLight) * .5;
	Y.ptr[1] = (X.ptr[3] - X.ptr[7] / speedOfLight) * .5;
	Y.ptr[2] = (X.ptr[1] - X.ptr[5] * speedOfLight) * .5;
	Y.ptr[3] = (X.ptr[4] * speedOfLight + X.ptr[2]) * .5;
	Y.ptr[4] = (X.ptr[5] * speedOfLight + X.ptr[1]) * .5;
	Y.ptr[5] = (X.ptr[2] - X.ptr[4] * speedOfLight) * .5;
	Y.ptr[6] = (X.ptr[7] / speedOfLight + X.ptr[3]) * .5;
	Y.ptr[7] = (X.ptr[6] * speedOfLight + X.ptr[0]) * .5;
   
   <? elseif side==1 then ?>
   
	Y.ptr[0] = (X.ptr[1] - X.ptr[6] * speedOfLight) * .5;
	Y.ptr[1] = (X.ptr[4] - X.ptr[7] / speedOfLight) * .5;
	Y.ptr[2] = (X.ptr[5] * speedOfLight + X.ptr[0]) * .5;
	Y.ptr[3] = (X.ptr[2] - X.ptr[3] * speedOfLight) * .5;
	Y.ptr[4] = (X.ptr[0] - X.ptr[5] * speedOfLight) * .5;
	Y.ptr[5] = (X.ptr[3] * speedOfLight + X.ptr[2]) * .5;
	Y.ptr[6] = (X.ptr[7] / speedOfLight + X.ptr[4]) * .5;
	Y.ptr[7] = (X.ptr[6] * speedOfLight + X.ptr[1]) * .5;
   
   <? elseif side==2 then ?>
   
	Y.ptr[0] = (X.ptr[2] - X.ptr[6] * speedOfLight) * .5;
	Y.ptr[1] = (X.ptr[5] - X.ptr[7] / speedOfLight) * .5;
	Y.ptr[2] = (X.ptr[0] - X.ptr[4] * speedOfLight) * .5;
	Y.ptr[3] = (X.ptr[3] * speedOfLight + X.ptr[1]) * .5;
	Y.ptr[4] = (X.ptr[4] * speedOfLight + X.ptr[0]) * .5;
	Y.ptr[5] = (X.ptr[1] - X.ptr[3] * speedOfLight) * .5;
	Y.ptr[6] = (X.ptr[7] / speedOfLight + X.ptr[5]) * .5;
	Y.ptr[7] = (X.ptr[6] * speedOfLight + X.ptr[2]) * .5;
	
	<? end ?>
	
	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;

	<? if side==0 then ?>
	
	Y.ptr[0] = X.ptr[7] + X.ptr[0];
	Y.ptr[1] = X.ptr[4] + X.ptr[2];
	Y.ptr[2] = X.ptr[5] + X.ptr[3];
	Y.ptr[3] = X.ptr[6] + X.ptr[1];
	Y.ptr[4] = (X.ptr[3] - X.ptr[5]) / speedOfLight;
	Y.ptr[5] = (X.ptr[4] - X.ptr[2]) / speedOfLight;
	Y.ptr[6] = (X.ptr[7] - X.ptr[0]) / speedOfLight;
	Y.ptr[7] = (X.ptr[6] - X.ptr[1]) * speedOfLight;
   
   <? elseif side==1 then ?>

	Y.ptr[0] = X.ptr[4] + X.ptr[2];
	Y.ptr[1] = X.ptr[7] + X.ptr[0];
	Y.ptr[2] = X.ptr[5] + X.ptr[3];
	Y.ptr[3] = (X.ptr[5] - X.ptr[3]) / speedOfLight;
	Y.ptr[4] = X.ptr[6] + X.ptr[1];
	Y.ptr[5] = (X.ptr[2] - X.ptr[4]) / speedOfLight;
	Y.ptr[6] = (X.ptr[7] - X.ptr[0]) / speedOfLight;
	Y.ptr[7] = (X.ptr[6] - X.ptr[1]) * speedOfLight;
   
   <? elseif side==2 then ?>

	Y.ptr[0] = X.ptr[4] + X.ptr[2];
	Y.ptr[1] = X.ptr[5] + X.ptr[3];
	Y.ptr[2] = X.ptr[7] + X.ptr[0];
	Y.ptr[3] = (X.ptr[3] - X.ptr[5]) / speedOfLight;
	Y.ptr[4] = (X.ptr[4] - X.ptr[2]) / speedOfLight;
	Y.ptr[5] = X.ptr[6] + X.ptr[1];
	Y.ptr[6] = (X.ptr[7] - X.ptr[0]) / speedOfLight;
	Y.ptr[7] = (X.ptr[6] - X.ptr[1]) * speedOfLight;
	
	<? end ?>
	
	for (int i = 8; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	return fluxFromCons_<?=side?>(X, x);
}

<? end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	real3 mu0_J = real3_scale(U->E, mu0 / U->conductivity);
	deriv->E = real3_sub(deriv->E, mu0_J);
	deriv->phi += U->rhoCharge * divPhiWavespeed / eps0;
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){};
}
<? end ?>
