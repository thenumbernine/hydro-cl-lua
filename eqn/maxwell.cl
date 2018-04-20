#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 B = U.B;
	real3 epsE = U.epsE;
	real mu = U.mu;
	real eps = U.eps;
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.epsE = _real3(0., B.z / mu, -B.y / mu),
		.B = _real3(0., -epsE.z / eps, epsE.y / eps),
	<? elseif side == 1 then ?>
		.epsE = _real3(-B.z / mu, 0., B.x / mu),
		.B = _real3(epsE.z / eps, 0., -epsE.x / eps),
	<? elseif side == 2 then ?>
		.epsE = _real3(B.y / mu, -B.x / mu, 0.),
		.B = _real3(-epsE.y / eps, epsE.x / eps, 0.),
	<? end ?>
		.BPot = 0.,
		.sigma = 0.,
		.eps = 0.,
		.mu = 0.,
	};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	real lambda = 1. / sqrt(U->eps * U->mu);
	return (range_t){-lambda, lambda};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_eps = sqrt(.5 * (UL->eps + UR->eps)),
		.sqrt_mu = sqrt(.5 * (UL->mu + UR->mu)),
	};
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

	const real ise = sqrt_1_2 / eig.sqrt_eps;
	const real isu = sqrt_1_2 / eig.sqrt_mu;

	<? if side==0 then ?>
	
	Y.ptr[0] = X.ptr[2] *  ise + X.ptr[4] * isu;
	Y.ptr[1] = X.ptr[1] * -ise + X.ptr[5] * isu;
	Y.ptr[2] = X.ptr[0] * -ise + X.ptr[3] * isu;
	Y.ptr[3] = X.ptr[0] *  ise + X.ptr[3] * isu;
	Y.ptr[4] = X.ptr[1] *  ise + X.ptr[5] * isu;
	Y.ptr[5] = X.ptr[2] * -ise + X.ptr[4] * isu;
	
	<? elseif side==1 then ?>
	
	Y.ptr[0] = X.ptr[0] *  ise + X.ptr[5] * isu;
	Y.ptr[1] = X.ptr[2] * -ise + X.ptr[3] * isu;
	Y.ptr[2] = X.ptr[1] * -ise + X.ptr[4] * isu;
	Y.ptr[3] = X.ptr[1] *  ise + X.ptr[4] * isu;
	Y.ptr[4] = X.ptr[2] *  ise + X.ptr[3] * isu;
	Y.ptr[5] = X.ptr[0] * -ise + X.ptr[5] * isu;
	
	<? elseif side==2 then ?>
	
	Y.ptr[0] = X.ptr[1] *  ise + X.ptr[3] * isu;
	Y.ptr[1] = X.ptr[0] * -ise + X.ptr[4] * isu;
	Y.ptr[2] = X.ptr[2] * -ise + X.ptr[5] * isu;
	Y.ptr[3] = X.ptr[2] *  ise + X.ptr[5] * isu;
	Y.ptr[4] = X.ptr[0] *  ise + X.ptr[4] * isu;
	Y.ptr[5] = X.ptr[1] * -ise + X.ptr[3] * isu;
	
	<? end ?>

	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;
	
	const real se = sqrt_1_2 * eig.sqrt_eps;
	const real su = sqrt_1_2 * eig.sqrt_mu;

	<? if side==0 then ?>
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
	Y.ptr[0] = se * (-X.ptr[2] + X.ptr[3]);
	Y.ptr[1] = se * (-X.ptr[1] + X.ptr[4]);
	Y.ptr[2] = se * (X.ptr[0] + -X.ptr[5]);
	Y.ptr[3] = su * (X.ptr[2] + X.ptr[3]);
	Y.ptr[4] = su * (X.ptr[0] + X.ptr[5]);
	Y.ptr[5] = su * (X.ptr[1] + X.ptr[4]);
	
	<? elseif side==1 then ?>

/*
x, -z, -y, y, z, -x
z,  x,  y, y, x,  z

1  0  0 0 0 -1
0  0 -1 1 0  0
0 -1  0 0 1  0
0  1  0 0 1  0
0  0  1 1 0  0
1  0  0 0 0  1
*/
	Y.ptr[0] = se * (X.ptr[0] - X.ptr[5]);
	Y.ptr[1] = se * (-X.ptr[2] + X.ptr[3]);
	Y.ptr[2] = se * (-X.ptr[1] + X.ptr[4]);
	Y.ptr[3] = su * (X.ptr[1] + X.ptr[4]);
	Y.ptr[4] = su * (X.ptr[2] + X.ptr[3]);
	Y.ptr[5] = su * (X.ptr[0] + X.ptr[5]);
	
	<? elseif side==2 then ?>

/*
y, -x, -z, z, x, -y
x,  y,  z, z,  y,  x

0 -1  0 0 1  0
1  0  0 0 0 -1
0  0 -1 1 0  0
1  0  0 0 0  1
0  1  0 0 1  0
0  0  1 1 0  0
*/
	Y.ptr[0] = se * (-X.ptr[1] + X.ptr[4]);
	Y.ptr[1] = se * (X.ptr[0] - X.ptr[5]);
	Y.ptr[2] = se * (-X.ptr[2] + X.ptr[3]);
	Y.ptr[3] = su * (X.ptr[0] + X.ptr[5]);
	Y.ptr[4] = su * (X.ptr[1] + X.ptr[4]);
	Y.ptr[5] = su * (X.ptr[2] + X.ptr[3]);
	
	<? end ?>
	
	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;

	real3 epsE = X.epsE;
	real3 B = X.B;
	real ieps = 1. / (eig.sqrt_eps * eig.sqrt_eps);
	real imu = 1. / (eig.sqrt_mu * eig.sqrt_mu);

	<? if side==0 then ?>
	
	Y.ptr[0] = 0;
	Y.ptr[1] = B.z * imu;
	Y.ptr[2] = -B.y * imu;
	Y.ptr[3] = 0;
	Y.ptr[4] = -epsE.z * ieps;
	Y.ptr[5] = epsE.y * ieps;

	<? elseif side==1 then ?>
		
	Y.ptr[0] = -B.z * imu;
	Y.ptr[1] = 0;
	Y.ptr[2] = B.x * imu;
	Y.ptr[3] = epsE.z * ieps;
	Y.ptr[4] = 0;
	Y.ptr[5] = -epsE.x * ieps;
		
	<? elseif side==2 then ?>
		
	Y.ptr[0] = B.y * imu;
	Y.ptr[1] = -B.x * imu;
	Y.ptr[2] = 0;
	Y.ptr[3] = -epsE.y * ieps;
	Y.ptr[4] = epsE.x * ieps;
	Y.ptr[5] = 0;
		
	<? end ?>
	
	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}
<? end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(U->epsE, 1. / U->eps * U->sigma));
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_eps = sqrt(U.eps),
		.sqrt_mu = sqrt(U.mu),
	};
}
<? end ?>

<?=eqn.cons_t?> apply_dU_dW(
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) { return W; }

<?=eqn.cons_t?> apply_dW_dU(
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) { return U; }
