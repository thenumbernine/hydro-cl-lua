#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

//#define SUPPORT_CURVILINEAR

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
#if !defined(SUPPORT_CURVILINEAR)	//cartesian
	real lambda = 1. / sqrt(U->eps * U->mu);
	return (range_t){-lambda, lambda};
#else	//grid metric
	real det_gamma = volume_at(x);
	real det_gamma2 = det_gamma * det_gamma;
	real det_gamma3 = det_gamma * det_gamma2;
	
	<? if side == 0 then ?>
	real detg_gUjj = coord_g11(x) * coord_g22(x) - coord_g12(x) * coord_g12(x);
	<? elseif side == 1 then ?>
	real detg_gUjj = coord_g00(x) * coord_g22(x) - coord_g02(x) * coord_g02(x);
	<? elseif side == 2 then ?>
	real detg_gUjj = coord_g00(x) * coord_g11(x) - coord_g01(x) * coord_g01(x);
	<? end ?>

	real lambda = sqrt(detg_gUjj / (det_gamma3 * U->eps * U->mu));
	return (range_t){-lambda, lambda};
#endif
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x
) {
#if !defined(SUPPORT_CURVILINEAR)
	return (<?=eqn.eigen_t?>){
		.sqrt_eps = sqrt(.5 * (UL.eps + UR.eps)),
		.sqrt_mu = sqrt(.5 * (UL.mu + UR.mu)),
	};
#else	
	real3 xR = x;
	xR.s<?=side?> += .5 * grid_dx<?=side?>;
	
	real3 xL = x;
	xL.s<?=side?> -= .5 * grid_dx<?=side?>;


	real det_gammaR = volume_at(xR);
	real det_gammaR2 = det_gammaR * det_gammaR;
	real det_gammaR3 = det_gammaR * det_gammaR2;
	
	real det_gammaL = volume_at(xL);
	real det_gammaL2 = det_gammaL * det_gammaL;
	real det_gammaL3 = det_gammaL * det_gammaL2;


	real eps = .5 * (UL.eps + UR.eps);
	real mu = .5 * (UL.mu + UR.mu);

	<? if side == 0 then ?>
	real detg_gUjj = .5 * (
		coord_g11(xL) * coord_g22(xL) - coord_g12(xL) * coord_g12(xL)
		+ coord_g11(xR) * coord_g22(xR) - coord_g12(xR) * coord_g12(xR)
	);
	<? elseif side == 1 then ?>
	real detg_gUjj = .5 * (
		coord_g00(xL) * coord_g22(xL) - coord_g02(xL) * coord_g02(xL)
		+ coord_g00(xR) * coord_g22(xR) - coord_g02(xR) * coord_g02(xR)
	);
	<? elseif side == 2 then ?>
	real detg_gUjj = .5 * (
		coord_g00(xL) * coord_g11(xL) - coord_g01(xL) * coord_g01(xL)
		+ coord_g00(xR) * coord_g11(xR) - coord_g01(xR) * coord_g01(xR)
	);
	<? end ?>
	
	real det_gamma3 = sqrt(det_gammaL3 * det_gammaR3);
	real lambda = 1. / sqrt(detg_gUjj / (det_gamma3 * eps * mu));
	
	return (<?=eqn.eigen_t?>){
		.lambda = lambda,
		.sqrt_eps = sqrt(eps),
		.sqrt_mu = sqrt(mu),
	};
#endif
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
		
		<?=solver:getULRCode()?>
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		int indexInt = side + dim * index;	
		eigenBuf[indexInt] = eigen_forSide_<?=side?>(*UL, *UR, xInt);
	}<? end ?>
}
		
/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/
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
	real3 x = cell_x(i);
	
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(U->epsE, 1. / U->eps * U->sigma));

<? if not require 'coord.cartesian'.is(solver.coord) then ?>
	//grid coordinate connection coefficient source terms for contravariant representation
	//epsE^i += 2 / (mu_0 sqrt(g)) [ijk] Gamma_jkl B^l
	//B^i -= 2 / (eps_0 sqrt(g)) [ijk] Gamma_jkl epsE^l
	real sqrt_det_g = sqrt_det_g_grid(x);
	_3sym3 conn = coord_conn(x);

	deriv->epsE.x += 2. / (U->mu * sqrt_det_g) * (
		(conn.y.xz - conn.z.xy) * U->B.x
		+ (conn.y.yz - conn.z.yy) * U->B.y
		+ (conn.y.zz - conn.z.yz) * U->B.z);
	deriv->epsE.y += 2. / (U->mu * sqrt_det_g) * (
		(conn.z.xx - conn.x.xz) * U->B.x
		+ (conn.z.xy - conn.x.yz) * U->B.y
		+ (conn.z.xz - conn.x.zz) * U->B.z);
	deriv->epsE.z += 2. / (U->mu * sqrt_det_g) * (
		(conn.x.xy - conn.y.xx) * U->B.x
		+ (conn.x.yy - conn.y.xy) * U->B.y
		+ (conn.x.yz - conn.y.xz) * U->B.z);
	
	deriv->B.x -= 2. / (U->eps * sqrt_det_g) * (
		(conn.y.xz - conn.z.xy) * U->epsE.x
		+ (conn.y.yz - conn.z.yy) * U->epsE.y
		+ (conn.y.zz - conn.z.yz) * U->epsE.z);
	deriv->B.y -= 2. / (U->eps * sqrt_det_g) * (
		(conn.z.xx - conn.x.xz) * U->epsE.x
		+ (conn.z.xy - conn.x.yz) * U->epsE.y
		+ (conn.z.xz - conn.x.zz) * U->epsE.z);
	deriv->B.z -= 2. / (U->eps * sqrt_det_g) * (
		(conn.x.xy - conn.y.xx) * U->epsE.x
		+ (conn.x.yy - conn.y.xy) * U->epsE.y
		+ (conn.x.yz - conn.y.xz) * U->epsE.z);
<? end ?>
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.eigen_t?> eig;
	eig.sqrt_eps = sqrt(U.eps);
	eig.sqrt_mu = sqrt(U.mu);
	eig.lambda = 1. / (eig.sqrt_eps * eig.sqrt_mu);
	return eig;
}
<? end ?>
