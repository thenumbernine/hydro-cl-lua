<? if moduleName == nil then ?>
<? elseif moduleName == "sqrt_1_2" then ?>

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? elseif moduleName == "eqn.common" then ?>

static inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) { return U; }
static inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) { return W; }

<? elseif moduleName == "applyInitCond" then ?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	global <?=eqn.cons_t?>* U = UBuf + index;

	//used
	real3 D = real3_zero;
	real3 B = real3_zero;
	real conductivity = 1.;
	
	//natural units say eps0 = 1/4pi, mu0 = 4pi
	//but waves don't make it to the opposite side...
	//mu0 eps0 = 1/c^2
	real permittivity = 1.; //1. / (4. * M_PI);
	real permeability = 1.; //4. * M_PI;
	
	//throw-away
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
	
	<?=code?>
	
	U->D = D;
	U->B = B;
	U->divBPot = 0;
	U->sigma = conductivity;
	U->eps = permittivity;
	U->mu = permeability;
}


<? elseif moduleName == "fluxFromCons" then ?>

#error this is out of date, should be made with normal_t, then moved to fluxFromCons
<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	
	real3 B_u = U.B;
	real3 D_u = U.D;

	//this only works for fluxFromCons
	real3 B_l = sym3_real3_mul(gamma, B_u);
	real3 D_l = sym3_real3_mul(gamma, D_u);
	
	real mu = U.mu;
	real eps = U.eps;

	//TODO update this and addSource to your worksheet
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.D = _real3(0., B_l.z / mu, -B_l.y / mu),
		.B = _real3(0., -D_l.z / eps, D_l.y / eps),
	<? elseif side == 1 then ?>
		.D = _real3(-B_l.z / mu, 0., B_l.x / mu),
		.B = _real3(D_l.z / eps, 0., -D_l.x / eps),
	<? elseif side == 2 then ?>
		.D = _real3(B_l.y / mu, -B_l.x / mu, 0.),
		.B = _real3(-D_l.y / eps, D_l.x / eps, 0.),
	<? end ?>
		.divBPot = 0.,
		.sigma = 0.,
		.eps = 0.,
		.mu = 0.,
	};
}

<? elseif moduleName == "fluxFromCons" then ?>

range_t calcCellMinMaxEigenvalues(
	const global <?=eqn.cons_t?>* U,
	real3 x<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>	
	real det_gamma = sym3_det(gamma);	
	real det_gamma2 = det_gamma * det_gamma;
	real det_gamma3 = det_gamma * det_gamma2;
	
	<? if side == 0 then ?>
	real detg_gUjj = gamma.yy * gamma.zz - gamma.yz * gamma.yz;
	<? elseif side == 1 then ?>
	real detg_gUjj = gamma.xx * gamma.zz - gamma.xz * gamma.xz;
	<? elseif side == 2 then ?>
	real detg_gUjj = gamma.xx * gamma.yy - gamma.xy * gamma.xy;
	<? end ?>

	real lambda = alpha * sqrt(detg_gUjj / (det_gamma3 * U->eps * U->mu));
	return (range_t){-lambda, lambda};
}

<? elseif moduleName == "calcEigenBasis" then ?>

//TODO HLL needs eigen_forInterface 
//but it would have to pass the extra ADM args into it
/*
<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forInterface(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	real3 n
) {
}
<? end ?>
*/

#error calcEigenBasis has been removed, and eigen_t structs are now calculated inline ... soooo ... convert this to something compatible
kernel void calcEigenBasis(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.eigen_t?>* eigenBuf,
	const global <?=solver.getULRArg?><?=
	solver:getADMArgs()?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;
	
	<?=solver:getADMVarCode{suffix='R'} --[[ produce alphaR, betaR, gammaR at indexR ]] ?>
	real det_gammaR = sym3_det(gammaR);
	real det_gammaR2 = det_gammaR * det_gammaR;
	real det_gammaR3 = det_gammaR * det_gammaR2;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - solver->stepsize.s<?=side?>;
		
		<?=solver:getULRCode()?>
		
		<?=solver:getADMVarCode{suffix='L'} --[[ produce alphaL, betaL, gammaL at indexL ]] ?>
		real det_gammaL = sym3_det(gammaL);
		real det_gammaL2 = det_gammaL * det_gammaL;
		real det_gammaL3 = det_gammaL * det_gammaL2;
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		//*eig = eigen_forInterface(*UL, *UR, xInt, normalForSide<?=side?>());
		eig->eps = .5 * (UL->eps + UR->eps);
		eig->mu = .5 * (UL->mu + UR->mu);
		real alpha = .5 * (alphaL + alphaR);
		real det_gamma = .5 * (det_gammaL + det_gammaR);

		<? if side == 0 then ?>
		real detg_gUjj = .5 * (
			gammaL.yy * gammaL.zz - gammaL.yz * gammaL.yz
			+ gammaR.yy * gammaR.zz - gammaR.yz * gammaR.yz
		);
		<? elseif side == 1 then ?>
		real detg_gUjj = .5 * (
			gammaL.xx * gammaL.zz - gammaL.xz * gammaL.xz
			+ gammaR.xx * gammaR.zz - gammaR.xz * gammaR.xz
		);
		<? elseif side == 2 then ?>
		real detg_gUjj = .5 * (
			gammaL.xx * gammaL.yy - gammaL.xy * gammaL.xy
			+ gammaR.xx * gammaR.yy - gammaR.xy * gammaR.xy
		);
		<? end ?>
	
		real det_gamma3 = sqrt(det_gammaL3 * det_gammaR3);
		eig->lambda = alpha / sqrt(detg_gUjj / (det_gamma3 * eig->eps * eig->mu));
	}<? end ?>
}

<? elseif moduleName == "eigen_left/rightTransform" then ?>

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/

<?=eqn.waves_t?> eigen_leftTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> UX,
	real3 x
) {
	const real ise = sqrt_1_2 / sqrt(eig.eps);
	const real isu = sqrt_1_2 / sqrt(eig.mu);

	<?=eqn.waves_t?> UY;
	real* X = UX.ptr;
	real* Y = UY.ptr;

	<? if side == 0 then ?>

	Y[0] = 								X[2] *  ise 				+ X[4] * isu;
	Y[1] = 				X[1] * -ise 												+ X[5] * isu;
	Y[2] = X[0] * -ise 								+ X[3] * isu;
	Y[3] = X[0] *  ise 								+ X[3] * isu;
	Y[4] = 				X[1] *  ise 												+ X[5] * isu;
	Y[5] = 								X[2] * -ise 				+ X[4] * isu;
	
	<? elseif side == 1 then ?>
	
	Y[0] = X[0] *  ise + X[5] * isu;
	Y[1] = X[2] * -ise + X[3] * isu;
	Y[2] = X[1] * -ise + X[4] * isu;
	Y[3] = X[1] *  ise + X[4] * isu;
	Y[4] = X[2] *  ise + X[3] * isu;
	Y[5] = X[0] * -ise + X[5] * isu;
	
	<? elseif side == 2 then ?>
	
	Y[0] = X[1] *  ise + X[3] * isu;
	Y[1] = X[0] * -ise + X[4] * isu;
	Y[2] = X[2] * -ise + X[5] * isu;
	Y[3] = X[2] *  ise + X[5] * isu;
	Y[4] = X[0] *  ise + X[4] * isu;
	Y[5] = X[1] * -ise + X[3] * isu;
	
	<? end ?>

	return UY;
}

<?=eqn.cons_t?> eigen_rightTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> UX,
	real3 x
) {
	const real se = sqrt_1_2 * sqrt(eig.eps * eig.detg_gUjj);
	const real su = sqrt_1_2 * sqrt(eig.mu * eig.detg_gUjj);

	<?=eqn.cons_t?> UY;
	real* X = UX.ptr;
	real* Y = UY.ptr;

	<? if side==0 then ?>
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
	Y[0] = se * (				-X[2] + X[3]					);
	Y[1] = se * (		-X[1] 					+ X[4]			);
	Y[2] = se * (X[0] 									+ -X[5]	);
	Y[3] = su * (				X[2] + X[3]						);
	Y[4] = su * (X[0] 									+ X[5]	);
	Y[5] = su * (		X[1] 					+ X[4]			);
	
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
	Y[0] = se * (X[0] - X[5]);
	Y[1] = se * (-X[2] + X[3]);
	Y[2] = se * (-X[1] + X[4]);
	Y[3] = su * (X[1] + X[4]);
	Y[4] = su * (X[2] + X[3]);
	Y[5] = su * (X[0] + X[5]);
	
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
	Y[0] = se * (-X[1] + X[4]);
	Y[1] = se * (X[0] - X[5]);
	Y[2] = se * (-X[2] + X[3]);
	Y[3] = su * (X[0] + X[5]);
	Y[4] = su * (X[1] + X[4]);
	Y[5] = su * (X[2] + X[3]);
	
	<? end ?>
	
	Y[6] = 0;	//divBPot

	return UY;
}

<? elseif moduleName == "eigen_fluxTransform" then ?>

<?=eqn.cons_t?> eigen_fluxTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> UX,
	real3 x
) {
	//swap input dim x<->side
	real3 D = UX.D;
	real3 B = UX.B;
	real ieps = 1. / eig.eps;
	real imu = 1. / eig.mu;

	<?=eqn.cons_t?> UY;
	real* X = UX.ptr;
	real* Y = UY.ptr;

	<? if side==0 then ?>
	
	Y[0] = 0;
	Y[1] = B.z * imu;
	Y[2] = -B.y * imu;
	Y[3] = 0;
	Y[4] = -D.z * ieps;
	Y[5] = D.y * ieps;

	<? elseif side==1 then ?>
		
	Y[0] = -B.z * imu;
	Y[1] = 0;
	Y[2] = B.x * imu;
	Y[3] = D.z * ieps;
	Y[4] = 0;
	Y[5] = -D.x * ieps;
		
	<? elseif side==2 then ?>
		
	Y[0] = B.y * imu;
	Y[1] = -B.x * imu;
	Y[2] = 0;
	Y[3] = -D.y * ieps;
	Y[4] = D.x * ieps;
	Y[5] = 0;
		
	<? end ?>

	Y[6] = 0.;

	return UY;
}

<? elseif moduleName == "addSource" then ?>

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	deriv->D = real3_sub(deriv->D, real3_real_mul(U->D, 1. / U->eps * U->sigma));
}

<? elseif moduleName == "eigen_forCell" then ?>

//used by PLM

//TODO FINISHME
<?=eqn.eigen_t?> eigen_forCell(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.eigen_t?> eig;
	eig.eps = U.eps;
	eig.mu = U.mu;
	//eig.lambda = eig.alpha / sqrt(eig.detg_gUjj / (det_gamma3 * eig.eps * eig.mu));
	return eig;
}

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
