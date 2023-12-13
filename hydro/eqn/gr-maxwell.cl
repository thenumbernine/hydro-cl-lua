//// MODULE_NAME: sqrt_1_2

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

//// MODULE_NAME: <?=eqn_common?>

static inline <?=prim_t?> primFromCons(<?=cons_t?> U, real3 x) { return U; }
static inline <?=cons_t?> consFromPrim(<?=prim_t?> W, real3 x) { return W; }

//// MODULE_NAME: <?=applyInitCondCell?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool const lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;

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


//// MODULE_NAME: <?=fluxFromCons?>

#error this is out of date, should be made with <?=normal_t?>, and should probably get rid of getADMArgs() somehow
<?=cons_t?> <?=fluxFromCons?>(
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const U,
	<?=cell_t?> const * const cell,
	<?=normal_t?> const n<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	
	real3 B_u = U.B;
	real3 D_u = U.D;

	//this only works for fluxFromCons
	real3 B_l = real3s3_real3_mul(gamma, B_u);
	real3 D_l = real3s3_real3_mul(gamma, D_u);
	
	real mu = U.mu;
	real eps = U.eps;

	//TODO update this and addSource to your worksheet
	return (<?=cons_t?>){
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

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>

<?=range_t?> <?=calcCellMinMaxEigenvalues?>(
	global <?=cons_t?> const * const U,
	real3 const x<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>	
	real det_gamma = real3s3_det(gamma);	
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
	return (<?=range_t?>){-lambda, lambda};
}

//// MODULE_NAME: calcEigenBasis

//TODO HLL needs eigen_forInterface 
//but it would have to pass the extra ADM args into it
#if 0
#define eigen_forInterface(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
)
#endif

#error calcEigenBasis has been removed, and eigen_t structs are now calculated inline ... soooo ... convert this to something compatible
kernel void calcEigenBasis(
	constant <?=solver_t?> const * const solver,
	global <?=eigen_t?> * const eigenBuf,
	global const <?=solver.getULRArg?><?=
	solver:getADMArgs()?>
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost - 1);
	real3 const x = cellBuf[index].pos;
	
	int const indexR = index;
	
	<?=solver:getADMVarCode{suffix='R'} --[[ produce alphaR, betaR, gammaR at indexR ]] ?>
	real det_gammaR = real3s3_det(gammaR);
	real det_gammaR2 = det_gammaR * det_gammaR;
	real det_gammaR3 = det_gammaR * det_gammaR2;

	<? for side=0,solver.dim-1 do ?>{
		int const side = <?=side?>;
		
		int const indexL = index - solver->stepsize.s<?=side?>;
		
		<?=solver:getULRCode()?>
		
		<?=solver:getADMVarCode{suffix='L'} --[[ produce alphaL, betaL, gammaL at indexL ]] ?>
		real det_gammaL = real3s3_det(gammaL);
		real det_gammaL2 = det_gammaL * det_gammaL;
		real det_gammaL3 = det_gammaL * det_gammaL2;
		
		int const indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		
		global <?=eigen_t?>* eig = eigenBuf + indexInt;
		//*eig = eigen_forInterface(eig, solver, UL, UR, cellL, cellR, xInt, normalForSide<?=side?>());
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

//// MODULE_NAME: eigen_leftTransform eigen_rightTransform

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/

<?=waves_t?> eigen_leftTransform(
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const eig,
	<?=cons_t?> const UX,
	real3 const x
) {
	real const ise = sqrt_1_2 / sqrt(eig.eps);
	real const isu = sqrt_1_2 / sqrt(eig.mu);

	<?=waves_t?> UY;
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

<?=cons_t?> eigen_rightTransform(
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const eig,
	<?=waves_t?> const UX,
	real3 const x
) {
	real const se = sqrt_1_2 * sqrt(eig.eps * eig.detg_gUjj);
	real const su = sqrt_1_2 * sqrt(eig.mu * eig.detg_gUjj);

	<?=cons_t?> UY;
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

//// MODULE_NAME: <?=eigen_fluxTransform?>

void <?=eigen_fluxTransform?>(
	<?=cons_t?> * const result,
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const eig,
	<?=cons_t?> const UX,
	<?=cell_t?> const * const cell,
	<?=normal_t?> const n
) {
	//swap input dim x<->side
	real3 const D = UX.D;
	real3 const B = UX.B;
	real const ieps = 1. / eig.eps;
	real const imu = 1. / eig.mu;

	real const * const X = UX.ptr;
	real * const Y = result->ptr;

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
}

//// MODULE_NAME: <?=addSource?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	deriv->D = real3_sub(deriv->D, real3_real_mul(U->D, 1. / U->eps * U->sigma));
}

//// MODULE_NAME: <?=eigen_forCell?>

//used by PLM

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*real3 const */n\
) {\
	(result)->eps = U.eps;\
	(result)->mu = U.mu;\
	/*(result)->lambda = (result)->alpha / sqrt((result)->detg_gUjj / (det_gamma3 * (result)->eps * (result)->mu));*/\
}
