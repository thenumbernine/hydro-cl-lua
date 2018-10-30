<?
local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym

local derivOrder = 2 * solver.numGhost
local makePartials = require 'eqn.makepartial'
local makePartial = function(...) return makePartials.makePartial(derivOrder, solver, ...) end
local makePartial2 = function(...) return makePartials.makePartial2(derivOrder, solver, ...) end
?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

kernel void calcDT(
	constant solver_t* solver,
	global real* dtBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const global cons_t* U = UBuf + index;
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma_ij) is only called once
	real f = calc_f(U->alpha);
	real det_gamma = sym3_det(U->gamma_ll);
	real sqrt_f = sqrt(f);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side==0 then ?>
		real gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
		<? elseif side==1 then ?>
		real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
		<? elseif side==2 then ?>
		real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
		<? end ?>	
		real lambdaLight = U->alpha * sqrt(gammaUjj);
		
		real lambdaGauge = lambdaLight * sqrt_f;
		real lambda = (real)max(lambdaGauge, lambdaLight);

		<? if eqn.useShift ~= 'none' then ?>
		real betaUi = U->beta_u.s<?=side?>;
		<? else ?>
		const real betaUi = 0.;
		<? end ?>
		
		real lambdaMin = (real)min((real)0., -betaUi - lambda);
		real lambdaMax = (real)max((real)0., -betaUi + lambda);

		dt = (real)min((real)dt, (real)(solver->grid_dx.s<?=side?> / (fabs(lambdaMax - lambdaMin) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt; 
}

//used by PLM
<? for side=0,solver.dim-1 do ?>
eigen_t eigen_forCell_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x 
) {
	eigen_t eig;
	eig.alpha = U.alpha;
	eig.sqrt_f = sqrt(calc_f(U.alpha));
	real det_gamma = sym3_det(U.gamma_ll);
	eig.gamma_uu = sym3_inv(U.gamma_ll, det_gamma);
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gamma_uu.xx), sqrt(eig.gamma_uu.yy), sqrt(eig.gamma_uu.zz));
	
	<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = U.beta_u;
	<? end ?>

	return eig;
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global cons_t* U,
	real3 x
) {
	real det_gamma = sym3_det(U->gamma_ll);
	
	<? if side==0 then ?>
	real gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
	<? elseif side==1 then ?>
	real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
	<? elseif side==2 then ?>
	real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
	<? end ?>
	
	real lambdaLight = U->alpha * sqrt(gammaUjj);
	
	real f = calc_f(U->alpha);
	real lambdaGauge = lambdaLight * sqrt(f);

	real lambdaMax = max(lambdaGauge, lambdaLight);
	//= lambdaLight * max(sqrt(f), 1)
	real lambdaMin = -lambdaMin;

	<? if eqn.useShift ~= 'none' then ?>
	lambdaMin -= U->beta_u.s<?=side?>;
	lambdaMax -= U->beta_u.s<?=side?>;
	<? end ?>

	return (range_t){
		.min = lambdaMin, 
		.max = lambdaMax,
	};
}
<? end ?>

//used for interface eigen basis
eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3 n
) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	sym3 avg_gamma = (sym3){
		.xx = .5 * (UL.gamma_ll.xx + UR.gamma_ll.xx),
		.xy = .5 * (UL.gamma_ll.xy + UR.gamma_ll.xy),
		.xz = .5 * (UL.gamma_ll.xz + UR.gamma_ll.xz),
		.yy = .5 * (UL.gamma_ll.yy + UR.gamma_ll.yy),
		.yz = .5 * (UL.gamma_ll.yz + UR.gamma_ll.yz),
		.zz = .5 * (UL.gamma_ll.zz + UR.gamma_ll.zz),
	};
	real det_avg_gamma = sym3_det(avg_gamma);

	eigen_t eig;
	eig.alpha = alpha;
	eig.sqrt_f = sqrt(calc_f(alpha));
	eig.gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);
	eig.sqrt_gammaUjj.x = sqrt(eig.gamma_uu.xx);
	eig.sqrt_gammaUjj.y = sqrt(eig.gamma_uu.yy);
	eig.sqrt_gammaUjj.z = sqrt(eig.gamma_uu.zz);
	
	<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = real3_real_mul(real3_add(UL.beta_u, UR.beta_u), .5);
	<? end ?>

	return eig;
}

<? for side=0,solver.dim-1 do ?>
waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x 
) {
	waves_t results;

	real3 a_l = real3_swap<?=side?>(inputU.a_l);							//0-2
	_3sym3 d_lll = (_3sym3){
		.x = sym3_swap<?=side?>(inputU.d_lll.v<?=side?>),					//3-8
		.y = sym3_swap<?=side?>(inputU.d_lll.v<?=side==1 and 0 or 1?>),		//9-14
		.z = sym3_swap<?=side?>(inputU.d_lll.v<?=side==2 and 0 or 2?>),		//15-20
	};
	sym3 K_ll = sym3_swap<?=side?>(inputU.K_ll);							//21-26
	real Theta = inputU.Theta;												//27
	real3 Z_l = real3_swap<?=side?>(inputU.Z_l);							//28-30
	
	sym3 gamma_ll = sym3_swap<?=side?>(eig.gamma_ll);
	sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);

	real sqrt_gammaUUxx = sqrt(gamma_uu.xx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;
	real gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;
	
	// but what is lambda?
	real lambda_1 = (2. - lambda) / (f - 1);
	real lambda_2 = (2. * f - lambda) / (f - 1);

	results.ptr[0] = (((lambda_2 * gammaUxx * Z_l.x) - (lambda_2 * gammaUxx * gammaUyy * d_lll.x.yy)) + ((lambda_2 * gammaUxx * gammaUyy * d_lll.y.xy) - (2. * lambda_2 * gammaUxx * gammaUyz * d_lll.x.yz)) + (lambda_2 * gammaUxx * gammaUyz * d_lll.y.xz) + ((lambda_2 * gammaUxx * gammaUyz * d_lll.z.xy) - (lambda_2 * gammaUxx * gammaUzz * d_lll.x.zz)) + (lambda_2 * gammaUxx * gammaUzz * d_lll.z.xz) + ((lambda_2 * (gammaUxy ^ 2.) * d_lll.x.yy) - (lambda_2 * (gammaUxy ^ 2.) * d_lll.y.xy)) + (lambda_2 * gammaUxy * Z_l.y) + ((((2. * lambda_2 * gammaUxy * gammaUxz * d_lll.x.yz) - (lambda_2 * gammaUxy * gammaUxz * d_lll.y.xz)) - (lambda_2 * gammaUxy * gammaUxz * d_lll.z.xy)) - (lambda_2 * gammaUxy * gammaUyz * d_lll.y.yz)) + ((lambda_2 * gammaUxy * gammaUyz * d_lll.z.yy) - (lambda_2 * gammaUxy * gammaUzz * d_lll.y.zz)) + (lambda_2 * gammaUxy * gammaUzz * d_lll.z.yz) + ((lambda_2 * (gammaUxz ^ 2.) * d_lll.x.zz) - (lambda_2 * (gammaUxz ^ 2.) * d_lll.z.xz)) + (lambda_2 * gammaUxz * Z_l.z) + ((lambda_2 * gammaUxz * gammaUyy * d_lll.y.yz) - (lambda_2 * gammaUxz * gammaUyy * d_lll.z.yy)) + ((lambda_2 * gammaUxz * gammaUyz * d_lll.y.zz) - (lambda_2 * gammaUxz * gammaUyz * d_lll.z.yz)) + (sqrt_f * lambda_1 * Theta) + (sqrt_f * gammaUxx * K_ll.xx) + (2. * sqrt_f * gammaUxy * K_ll.xy) + (2. * sqrt_f * gammaUxz * K_ll.xz) + (sqrt_f * gammaUyy * K_ll.yy) + (2. * sqrt_f * gammaUyz * K_ll.yz) + ((((sqrt_f * gammaUzz * K_ll.zz) - (gammaUxx * a.x)) - (gammaUxy * a.y)) - (gammaUxz * a.z)));
	results.ptr[1] = (Theta + ((gammaUxx * Z_l.x) - (gammaUxx * gammaUyy * d_lll.x.yy)) + ((gammaUxx * gammaUyy * d_lll.y.xy) - (2. * gammaUxx * gammaUyz * d_lll.x.yz)) + (gammaUxx * gammaUyz * d_lll.y.xz) + ((gammaUxx * gammaUyz * d_lll.z.xy) - (gammaUxx * gammaUzz * d_lll.x.zz)) + (gammaUxx * gammaUzz * d_lll.z.xz) + (((gammaUxy ^ 2.) * d_lll.x.yy) - ((gammaUxy ^ 2.) * d_lll.y.xy)) + (gammaUxy * Z_l.y) + ((((2. * gammaUxy * gammaUxz * d_lll.x.yz) - (gammaUxy * gammaUxz * d_lll.y.xz)) - (gammaUxy * gammaUxz * d_lll.z.xy)) - (gammaUxy * gammaUyz * d_lll.y.yz)) + ((gammaUxy * gammaUyz * d_lll.z.yy) - (gammaUxy * gammaUzz * d_lll.y.zz)) + (gammaUxy * gammaUzz * d_lll.z.yz) + (((gammaUxz ^ 2.) * d_lll.x.zz) - ((gammaUxz ^ 2.) * d_lll.z.xz)) + (gammaUxz * Z_l.z) + ((gammaUxz * gammaUyy * d_lll.y.yz) - (gammaUxz * gammaUyy * d_lll.z.yy)) + ((gammaUxz * gammaUyz * d_lll.y.zz) - (gammaUxz * gammaUyz * d_lll.z.yz)));
	results.ptr[2] = ((-(((a.y - (2. * K_ll.xy)) - (2. * Z_l.y)) + ((gammaUxx * d_lll.y.xx) - (2. * gammaUxy * d_lll.x.yy)) + ((2. * gammaUxy * d_lll.y.xy) - (2. * gammaUxz * d_lll.x.yz)) + (((2. * gammaUxz * d_lll.y.xz) - (gammaUyy * d_lll.y.yy)) - (2. * gammaUyz * d_lll.z.yy)) + ((gammaUzz * d_lll.y.zz) - (2. * gammaUzz * d_lll.z.yz)))) / 2.);
	results.ptr[3] = ((-(((a.z - (2. * K_ll.xz)) - (2. * Z_l.z)) + ((gammaUxx * d_lll.z.xx) - (2. * gammaUxy * d_lll.x.yz)) + ((2. * gammaUxy * d_lll.z.xy) - (2. * gammaUxz * d_lll.x.zz)) + ((2. * gammaUxz * d_lll.z.xz) - (2. * gammaUyy * d_lll.y.yz)) + (((gammaUyy * d_lll.z.yy) - (2. * gammaUyz * d_lll.y.zz)) - (gammaUzz * d_lll.z.zz)))) / 2.);
	results.ptr[4] = (((K_ll.yy - (gammaUxx * d_lll.x.yy)) - (gammaUxy * d_lll.y.yy)) - (gammaUxz * d_lll.z.yy));
	results.ptr[5] = (((K_ll.yz - (gammaUxx * d_lll.x.yz)) - (gammaUxy * d_lll.y.yz)) - (gammaUxz * d_lll.z.yz));
	results.ptr[6] = (((K_ll.zz - (gammaUxx * d_lll.x.zz)) - (gammaUxy * d_lll.y.zz)) - (gammaUxz * d_lll.z.zz));
	results.ptr[7] = a.y;
	results.ptr[8] = a.z;
	results.ptr[9] = d_lll.y.xx;
	results.ptr[10] = d_lll.y.xy;
	results.ptr[11] = d_lll.y.xz;
	results.ptr[12] = d_lll.y.yy;
	results.ptr[13] = d_lll.y.yz;
	results.ptr[14] = d_lll.y.zz;
	results.ptr[15] = d_lll.z.xx;
	results.ptr[16] = d_lll.z.xy;
	results.ptr[17] = d_lll.z.xz;
	results.ptr[18] = d_lll.z.yy;
	results.ptr[19] = d_lll.z.yz;
	results.ptr[20] = d_lll.z.zz;
	results.ptr[21] = (a.x + ((Z_l.x - (f * gammaUxx * d_lll.x.xx)) - (gammaUxy * d_lll.x.xy)) + (((gammaUxy * d_lll.y.xx) - (2. * gammaUxy * f * d_lll.x.xy)) - (gammaUxz * d_lll.x.xz)) + (((gammaUxz * d_lll.z.xx) - (2. * gammaUxz * f * d_lll.x.xz)) - (gammaUyy * d_lll.x.yy)) + (((gammaUyy * d_lll.y.xy) - (gammaUyy * f * d_lll.x.yy)) - (2. * gammaUyz * d_lll.x.yz)) + (gammaUyz * d_lll.y.xz) + (((gammaUyz * d_lll.z.xy) - (2. * gammaUyz * f * d_lll.x.yz)) - (gammaUzz * d_lll.x.zz)) + ((gammaUzz * d_lll.z.xz) - (gammaUzz * f * d_lll.x.zz)));
	results.ptr[22] = (a.y + (Z_l.y - (f * gammaUyy * d_lll.y.yy)) + (((gammaUxx * d_lll.x.xy) - (gammaUxx * d_lll.y.xx)) - (gammaUxx * f * d_lll.y.xx)) + (((gammaUxy * d_lll.x.yy) - (gammaUxy * d_lll.y.xy)) - (2. * gammaUxy * f * d_lll.y.xy)) + ((gammaUxz * d_lll.x.yz) - (2. * gammaUxz * d_lll.y.xz)) + (((gammaUxz * d_lll.z.xy) - (2. * gammaUxz * f * d_lll.y.xz)) - (gammaUyz * d_lll.y.yz)) + (((gammaUyz * d_lll.z.yy) - (2. * gammaUyz * f * d_lll.y.yz)) - (gammaUzz * d_lll.y.zz)) + ((gammaUzz * d_lll.z.yz) - (gammaUzz * f * d_lll.y.zz)));
	results.ptr[23] = (a.z + (Z_l.z - (f * gammaUzz * d_lll.z.zz)) + (((gammaUxx * d_lll.x.xz) - (gammaUxx * d_lll.z.xx)) - (gammaUxx * f * d_lll.z.xx)) + (gammaUxy * d_lll.x.yz) + (((gammaUxy * d_lll.y.xz) - (2. * gammaUxy * d_lll.z.xy)) - (2. * gammaUxy * f * d_lll.z.xy)) + (((gammaUxz * d_lll.x.zz) - (gammaUxz * d_lll.z.xz)) - (2. * gammaUxz * f * d_lll.z.xz)) + (((gammaUyy * d_lll.y.yz) - (gammaUyy * d_lll.z.yy)) - (gammaUyy * f * d_lll.z.yy)) + (((gammaUyz * d_lll.y.zz) - (gammaUyz * d_lll.z.yz)) - (2. * gammaUyz * f * d_lll.z.yz)));
	results.ptr[24] = ((a.y + ((2. * K_ll.xy) - (2. * Z_l.y)) + ((gammaUxx * d_lll.y.xx) - (2. * gammaUxy * d_lll.x.yy)) + ((2. * gammaUxy * d_lll.y.xy) - (2. * gammaUxz * d_lll.x.yz)) + (((2. * gammaUxz * d_lll.y.xz) - (gammaUyy * d_lll.y.yy)) - (2. * gammaUyz * d_lll.z.yy)) + ((gammaUzz * d_lll.y.zz) - (2. * gammaUzz * d_lll.z.yz))) / 2.);
	results.ptr[25] = ((a.z + ((2. * K_ll.xz) - (2. * Z_l.z)) + ((gammaUxx * d_lll.z.xx) - (2. * gammaUxy * d_lll.x.yz)) + ((2. * gammaUxy * d_lll.z.xy) - (2. * gammaUxz * d_lll.x.zz)) + ((2. * gammaUxz * d_lll.z.xz) - (2. * gammaUyy * d_lll.y.yz)) + (((gammaUyy * d_lll.z.yy) - (2. * gammaUyz * d_lll.y.zz)) - (gammaUzz * d_lll.z.zz))) / 2.);
	results.ptr[26] = (K_ll.yy + (gammaUxx * d_lll.x.yy) + (gammaUxy * d_lll.y.yy) + (gammaUxz * d_lll.z.yy));
	results.ptr[27] = (K_ll.yz + (gammaUxx * d_lll.x.yz) + (gammaUxy * d_lll.y.yz) + (gammaUxz * d_lll.z.yz));
	results.ptr[28] = (K_ll.zz + (gammaUxx * d_lll.x.zz) + (gammaUxy * d_lll.y.zz) + (gammaUxz * d_lll.z.zz));
	results.ptr[29] = ((Theta - (gammaUxx * Z_l.x)) + ((gammaUxx * gammaUyy * d_lll.x.yy) - (gammaUxx * gammaUyy * d_lll.y.xy)) + (((2. * gammaUxx * gammaUyz * d_lll.x.yz) - (gammaUxx * gammaUyz * d_lll.y.xz)) - (gammaUxx * gammaUyz * d_lll.z.xy)) + (((gammaUxx * gammaUzz * d_lll.x.zz) - (gammaUxx * gammaUzz * d_lll.z.xz)) - ((gammaUxy ^ 2.) * d_lll.x.yy)) + ((((gammaUxy ^ 2.) * d_lll.y.xy) - (gammaUxy * Z_l.y)) - (2. * gammaUxy * gammaUxz * d_lll.x.yz)) + (gammaUxy * gammaUxz * d_lll.y.xz) + (gammaUxy * gammaUxz * d_lll.z.xy) + ((gammaUxy * gammaUyz * d_lll.y.yz) - (gammaUxy * gammaUyz * d_lll.z.yy)) + (((gammaUxy * gammaUzz * d_lll.y.zz) - (gammaUxy * gammaUzz * d_lll.z.yz)) - ((gammaUxz ^ 2.) * d_lll.x.zz)) + ((((gammaUxz ^ 2.) * d_lll.z.xz) - (gammaUxz * Z_l.z)) - (gammaUxz * gammaUyy * d_lll.y.yz)) + ((gammaUxz * gammaUyy * d_lll.z.yy) - (gammaUxz * gammaUyz * d_lll.y.zz)) + (gammaUxz * gammaUyz * d_lll.z.yz));
	results.ptr[30] = (-(((lambda_2 * gammaUxx * Z_l.x) - (lambda_2 * gammaUxx * gammaUyy * d_lll.x.yy)) + ((lambda_2 * gammaUxx * gammaUyy * d_lll.y.xy) - (2. * lambda_2 * gammaUxx * gammaUyz * d_lll.x.yz)) + (lambda_2 * gammaUxx * gammaUyz * d_lll.y.xz) + ((lambda_2 * gammaUxx * gammaUyz * d_lll.z.xy) - (lambda_2 * gammaUxx * gammaUzz * d_lll.x.zz)) + (lambda_2 * gammaUxx * gammaUzz * d_lll.z.xz) + ((lambda_2 * (gammaUxy ^ 2.) * d_lll.x.yy) - (lambda_2 * (gammaUxy ^ 2.) * d_lll.y.xy)) + (lambda_2 * gammaUxy * Z_l.y) + ((((2. * lambda_2 * gammaUxy * gammaUxz * d_lll.x.yz) - (lambda_2 * gammaUxy * gammaUxz * d_lll.y.xz)) - (lambda_2 * gammaUxy * gammaUxz * d_lll.z.xy)) - (lambda_2 * gammaUxy * gammaUyz * d_lll.y.yz)) + ((lambda_2 * gammaUxy * gammaUyz * d_lll.z.yy) - (lambda_2 * gammaUxy * gammaUzz * d_lll.y.zz)) + (lambda_2 * gammaUxy * gammaUzz * d_lll.z.yz) + ((lambda_2 * (gammaUxz ^ 2.) * d_lll.x.zz) - (lambda_2 * (gammaUxz ^ 2.) * d_lll.z.xz)) + (lambda_2 * gammaUxz * Z_l.z) + ((lambda_2 * gammaUxz * gammaUyy * d_lll.y.yz) - (lambda_2 * gammaUxz * gammaUyy * d_lll.z.yy)) + ((((((((((((lambda_2 * gammaUxz * gammaUyz * d_lll.y.zz) - (lambda_2 * gammaUxz * gammaUyz * d_lll.z.yz)) - (sqrt_f * lambda_1 * Theta)) - (sqrt_f * gammaUxx * K_ll.xx)) - (2. * sqrt_f * gammaUxy * K_ll.xy)) - (2. * sqrt_f * gammaUxz * K_ll.xz)) - (sqrt_f * gammaUyy * K_ll.yy)) - (2. * sqrt_f * gammaUyz * K_ll.yz)) - (sqrt_f * gammaUzz * K_ll.zz)) - (gammaUxx * a.x)) - (gammaUxy * a.y)) - (gammaUxz * a.z))));


	return results;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t input,
	real3 x
) {
	cons_t resultU;
	for (int j = 0; j < numStates; ++j) {
		resultU.ptr[j] = 0;
	}
	
	sym3 gamma_ll = sym3_swap<?=side?>(eig.gamma_ll);
	sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);
		
	real sqrt_gammaUUxx = sqrt(gamma_uu.xx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;
	real gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;
	
	results[0] = ((((input[0] - input[30]) - (\lambda_2 * input[1])) + (\lambda_2 * input[29]) + (2. * gammaUxy * input[7]) + (2. * gammaUxz * input[8])) / (-(2. * gammaUxx)));
	results[1] = input[7];
	results[2] = input[8];
	results[3] = ((((2. * input[0]) - (2. * input[1])) + (4. * input[21] * gammaUxx) + (((2. * input[29]) - (2. * input[30])) - (2. * \lambda_2 * input[1])) + (((2. * \lambda_2 * input[29]) - (4. * gammaUxy * input[2] * f)) - (12. * gammaUxy * input[7] * f)) + (8. * gammaUxy * input[9] * gammaUxx * f) + (8. * (gammaUxy ^ 2.) * input[10] * f) + (4. * gammaUxy * input[22]) + (4. * gammaUxy * input[24] * f) + (8. * gammaUxy * (f ^ 2.) * input[9] * gammaUxx) + (16. * (gammaUxy ^ 2.) * (f ^ 2.) * input[10]) + (8. * gammaUxy * f * input[22]) + (8. * gammaUxy * gammaUxz * input[11] * f) + (8. * gammaUxy * gammaUxz * input[16] * f) + (16. * gammaUxy * gammaUxz * (f ^ 2.) * input[11]) + (16. * gammaUxy * gammaUxz * (f ^ 2.) * input[16]) + (8. * gammaUxy * gammaUyz * input[13] * f) + (((16. * gammaUxy * gammaUyz * (f ^ 2.) * input[13]) - (4. * gammaUxz * input[3] * f)) - (12. * gammaUxz * input[8] * f)) + (8. * gammaUxz * input[15] * gammaUxx * f) + (8. * (gammaUxz ^ 2.) * input[17] * f) + (4. * gammaUxz * input[23]) + (4. * gammaUxz * input[25] * f) + (8. * gammaUxz * (f ^ 2.) * input[15] * gammaUxx) + (16. * (gammaUxz ^ 2.) * (f ^ 2.) * input[17]) + (8. * gammaUxz * f * input[23]) + (8. * gammaUxz * gammaUyz * input[19] * f) + ((16. * gammaUxz * gammaUyz * (f ^ 2.) * input[19]) - (2. * gammaUyy * input[4] * f)) + (2. * gammaUyy * input[26] * f) + (4. * gammaUyy * gammaUxy * input[12] * f) + (8. * gammaUyy * gammaUxy * (f ^ 2.) * input[12]) + (4. * gammaUyy * gammaUxz * input[18] * f) + ((8. * gammaUyy * gammaUxz * (f ^ 2.) * input[18]) - (4. * gammaUyz * input[5] * f)) + ((4. * gammaUyz * input[27] * f) - (2. * gammaUzz * input[6] * f)) + (2. * gammaUzz * input[28] * f) + (4. * gammaUzz * gammaUxy * input[14] * f) + (8. * gammaUzz * gammaUxy * (f ^ 2.) * input[14]) + (4. * gammaUzz * gammaUxz * input[20] * f) + (8. * gammaUzz * gammaUxz * (f ^ 2.) * input[20])) / (-(4. * gammaUxxSq * f)));
	results[4] = ((-(input[2] + (((((((3. * input[7]) - (input[9] * gammaUxx)) - (2. * input[22])) - input[24]) - (2. * f * input[9] * gammaUxx)) - (4. * gammaUxy * f * input[10])) - (2. * gammaUxz * input[11])) + ((((((((2. * gammaUxz * input[16]) - (4. * gammaUxz * f * input[11])) - (gammaUyy * input[12])) - (2. * gammaUyy * f * input[12])) - (2. * gammaUyz * input[13])) - (4. * gammaUyz * f * input[13])) - (gammaUzz * input[14])) - (2. * gammaUzz * f * input[14])))) / (2. * gammaUxx));
	results[5] = ((-(input[3] + (((((3. * input[8]) - (input[15] * gammaUxx)) - (2. * input[23])) - input[25]) - (2. * f * input[15] * gammaUxx)) + ((((((((((2. * gammaUxy * input[11]) - (2. * gammaUxy * input[16])) - (4. * gammaUxy * f * input[16])) - (4. * gammaUxz * f * input[17])) - (gammaUyy * input[18])) - (2. * gammaUyy * f * input[18])) - (2. * gammaUyz * input[19])) - (4. * gammaUyz * f * input[19])) - (gammaUzz * input[20])) - (2. * gammaUzz * f * input[20])))) / (2. * gammaUxx));
	results[6] = ((-((input[4] - input[26]) + (2. * gammaUxy * input[12]) + (2. * gammaUxz * input[18]))) / (2. * gammaUxx));
	results[7] = ((-((input[5] - input[27]) + (2. * gammaUxy * input[13]) + (2. * gammaUxz * input[19]))) / (2. * gammaUxx));
	results[8] = ((-((input[6] - input[28]) + (2. * gammaUxy * input[14]) + (2. * gammaUxz * input[20]))) / (2. * gammaUxx));
	results[9] = input[9];
	results[10] = input[10];
	results[11] = input[11];
	results[12] = input[12];
	results[13] = input[13];
	results[14] = input[14];
	results[15] = input[15];
	results[16] = input[16];
	results[17] = input[17];
	results[18] = input[18];
	results[19] = input[19];
	results[20] = input[20];
	results[21] = ((input[0] + ((((((((((((input[30] - (\lambda_1 * input[1] * sqrt_f)) - (\lambda_1 * input[29] * sqrt_f)) - (2. * gammaUxy * input[2] * sqrt_f)) - (2. * gammaUxy * input[24] * sqrt_f)) - (2. * gammaUxz * input[3] * sqrt_f)) - (2. * gammaUxz * input[25] * sqrt_f)) - (gammaUyy * input[4] * sqrt_f)) - (gammaUyy * input[26] * sqrt_f)) - (2. * gammaUyz * input[5] * sqrt_f)) - (2. * gammaUyz * input[27] * sqrt_f)) - (gammaUzz * input[6] * sqrt_f)) - (gammaUzz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUxx));
	results[22] = ((input[2] + input[24]) / 2.);
	results[23] = ((input[3] + input[25]) / 2.);
	results[24] = ((input[4] + input[26]) / 2.);
	results[25] = ((input[5] + input[27]) / 2.);
	results[26] = ((input[6] + input[28]) / 2.);
	results[27] = ((input[1] + input[29]) / 2.);
	results[28] = ((((((input[1] - input[29]) - (gammaUxy * input[2])) - (gammaUxy * input[7])) - (gammaUxy * input[9] * gammaUxx)) + (((((gammaUxy * input[24]) - (2. * gammaUxy * gammaUyz * input[13])) - (gammaUxz * input[3])) - (gammaUxz * input[8])) - (gammaUxz * input[15] * gammaUxx)) + (((gammaUxz * input[25]) - (gammaUyy * input[4])) - (2. * gammaUyy * input[10] * gammaUxx)) + ((((((gammaUyy * input[26]) - (gammaUyy * gammaUxy * input[12])) - (gammaUyy * gammaUxz * input[18])) - (2. * gammaUyz * input[5])) - (2. * gammaUyz * input[11] * gammaUxx)) - (2. * gammaUyz * input[16] * gammaUxx)) + ((((2. * gammaUyz * input[27]) - (2. * gammaUyz * gammaUxz * input[19])) - (gammaUzz * input[6])) - (2. * gammaUzz * input[17] * gammaUxx)) + (((gammaUzz * input[28]) - (gammaUzz * gammaUxy * input[14])) - (gammaUzz * gammaUxz * input[20]))) / (2. * gammaUxx));
	results[29] = (((input[2] * gammaUxx) + ((input[7] * gammaUxx) - (input[24] * gammaUxx)) + (((gammaUxxSq * input[9]) - (gammaUxx * gammaUyy * input[12])) - (2. * gammaUxx * gammaUyz * input[18])) + (gammaUxy * input[4]) + (2. * gammaUxy * input[10] * gammaUxx) + ((2. * (gammaUxy ^ 2.) * input[12]) - (gammaUxy * input[26])) + (2. * gammaUxy * gammaUxz * input[13]) + (gammaUxz * input[5]) + (2. * gammaUxz * input[11] * gammaUxx) + ((2. * (gammaUxz ^ 2.) * input[19]) - (gammaUxz * input[27])) + (2. * gammaUxz * gammaUxy * input[18]) + ((gammaUzz * input[14] * gammaUxx) - (2. * gammaUzz * gammaUxx * input[19]))) / (2. * gammaUxx));
	results[30] = (((input[3] * gammaUxx) + ((input[8] * gammaUxx) - (input[25] * gammaUxx)) + ((((gammaUxxSq * input[15]) - (2. * gammaUxx * gammaUyy * input[13])) - (2. * gammaUxx * gammaUyz * input[14])) - (gammaUxx * gammaUzz * input[20])) + (gammaUxy * input[5]) + (2. * (gammaUxy ^ 2.) * input[13]) + ((2. * gammaUxy * input[16] * gammaUxx) - (gammaUxy * input[27])) + (2. * gammaUxy * gammaUxz * input[19]) + (gammaUxz * input[6]) + (2. * gammaUxz * input[17] * gammaUxx) + ((2. * (gammaUxz ^ 2.) * input[20]) - (gammaUxz * input[28])) + (2. * gammaUxz * gammaUxy * input[14]) + (gammaUyy * input[18] * gammaUxx)) / (2. * gammaUxx));

	return resultU;
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x
) {
}
<? end ?>

//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf)
{
	SETBOUNDS_NOGHOST();
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real f = calc_f(U->alpha);	//could be based on alpha...

<? if eqn.useStressEnergyTerms then ?> 
	real3 S_l = sym3_real3_mul(U->gamma_ll, U->S_u);
	sym3 S_ll = U->S_ll;
	real S = sym3_dot(S_ll, gamma_uu);
<? else ?>
	real3 S_l = real3_zero;
	sym3 S_ll = sym3_zero;
	real S = 0.;
<? end ?>

	// source terms
	
	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real tr_K = real3x3_trace(K_ul);								//K^k_k
	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		//KSq_ij = K_ik K^k_j

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};

	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//e_i = d^j_ji
	real3 e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end	?>,
<? end
?>	};

	real3 V_l = real3_sub(d_l, e_l);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d_lll.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	sym3 R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>
			+ conn_ull.<?=xk?>.<?=xij?> * (V_l.<?=xk?> - e_l.<?=xk?>)

<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(j,l)?>
			- 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=l-1?>].<?=xj?>.<?=xk?>
			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=j-1?>].<?=xl?>.<?=xk?>
			+ 2. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=k-1?>].<?=xj?>.<?=xl?>
			- 3. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=j-1?>].<?=xk?>.<?=xl?>
<? 		end
	end
?>		,
<? end
?>	};

	//srcK_ij = 
	
	//lapse:
	//	(-a_i a_j 
	//		+ conn^k_ij a_k

	//Riemann curvature:
	//		+ conn^k_ij (V_k - d^l_lk) 
	//		+ 2 d_ki^l d^k_jl
	//		- 2 d_ki^l d_lj^k
	//		+ 2 d_ki^l d_jl^k
	//		+ 2 d_il^k d_kj^l
	//		- 3 d_il^k d_jk^l
	
	//extrinsic curvature:
	//		+ K K_ij - 2 K_ik K^k_j

	//matter terms:
	//		- 4 pi (2 S_ij - gamma_ij (S - rho))
	
	//...and the shift comes later ...
	
	sym3 srcK_ll_over_alpha = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 
			- U->a_l.<?=xi?> * U->a_l.<?=xj?>
<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * U->a_l.<?=xk?>
<?	end
?>			+ R_ll.<?=xij?>
			+ tr_K * U->K_ll.<?=xij?>
			- 2. * KSq_ll.<?=xij?>
			- (8. * M_PI * S_ll.<?=xij?> - 4. * M_PI * U->gamma_ll.<?=xij?> * (S - U->rho))
		,
<? end
?>	};

	//d_i = d_ij^j
	real3 d_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_trace(d_llu[<?=i-1?>]),
<? end
?>	};

	//alpha_,t = shift terms - alpha^2 f gamma^ij K_ij
	deriv->alpha += -U->alpha * U->alpha * f * tr_K;
	
	//gamma_ij,t = shift terms - 2 alpha K_ij
	sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));
	
	//K_ij,t = shift terms + alpha srcK_ij
	sym3_add(deriv->K_ll, sym3_real_mul(srcK_ll_over_alpha, U->alpha));

	// and now for the first-order constraints

	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	<? for i,xi in ipairs(xNames) do ?>{
		<? if i <= solver.dim then ?>
		real di_alpha = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (2. * solver->grid_dx.s<?=i-1?>);
		<? else ?>
		real di_alpha = 0.;
		<? end ?>
		deriv->a_l.<?=xi?> += solver->a_convCoeff * (di_alpha / U->alpha - U->a_l.<?=xi?>);
	}<? end ?>	
	
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	<? 
for i,xi in ipairs(xNames) do 
	for jk,xjk in ipairs(symNames) do ?>{
		<? if i <= solver.dim then ?>
		real di_gamma_jk = (U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> - U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>) / (2. * solver->grid_dx.s<?=i-1?>);
		<? else ?>
		real di_gamma_jk = 0;
		<? end ?>
		deriv->d_lll.<?=xi?>.<?=xjk?> += solver->d_convCoeff * (.5 * di_gamma_jk - U->d_lll.<?=xi?>.<?=xjk?>);
	}<? 
	end
end ?>
}

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf
) {
	SETBOUNDS(numGhost,numGhost);		
	global cons_t* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real tr_K = real3x3_trace(K_ul);								//K^k_k
	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		//KSq_ij = K_ik K^k_j

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};
	
	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//e_i = d^j_ji
	real3 e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end	?>,
<? end
?>	};
	
	real3 V_l = real3_sub(d_l, e_l);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d_lll.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	sym3 R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>
			+ conn_ull.<?=xk?>.<?=xij?> * (V_l.<?=xk?> - e_l.<?=xk?>)

<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(j,l)?>
			- 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=l-1?>].<?=xj?>.<?=xk?>
			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=j-1?>].<?=xl?>.<?=xk?>
			+ 2. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=k-1?>].<?=xj?>.<?=xl?>
			- 3. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=j-1?>].<?=xk?>.<?=xl?>
<? 		end
	end
?>		,
<? end
?>	};

	//scaled down by 1/8 to match B&S BSSNOK equations ... maybe I'll scale theirs up by 8 ...
	//B&S eqn 2.125 ... divded by two
	//Alcubierre eqn 2.5.9
	//H = 1/2 (R + K^2 - K_ij K^ij) - 8 pi rho
	real R = sym3_dot(R_ll, gamma_uu);
	real tr_KSq = sym3_dot(KSq_ll, gamma_uu);
	U->H = .5 * (R + tr_K * tr_K - tr_KSq) <? 
if eqn.useStressEnergyTerms then ?>
	- 8. * M_PI * U->rho <? 
end ?>;
	//momentum constraint
}
