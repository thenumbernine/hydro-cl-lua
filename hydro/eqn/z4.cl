<?
local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym

local derivOrder = 2 * solver.numGhost
local makePartials = require 'hydro.eqn.makepartial'
local makePartial1 = function(...) return makePartials.makePartial1(derivOrder, solver, ...) end
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
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
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
	
	real sqrt_f = eig.sqrt_f;
	real f = sqrt_f * sqrt_f;
	real lambda_1 = (2. - solver->m) / (f - 1);
	real lambda_2 = (2. * f - solver->m) / (f - 1);

	results.ptr[0] = (((lambda_2 * gamma_uu.xx * Z_l.x) - (lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (lambda_2 * gamma_uu.xy * Z_l.y) + ((((2. * lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (lambda_2 * gamma_uu.xz * Z_l.z) + ((lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)) + (sqrt_f * lambda_1 * Theta) + (sqrt_f * gamma_uu.xx * K_ll.xx) + (2. * sqrt_f * gamma_uu.xy * K_ll.xy) + (2. * sqrt_f * gamma_uu.xz * K_ll.xz) + (sqrt_f * gamma_uu.yy * K_ll.yy) + (2. * sqrt_f * gamma_uu.yz * K_ll.yz) + ((((sqrt_f * gamma_uu.zz * K_ll.zz) - (gamma_uu.xx * a_l.x)) - (gamma_uu.xy * a_l.y)) - (gamma_uu.xz * a_l.z)));
	results.ptr[1] = (Theta + ((gamma_uu.xx * Z_l.x) - (gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (gamma_uu.xy * Z_l.y) + ((((2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (gamma_uu.xz * Z_l.z) + ((gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)));
	results.ptr[2] = ((-(((a_l.y - (2. * K_ll.xy)) - (2. * Z_l.y)) + ((gamma_uu.xx * d_lll.y.xx) - (2. * gamma_uu.xy * d_lll.x.yy)) + ((2. * gamma_uu.xy * d_lll.y.xy) - (2. * gamma_uu.xz * d_lll.x.yz)) + (((2. * gamma_uu.xz * d_lll.y.xz) - (gamma_uu.yy * d_lll.y.yy)) - (2. * gamma_uu.yz * d_lll.z.yy)) + ((gamma_uu.zz * d_lll.y.zz) - (2. * gamma_uu.zz * d_lll.z.yz)))) / 2.);
	results.ptr[3] = ((-(((a_l.z - (2. * K_ll.xz)) - (2. * Z_l.z)) + ((gamma_uu.xx * d_lll.z.xx) - (2. * gamma_uu.xy * d_lll.x.yz)) + ((2. * gamma_uu.xy * d_lll.z.xy) - (2. * gamma_uu.xz * d_lll.x.zz)) + ((2. * gamma_uu.xz * d_lll.z.xz) - (2. * gamma_uu.yy * d_lll.y.yz)) + (((gamma_uu.yy * d_lll.z.yy) - (2. * gamma_uu.yz * d_lll.y.zz)) - (gamma_uu.zz * d_lll.z.zz)))) / 2.);
	results.ptr[4] = (((K_ll.yy - (gamma_uu.xx * d_lll.x.yy)) - (gamma_uu.xy * d_lll.y.yy)) - (gamma_uu.xz * d_lll.z.yy));
	results.ptr[5] = (((K_ll.yz - (gamma_uu.xx * d_lll.x.yz)) - (gamma_uu.xy * d_lll.y.yz)) - (gamma_uu.xz * d_lll.z.yz));
	results.ptr[6] = (((K_ll.zz - (gamma_uu.xx * d_lll.x.zz)) - (gamma_uu.xy * d_lll.y.zz)) - (gamma_uu.xz * d_lll.z.zz));
	results.ptr[7] = a_l.y;
	results.ptr[8] = a_l.z;
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
	results.ptr[21] = (a_l.x + ((Z_l.x - (f * gamma_uu.xx * d_lll.x.xx)) - (gamma_uu.xy * d_lll.x.xy)) + (((gamma_uu.xy * d_lll.y.xx) - (2. * gamma_uu.xy * f * d_lll.x.xy)) - (gamma_uu.xz * d_lll.x.xz)) + (((gamma_uu.xz * d_lll.z.xx) - (2. * gamma_uu.xz * f * d_lll.x.xz)) - (gamma_uu.yy * d_lll.x.yy)) + (((gamma_uu.yy * d_lll.y.xy) - (gamma_uu.yy * f * d_lll.x.yy)) - (2. * gamma_uu.yz * d_lll.x.yz)) + (gamma_uu.yz * d_lll.y.xz) + (((gamma_uu.yz * d_lll.z.xy) - (2. * gamma_uu.yz * f * d_lll.x.yz)) - (gamma_uu.zz * d_lll.x.zz)) + ((gamma_uu.zz * d_lll.z.xz) - (gamma_uu.zz * f * d_lll.x.zz)));
	results.ptr[22] = (a_l.y + (Z_l.y - (f * gamma_uu.yy * d_lll.y.yy)) + (((gamma_uu.xx * d_lll.x.xy) - (gamma_uu.xx * d_lll.y.xx)) - (gamma_uu.xx * f * d_lll.y.xx)) + (((gamma_uu.xy * d_lll.x.yy) - (gamma_uu.xy * d_lll.y.xy)) - (2. * gamma_uu.xy * f * d_lll.y.xy)) + ((gamma_uu.xz * d_lll.x.yz) - (2. * gamma_uu.xz * d_lll.y.xz)) + (((gamma_uu.xz * d_lll.z.xy) - (2. * gamma_uu.xz * f * d_lll.y.xz)) - (gamma_uu.yz * d_lll.y.yz)) + (((gamma_uu.yz * d_lll.z.yy) - (2. * gamma_uu.yz * f * d_lll.y.yz)) - (gamma_uu.zz * d_lll.y.zz)) + ((gamma_uu.zz * d_lll.z.yz) - (gamma_uu.zz * f * d_lll.y.zz)));
	results.ptr[23] = (a_l.z + (Z_l.z - (f * gamma_uu.zz * d_lll.z.zz)) + (((gamma_uu.xx * d_lll.x.xz) - (gamma_uu.xx * d_lll.z.xx)) - (gamma_uu.xx * f * d_lll.z.xx)) + (gamma_uu.xy * d_lll.x.yz) + (((gamma_uu.xy * d_lll.y.xz) - (2. * gamma_uu.xy * d_lll.z.xy)) - (2. * gamma_uu.xy * f * d_lll.z.xy)) + (((gamma_uu.xz * d_lll.x.zz) - (gamma_uu.xz * d_lll.z.xz)) - (2. * gamma_uu.xz * f * d_lll.z.xz)) + (((gamma_uu.yy * d_lll.y.yz) - (gamma_uu.yy * d_lll.z.yy)) - (gamma_uu.yy * f * d_lll.z.yy)) + (((gamma_uu.yz * d_lll.y.zz) - (gamma_uu.yz * d_lll.z.yz)) - (2. * gamma_uu.yz * f * d_lll.z.yz)));
	results.ptr[24] = ((a_l.y + ((2. * K_ll.xy) - (2. * Z_l.y)) + ((gamma_uu.xx * d_lll.y.xx) - (2. * gamma_uu.xy * d_lll.x.yy)) + ((2. * gamma_uu.xy * d_lll.y.xy) - (2. * gamma_uu.xz * d_lll.x.yz)) + (((2. * gamma_uu.xz * d_lll.y.xz) - (gamma_uu.yy * d_lll.y.yy)) - (2. * gamma_uu.yz * d_lll.z.yy)) + ((gamma_uu.zz * d_lll.y.zz) - (2. * gamma_uu.zz * d_lll.z.yz))) / 2.);
	results.ptr[25] = ((a_l.z + ((2. * K_ll.xz) - (2. * Z_l.z)) + ((gamma_uu.xx * d_lll.z.xx) - (2. * gamma_uu.xy * d_lll.x.yz)) + ((2. * gamma_uu.xy * d_lll.z.xy) - (2. * gamma_uu.xz * d_lll.x.zz)) + ((2. * gamma_uu.xz * d_lll.z.xz) - (2. * gamma_uu.yy * d_lll.y.yz)) + (((gamma_uu.yy * d_lll.z.yy) - (2. * gamma_uu.yz * d_lll.y.zz)) - (gamma_uu.zz * d_lll.z.zz))) / 2.);
	results.ptr[26] = (K_ll.yy + (gamma_uu.xx * d_lll.x.yy) + (gamma_uu.xy * d_lll.y.yy) + (gamma_uu.xz * d_lll.z.yy));
	results.ptr[27] = (K_ll.yz + (gamma_uu.xx * d_lll.x.yz) + (gamma_uu.xy * d_lll.y.yz) + (gamma_uu.xz * d_lll.z.yz));
	results.ptr[28] = (K_ll.zz + (gamma_uu.xx * d_lll.x.zz) + (gamma_uu.xy * d_lll.y.zz) + (gamma_uu.xz * d_lll.z.zz));
	results.ptr[29] = ((Theta - (gamma_uu.xx * Z_l.x)) + ((gamma_uu.xx * gamma_uu.yy * d_lll.x.yy) - (gamma_uu.xx * gamma_uu.yy * d_lll.y.xy)) + (((2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz) - (gamma_uu.xx * gamma_uu.yz * d_lll.y.xz)) - (gamma_uu.xx * gamma_uu.yz * d_lll.z.xy)) + (((gamma_uu.xx * gamma_uu.zz * d_lll.x.zz) - (gamma_uu.xx * gamma_uu.zz * d_lll.z.xz)) - (gamma_uu.xy * gamma_uu.xy * d_lll.x.yy)) + (((gamma_uu.xy * gamma_uu.xy * d_lll.y.xy) - (gamma_uu.xy * Z_l.y)) - (2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz)) + (gamma_uu.xy * gamma_uu.xz * d_lll.y.xz) + (gamma_uu.xy * gamma_uu.xz * d_lll.z.xy) + ((gamma_uu.xy * gamma_uu.yz * d_lll.y.yz) - (gamma_uu.xy * gamma_uu.yz * d_lll.z.yy)) + (((gamma_uu.xy * gamma_uu.zz * d_lll.y.zz) - (gamma_uu.xy * gamma_uu.zz * d_lll.z.yz)) - (gamma_uu.xz * gamma_uu.xz * d_lll.x.zz)) + (((gamma_uu.xz * gamma_uu.xz * d_lll.z.xz) - (gamma_uu.xz * Z_l.z)) - (gamma_uu.xz * gamma_uu.yy * d_lll.y.yz)) + ((gamma_uu.xz * gamma_uu.yy * d_lll.z.yy) - (gamma_uu.xz * gamma_uu.yz * d_lll.y.zz)) + (gamma_uu.xz * gamma_uu.yz * d_lll.z.yz));
	results.ptr[30] = (-(((lambda_2 * gamma_uu.xx * Z_l.x) - (lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (lambda_2 * gamma_uu.xy * Z_l.y) + ((((2. * lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (lambda_2 * gamma_uu.xz * Z_l.z) + ((lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((((((((((((lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)) - (sqrt_f * lambda_1 * Theta)) - (sqrt_f * gamma_uu.xx * K_ll.xx)) - (2. * sqrt_f * gamma_uu.xy * K_ll.xy)) - (2. * sqrt_f * gamma_uu.xz * K_ll.xz)) - (sqrt_f * gamma_uu.yy * K_ll.yy)) - (2. * sqrt_f * gamma_uu.yz * K_ll.yz)) - (sqrt_f * gamma_uu.zz * K_ll.zz)) - (gamma_uu.xx * a_l.x)) - (gamma_uu.xy * a_l.y)) - (gamma_uu.xz * a_l.z))));


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

	real sqrt_f = eig.sqrt_f;
	real f = sqrt_f * sqrt_f;
	real fSq = f * f;
	real lambda_1 = (2. - solver->m) / (f - 1);
	real lambda_2 = (2. * f - solver->m) / (f - 1);
	
	resultU.ptr[0] = ((((input.ptr[0] - input.ptr[30]) - (lambda_2 * input.ptr[1])) + (lambda_2 * input.ptr[29]) + (2. * gamma_uu.xy * input.ptr[7]) + (2. * gamma_uu.xz * input.ptr[8])) / (-(2. * gamma_uu.xx)));
	resultU.ptr[1] = input.ptr[7];
	resultU.ptr[2] = input.ptr[8];
	resultU.ptr[3] = ((((2. * input.ptr[0]) - (2. * input.ptr[1])) + (4. * input.ptr[21] * gamma_uu.xx) + (((2. * input.ptr[29]) - (2. * input.ptr[30])) - (2. * lambda_2 * input.ptr[1])) + (((2. * lambda_2 * input.ptr[29]) - (4. * gamma_uu.xy * input.ptr[2] * f)) - (12. * gamma_uu.xy * input.ptr[7] * f)) + (8. * gamma_uu.xy * input.ptr[9] * gamma_uu.xx * f) + (8. * gamma_uu.xy * gamma_uu.xy * input.ptr[10] * f) + (4. * gamma_uu.xy * input.ptr[22]) + (4. * gamma_uu.xy * input.ptr[24] * f) + (8. * gamma_uu.xy * fSq * input.ptr[9] * gamma_uu.xx) + (16. * gamma_uu.xy * gamma_uu.xy * fSq * input.ptr[10]) + (8. * gamma_uu.xy * f * input.ptr[22]) + (8. * gamma_uu.xy * gamma_uu.xz * input.ptr[11] * f) + (8. * gamma_uu.xy * gamma_uu.xz * input.ptr[16] * f) + (16. * gamma_uu.xy * gamma_uu.xz * fSq * input.ptr[11]) + (16. * gamma_uu.xy * gamma_uu.xz * fSq * input.ptr[16]) + (8. * gamma_uu.xy * gamma_uu.yz * input.ptr[13] * f) + (((16. * gamma_uu.xy * gamma_uu.yz * fSq * input.ptr[13]) - (4. * gamma_uu.xz * input.ptr[3] * f)) - (12. * gamma_uu.xz * input.ptr[8] * f)) + (8. * gamma_uu.xz * input.ptr[15] * gamma_uu.xx * f) + (8. * gamma_uu.xz * gamma_uu.xz * input.ptr[17] * f) + (4. * gamma_uu.xz * input.ptr[23]) + (4. * gamma_uu.xz * input.ptr[25] * f) + (8. * gamma_uu.xz * fSq * input.ptr[15] * gamma_uu.xx) + (16. * gamma_uu.xz * gamma_uu.xz * fSq * input.ptr[17]) + (8. * gamma_uu.xz * f * input.ptr[23]) + (8. * gamma_uu.xz * gamma_uu.yz * input.ptr[19] * f) + ((16. * gamma_uu.xz * gamma_uu.yz * fSq * input.ptr[19]) - (2. * gamma_uu.yy * input.ptr[4] * f)) + (2. * gamma_uu.yy * input.ptr[26] * f) + (4. * gamma_uu.yy * gamma_uu.xy * input.ptr[12] * f) + (8. * gamma_uu.yy * gamma_uu.xy * fSq * input.ptr[12]) + (4. * gamma_uu.yy * gamma_uu.xz * input.ptr[18] * f) + ((8. * gamma_uu.yy * gamma_uu.xz * fSq * input.ptr[18]) - (4. * gamma_uu.yz * input.ptr[5] * f)) + ((4. * gamma_uu.yz * input.ptr[27] * f) - (2. * gamma_uu.zz * input.ptr[6] * f)) + (2. * gamma_uu.zz * input.ptr[28] * f) + (4. * gamma_uu.zz * gamma_uu.xy * input.ptr[14] * f) + (8. * gamma_uu.zz * gamma_uu.xy * fSq * input.ptr[14]) + (4. * gamma_uu.zz * gamma_uu.xz * input.ptr[20] * f) + (8. * gamma_uu.zz * gamma_uu.xz * fSq * input.ptr[20])) / (-(4. * gammaUUxxSq * f)));
	resultU.ptr[4] = ((-(input.ptr[2] + (((((((3. * input.ptr[7]) - (input.ptr[9] * gamma_uu.xx)) - (2. * input.ptr[22])) - input.ptr[24]) - (2. * f * input.ptr[9] * gamma_uu.xx)) - (4. * gamma_uu.xy * f * input.ptr[10])) - (2. * gamma_uu.xz * input.ptr[11])) + ((((((((2. * gamma_uu.xz * input.ptr[16]) - (4. * gamma_uu.xz * f * input.ptr[11])) - (gamma_uu.yy * input.ptr[12])) - (2. * gamma_uu.yy * f * input.ptr[12])) - (2. * gamma_uu.yz * input.ptr[13])) - (4. * gamma_uu.yz * f * input.ptr[13])) - (gamma_uu.zz * input.ptr[14])) - (2. * gamma_uu.zz * f * input.ptr[14])))) / (2. * gamma_uu.xx));
	resultU.ptr[5] = ((-(input.ptr[3] + (((((3. * input.ptr[8]) - (input.ptr[15] * gamma_uu.xx)) - (2. * input.ptr[23])) - input.ptr[25]) - (2. * f * input.ptr[15] * gamma_uu.xx)) + ((((((((((2. * gamma_uu.xy * input.ptr[11]) - (2. * gamma_uu.xy * input.ptr[16])) - (4. * gamma_uu.xy * f * input.ptr[16])) - (4. * gamma_uu.xz * f * input.ptr[17])) - (gamma_uu.yy * input.ptr[18])) - (2. * gamma_uu.yy * f * input.ptr[18])) - (2. * gamma_uu.yz * input.ptr[19])) - (4. * gamma_uu.yz * f * input.ptr[19])) - (gamma_uu.zz * input.ptr[20])) - (2. * gamma_uu.zz * f * input.ptr[20])))) / (2. * gamma_uu.xx));
	resultU.ptr[6] = ((-((input.ptr[4] - input.ptr[26]) + (2. * gamma_uu.xy * input.ptr[12]) + (2. * gamma_uu.xz * input.ptr[18]))) / (2. * gamma_uu.xx));
	resultU.ptr[7] = ((-((input.ptr[5] - input.ptr[27]) + (2. * gamma_uu.xy * input.ptr[13]) + (2. * gamma_uu.xz * input.ptr[19]))) / (2. * gamma_uu.xx));
	resultU.ptr[8] = ((-((input.ptr[6] - input.ptr[28]) + (2. * gamma_uu.xy * input.ptr[14]) + (2. * gamma_uu.xz * input.ptr[20]))) / (2. * gamma_uu.xx));
	resultU.ptr[9] = input.ptr[9];
	resultU.ptr[10] = input.ptr[10];
	resultU.ptr[11] = input.ptr[11];
	resultU.ptr[12] = input.ptr[12];
	resultU.ptr[13] = input.ptr[13];
	resultU.ptr[14] = input.ptr[14];
	resultU.ptr[15] = input.ptr[15];
	resultU.ptr[16] = input.ptr[16];
	resultU.ptr[17] = input.ptr[17];
	resultU.ptr[18] = input.ptr[18];
	resultU.ptr[19] = input.ptr[19];
	resultU.ptr[20] = input.ptr[20];
	resultU.ptr[21] = ((input.ptr[0] + ((((((((((((input.ptr[30] - (lambda_1 * input.ptr[1] * sqrt_f)) - (lambda_1 * input.ptr[29] * sqrt_f)) - (2. * gamma_uu.xy * input.ptr[2] * sqrt_f)) - (2. * gamma_uu.xy * input.ptr[24] * sqrt_f)) - (2. * gamma_uu.xz * input.ptr[3] * sqrt_f)) - (2. * gamma_uu.xz * input.ptr[25] * sqrt_f)) - (gamma_uu.yy * input.ptr[4] * sqrt_f)) - (gamma_uu.yy * input.ptr[26] * sqrt_f)) - (2. * gamma_uu.yz * input.ptr[5] * sqrt_f)) - (2. * gamma_uu.yz * input.ptr[27] * sqrt_f)) - (gamma_uu.zz * input.ptr[6] * sqrt_f)) - (gamma_uu.zz * input.ptr[28] * sqrt_f))) / (2. * sqrt_f * gamma_uu.xx));
	resultU.ptr[22] = ((input.ptr[2] + input.ptr[24]) / 2.);
	resultU.ptr[23] = ((input.ptr[3] + input.ptr[25]) / 2.);
	resultU.ptr[24] = ((input.ptr[4] + input.ptr[26]) / 2.);
	resultU.ptr[25] = ((input.ptr[5] + input.ptr[27]) / 2.);
	resultU.ptr[26] = ((input.ptr[6] + input.ptr[28]) / 2.);
	resultU.ptr[27] = ((input.ptr[1] + input.ptr[29]) / 2.);
	resultU.ptr[28] = ((((((input.ptr[1] - input.ptr[29]) - (gamma_uu.xy * input.ptr[2])) - (gamma_uu.xy * input.ptr[7])) - (gamma_uu.xy * input.ptr[9] * gamma_uu.xx)) + (((((gamma_uu.xy * input.ptr[24]) - (2. * gamma_uu.xy * gamma_uu.yz * input.ptr[13])) - (gamma_uu.xz * input.ptr[3])) - (gamma_uu.xz * input.ptr[8])) - (gamma_uu.xz * input.ptr[15] * gamma_uu.xx)) + (((gamma_uu.xz * input.ptr[25]) - (gamma_uu.yy * input.ptr[4])) - (2. * gamma_uu.yy * input.ptr[10] * gamma_uu.xx)) + ((((((gamma_uu.yy * input.ptr[26]) - (gamma_uu.yy * gamma_uu.xy * input.ptr[12])) - (gamma_uu.yy * gamma_uu.xz * input.ptr[18])) - (2. * gamma_uu.yz * input.ptr[5])) - (2. * gamma_uu.yz * input.ptr[11] * gamma_uu.xx)) - (2. * gamma_uu.yz * input.ptr[16] * gamma_uu.xx)) + ((((2. * gamma_uu.yz * input.ptr[27]) - (2. * gamma_uu.yz * gamma_uu.xz * input.ptr[19])) - (gamma_uu.zz * input.ptr[6])) - (2. * gamma_uu.zz * input.ptr[17] * gamma_uu.xx)) + (((gamma_uu.zz * input.ptr[28]) - (gamma_uu.zz * gamma_uu.xy * input.ptr[14])) - (gamma_uu.zz * gamma_uu.xz * input.ptr[20]))) / (2. * gamma_uu.xx));
	resultU.ptr[29] = (((input.ptr[2] * gamma_uu.xx) + ((input.ptr[7] * gamma_uu.xx) - (input.ptr[24] * gamma_uu.xx)) + (((gammaUUxxSq * input.ptr[9]) - (gamma_uu.xx * gamma_uu.yy * input.ptr[12])) - (2. * gamma_uu.xx * gamma_uu.yz * input.ptr[18])) + (gamma_uu.xy * input.ptr[4]) + (2. * gamma_uu.xy * input.ptr[10] * gamma_uu.xx) + ((2. * gamma_uu.xy * gamma_uu.xy * input.ptr[12]) - (gamma_uu.xy * input.ptr[26])) + (2. * gamma_uu.xy * gamma_uu.xz * input.ptr[13]) + (gamma_uu.xz * input.ptr[5]) + (2. * gamma_uu.xz * input.ptr[11] * gamma_uu.xx) + ((2. * gamma_uu.xz * gamma_uu.xz * input.ptr[19]) - (gamma_uu.xz * input.ptr[27])) + (2. * gamma_uu.xz * gamma_uu.xy * input.ptr[18]) + ((gamma_uu.zz * input.ptr[14] * gamma_uu.xx) - (2. * gamma_uu.zz * gamma_uu.xx * input.ptr[19]))) / (2. * gamma_uu.xx));
	resultU.ptr[30] = (((input.ptr[3] * gamma_uu.xx) + ((input.ptr[8] * gamma_uu.xx) - (input.ptr[25] * gamma_uu.xx)) + ((((gammaUUxxSq * input.ptr[15]) - (2. * gamma_uu.xx * gamma_uu.yy * input.ptr[13])) - (2. * gamma_uu.xx * gamma_uu.yz * input.ptr[14])) - (gamma_uu.xx * gamma_uu.zz * input.ptr[20])) + (gamma_uu.xy * input.ptr[5]) + (2. * gamma_uu.xy * gamma_uu.xy * input.ptr[13]) + ((2. * gamma_uu.xy * input.ptr[16] * gamma_uu.xx) - (gamma_uu.xy * input.ptr[27])) + (2. * gamma_uu.xy * gamma_uu.xz * input.ptr[19]) + (gamma_uu.xz * input.ptr[6]) + (2. * gamma_uu.xz * input.ptr[17] * gamma_uu.xx) + ((2. * gamma_uu.xz * gamma_uu.xz * input.ptr[20]) - (gamma_uu.xz * input.ptr[28])) + (2. * gamma_uu.xz * gamma_uu.xy * input.ptr[14]) + (gamma_uu.yy * input.ptr[18] * gamma_uu.xx)) / (2. * gamma_uu.xx));

	resultU.a_l = real3_swap<?=side?>(resultU.a_l);							//0-2
	resultU.d_lll = (_3sym3){
		.x = sym3_swap<?=side?>(resultU.d_lll.v<?=side?>),					//3-8
		.y = sym3_swap<?=side?>(resultU.d_lll.v<?=side==1 and 0 or 1?>),		//9-14
		.z = sym3_swap<?=side?>(resultU.d_lll.v<?=side==2 and 0 or 2?>),		//15-20
	};
	resultU.K_ll = sym3_swap<?=side?>(resultU.K_ll);							//21-26
	resultU.Theta = resultU.Theta;												//27
	resultU.Z_l = real3_swap<?=side?>(resultU.Z_l);							//28-30

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

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf)
{
	SETBOUNDS_NOGHOST();
	global cons_t* U = UBuf + index;
	global cons_t* deriv = derivBuf + index;

	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real f = calc_f(U->alpha);	//could be based on alpha...

	real3 S_l = real3_zero;
	sym3 S_ll = sym3_zero;
	real S = 0.;
	real rho = 0.;

	// source terms
	
	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real trK = real3x3_trace(K_ul);								//K^k_k
	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		//KSq_ij = K_ik K^k_j

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};

	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//e_l = d^j_ji
	real3 e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end	?>,
<? end
?>	};

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


	//d_l = d_ij^j
	real3 d_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3x3_trace(d_llu[<?=i-1?>]),
<? end
?>	};
	
	real3 d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 Z_u = sym3_real3_mul(gamma_uu, U->Z_l);

	real3 V_l = real3_sub(d_l, e_l);

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


	//d_luu
	_3sym3 d_luu = (_3sym3){
<? for i,xi in ipairs(xNames) do		
?>		.<?=xi?> = sym3_real3x3_to_sym3_mul(gamma_uu, d_llu[<?=i-1?>]),
<? end
?>	};

	//alpha_,t = shift terms - alpha^2 f gamma^ij K_ij
	deriv->alpha += -U->alpha * U->alpha * f * trK;
	
	//gamma_ij,t = shift terms - 2 alpha K_ij
	sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));

	//2005 Bona et al A.1
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->K_ll.<?=xij?> += U->alpha * (
		.5 * (U->a_l.<?=xi?> * d_l.<?=xj?> + U->a_l.<?=xj?> * d_l.<?=xi?>)
		- .25 * (U->a_l.<?=xj?> * (2. * e_l.<?=xi?> - d_l.<?=xi?>) + U->a_l.<?=xi?> * (2. * e_l.<?=xj?> - d_l.<?=xj?>))
		- U->a_l.<?=xi?> * U->Z_l.<?=xj?> - U->a_l.<?=xj?> * U->Z_l.<?=xi?>
		+ (trK - 2. * U->Theta) * U->K_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- .5 * U->a_l.<?=xk?> * conn_ull.<?=xk?>.<?=xij?>
		+ .5 * U->a_l.<?=xk?> * d_ull.<?=xk?>.<?=xij?>
		- 2. * e_l.<?=xk?> * (d_llu[<?=i-1?>].<?=xj?>.<?=xk?> + d_llu[<?=j-1?>].<?=xi?>.<?=xk?>)
		+ (d_l.<?=xk?> + U->a_l.<?=xk?> - 2. * U->Z_l.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		- 2. * K_ul.<?=xk?>.<?=xi?> * U->K_ll.<?=sym(k,j)?>
<?		for l,xl in ipairs(xNames) do
?>		+ 2. * (d_llu[<?=i-1?>].<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,j)?> + d_llu[<?=j-1?>].<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,i)?>)
		- conn_ull.<?=xk?>.<?=sym(l,j)?> * conn_ull.<?=xl?>.<?=sym(k,i)?>
<?		end
	end
?> - 8. * M_PI * (S_ll.<?=xij?> - .5 * U->gamma_ll.<?=xij?> * (S - rho)));
<? end
?>

	//2005 Bona et al A.2
<? for i,xi in ipairs(xNames) do
?>	deriv->Z_l.<?=xi?> += U->alpha * (
		U->a_l.<?=xi?> * (trK - 2. * U->Theta)
<?	for k,xk in ipairs(xNames) do
?>		- U->a_l.<?=xk?> * K_ul.<?=xk?>.<?=xi?>
		+ K_ul.<?=xk?>.<?=xi?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?		for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * conn_ull.<?=xr?>.<?=sym(k,i)?>
<?		end
	end
?>
	) - 8 * M_PI * U->alpha * S_l.<?=xi?>;
<? end 
?>
	//2005 Bona et al A.3
	deriv->Theta += U->alpha * .5 * ( trK * (trK - 2. * U->Theta)
<? 
for k,xk in ipairs(xNames) do 
?>		+ 2. * U->a_l.<?=xk?> * (d_u.<?=xk?> - e_u.<?=xk?> - 2. * Z_u.<?=xk?>) 
		- d_u.<?=xk?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?	for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * K_ul.<?=xr?>.<?=xk?>
<?		for s,xs in ipairs(xNames) do
?>		+ d_luu.<?=xk?>.<?=sym(r,k)?> * conn_ull.<?=xk?>.<?=sym(r,s)?>
<?		end
	end
end?>
	) - 8. * M_PI * U->alpha * rho;
}
