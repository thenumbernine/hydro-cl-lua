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

	real sqrt_gammaUUxx = sqrt(gammaUUxx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gammaUUxx;
	real gammaUUxxSq = gammaUUxx * gammaUUxx;

	results[0] = (((\lambda_2 * gammaUxx * input[28]) - (\lambda_2 * gammaUxx * gammaUyy * input[6])) + ((\lambda_2 * gammaUxx * gammaUyy * input[10]) - (2. * \lambda_2 * gammaUxx * gammaUyz * input[7])) + (\lambda_2 * gammaUxx * gammaUyz * input[11]) + ((\lambda_2 * gammaUxx * gammaUyz * input[16]) - (\lambda_2 * gammaUxx * gammaUzz * input[8])) + (\lambda_2 * gammaUxx * gammaUzz * input[17]) + ((\lambda_2 * (gammaUxy ^ 2.) * input[6]) - (\lambda_2 * (gammaUxy ^ 2.) * input[10])) + (\lambda_2 * gammaUxy * input[29]) + ((((2. * \lambda_2 * gammaUxy * gammaUxz * input[7]) - (\lambda_2 * gammaUxy * gammaUxz * input[11])) - (\lambda_2 * gammaUxy * gammaUxz * input[16])) - (\lambda_2 * gammaUxy * gammaUyz * input[13])) + ((\lambda_2 * gammaUxy * gammaUyz * input[18]) - (\lambda_2 * gammaUxy * gammaUzz * input[14])) + (\lambda_2 * gammaUxy * gammaUzz * input[19]) + ((\lambda_2 * (gammaUxz ^ 2.) * input[8]) - (\lambda_2 * (gammaUxz ^ 2.) * input[17])) + (\lambda_2 * gammaUxz * input[30]) + ((\lambda_2 * gammaUxz * gammaUyy * input[13]) - (\lambda_2 * gammaUxz * gammaUyy * input[18])) + ((\lambda_2 * gammaUxz * gammaUyz * input[14]) - (\lambda_2 * gammaUxz * gammaUyz * input[19])) + (sqrt_f * \lambda_1 * input[27]) + (sqrt_f * gammaUxx * input[21]) + (2. * sqrt_f * gammaUxy * input[22]) + (2. * sqrt_f * gammaUxz * input[23]) + (sqrt_f * gammaUyy * input[24]) + (2. * sqrt_f * gammaUyz * input[25]) + ((((sqrt_f * gammaUzz * input[26]) - (gammaUxx * input[0])) - (gammaUxy * input[1])) - (gammaUxz * input[2])));
	results[1] = (input[27] + ((gammaUxx * input[28]) - (gammaUxx * gammaUyy * input[6])) + ((gammaUxx * gammaUyy * input[10]) - (2. * gammaUxx * gammaUyz * input[7])) + (gammaUxx * gammaUyz * input[11]) + ((gammaUxx * gammaUyz * input[16]) - (gammaUxx * gammaUzz * input[8])) + (gammaUxx * gammaUzz * input[17]) + (((gammaUxy ^ 2.) * input[6]) - ((gammaUxy ^ 2.) * input[10])) + (gammaUxy * input[29]) + ((((2. * gammaUxy * gammaUxz * input[7]) - (gammaUxy * gammaUxz * input[11])) - (gammaUxy * gammaUxz * input[16])) - (gammaUxy * gammaUyz * input[13])) + ((gammaUxy * gammaUyz * input[18]) - (gammaUxy * gammaUzz * input[14])) + (gammaUxy * gammaUzz * input[19]) + (((gammaUxz ^ 2.) * input[8]) - ((gammaUxz ^ 2.) * input[17])) + (gammaUxz * input[30]) + ((gammaUxz * gammaUyy * input[13]) - (gammaUxz * gammaUyy * input[18])) + ((gammaUxz * gammaUyz * input[14]) - (gammaUxz * gammaUyz * input[19])));
	results[2] = ((-(((input[1] - (2. * input[22])) - (2. * input[29])) + ((gammaUxx * input[9]) - (2. * gammaUxy * input[6])) + ((2. * gammaUxy * input[10]) - (2. * gammaUxz * input[7])) + (((2. * gammaUxz * input[11]) - (gammaUyy * input[12])) - (2. * gammaUyz * input[18])) + ((gammaUzz * input[14]) - (2. * gammaUzz * input[19])))) / 2.);
	results[3] = ((-(((input[2] - (2. * input[23])) - (2. * input[30])) + ((gammaUxx * input[15]) - (2. * gammaUxy * input[7])) + ((2. * gammaUxy * input[16]) - (2. * gammaUxz * input[8])) + ((2. * gammaUxz * input[17]) - (2. * gammaUyy * input[13])) + (((gammaUyy * input[18]) - (2. * gammaUyz * input[14])) - (gammaUzz * input[20])))) / 2.);
	results[4] = (((input[24] - (gammaUxx * input[6])) - (gammaUxy * input[12])) - (gammaUxz * input[18]));
	results[5] = (((input[25] - (gammaUxx * input[7])) - (gammaUxy * input[13])) - (gammaUxz * input[19]));
	results[6] = (((input[26] - (gammaUxx * input[8])) - (gammaUxy * input[14])) - (gammaUxz * input[20]));
	results[7] = input[1];
	results[8] = input[2];
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
	results[21] = (input[0] + ((input[28] - (f * gammaUxx * input[3])) - (gammaUxy * input[4])) + (((gammaUxy * input[9]) - (2. * gammaUxy * f * input[4])) - (gammaUxz * input[5])) + (((gammaUxz * input[15]) - (2. * gammaUxz * f * input[5])) - (gammaUyy * input[6])) + (((gammaUyy * input[10]) - (gammaUyy * f * input[6])) - (2. * gammaUyz * input[7])) + (gammaUyz * input[11]) + (((gammaUyz * input[16]) - (2. * gammaUyz * f * input[7])) - (gammaUzz * input[8])) + ((gammaUzz * input[17]) - (gammaUzz * f * input[8])));
	results[22] = (input[1] + (input[29] - (f * gammaUyy * input[12])) + (((gammaUxx * input[4]) - (gammaUxx * input[9])) - (gammaUxx * f * input[9])) + (((gammaUxy * input[6]) - (gammaUxy * input[10])) - (2. * gammaUxy * f * input[10])) + ((gammaUxz * input[7]) - (2. * gammaUxz * input[11])) + (((gammaUxz * input[16]) - (2. * gammaUxz * f * input[11])) - (gammaUyz * input[13])) + (((gammaUyz * input[18]) - (2. * gammaUyz * f * input[13])) - (gammaUzz * input[14])) + ((gammaUzz * input[19]) - (gammaUzz * f * input[14])));
	results[23] = (input[2] + (input[30] - (f * gammaUzz * input[20])) + (((gammaUxx * input[5]) - (gammaUxx * input[15])) - (gammaUxx * f * input[15])) + (gammaUxy * input[7]) + (((gammaUxy * input[11]) - (2. * gammaUxy * input[16])) - (2. * gammaUxy * f * input[16])) + (((gammaUxz * input[8]) - (gammaUxz * input[17])) - (2. * gammaUxz * f * input[17])) + (((gammaUyy * input[13]) - (gammaUyy * input[18])) - (gammaUyy * f * input[18])) + (((gammaUyz * input[14]) - (gammaUyz * input[19])) - (2. * gammaUyz * f * input[19])));
	results[24] = ((input[1] + ((2. * input[22]) - (2. * input[29])) + ((gammaUxx * input[9]) - (2. * gammaUxy * input[6])) + ((2. * gammaUxy * input[10]) - (2. * gammaUxz * input[7])) + (((2. * gammaUxz * input[11]) - (gammaUyy * input[12])) - (2. * gammaUyz * input[18])) + ((gammaUzz * input[14]) - (2. * gammaUzz * input[19]))) / 2.);
	results[25] = ((input[2] + ((2. * input[23]) - (2. * input[30])) + ((gammaUxx * input[15]) - (2. * gammaUxy * input[7])) + ((2. * gammaUxy * input[16]) - (2. * gammaUxz * input[8])) + ((2. * gammaUxz * input[17]) - (2. * gammaUyy * input[13])) + (((gammaUyy * input[18]) - (2. * gammaUyz * input[14])) - (gammaUzz * input[20]))) / 2.);
	results[26] = (input[24] + (gammaUxx * input[6]) + (gammaUxy * input[12]) + (gammaUxz * input[18]));
	results[27] = (input[25] + (gammaUxx * input[7]) + (gammaUxy * input[13]) + (gammaUxz * input[19]));
	results[28] = (input[26] + (gammaUxx * input[8]) + (gammaUxy * input[14]) + (gammaUxz * input[20]));
	results[29] = ((input[27] - (gammaUxx * input[28])) + ((gammaUxx * gammaUyy * input[6]) - (gammaUxx * gammaUyy * input[10])) + (((2. * gammaUxx * gammaUyz * input[7]) - (gammaUxx * gammaUyz * input[11])) - (gammaUxx * gammaUyz * input[16])) + (((gammaUxx * gammaUzz * input[8]) - (gammaUxx * gammaUzz * input[17])) - ((gammaUxy ^ 2.) * input[6])) + ((((gammaUxy ^ 2.) * input[10]) - (gammaUxy * input[29])) - (2. * gammaUxy * gammaUxz * input[7])) + (gammaUxy * gammaUxz * input[11]) + (gammaUxy * gammaUxz * input[16]) + ((gammaUxy * gammaUyz * input[13]) - (gammaUxy * gammaUyz * input[18])) + (((gammaUxy * gammaUzz * input[14]) - (gammaUxy * gammaUzz * input[19])) - ((gammaUxz ^ 2.) * input[8])) + ((((gammaUxz ^ 2.) * input[17]) - (gammaUxz * input[30])) - (gammaUxz * gammaUyy * input[13])) + ((gammaUxz * gammaUyy * input[18]) - (gammaUxz * gammaUyz * input[14])) + (gammaUxz * gammaUyz * input[19]));
	results[30] = (-(((\lambda_2 * gammaUxx * input[28]) - (\lambda_2 * gammaUxx * gammaUyy * input[6])) + ((\lambda_2 * gammaUxx * gammaUyy * input[10]) - (2. * \lambda_2 * gammaUxx * gammaUyz * input[7])) + (\lambda_2 * gammaUxx * gammaUyz * input[11]) + ((\lambda_2 * gammaUxx * gammaUyz * input[16]) - (\lambda_2 * gammaUxx * gammaUzz * input[8])) + (\lambda_2 * gammaUxx * gammaUzz * input[17]) + ((\lambda_2 * (gammaUxy ^ 2.) * input[6]) - (\lambda_2 * (gammaUxy ^ 2.) * input[10])) + (\lambda_2 * gammaUxy * input[29]) + ((((2. * \lambda_2 * gammaUxy * gammaUxz * input[7]) - (\lambda_2 * gammaUxy * gammaUxz * input[11])) - (\lambda_2 * gammaUxy * gammaUxz * input[16])) - (\lambda_2 * gammaUxy * gammaUyz * input[13])) + ((\lambda_2 * gammaUxy * gammaUyz * input[18]) - (\lambda_2 * gammaUxy * gammaUzz * input[14])) + (\lambda_2 * gammaUxy * gammaUzz * input[19]) + ((\lambda_2 * (gammaUxz ^ 2.) * input[8]) - (\lambda_2 * (gammaUxz ^ 2.) * input[17])) + (\lambda_2 * gammaUxz * input[30]) + ((\lambda_2 * gammaUxz * gammaUyy * input[13]) - (\lambda_2 * gammaUxz * gammaUyy * input[18])) + ((((((((((((\lambda_2 * gammaUxz * gammaUyz * input[14]) - (\lambda_2 * gammaUxz * gammaUyz * input[19])) - (sqrt_f * \lambda_1 * input[27])) - (sqrt_f * gammaUxx * input[21])) - (2. * sqrt_f * gammaUxy * input[22])) - (2. * sqrt_f * gammaUxz * input[23])) - (sqrt_f * gammaUyy * input[24])) - (2. * sqrt_f * gammaUyz * input[25])) - (sqrt_f * gammaUzz * input[26])) - (gammaUxx * input[0])) - (gammaUxy * input[1])) - (gammaUxz * input[2]))));


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
		
	real sqrt_gammaUUxx = sqrt(gammaUUxx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gammaUUxx;
	real gammaUUxxSq = gammaUUxx * gammaUUxx;
	
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
			+ conn_ull.<?=xk?>.<?=xij?> * (U->V_l.<?=xk?> - e_l.<?=xk?>)

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

	//srcV_k = -a_k K 
	//			+ 8 pi j_k
	//			+ a^j K_jk 
	//			- 4 K d^j_jk 
	//			+ 2 K^i_k d^j_ji 
	//			+ K^i_k d_ij^j 
	//			- K^i_j d_ki^j
	//			+ 2 K^i_j d_ik^j 
	real3 srcV_l = (real3){
<? for k,xk in ipairs(xNames) do
?>		.<?=xk?> = -tr_K * U->a_l.<?=xk?>
			- 8. * M_PI * S_l.<?=xk?>
<?	for j,xj in ipairs(xNames) do
?>			+ U->a_l.<?=xj?> * K_ul.<?=xj?>.<?=xk?>
			- 4. * tr_K * d_ull.<?=xj?>.<?=sym(j,k)?>
			+ 2. * K_ul.<?=xj?>.<?=xk?> * e_l.<?=xj?>
			+ 2. * K_ul.<?=xj?>.<?=xk?> * d_l.<?=xj?>
<?		for i,xi in ipairs(xNames) do
?>			- K_ul.<?=xi?>.<?=xj?> * d_llu[<?=k-1?>].<?=xi?>.<?=xj?>
			+ 2. * K_ul.<?=xi?>.<?=xj?> * d_llu[<?=i-1?>].<?=xk?>.<?=xj?>
<? 		end
	end ?>,
<? end
?>	};


	//alpha_,t = shift terms - alpha^2 f gamma^ij K_ij
	deriv->alpha += -U->alpha * U->alpha * f * tr_K;
	
	//gamma_ij,t = shift terms - 2 alpha K_ij
	sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));
	
	//K_ij,t = shift terms + alpha srcK_ij
	sym3_add(deriv->K_ll, sym3_real_mul(srcK_ll_over_alpha, U->alpha));

	//V_k,t = shift terms + alpha srcV_k
	real3_add(deriv->V_l, real3_real_mul(srcV_l, U->alpha));

<? if eqn.useShift ~= 'none' then ?>


	<?
	if eqn.useShift == 'HarmonicShiftCondition-FiniteDifference'
	or eqn.useShift == 'MinimalDistortionElliptic'
	or eqn.useShift == 'MinimalDistortionEllipticEvolve'
	then 
	?>

	//make each useShift option into an object
	//and give it a method for producing the partial_beta_ul
	// so the finite difference ones can use this
	// and the state ones can just assign a variable
	//then all this can be moved into the "if eqn.useShift ~= 'none' then" block
	//but then again, any FVS-based method will already have the U_,k beta^k terms incorporated into its eigensystem ...
	//this means the useShift will also determine whether the flux has any shift terms
	// for any '-FiniteDifference' shift, the flux won't have any shift terms.

	//partial_beta_ul[j].i := beta^i_,j
<?=makePartial('beta_u', 'real3')?>	

	//= gamma_ik beta^k_,j
	real3x3 partial_beta_u_ll = (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_real3_mul(gamma_uu, partial_beta_ul[<?=i-1?>]),
<? end
?>	};

	//first add the Lie derivative beta terms to the time derivatives

	//alpha_,t = alpha_,i beta^i + ...
	// = alpha a_i beta^i + ..
	deriv->alpha += real3_dot(U->beta_u, U->a_l) * U->alpha;

	//Alcubierre 2.3.11
	//gamma_ij,t = gamma_ij,k beta^k + gamma_kj beta^k_,i + gamma_ik beta^k_,j
	// = 2 d_kij beta^k + gamma_kj beta^k_,i + gamma_ik beta^k_,j
	//
	//Alcubierre eqn 2.3.11 and B&S eqn 2.128 rewrite it as:
	//gamma_ij,t = -2 alpha K_ij + D_i beta_j + D_j beta_i
	//	= -2 alpha K_ij + D_i beta^k gamma_kj + D_j beta^k gamma_ki
	//	= -2 alpha K_ij + (beta^k_,i + Gamma^k_li beta^l) gamma_kj + (beta^k_,j + Gamma^k_lj beta^l) gamma_ki
	//	= -2 alpha K_ij + beta^k_,i gamma_kj + beta^k_,j gamma_ki + (Gamma_jki + Gamma_ikj) beta^k
	//	= -2 alpha K_ij + beta^k_,i gamma_kj + beta^k_,j gamma_ki + (Gamma_jik + Gamma_ijk) beta^k
	//Gamma_ijk = -d_ijk + d_jik + d_kij
	//Gamma_jik = -d_jik + d_ijk + d_kji
	//Gamma_ijk + Gamma_jik = 2 d_kij
	//gamma_ij,t = -2 alpha K_ij + beta^k_,i gamma_kj + beta^k_,j gamma_ki + 2 d_kij beta^k
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	deriv->gamma_ll.<?=xij?> += partial_beta_u_ll.<?=xi?>.<?=xj?>
		+ partial_beta_u_ll.<?=xj?>.<?=xi?>
<?	for k,xk in ipairs(xNames) do
?> 		+ 2. * U->d_lll.<?=xk?>.<?=xij?> * U->beta_u.<?=xk?>
<?	end
?>	;
<? end ?>

#if 0	//the hyperbolic vars should get shifted via the wavespeeds
		//I should only need to manually shift the source-only vars


	//partial_a_l[j].i = a_i,j
<?=makePartial('a', 'real3')?>

	//a_i,t = a_i,j beta^j + a_j beta^j_,i
<? for i,xi in ipairs(xNames) do
?>	deriv->a.<?=xi?> += 0.<?
	for j,xj in ipairs(xNames) do
?>		+ partial_a_l[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
		+ U->a.<?=xj?> * partial_beta_ul[<?=i-1?>].<?=xj?>
<?	end
?>	;
<? end
?>

	//technically this can be represented as a sym3sym3 since d_ijk,l is symmetric with jk and il
	//partial_d_llll[k][l].ij = d_kij,l
<?=makePartial('d_lll', '_3sym3')?>

	//d_kij,t = d_kij,l beta^l + d_lij beta^l_,k + d_klj beta^l_,i + d_kil beta^l_,j
<? for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>	deriv->d_lll.<?=xk?>.<?=xij?> += 0.
<?		for l,xl in ipairs(xNames) do
?>			+ partial_d_llll[<?=l-1?>].<?=xk?>.<?=xij?>
<?		end
?>	;
<?	end
end ?>

	//partial_K_l[k].ij = K_ij,k
<?=makePartial('K', 'sym3')?>

	//K_ij,t = K_ij,k beta^k + K_kj beta^k_,i + K_ik beta^k_,j
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>
	deriv->K.<?=xij?> += 0.
<?	for k,xk in ipairs(xNames) do
?>		+ partial_K_l[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>
		+ U->K.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
		+ U->K.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
<? 	end 
?>	;
<? end ?>

	//partial_V_l[j].i = V_i,j
<?=makePartial('V_l', 'real3')?>

	//V_i,t = V_i,j beta^j + V_j beta^i_,j
<? for i,xi in ipairs(xNames) do
?>	deriv->V_l.<?=xi?> += 0. <?
	for j,xj in ipairs(xNames) do
?>		+ partial_V_l[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
		+ U->V_l.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?	end
?>	;
<? end
?>

#endif


	<? 
	end	-- using partial_beta_ul to integrate alpha, gamma_ll, 
	if eqn.useShift == '2005 Bona / 2008 Yano' then ?>
	
	
	
	<? elseif eqn.useShift == 'HarmonicShiftCondition-FiniteDifference' then 
	?>

	//next add the source for the particular useShift

	// 2008 Alcubierre 4.3.37
	//beta^i_,t = beta^j beta^i_,j - alpha alpha^,i + alpha^2 conn^i + beta^i / alpha (alpha_,t - beta^j alpha_,j + alpha^2 K)
	//using alpha_,t = beta^i alpha_,i - alpha^2 f K
	//= beta^j beta^i_,j + alpha^2 (conn^i - a^i) + beta^i (beta^j a_j - alpha f K - beta^j a_j + alpha K)
	//= beta^j beta^i_,j + alpha^2 (conn^i - a^i) + alpha K beta^i (1 - f)

	real3 dbeta_beta = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + partial_beta_ul[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?><?
	end ?>,
<? end
?>	};
	
	real3 a_u = sym3_real3_mul(gamma_uu, U->a_l);

	//conn^i = conn^i_jk gamma^jk
	real3 conn_u = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(conn_ull.<?=xi?>, gamma_uu),
<? end
?>	};

	deriv->beta_u = 
	real3_add(
		deriv->beta_u,
		real3_add(
			real3_add(
				dbeta_beta,
				real3_real_mul(
					U->beta_u,
					U->alpha * tr_K * (1. - f)
				)
			),
			real3_real_mul(
				real3_sub(
					conn_u,
					a_u
				),
				U->alpha * U->alpha
			)
		)
	);

	<?
	elseif eqn.useShift == 'MinimalDistortionEllipticEvolve' then
	?>

	//beta^i_,t = epsilon (D^2 beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j - D_j (2 alpha (K^ij - 1/3 K gamma^ij))

	<?
	elseif eqn.useShift == 'LagrangianCoordinates' 
	or eqn.useShift == 'MinimalDistortionElliptic' 
	then
		-- do nothing
	else
		error("I don't have any source terms implemented for this particular useShift")
	end ?>

<? end	-- eqn.useShift ?>

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

	//V_i = d_ik^k - d^k_ki <=> V_i += eta (d_ik^k - d^k_ki - V_i)
	deriv->V_l = real3_add(
		deriv->V_l,
		real3_real_mul(
			real3_sub(real3_sub(d_l, e_l), U->V_l),
			solver->V_convCoeff));

	//Kreiss-Oligar diffusion, for stability's sake?
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

<?
local constrainVGuiVar = eqn.guiVars['constrain V']
local constrainV = constrainVGuiVar.options[constrainVGuiVar.value]  
if constrainV ~= 'none' then 
?>	//gravitational_wave_sim uses this (for 1D), HydroGPU doesn't (for 2D/3D)

	real3 delta;
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d_lll.<?=xi?>, gamma_uu);
		real d2 = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + U->d_lll.<?=xj?>.<?=sym(k,i)?> * gamma_uu.<?=sym(j,k)?><?
		end
	end ?>;
		delta.<?=xi?> = U->V_l.<?=xi?> - (d1 - d2);
	}<? end ?>

<?
	if constrainV == 'replace V' then
?>
	//directly assign to V
	U->V_l = real3_sub(U->V_l, delta);
<?
	elseif constrainV == 'average' then
?>
	//average between V and d

/*
interpolation between V and d
weight = 1 means only adjust V (gives a 1e-20 error in the constraint eqns)
weight = 0 means only adjust d (hmm, is only at 1e-2 ...)

V_i - (d_im^m - d^m_mi) = Q_i ...
V'_i = V_i - alpha Q_i
(d_im^m - d^m_mi)' = (d_im^m - d^m_mi) + (1-alpha) Q_i
so (V_i - (d_im^m - d^m_mi))' = V'_i - (d_im^m - d^m_mi)' = V_i - alpha Q_i - (d_im^m - d^m_mi) - (1-alpha) Q_i = Q_i - Q_i = 0
works for d'_ijk = d_ijk + (1-alpha)/4 (Q_i gamma_jk - Q_k gamma_ij)
so d'_im^m = d'_ijk gamma^jk = d_im^m + (1-alpha)/4 (Q_i gamma_jk - Q_k gamma_ij) gamma^jk = d_im^m + (1-alpha)/2 Q_i
and d'^m_mi = d'_kji gamma^jk = d^m_mi + (1-alpha)/4 (Q_k gamma_ij - Q_i gamma_jk) gamma^jk = d^m_mi - (1-alpha)/2 Q_i
therefore (d_im^m - d^m_mi)' = d_im^m - d^m_mi + (1-alpha) Q_i

...but how come it's not working?
...probably because neighbors are influencing one another
so to solve this, you must solve a giant linear system of all variables
*/
	const real weight = .5;
	const real v_weight = weight;
	const real d_weight = (1. - weight) / 4.;

	U->V_l = real3_sub(U->V_l, real3_real_mul(delta, v_weight));
<? 
		for i,xi in ipairs(xNames) do 
			for jk,xjk in ipairs(symNames) do
				local j,k = from6to3x3(jk)
				local xk = xNames[k]
?>	U->d_lll.<?=xi?>.<?=xjk?> += (
		delta.<?=xi?> * U->gamma_ll.<?=xjk?> 
		- delta.<?=xk?> * U->gamma_ll.<?=sym(i,j)?>
	) * d_weight;
<?	
			end
		end 
	end
?>

//...or linearly project out the [V_i, U->d_ijk] vector
//...or do a single gradient descent step
<?
end	-- constrain V
?>

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
			+ conn_ull.<?=xk?>.<?=xij?> * (U->V_l.<?=xk?> - e_l.<?=xk?>)

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
