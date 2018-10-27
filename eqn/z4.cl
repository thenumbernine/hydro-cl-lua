<?
local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym
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
	real det_gamma = sym3_det(U->gamma_ll);
	real f = calc_f(U->alpha);
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) is only called once
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
		
		real lambdaMin = (real)min((real)0., -lambda);
		real lambdaMax = (real)max((real)0., lambda);
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
	eig.gamma_ll = U.gamma_ll;
	real det_gamma = sym3_det(U.gamma_ll);
	eig.gamma_uu = sym3_inv(U.gamma_ll, det_gamma);
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gamma_uu.xx), sqrt(eig.gamma_uu.yy), sqrt(eig.gamma_uu.zz));
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

	return (range_t){
		.min = -lambdaMax, 
		.max = lambdaMax,
	};
}
<? end ?>

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
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gamma_uu.xx), sqrt(eig.gamma_uu.yy), sqrt(eig.gamma_uu.zz));
	return eig;
}

<? for side=0,solver.dim-1 do ?>
waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 pt
) {
	waves_t results;

	//input
	real3 a_l = real3_swap<?=side?>(inputU.a_l);							//0-2
	_3sym3 d_lll = (_3sym3){
		.x = sym3_swap<?=side?>(inputU.d_lll.v<?=side?>),					//3-8
		.y = sym3_swap<?=side?>(inputU.d_lll.v<?=side==1 and 0 or 1?>),		//9-14
		.z = sym3_swap<?=side?>(inputU.d_lll.v<?=side==2 and 0 or 2?>),		//15-20
	};
	sym3 K_ll = sym3_swap<?=side?>(inputU.K_ll);							//21-26
	real Theta = inputU.Theta;												//27
	real3 Z_l = real3_swap<?=side?>(inputU.Z_l);							//28-30

	//eig
	real sqrt_f = eig.sqrt_f;
	real f = sqrt_f * sqrt_f;	
	real _1_sqrt_f = 1. / sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f;

	real f_minus_1 = f - 1.;
	real f_minus_1_sq = f_minus_1 * f_minus_1;

	sym3 gamma_ll = sym3_swap<?=side?>(eig.gamma_ll);
	sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);
	
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, d_lll);
	
	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(d_lll.<?=xi?>, gamma_uu),
<? end
?>	};

	real3x3 K_ul = sym3_sym3_mul(gamma_uu, K_ll);

	real tr_K = real3x3_trace(K_ul);

	real sqrt_gUxx = eig.sqrt_gammaUjj.s<?=side?>;
	real _1_sqrt_gUxx = 1. / sqrt_gUxx;
	real _1_gammaUxx = _1_sqrt_gUxx * _1_sqrt_gUxx;

	real a_u_x = a_l.x * gamma_uu.xx + a_l.y * gamma_uu.xy + a_l.z * gamma_uu.xz;
	real Z_u_x = Z_l.x * gamma_uu.xx + Z_l.y * gamma_uu.xy + Z_l.z * gamma_uu.xz;

	//d_i = d_ijk gamma^jk
	real3 d_l = _real3(
		sym3_dot(d_lll.x, gamma_uu),
		sym3_dot(d_lll.y, gamma_uu),
		sym3_dot(d_lll.z, gamma_uu));
	
	real d_u_x = d_l.x * gamma_uu.xx + d_l.y * gamma_uu.xy + d_l.z * gamma_uu.xz;

	//e_i = d_jki gamma^jk
	real3 e_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0. <? 
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + d_lll.<?=xj?>.<?=sym(k,i)?> * gamma_uu.<?=sym(j,k)?><?
		end
	end
?>,
<? end
?>	};
	
	real e_u_x = e_l.x * gamma_uu.xx + e_l.y * gamma_uu.xy + e_l.z * gamma_uu.xz;

	real tmp_0_30_a = .5 * _1_gammaUxx * tr_K;
	real tmp_0_30_b = .5 * _1_gammaUxx * (
			a_u_x * _1_sqrt_f * _1_sqrt_gUxx
			+ (
				+ f * (solver->m - 2.) * (Z_u_x + e_u_x - d_u_x) * _1_sqrt_f * _1_sqrt_gUxx
				- (solver->m * f - 2.) * Theta
			) / f_minus_1
		);
	results.ptr[0] = tmp_0_30_a - tmp_0_30_b;
	results.ptr[30] = tmp_0_30_a + tmp_0_30_b;

	real tmp_1_24_a = .25 * (
			+ gamma_uu.xx * (d_llu[1].z.y - d_llu[2].y.y)
			+ gamma_uu.xy * (d_llu[2].y.x - d_llu[1].z.x)
		) * _1_gammaUxx
		- .25 * a_l.z
		+ .5 * Z_l.z;
	real tmp_1_24_b = .5 * K_ul.x.z * _1_sqrt_gUxx;
	results.ptr[1] = tmp_1_24_a + tmp_1_24_b;
	results.ptr[24] = tmp_1_24_a - tmp_1_24_b;

	real tmp_2_25_a = .25 * (
			+ gamma_uu.xx * (d_llu[2].y.z - d_llu[1].z.z)
			+ gamma_uu.xz * (d_llu[1].z.x - d_llu[2].y.x)
		) * _1_gammaUxx
		- .25 * a_l.y
		+ .5 * Z_l.y;
	real tmp_2_25_b = .5 * K_ul.x.y * _1_sqrt_gUxx;
	results.ptr[2] = tmp_2_25_a + tmp_2_25_b;
	results.ptr[25] = tmp_2_25_a - tmp_2_25_b;

	real tmp_3_26_a = 
		-.25 * a_l.x
		+ .5 * Z_l.x
		+ .25 * (
			+ a_u_x
			+ gamma_uu.xy * (d_llu[1].z.z - d_llu[2].y.z)
			+ gamma_uu.xz * (d_llu[2].y.y - d_llu[1].z.y)
		) * _1_gammaUxx;
	real tmp_3_26_b = 2. * _1_sqrt_gUxx * (K_ul.x.x - tr_K + Theta);
	results.ptr[3] = tmp_3_26_a + tmp_3_26_b;
	results.ptr[26] = tmp_3_26_a - tmp_3_26_b;

	real tmp_4_27_a = .5 * Theta;
	real tmp_4_27_b = .5 * (
			+ e_u_x
			- d_u_x
			+ Z_u_x
		) * _1_sqrt_gUxx;
	results.ptr[4] = tmp_4_27_a + tmp_4_27_b;
	results.ptr[27] = tmp_4_27_a - tmp_4_27_b;

	real tmp_5_28_a = .5 * K_ll.zz;
	real tmp_5_28_b = .5 * (
			+ d_llu[2].z.x
			- d_ull.x.zz
		) * _1_sqrt_gUxx;
	results.ptr[5] = tmp_5_28_a + tmp_5_28_b;
	results.ptr[28] = tmp_5_28_a - tmp_5_28_b;

	real tmp_6_29_a = .5 * K_ll.yz;
	real tmp_6_29_b = .25 * (
			+ d_llu[2].y.x 
			+ d_llu[1].z.x 
			- 2. * d_ull.x.yz
		) * _1_sqrt_gUxx;
	results.ptr[6] = tmp_6_29_a + tmp_6_29_b;
	results.ptr[29] = tmp_6_29_a - tmp_6_29_b;

	results.ptr[8] = -.5 * (
		+ d_l.y
		- e_l.y
		
		+ d_ull.x.xy
		+ d_ull.x.xy
		
		- d_llu[0].y.x
		- d_llu[1].x.x
		
		+ a_l.y
		- 2. * Z_l.y
	) * _1_gammaUxx
	+ d_lll.x.xy;

	results.ptr[9] = -.5 * (
		+ d_l.z
		- e_l.z
		
		+ d_ull.x.xz
		+ d_ull.x.xz
		
		- d_llu[2].x.x
		- d_llu[0].z.x
		
		+ a_l.z
		- 2. * Z_l.z
	) * _1_gammaUxx
	+ d_lll.x.xz;
	
	results.ptr[22] = .5 * (
		- gamma_uu.xx * (d_llu[2].y.z - d_llu[1].z.z)
		- gamma_uu.xz * (d_llu[1].z.x - d_llu[2].y.x)
	) * _1_gammaUxx
	+ .5 * a_l.y;

	results.ptr[23] = .5 * (
		- gamma_uu.xx * (d_llu[1].z.y - d_llu[2].y.y)
		- gamma_uu.xy * (d_llu[2].y.x - d_llu[1].z.x)
	) * _1_gammaUxx
	+ .5 * a_l.z;

	results.ptr[7] = (
		(solver->m - 1.) * (e_u_x - d_u_x)
		
		+ gamma_uu.xx * gamma_uu.xx * d_lll.x.xx
		+ gamma_uu.xy * gamma_uu.xy * d_lll.y.xy
		+ gamma_uu.xy * gamma_uu.yz * d_lll.y.yz
		+ gamma_uu.xy * gamma_uu.zz * d_lll.y.zz
		+ gamma_uu.xy * gamma_uu.xz * d_lll.y.xz
		+ gamma_uu.xz * gamma_uu.xy * d_lll.z.xy
		+ gamma_uu.xz * gamma_uu.xz * d_lll.z.xz
		+ gamma_uu.xz * gamma_uu.yy * d_lll.z.yy
		+ gamma_uu.xz * gamma_uu.yz * d_lll.z.yz
		+ gamma_uu.xz * gamma_uu.zz * d_lll.z.zz
		
		- gamma_uu.xy * gamma_uu.xy * d_lll.x.yy
		- gamma_uu.xy * gamma_uu.xy * d_lll.x.yy
		- gamma_uu.xy * gamma_uu.xz * d_lll.x.yz
		- gamma_uu.xy * gamma_uu.xz * d_lll.x.yz
		- gamma_uu.xz * gamma_uu.xz * d_lll.x.zz
		- gamma_uu.yy * gamma_uu.xz * d_lll.y.yz
		- gamma_uu.yz * gamma_uu.xz * d_lll.y.zz
		- gamma_uu.yz * gamma_uu.xy * d_lll.z.yy
		- gamma_uu.zz * gamma_uu.xy * d_lll.z.yz
		- gamma_uu.zz * gamma_uu.xz * d_lll.z.zz
		
		- gamma_uu.xx * a_l.x

		+ a_u_x * f_minus_1 / f
		+ 2 * gamma_uu.xx * Z_l.x
		+ (solver->m - 2.) * Z_u_x
	
	) * _1_gammaUxx * _1_gammaUxx;

	results.ptr[10] = d_lll.y.xx;
	results.ptr[11] = d_lll.y.xy;
	results.ptr[12] = d_lll.y.xz;
	results.ptr[13] = d_lll.y.yy;
	results.ptr[14] = d_lll.y.yz;
	results.ptr[15] = d_lll.y.zz;
	results.ptr[16] = d_lll.z.xx;
	results.ptr[17] = d_lll.z.xy;
	results.ptr[18] = d_lll.z.xz;
	results.ptr[19] = d_lll.z.yy;
	results.ptr[20] = d_lll.z.yz;
	results.ptr[21] = d_lll.z.zz;

	return results;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t input,
	real3 unused
) {
	cons_t resultU;
	for (int j = 0; j < numStates; ++j) {
		resultU.ptr[j] = 0;
	}
	
	sym3 gamma_ll = sym3_swap<?=side?>(eig.gamma_ll);
	sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);

	real sqrt_f = eig.sqrt_f;
	real f = sqrt_f * sqrt_f;
	real fSq = f * f;
	real f_3_2 = f * sqrt_f;

	real f_minus_1 = f - 1.;
	real sqrt_gUxx = sqrt(gamma_uu.xx);
	real _1_sqrt_gUxx = 1. / sqrt_gUxx;
	real _1_gammaUxx = _1_sqrt_gUxx * _1_sqrt_gUxx;

	resultU.a_l.x = -(
		+ (	
			- gamma_uu.xy * gamma_uu.yz * input.ptr[14] 
			- gamma_uu.xy * gamma_uu.zz * input.ptr[15] 
			+ gamma_uu.xy * gamma_uu.yz * input.ptr[19] 
			+ gamma_uu.xy * gamma_uu.zz * input.ptr[20] 
			+ gamma_uu.xy * .2 * input.ptr[22] 
			
			+ gamma_uu.xz * gamma_uu.yy * input.ptr[14] 
			+ gamma_uu.xz * gamma_uu.yz * input.ptr[15] 
			- gamma_uu.xz * gamma_uu.yy * input.ptr[19] 
			- gamma_uu.xz * gamma_uu.yz * input.ptr[20] 
			+ gamma_uu.xz * .2 * input.ptr[23] 
		
		) * _1_sqrt_gUxx

		+ sqrt_f * gamma_uu.xx * (input.ptr[0] - input.ptr[30])
		+ (solver->m - 2.) * (input.ptr[4] - input.ptr[27]) * f / f_minus_1
	
	) * _1_sqrt_gUxx;
	
	resultU.a_l.y = (
		+ gamma_uu.xz * gamma_uu.xy * input.ptr[14] 
		+ gamma_uu.xz * gamma_uu.xz * input.ptr[15] 
		- gamma_uu.xz * gamma_uu.xy * input.ptr[19] 
		- gamma_uu.xz * gamma_uu.xz * input.ptr[20] 
	) * _1_gammaUxx
	- gamma_uu.yz * input.ptr[14] 
	- gamma_uu.zz * input.ptr[15] 
	+ gamma_uu.yz * input.ptr[19] 
	+ gamma_uu.zz * input.ptr[20] 
	+ 2. * input.ptr[22]
	;
	
	resultU.a_l.z = -(
		(
			+ gamma_uu.xy * gamma_uu.xy * input.ptr[14] 
			+ gamma_uu.xy * gamma_uu.xz * input.ptr[15] 
			- gamma_uu.xy * gamma_uu.xy * input.ptr[19] 
			- gamma_uu.xy * gamma_uu.xz * input.ptr[20] 
		) * _1_gammaUxx
			
		- gamma_uu.yy * input.ptr[14] 
		- gamma_uu.yz * input.ptr[15] 
		+ gamma_uu.yy * input.ptr[19] 
		+ gamma_uu.yz * input.ptr[20] 
		- 2. * input.ptr[23]
	);
	
	resultU.d_lll.x.xx = (
		- (gamma_uu.xx * gamma_uu.yy - gamma_uu.xy * gamma_uu.xy) * (gamma_uu.xx * gamma_uu.yy - gamma_uu.xy * gamma_uu.xy) * input.ptr[0] * sqrt_gUxx * f_minus_1
		+ (gamma_uu.xx * gamma_uu.yy - gamma_uu.xy * gamma_uu.xy) * gamma_uu.xz * gamma_uu.yy * input.ptr[1] * sqrt_f * f_minus_1
		+ (gamma_uu.xx * gamma_uu.yy - gamma_uu.xy * gamma_uu.xy) * gamma_uu.xy * gamma_uu.yy * input.ptr[2] * sqrt_f * f_minus_1
	
		+ (
			+ 3. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx
			- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		) * input.ptr[3] * sqrt_f * f_minus_1

		+ (
			+ f * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy 
			- f * 3. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.yy 
			+ f * 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy
			+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy 
		) * (input.ptr[4] - input.ptr[27]) * sqrt_f * _1_sqrt_gUxx

		+ (gamma_uu.xx * gamma_uu.yy - gamma_uu.xy * gamma_uu.xy)
			* (gamma_uu.xx * gamma_uu.yy - gamma_uu.xy * gamma_uu.xy)
			* (input.ptr[4] - input.ptr[27])
			* solver->m * f_3_2 * _1_sqrt_gUxx
	

		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * sqrt_gUxx * gamma_uu.yy * sqrt_f 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * 1. / sqrt_gUxx * sqrt_f 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * f_3_2 * 1. / sqrt_gUxx 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * f_3_2 * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[5] * 1. / sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * sqrt_f 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[5] * sqrt_gUxx * sqrt_f 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[5] * f_3_2 * sqrt_gUxx 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[5] * f_3_2 * 1. / sqrt_gUxx * gamma_uu.xy * gamma_uu.xy 

		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[6] * 1. / sqrt_gUxx * sqrt_f 
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[6] * sqrt_gUxx * sqrt_f 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[6] * f_3_2 * 1. / sqrt_gUxx 
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[6] * f_3_2 * sqrt_gUxx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * sqrt_gUxx * gamma_uu.yy * sqrt_f 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * 1. / sqrt_gUxx * sqrt_f 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * f_3_2 * 1. / sqrt_gUxx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * f_3_2 * sqrt_gUxx * gamma_uu.yy 

		- input.ptr[7] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * sqrt_f 
		+ input.ptr[7] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * f_3_2 * gamma_uu.yy * gamma_uu.yy 
		+ 2. * input.ptr[7] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * sqrt_f * gamma_uu.xx * gamma_uu.xx 
		- input.ptr[7] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_f * gamma_uu.xx 
		+ input.ptr[7] * f_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		- 2. * input.ptr[7] * f_3_2 * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx
		
		+ gamma_uu.xz * gamma_uu.yy * input.ptr[24] * gamma_uu.xy * gamma_uu.xy * sqrt_f 
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[24] * f_3_2 * gamma_uu.xx 
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[24] * sqrt_f * gamma_uu.xx 
		- gamma_uu.xz * gamma_uu.yy * input.ptr[24] * f_3_2 * gamma_uu.xy * gamma_uu.xy 
	
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * input.ptr[25] * sqrt_f 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * input.ptr[25] * f_3_2 
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * input.ptr[25] * f_3_2 * gamma_uu.xx 
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * input.ptr[25] * sqrt_f * gamma_uu.xx 

		- 3. * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * gamma_uu.yy * sqrt_f * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * sqrt_f 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * f_3_2 
		+ 3. * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * f_3_2 * gamma_uu.yy * gamma_uu.xx 
		+ gamma_uu.yy * gamma_uu.yy * input.ptr[26] * sqrt_f * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.yy * gamma_uu.yy * input.ptr[26] * f_3_2 * gamma_uu.xx * gamma_uu.xx 
	
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * sqrt_gUxx * gamma_uu.yy * sqrt_f 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * 1. / sqrt_gUxx * sqrt_f 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * f_3_2 * 1. / sqrt_gUxx 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * f_3_2 * sqrt_gUxx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[28] * 1. / sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * sqrt_f 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[28] * sqrt_gUxx * sqrt_f 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[28] * f_3_2 * sqrt_gUxx 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[28] * f_3_2 * 1. / sqrt_gUxx * gamma_uu.xy * gamma_uu.xy 

		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[29] * sqrt_gUxx * sqrt_f 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[29] * 1. / sqrt_gUxx * sqrt_f 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[29] * f_3_2 * 1. / sqrt_gUxx 
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[29] * f_3_2 * sqrt_gUxx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * sqrt_gUxx * gamma_uu.yy * sqrt_f 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * 1. / sqrt_gUxx * sqrt_f 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * f_3_2 * 1. / sqrt_gUxx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * f_3_2 * sqrt_gUxx * gamma_uu.yy 
			
		- input.ptr[30] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		+ 2. * input.ptr[30] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
		- input.ptr[30] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy 
		+ input.ptr[30] * sqrt_gUxx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * input.ptr[30] * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
		+ input.ptr[30] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.yy * gamma_uu.yy 
	) / (
		-gamma_uu.xx * sqrt_f * (
			+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
			+ 2. * gamma_uu.xx * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
			+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy 
			- gamma_uu.xx * gamma_uu.xx * f * gamma_uu.yy * gamma_uu.yy
		)
	);
	
	resultU.d_lll.x.xy = (
		+ 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] 
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * gamma_uu.yy 
		- 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] 
		+ 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * gamma_uu.yy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * input.ptr[28] 
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * input.ptr[28] * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * input.ptr[5] 
		+ gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * input.ptr[5] * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[1] 
		+ gamma_uu.xy * gamma_uu.xz * input.ptr[1] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[24] 
		+ gamma_uu.xy * gamma_uu.xz * input.ptr[24] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[28] * 1. / sqrt_gUxx 
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[28] * sqrt_gUxx * gamma_uu.yy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[5] * 1. / sqrt_gUxx 
		- gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[5] * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * gamma_uu.xx 
		+ gamma_uu.xy * input.ptr[26] * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[27] * sqrt_gUxx 
		+ gamma_uu.xy * input.ptr[27] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * gamma_uu.xx 
		+ gamma_uu.xy * input.ptr[3] * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[4] * sqrt_gUxx 
		- gamma_uu.xy * input.ptr[4] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * input.ptr[29] 
		- gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * input.ptr[6] 
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[29] * 1. / sqrt_gUxx 
		+ gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[6] * 1. / sqrt_gUxx 
		+ gamma_uu.yy * gamma_uu.yy * input.ptr[25] * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.yy * input.ptr[25] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		+ gamma_uu.yy * gamma_uu.yy * input.ptr[2] * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.yy * input.ptr[2] * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy 
		- input.ptr[8] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy 
		- input.ptr[8] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		+ 2. * input.ptr[8] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy
	) / (
		-gamma_uu.xx * (
			+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
			+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy
		)
	);
	
	resultU.d_lll.x.xz = -(
		+ gamma_uu.xy * input.ptr[29] * 1. / sqrt_gUxx 
		- gamma_uu.xy * input.ptr[6] * 1. / sqrt_gUxx 
		+ gamma_uu.xz * input.ptr[28] * 1. / sqrt_gUxx 
		- gamma_uu.xz * input.ptr[5] * 1. / sqrt_gUxx 
		+ input.ptr[1] 
		+ input.ptr[24] 
		- input.ptr[9] * gamma_uu.xx
	) / gamma_uu.xx;
	
	resultU.d_lll.x.yy = (
		+ 2. * gamma_uu.xx * gamma_uu.yz * input.ptr[29] * gamma_uu.xy * gamma_uu.xy 
		- 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * input.ptr[29] * gamma_uu.yy 
		- 2. * gamma_uu.xx * gamma_uu.yz * input.ptr[6] * gamma_uu.xy * gamma_uu.xy 
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * input.ptr[6] * gamma_uu.yy 
		+ gamma_uu.xx * gamma_uu.zz * input.ptr[28] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * input.ptr[28] * gamma_uu.yy 
		- gamma_uu.xx * gamma_uu.zz * input.ptr[5] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * input.ptr[5] * gamma_uu.yy 
		- gamma_uu.xx * sqrt_gUxx * input.ptr[26] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * input.ptr[26] * gamma_uu.yy 
		- gamma_uu.xx * input.ptr[27] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * input.ptr[27] * gamma_uu.yy 
		- gamma_uu.xx * sqrt_gUxx * input.ptr[3] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * input.ptr[3] * gamma_uu.yy 
		+ gamma_uu.xx * input.ptr[4] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * gamma_uu.xx * input.ptr[4] * gamma_uu.yy 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[29] 
		+ 2. * gamma_uu.xy * gamma_uu.xz * input.ptr[29] * gamma_uu.xx * gamma_uu.yy 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[6] 
		- 2. * gamma_uu.xy * gamma_uu.xz * input.ptr[6] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[25] * sqrt_gUxx 
		+ gamma_uu.xy * input.ptr[25] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[2] * sqrt_gUxx 
		+ gamma_uu.xy * input.ptr[2] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		+ gamma_uu.xz * input.ptr[14] * 1. / sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * gamma_uu.xz * input.ptr[14] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
		+ gamma_uu.xz * input.ptr[14] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy 
		- gamma_uu.xz * input.ptr[19] * 1. / sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		+ 2. * gamma_uu.xz * input.ptr[19] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
		- gamma_uu.xz * input.ptr[19] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy 
		+ gamma_uu.xz * input.ptr[1] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xz * input.ptr[1] * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx 
		- gamma_uu.xz * input.ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xz * input.ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xz * input.ptr[28] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xz * gamma_uu.xz * input.ptr[28] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xz * gamma_uu.xz * input.ptr[5] * gamma_uu.xx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xz * input.ptr[5] * gamma_uu.xy * gamma_uu.xy 
		+ input.ptr[11] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * input.ptr[11] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
		+ input.ptr[11] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy
	) / (
		sqrt_gUxx * (
			+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
			+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy
		)
	);
	
	resultU.d_lll.x.yz = -.5 * (
		+ gamma_uu.xy * input.ptr[14] * 1. / sqrt_gUxx 
		+ gamma_uu.xz * input.ptr[20] * 1. / sqrt_gUxx 
		- gamma_uu.xy * input.ptr[19] * 1. / sqrt_gUxx 
		- gamma_uu.xz * input.ptr[15] * 1. / sqrt_gUxx 
		
		- input.ptr[12] * sqrt_gUxx 
		- input.ptr[17] * sqrt_gUxx 
		- 2. * input.ptr[29] 
		+ 2. * input.ptr[6]
	) * _1_sqrt_gUxx;
	
	resultU.d_lll.x.zz = -(
		+ gamma_uu.xy * input.ptr[15] * 1. / sqrt_gUxx 
		- gamma_uu.xy * input.ptr[20] * 1. / sqrt_gUxx 
		- input.ptr[18] * sqrt_gUxx 
		- input.ptr[28] 
		+ input.ptr[5]
	) * _1_sqrt_gUxx;
	
	resultU.d_lll.y.xx = input.ptr[10];
	resultU.d_lll.y.xy = input.ptr[11];
	resultU.d_lll.y.xz = input.ptr[12];
	resultU.d_lll.y.yy = input.ptr[13];
	resultU.d_lll.y.yz = input.ptr[14];
	resultU.d_lll.y.zz = input.ptr[15];
	resultU.d_lll.z.xx = input.ptr[16];
	resultU.d_lll.z.xy = input.ptr[17];
	resultU.d_lll.z.xz = input.ptr[18];
	resultU.d_lll.z.yy = input.ptr[19];
	resultU.d_lll.z.yz = input.ptr[20];
	resultU.d_lll.z.zz = input.ptr[21];
	
	resultU.K_ll.xx = (
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[26] 
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[3] 
		
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[26] * f
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[3] * f 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[29] * f 
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[29] * gamma_uu.xx * gamma_uu.xx * f
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[29] * gamma_uu.xx * f
		
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[29] 
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[29] * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[29] * gamma_uu.xx * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[6] 
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[6] * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[6] * gamma_uu.xx * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * input.ptr[6] * f 
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[6] * f * gamma_uu.xx * gamma_uu.xx 
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[6] * f * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * input.ptr[25] * gamma_uu.xx * sqrt_gUxx 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * input.ptr[25] * sqrt_gUxx 
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[25] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * input.ptr[25] * sqrt_gUxx * f 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * input.ptr[25] * gamma_uu.xx * sqrt_gUxx * f 
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[25] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f 
		+ gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[2] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * input.ptr[2] * gamma_uu.xx * sqrt_gUxx 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * input.ptr[2] * sqrt_gUxx 
		- gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[2] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * input.ptr[2] * f * sqrt_gUxx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * input.ptr[2] * f * gamma_uu.xx * sqrt_gUxx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * gamma_uu.yy * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * gamma_uu.yy * gamma_uu.yy * f * gamma_uu.xx * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * f 
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * f * gamma_uu.yy * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * gamma_uu.yy * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * gamma_uu.yy * gamma_uu.yy * f * gamma_uu.xx * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * f 
		+ 4. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * f * gamma_uu.yy * gamma_uu.xx 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * gamma_uu.yy * gamma_uu.xx 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * f 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[28] * f * gamma_uu.yy * gamma_uu.xx 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * gamma_uu.yy * gamma_uu.xx 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * f 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.zz * input.ptr[5] * f * gamma_uu.yy * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * input.ptr[27] 
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * input.ptr[27] * gamma_uu.yy * gamma_uu.xx 
		+ 4. * gamma_uu.xy * gamma_uu.xy * f * input.ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * input.ptr[4] 
		+ 4. * gamma_uu.xy * gamma_uu.xy * f * input.ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * input.ptr[4] * gamma_uu.yy * gamma_uu.xx 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * sqrt_gUxx 
		- 4. * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy 
		+ 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		+ 4. * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * f 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * sqrt_gUxx * f 
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] * gamma_uu.xx * sqrt_gUxx * f * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[27] * gamma_uu.yy * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * input.ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * sqrt_gUxx 
		- 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		+ 4. * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * f * sqrt_gUxx 
		- 4. * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy 
		+ 5. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[3] * f * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[4] * gamma_uu.yy * gamma_uu.xx 
		+ 2. * gamma_uu.xy * gamma_uu.xy * input.ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[1] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx 
		- 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[1] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * sqrt_gUxx 
		+ gamma_uu.xz * gamma_uu.yy * input.ptr[1] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx 
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[1] * f * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx 
		+ 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[1] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * sqrt_gUxx 
		- gamma_uu.xz * gamma_uu.yy * input.ptr[1] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx 
		- gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[24] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx 
		+ 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xz * gamma_uu.yy * input.ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[24] * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * f 
		+ gamma_uu.xz * gamma_uu.yy * input.ptr[24] * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f 
		+ gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[24] * gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * f 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[28] * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[28] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		+ 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[28] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[28] * f * gamma_uu.xx * gamma_uu.xx 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[28] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[28] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[5] * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[5] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		+ 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[5] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[5] * f * gamma_uu.xx * gamma_uu.xx 
		+ gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * input.ptr[5] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * gamma_uu.xz * gamma_uu.xz * gamma_uu.yy * gamma_uu.yy * input.ptr[5] * f * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		- gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- f * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- f * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- f * solver->m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[27] 
		- 3. * f * solver->m * gamma_uu.xy * gamma_uu.xy * input.ptr[27] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ 3. * f * solver->m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[27] * gamma_uu.yy * gamma_uu.xx 
		- f * solver->m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[4] 
		- 3. * f * solver->m * gamma_uu.xy * gamma_uu.xy * input.ptr[4] * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		+ 3. * f * solver->m * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[4] * gamma_uu.yy * gamma_uu.xx 
		+ f * solver->m * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[27] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		+ f * solver->m * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * input.ptr[4] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- input.ptr[0] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy 
		+ input.ptr[0] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f 
		+ input.ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		+ 3. * input.ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- 3. * input.ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- input.ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.xx 
		+ 3. * input.ptr[0] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- 3. * input.ptr[0] * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- input.ptr[30] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy 
		+ input.ptr[30] * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f 
		+ input.ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx 
		+ 3. * input.ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		- 3. * input.ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx 
		- input.ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.xx 
		- 3. * input.ptr[30] * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx * gamma_uu.xx 
		+ 3. * input.ptr[30] * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.xx * gamma_uu.xx
	) / (
		+ gamma_uu.xx * (
			+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- 3. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
			+ 3. * gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy * gamma_uu.yy 
			- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f 
			+ 3. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy 
			- 3. * gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * f * gamma_uu.yy * gamma_uu.yy 
			- gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy 
			+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * gamma_uu.yy * f
		)
	);
		
	resultU.K_ll.xy = (2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[29] 
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * input.ptr[29] * gamma_uu.yy 
		+ 2. * sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.yz * input.ptr[6] 
		- 2. * gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * gamma_uu.yz * input.ptr[6] * gamma_uu.yy 
		- gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[26] 
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * input.ptr[26] * gamma_uu.yy 
		+ gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[3] 
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.xy * input.ptr[3] * gamma_uu.yy 
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * input.ptr[25] 
		- gamma_uu.xx * gamma_uu.yy * input.ptr[25] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy * input.ptr[2] 
		+ gamma_uu.xx * gamma_uu.yy * input.ptr[2] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * input.ptr[28] 
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * input.ptr[28] * gamma_uu.yy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * sqrt_gUxx * gamma_uu.zz * input.ptr[5] 
		- gamma_uu.xy * gamma_uu.xx * sqrt_gUxx * gamma_uu.zz * input.ptr[5] * gamma_uu.yy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[1] 
		- gamma_uu.xy * gamma_uu.xz * input.ptr[1] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[24] 
		+ gamma_uu.xy * gamma_uu.xz * input.ptr[24] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[28] * 1. / sqrt_gUxx 
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[28] * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[5] * 1. / sqrt_gUxx 
		+ gamma_uu.xy * gamma_uu.xz * gamma_uu.xz * input.ptr[5] * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[27] * sqrt_gUxx 
		+ gamma_uu.xy * input.ptr[27] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		- gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[4] * sqrt_gUxx 
		+ gamma_uu.xy * input.ptr[4] * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * input.ptr[29] 
		+ gamma_uu.xz * gamma_uu.xx * sqrt_gUxx * gamma_uu.yy * gamma_uu.yy * input.ptr[6] 
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[29] * 1. / sqrt_gUxx 
		- gamma_uu.xz * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[6] * 1. / sqrt_gUxx
	) / (
		-sqrt_gUxx * (
			+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
			- 2. * gamma_uu.xx * gamma_uu.xy * gamma_uu.xy * gamma_uu.yy 
			+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy
		)
	);
	
	resultU.K_ll.xz = -(
		+ (
			+ gamma_uu.xy * input.ptr[29]
			+ gamma_uu.xy * input.ptr[6]
			+ gamma_uu.xz * input.ptr[28]
			+ gamma_uu.xz * input.ptr[5]
		) * _1_sqrt_gUxx
		+ input.ptr[24]
		- input.ptr[1] 
	) * _1_sqrt_gUxx;
	
	resultU.K_ll.yy = -(
		+ sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[25] 
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * input.ptr[25] * gamma_uu.yy 
		- sqrt_gUxx * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * input.ptr[2] 
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.xy * input.ptr[2] * gamma_uu.yy 
		- sqrt_gUxx * gamma_uu.xz * input.ptr[1] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * sqrt_gUxx * gamma_uu.xz * input.ptr[1] * gamma_uu.yy 
		+ sqrt_gUxx * gamma_uu.xz * input.ptr[24] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * sqrt_gUxx * gamma_uu.xz * input.ptr[24] * gamma_uu.yy 
		- 2. * gamma_uu.xx * gamma_uu.yz * input.ptr[29] * gamma_uu.xy * gamma_uu.xy 
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * input.ptr[29] * gamma_uu.yy 
		- 2. * gamma_uu.xx * gamma_uu.yz * input.ptr[6] * gamma_uu.xy * gamma_uu.xy 
		+ 2. * gamma_uu.xx * gamma_uu.xx * gamma_uu.yz * input.ptr[6] * gamma_uu.yy 
		- gamma_uu.xx * gamma_uu.zz * input.ptr[28] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * input.ptr[28] * gamma_uu.yy 
		- gamma_uu.xx * gamma_uu.zz * input.ptr[5] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.zz * input.ptr[5] * gamma_uu.yy 
		+ gamma_uu.xx * sqrt_gUxx * input.ptr[26] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * input.ptr[26] * gamma_uu.yy 
		+ gamma_uu.xx * input.ptr[27] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * gamma_uu.xx * input.ptr[27] * gamma_uu.yy 
		- gamma_uu.xx * sqrt_gUxx * input.ptr[3] * gamma_uu.xy * gamma_uu.xy 
		+ gamma_uu.xx * gamma_uu.xx * sqrt_gUxx * input.ptr[3] * gamma_uu.yy 
		+ gamma_uu.xx * input.ptr[4] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xx * gamma_uu.xx * input.ptr[4] * gamma_uu.yy 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[29] 
		- 2. * gamma_uu.xy * gamma_uu.xz * input.ptr[29] * gamma_uu.xx * gamma_uu.yy 
		+ 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xz * input.ptr[6] 
		- 2. * gamma_uu.xy * gamma_uu.xz * input.ptr[6] * gamma_uu.xx * gamma_uu.yy 
		- gamma_uu.xz * gamma_uu.xz * input.ptr[28] * gamma_uu.xx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xz * input.ptr[28] * gamma_uu.xy * gamma_uu.xy 
		- gamma_uu.xz * gamma_uu.xz * input.ptr[5] * gamma_uu.xx * gamma_uu.yy 
		+ gamma_uu.xz * gamma_uu.xz * input.ptr[5] * gamma_uu.xy * gamma_uu.xy
	) / (
		+ gamma_uu.xx * gamma_uu.xx * gamma_uu.yy * gamma_uu.yy 
		+ gamma_uu.xy * gamma_uu.xy * gamma_uu.xy * gamma_uu.xy 
		- 2. * gamma_uu.xy * gamma_uu.xy * gamma_uu.xx * gamma_uu.yy
	);
	
	resultU.K_ll.yz = input.ptr[29] + input.ptr[6];
	resultU.K_ll.zz = input.ptr[28] + input.ptr[5];
	resultU.Theta = input.ptr[27] + input.ptr[4];
	
	resultU.Z_l.x = (
		- input.ptr[22] * gamma_uu.xy 
		- input.ptr[23] * gamma_uu.xz 
		+ input.ptr[26] * gamma_uu.xx 
		+ input.ptr[3] * gamma_uu.xx
	) * _1_gammaUxx;
	
	resultU.Z_l.y = input.ptr[2] + input.ptr[22] + input.ptr[25];
	resultU.Z_l.z = input.ptr[1] + input.ptr[23] + input.ptr[24];

	return resultU;
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t input,
	real3 pt 
) {
	cons_t results;
	for (int i = 0; i < numStates; ++i) {
		results.ptr[i] = 0;
	}
	
	sym3 gamma_ll = sym3_swap<?=side?>(eig.gamma_ll);
	sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);
	
	real sqrt_f = eig.sqrt_f;
	real f = sqrt_f * sqrt_f;
	
	results.ptr[0] = f * (gamma_uu.xx * input.ptr[21] 
		+ 2. * gamma_uu.xy * input.ptr[22] 
		+ 2. * gamma_uu.xz * input.ptr[23] 
		+ gamma_uu.yy * input.ptr[24] 
		+ 2. * gamma_uu.yz * input.ptr[25] 
		+ gamma_uu.zz * input.ptr[26] 
		- solver->m * input.ptr[27]);
	results.ptr[1] = 0.;
	results.ptr[2] = 0.;
	results.ptr[3] = input.ptr[21];
	results.ptr[4] = input.ptr[22];
	results.ptr[5] = input.ptr[23];
	results.ptr[6] = input.ptr[24];
	results.ptr[7] = input.ptr[25];
	results.ptr[8] = input.ptr[26];
	results.ptr[9] = 0.;
	results.ptr[10] = 0.;
	results.ptr[11] = 0.;
	results.ptr[12] = 0.;
	results.ptr[13] = 0.;
	results.ptr[14] = 0.;
	results.ptr[15] = 0.;
	results.ptr[16] = 0.;
	results.ptr[17] = 0.;
	results.ptr[18] = 0.;
	results.ptr[19] = 0.;
	results.ptr[20] = 0.;
	results.ptr[21] = 
		-(gamma_uu.yy * input.ptr[10] 
		- gamma_uu.yy * input.ptr[6] 
		+ gamma_uu.yz * input.ptr[11] 
		+ gamma_uu.yz * input.ptr[16] 
		- 2. * gamma_uu.yz * input.ptr[7] 
		+ gamma_uu.zz * input.ptr[17] 
		- gamma_uu.zz * input.ptr[8] 
		- input.ptr[0] 
		+ 2. * input.ptr[28]);
	results.ptr[22] = (2. * gamma_uu.xy * input.ptr[10] 
		- 2. * gamma_uu.xy * input.ptr[6] 
		+ gamma_uu.xz * input.ptr[11] 
		+ gamma_uu.xz * input.ptr[16] 
		- 2. * gamma_uu.xz * input.ptr[7] 
		+ gamma_uu.yz * input.ptr[13] 
		- gamma_uu.yz * input.ptr[18] 
		+ gamma_uu.zz * input.ptr[14] 
		- gamma_uu.zz * input.ptr[19] 
		+ input.ptr[1] 
		- 2. * input.ptr[29]) / 2.;
	results.ptr[23] = (gamma_uu.xy * input.ptr[11] 
		+ gamma_uu.xy * input.ptr[16] 
		- 2. * gamma_uu.xy * input.ptr[7] 
		+ 2. * gamma_uu.xz * input.ptr[17] 
		- 2. * gamma_uu.xz * input.ptr[8] 
		- gamma_uu.yy * input.ptr[13] 
		+ gamma_uu.yy * input.ptr[18] 
		- gamma_uu.yz * input.ptr[14] 
		+ gamma_uu.yz * input.ptr[19] 
		+ input.ptr[2] 
		- 2. * input.ptr[30]) / 2.;
	results.ptr[24] = 
		-(gamma_uu.xx * input.ptr[10] 
		- gamma_uu.xx * input.ptr[6] 
		+ gamma_uu.xz * input.ptr[13] 
		- gamma_uu.xz * input.ptr[18]);
	results.ptr[25] = (
		-(gamma_uu.xx * input.ptr[11] 
		+ gamma_uu.xx * input.ptr[16] 
		- 2. * gamma_uu.xx * input.ptr[7] 
		- gamma_uu.xy * input.ptr[13] 
		+ gamma_uu.xy * input.ptr[18] 
		+ gamma_uu.xz * input.ptr[14] 
		- gamma_uu.xz * input.ptr[19])) / 2.;
	results.ptr[26] = 
		-(gamma_uu.xx * input.ptr[17] 
		- gamma_uu.xx * input.ptr[8] 
		- gamma_uu.xy * input.ptr[14] 
		+ gamma_uu.xy * input.ptr[19]);
	results.ptr[27] = 
		-(gamma_uu.xx * gamma_uu.yy * input.ptr[10] 
		- gamma_uu.xx * gamma_uu.yy * input.ptr[6] 
		+ gamma_uu.xx * gamma_uu.yz * input.ptr[11] 
		+ gamma_uu.xx * gamma_uu.yz * input.ptr[16] 
		- 2. * gamma_uu.xx * gamma_uu.yz * input.ptr[7] 
		+ gamma_uu.xx * gamma_uu.zz * input.ptr[17] 
		- gamma_uu.xx * gamma_uu.zz * input.ptr[8] 
		+ gamma_uu.xx * input.ptr[28] 
		- gamma_uu.xy * gamma_uu.xz * input.ptr[11] 
		- gamma_uu.xy * gamma_uu.xz * input.ptr[16] 
		+ 2. * gamma_uu.xy * gamma_uu.xz * input.ptr[7] 
		- gamma_uu.xy * gamma_uu.yz * input.ptr[13] 
		+ gamma_uu.xy * gamma_uu.yz * input.ptr[18] 
		- gamma_uu.xy * gamma_uu.zz * input.ptr[14] 
		+ gamma_uu.xy * gamma_uu.zz * input.ptr[19] 
		- gamma_uu.xy * gamma_uu.xy * input.ptr[10] 
		+ gamma_uu.xy * input.ptr[29] 
		+ gamma_uu.xy * gamma_uu.xy * input.ptr[6] 
		+ gamma_uu.xz * gamma_uu.yy * input.ptr[13] 
		- gamma_uu.xz * gamma_uu.yy * input.ptr[18] 
		+ gamma_uu.xz * gamma_uu.yz * input.ptr[14] 
		- gamma_uu.xz * gamma_uu.yz * input.ptr[19] 
		- gamma_uu.xz * gamma_uu.xz * input.ptr[17] 
		+ gamma_uu.xz * input.ptr[30] 
		+ gamma_uu.xz * gamma_uu.xz * input.ptr[8]);
	results.ptr[28] = gamma_uu.xy * input.ptr[22] 
		+ gamma_uu.xz * input.ptr[23] 
		+ gamma_uu.yy * input.ptr[24] 
		+ 2. * gamma_uu.yz * input.ptr[25] 
		+ gamma_uu.zz * input.ptr[26] 
		- input.ptr[27];
	results.ptr[29] = 
		-(gamma_uu.xx * input.ptr[22] 
		+ gamma_uu.xy * input.ptr[24] 
		+ gamma_uu.xz * input.ptr[25]);
	results.ptr[30] = 
		-(gamma_uu.xx * input.ptr[23] 
		+ gamma_uu.xy * input.ptr[25] 
		+ gamma_uu.xz * input.ptr[26]);

	return results;
}
<? end ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf)
{
#if 0
	SETBOUNDS_NOGHOST();
	const global cons_t* U = UBuf + index;
	global cons_t* deriv = derivBuf + index;

	const real xi = 1.;	//which is which parameter?  I alwasy forget .. m .. lambda ... xi ... always changing names

	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real f = calc_f(U->alpha);	//could be based on alpha...

	//K^i_j = gamma^ik K_kj
	real3x3 K_ul = (real3x3){
<? 
for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<? 	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
		?> + gamma_uu.<?=sym(i,k)?> * U->K_ll.<?=sym(k,j)?><? 
		end ?>,
<? 	end
?>		},
<? 
end
?>	};

	real tr_K = real3x3_trace(K_ul);

	//TODO correct source terms
	//I'm taking this from the 2008 Yano et al "Flux-Vector Splitting..."
	//which itself references (for the source terms) the 2005 Bona et al "Geometrically Motivated ..." 
	//but notice, 2005 Bona paper shows the flux as densitized, so the eigenvalues without any gamma_ll influence

	//alpha,t + ... = alpha beta^k a_k - alpha^2 f (K_ll - 2 Theta)
	deriv->alpha += -U->alpha * U->alpha * f * (tr_K - 2. * U->Theta);

	//a_i,t + ... = b_i^k a_k - b_k^k a_i

	//beta^i_,t + ... = beta^k b_k^i - alpha Q^i
	//Q^i = alpha (a^i - d_lll^i + 2 V^i)
	//V^i = d_lll^i - e^i - Z^i
	
	//gamma_ij,t + ... = 2 beta^k d_kij + b_ji + b_ij - 2 alpha K_ij
	deriv->gamma_ll = sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));

	//d_kij,t + ... = b_k^l d_lij - b^l d_kij

	//conn_ijk = .5 * (g_ij,k + g_ik,j - g_jk,i) 
	//= d_kij + d_jik - d_ijk
	_3sym3 conn_lll = {
<? 
for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<? 	for jk,xjk in ipairs(symNames) do 
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
?>			.<?=xjk?> = U->d_lll.<?=xk?>.<?=sym(i,j)?> + U->d_lll.<?=xj?>.<?=sym(i,k)?> - U->d_lll.<?=xi?>.<?=sym(j,k)?>,
<? 	end 
?>		},
<? 
end
?>	};

	_3sym3 conn_ull = sym3_3sym3_mul(gamma_uu, conn_lll);
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//conn^ij_k = gamma^jl conn^i_lk
	_3sym3 conn_uul = {
<? 
	for i,xi in ipairs(xNames) do
?>
		.<?=xi?> = {
<? 
		for jk,xjk in ipairs(symNames) do
			local j,k = from6to3x3(jk)
?>			.<?=xjk?> = 0. <?
			for l,xl in ipairs(xNames) do
				?> + gamma_uu.<?=sym(l,j)?> * conn_ull.<?=xi?>.<?=sym(l,k)?><?
			end	?>,
<? end		
?>		},
<? 
	end
?>	};

	//d_ijk conn^ij_l gamma^kl
	real d_dot_conn = 0.<?
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
?> + U->d_lll.<?=xi?>.<?=sym(j,k)?> * conn_uul.<?=xi?>.<?=sym(j,k)?> * gamma_uu.<?=sym(k,l)?><?
			end
		end
	end
end
?>;

	//d_i = d_ik^k
	real3 d_lll = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(U->d_lll.<?=xi?>, gamma_uu),
<? end
?>	};
	
	//e_i = d_lll^k_ki
	real3 e = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0. <?
	for j,xj in ipairs(xNames) do
?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end
?>,
<? end
?>	};

	real3 d_u = sym3_real3_mul(gamma_uu, d_lll);
	real3 e_u = sym3_real3_mul(gamma_uu, e);
	real3 a_u = sym3_real3_mul(gamma_uu, U->a_l);
	real3 Z_u = sym3_real3_mul(gamma_uu, U->Z_l);

	//connSq_ij = conn^k_li conn^l_kj
	real connSq_ll[3][3] = {
<? for i,xi in ipairs(xNames) do 
?>		{
<?	for j,xj in ipairs(xNames) do
?>
			0. <?
		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				?> + conn_ull.<?=xk?>.<?=sym(l,i)?> * conn_ull.<?=xl?>.<?=sym(k,j)?><?
			end
		end
	?>,
<? end 
?>		},
<? end
?>	};

	//K_ik K^k_j
	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);

	real tr_KSq = sym3_dot(KSq_ll, gamma_uu);

	//dsq_ij = d_ikl gamma^lm d_lll^k_mj
	real3x3 dsq = {
<? 
for i,xi in ipairs(xNames) do
?>
		.<?=xi?> = {
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 
				0.<?
		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
					?> + U->d_lll.<?=xi?>.<?=sym(j,l)?> * gamma_uu.<?=sym(l,m)?> * d_ull.<?=xk?>.<?=sym(m,j)?>
<?
				end
			end
		end	
?>			,
<? 	end
?>		},
<? 
end
?>	};

	/*
	K_ij,t + ... = -K_ij b_k^k + K_ik b_j^k + K_jk b_i^k 
		+ alpha ( 1/2 (1 + xi) ( - a_k Gamma^k_ij + 1/2 (a_i d_j + a_j d_i))
			+ 1/2 (1 - xi) (a_k d_lll^k_ij - 1/2 (a_j (2 e_i - d_i) + a_i (2 e_j - d_j))
				+ 2 (d_ir^m d_lll^r_mj + d_jr^m d_lll^r_mi) - 2 e_k (d_ij^k + d_ji^k)
			)
			+ (d_k + a_k - 2 Z_k) Gamma^k_ij - Gamma^k_mj Gamma^m_ki - (a_i Z_j + a_j Z_i)
			- 2 K^k_i K_kj + (K_ll - 2 Theta) K_ij
		) - 8 pi alpha (S_ij - 1/2 (S - tau) gamma_ij)
	*/
<? for ij,xij in ipairs(symNames) do 
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i], xNames[j]
?>
	deriv->K_ll.<?=xij?> += 
		//TODO shift terms
	U->alpha * (
		.5 * (1. + xi) * (
		<? for k,xk in ipairs(xNames) do ?> 
			-U->a_l.<?=xk?> * conn_ull.<?=xk?>.<?=xij?>
		<? end ?>
			+ .5 * (U->a_l.<?=xi?> * inputU.d_lll.<?=xj?> + U->a_l.<?=xj?> * inputU.d_lll.<?=xi?>)
		)
		+ .5 * (1. - xi) * (0.
		<? for k,xk in ipairs(xNames) do ?> 
			+ a_u.<?=xk?> * U->d_lll.<?=xk?>.<?=xij?>
		<? end ?>
			- .5 * ( 
				U->a_l.<?=xj?> * (2. * e.<?=xi?> - inputU.d_lll.<?=xi?>) 
				+ U->a_l.<?=xi?> * (2. * e.<?=xj?> - inputU.d_lll.<?=xj?>)
			)
			+ 2. * (dsq.v[<?=i-1?>].s[<?=j-1?>] + dsq.v[<?=j-1?>].s[<?=i-1?>])
		<? for k,xk in ipairs(xNames) do ?> 
			- 2. * e_u.<?=xk?> * (U->d_lll.<?=xi?>.<?=sym(j,k)?> + U->d_lll.<?=xj?>.<?=sym(i,k)?>)
		<? end ?>
		)
		<? for k,xk in ipairs(xNames) do ?> 
		+ (d_lll.<?=xk?> + U->a_l.<?=xk?> - 2. * U->Z_l.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		<? end ?>
		- connSq_ll[<?=j-1?>][<?=i-1?>]
		- U->a_l.<?=xi?> * U->Z_l.<?=xj?>
		- U->a_l.<?=xj?> * U->Z_l.<?=xi?>
		- 2. * KSq_ll.<?=xij?>
		+ U->K_ll.<?=xij?> * (tr_K - 2. * U->Theta)
	)
		// TODO source terms
	;
<? end
?>
		
	/*
	Theta,t + ... = 
		-Theta b_k^k 
		+ 1/2 alpha (
			2 a_k (d^k - e^k - 2 Z^k) 
			+ d_krs Gamma^krs 
			- d^k (d_k - 2 Z_k)
			- K^k_r K^r_k 
			+ K (K - 2 Theta)
		)
		- 8 pi alpha tau
	
	d_ijk conn^ijk
	= d_ijk (d^kij + d^jik - d^ijk)
	*/
	deriv->Theta += .5 * U->alpha * (0. 
		<? for k,xk in ipairs(xNames) do ?> 
		+ 2. * U->a_l.<?=xk?> * (d_u.<?=xk?> - e_u.<?=xk?> - 2. * Z_u.<?=xk?>)
		- d_u.<?=xk?> * (d_lll.<?=xk?> - 2. * U->Z_l.<?=xk?>)
		<? end ?>	
		+ d_dot_conn
		- tr_KSq
		+ tr_K * (tr_K - 2. * U->Theta)
	);

	// K^k_j Gamma^j_ki
	real3 K_times_conn = (real3){
<? 
for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + K_ul.v[<?=k-1?>].s[<?=j-1?>] * conn_ull.<?=xj?>.<?=sym(k,i)?><?
		end
	end	?>,
<? 
end
?>	};

	/*
	Z_i,t + ... = 
		-Z_i b_k^k 
		+ Z_k b_i^k 
		- 8 pi alpha S_i
		+ alpha (
			a_i (K_ll - 2 Theta) 
			- a_k K^k_i 
			- K^k_r Gamma^r_ki 
			+ K^k_i (d_k - 2 Z_k)
		)
	*/
	<? for i,xi in ipairs(xNames) do ?> 
	deriv->Z_l.<?=xi?> += U->alpha * (
		U->a_l.<?=xi?> * (tr_K - 2. * U->Theta)
		<? for k,xk in ipairs(xNames) do ?> 
		- a_u.<?=xk?> * U->K_ll.<?=sym(i,k)?>
		+ U->K_ll.<?=sym(i,k)?> * (d_u.<?=xk?> - 2. * Z_u.<?=xk?>)
		<? end ?>
		- K_times_conn.<?=xi?>
	);
	<? end ?>

#endif

}
