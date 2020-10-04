<?
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

//used by PLM
eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	
	//This is interesting, because normal_t varies based on our vector components.
	//However in my GR solvers the components are irrespective of the grid -- instead they are based on the metric of the state variables.
	normal_t n
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

range_t calcCellMinMaxEigenvalues(
	const global cons_t* U,
	real3 x,
	normal_t n
) {
	real det_gamma = sym3_det(U->gamma_ll);

	real gammaUjj;
	if (n.side == 0) {
		gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
	} else if (n.side == 1) {
		gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
	} else if (n.side == 2) {
		gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
	}
	
	real lambdaLight = U->alpha * sqrt(gammaUjj);
	
	real f = calc_f(U->alpha);
	real lambdaGauge = lambdaLight * sqrt(f);

	real lambdaMax = max(lambdaGauge, lambdaLight);
	//= lambdaLight * max(sqrt(f), 1)
	real lambdaMin = -lambdaMin;

	<? if eqn.useShift ~= 'none' then ?>
	lambdaMin -= normal_vecDotN1(n, U->beta_u);
	lambdaMax -= normal_vecDotN1(n, U->beta_u);
	<? end ?>

	return (range_t){
		.min = lambdaMin, 
		.max = lambdaMax,
	};
}

//used for interface eigen basis
eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	normal_t n
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

waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x,
	normal_t n
) {
	waves_t results;
<? if not eqn.noZeroRowsInFlux then ?>

	if (n.side == 0) {

		//a_y, a_z
		results.ptr[6] = inputU.a_l.y;
		results.ptr[7] = inputU.a_l.z;

		//d_yij
		results.ptr[8] = inputU.d_lll.y.xx;
		results.ptr[9] = inputU.d_lll.y.xy;
		results.ptr[10] = inputU.d_lll.y.xz;
		results.ptr[11] = inputU.d_lll.y.yy;
		results.ptr[12] = inputU.d_lll.y.yz;
		results.ptr[13] = inputU.d_lll.y.zz;

		//d_zij
		results.ptr[14] = inputU.d_lll.z.xx;
		results.ptr[15] = inputU.d_lll.z.xy;
		results.ptr[16] = inputU.d_lll.z.xz;
		results.ptr[17] = inputU.d_lll.z.yy;
		results.ptr[18] = inputU.d_lll.z.yz;
		results.ptr[19] = inputU.d_lll.z.zz;

		//V_j
		results.ptr[20] = inputU.V_l.x;
		results.ptr[21] = inputU.V_l.y;
		results.ptr[22] = inputU.V_l.z;
		
		sym3 K_sqrt_gammaUxx = sym3_real_mul(inputU.K_ll, eig.sqrt_gammaUjj.x);

		//a^x - f d^xj_j

		real f = eig.sqrt_f * eig.sqrt_f;
		real d_x_input = sym3_dot(eig.gamma_uu, inputU.d_lll.x);
		results.ptr[23] = inputU.a_l.x - f * d_x_input;

		//gauge:
		//sqrt(f gamma^xx) K +- (a^x + 2 V^x)

		real ev0a = eig.sqrt_f * sym3_dot(eig.gamma_uu, K_sqrt_gammaUxx);
		real ev0b = eig.gamma_uu.xx * (inputU.a_l.x + 2. * inputU.V_l.x) 
					+ eig.gamma_uu.xy * (inputU.a_l.y + 2. * inputU.V_l.y)
					+ eig.gamma_uu.xz * (inputU.a_l.z + 2. * inputU.V_l.z);
		results.ptr[0] = ev0a - ev0b;
		results.ptr[29] = ev0a + ev0b;

		//light:
		//sqrt(gamma^xx) K_xy +- (d^x_xy + .5 (a_y - d_yj^j) + V_y)

		real d_y_input = sym3_dot(eig.gamma_uu, inputU.d_lll.y);
		real dUx_xy_input = eig.gamma_uu.xx * inputU.d_lll.x.xy + eig.gamma_uu.xy * inputU.d_lll.y.xy + eig.gamma_uu.xz * inputU.d_lll.z.xy;
		real ev1b = .5 * (inputU.a_l.y - d_y_input) + inputU.V_l.y + dUx_xy_input;
		results.ptr[1] = K_sqrt_gammaUxx.xy - ev1b;
		results.ptr[24] = K_sqrt_gammaUxx.xy + ev1b;

		//light:
		//sqrt(gamma^xx) K_xz +- (d^x_xz + .5 (a_z - d_zj^j) + V_z)

		real d_z_input = sym3_dot(eig.gamma_uu, inputU.d_lll.z);
		real dUx_xz_input = eig.gamma_uu.xx * inputU.d_lll.x.xz + eig.gamma_uu.xy * inputU.d_lll.y.xz + eig.gamma_uu.xz * inputU.d_lll.z.xz;
		real ev2b = .5 * (inputU.a_l.z - d_z_input) + inputU.V_l.z + dUx_xz_input;
		results.ptr[2] = K_sqrt_gammaUxx.xz - ev2b;
		results.ptr[25] = K_sqrt_gammaUxx.xz + ev2b;

		//light:
		//sqrt(gamma^xx) K_yy +- d^x_yy

		real dUx_yy_input = eig.gamma_uu.xx * inputU.d_lll.x.yy + eig.gamma_uu.xy * inputU.d_lll.y.yy + eig.gamma_uu.xz * inputU.d_lll.z.yy;
		results.ptr[3] = K_sqrt_gammaUxx.yy - dUx_yy_input;
		results.ptr[26] = K_sqrt_gammaUxx.yy + dUx_yy_input;

		//light:
		//sqrt(gamma^xx) K_yz +- d^x_yz

		real dUx_yz_input = eig.gamma_uu.xx * inputU.d_lll.x.yz + eig.gamma_uu.xy * inputU.d_lll.y.yz + eig.gamma_uu.xz * inputU.d_lll.z.yz;
		results.ptr[4] = K_sqrt_gammaUxx.yz - dUx_yz_input; 
		results.ptr[27] = K_sqrt_gammaUxx.yz + dUx_yz_input;

		//light:
		//sqrt(gamma^xx) K_zz +- d^x_zz

		real dUx_zz_input = eig.gamma_uu.xx * inputU.d_lll.x.zz + eig.gamma_uu.xy * inputU.d_lll.y.zz + eig.gamma_uu.xz * inputU.d_lll.z.zz;
		results.ptr[5] = K_sqrt_gammaUxx.zz - dUx_zz_input;
		results.ptr[28] = K_sqrt_gammaUxx.zz + dUx_zz_input;

	} else if (n.side == 1) {

		//a_x, a_z
		results.ptr[6] = inputU.a_l.x;
		results.ptr[7] = inputU.a_l.z;

		//d_xij
		results.ptr[8] = inputU.d_lll.x.xx;
		results.ptr[9] = inputU.d_lll.x.xy;
		results.ptr[10] = inputU.d_lll.x.xz;
		results.ptr[11] = inputU.d_lll.x.yy;
		results.ptr[12] = inputU.d_lll.x.yz;
		results.ptr[13] = inputU.d_lll.x.zz;
		
		//d_zij
		results.ptr[14] = inputU.d_lll.z.xx;
		results.ptr[15] = inputU.d_lll.z.xy;
		results.ptr[16] = inputU.d_lll.z.xz;
		results.ptr[17] = inputU.d_lll.z.yy;
		results.ptr[18] = inputU.d_lll.z.yz;
		results.ptr[19] = inputU.d_lll.z.zz;
		
		//V_j
		results.ptr[20] = inputU.V_l.x;
		results.ptr[21] = inputU.V_l.y;
		results.ptr[22] = inputU.V_l.z;
		
		sym3 K_sqrt_gammaUyy = sym3_real_mul(inputU.K_ll, eig.sqrt_gammaUjj.y);

		//a^y - f d^yj_j

		real f = eig.sqrt_f * eig.sqrt_f;
		real d_y_input = sym3_dot(eig.gamma_uu, inputU.d_lll.y);
		results.ptr[23] = inputU.a_l.y - f * d_y_input;
		
		//gauge:
		//sqrt(f gamma^yy) K +- (a^y + 2 V^y)

		real ev0a = eig.sqrt_f * sym3_dot(eig.gamma_uu, K_sqrt_gammaUyy);
		real ev0b = eig.gamma_uu.xy * (inputU.a_l.x + 2. * inputU.V_l.x)
					+ eig.gamma_uu.yy * (inputU.a_l.y + 2. * inputU.V_l.y)
					+ eig.gamma_uu.yz * (inputU.a_l.z + 2. * inputU.V_l.z);
		results.ptr[0] = ev0a - ev0b;
		results.ptr[29] = ev0a + ev0b;

		//light:
		//sqrt(gamma^yy) K_xx +- d^y_xx

		real dUy_xx_input = eig.gamma_uu.xy * inputU.d_lll.x.xx + eig.gamma_uu.yy * inputU.d_lll.y.xx + eig.gamma_uu.yz * inputU.d_lll.z.xx;
		results.ptr[1] = K_sqrt_gammaUyy.xx - dUy_xx_input;
		results.ptr[24] = K_sqrt_gammaUyy.xx + dUy_xx_input;

		//light:
		//sqrt(gamma^yy) K_xy +- (d^y_xy + .5 (a_x - d_xj^j) + V_x)

		real d_x_input = sym3_dot(eig.gamma_uu, inputU.d_lll.x);
		real dUy_xy_input = eig.gamma_uu.xy * inputU.d_lll.x.xy + eig.gamma_uu.yy * inputU.d_lll.y.xy + eig.gamma_uu.yz * inputU.d_lll.z.xy;
		real ev2b = dUy_xy_input + .5 * (inputU.a_l.x - d_x_input) + inputU.V_l.x;
		results.ptr[2] = K_sqrt_gammaUyy.xy - ev2b;
		results.ptr[25] = K_sqrt_gammaUyy.xy + ev2b;

		//light:
		//sqrt(gamma^yy) K_xz +- d^y_xz

		real dUy_xz_input = eig.gamma_uu.xy * inputU.d_lll.x.xz + eig.gamma_uu.yy * inputU.d_lll.y.xz + eig.gamma_uu.yz * inputU.d_lll.z.xz;
		results.ptr[3] = K_sqrt_gammaUyy.xz - dUy_xz_input;
		results.ptr[26] = K_sqrt_gammaUyy.xz + dUy_xz_input;

		//light:
		//sqrt(gamma^yy) K_yz +- (d^y_yz + .5 (a_z - d_zj^j) + V_z)

		real dUy_yz_input = eig.gamma_uu.xy * inputU.d_lll.x.yz + eig.gamma_uu.yy * inputU.d_lll.y.yz + eig.gamma_uu.yz * inputU.d_lll.z.yz;
		real d_z_input = sym3_dot(eig.gamma_uu, inputU.d_lll.z);
		real ev4b = dUy_yz_input + .5 * (inputU.a_l.z - d_z_input) + inputU.V_l.z;
		results.ptr[4] = K_sqrt_gammaUyy.yz - ev4b;
		results.ptr[27] = K_sqrt_gammaUyy.yz + ev4b;

		//light:
		//sqrt(gamma^yy) K_zz +- d^y_zz

		real dUy_zz_input = eig.gamma_uu.xy * inputU.d_lll.x.zz + eig.gamma_uu.yy * inputU.d_lll.y.zz + eig.gamma_uu.yz * inputU.d_lll.z.zz;
		results.ptr[5] = K_sqrt_gammaUyy.zz - dUy_zz_input;
		results.ptr[28] = K_sqrt_gammaUyy.zz - dUy_zz_input;
	
	} else if (n.side == 2) {

		//a_x, a_y
		results.ptr[6] = inputU.a_l.x;
		results.ptr[7] = inputU.a_l.y;
		
		//d_xij
		results.ptr[8] =  inputU.d_lll.x.xx;
		results.ptr[9] =  inputU.d_lll.x.xy;
		results.ptr[10] = inputU.d_lll.x.xz;
		results.ptr[11] = inputU.d_lll.x.yy;
		results.ptr[12] = inputU.d_lll.x.yz;
		results.ptr[13] = inputU.d_lll.x.zz;
		
		//d_yij
		results.ptr[14] = inputU.d_lll.y.xx;
		results.ptr[15] = inputU.d_lll.y.xy;
		results.ptr[16] = inputU.d_lll.y.xz;
		results.ptr[17] = inputU.d_lll.y.yy;
		results.ptr[18] = inputU.d_lll.y.yz;
		results.ptr[19] = inputU.d_lll.y.zz;
		
		//V_j
		results.ptr[20] = inputU.V_l.x;
		results.ptr[21] = inputU.V_l.y;
		results.ptr[22] = inputU.V_l.z;

		sym3 K_sqrt_gammaUzz = sym3_real_mul(inputU.K_ll, eig.sqrt_gammaUjj.z);

		//a^z - f d^zj_j

		real f = eig.sqrt_f * eig.sqrt_f;
		real d_z_input = sym3_dot(eig.gamma_uu, inputU.d_lll.z);
		results.ptr[23] = inputU.a_l.z - f * d_z_input;

		//gauge:
		//sqrt(f gamma^zz) K +- (a^z + 2 V^z)

		real ev0a = eig.sqrt_f * sym3_dot(eig.gamma_uu, K_sqrt_gammaUzz);
		real ev0b = eig.gamma_uu.xz * (inputU.a_l.x + 2. * inputU.V_l.x)
					+ eig.gamma_uu.yz * (inputU.a_l.y + 2. * inputU.V_l.y)
					+ eig.gamma_uu.zz * (inputU.a_l.z + 2. * inputU.V_l.z);
		results.ptr[0] = ev0a - ev0b;
		results.ptr[29] = ev0a + ev0b;

		//light:
		//sqrt(gamma^zz) K_xx +- d^z_xx
		
		real dUz_xx_input = eig.gamma_uu.xz * inputU.d_lll.x.xx + eig.gamma_uu.yz * inputU.d_lll.y.xx + eig.gamma_uu.zz * inputU.d_lll.z.xx;
		results.ptr[1] = K_sqrt_gammaUzz.xx - dUz_xx_input;
		results.ptr[24] = K_sqrt_gammaUzz.xx + dUz_xx_input;

		//light:
		//sqrt(gamma^zz) K_xy +- d^z_xy

		real dUz_xy_input = eig.gamma_uu.xz * inputU.d_lll.x.xy + eig.gamma_uu.yz * inputU.d_lll.y.xy + eig.gamma_uu.zz * inputU.d_lll.z.xy;
		results.ptr[2] = K_sqrt_gammaUzz.xy - dUz_xy_input;
		results.ptr[25] = K_sqrt_gammaUzz.xy + dUz_xy_input;

		//light:
		//sqrt(gamma^zz) K_xz +- (d^z_xz + .5 (a_x - d_xj^j) + V_x)
		
		real d_x_input = sym3_dot(eig.gamma_uu, inputU.d_lll.x);
		real dUz_xz_input = eig.gamma_uu.xz * inputU.d_lll.x.xz + eig.gamma_uu.yz * inputU.d_lll.y.xz + eig.gamma_uu.zz * inputU.d_lll.z.xz;
		real ev3b = .5 * (inputU.a_l.x - d_x_input) + inputU.V_l.x + dUz_xz_input;
		results.ptr[3] = K_sqrt_gammaUzz.xz - ev3b;
		results.ptr[26] = K_sqrt_gammaUzz.xz + ev3b;

		//light:
		//sqrt(gamma^zz) K_yy +- d^z_yy

		real dUz_yy_input = eig.gamma_uu.xz * inputU.d_lll.x.yy + eig.gamma_uu.yz * inputU.d_lll.y.yy + eig.gamma_uu.zz * inputU.d_lll.z.yy;
		results.ptr[4] = K_sqrt_gammaUzz.yy - dUz_yy_input;
		results.ptr[27] = K_sqrt_gammaUzz.yy + dUz_yy_input;
		
		//light:
		//sqrt(gamma^zz) K_yz

		real d_y_input = sym3_dot(eig.gamma_uu, inputU.d_lll.y);
		real dUz_yz_input = eig.gamma_uu.xz * inputU.d_lll.x.yz + eig.gamma_uu.yz * inputU.d_lll.y.yz + eig.gamma_uu.zz * inputU.d_lll.z.yz;
		real ev5b = .5 * (inputU.a_l.y - d_y_input) + inputU.V_l.y + dUz_yz_input;
		results.ptr[5] = K_sqrt_gammaUzz.yz - ev5b;
		results.ptr[28] = K_sqrt_gammaUzz.yz + ev5b;

	}

<? else -- eqn.noZeroRowsInFlux ?>

	real _1_sqrt_f = 1. / eig.sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f; 

	real sqrt_gammaUjj, _1_gammaUjj, a_j;
	sym3 d_lll, K_ll, gamma_uu;
	//now swap x and side on the sym3's
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		sqrt_gammaUjj = eig.sqrt_gammaUjj.s<?=side?>;
		_1_gammaUjj = 1. / eig.gamma_uu.s<?=side?><?=side?>;

		a_j = inputU.a_l.s<?=side?>;
		
		d_lll = sym3_swap<?=side?>(inputU.d_lll.v<?=side?>);
		K_ll = sym3_swap<?=side?>(inputU.K_ll);
		gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);
	}
	<? end ?>

	real K_dot_eig_gamma = sym3_dot(K_ll, gamma_uu);
	real dj_dot_eig_gamma = sym3_dot(d_lll, gamma_uu);

	results.ptr[0] = (a_j * -sqrt_gammaUjj * _1_sqrt_f + K_dot_eig_gamma) * .5 * _1_gammaUjj;

		<? for i=1,5 do ?>
	results.ptr[<?=i?>] = .5 * (-sqrt_gammaUjj * d_lll.s[<?=i?>] + K_ll.s[<?=i?>]);
		<? end ?>

	results.ptr[6] = (-a_j * _1_f + dj_dot_eig_gamma) * _1_gammaUjj;

		<? for i=1,5 do ?>
	results.ptr[<?=6+i?>] = .5 * (sqrt_gammaUjj * d_lll.s[<?=i?>] + K_ll.s[<?=i?>]);
		<? end ?>

	results.ptr[12] = (a_j * sqrt_gammaUjj * _1_sqrt_f + K_dot_eig_gamma) * .5 * _1_gammaUjj;
	
<? end -- eqn.noZeroRowsInFlux ?>
	return results;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t input,
	real3 x,
	normal_t n
) {
	cons_t resultU;
	for (int j = 0; j < numStates; ++j) {
		resultU.ptr[j] = 0;
	}

<? if not eqn.noZeroRowsInFlux then ?>
	
	if (n.side == 0) {

		//write zeros to the alpha and gammaLL terms
		resultU.alpha = 0;
		resultU.gamma_ll = sym3_zero;

		resultU.a_l.y = input.ptr[6];
		resultU.a_l.z = input.ptr[7];
		
		resultU.d_lll.y.xx = input.ptr[8];
		resultU.d_lll.y.xy = input.ptr[9];
		resultU.d_lll.y.xz = input.ptr[10];
		resultU.d_lll.y.yy = input.ptr[11];
		resultU.d_lll.y.yz = input.ptr[12];
		resultU.d_lll.y.zz = input.ptr[13];
		
		resultU.d_lll.z.xx = input.ptr[14];
		resultU.d_lll.z.xy = input.ptr[15];
		resultU.d_lll.z.xz = input.ptr[16];
		resultU.d_lll.z.yy = input.ptr[17];
		resultU.d_lll.z.yz = input.ptr[18];
		resultU.d_lll.z.zz = input.ptr[19];
		
		resultU.V_l.x = input.ptr[20];
		resultU.V_l.y = input.ptr[21];
		resultU.V_l.z = input.ptr[22];

		real sqrt_gammaUxx = eig.sqrt_gammaUjj.x;
		real _1_sqrt_gammaUxx = 1. / sqrt_gammaUxx;
		real _1_gammaUxx = _1_sqrt_gammaUxx * _1_sqrt_gammaUxx;
		real invDenom = .5 * _1_gammaUxx;

		real VUx_input = eig.gamma_uu.xx * input.ptr[20]		//V_x
						+ eig.gamma_uu.xy * input.ptr[21]	//V_y
						+ eig.gamma_uu.xz * input.ptr[22];	//V_z

		real gauge_diff = input.ptr[0] - input.ptr[29];
		real gauge_sum = input.ptr[0] + input.ptr[29];

		resultU.a_l.x = -(
				gauge_diff
				+ 2. * (
					eig.gamma_uu.xy * input.ptr[6]	//a_y
					+ eig.gamma_uu.xz * input.ptr[7]	//a_z
				)	
				+ 4. * VUx_input
			) * invDenom;
		
		real K_input_minus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[1]
				+ eig.gamma_uu.xz * input.ptr[2]
				+ eig.gamma_uu.yz * input.ptr[4]
			) + eig.gamma_uu.yy * input.ptr[3]
			+ eig.gamma_uu.zz * input.ptr[5];

		real K_input_plus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[24]
				+ eig.gamma_uu.xz * input.ptr[25]
				+ eig.gamma_uu.yz * input.ptr[27]
			) + eig.gamma_uu.yy * input.ptr[26]
			+ eig.gamma_uu.zz * input.ptr[28];

		real _1_sqrt_f = 1. / eig.sqrt_f;
		real _1_f = _1_sqrt_f * _1_sqrt_f;

		resultU.d_lll.x.xx = -(
				- K_input_minus
				+ K_input_plus

				+ 2. * eig.gamma_uu.xy * (
					eig.gamma_uu.xx * input.ptr[8]	//d_yxx
					- input.ptr[6]		//a_y
					- 2. * input.ptr[21]	//V_y
				)
				+ 2. * eig.gamma_uu.xz * (
					eig.gamma_uu.xx * input.ptr[14]	//d_zxx
					- input.ptr[7]		//a_z
					- 2. * input.ptr[22]	//V_z
				)
				//(a^x + 2 V^x) / f - 2 gamma^xx d^xj_j
				+ (
					+ gauge_diff	//-ev0b = -a^x - 2 V^x
					+ 2. * (
						+ eig.gamma_uu.xx * input.ptr[23]	//a_x - f d^xj_j
						+ eig.gamma_uu.xy * input.ptr[6]	//a_y
						+ eig.gamma_uu.xz * input.ptr[7]	//a_z
						+ 2. * VUx_input
					)
				) * _1_f
			) * invDenom * _1_gammaUxx;

		//d_yj^j
		real d_y_input = 
			eig.gamma_uu.xx * input.ptr[8]
			+ eig.gamma_uu.yy * input.ptr[11]
			+ eig.gamma_uu.zz * input.ptr[13]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[9]
				+ eig.gamma_uu.xz * input.ptr[10]
				+ eig.gamma_uu.yz * input.ptr[12]
			);
		
		resultU.d_lll.x.xy = -(
				input.ptr[1]
				+ input.ptr[6]
				- d_y_input
				+ 2. * eig.gamma_uu.xy * input.ptr[9]
				+ 2. * eig.gamma_uu.xz * input.ptr[15]
				+ 2. * input.ptr[21]
				- input.ptr[24]
			) * invDenom;
		
		//d_zj^j
		real d_z_input = 
			eig.gamma_uu.xx * input.ptr[14]
			+ eig.gamma_uu.yy * input.ptr[17]
			+ eig.gamma_uu.zz * input.ptr[19]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[15]
				+ eig.gamma_uu.xz * input.ptr[16]
				+ eig.gamma_uu.yz * input.ptr[18]
			);

		resultU.d_lll.x.xz = -(
				input.ptr[2]
				+ input.ptr[7]
				- d_z_input
				+ 2. * eig.gamma_uu.xy * input.ptr[10]
				+ 2. * eig.gamma_uu.xz * input.ptr[16]
				+ 2. * input.ptr[22]
				- input.ptr[25]
			) * invDenom;
		resultU.d_lll.x.yy = -(
				input.ptr[3]
				+ 2. * eig.gamma_uu.xy * input.ptr[11]
				+ 2. * eig.gamma_uu.xz * input.ptr[17]
				- input.ptr[26]
			) * invDenom;
		resultU.d_lll.x.yz = -(
				input.ptr[4]
				+ 2. * eig.gamma_uu.xy * input.ptr[12]
				+ 2. * eig.gamma_uu.xz * input.ptr[18]
				- input.ptr[27]
			) * invDenom;
		resultU.d_lll.x.zz = -(
				input.ptr[5]
				+ 2. * eig.gamma_uu.xy * input.ptr[13]
				+ 2. * eig.gamma_uu.xz * input.ptr[19]
				- input.ptr[28]
			) * invDenom;

		resultU.K_ll.xx = (
				- K_input_minus
				- K_input_plus	
				+ gauge_sum * _1_sqrt_f
			) * invDenom * _1_sqrt_gammaUxx;

		real tmp = .5 * _1_sqrt_gammaUxx;
		resultU.K_ll.xy = (input.ptr[1] + input.ptr[24]) * tmp;
		resultU.K_ll.xz = (input.ptr[2] + input.ptr[25]) * tmp;
		resultU.K_ll.yy = (input.ptr[3] + input.ptr[26]) * tmp;
		resultU.K_ll.yz = (input.ptr[4] + input.ptr[27]) * tmp;
		resultU.K_ll.zz = (input.ptr[5] + input.ptr[28]) * tmp;

	} else if (n.side == 1) {
	
		//write zeros to the alpha and gammaLL terms
		resultU.alpha = 0;
		resultU.gamma_ll = sym3_zero;
		
		resultU.a_l.x = input.ptr[6];
		resultU.a_l.z = input.ptr[7];
		
		resultU.d_lll.x.xx = input.ptr[8];
		resultU.d_lll.x.xy = input.ptr[9];
		resultU.d_lll.x.xz = input.ptr[10];
		resultU.d_lll.x.yy = input.ptr[11];
		resultU.d_lll.x.yz = input.ptr[12];
		resultU.d_lll.x.zz = input.ptr[13];
		
		resultU.d_lll.z.xx = input.ptr[14];
		resultU.d_lll.z.xy = input.ptr[15];
		resultU.d_lll.z.xz = input.ptr[16];
		resultU.d_lll.z.yy = input.ptr[17];
		resultU.d_lll.z.yz = input.ptr[18];
		resultU.d_lll.z.zz = input.ptr[19];
		
		resultU.V_l.x = input.ptr[20];
		resultU.V_l.y = input.ptr[21];
		resultU.V_l.z = input.ptr[22];

		real sqrt_gammaUyy = eig.sqrt_gammaUjj.y;
		real _1_sqrt_gammaUyy = 1. / sqrt_gammaUyy;
		real inv_gammaUyy = _1_sqrt_gammaUyy * _1_sqrt_gammaUyy;
		real invDenom = .5 * inv_gammaUyy;

		real VUy_input = eig.gamma_uu.xy * input.ptr[20]
						+ eig.gamma_uu.yy * input.ptr[21]
						+ eig.gamma_uu.yz * input.ptr[22];
		
		real gauge_diff = input.ptr[0] - input.ptr[29];
		real gauge_sum = input.ptr[0] + input.ptr[29];

		resultU.a_l.y = -(
				gauge_diff
				+ 2. * (
					eig.gamma_uu.xy * input.ptr[6]
					+ eig.gamma_uu.yz * input.ptr[7]
				)
				+ 4. * VUy_input
			) * invDenom;

		resultU.d_lll.y.xx = -(
				+ input.ptr[1]
				- input.ptr[24]
				+ 2. * eig.gamma_uu.xy * input.ptr[8]
				+ 2. * eig.gamma_uu.yz * input.ptr[14]
			) * invDenom;

		//d_xj^j
		real d_x_input = 
			eig.gamma_uu.xx * input.ptr[8]
			+ eig.gamma_uu.yy * input.ptr[11] 
			+ eig.gamma_uu.zz * input.ptr[13]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[9]
				+ eig.gamma_uu.xz * input.ptr[10]
				+ eig.gamma_uu.yz * input.ptr[12]
			);
		
		resultU.d_lll.y.xy = -(
				+ input.ptr[2]
				+ input.ptr[6]
				- d_x_input			
				+ 2. * eig.gamma_uu.xy * input.ptr[9]
				+ 2. * eig.gamma_uu.yz * input.ptr[15]
				+ 2. * input.ptr[20]
				- input.ptr[25]
			) * invDenom;
		
		resultU.d_lll.y.xz = -(
				+ input.ptr[3]
				+ 2. * eig.gamma_uu.xy * input.ptr[10]
				+ 2. * eig.gamma_uu.yz * input.ptr[16]
				- input.ptr[26]
			) * invDenom;
		
		real K_input_minus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[2]
				+ eig.gamma_uu.xz * input.ptr[3]
				+ eig.gamma_uu.yz * input.ptr[4]
			) + eig.gamma_uu.xx * input.ptr[1]
			+ eig.gamma_uu.zz * input.ptr[5];

		real K_input_plus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[25]
				+ eig.gamma_uu.xz * input.ptr[26]
				+ eig.gamma_uu.yz * input.ptr[27]
			) + eig.gamma_uu.xx * input.ptr[24]
			+ eig.gamma_uu.zz * input.ptr[28];

		real _1_sqrt_f = 1. / eig.sqrt_f;
		real _1_f = _1_sqrt_f  * _1_sqrt_f;

		resultU.d_lll.y.yy = -(
				- K_input_minus	
				+ K_input_plus	
				
				+ 2. * eig.gamma_uu.xy * (
					input.ptr[11] * eig.gamma_uu.yy
					- input.ptr[6]
					- 2. * input.ptr[20]
				)	
				+ 2. * eig.gamma_uu.yz * (
					input.ptr[17] * eig.gamma_uu.yy
					- input.ptr[7]
					- 2. * input.ptr[22]
				)
				+ (
					+ gauge_diff
					+ 2. * (
						+ eig.gamma_uu.yy * input.ptr[23]
						+ eig.gamma_uu.xy * input.ptr[6]
						+ eig.gamma_uu.yz * input.ptr[7]
						+ 2. * VUy_input
					)
				) * _1_f
			) * invDenom * inv_gammaUyy;

		//gamma_zj^j
		real d_z_input = eig.gamma_uu.xx * input.ptr[14]
			+ eig.gamma_uu.yy * input.ptr[17]
			+ eig.gamma_uu.zz * input.ptr[19]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[15]
				+ eig.gamma_uu.xz * input.ptr[16]
				+ eig.gamma_uu.yz * input.ptr[18]);

		resultU.d_lll.y.yz = -(
				+ input.ptr[4]
				+ input.ptr[7]
				- d_z_input	
				+ 2. * eig.gamma_uu.xy * input.ptr[12]
				+ 2. * eig.gamma_uu.yz * input.ptr[18]
				+ 2. * input.ptr[22]
				- input.ptr[27]
			) * invDenom;
		
		resultU.d_lll.y.zz = -(
				+ input.ptr[5]
				- input.ptr[28]
				+ 2. * eig.gamma_uu.xy * input.ptr[13]
				+ 2. * eig.gamma_uu.yz * input.ptr[19]
			) * invDenom;

		resultU.K_ll.yy = (
				- K_input_minus
				- K_input_plus
				+ gauge_sum * _1_sqrt_f
			) * invDenom * _1_sqrt_gammaUyy;

		real tmp = .5 * _1_sqrt_gammaUyy;
		resultU.K_ll.xx = (input.ptr[1] + input.ptr[24]) * tmp;
		resultU.K_ll.xy = (input.ptr[2] + input.ptr[25]) * tmp;		//once this gets enabled, compiling crashes
		resultU.K_ll.xz = (input.ptr[3] + input.ptr[26]) * tmp;
		resultU.K_ll.yz = (input.ptr[4] + input.ptr[27]) * tmp;
		resultU.K_ll.zz = (input.ptr[5] + input.ptr[28]) * tmp;
	
	} else if (n.side == 2) {

		//write zeros to the alpha and gammaLL terms
		resultU.alpha = 0;
		resultU.gamma_ll = sym3_zero;
		
		resultU.a_l.x = input.ptr[6];
		resultU.a_l.y = input.ptr[7];
		
		resultU.d_lll.x.xx = input.ptr[8];
		resultU.d_lll.x.xy = input.ptr[9];
		resultU.d_lll.x.xz = input.ptr[10];
		resultU.d_lll.x.yy = input.ptr[11];
		resultU.d_lll.x.yz = input.ptr[12];
		resultU.d_lll.x.zz = input.ptr[13];
		
		resultU.d_lll.y.xx = input.ptr[14];
		resultU.d_lll.y.xy = input.ptr[15];
		resultU.d_lll.y.xz = input.ptr[16];
		resultU.d_lll.y.yy = input.ptr[17];
		resultU.d_lll.y.yz = input.ptr[18];
		resultU.d_lll.y.zz = input.ptr[19];

		resultU.V_l.x = input.ptr[20];
		resultU.V_l.y = input.ptr[21];
		resultU.V_l.z = input.ptr[22];

		real sqrt_gammaUzz = eig.sqrt_gammaUjj.z;
		real _1_sqrt_gammaUzz = 1. / sqrt_gammaUzz;
		real inv_gammaUzz = _1_sqrt_gammaUzz * _1_sqrt_gammaUzz;
		real invDenom = .5 * inv_gammaUzz;

		real VUz_input = eig.gamma_uu.xz * input.ptr[20]
						+ eig.gamma_uu.yz * input.ptr[21]
						+ eig.gamma_uu.zz * input.ptr[22];

		real gauge_diff = input.ptr[0] - input.ptr[29];
		real gauge_sum = input.ptr[0] + input.ptr[29];
		
		resultU.a_l.z = -(
				gauge_diff
				+ 2. * (
					eig.gamma_uu.xz * input.ptr[6]
					+ eig.gamma_uu.yz * input.ptr[7]
				)
				+ 4. * VUz_input
			) * invDenom;
		
		resultU.d_lll.z.xx = -(
				+ input.ptr[1]
				- input.ptr[24]
				+ 2. * eig.gamma_uu.xz * input.ptr[8]
				+ 2. * eig.gamma_uu.yz * input.ptr[14]
			) * invDenom;
		
		resultU.d_lll.z.xy = -(
				+ input.ptr[2]
				- input.ptr[25]
				+ 2. * eig.gamma_uu.xz * input.ptr[9]
				+ 2. * eig.gamma_uu.yz * input.ptr[15]
			) * invDenom;

		//d_xj^j
		real d_x_input = 
			eig.gamma_uu.xx * input.ptr[8]
			+ eig.gamma_uu.yy * input.ptr[11]
			+ eig.gamma_uu.zz * input.ptr[13]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[9]
				+ eig.gamma_uu.xz * input.ptr[10]
				+ eig.gamma_uu.yz * input.ptr[12]);

		resultU.d_lll.z.xz = -(
				+ input.ptr[3]
				+ input.ptr[6]
				- d_x_input
				+ 2. * eig.gamma_uu.xz * input.ptr[10]
				+ 2. * eig.gamma_uu.yz * input.ptr[16]
				+ 2. * input.ptr[20]
				- input.ptr[26]
			) * invDenom;
		
		resultU.d_lll.z.yy = -(
				+ input.ptr[4]
				- input.ptr[27]
				+ 2. * eig.gamma_uu.xz * input.ptr[11]
				+ 2. * eig.gamma_uu.yz * input.ptr[17]
			) * invDenom;

		//d_yj^j
		real d_y_input = 
			eig.gamma_uu.xx * input.ptr[14]
			+ eig.gamma_uu.yy * input.ptr[17]
			+ eig.gamma_uu.zz * input.ptr[19]
			+ 2. * (
				eig.gamma_uu.xy * input.ptr[15]
				+ eig.gamma_uu.xz * input.ptr[16]
				+ eig.gamma_uu.yz * input.ptr[18]);

		resultU.d_lll.z.yz = -(
				+ input.ptr[5]
				+ input.ptr[7]
				- d_y_input	
				+ 2. * eig.gamma_uu.xz * input.ptr[12]
				+ 2. * eig.gamma_uu.yz * input.ptr[18]
				+ 2. * input.ptr[21]
				- input.ptr[28]
			) * invDenom;

		real K_input_minus = 
			2. * (
				eig.gamma_uu.xy * input.ptr[2]
				+ eig.gamma_uu.xz * input.ptr[3]
				+ eig.gamma_uu.yz * input.ptr[5]
			) + eig.gamma_uu.xx * input.ptr[1]
			+ eig.gamma_uu.yy * input.ptr[4];

		real K_input_plus = 
			2. * (eig.gamma_uu.xy * input.ptr[25]
				+ eig.gamma_uu.xz * input.ptr[26]
				+ eig.gamma_uu.yz * input.ptr[28]
			) + eig.gamma_uu.xx * input.ptr[24]
			+ eig.gamma_uu.yy * input.ptr[27];

		real _1_sqrt_f = 1. / eig.sqrt_f;
		real _1_f = _1_sqrt_f * _1_sqrt_f;

		resultU.d_lll.z.zz = -(
				- K_input_minus
				+ K_input_plus
				
				+ 2. * eig.gamma_uu.xz * (
					input.ptr[13] * eig.gamma_uu.zz
					- input.ptr[6]
					- 2. * input.ptr[20]
				)
				+ 2. * eig.gamma_uu.yz * (
					eig.gamma_uu.zz * input.ptr[19]
					- input.ptr[7]
					- 2. * input.ptr[21]
				)	
				+ (	
					+ gauge_diff
					+ 2. * (
						eig.gamma_uu.zz * input.ptr[23]
						+ eig.gamma_uu.xz * input.ptr[6]
						+ eig.gamma_uu.yz * input.ptr[7]
						+ 2. * VUz_input
					)
				) * _1_f
			) * invDenom * inv_gammaUzz;

		resultU.K_ll.zz = (
				+ gauge_sum * _1_sqrt_f
				- K_input_minus
				- K_input_plus
			) * invDenom * _1_sqrt_gammaUzz;
		
		real tmp = .5 * _1_sqrt_gammaUzz;
		resultU.K_ll.xx = (input.ptr[1] + input.ptr[24]) * tmp;
		resultU.K_ll.xy = (input.ptr[2] + input.ptr[25]) * tmp;
		resultU.K_ll.xz = (input.ptr[3] + input.ptr[26]) * tmp;
		resultU.K_ll.yy = (input.ptr[4] + input.ptr[27]) * tmp;
		resultU.K_ll.yz = (input.ptr[5] + input.ptr[28]) * tmp;
		
	}

<? else -- eqn.noZeroRowsInFlux ?>
	
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		//TODO swap size inside eigen_t structure
		//instead of doing it here
		sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);

		real input1_dot_gammaU = input.ptr[1] * 2. * gamma_uu.xy
			+ input.ptr[2] * 2. * gamma_uu.xz
			+ input.ptr[3] * gamma_uu.yy
			+ input.ptr[4] * 2. * gamma_uu.yz
			+ input.ptr[5] * gamma_uu.zz;
		real input7_dot_gammaU = input.ptr[7] * 2. * gamma_uu.xy
			+ input.ptr[8] * 2. * gamma_uu.xz
			+ input.ptr[9] * gamma_uu.yy
			+ input.ptr[10] * 2. * gamma_uu.yz
			+ input.ptr[11] * gamma_uu.zz;

		real _1_sqrt_f = 1. / eig.sqrt_f;
		real sqrt_gammaUjj = eig.sqrt_gammaUjj.s<?=side?>;
		real _1_sqrt_gammaUjj = 1. / sqrt_gammaUjj;
		real _1_gammaUjj = _1_sqrt_gammaUjj * _1_sqrt_gammaUjj; 

		resultU.a_l.s<?=side?> = eig.sqrt_f * sqrt_gammaUjj * (input.ptr[12] - input.ptr[0]);

		sym3 d_lll, K_ll;
		d_lll.xx = (
			(input.ptr[12] - input.ptr[0]) * _1_sqrt_f
			+ (input1_dot_gammaU - input7_dot_gammaU) * _1_gammaUjj
		) * _1_sqrt_gammaUjj + input.ptr[6];

		<? for i=1,5 do ?>
		d_lll.s[<?=i?>] = (input.ptr[<?=i+6?>] - input.ptr[<?=i?>]) * _1_sqrt_gammaUjj;
		<? end ?>

		K_ll.xx = input.ptr[0] + input.ptr[12] - (input1_dot_gammaU + input7_dot_gammaU) * _1_gammaUjj;

		<? for i=1,5 do ?>
		K_ll.s[<?=i?>] = input.ptr[<?=i?>] + input.ptr[<?=i+6?>];
		<? end ?>

		//now swap x and side on the sym3's
		resultU.d_lll.v<?=side?> = sym3_swap<?=side?>(d_lll);
		resultU.K_ll = sym3_swap<?=side?>(K_ll);
	}
	<? end ?>

<? end 	-- eqn.noZeroRowsInFlux ?>
	return resultU;
}

cons_t eigen_fluxTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x,
	normal_t n
) {
<? if not eqn.noZeroRowsInFlux then ?>

	// TODO make this a default implementation somewhere.
	// If no one has provided one then just fall back on left / wave / right transform.
	// TODO use that static function for the calc waves as well
	

	waves_t waves = eigen_leftTransform(solver, eig, inputU, x);

	<?=eqn:eigenWaveCodePrefix(n, 'eig', 'x')?>

<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode(n, 'eig', 'x', j)?>;
<? end 
?>
	return eigen_rightTransform(solver, eig, waves, x);

<? else -- noZeroRowsInFlux ?>
<? if false then 	-- by-hand ?>

	cons_t resultU;
	for (int i = 0; i < numStates; ++i) {
		resultU.ptr[i] = 0;
	}

	real f = eig.sqrt_f * eig.sqrt_f;

	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {

		//now swap x and side on the sym3's
		sym3 input_d = sym3_swap<?=side?>(inputU.d_lll.v<?=side?>);
		sym3 input_K = sym3_swap<?=side?>(inputU.K_ll);
		sym3 gamma_uu = sym3_swap<?=side?>(eig.gamma_uu);

		resultU.a_l.s<?=side?> = sym3_dot(input_K, gamma_uu) * eig.alpha * f;
		sym3 result_d = sym3_real_mul(input_K, eig.alpha);
		sym3 result_K = sym3_real_mul(input_d, eig.alpha * gamma_uu.xx);
		result_K.xx += (inputU.a_l.s<?=side?> - sym3_dot(input_d, gamma_uu)) * eig.alpha;

		//now swap x and side on the sym3's
		resultU.d_lll.v<?=side?> = sym3_swap<?=side?>(result_d);
		resultU.K_ll = sym3_swap<?=side?>(result_K);

	}
	<? end ?>

	return resultU;
<? else	-- codegen ?>
	real f = eig.sqrt_f * eig.sqrt_f;
	
	sym3 gamma_uu = eig.gamma_uu;
	
	real alpha = inputU.alpha;
	real3 V_l = inputU.V_l;
	real3 a_l = inputU.a_l;
	sym3 gamma_ll = inputU.gamma_ll;
	sym3 K_ll = inputU.K_ll;
	_3sym3 d_lll = inputU.d_lll;
	
	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		V_l = real3_swap<?=side?>(V_l);
		a_l = real3_swap<?=side?>(a_l);
		gamma_ll = sym3_swap<?=side?>(gamma_ll);
		K_ll = sym3_swap<?=side?>(K_ll);
		d_lll = _3sym3_swap<?=side?>(d_lll);
		gamma_uu = sym3_swap<?=side?>(gamma_uu);
	}
	<? end ?>

	<?=eqn.cons_t?> resultU = {.ptr={0}};

	// BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/adm_noZeroRows.html
	resultU.a_l.x = eig.alpha * f * (2. * K_ll.xy * gamma_uu.xy + 2. * K_ll.xz * gamma_uu.xz + 2. * K_ll.yz * gamma_uu.yz + K_ll.xx * gamma_uu.xx + K_ll.yy * gamma_uu.yy + K_ll.zz * gamma_uu.zz);
	resultU.d_lll.x.xx = K_ll.xx * eig.alpha;
	resultU.d_lll.x.xy = K_ll.xy * eig.alpha;
	resultU.d_lll.x.xz = K_ll.xz * eig.alpha;
	resultU.d_lll.x.yy = K_ll.yy * eig.alpha;
	resultU.d_lll.x.yz = K_ll.yz * eig.alpha;
	resultU.d_lll.x.zz = K_ll.zz * eig.alpha;
	resultU.K_ll.xx = eig.alpha * (a_l.x + d_lll.x.yy * gamma_uu.yy + 2. * d_lll.x.yz * gamma_uu.yz + d_lll.x.zz * gamma_uu.zz - d_lll.y.xx * gamma_uu.xy - 2. * d_lll.y.xy * gamma_uu.yy - 2. * d_lll.y.xz * gamma_uu.yz - d_lll.z.xx * gamma_uu.xz - 2. * d_lll.z.xy * gamma_uu.yz - 2. * d_lll.z.xz * gamma_uu.zz);
	resultU.K_ll.xy = (eig.alpha * (a_l.y - 2. * d_lll.x.yy * gamma_uu.xy - 2. * d_lll.x.yz * gamma_uu.xz - 2. * d_lll.y.yy * gamma_uu.yy - 2. * d_lll.y.yz * gamma_uu.yz - 2. * d_lll.z.yy * gamma_uu.yz - 2. * d_lll.z.yz * gamma_uu.zz)) / 2.;
	resultU.K_ll.xz = (eig.alpha * (a_l.z - 2. * d_lll.x.yz * gamma_uu.xy - 2. * d_lll.x.zz * gamma_uu.xz - 2. * d_lll.y.yz * gamma_uu.yy - 2. * d_lll.y.zz * gamma_uu.yz - 2. * d_lll.z.yz * gamma_uu.yz - 2. * d_lll.z.zz * gamma_uu.zz)) / 2.;
	resultU.K_ll.yy = eig.alpha * (d_lll.x.yy * gamma_uu.xx + d_lll.y.yy * gamma_uu.xy + d_lll.z.yy * gamma_uu.xz);
	resultU.K_ll.yz = eig.alpha * (d_lll.x.yz * gamma_uu.xx + d_lll.y.yz * gamma_uu.xy + d_lll.z.yz * gamma_uu.xz);
	resultU.K_ll.zz = eig.alpha * (d_lll.x.zz * gamma_uu.xx + d_lll.y.zz * gamma_uu.xy + d_lll.z.zz * gamma_uu.xz);
	// END CUT

	<? for side=0,solver.dim-1 do ?>
	if (n.side == <?=side?>) {
		resultU.V_l = real3_swap<?=side?>(resultU.V_l);
		resultU.a_l = real3_swap<?=side?>(resultU.a_l);
		resultU.gamma_ll = sym3_swap<?=side?>(resultU.gamma_ll);
		resultU.K_ll = sym3_swap<?=side?>(resultU.K_ll);
		resultU.d_lll = _3sym3_swap<?=side?>(resultU.d_lll);
	}
	<? end ?>

	return resultU;
<? end ?>
<? end -- noZeroRowsInFlux ?>
}

/*
this should just be 

alpha_,t + F^i^alpha_,i = -f alpha gamma^ij K_ij
*/


//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
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


#if 1	//hand-rolled
	// source terms

	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real tr_K = real3x3_trace(K_ul);							//K^k_k
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
	//		+ K K_ij 
	//		- 2 K_ik K^k_j

	//matter terms:
	//		- 8 pi S_ij 
	//		+ 4 pi  gamma_ij (S - rho)
	
	//...and the shift comes later ...

	sym3 stressConstraint_ll = (sym3){
<? for ij,xij in ipairs(symNames) do	
?>		.<?=xij?> = 
			R_ll.<?=xij?> 
			+ tr_K * U->K_ll.<?=xij?> 
			- KSq_ll.<?=xij?>
			- 8. * M_PI * S_ll.<?=xij?> 
			+ 4. * M_PI * U->gamma_ll.<?=xij?> * (S - U->rho)
		,
<? end
?>	};

	real HamiltonianConstraint = sym3_dot(stressConstraint_ll, gamma_uu);


	sym3 srcK_ll_over_alpha = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 
			- U->a_l.<?=xi?> * U->a_l.<?=xj?>
<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * U->a_l.<?=xk?>
<?	end
?>			- KSq_ll.<?=xij?>
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
<? for ij,xij in ipairs(symNames) do
?>	deriv->gamma_ll.<?=xij?> += -2. * U->alpha * U->K_ll.<?=xij?>;
<? end
?>

	//K_ij,t = shift terms + alpha srcK_ij
<? for ij,xij in ipairs(symNames) do
?>	deriv->K_ll.<?=xij?> += U->alpha * srcK_ll_over_alpha.<?=xij?>;
<? end
?>
	//V_k,t = shift terms + alpha srcV_k
<? for i,xi in ipairs(xNames) do
?>	deriv->V_l.<?=xi?> += U->alpha * srcV_l.<?=xi?>;
<? end
?>
#else	//code-generated

// TODO 

#endif

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

	deriv->beta_u = real3_add(
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

	/*
	beta^i_,t = epsilon (gamma^jk D_j D_k beta^i + 1/3 D^i D_j beta^j + R^i_j beta^j - D_j (2 alpha (K^ij - 1/3 K gamma^ij))
	= epsilon (
		gamma^jk D_j (beta^i_,k + Gamma^i_lk beta^l) 
		+ 1/3 gamma^ik D_k (beta^j_,j + Gamma^j_lj beta^l) 
		+ R^i_j beta^j 
		- (2 alpha_,j (K^ij - 1/3 K gamma^ij)
		- (2 alpha (D_j K^ij - 1/3 K_,j gamma^ij)
	)
	= epsilon (
		- 2 alpha_,j K^ij 
		+ 2/3 alpha_,j K gamma^ij
		+ gamma^jk Gamma^i_lk,j beta^l 
		+ 2 Gamma^ik_j beta^j_,k 
		- Gamma^kj_j beta^i_,k 
		+ 1/3 gamma^ik Gamma^j_lj beta^l_,k 
		+ gamma^jk beta^i_,kj 
		+ 1/3 gamma^ik beta^j_,jk 
		+ 1/3 gamma^ik Gamma^j_lj,k beta^l 
		- 2 alpha K^ij_,j 
		+ 2/3 alpha K_,j gamma^ij
		+ Gamma^i_jm Gamma^mj_l beta^l 
		- Gamma^m_jl Gamma^ij_m beta^l 
		- Gamma^mj_j Gamma^i_lm beta^l
		+ Gamma^ij_l Gamma^l_kj beta^k
		+ 1/3 gamma^ik Gamma^j_lj Gamma^l_mk beta^m
		+ R^i_j beta^j 
		- 2 alpha Gamma^i_kj K^ki 
		- 2 alpha Gamma^j_kj K^ik 
	)
	*/

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


	//add stress-energy terms and H damping:
	//TODO I see it's a term in the K_ij,t ADM formalism
	//but where is the idea of adding this to gamma_ij,t from?

<? for ij,xij in ipairs(symNames) do
?>	deriv->gamma_ll.<?=xij?> += U->alpha * (
		+ solver->gamma_ll_srcStressCoeff * stressConstraint_ll.<?=xij?>
		+ solver->gamma_ll_srcHCoeff * U->gamma_ll.<?=xij?> * HamiltonianConstraint
	);
<? end 
?>

<? for ij,xij in ipairs(symNames) do
?>	deriv->K_ll.<?=xij?> += U->alpha * (
			+ solver->K_ll_srcStressCoeff * stressConstraint_ll.<?=xij?>
			+ solver->K_ll_srcHCoeff * U->gamma_ll.<?=xij?> * HamiltonianConstraint
	);
<? end 
?>

	//converge V_i to its constrained value:


	//V_i = d_ik^k - d^k_ki <=> V_i += eta (d_ik^k - d^k_ki - V_i)
	deriv->V_l = real3_add(
		deriv->V_l,
		real3_real_mul(
			real3_sub(real3_sub(d_l, e_l), U->V_l),
			solver->V_convCoeff));

	//Kreiss-Oligar diffusion, for stability's sake?
	//isn't that just a hack to improve the stability of finite-difference, which diverges by nature?
}

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost);		
	global cons_t* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real tr_K = real3x3_trace(K_ul);							//K^k_k
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
