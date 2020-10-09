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

//used by PLM, and by the default fluxFromCons (used by hll, or roe when roeUseFluxFromCons is set)
eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	normal_t n
) {
	eigen_t eig;
	eig.alpha = U.alpha;
	eig.alpha_sqrt_f = sqrt(calc_f_alphaSq(U.alpha));
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

	real gammaUjj = 0./0.;
	if (n.side == 0) {
		gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
	} else if (n.side == 1) {
		gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
	} else if (n.side == 2) {
		gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
	}
	real sqrt_gammaUjj = sqrt(gammaUjj);
	real lambdaLight = U->alpha * sqrt_gammaUjj;
	
	real f_alphaSq = calc_f_alphaSq(U->alpha);
	real lambdaGauge = sqrt(f_alphaSq) * sqrt_gammaUjj;

	real lambdaMax = max(lambdaGauge, lambdaLight);
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
	eig.alpha_sqrt_f = sqrt(calc_f_alphaSq(alpha));
	eig.gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);
	eig.sqrt_gammaUjj.x = sqrt(eig.gamma_uu.xx);
	eig.sqrt_gammaUjj.y = sqrt(eig.gamma_uu.yy);
	eig.sqrt_gammaUjj.z = sqrt(eig.gamma_uu.zz);
	
	<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = real3_real_mul(real3_add(UL.beta_u, UR.beta_u), .5);
	<? end ?>

	return eig;
}

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x,
	normal_t n
) {
	waves_t results = (waves_t){.ptr={0. / 0.}};
	//for (int i = 0; i < numWaves; ++i) {
	//	results.ptr[i] = 0;
	//}

#if 0	//don't enable this.  it's made for > waves than I'm using, so it will cause buffer corruption

	real3 a_l = real3_swap(inputU.a_l, n.side);							//0-2
	_3sym3 d_lll = _3sym3_swap(inputU.d_lll, n.side);					//3-20 ... .x = 3-8, .y = 9-14, .z = 15-20
	sym3 K_ll = sym3_swap(inputU.K_ll, n.side);							//21-26
	real Theta = inputU.Theta;												//27
	real3 Z_l = real3_swap(inputU.Z_l, n.side);							//28-30
	
	sym3 gamma_ll = sym3_swap(eig.gamma_ll, n.side);
	sym3 gamma_uu = sym3_swap(eig.gamma_uu, n.side);

	real sqrt_gammaUUxx = sqrt(gamma_uu.xx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;
	real gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;
	
	real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
	real f = sqrt_f * sqrt_f;
	// from 2004 Bona et al, "A symmetry breaking..." eqn A.20
	// mind you the 'm' in that form is for alpha_,t = -alpha^2 (f K - m Theta)
	// while in more modern Z4 papers it is alpha_,t = -alpha^2 f (K - m Theta)
	// the difference means that this eigendecomposition is probably incorrect.
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

#endif
	return results;
}

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t input,
	real3 x,
	normal_t n
) {
	cons_t resultU = (cons_t){.ptr={0. / 0.}};
	//for (int j = 0; j < numStates; ++j) {
	//	resultU.ptr[j] = 0;
	//}
	
#if 0	
	sym3 gamma_ll = sym3_swap(eig.gamma_ll, n.side);
	sym3 gamma_uu = sym3_swap(eig.gamma_uu, n.side);
	
	real sqrt_gammaUUxx = sqrt(gamma_uu.xx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;
	real gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;

	real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
	real f = sqrt_f * sqrt_f;
	real fSq = f * f;
	// from 2004 Bona et al, "A symmetry breaking..." eqn A.20
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

	resultU.a_l = real3_swap(resultU.a_l, n.side);			//0-2
	resultU.d_lll = _3sym3_swap(resultU.d_lll, n.side);		//3-20
	resultU.K_ll = sym3_swap(resultU.K_ll, n.side);			//21-26
	resultU.Theta = resultU.Theta;							//27
	resultU.Z_l = real3_swap(resultU.Z_l, n.side);			//28-30
#endif
	return resultU;
}

//so long as roeUseFluxFromCons isn't set for the roe solver, 
// and fluxFromCons is provided/unused,
// eigen_fluxTransform isn't needed.
// but some solvers do use a boilerplate right(lambda(left(U)))
//however if you want to use the HLL solver then fluxFromCons is needed
//...however fluxFromCons is not provided by this eqn.
cons_t eigen_fluxTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t inputU,
	real3 x,
	normal_t n
) {
#if 0	
	//default
	waves_t waves = eigen_leftTransform(solver, eig, inputU, x, n);
	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'x')?>
<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode('n', 'eig', 'x', j)?>;
<? end 
?>	return eigen_rightTransform(solver, eig, waves, x, n);
#else
	cons_t F = {.ptr={0. / 0.}};
	return F;
#endif
}

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	global cons_t* U = UBuf + index;
	global cons_t* deriv = derivBuf + index;

	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real3 S_l = real3_zero;
	sym3 S_ll = sym3_zero;
	real S = 0.;
	real rho = 0.;

#if 0	//hand-rolled
	// source terms
	
	real3x3 K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real trK = real3x3_trace(K_ul);								//K^k_k
//	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);		//KSq_ij = K_ik K^k_j

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3x3 d_llu = _3sym3_sym3_mul(U->d_lll, gamma_uu);
	
	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);
	
	//e_l = e_i = d^j_ji
	real3 e_l = _3sym3_tr12(d_ull);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?> - d_ull.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	//d_l = d_i = d_ij^j
	real3 d_l = real3x3x3_tr23(d_llu);
	
	real3 d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 Z_u = sym3_real3_mul(gamma_uu, U->Z_l);

	//d_luu = d_i^jk = gamma^jl d_il^k
	_3sym3 d_luu = (_3sym3){
<? for i,xi in ipairs(xNames) do		
?>		.<?=xi?> = sym3_real3x3_to_sym3_mul(gamma_uu, d_llu.<?=xi?>),
<? end
?>	};

	//alpha_,t = shift terms - alpha^2 f (K - m Theta)
	real f_alphaSq = calc_f_alphaSq(U->alpha);
	deriv->alpha += -f_alphaSq * (trK - solver->m * U->Theta);
	
	//gamma_ij,t = shift terms - 2 alpha K_ij
	deriv->gamma_ll = sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));

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
		- 2. * e_l.<?=xk?> * (d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?>)
		+ (d_l.<?=xk?> + U->a_l.<?=xk?> - 2. * U->Z_l.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		- 2. * K_ul.<?=xk?>.<?=xi?> * U->K_ll.<?=sym(k,j)?>
<?		for l,xl in ipairs(xNames) do
?>		+ 2. * (d_llu.<?=xi?>.<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,j)?> + d_llu.<?=xj?>.<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,i)?>)
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
#else	//code-generated

	real alpha = U->alpha;
	real Theta = U->Theta;
	real3 Z_l = U->Z_l;
	real3 a_l = U->a_l;
	sym3 gamma_ll = U->gamma_ll;
	sym3 K_ll = U->K_ll;
	_3sym3 d_lll = U->d_lll;

	//for 1 + log(alpha) slicing:
	//f = 2 / alpha
	//f alpha = 2
	//df/dalpha = -2 / alpha^2
	//df/dalpha alpha = -2 / alpha
	//df/dalpha alpha^2 = -2
	real f_alpha = calc_f_alpha(alpha);	
	real f_alphaSq = calc_f_alphaSq(alpha);	
	real alphaSq_dalpha_f = calc_alphaSq_dalpha_f(alpha);

	//TODO:
	//alpha_,t = f alpha^2 ( ... )
	//a_k,t = f alpha ( ... ) + ...
	//
	// BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/adm_noZeroRows.html
	real tmp1 = K_ll.xy * gamma_uu.xy;
	real tmp2 = K_ll.xz * gamma_uu.xz;
	real tmp3 = K_ll.yz * gamma_uu.yz;
	real tmp4 = K_ll.zz * gamma_uu.zz;
	real tmp8 = K_ll.yy * gamma_uu.yy;
	real tmp37 = gamma_uu.xx * gamma_uu.xx;
	real tmp43 = gamma_uu.xx * gamma_uu.xy;
	real tmp49 = gamma_uu.xx * gamma_uu.xz;
	real tmp55 = gamma_uu.xy * gamma_uu.xy;
	real tmp61 = gamma_uu.xy * gamma_uu.xz;
	real tmp67 = gamma_uu.xz * gamma_uu.xz;
	real tmp86 = gamma_uu.xx * gamma_uu.yy;
	real tmp98 = gamma_uu.xx * gamma_uu.yz;
	real tmp110 = gamma_uu.xy * gamma_uu.yy;
	real tmp116 = gamma_uu.xy * gamma_uu.yz;
	real tmp122 = gamma_uu.xz * gamma_uu.yy;
	real tmp128 = gamma_uu.xz * gamma_uu.yz;
	real tmp159 = gamma_uu.xx * gamma_uu.zz;
	real tmp177 = gamma_uu.xy * gamma_uu.zz;
	real tmp189 = gamma_uu.xz * gamma_uu.zz;
	real tmp218 = gamma_uu.yy * gamma_uu.yy;
	real tmp224 = gamma_uu.yy * gamma_uu.yz;
	real tmp230 = gamma_uu.yz * gamma_uu.yz;
	real tmp279 = gamma_uu.yy * gamma_uu.zz;
	real tmp291 = gamma_uu.yz * gamma_uu.zz;
	real tmp332 = gamma_uu.zz * gamma_uu.zz;
	real tmp1049 = a_l.x * alpha;
	real tmp1061 = a_l.y * alpha;
	real tmp1073 = a_l.z * alpha;
	real tmp1109 = d_lll.x.xx * gamma_uu.xx;
	real tmp1111 = d_lll.x.xy * gamma_uu.xy;
	real tmp1113 = d_lll.x.xz * gamma_uu.xz;
	real tmp1115 = d_lll.y.xx * gamma_uu.xy;
	real tmp1119 = d_lll.z.xx * gamma_uu.xz;
	real tmp1123 = d_lll.x.xx * gamma_uu.xy;
	real tmp1125 = d_lll.x.xy * gamma_uu.yy;
	real tmp1127 = d_lll.x.xz * gamma_uu.yz;
	real tmp1129 = d_lll.y.xx * gamma_uu.yy;
	real tmp1133 = d_lll.z.xx * gamma_uu.yz;
	real tmp1137 = d_lll.x.xx * gamma_uu.xz;
	real tmp1139 = d_lll.x.xy * gamma_uu.yz;
	real tmp1141 = d_lll.x.xz * gamma_uu.zz;
	real tmp1143 = d_lll.y.xx * gamma_uu.yz;
	real tmp1147 = d_lll.z.xx * gamma_uu.zz;
	real tmp1191 = d_lll.x.yy * tmp86;
	real tmp1195 = d_lll.x.yz * tmp98;
	real tmp1200 = d_lll.x.zz * tmp159;
	real tmp1204 = d_lll.y.xy * tmp86;
	real tmp1207 = d_lll.y.xz * tmp98;
	real tmp1210 = d_lll.y.yy * tmp110;
	real tmp1212 = d_lll.y.yz * tmp122;
	real tmp1215 = d_lll.y.zz * tmp177;
	real tmp1219 = d_lll.y.zz * tmp128;
	real tmp1222 = d_lll.z.xy * tmp98;
	real tmp1225 = d_lll.z.xz * tmp159;
	real tmp1228 = d_lll.z.yy * tmp116;
	real tmp1231 = d_lll.z.yy * tmp122;
	real tmp1235 = d_lll.z.yz * tmp177;
	real tmp1238 = d_lll.z.zz * tmp189;
	real tmp1240 = d_lll.x.yy * tmp110;
	real tmp1245 = d_lll.x.yz * tmp116;
	real tmp1250 = d_lll.x.zz * tmp177;
	real tmp1255 = d_lll.y.xy * tmp110;
	real tmp1258 = d_lll.y.xz * tmp116;
	real tmp1261 = d_lll.y.yy * tmp218;
	real tmp1264 = d_lll.y.yz * tmp224;
	real tmp1267 = d_lll.y.zz * tmp279;
	real tmp1272 = d_lll.y.zz * tmp230;
	real tmp1275 = d_lll.z.xy * tmp116;
	real tmp1278 = d_lll.z.xz * tmp177;
	real tmp1281 = d_lll.z.yy * tmp224;
	real tmp1284 = d_lll.z.yz * tmp279;
	real tmp1287 = d_lll.z.zz * tmp291;
	real tmp1290 = d_lll.x.yy * tmp122;
	real tmp1295 = d_lll.x.yz * tmp128;
	real tmp1300 = d_lll.x.zz * tmp189;
	real tmp1305 = d_lll.y.xy * tmp122;
	real tmp1308 = d_lll.y.xz * tmp128;
	real tmp1311 = d_lll.y.yy * tmp224;
	real tmp1314 = d_lll.y.yz * tmp279;
	real tmp1317 = d_lll.y.zz * tmp291;
	real tmp1320 = d_lll.z.xy * tmp128;
	real tmp1323 = d_lll.z.xz * tmp189;
	real tmp1326 = d_lll.z.yy * tmp279;
	real tmp1331 = d_lll.z.yy * tmp230;
	real tmp1334 = d_lll.z.yz * tmp291;
	real tmp1337 = d_lll.z.zz * tmp332;
	real tmp1340 = d_lll.x.yz * tmp224;
	real tmp1345 = d_lll.x.zz * tmp230;
	real tmp1359 = d_lll.z.xx * tmp122;
	real tmp1360 = d_lll.x.yy * d_lll.x.yy;
	real tmp1363 = tmp1360 * tmp218;
	real tmp1365 = d_lll.x.zz * tmp291;
	real tmp1375 = d_lll.z.xx * tmp177;
	real tmp1379 = d_lll.x.yz * d_lll.x.yz;
	real tmp1381 = tmp1379 * tmp230;
	real tmp1395 = d_lll.x.zz * d_lll.x.zz;
	real tmp1398 = tmp1395 * tmp332;
	real tmp1405 = d_lll.y.xz * tmp122;
	real tmp1426 = d_lll.z.xx * tmp98;
	real tmp1436 = d_lll.z.xy * tmp122;
	real tmp1465 = d_lll.z.xx * tmp128;
	real tmp1470 = d_lll.z.xy * tmp279;
	real tmp1473 = d_lll.z.xy * tmp230;
	real tmp1477 = d_lll.y.xz * d_lll.y.xz;
	real tmp1479 = tmp1477 * tmp230;
	real tmp1494 = d_lll.z.xy * tmp177;
	real tmp1519 = d_lll.z.xy * d_lll.z.xy;
	real tmp1521 = tmp1519 * tmp230;
	real tmp1523 = d_lll.y.xx * d_lll.y.xx;
	real tmp1527 = d_lll.z.xx * d_lll.z.xx;
	real tmp1531 = K_ll.xy * K_ll.xy;
	real tmp1532 = gamma_uu.yy * tmp1531;
	real tmp1534 = gamma_uu.zz * tmp1379;
	real tmp1535 = gamma_uu.yy * tmp1534;
	real tmp1539 = gamma_uu.zz * tmp1477;
	real tmp1540 = gamma_uu.yy * tmp1539;
	real tmp1544 = gamma_uu.zz * tmp1519;
	real tmp1545 = gamma_uu.yy * tmp1544;
	real tmp1548 = K_ll.xz * K_ll.xz;
	real tmp1549 = gamma_uu.zz * tmp1548;
	real tmp1679 = K_ll.yy * gamma_uu.xy;
	real tmp1681 = K_ll.yz * gamma_uu.xz;
	real tmp1699 = d_lll.x.yy * gamma_uu.xy;
	real tmp1701 = d_lll.x.yz * gamma_uu.xz;
	real tmp1703 = d_lll.y.xx * gamma_uu.xx;
	real tmp1705 = d_lll.y.xz * gamma_uu.xz;
	real tmp1707 = d_lll.z.xy * gamma_uu.xz;
	real tmp1711 = d_lll.x.yy * gamma_uu.yy;
	real tmp1713 = d_lll.x.yz * gamma_uu.yz;
	real tmp1717 = d_lll.y.xz * gamma_uu.yz;
	real tmp1719 = d_lll.z.xy * gamma_uu.yz;
	real tmp1723 = d_lll.x.yy * gamma_uu.yz;
	real tmp1725 = d_lll.x.yz * gamma_uu.zz;
	real tmp1727 = d_lll.y.xx * gamma_uu.xz;
	real tmp1729 = d_lll.y.xz * gamma_uu.zz;
	real tmp1731 = d_lll.z.xy * gamma_uu.zz;
	real tmp1775 = d_lll.x.yy * tmp43;
	real tmp1777 = d_lll.x.yz * tmp49;
	real tmp1779 = d_lll.y.xy * tmp43;
	real tmp1784 = d_lll.y.xz * tmp49;
	real tmp1788 = d_lll.y.yy * tmp55;
	real tmp1792 = d_lll.y.yz * tmp61;
	real tmp1797 = d_lll.y.zz * tmp67;
	real tmp1801 = d_lll.z.xy * tmp49;
	real tmp1805 = d_lll.x.yy * tmp55;
	real tmp1808 = d_lll.x.yz * tmp61;
	real tmp1811 = d_lll.y.xy * tmp55;
	real tmp1816 = d_lll.y.xz * tmp61;
	real tmp1826 = d_lll.y.yz * tmp116;
	real tmp1836 = d_lll.z.xy * tmp61;
	real tmp1855 = d_lll.y.xy * tmp61;
	real tmp1860 = d_lll.y.xz * tmp67;
	real tmp1865 = d_lll.y.yy * tmp116;
	real tmp1870 = d_lll.y.yz * tmp177;
	real tmp1875 = d_lll.y.zz * tmp189;
	real tmp1880 = d_lll.z.xy * tmp67;
	real tmp1885 = d_lll.z.yy * tmp177;
	real tmp1888 = d_lll.z.yy * tmp128;
	real tmp1896 = d_lll.x.yz * tmp122;
	real tmp1902 = d_lll.x.zz * tmp128;
	real tmp1905 = d_lll.y.xx * tmp55;
	real tmp1936 = d_lll.y.xx * tmp98;
	real tmp1940 = d_lll.y.xx * tmp61;
	real tmp1943 = d_lll.y.xy * tmp116;
	real tmp1952 = d_lll.y.yz * tmp230;
	real tmp1961 = d_lll.z.xx * tmp159;
	real tmp1990 = d_lll.y.xx * tmp67;
	real tmp1992 = d_lll.y.xy * tmp128;
	real tmp1997 = d_lll.y.xz * tmp189;
	real tmp2001 = d_lll.y.yy * tmp230;
	real tmp2005 = d_lll.y.yz * tmp291;
	real tmp2010 = d_lll.y.zz * tmp332;
	real tmp2038 = d_lll.z.xx * tmp49;
	real tmp2060 = d_lll.z.xx * tmp61;
	real tmp2079 = d_lll.z.xx * tmp67;
	real tmp2100 = d_lll.z.xy * tmp224;
	real tmp2121 = d_lll.z.xy * tmp159;
	real tmp2153 = gamma_uu.yy * tmp1360;
	real tmp2161 = gamma_uu.yz * tmp1519;
	real tmp2162 = gamma_uu.xz * tmp2161;
	real tmp2331 = K_ll.xz * gamma_uu.xx;
	real tmp2332 = K_ll.yz * gamma_uu.xy;
	real tmp2334 = K_ll.zz * gamma_uu.xz;
	real tmp2336 = K_ll.yz * gamma_uu.yy;
	real tmp2338 = K_ll.zz * gamma_uu.yz;
	real tmp2352 = d_lll.x.yz * gamma_uu.xy;
	real tmp2354 = d_lll.x.zz * gamma_uu.xz;
	real tmp2356 = d_lll.y.xz * gamma_uu.xy;
	real tmp2360 = d_lll.z.xx * gamma_uu.xx;
	real tmp2362 = d_lll.z.xy * gamma_uu.xy;
	real tmp2364 = d_lll.x.yz * gamma_uu.yy;
	real tmp2366 = d_lll.x.zz * gamma_uu.yz;
	real tmp2368 = d_lll.y.xz * gamma_uu.yy;
	real tmp2372 = d_lll.z.xx * gamma_uu.xy;
	real tmp2374 = d_lll.z.xy * gamma_uu.yy;
	real tmp2378 = d_lll.x.zz * gamma_uu.zz;
	real tmp2428 = d_lll.x.yz * tmp43;
	real tmp2430 = d_lll.x.zz * tmp49;
	real tmp2432 = d_lll.y.xz * tmp43;
	real tmp2436 = d_lll.z.xy * tmp43;
	real tmp2440 = d_lll.z.xz * tmp49;
	real tmp2445 = d_lll.z.yy * tmp55;
	real tmp2449 = d_lll.z.yz * tmp61;
	real tmp2454 = d_lll.z.zz * tmp67;
	real tmp2461 = d_lll.x.zz * tmp61;
	real tmp2464 = d_lll.y.xz * tmp55;
	real tmp2469 = d_lll.y.zz * tmp116;
	real tmp2474 = d_lll.y.zz * tmp122;
	real tmp2477 = d_lll.z.xy * tmp55;
	real tmp2482 = d_lll.z.xz * tmp61;
	real tmp2487 = d_lll.z.yy * tmp110;
	real tmp2492 = d_lll.z.yz * tmp122;
	real tmp2497 = d_lll.z.zz * tmp128;
	real tmp2505 = d_lll.x.zz * tmp67;
	real tmp2526 = d_lll.z.xz * tmp67;
	real tmp2536 = d_lll.z.yz * tmp128;
	real tmp2555 = d_lll.y.xz * tmp110;
	real tmp2559 = d_lll.z.xx * tmp86;
	real tmp2563 = d_lll.z.xx * tmp55;
	real tmp2565 = d_lll.z.xy * tmp110;
	real tmp2569 = d_lll.z.xz * tmp116;
	real tmp2574 = d_lll.z.yy * tmp218;
	real tmp2578 = d_lll.z.yz * tmp224;
	real tmp2583 = d_lll.z.zz * tmp230;
	real tmp2592 = d_lll.y.xx * tmp86;
	real tmp2624 = d_lll.z.xz * tmp128;
	real tmp2633 = d_lll.z.yz * tmp230;
	real tmp2647 = d_lll.y.xz * tmp177;
	real tmp2671 = d_lll.y.xz * tmp86;
	real tmp2683 = d_lll.z.xx * tmp43;
	real tmp2685 = d_lll.z.xy * tmp86;
	real tmp2692 = d_lll.z.xz * tmp98;
	real tmp2700 = d_lll.z.yz * tmp116;
	real tmp2806 = gamma_uu.yz * tmp1477;
	real tmp2811 = gamma_uu.yy * tmp1379;
	real tmp2814 = gamma_uu.yy * tmp1519;
	real tmp2817 = gamma_uu.zz * tmp1395;
	real tmp3012 = d_lll.x.yy * gamma_uu.xx;
	real tmp3014 = d_lll.y.xy * gamma_uu.xx;
	real tmp3018 = d_lll.y.yy * gamma_uu.xy;
	real tmp3022 = d_lll.y.yz * gamma_uu.xz;
	real tmp3026 = d_lll.z.yy * gamma_uu.xz;
	real tmp3030 = d_lll.y.xy * gamma_uu.xy;
	real tmp3034 = d_lll.y.yy * gamma_uu.yy;
	real tmp3038 = d_lll.y.yz * gamma_uu.yz;
	real tmp3042 = d_lll.z.yy * gamma_uu.yz;
	real tmp3044 = d_lll.x.yy * gamma_uu.xz;
	real tmp3046 = d_lll.y.xy * gamma_uu.xz;
	real tmp3050 = d_lll.y.yy * gamma_uu.yz;
	real tmp3054 = d_lll.y.yz * gamma_uu.zz;
	real tmp3058 = d_lll.z.yy * gamma_uu.zz;
	real tmp3122 = d_lll.y.yy * tmp86;
	real tmp3127 = d_lll.y.yz * tmp98;
	real tmp3132 = d_lll.z.yy * tmp98;
	real tmp3153 = d_lll.z.yy * tmp159;
	real tmp3166 = d_lll.y.xx * tmp43;
	real tmp3204 = d_lll.y.yy * tmp122;
	real tmp3209 = d_lll.y.yz * tmp128;
	real tmp3231 = tmp1379 * tmp67;
	real tmp3255 = d_lll.z.yy * tmp189;
	real tmp3275 = d_lll.z.yy * tmp61;
	real tmp3324 = tmp1477 * tmp67;
	real tmp3386 = d_lll.y.zz * d_lll.y.zz;
	real tmp3390 = tmp3386 * tmp332;
	real tmp3400 = tmp1519 * tmp67;
	real tmp3429 = d_lll.z.yy * d_lll.z.yy;
	real tmp3431 = K_ll.yz * K_ll.yz;
	real tmp3432 = gamma_uu.zz * tmp3431;
	real tmp3479 = tmp1523 * tmp37;
	real tmp3611 = d_lll.x.yz * gamma_uu.xx;
	real tmp3613 = d_lll.y.xz * gamma_uu.xx;
	real tmp3617 = d_lll.y.zz * gamma_uu.xz;
	real tmp3621 = d_lll.z.xy * gamma_uu.xx;
	real tmp3625 = d_lll.z.yy * gamma_uu.xy;
	real tmp3635 = d_lll.y.zz * gamma_uu.yz;
	real tmp3643 = d_lll.z.yy * gamma_uu.yy;
	real tmp3653 = d_lll.y.zz * gamma_uu.zz;
	real tmp3715 = d_lll.y.zz * tmp98;
	real tmp3725 = d_lll.z.yy * tmp86;
	real tmp3755 = d_lll.x.zz * tmp98;
	real tmp3844 = d_lll.y.xy * tmp98;
	real tmp3852 = d_lll.y.xz * tmp159;
	real tmp3890 = d_lll.z.xx * tmp37;
	real tmp3968 = d_lll.y.zz * tmp224;
	real tmp4086 = gamma_uu.zz * tmp3386;
	real tmp4256 = d_lll.x.zz * gamma_uu.xx;
	real tmp4258 = d_lll.y.zz * gamma_uu.xy;
	real tmp4260 = d_lll.z.xz * gamma_uu.xx;
	real tmp4264 = d_lll.z.yz * gamma_uu.xy;
	real tmp4268 = d_lll.z.zz * gamma_uu.xz;
	real tmp4272 = d_lll.x.zz * gamma_uu.xy;
	real tmp4274 = d_lll.y.zz * gamma_uu.yy;
	real tmp4276 = d_lll.z.xz * gamma_uu.xy;
	real tmp4280 = d_lll.z.yz * gamma_uu.yy;
	real tmp4284 = d_lll.z.zz * gamma_uu.yz;
	real tmp4292 = d_lll.z.xz * gamma_uu.xz;
	real tmp4296 = d_lll.z.yz * gamma_uu.yz;
	real tmp4300 = d_lll.z.zz * gamma_uu.zz;
	real tmp4361 = d_lll.y.zz * tmp86;
	real tmp4364 = d_lll.z.xz * tmp43;
	real tmp4369 = d_lll.z.yz * tmp86;
	real tmp4374 = d_lll.z.zz * tmp98;
	real tmp4390 = d_lll.z.yz * tmp98;
	real tmp4395 = d_lll.z.zz * tmp159;
	real tmp4407 = d_lll.y.zz * tmp110;
	real tmp4409 = d_lll.z.xz * tmp86;
	real tmp4417 = d_lll.z.yz * tmp110;
	real tmp4427 = d_lll.z.zz * tmp122;
	real tmp4458 = d_lll.z.zz * tmp177;
	real tmp4464 = tmp1379 * tmp55;
	real tmp4565 = tmp1477 * tmp55;
	real tmp4598 = d_lll.z.zz * tmp279;
	real tmp4644 = tmp1519 * tmp55;
	real tmp4657 = tmp3429 * tmp218;
	real tmp4705 = tmp1527 * tmp37;
	real tmp4920 = d_lll.x.xx * tmp37;
	real tmp4924 = d_lll.x.xy * tmp43;
	real tmp4929 = d_lll.x.xz * tmp49;
	real tmp5025 = d_lll.x.xx * tmp43;
	real tmp5029 = d_lll.x.xy * tmp86;
	real tmp5034 = d_lll.x.xz * tmp98;
	real tmp5130 = d_lll.x.xx * tmp49;
	real tmp5134 = d_lll.x.xy * tmp98;
	real tmp5139 = d_lll.x.xz * tmp159;
	real tmp5144 = d_lll.x.yy * tmp116;
	real tmp5151 = d_lll.x.yz * tmp177;
	real tmp5232 = gamma_uu.xx * tmp55;
	real tmp5237 = gamma_uu.yy * tmp37;
	real tmp5240 = gamma_uu.xx * tmp61;
	real tmp5246 = gamma_uu.yz * tmp37;
	real tmp5250 = gamma_uu.xx * tmp67;
	real tmp5255 = gamma_uu.zz * tmp37;
	real tmp5278 = gamma_uu.xx * tmp110;
	real tmp5283 = gamma_uu.xy * tmp55;
	real tmp5286 = gamma_uu.xx * tmp122;
	real tmp5292 = gamma_uu.xz * tmp55;
	real tmp5296 = gamma_uu.xx * tmp177;
	real tmp5299 = gamma_uu.xx * tmp128;
	real tmp5305 = gamma_uu.xy * tmp67;
	real tmp5328 = gamma_uu.xx * tmp116;
	real tmp5350 = gamma_uu.xx * tmp189;
	real tmp5355 = gamma_uu.xz * tmp67;
	real tmp5389 = d_lll.y.xy * tmp5278;
	real tmp5395 = d_lll.y.xy * tmp5283;
	real tmp5399 = d_lll.y.xz * tmp5328;
	real tmp5405 = d_lll.y.xz * tmp5292;
	real tmp5408 = gamma_uu.xx * tmp218;
	real tmp5409 = d_lll.y.yy * tmp5408;
	real tmp5414 = gamma_uu.yy * tmp55;
	real tmp5415 = d_lll.y.yy * tmp5414;
	real tmp5418 = gamma_uu.xx * tmp224;
	real tmp5419 = d_lll.y.yz * tmp5418;
	real tmp5424 = gamma_uu.yz * tmp55;
	real tmp5425 = d_lll.y.yz * tmp5424;
	real tmp5428 = gamma_uu.xx * tmp279;
	real tmp5429 = d_lll.y.zz * tmp5428;
	real tmp5432 = gamma_uu.xx * tmp230;
	real tmp5433 = d_lll.y.zz * tmp5432;
	real tmp5438 = gamma_uu.xy * tmp128;
	real tmp5439 = d_lll.y.zz * tmp5438;
	real tmp5442 = gamma_uu.yy * tmp67;
	real tmp5443 = d_lll.y.zz * tmp5442;
	real tmp5449 = d_lll.z.xy * tmp5328;
	real tmp5455 = d_lll.z.xy * tmp5292;
	real tmp5459 = d_lll.z.xz * tmp5296;
	real tmp5465 = d_lll.z.xz * tmp5305;
	real tmp5469 = d_lll.z.yy * tmp5418;
	real tmp5474 = gamma_uu.xy * tmp122;
	real tmp5479 = d_lll.z.yy * tmp5424;
	real tmp5485 = d_lll.z.yz * tmp5428;
	real tmp5491 = d_lll.z.yz * tmp5442;
	real tmp5494 = gamma_uu.xx * tmp291;
	real tmp5495 = d_lll.z.zz * tmp5494;
	real tmp5500 = gamma_uu.yz * tmp67;
	real tmp5501 = d_lll.z.zz * tmp5500;
	real tmp5570 = gamma_uu.zz * tmp55;
	real tmp5580 = gamma_uu.xy * tmp189;
	real tmp5591 = d_lll.z.xy * tmp5299;
	real tmp5597 = d_lll.z.xy * tmp5305;
	real tmp5601 = d_lll.z.xz * tmp5350;
	real tmp5607 = d_lll.z.xz * tmp5355;
	real tmp5611 = d_lll.z.yy * tmp5428;
	real tmp5615 = d_lll.z.yy * tmp5432;
	real tmp5621 = d_lll.z.yy * tmp5438;
	real tmp5625 = d_lll.z.yy * tmp5570;
	real tmp5631 = d_lll.z.yz * tmp5494;
	real tmp5637 = d_lll.z.yz * tmp5500;
	real tmp5640 = gamma_uu.xx * tmp332;
	real tmp5641 = d_lll.z.zz * tmp5640;
	real tmp5646 = gamma_uu.zz * tmp67;
	real tmp5647 = d_lll.z.zz * tmp5646;
	real tmp5698 = d_lll.y.xz * tmp5474;
	real tmp5702 = d_lll.y.xz * tmp5424;
	real tmp5707 = gamma_uu.xy * tmp279;
	real tmp5708 = d_lll.y.zz * tmp5707;
	real tmp5710 = gamma_uu.xy * tmp230;
	real tmp5711 = d_lll.y.zz * tmp5710;
	real tmp5716 = d_lll.z.xx * tmp5286;
	real tmp5719 = d_lll.z.xx * tmp5292;
	real tmp5724 = d_lll.z.xy * tmp5474;
	real tmp5728 = d_lll.z.xy * tmp5424;
	real tmp5734 = d_lll.z.xz * tmp5428;
	real tmp5738 = d_lll.z.xz * tmp5432;
	real tmp5753 = gamma_uu.xy * tmp224;
	real tmp5759 = gamma_uu.xz * tmp218;
	real tmp5764 = d_lll.z.yz * tmp5707;
	real tmp5770 = d_lll.z.yz * tmp5710;
	real tmp5775 = gamma_uu.xz * tmp224;
	real tmp5779 = gamma_uu.xy * tmp291;
	real tmp5785 = gamma_uu.xz * tmp279;
	real tmp5786 = d_lll.z.zz * tmp5785;
	real tmp5788 = gamma_uu.xz * tmp230;
	real tmp5789 = d_lll.z.zz * tmp5788;
	real tmp5852 = d_lll.y.yy * tmp5753;
	real tmp5856 = d_lll.y.yy * tmp5759;
	real tmp5862 = d_lll.y.yz * tmp5710;
	real tmp5866 = d_lll.y.yz * tmp5775;
	real tmp5872 = d_lll.y.zz * tmp5779;
	real tmp5876 = d_lll.y.zz * tmp5788;
	real tmp5882 = d_lll.z.xx * tmp5299;
	real tmp5886 = d_lll.z.xx * tmp5305;
	real tmp5892 = d_lll.z.xy * tmp5428;
	real tmp5898 = d_lll.z.xy * tmp5432;
	real tmp5902 = d_lll.z.xy * tmp5438;
	real tmp5908 = d_lll.z.xy * tmp5442;
	real tmp5912 = d_lll.z.xy * tmp5570;
	real tmp5916 = d_lll.z.xz * tmp5580;
	real tmp5922 = d_lll.z.xz * tmp5500;
	real tmp5926 = d_lll.z.yy * tmp5710;
	real tmp5932 = d_lll.z.yy * tmp5775;
	real tmp5936 = d_lll.z.yz * tmp5779;
	real tmp5942 = d_lll.z.yz * tmp5788;
	real tmp5945 = gamma_uu.xy * tmp332;
	real tmp5946 = d_lll.z.zz * tmp5945;
	real tmp5951 = gamma_uu.xz * tmp291;
	real tmp5952 = d_lll.z.zz * tmp5951;
	real tmp6050 = d_lll.z.yy * tmp5785;
	real tmp6053 = d_lll.z.yy * tmp5788;
	real tmp6386 = gamma_uu.yy * tmp230;
	real tmp6391 = gamma_uu.zz * tmp218;
	real tmp6436 = gamma_uu.yy * tmp291;
	real tmp6441 = gamma_uu.yz * tmp230;
	real tmp6495 = d_lll.z.yz * tmp6436;
	real tmp6501 = d_lll.z.yz * tmp6441;
	real tmp6504 = gamma_uu.yy * tmp332;
	real tmp6505 = d_lll.z.zz * tmp6504;
	real tmp6510 = gamma_uu.zz * tmp230;
	real tmp6511 = d_lll.z.zz * tmp6510;
	real tmp7192 = Z_l.x * gamma_uu.xx;
	real tmp7194 = Z_l.y * gamma_uu.xy;
	real tmp7196 = Z_l.z * gamma_uu.xz;
	real tmp7258 = Z_l.x * gamma_uu.xy;
	real tmp7260 = Z_l.y * gamma_uu.yy;
	real tmp7262 = Z_l.z * gamma_uu.yz;
	real tmp7341 = Z_l.x * gamma_uu.xz;
	real tmp7343 = Z_l.y * gamma_uu.yz;
	real tmp7345 = Z_l.z * gamma_uu.zz;
	deriv->alpha += -f_alphaSq * (tmp8
		+ tmp4
		- 2. * Theta
		+ 2. * tmp3
		+ 2. * tmp2
		+ 2. * tmp1
		+ K_ll.xx * gamma_uu.xx);
	deriv->gamma_ll.xx += -2. * K_ll.xx * alpha;
	deriv->gamma_ll.xy += -2. * K_ll.xy * alpha;
	deriv->gamma_ll.xz += -2. * K_ll.xz * alpha;
	deriv->gamma_ll.yy += -2. * K_ll.yy * alpha;
	deriv->gamma_ll.yz += -2. * K_ll.yz * alpha;
	deriv->gamma_ll.zz += -2. * K_ll.zz * alpha;
	deriv->a_l.x += -alphaSq_dalpha_f * (
		K_ll.zz * a_l.x * gamma_uu.zz
		- 2. * Theta * a_l.x
		+ 2. * K_ll.yz * a_l.x * gamma_uu.yz
		+ K_ll.yy * a_l.x * gamma_uu.yy
		+ 2. * K_ll.xz * a_l.x * gamma_uu.xz
		+ 2. * K_ll.xy * a_l.x * gamma_uu.xy
		+ K_ll.xx * a_l.x * gamma_uu.xx
	)
	- f_alpha * (
		- 2. * K_ll.zz * d_lll.x.xx * tmp67
		- 4. * K_ll.zz * d_lll.x.xy * tmp128
		- 4. * K_ll.zz * d_lll.x.xz * tmp189
		- 2. * K_ll.zz * d_lll.x.yy * tmp230
		- 4. * K_ll.zz * d_lll.x.yz * tmp291
		- 2. * K_ll.zz * d_lll.x.zz * tmp332
		- 2. * Theta * a_l.x
		+ K_ll.zz * a_l.x * gamma_uu.zz
		- 4. * K_ll.yz * d_lll.x.xx * tmp61
		- 4. * K_ll.yz * d_lll.x.xy * tmp116
		- 4. * K_ll.yz * d_lll.x.xy * tmp122
		- 4. * K_ll.yz * d_lll.x.xz * tmp177
		- 4. * K_ll.yz * d_lll.x.xz * tmp128
		- 4. * K_ll.yz * d_lll.x.yy * tmp224
		- 4. * K_ll.yz * d_lll.x.yz * tmp279
		- 4. * K_ll.yz * d_lll.x.yz * tmp230
		- 4. * K_ll.yz * d_lll.x.zz * tmp291
		+ 2. * K_ll.yz * a_l.x * gamma_uu.yz
		- 2. * K_ll.yy * d_lll.x.xx * tmp55
		- 4. * K_ll.yy * d_lll.x.xy * tmp110
		- 4. * K_ll.yy * d_lll.x.xz * tmp116
		- 2. * K_ll.yy * d_lll.x.yy * tmp218
		- 4. * K_ll.yy * d_lll.x.yz * tmp224
		- 2. * K_ll.yy * d_lll.x.zz * tmp230
		+ K_ll.yy * a_l.x * gamma_uu.yy
		- 4. * K_ll.xz * d_lll.x.xx * tmp49
		- 4. * K_ll.xz * d_lll.x.xy * tmp98
		- 4. * K_ll.xz * d_lll.x.xy * tmp61
		- 4. * K_ll.xz * d_lll.x.xz * tmp159
		- 4. * K_ll.xz * d_lll.x.xz * tmp67
		- 4. * K_ll.xz * d_lll.x.yy * tmp116
		- 4. * K_ll.xz * d_lll.x.yz * tmp177
		- 4. * K_ll.xz * d_lll.x.yz * tmp128
		- 4. * K_ll.xz * d_lll.x.zz * tmp189
		+ 2. * K_ll.xz * a_l.x * gamma_uu.xz
		- 4. * K_ll.xy * d_lll.x.xx * tmp43
		- 4. * K_ll.xy * d_lll.x.xy * tmp86
		- 4. * K_ll.xy * d_lll.x.xy * tmp55
		- 4. * K_ll.xy * d_lll.x.xz * tmp98
		- 4. * K_ll.xy * d_lll.x.xz * tmp61
		- 4. * K_ll.xy * d_lll.x.yy * tmp110
		- 4. * K_ll.xy * d_lll.x.yz * tmp116
		- 4. * K_ll.xy * d_lll.x.yz * tmp122
		- 4. * K_ll.xy * d_lll.x.zz * tmp128
		+ 2. * K_ll.xy * a_l.x * gamma_uu.xy
		- 2. * K_ll.xx * d_lll.x.xx * tmp37
		- 4. * K_ll.xx * d_lll.x.xy * tmp43
		- 4. * K_ll.xx * d_lll.x.xz * tmp49
		- 2. * K_ll.xx * d_lll.x.yy * tmp55
		- 4. * K_ll.xx * d_lll.x.yz * tmp61
		- 2. * K_ll.xx * d_lll.x.zz * tmp67
		+ K_ll.xx * a_l.x * gamma_uu.xx);
	deriv->a_l.y += -alphaSq_dalpha_f * (
		K_ll.zz * a_l.y * gamma_uu.zz
		- 2. * Theta * a_l.y
		+ 2. * K_ll.yz * a_l.y * gamma_uu.yz
		+ K_ll.yy * a_l.y * gamma_uu.yy
		+ 2. * K_ll.xz * a_l.y * gamma_uu.xz
		+ 2. * K_ll.xy * a_l.y * gamma_uu.xy
		+ K_ll.xx * a_l.y * gamma_uu.xx
	)
	- f_alpha * (
		- 2. * K_ll.zz * d_lll.y.xx * tmp67
		- 4. * K_ll.zz * d_lll.y.xy * tmp128
		- 4. * K_ll.zz * d_lll.y.xz * tmp189
		- 2. * K_ll.zz * d_lll.y.yy * tmp230
		- 4. * K_ll.zz * d_lll.y.yz * tmp291
		- 2. * K_ll.zz * d_lll.y.zz * tmp332
		- 2. * Theta * a_l.y
		+ K_ll.zz * a_l.y * gamma_uu.zz
		- 4. * K_ll.yz * d_lll.y.xx * tmp61
		- 4. * K_ll.yz * d_lll.y.xy * tmp116
		- 4. * K_ll.yz * d_lll.y.xy * tmp122
		- 4. * K_ll.yz * d_lll.y.xz * tmp177
		- 4. * K_ll.yz * d_lll.y.xz * tmp128
		- 4. * K_ll.yz * d_lll.y.yy * tmp224
		- 4. * K_ll.yz * d_lll.y.yz * tmp279
		- 4. * K_ll.yz * d_lll.y.yz * tmp230
		- 4. * K_ll.yz * d_lll.y.zz * tmp291
		+ 2. * K_ll.yz * a_l.y * gamma_uu.yz
		- 2. * K_ll.yy * d_lll.y.xx * tmp55
		- 4. * K_ll.yy * d_lll.y.xy * tmp110
		- 4. * K_ll.yy * d_lll.y.xz * tmp116
		- 2. * K_ll.yy * d_lll.y.yy * tmp218
		- 4. * K_ll.yy * d_lll.y.yz * tmp224
		- 2. * K_ll.yy * d_lll.y.zz * tmp230
		+ K_ll.yy * a_l.y * gamma_uu.yy
		- 4. * K_ll.xz * d_lll.y.xx * tmp49
		- 4. * K_ll.xz * d_lll.y.xy * tmp98
		- 4. * K_ll.xz * d_lll.y.xy * tmp61
		- 4. * K_ll.xz * d_lll.y.xz * tmp159
		- 4. * K_ll.xz * d_lll.y.xz * tmp67
		- 4. * K_ll.xz * d_lll.y.yy * tmp116
		- 4. * K_ll.xz * d_lll.y.yz * tmp177
		- 4. * K_ll.xz * d_lll.y.yz * tmp128
		- 4. * K_ll.xz * d_lll.y.zz * tmp189
		+ 2. * K_ll.xz * a_l.y * gamma_uu.xz
		- 4. * K_ll.xy * d_lll.y.xx * tmp43
		- 4. * K_ll.xy * d_lll.y.xy * tmp86
		- 4. * K_ll.xy * d_lll.y.xy * tmp55
		- 4. * K_ll.xy * d_lll.y.xz * tmp98
		- 4. * K_ll.xy * d_lll.y.xz * tmp61
		- 4. * K_ll.xy * d_lll.y.yy * tmp110
		- 4. * K_ll.xy * d_lll.y.yz * tmp116
		- 4. * K_ll.xy * d_lll.y.yz * tmp122
		- 4. * K_ll.xy * d_lll.y.zz * tmp128
		+ 2. * K_ll.xy * a_l.y * gamma_uu.xy
		- 2. * K_ll.xx * d_lll.y.xx * tmp37
		- 4. * K_ll.xx * d_lll.y.xy * tmp43
		- 4. * K_ll.xx * d_lll.y.xz * tmp49
		- 2. * K_ll.xx * d_lll.y.yy * tmp55
		- 4. * K_ll.xx * d_lll.y.yz * tmp61
		- 2. * K_ll.xx * d_lll.y.zz * tmp67
		+ K_ll.xx * a_l.y * gamma_uu.xx);
	deriv->a_l.z += -alphaSq_dalpha_f * (
		K_ll.zz * a_l.z * gamma_uu.zz
		- 2. * Theta * a_l.z
		+ 2. * K_ll.yz * a_l.z * gamma_uu.yz
		+ K_ll.yy * a_l.z * gamma_uu.yy
		+ 2. * K_ll.xz * a_l.z * gamma_uu.xz
		+ 2. * K_ll.xy * a_l.z * gamma_uu.xy
		+ K_ll.xx * a_l.z * gamma_uu.xx
	)
	- f_alpha * (
		- 2. * K_ll.zz * d_lll.z.xx * tmp67
		- 4. * K_ll.zz * d_lll.z.xy * tmp128
		- 4. * K_ll.zz * d_lll.z.xz * tmp189
		- 2. * K_ll.zz * d_lll.z.yy * tmp230
		- 4. * K_ll.zz * d_lll.z.yz * tmp291
		- 2. * K_ll.zz * d_lll.z.zz * tmp332
		- 2. * Theta * a_l.z
		+ K_ll.zz * a_l.z * gamma_uu.zz
		- 4. * K_ll.yz * d_lll.z.xx * tmp61
		- 4. * K_ll.yz * d_lll.z.xy * tmp116
		- 4. * K_ll.yz * d_lll.z.xy * tmp122
		- 4. * K_ll.yz * d_lll.z.xz * tmp177
		- 4. * K_ll.yz * d_lll.z.xz * tmp128
		- 4. * K_ll.yz * d_lll.z.yy * tmp224
		- 4. * K_ll.yz * d_lll.z.yz * tmp279
		- 4. * K_ll.yz * d_lll.z.yz * tmp230
		- 4. * K_ll.yz * d_lll.z.zz * tmp291
		+ 2. * K_ll.yz * a_l.z * gamma_uu.yz
		- 2. * K_ll.yy * d_lll.z.xx * tmp55
		- 4. * K_ll.yy * d_lll.z.xy * tmp110
		- 4. * K_ll.yy * d_lll.z.xz * tmp116
		- 2. * K_ll.yy * d_lll.z.yy * tmp218
		- 4. * K_ll.yy * d_lll.z.yz * tmp224
		- 2. * K_ll.yy * d_lll.z.zz * tmp230
		+ K_ll.yy * a_l.z * gamma_uu.yy
		- 4. * K_ll.xz * d_lll.z.xx * tmp49
		- 4. * K_ll.xz * d_lll.z.xy * tmp98
		- 4. * K_ll.xz * d_lll.z.xy * tmp61
		- 4. * K_ll.xz * d_lll.z.xz * tmp159
		- 4. * K_ll.xz * d_lll.z.xz * tmp67
		- 4. * K_ll.xz * d_lll.z.yy * tmp116
		- 4. * K_ll.xz * d_lll.z.yz * tmp177
		- 4. * K_ll.xz * d_lll.z.yz * tmp128
		- 4. * K_ll.xz * d_lll.z.zz * tmp189
		+ 2. * K_ll.xz * a_l.z * gamma_uu.xz
		- 4. * K_ll.xy * d_lll.z.xx * tmp43
		- 4. * K_ll.xy * d_lll.z.xy * tmp86
		- 4. * K_ll.xy * d_lll.z.xy * tmp55
		- 4. * K_ll.xy * d_lll.z.xz * tmp98
		- 4. * K_ll.xy * d_lll.z.xz * tmp61
		- 4. * K_ll.xy * d_lll.z.yy * tmp110
		- 4. * K_ll.xy * d_lll.z.yz * tmp116
		- 4. * K_ll.xy * d_lll.z.yz * tmp122
		- 4. * K_ll.xy * d_lll.z.zz * tmp128
		+ 2. * K_ll.xy * a_l.z * gamma_uu.xy
		- 2. * K_ll.xx * d_lll.z.xx * tmp37
		- 4. * K_ll.xx * d_lll.z.xy * tmp43
		- 4. * K_ll.xx * d_lll.z.xz * tmp49
		- 2. * K_ll.xx * d_lll.z.yy * tmp55
		- 4. * K_ll.xx * d_lll.z.yz * tmp61
		- 2. * K_ll.xx * d_lll.z.zz * tmp67
		+ K_ll.xx * a_l.z * gamma_uu.xx);
	deriv->d_lll.x.xx += -K_ll.xx * tmp1049;
	deriv->d_lll.x.xy += -K_ll.xy * tmp1049;
	deriv->d_lll.x.xz += -K_ll.xz * tmp1049;
	deriv->d_lll.x.yy += -K_ll.yy * tmp1049;
	deriv->d_lll.x.yz += -K_ll.yz * tmp1049;
	deriv->d_lll.x.zz += -K_ll.zz * tmp1049;
	deriv->d_lll.y.xx += -K_ll.xx * tmp1061;
	deriv->d_lll.y.xy += -K_ll.xy * tmp1061;
	deriv->d_lll.y.xz += -K_ll.xz * tmp1061;
	deriv->d_lll.y.yy += -K_ll.yy * tmp1061;
	deriv->d_lll.y.yz += -K_ll.yz * tmp1061;
	deriv->d_lll.y.zz += -K_ll.zz * tmp1061;
	deriv->d_lll.z.xx += -K_ll.xx * tmp1073;
	deriv->d_lll.z.xy += -K_ll.xy * tmp1073;
	deriv->d_lll.z.xz += -K_ll.xz * tmp1073;
	deriv->d_lll.z.yy += -K_ll.yy * tmp1073;
	deriv->d_lll.z.yz += -K_ll.yz * tmp1073;
	deriv->d_lll.z.zz += -K_ll.zz * tmp1073;
	deriv->K_ll.xx += -alpha * (2. * tmp1549
		+ 2. * tmp1532
		- 2. * tmp1535
		- 2. * tmp1540
		- 2. * tmp1545
		+ gamma_uu.xx * K_ll.xx * K_ll.xx
		- gamma_uu.xx * gamma_uu.yy * tmp1523
		- gamma_uu.xx * gamma_uu.zz * tmp1527
		+ 2. * tmp1521
		+ d_lll.z.xx * tmp1326
		- 2. * d_lll.z.xx * tmp1331
		- 2. * d_lll.z.xx * tmp1334
		- d_lll.z.xx * tmp1337
		+ 2. * tmp1479
		- d_lll.y.yy * d_lll.z.xx * tmp224
		- 2. * d_lll.y.yz * d_lll.z.xx * tmp279
		- d_lll.y.zz * d_lll.z.xx * tmp291
		- 2. * d_lll.z.xx * tmp1494
		- 2. * d_lll.z.xx * tmp1323
		+ 4. * d_lll.y.xz * tmp1470
		- 4. * d_lll.y.xz * tmp1473
		+ 2. * d_lll.y.xz * tmp1375
		- 4. * d_lll.y.xz * tmp1465
		+ 2. * d_lll.y.xx * tmp1436
		- 2. * d_lll.y.xx * tmp1278
		- d_lll.y.xx * tmp1281
		- 2. * d_lll.y.xx * tmp1284
		- d_lll.y.xx * tmp1287
		- 2. * d_lll.y.xy * tmp1359
		+ d_lll.y.xx * tmp1267
		- 2. * d_lll.y.xx * tmp1272
		- 2. * d_lll.y.xx * tmp1426
		- 4. * d_lll.y.xx * tmp1275
		+ d_lll.x.zz * d_lll.y.xx * tmp177
		- 2. * d_lll.x.zz * d_lll.y.xx * tmp128
		- d_lll.x.zz * d_lll.z.xx * tmp189
		- tmp1398
		- 2. * d_lll.y.xx * tmp1255
		- 2. * d_lll.y.xx * tmp1405
		- d_lll.y.xx * tmp1261
		- 2. * d_lll.y.xx * tmp1264
		+ d_lll.x.yy * tmp1359
		- tmp1363
		- 4. * d_lll.x.yz * tmp1365
		- 2. * d_lll.x.yz * d_lll.y.xx * tmp122
		- 2. * d_lll.x.yz * tmp1375
		- 2. * tmp1381
		+ 2. * d_lll.x.xz * tmp1337
		- 4. * d_lll.x.yy * tmp1340
		- 2. * d_lll.x.yy * tmp1345
		- d_lll.x.yy * d_lll.y.xx * tmp110
		- 2. * d_lll.x.yy * d_lll.z.xx * tmp116
		+ 4. * d_lll.x.xz * tmp1334
		+ 4. * d_lll.x.xz * tmp1331
		+ 4. * d_lll.x.xz * tmp1323
		- 2. * d_lll.x.xz * tmp1326
		+ 4. * d_lll.x.xz * tmp1320
		+ 2. * d_lll.x.xz * tmp1317
		+ 4. * d_lll.x.xz * tmp1314
		+ 2. * d_lll.x.xz * tmp1311
		+ 4. * d_lll.x.xz * tmp1308
		+ 4. * d_lll.x.xz * tmp1305
		+ 2. * d_lll.x.xy * tmp1287
		- 2. * d_lll.x.xz * tmp1290
		- 4. * d_lll.x.xz * tmp1295
		- 2. * d_lll.x.xz * tmp1300
		+ 4. * d_lll.x.xy * tmp1284
		+ 2. * d_lll.x.xy * tmp1281
		+ 4. * d_lll.x.xy * tmp1278
		+ 4. * d_lll.x.xy * tmp1275
		+ 4. * d_lll.x.xy * tmp1272
		+ 4. * d_lll.x.xy * tmp1264
		- 2. * d_lll.x.xy * tmp1267
		+ 2. * d_lll.x.xy * tmp1261
		+ 4. * d_lll.x.xy * tmp1258
		+ 4. * d_lll.x.xy * tmp1255
		+ d_lll.x.xx * tmp1238
		- 2. * d_lll.x.xy * tmp1240
		- 4. * d_lll.x.xy * tmp1245
		- 2. * d_lll.x.xy * tmp1250
		+ 2. * d_lll.x.xx * tmp1235
		+ 2. * d_lll.x.xx * tmp1228
		- d_lll.x.xx * tmp1231
		+ 2. * d_lll.x.xx * tmp1225
		+ 2. * d_lll.x.xx * tmp1222
		+ 2. * d_lll.x.xx * tmp1219
		+ 2. * d_lll.x.xx * tmp1212
		- d_lll.x.xx * tmp1215
		+ d_lll.x.xx * tmp1210
		+ 2. * d_lll.x.xx * tmp1207
		+ 2. * d_lll.x.xx * tmp1204
		+ a_l.x * a_l.x
		- d_lll.x.xx * tmp1191
		- 2. * d_lll.x.xx * tmp1195
		- d_lll.x.xx * tmp1200
		+ a_l.z * tmp1147
		+ a_l.z * tmp1143
		+ a_l.y * tmp1133
		- a_l.z * tmp1137
		- 2. * a_l.z * tmp1139
		- 2. * a_l.z * tmp1141
		+ a_l.y * tmp1129
		+ a_l.x * tmp1119
		- a_l.y * tmp1123
		- 2. * a_l.y * tmp1125
		- 2. * a_l.y * tmp1127
		+ a_l.x * tmp1115
		+ 4. * Z_l.z * tmp1141
		- 2. * Z_l.z * tmp1143
		- 2. * Z_l.z * tmp1147
		- a_l.x * tmp1109
		- 2. * a_l.x * tmp1111
		- 2. * a_l.x * tmp1113
		+ 4. * Z_l.z * tmp1139
		+ 2. * Z_l.z * tmp1137
		+ 4. * Z_l.y * tmp1127
		- 2. * Z_l.y * tmp1129
		- 2. * Z_l.y * tmp1133
		+ 4. * Z_l.y * tmp1125
		+ 2. * Z_l.y * tmp1123
		+ 4. * Z_l.x * tmp1113
		- 2. * Z_l.x * tmp1115
		- 2. * Z_l.x * tmp1119
		+ 4. * Z_l.x * tmp1111
		+ 2. * Z_l.x * tmp1109
		+ 4. * M_PI * gamma_ll.xx * rho
		+ 8. * M_PI * S_ll.xx
		+ 4. * K_ll.xy * K_ll.xz * gamma_uu.yz
		- 4. * M_PI * S * gamma_ll.xx
		+ 2. * K_ll.xx * Theta
		+ 2. * K_ll.xx * tmp2
		- K_ll.xx * tmp8
		- 2. * K_ll.xx * tmp3
		- K_ll.xx * tmp4
		+ 2. * K_ll.xx * tmp1);
	deriv->K_ll.xy += -alpha * (2. * gamma_uu.xy * tmp1539
		- 2. * tmp2162
		+ 2. * gamma_uu.xy * tmp1534
		+ gamma_uu.xy * tmp2153
		+ gamma_uu.xx * gamma_uu.xy * tmp1523
		+ 2. * d_lll.z.xx * tmp1888
		- 2. * d_lll.z.xy * tmp1323
		- d_lll.z.xy * tmp1326
		- 2. * d_lll.z.xy * tmp1334
		- d_lll.z.xy * tmp1337
		+ 2. * d_lll.y.yz * tmp1375
		- 2. * d_lll.y.yz * tmp1465
		- 2. * d_lll.y.yz * tmp1473
		- d_lll.y.zz * d_lll.z.xy * tmp291
		- d_lll.z.xx * tmp2121
		- 2. * d_lll.z.xx * tmp1885
		+ d_lll.y.xz * tmp1337
		- d_lll.y.yy * tmp2100
		+ 2. * d_lll.y.xz * tmp1334
		+ d_lll.y.xz * tmp1326
		+ 2. * d_lll.y.xz * tmp1323
		+ 2. * d_lll.y.xz * tmp1320
		+ 2. * d_lll.y.xz * tmp2079
		- 2. * d_lll.y.xz * tmp1494
		+ d_lll.y.xz * tmp1317
		- d_lll.y.xz * tmp1961
		+ 2. * d_lll.y.xz * tmp1952
		+ d_lll.y.xz * tmp1311
		+ 2. * d_lll.y.xy * tmp2060
		- 2. * d_lll.y.xy * tmp1275
		+ 2. * d_lll.y.xy * tmp1258
		- 2. * d_lll.y.xy * tmp1426
		+ d_lll.y.xx * tmp1238
		+ 2. * d_lll.y.xx * tmp1235
		+ d_lll.y.xx * tmp1231
		+ 2. * d_lll.y.xx * tmp1225
		+ d_lll.y.xx * tmp1222
		+ d_lll.y.xx * tmp2038
		+ 2. * d_lll.y.xx * tmp1219
		+ 2. * d_lll.y.xx * tmp1826
		- d_lll.y.xx * tmp1215
		+ d_lll.y.xx * tmp1210
		+ 2. * d_lll.y.xx * tmp1816
		+ d_lll.y.xx * tmp1207
		+ 2. * d_lll.y.xx * tmp1811
		+ d_lll.x.zz * tmp1990
		- 2. * d_lll.x.zz * tmp1992
		- d_lll.x.zz * tmp1997
		- d_lll.x.zz * tmp2001
		- 2. * d_lll.x.zz * tmp2005
		- d_lll.x.zz * tmp2010
		- d_lll.x.zz * d_lll.z.xy * tmp189
		+ d_lll.x.yz * tmp1337
		- d_lll.x.zz * d_lll.y.xx * tmp159
		+ 2. * d_lll.x.yz * tmp1334
		+ 2. * d_lll.x.yz * tmp1331
		+ 2. * d_lll.x.yz * tmp1323
		- d_lll.x.yz * tmp1326
		+ 2. * d_lll.x.yz * tmp1320
		+ d_lll.x.yz * tmp1961
		- 2. * d_lll.x.yz * tmp1494
		+ 2. * d_lll.x.yz * tmp1940
		- 2. * d_lll.x.yz * tmp1943
		- d_lll.x.yz * tmp1311
		- 2. * d_lll.x.yz * tmp1952
		- d_lll.x.yz * tmp1317
		+ d_lll.x.yz * tmp1300
		- d_lll.x.yz * tmp1936
		+ d_lll.x.yy * tmp1287
		+ 2. * d_lll.x.yy * tmp1284
		+ d_lll.x.yy * tmp1281
		+ 2. * d_lll.x.yy * tmp1278
		+ d_lll.x.yy * tmp1436
		+ d_lll.x.yy * tmp1426
		+ d_lll.x.yy * tmp1272
		+ 2. * d_lll.x.yy * tmp1258
		- d_lll.x.yy * tmp1405
		- d_lll.x.yy * tmp1267
		+ d_lll.x.yy * tmp1905
		+ 2. * d_lll.x.yy * tmp1902
		+ d_lll.x.yy * tmp1896
		- d_lll.x.yy * tmp1250
		+ 2. * d_lll.x.yy * tmp1245
		+ 2. * d_lll.x.xz * tmp1885
		- 2. * d_lll.x.xz * tmp1888
		+ 2. * d_lll.x.xz * d_lll.x.yz * tmp67
		- 4. * d_lll.x.xz * tmp1855
		- 2. * d_lll.x.xz * tmp1860
		- 2. * d_lll.x.xz * tmp1865
		- 4. * d_lll.x.xz * tmp1870
		- 2. * d_lll.x.xz * tmp1875
		- 2. * d_lll.x.xz * tmp1880
		+ 2. * d_lll.x.xz * d_lll.x.yy * tmp61
		+ 2. * d_lll.x.xy * tmp1228
		- 2. * d_lll.x.xy * tmp1231
		+ 2. * d_lll.x.xy * tmp1808
		- 4. * d_lll.x.xy * tmp1811
		- 2. * d_lll.x.xy * tmp1816
		- 2. * d_lll.x.xy * tmp1210
		- 4. * d_lll.x.xy * tmp1826
		- 2. * d_lll.x.xy * tmp1219
		- 2. * d_lll.x.xy * tmp1836
		+ 2. * d_lll.x.xy * tmp1805
		+ d_lll.x.xx * tmp1777
		- 2. * d_lll.x.xx * tmp1779
		- d_lll.x.xx * tmp1784
		- d_lll.x.xx * tmp1788
		- 2. * d_lll.x.xx * tmp1792
		- d_lll.x.xx * tmp1797
		- d_lll.x.xx * tmp1801
		+ d_lll.x.xx * tmp1775
		+ a_l.z * tmp1731
		+ a_l.y * tmp1719
		- a_l.z * tmp1723
		- a_l.z * tmp1725
		- a_l.z * tmp1727
		- a_l.z * tmp1729
		+ a_l.x * tmp1707
		- a_l.y * tmp1711
		- a_l.y * tmp1713
		- a_l.y * tmp1115
		- a_l.y * tmp1717
		+ a_l.x * a_l.y
		- a_l.x * tmp1699
		- a_l.x * tmp1701
		- a_l.x * tmp1703
		- a_l.x * tmp1705
		+ 2. * Z_l.z * tmp1729
		- 2. * Z_l.z * tmp1731
		+ 2. * Z_l.z * tmp1727
		+ 2. * Z_l.z * tmp1725
		+ 2. * Z_l.z * tmp1723
		+ 2. * Z_l.y * tmp1717
		- 2. * Z_l.y * tmp1719
		+ 2. * Z_l.y * tmp1115
		+ 2. * Z_l.y * tmp1713
		+ 2. * Z_l.y * tmp1711
		+ 2. * Z_l.x * tmp1705
		- 2. * Z_l.x * tmp1707
		+ 2. * Z_l.x * tmp1703
		+ 2. * Z_l.x * tmp1701
		+ 2. * Z_l.x * tmp1699
		+ 4. * M_PI * gamma_ll.xy * rho
		+ 8. * M_PI * S_ll.xy
		+ 2. * K_ll.xz * K_ll.yz * gamma_uu.zz
		- 4. * M_PI * S * gamma_ll.xy
		+ 2. * K_ll.xz * K_ll.yy * gamma_uu.yz
		+ 2. * K_ll.xy * Theta
		+ K_ll.xy * tmp8
		- K_ll.xy * tmp4
		+ 2. * K_ll.xx * tmp1681
		+ 2. * K_ll.xx * tmp1679
		+ K_ll.xx * K_ll.xy * gamma_uu.xx);
	deriv->K_ll.xz += -alpha * (gamma_uu.xz * tmp2817
		+ 2. * gamma_uu.xz * tmp2814
		+ 2. * gamma_uu.xz * tmp2811
		+ gamma_uu.xx * gamma_uu.xz * tmp1527
		- 2. * gamma_uu.xy * tmp2806
		+ d_lll.z.xy * tmp1287
		+ 2. * d_lll.z.xy * tmp2633
		+ d_lll.z.xy * tmp1281
		+ 2. * d_lll.z.xy * tmp2624
		+ d_lll.z.xx * tmp1238
		+ 2. * d_lll.z.xx * tmp2536
		+ 2. * d_lll.z.xx * tmp1228
		- d_lll.z.xx * tmp1231
		+ 2. * d_lll.z.xx * tmp2526
		+ 2. * d_lll.z.xx * tmp1836
		+ d_lll.z.xx * tmp1222
		+ d_lll.y.zz * tmp1470
		+ d_lll.y.zz * tmp1375
		+ 2. * d_lll.y.yz * tmp2100
		+ 2. * d_lll.y.yz * tmp1359
		+ d_lll.y.yy * d_lll.z.xy * tmp218
		+ d_lll.y.yy * d_lll.z.xx * tmp110
		+ 2. * d_lll.y.xz * tmp1275
		- 2. * d_lll.y.xz * tmp1436
		- 2. * d_lll.y.xz * tmp2624
		- d_lll.y.xz * tmp1281
		- 2. * d_lll.y.xz * tmp2633
		- d_lll.y.xz * tmp1287
		+ d_lll.y.xz * tmp1426
		+ 2. * d_lll.y.xy * tmp2565
		- d_lll.y.xz * tmp1261
		- 2. * d_lll.y.xz * tmp1264
		- d_lll.y.xz * tmp1267
		+ 2. * d_lll.y.xy * tmp2559
		+ 2. * d_lll.y.xx * tmp2492
		- 2. * d_lll.y.xy * tmp2555
		+ 2. * d_lll.y.xx * tmp2482
		- 2. * d_lll.y.xx * tmp2700
		+ 2. * d_lll.y.xx * tmp2477
		- 2. * d_lll.y.xx * tmp2692
		+ d_lll.y.xx * tmp2683
		- d_lll.y.xx * tmp2685
		+ 2. * d_lll.y.xx * tmp2469
		- 2. * d_lll.y.xx * tmp2474
		+ d_lll.x.zz * tmp1331
		- d_lll.y.xx * tmp2671
		+ 2. * d_lll.x.zz * tmp1320
		- d_lll.x.zz * tmp1326
		+ d_lll.x.zz * tmp2079
		- d_lll.x.zz * tmp1494
		+ d_lll.x.zz * tmp1317
		+ 2. * d_lll.x.zz * tmp1314
		+ d_lll.x.zz * tmp1311
		+ d_lll.x.zz * tmp2647
		+ 2. * d_lll.x.zz * tmp1305
		+ d_lll.x.zz * tmp1936
		+ 2. * d_lll.x.yz * tmp2060
		- 2. * d_lll.x.yz * tmp2624
		- d_lll.x.yz * tmp1281
		- 2. * d_lll.x.yz * tmp2633
		- d_lll.x.yz * tmp1287
		+ 2. * d_lll.x.yz * tmp1272
		- d_lll.x.yz * tmp1426
		+ 2. * d_lll.x.yz * tmp1264
		- d_lll.x.yz * tmp1267
		+ d_lll.x.yz * tmp1261
		+ 2. * d_lll.x.yz * tmp1258
		- 2. * d_lll.x.yz * tmp1405
		+ 2. * d_lll.x.yz * tmp1255
		+ d_lll.x.yz * tmp2592
		+ 2. * d_lll.x.yz * tmp1902
		+ d_lll.x.yz * tmp1250
		+ d_lll.x.yy * tmp2563
		- d_lll.x.yy * tmp2565
		- 2. * d_lll.x.yy * tmp2569
		- d_lll.x.yy * tmp2574
		- 2. * d_lll.x.yy * tmp2578
		- d_lll.x.yy * tmp2583
		+ 2. * d_lll.x.yy * d_lll.x.zz * tmp116
		- d_lll.x.yy * d_lll.x.zz * tmp122
		- d_lll.x.yy * tmp2555
		- d_lll.x.yy * tmp2559
		+ d_lll.x.yy * d_lll.x.yz * tmp110
		+ 2. * d_lll.x.xz * tmp1219
		- 2. * d_lll.x.xz * tmp1836
		- 4. * d_lll.x.xz * tmp2526
		- 2. * d_lll.x.xz * tmp1228
		- 4. * d_lll.x.xz * tmp2536
		- 2. * d_lll.x.xz * tmp1238
		+ 2. * d_lll.x.xz * tmp2505
		- 2. * d_lll.x.xz * tmp1816
		- 2. * d_lll.x.xz * tmp1215
		+ 2. * d_lll.x.xz * tmp1808
		+ 2. * d_lll.x.xy * tmp2474
		- 2. * d_lll.x.xy * tmp2477
		- 4. * d_lll.x.xy * tmp2482
		- 2. * d_lll.x.xy * tmp2487
		- 4. * d_lll.x.xy * tmp2492
		- 2. * d_lll.x.xy * tmp2497
		+ 2. * d_lll.x.xy * tmp2461
		- 2. * d_lll.x.xy * tmp2464
		- 2. * d_lll.x.xy * tmp2469
		+ 2. * d_lll.x.xy * d_lll.x.yz * tmp55
		+ d_lll.x.xx * tmp2430
		- d_lll.x.xx * tmp2432
		- d_lll.x.xx * tmp2436
		- 2. * d_lll.x.xx * tmp2440
		- d_lll.x.xx * tmp2445
		- 2. * d_lll.x.xx * tmp2449
		- d_lll.x.xx * tmp2454
		+ d_lll.x.xx * tmp2428
		+ a_l.z * tmp1717
		- a_l.z * tmp1119
		- a_l.z * tmp1719
		+ a_l.y * tmp2368
		- a_l.y * tmp2372
		- a_l.y * tmp2374
		- a_l.z * tmp1713
		- a_l.z * tmp2378
		+ a_l.x * tmp2356
		- a_l.x * tmp2360
		- a_l.x * tmp2362
		- a_l.y * tmp2364
		- a_l.y * tmp2366
		+ a_l.x * a_l.z
		- a_l.x * tmp2352
		- a_l.x * tmp2354
		+ 2. * Z_l.z * tmp1719
		+ 2. * Z_l.z * tmp1119
		+ 2. * Z_l.z * tmp2378
		- 2. * Z_l.z * tmp1717
		+ 2. * Z_l.z * tmp1713
		+ 2. * Z_l.y * tmp2374
		+ 2. * Z_l.y * tmp2372
		+ 2. * Z_l.y * tmp2366
		- 2. * Z_l.y * tmp2368
		+ 2. * Z_l.y * tmp2364
		+ 2. * Z_l.x * tmp2362
		+ 2. * Z_l.x * tmp2360
		+ 2. * Z_l.x * tmp2354
		- 2. * Z_l.x * tmp2356
		+ 2. * Z_l.x * tmp2352
		+ 4. * M_PI * gamma_ll.xz * rho
		+ 8. * M_PI * S_ll.xz
		+ 2. * K_ll.xz * Theta
		- 4. * M_PI * S * gamma_ll.xz
		+ K_ll.xz * tmp4
		+ 2. * K_ll.xy * tmp2338
		- K_ll.xz * tmp8
		+ 2. * K_ll.xy * tmp2336
		+ 2. * K_ll.xx * tmp2334
		+ 2. * K_ll.xx * tmp2332
		+ K_ll.xx * tmp2331);
	deriv->K_ll.yy += alpha * (tmp3479
		+ gamma_uu.yy * gamma_uu.zz * tmp3429
		- 2. * tmp3432
		+ 2. * gamma_uu.xx * tmp1544
		- gamma_uu.yy * K_ll.yy * K_ll.yy
		+ 2. * gamma_uu.xx * tmp1539
		+ 2. * gamma_uu.xx * tmp1534
		+ gamma_uu.xx * tmp2153
		+ d_lll.z.yy * tmp1337
		- 2. * gamma_uu.xx * tmp1531
		+ 2. * d_lll.z.yy * tmp1334
		+ 2. * d_lll.z.xz * tmp3255
		+ 2. * d_lll.z.xy * tmp1885
		- 2. * tmp3400
		+ 2. * d_lll.z.xx * d_lll.z.yy * tmp67
		+ tmp3390
		- d_lll.z.xx * tmp3153
		+ d_lll.y.zz * d_lll.z.yy * tmp291
		+ 2. * d_lll.y.yz * tmp1961
		- 4. * d_lll.y.yz * tmp2079
		- 4. * d_lll.y.yz * tmp1320
		- 4. * d_lll.y.yz * tmp1323
		- 4. * d_lll.y.yz * tmp1334
		- 2. * d_lll.y.yz * tmp1337
		+ 2. * d_lll.y.yz * tmp1317
		+ d_lll.y.yy * tmp1426
		- 2. * d_lll.y.yy * tmp2060
		- 2. * d_lll.y.yy * tmp1436
		- 2. * d_lll.y.yy * tmp1278
		- 2. * d_lll.y.yy * tmp1284
		- d_lll.y.yy * tmp1287
		+ d_lll.y.yy * tmp1267
		+ 2. * tmp3324
		+ 2. * d_lll.y.xz * tmp1885
		+ 4. * d_lll.y.xz * tmp1875
		+ 4. * d_lll.y.xz * tmp3209
		+ 2. * d_lll.y.xz * tmp3204
		+ 2. * d_lll.y.xy * tmp1215
		- 2. * d_lll.y.xy * tmp2038
		- 4. * d_lll.y.xy * tmp1836
		- 4. * d_lll.y.xy * tmp1225
		- 4. * d_lll.y.xy * tmp1235
		- 2. * d_lll.y.xy * tmp1238
		+ 4. * d_lll.y.xy * tmp1816
		+ 2. * d_lll.y.xx * tmp3275
		+ 2. * d_lll.y.xx * tmp1797
		- d_lll.y.xx * tmp3132
		+ 2. * d_lll.y.xx * tmp3127
		+ d_lll.y.xx * tmp3122
		+ 4. * d_lll.y.xx * tmp1784
		+ 2. * d_lll.y.xx * tmp1779
		+ d_lll.x.zz * tmp3255
		+ d_lll.x.zz * d_lll.y.yy * tmp177
		- 2. * d_lll.x.zz * d_lll.y.yy * tmp128
		- 2. * d_lll.x.zz * d_lll.y.yz * tmp189
		+ 2. * d_lll.x.zz * d_lll.y.xy * tmp159
		- 4. * d_lll.x.zz * d_lll.y.xy * tmp67
		+ 4. * d_lll.x.yz * tmp1888
		- 2. * tmp3231
		+ 4. * d_lll.x.yz * tmp1880
		- 2. * d_lll.x.yz * tmp1885
		+ d_lll.x.yy * tmp1238
		- 4. * d_lll.x.yz * tmp1855
		- 2. * d_lll.x.yz * tmp3204
		- 4. * d_lll.x.yz * tmp3209
		- 4. * d_lll.x.yz * tmp2121
		+ 2. * d_lll.x.yy * tmp1235
		+ 2. * d_lll.x.yy * tmp1231
		+ 2. * d_lll.x.yy * tmp1225
		+ 4. * d_lll.x.yy * tmp1836
		+ d_lll.x.yy * tmp2038
		- 2. * d_lll.x.yy * tmp1222
		+ 2. * d_lll.x.yy * tmp1219
		+ 2. * d_lll.x.yy * tmp1207
		- d_lll.x.yy * tmp1215
		+ d_lll.x.yy * tmp3166
		+ 2. * d_lll.x.yy * tmp2505
		+ 2. * d_lll.x.yy * tmp1195
		- d_lll.x.yy * tmp1200
		+ 2. * d_lll.x.xz * tmp3153
		+ 2. * d_lll.x.xz * d_lll.x.yy * tmp49
		- 4. * d_lll.x.xz * d_lll.y.xy * tmp49
		- 2. * d_lll.x.xz * d_lll.y.yy * tmp98
		- 4. * d_lll.x.xz * d_lll.y.yz * tmp159
		+ 2. * d_lll.x.xy * tmp3132
		+ 2. * d_lll.x.xy * tmp1775
		- 4. * d_lll.x.xy * tmp1779
		- 2. * d_lll.x.xy * tmp3122
		- 4. * d_lll.x.xy * tmp3127
		+ d_lll.x.xx * d_lll.z.yy * tmp49
		+ d_lll.x.xx * d_lll.x.yy * tmp37
		- 2. * d_lll.x.xx * d_lll.y.xy * tmp37
		- d_lll.x.xx * d_lll.y.yy * tmp43
		- 2. * d_lll.x.xx * d_lll.y.yz * tmp49
		+ 2. * a_l.z * tmp3054
		- a_l.z * tmp3058
		- a_l.y * a_l.y
		+ a_l.z * tmp3050
		+ 2. * a_l.z * tmp3046
		+ 2. * a_l.y * tmp3038
		- a_l.y * tmp3042
		- a_l.z * tmp3044
		+ a_l.y * tmp3034
		+ 2. * a_l.y * tmp3030
		+ 2. * a_l.x * tmp3022
		- a_l.x * tmp3026
		- a_l.y * tmp1699
		+ a_l.x * tmp3018
		+ 2. * a_l.x * tmp3014
		+ 2. * Z_l.z * tmp3058
		- a_l.x * tmp3012
		+ 2. * Z_l.z * tmp3044
		- 4. * Z_l.z * tmp3046
		- 2. * Z_l.z * tmp3050
		- 4. * Z_l.z * tmp3054
		+ 2. * Z_l.y * tmp3042
		+ 2. * Z_l.y * tmp1699
		- 4. * Z_l.y * tmp3030
		- 2. * Z_l.y * tmp3034
		- 4. * Z_l.y * tmp3038
		+ 2. * Z_l.x * tmp3026
		+ 2. * Z_l.x * tmp3012
		- 4. * Z_l.x * tmp3014
		- 2. * Z_l.x * tmp3018
		- 4. * Z_l.x * tmp3022
		+ 4. * M_PI * S * gamma_ll.yy
		- 8. * M_PI * S_ll.yy
		- 4. * M_PI * gamma_ll.yy * rho
		+ K_ll.yy * tmp4
		- 2. * K_ll.yy * Theta
		+ 2. * K_ll.xz * K_ll.yy * gamma_uu.xz
		- 2. * K_ll.yy * tmp3
		+ K_ll.xx * K_ll.yy * gamma_uu.xx
		- 2. * K_ll.xy * tmp1679
		- 4. * K_ll.xy * tmp1681);
	deriv->K_ll.yz += alpha * (2. * gamma_uu.xy * gamma_uu.xz * tmp1379
		- gamma_uu.yy * gamma_uu.yz * tmp3429
		- gamma_uu.yz * tmp4086
		+ d_lll.z.xx * tmp3132
		- 2. * d_lll.z.xx * tmp3275
		- 2. * d_lll.z.xy * tmp2526
		- 2. * d_lll.z.xy * tmp1228
		- d_lll.z.xy * tmp1231
		- 2. * d_lll.z.xy * tmp2536
		- d_lll.z.xy * tmp1238
		- 2. * d_lll.z.xz * tmp1888
		- 2. * d_lll.z.yy * tmp2633
		- d_lll.z.yy * tmp1287
		- 2. * gamma_uu.xx * tmp2806
		- 2. * gamma_uu.xx * tmp2161
		+ d_lll.y.zz * tmp1494
		- 2. * d_lll.y.zz * tmp1320
		- d_lll.y.zz * tmp1331
		- d_lll.z.xx * tmp1801
		+ d_lll.y.zz * tmp1961
		- d_lll.y.zz * tmp2079
		+ 2. * d_lll.y.yz * tmp1287
		+ 4. * d_lll.y.yz * tmp2633
		+ 4. * d_lll.y.yz * tmp2624
		+ 2. * d_lll.y.yz * tmp1275
		+ 2. * d_lll.y.yz * tmp2060
		+ d_lll.y.yy * tmp2583
		- 2. * d_lll.y.yz * tmp1272
		+ 2. * d_lll.y.yy * tmp2578
		+ 2. * d_lll.y.yy * tmp2569
		+ d_lll.y.yy * tmp2565
		+ d_lll.y.yy * tmp2563
		+ d_lll.y.xz * tmp1238
		- d_lll.y.yy * tmp3968
		+ 2. * d_lll.y.xz * tmp2536
		+ d_lll.y.xz * tmp1231
		+ 2. * d_lll.y.xz * tmp2526
		- 2. * d_lll.y.xz * tmp1228
		+ d_lll.y.xz * tmp2038
		+ 2. * d_lll.y.xy * tmp2497
		- d_lll.y.xz * tmp1210
		- 2. * d_lll.y.xz * tmp1826
		- d_lll.y.xz * tmp1215
		- 2. * d_lll.y.xz * tmp1219
		+ 4. * d_lll.y.xy * tmp2700
		+ 4. * d_lll.y.xy * tmp2692
		+ 2. * d_lll.y.xy * tmp2477
		+ 2. * d_lll.y.xy * tmp2683
		+ d_lll.y.xx * tmp2454
		- 2. * d_lll.y.xy * tmp2464
		- 2. * d_lll.y.xy * tmp2469
		+ 2. * d_lll.y.xx * tmp2449
		+ d_lll.y.xx * tmp3725
		- d_lll.y.xx * tmp2445
		+ 2. * d_lll.y.xx * tmp2440
		+ d_lll.y.xx * tmp2436
		+ d_lll.y.xx * tmp3890
		+ d_lll.y.xx * tmp3715
		- 2. * d_lll.y.xx * d_lll.y.zz * tmp61
		+ 2. * d_lll.x.zz * tmp1870
		- 2. * d_lll.x.zz * tmp3209
		- d_lll.x.zz * tmp1875
		- d_lll.x.zz * tmp2121
		- d_lll.x.zz * tmp1885
		- d_lll.y.xx * tmp2432
		+ d_lll.x.zz * tmp3852
		- 2. * d_lll.x.zz * tmp1860
		+ 2. * d_lll.x.zz * tmp1855
		+ d_lll.x.yz * tmp1238
		- 2. * d_lll.x.zz * tmp3844
		+ 2. * d_lll.x.yz * tmp2536
		+ 2. * d_lll.x.yz * tmp2526
		- d_lll.x.yz * tmp1231
		+ 2. * d_lll.x.yz * tmp1222
		- 2. * d_lll.x.yz * tmp1836
		+ d_lll.x.yz * tmp2038
		+ 2. * d_lll.x.yz * tmp1826
		- d_lll.x.yz * tmp1215
		+ d_lll.x.yz * tmp1210
		+ 2. * d_lll.x.yz * tmp1207
		- 2. * d_lll.x.yz * tmp1816
		+ 2. * d_lll.x.yz * tmp1811
		+ d_lll.x.yz * tmp3166
		+ d_lll.x.yz * tmp1200
		+ 2. * d_lll.x.yy * tmp2492
		+ 2. * d_lll.x.yy * tmp2482
		- d_lll.x.yy * tmp2487
		- 2. * d_lll.x.yy * tmp2700
		+ d_lll.x.yy * tmp2685
		- 2. * d_lll.x.yy * tmp2477
		- 2. * d_lll.x.yy * tmp2692
		+ 2. * d_lll.x.yy * tmp3755
		- 2. * d_lll.x.yy * tmp2461
		- d_lll.x.yy * tmp2671
		- d_lll.x.yy * tmp2474
		+ d_lll.x.yy * d_lll.x.yz * tmp86
		+ 2. * d_lll.x.xz * tmp1777
		- 2. * d_lll.x.xz * tmp1784
		- 2. * d_lll.x.xz * d_lll.y.zz * tmp159
		- 2. * d_lll.x.xz * tmp1801
		- 2. * d_lll.x.xz * tmp3132
		+ 2. * d_lll.x.xy * tmp2428
		- 2. * d_lll.x.xy * tmp2432
		- 2. * d_lll.x.xy * tmp3715
		- 2. * d_lll.x.xy * tmp2436
		- 2. * d_lll.x.xy * tmp3725
		+ d_lll.x.xx * d_lll.x.yz * tmp37
		- d_lll.x.xx * d_lll.y.xz * tmp37
		- d_lll.x.xx * d_lll.y.zz * tmp49
		- d_lll.x.xx * d_lll.z.xy * tmp37
		- d_lll.x.xx * d_lll.z.yy * tmp43
		+ a_l.z * tmp3042
		+ a_l.z * tmp1707
		+ a_l.z * tmp3653
		+ a_l.z * tmp1705
		+ a_l.y * tmp3643
		- a_l.z * tmp1701
		+ a_l.y * tmp2362
		+ a_l.y * tmp3635
		+ a_l.y * tmp2356
		+ a_l.x * tmp3625
		- a_l.y * a_l.z
		- a_l.y * tmp2352
		+ a_l.x * tmp3621
		+ a_l.x * tmp3617
		+ a_l.x * tmp3613
		+ 2. * Z_l.z * tmp1701
		- 2. * Z_l.z * tmp1705
		- 2. * Z_l.z * tmp3653
		- 2. * Z_l.z * tmp1707
		- 2. * Z_l.z * tmp3042
		- a_l.x * tmp3611
		+ 2. * Z_l.y * tmp2352
		- 2. * Z_l.y * tmp2356
		- 2. * Z_l.y * tmp3635
		- 2. * Z_l.y * tmp2362
		- 2. * Z_l.y * tmp3643
		+ 2. * Z_l.x * tmp3611
		- 2. * Z_l.x * tmp3613
		- 2. * Z_l.x * tmp3617
		- 2. * Z_l.x * tmp3621
		- 2. * Z_l.x * tmp3625
		+ 4. * M_PI * S * gamma_ll.yz
		- 8. * M_PI * S_ll.yz
		- 4. * M_PI * gamma_ll.yz * rho
		+ K_ll.xx * K_ll.yz * gamma_uu.xx
		- 2. * K_ll.xy * tmp2331
		- 2. * K_ll.xy * tmp2334
		- 2. * K_ll.xz * tmp1679
		- K_ll.yy * tmp2336
		- 2. * K_ll.yy * tmp2338
		- K_ll.yz * tmp4
		- 2. * K_ll.yz * Theta);
	deriv->K_ll.zz += alpha * (tmp4705
		+ gamma_uu.yy * tmp4086
		- gamma_uu.zz * K_ll.zz * K_ll.zz
		+ gamma_uu.xx * tmp2817
		- 2. * gamma_uu.yy * tmp3431
		+ 2. * gamma_uu.xx * tmp2814
		+ 2. * gamma_uu.xx * gamma_uu.yy * tmp1477
		+ 2. * gamma_uu.xx * tmp2811
		+ tmp4657
		- 2. * gamma_uu.xx * tmp1548
		+ d_lll.z.yy * tmp4598
		+ 2. * d_lll.z.yy * tmp2578
		+ 2. * d_lll.z.xz * tmp1231
		+ 2. * tmp4644
		+ 2. * d_lll.z.xy * tmp4458
		+ 4. * d_lll.z.xy * tmp2700
		+ 4. * d_lll.z.xy * tmp2487
		+ 4. * d_lll.z.xy * tmp2482
		+ d_lll.z.xx * tmp4395
		+ 2. * d_lll.z.xx * tmp4390
		+ 2. * d_lll.z.xx * tmp2445
		+ 2. * d_lll.z.xx * tmp2440
		+ 4. * d_lll.z.xx * tmp2436
		+ d_lll.y.zz * tmp1281
		+ 2. * d_lll.y.zz * tmp1436
		+ 2. * d_lll.y.zz * tmp2060
		+ 2. * d_lll.y.yz * tmp3968
		- 4. * d_lll.y.yz * d_lll.z.xz * tmp122
		- 4. * d_lll.y.yz * tmp2578
		- 2. * d_lll.y.yz * tmp4598
		- d_lll.y.zz * tmp1426
		+ d_lll.y.yy * d_lll.y.zz * tmp218
		- 2. * d_lll.y.yy * d_lll.z.xz * tmp110
		- 2. * d_lll.y.yy * d_lll.z.yz * tmp218
		- d_lll.y.yy * d_lll.z.zz * tmp224
		+ 2. * d_lll.y.xz * tmp2474
		- 4. * d_lll.y.xz * tmp2482
		- 4. * d_lll.y.xz * tmp2700
		- 2. * d_lll.y.xz * tmp4458
		- 2. * tmp4565
		+ 2. * d_lll.y.xy * tmp4407
		- 4. * d_lll.y.xy * tmp4409
		- 4. * d_lll.y.xy * tmp4417
		- 2. * d_lll.y.xy * tmp4427
		+ d_lll.y.xx * tmp4374
		- 2. * d_lll.y.xx * d_lll.z.zz * tmp61
		+ 2. * d_lll.y.xx * tmp4369
		- 4. * d_lll.y.xx * d_lll.z.yz * tmp55
		+ 2. * d_lll.y.xx * d_lll.y.zz * tmp55
		- 2. * d_lll.y.xx * tmp4364
		+ 2. * d_lll.x.zz * tmp1228
		- d_lll.x.zz * tmp1231
		- d_lll.y.xx * tmp4361
		+ 2. * d_lll.x.zz * tmp1222
		+ d_lll.x.zz * tmp2038
		+ 2. * d_lll.x.zz * tmp1215
		+ 2. * d_lll.x.zz * tmp1212
		+ d_lll.x.zz * tmp1210
		+ 4. * d_lll.x.zz * tmp1816
		+ 2. * d_lll.x.zz * tmp1204
		- 2. * d_lll.x.zz * tmp1207
		+ d_lll.x.zz * tmp3166
		+ 4. * d_lll.x.yz * tmp2469
		- 2. * d_lll.x.yz * tmp2474
		- 4. * d_lll.x.yz * tmp2482
		- 4. * d_lll.x.yz * tmp2700
		- 2. * d_lll.x.yz * tmp4458
		- 2. * tmp4464
		+ 4. * d_lll.x.yz * tmp2464
		+ 2. * d_lll.x.yz * tmp3755
		- 4. * d_lll.x.yz * tmp2671
		+ d_lll.x.yy * tmp4427
		+ 2. * d_lll.x.yy * tmp4409
		- 4. * d_lll.x.yy * d_lll.z.xz * tmp55
		- 2. * d_lll.x.yy * tmp4417
		- 2. * d_lll.x.yy * d_lll.z.zz * tmp116
		+ d_lll.x.yy * tmp4407
		+ 2. * d_lll.x.yy * d_lll.x.zz * tmp55
		+ 2. * d_lll.x.xz * tmp3715
		- 4. * d_lll.x.xz * tmp2440
		- 4. * d_lll.x.xz * tmp4390
		- 2. * d_lll.x.xz * tmp4395
		- d_lll.x.yy * d_lll.x.zz * tmp86
		+ 2. * d_lll.x.xz * tmp2430
		+ 2. * d_lll.x.xy * tmp4361
		- 4. * d_lll.x.xy * tmp4364
		- 4. * d_lll.x.xy * tmp4369
		- 2. * d_lll.x.xy * tmp4374
		+ 2. * d_lll.x.xy * d_lll.x.zz * tmp43
		+ d_lll.x.xx * d_lll.y.zz * tmp43
		- 2. * d_lll.x.xx * d_lll.z.xz * tmp37
		- 2. * d_lll.x.xx * d_lll.z.yz * tmp43
		- d_lll.x.xx * d_lll.z.zz * tmp49
		+ d_lll.x.xx * d_lll.x.zz * tmp37
		+ a_l.z * tmp4300
		- a_l.z * a_l.z
		+ 2. * a_l.z * tmp4296
		+ 2. * a_l.z * tmp4292
		+ a_l.y * tmp4284
		- a_l.z * tmp2354
		- a_l.z * tmp3635
		+ 2. * a_l.y * tmp4280
		+ 2. * a_l.y * tmp4276
		+ a_l.x * tmp4268
		- a_l.y * tmp4272
		- a_l.y * tmp4274
		+ 2. * a_l.x * tmp4264
		+ 2. * a_l.x * tmp4260
		+ 2. * Z_l.z * tmp3635
		- 4. * Z_l.z * tmp4292
		- 4. * Z_l.z * tmp4296
		- 2. * Z_l.z * tmp4300
		- a_l.x * tmp4256
		- a_l.x * tmp4258
		+ 2. * Z_l.z * tmp2354
		+ 2. * Z_l.y * tmp4274
		- 4. * Z_l.y * tmp4276
		- 4. * Z_l.y * tmp4280
		- 2. * Z_l.y * tmp4284
		+ 2. * Z_l.y * tmp4272
		+ 2. * Z_l.x * tmp4258
		- 4. * Z_l.x * tmp4260
		- 4. * Z_l.x * tmp4264
		- 2. * Z_l.x * tmp4268
		+ 2. * Z_l.x * tmp4256
		+ 4. * M_PI * S * gamma_ll.zz
		- 8. * M_PI * S_ll.zz
		- 4. * M_PI * gamma_ll.zz * rho
		+ K_ll.yy * K_ll.zz * gamma_uu.yy
		- 2. * K_ll.yz * tmp2338
		- 2. * K_ll.zz * Theta
		+ 2. * K_ll.xy * K_ll.zz * gamma_uu.xy
		- 4. * K_ll.xz * tmp2332
		- 2. * K_ll.xz * tmp2334
		+ K_ll.xx * K_ll.zz * gamma_uu.xx);
	deriv->Theta += alpha * (gamma_uu.zz * tmp4657
		+ gamma_uu.zz * tmp4644
		+ gamma_uu.zz * tmp4705
		+ gamma_uu.yy * tmp3390
		- 3. * gamma_uu.yy * tmp3400
		- gamma_uu.yy * tmp3429 * tmp230
		- gamma_uu.yy * tmp3432
		- 3. * gamma_uu.zz * tmp4464
		- gamma_uu.zz * tmp1395 * tmp67
		- 3. * gamma_uu.zz * tmp4565
		- gamma_uu.zz * tmp3386 * tmp230
		+ gamma_uu.yy * tmp3324
		+ gamma_uu.yy * tmp3479
		+ 2. * gamma_uu.xy * tmp2162
		- gamma_uu.yy * tmp1360 * tmp55
		- 3. * gamma_uu.yy * tmp3231
		+ 2. * gamma_uu.xy * gamma_uu.xz * tmp2806
		+ 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yz * tmp1379
		+ 3. * gamma_uu.xx * tmp1545
		- gamma_uu.xx * tmp1549
		+ 3. * gamma_uu.xx * tmp1540
		+ 3. * gamma_uu.xx * tmp1535
		+ gamma_uu.xx * tmp1398
		- gamma_uu.xx * tmp1523 * tmp55
		- 3. * gamma_uu.xx * tmp1479
		- gamma_uu.xx * tmp1527 * tmp67
		- 3. * gamma_uu.xx * tmp1521
		- gamma_uu.xx * tmp1532
		+ gamma_uu.xx * tmp1381
		+ gamma_uu.xx * tmp1363
		+ d_lll.z.yy * tmp6505
		- d_lll.z.yy * tmp6511
		+ 2. * d_lll.z.yy * tmp6495
		- 2. * d_lll.z.yy * tmp6501
		+ 2. * d_lll.z.xz * tmp6050
		- 2. * d_lll.z.xz * tmp6053
		+ 2. * d_lll.z.xy * tmp5946
		- 2. * d_lll.z.xy * tmp5952
		+ 4. * d_lll.z.xy * tmp5936
		- 4. * d_lll.z.xy * tmp5942
		+ 4. * d_lll.z.xy * d_lll.z.yy * tmp5707
		- 2. * d_lll.z.xy * tmp5926
		- 2. * d_lll.z.xy * tmp5932
		+ 4. * d_lll.z.xy * tmp5916
		- 4. * d_lll.z.xy * tmp5922
		+ d_lll.z.xx * tmp5641
		- d_lll.z.xx * tmp5647
		+ 2. * d_lll.z.xx * tmp5631
		- 2. * d_lll.z.xx * tmp5637
		+ 3. * d_lll.z.xx * tmp5625
		+ 2. * d_lll.z.xx * d_lll.z.yy * tmp5442
		+ 2. * d_lll.z.xx * tmp5615
		- 6. * d_lll.z.xx * tmp5621
		+ 2. * d_lll.z.xx * tmp5601
		- 2. * d_lll.z.xx * tmp5607
		- d_lll.z.xx * tmp5611
		+ 4. * d_lll.z.xx * d_lll.z.xy * tmp5296
		- 2. * d_lll.z.xx * tmp5591
		- 2. * d_lll.z.xx * tmp5597
		+ d_lll.y.zz * d_lll.z.yy * tmp6436
		- d_lll.y.zz * d_lll.z.yy * tmp6441
		+ 2. * d_lll.y.zz * d_lll.z.xy * tmp5779
		- 2. * d_lll.y.zz * d_lll.z.xy * tmp5788
		+ d_lll.y.zz * d_lll.z.xx * tmp5494
		- d_lll.y.zz * d_lll.z.xx * tmp5500
		+ 2. * d_lll.y.yz * tmp6511
		+ 4. * d_lll.y.yz * tmp6501
		- 2. * d_lll.y.yz * tmp6505
		+ 4. * d_lll.y.yz * d_lll.z.xz * tmp5788
		- 4. * d_lll.y.yz * tmp6495
		+ 4. * d_lll.y.yz * d_lll.z.xy * tmp5710
		- 4. * d_lll.y.yz * d_lll.z.xy * tmp5775
		- 4. * d_lll.y.yz * d_lll.z.xz * tmp5785
		+ 4. * d_lll.y.yz * d_lll.z.xx * tmp5438
		- 4. * d_lll.y.yz * d_lll.z.xx * tmp5442
		- 2. * d_lll.y.yz * d_lll.z.xx * tmp5570
		+ 2. * d_lll.y.yz * d_lll.z.xx * tmp5428
		+ 2. * d_lll.y.yz * d_lll.y.zz * tmp6436
		- 2. * d_lll.y.yz * d_lll.y.zz * tmp6441
		+ d_lll.y.yy * d_lll.z.zz * tmp6441
		+ 2. * d_lll.y.yy * d_lll.z.yz * tmp6386
		- 2. * d_lll.y.yy * d_lll.z.yz * tmp6391
		- d_lll.y.yy * d_lll.z.zz * tmp6436
		+ 2. * d_lll.y.yy * d_lll.z.xz * tmp5710
		+ 2. * d_lll.y.yy * d_lll.z.xy * tmp5753
		- 2. * d_lll.y.yy * d_lll.z.xy * tmp5759
		- 2. * d_lll.y.yy * d_lll.z.xz * tmp5707
		+ d_lll.y.yy * d_lll.z.xx * tmp5424
		+ d_lll.y.yy * d_lll.z.xx * tmp5418
		- 2. * d_lll.y.yy * d_lll.z.xx * tmp5474
		+ d_lll.y.yy * d_lll.y.zz * tmp6391
		+ 2. * d_lll.y.xz * tmp5952
		- d_lll.y.yy * d_lll.y.zz * tmp6386
		+ 4. * d_lll.y.xz * tmp5942
		- 2. * d_lll.y.xz * tmp5946
		+ 2. * d_lll.y.xz * tmp5932
		- 4. * d_lll.y.xz * tmp5936
		+ 4. * d_lll.y.xz * tmp5922
		- 2. * d_lll.y.xz * tmp5926
		+ 2. * d_lll.y.xz * tmp5912
		- 4. * d_lll.y.xz * tmp5916
		+ 2. * d_lll.y.xz * tmp5908
		+ 2. * d_lll.y.xz * tmp5898
		- 4. * d_lll.y.xz * tmp5902
		+ 2. * d_lll.y.xz * tmp5882
		- 2. * d_lll.y.xz * tmp5886
		- 2. * d_lll.y.xz * tmp5892
		+ 4. * d_lll.y.xz * d_lll.y.zz * tmp5785
		- 2. * d_lll.y.xz * tmp5876
		+ 4. * d_lll.y.xz * tmp5866
		- 2. * d_lll.y.xz * tmp5872
		+ 2. * d_lll.y.xz * tmp5856
		- 4. * d_lll.y.xz * tmp5862
		+ 2. * d_lll.y.xy * tmp5789
		- 2. * d_lll.y.xz * tmp5852
		+ 4. * d_lll.y.xy * tmp5770
		- 2. * d_lll.y.xy * tmp5786
		+ 4. * d_lll.y.xy * tmp5738
		- 4. * d_lll.y.xy * tmp5764
		+ 4. * d_lll.y.xy * tmp5728
		- 4. * d_lll.y.xy * tmp5734
		+ 4. * d_lll.y.xy * d_lll.z.xx * tmp5328
		- 2. * d_lll.y.xy * tmp5716
		- 2. * d_lll.y.xy * tmp5719
		- 4. * d_lll.y.xy * tmp5724
		+ 2. * d_lll.y.xy * tmp5708
		- 2. * d_lll.y.xy * tmp5711
		+ 4. * d_lll.y.xy * tmp5698
		- 4. * d_lll.y.xy * tmp5702
		+ d_lll.y.xx * tmp5501
		+ d_lll.y.xx * tmp5495
		- 2. * d_lll.y.xx * d_lll.z.zz * tmp5580
		+ 4. * d_lll.y.xx * d_lll.z.yz * tmp5438
		- 2. * d_lll.y.xx * tmp5491
		- 4. * d_lll.y.xx * d_lll.z.yz * tmp5570
		+ 2. * d_lll.y.xx * tmp5485
		+ d_lll.y.xx * tmp5469
		- d_lll.y.xx * tmp5479
		+ 4. * d_lll.y.xx * d_lll.z.xz * tmp5299
		- 2. * d_lll.y.xx * tmp5465
		+ 2. * d_lll.y.xx * tmp5449
		- 2. * d_lll.y.xx * tmp5455
		- 2. * d_lll.y.xx * tmp5459
		+ 2. * d_lll.y.xx * d_lll.z.xx * tmp5246
		+ 2. * d_lll.y.xx * d_lll.y.zz * tmp5570
		- 2. * d_lll.y.xx * d_lll.z.xx * tmp5240
		+ 3. * d_lll.y.xx * tmp5443
		+ 2. * d_lll.y.xx * tmp5433
		- 6. * d_lll.y.xx * tmp5439
		+ 2. * d_lll.y.xx * tmp5419
		- 2. * d_lll.y.xx * tmp5425
		- d_lll.y.xx * tmp5429
		+ d_lll.y.xx * tmp5409
		- d_lll.y.xx * tmp5415
		+ 4. * d_lll.y.xx * d_lll.y.xz * tmp5286
		- 2. * d_lll.y.xx * tmp5405
		+ 2. * d_lll.y.xx * tmp5389
		- 2. * d_lll.y.xx * tmp5395
		- 2. * d_lll.y.xx * tmp5399
		+ d_lll.x.zz * tmp6050
		- d_lll.x.zz * tmp6053
		+ 2. * d_lll.x.zz * d_lll.z.xy * tmp5580
		- 2. * d_lll.x.zz * d_lll.z.xy * tmp5500
		+ d_lll.x.zz * d_lll.z.xx * tmp5350
		- d_lll.x.zz * d_lll.z.xx * tmp5355
		+ 2. * d_lll.x.zz * d_lll.y.zz * tmp5945
		- 2. * d_lll.x.zz * d_lll.y.zz * tmp5951
		+ 4. * d_lll.x.zz * d_lll.y.yz * tmp5779
		- 2. * d_lll.x.zz * d_lll.y.yz * tmp5785
		- 2. * d_lll.x.zz * d_lll.y.yz * tmp5788
		+ d_lll.x.zz * d_lll.y.yy * tmp5710
		- 2. * d_lll.x.zz * d_lll.y.yy * tmp5775
		+ d_lll.x.zz * d_lll.y.yy * tmp5707
		+ 2. * d_lll.x.zz * d_lll.y.xz * tmp5580
		- 2. * d_lll.x.zz * d_lll.y.xz * tmp5500
		+ 4. * d_lll.x.zz * d_lll.y.xy * tmp5438
		- 4. * d_lll.x.zz * d_lll.y.xy * tmp5442
		+ 2. * d_lll.x.zz * d_lll.y.xy * tmp5428
		- 2. * d_lll.x.zz * d_lll.y.xy * tmp5432
		+ d_lll.x.zz * d_lll.y.xx * tmp5296
		- d_lll.x.zz * d_lll.y.xx * tmp5305
		+ 2. * d_lll.x.yz * tmp5952
		+ 4. * d_lll.x.yz * tmp5942
		- 2. * d_lll.x.yz * tmp5946
		+ 2. * d_lll.x.yz * tmp5932
		- 4. * d_lll.x.yz * tmp5936
		+ 4. * d_lll.x.yz * tmp5922
		- 2. * d_lll.x.yz * tmp5926
		+ 2. * d_lll.x.yz * tmp5912
		- 4. * d_lll.x.yz * tmp5916
		+ 2. * d_lll.x.yz * tmp5908
		+ 2. * d_lll.x.yz * tmp5898
		- 4. * d_lll.x.yz * tmp5902
		+ 2. * d_lll.x.yz * tmp5882
		- 2. * d_lll.x.yz * tmp5886
		- 2. * d_lll.x.yz * tmp5892
		+ 2. * d_lll.x.yz * tmp5872
		- 2. * d_lll.x.yz * tmp5876
		+ 4. * d_lll.x.yz * tmp5862
		- 4. * d_lll.x.yz * tmp5866
		+ 2. * d_lll.x.yz * tmp5852
		- 2. * d_lll.x.yz * tmp5856
		+ 2. * d_lll.x.yz * d_lll.y.xz * tmp5570
		+ 2. * d_lll.x.yz * d_lll.y.xz * tmp5442
		+ 2. * d_lll.x.yz * d_lll.y.xz * tmp5432
		- 4. * d_lll.x.yz * d_lll.y.xz * tmp5438
		+ 4. * d_lll.x.yz * d_lll.y.xy * tmp5424
		- 2. * d_lll.x.yz * d_lll.y.xz * tmp5428
		+ 2. * d_lll.x.yz * d_lll.y.xx * tmp5328
		- 2. * d_lll.x.yz * d_lll.y.xx * tmp5292
		- 4. * d_lll.x.yz * d_lll.y.xy * tmp5474
		+ 4. * d_lll.x.yz * d_lll.x.zz * tmp5494
		- 2. * d_lll.x.yz * d_lll.x.zz * tmp5580
		- 2. * d_lll.x.yz * d_lll.x.zz * tmp5500
		+ d_lll.x.yy * tmp5789
		+ d_lll.x.yy * tmp5786
		+ 4. * d_lll.x.yy * d_lll.z.yz * tmp5775
		- 2. * d_lll.x.yy * d_lll.z.zz * tmp5779
		+ 2. * d_lll.x.yy * d_lll.z.yy * tmp5759
		- 2. * d_lll.x.yy * tmp5764
		- 2. * d_lll.x.yy * tmp5770
		+ 4. * d_lll.x.yy * d_lll.z.xz * tmp5438
		- 4. * d_lll.x.yy * d_lll.z.xz * tmp5570
		- 2. * d_lll.x.yy * d_lll.z.yy * tmp5753
		+ 2. * d_lll.x.yy * tmp5734
		- 2. * d_lll.x.yy * tmp5738
		+ 2. * d_lll.x.yy * tmp5724
		- 2. * d_lll.x.yy * tmp5728
		+ d_lll.x.yy * tmp5716
		- d_lll.x.yy * tmp5719
		+ d_lll.x.yy * tmp5708
		- d_lll.x.yy * tmp5711
		+ 2. * d_lll.x.yy * tmp5698
		- 2. * d_lll.x.yy * tmp5702
		+ d_lll.x.yy * d_lll.y.xx * tmp5278
		- d_lll.x.yy * d_lll.y.xx * tmp5283
		+ 2. * d_lll.x.yy * d_lll.x.zz * tmp5570
		+ 2. * d_lll.x.yy * d_lll.x.zz * tmp5442
		+ 3. * d_lll.x.yy * d_lll.x.zz * tmp5432
		- 6. * d_lll.x.yy * d_lll.x.zz * tmp5438
		+ 4. * d_lll.x.yy * d_lll.x.yz * tmp5418
		- 2. * d_lll.x.yy * d_lll.x.yz * tmp5474
		- 2. * d_lll.x.yy * d_lll.x.yz * tmp5424
		- d_lll.x.yy * d_lll.x.zz * tmp5428
		+ 2. * d_lll.x.xz * tmp5647
		+ 4. * d_lll.x.xz * tmp5637
		- 2. * d_lll.x.xz * tmp5641
		+ 4. * d_lll.x.xz * tmp5621
		- 2. * d_lll.x.xz * tmp5625
		- 4. * d_lll.x.xz * tmp5631
		+ 2. * d_lll.x.xz * tmp5611
		- 4. * d_lll.x.xz * tmp5615
		+ 4. * d_lll.x.xz * tmp5607
		+ 4. * d_lll.x.xz * tmp5597
		- 4. * d_lll.x.xz * tmp5601
		+ 4. * d_lll.x.xz * d_lll.y.zz * tmp5580
		- 2. * d_lll.x.xz * d_lll.y.zz * tmp5500
		- 4. * d_lll.x.xz * tmp5591
		+ 4. * d_lll.x.xz * d_lll.y.yz * tmp5570
		- 2. * d_lll.x.xz * d_lll.y.zz * tmp5494
		+ 2. * d_lll.x.xz * d_lll.y.yy * tmp5424
		- 4. * d_lll.x.xz * d_lll.y.yz * tmp5428
		+ 4. * d_lll.x.xz * d_lll.y.xz * tmp5305
		- 2. * d_lll.x.xz * d_lll.y.yy * tmp5418
		+ 4. * d_lll.x.xz * d_lll.y.xy * tmp5292
		- 4. * d_lll.x.xz * d_lll.y.xz * tmp5299
		+ 2. * d_lll.x.xz * d_lll.x.zz * tmp5350
		- 2. * d_lll.x.xz * d_lll.x.zz * tmp5355
		- 4. * d_lll.x.xz * d_lll.y.xy * tmp5286
		+ 4. * d_lll.x.xz * d_lll.x.yz * tmp5299
		- 4. * d_lll.x.xz * d_lll.x.yz * tmp5305
		+ 2. * d_lll.x.xz * d_lll.x.yy * tmp5286
		- 2. * d_lll.x.xz * d_lll.x.yy * tmp5292
		+ 2. * d_lll.x.xy * tmp5501
		+ 4. * d_lll.x.xy * tmp5491
		- 2. * d_lll.x.xy * tmp5495
		+ 4. * d_lll.x.xy * d_lll.z.yy * tmp5474
		- 2. * d_lll.x.xy * tmp5479
		- 4. * d_lll.x.xy * tmp5485
		+ 4. * d_lll.x.xy * tmp5465
		- 2. * d_lll.x.xy * tmp5469
		+ 4. * d_lll.x.xy * tmp5455
		- 4. * d_lll.x.xy * tmp5459
		+ 4. * d_lll.x.xy * tmp5439
		- 2. * d_lll.x.xy * tmp5443
		- 4. * d_lll.x.xy * tmp5449
		+ 2. * d_lll.x.xy * tmp5429
		- 4. * d_lll.x.xy * tmp5433
		+ 4. * d_lll.x.xy * tmp5425
		+ 2. * d_lll.x.xy * tmp5415
		- 4. * d_lll.x.xy * tmp5419
		+ 4. * d_lll.x.xy * tmp5405
		- 2. * d_lll.x.xy * tmp5409
		+ 4. * d_lll.x.xy * tmp5395
		- 4. * d_lll.x.xy * tmp5399
		+ 2. * d_lll.x.xy * d_lll.x.zz * tmp5296
		- 2. * d_lll.x.xy * d_lll.x.zz * tmp5305
		- 4. * d_lll.x.xy * tmp5389
		+ 4. * d_lll.x.xy * d_lll.x.yz * tmp5328
		- 4. * d_lll.x.xy * d_lll.x.yz * tmp5292
		+ 2. * d_lll.x.xy * d_lll.x.yy * tmp5278
		- 2. * d_lll.x.xy * d_lll.x.yy * tmp5283
		+ d_lll.x.xx * d_lll.z.zz * tmp5355
		+ 2. * d_lll.x.xx * d_lll.z.yz * tmp5305
		- d_lll.x.xx * d_lll.z.zz * tmp5350
		+ d_lll.x.xx * d_lll.z.yy * tmp5292
		- 2. * d_lll.x.xx * d_lll.z.yz * tmp5296
		+ d_lll.x.xx * d_lll.z.yy * tmp5286
		+ 2. * d_lll.x.xx * d_lll.z.xz * tmp5250
		- 2. * d_lll.x.xx * d_lll.z.xz * tmp5255
		- 2. * d_lll.x.xx * d_lll.z.yy * tmp5328
		+ 2. * d_lll.x.xx * d_lll.z.xy * tmp5240
		- 2. * d_lll.x.xx * d_lll.z.xy * tmp5246
		+ d_lll.x.xx * d_lll.y.zz * tmp5305
		+ d_lll.x.xx * d_lll.y.zz * tmp5296
		- 2. * d_lll.x.xx * d_lll.y.zz * tmp5299
		+ 2. * d_lll.x.xx * d_lll.y.yz * tmp5292
		+ d_lll.x.xx * d_lll.y.yy * tmp5283
		- 2. * d_lll.x.xx * d_lll.y.yz * tmp5286
		+ 2. * d_lll.x.xx * d_lll.y.xz * tmp5240
		- 2. * d_lll.x.xx * d_lll.y.xz * tmp5246
		- d_lll.x.xx * d_lll.y.yy * tmp5278
		+ 2. * d_lll.x.xx * d_lll.y.xy * tmp5232
		- 2. * d_lll.x.xx * d_lll.y.xy * tmp5237
		+ d_lll.x.xx * d_lll.x.zz * tmp5255
		+ 2. * d_lll.x.xx * d_lll.x.yz * tmp5246
		- d_lll.x.xx * d_lll.x.zz * tmp5250
		+ d_lll.x.xx * d_lll.x.yy * tmp5237
		- 2. * d_lll.x.xx * d_lll.x.yz * tmp5240
		+ tmp3431 * tmp230
		- d_lll.x.xx * d_lll.x.yy * tmp5232
		+ tmp1548 * tmp67
		+ tmp1531 * tmp55
		+ Z_l.z * tmp1326
		- 2. * Z_l.z * tmp1331
		- 2. * Z_l.z * tmp1334
		- Z_l.z * tmp1337
		+ 2. * Z_l.z * tmp1494
		- 4. * Z_l.z * tmp1320
		- 2. * Z_l.z * tmp1323
		+ Z_l.z * tmp1961
		- 2. * Z_l.z * tmp2079
		+ Z_l.z * tmp1936
		- 2. * Z_l.z * tmp1940
		- 2. * Z_l.z * tmp1305
		- 2. * Z_l.z * tmp2647
		- Z_l.z * tmp1311
		- 2. * Z_l.z * tmp1314
		- Z_l.z * tmp1317
		+ Z_l.z * tmp1290
		- 2. * Z_l.z * tmp5151
		- Z_l.z * tmp1300
		+ Z_l.y * tmp1426
		- 2. * Z_l.y * tmp2060
		- 2. * Z_l.y * tmp1436
		- 2. * Z_l.y * tmp1278
		- Z_l.y * tmp1281
		- 2. * Z_l.y * tmp1284
		- Z_l.y * tmp1287
		- Z_l.z * a_l.x * gamma_uu.xz
		- Z_l.z * a_l.y * gamma_uu.yz
		- Z_l.z * a_l.z * gamma_uu.zz
		- Z_l.z * tmp5130
		- 2. * Z_l.z * tmp5134
		- 2. * Z_l.z * tmp5139
		- 2. * Z_l.z * tmp5144
		+ Z_l.y * tmp1267
		- 2. * Z_l.y * tmp1272
		+ 2. * Z_l.y * tmp1405
		- Z_l.y * tmp1261
		- 2. * Z_l.y * tmp1264
		+ Z_l.y * tmp2592
		- 2. * Z_l.y * tmp1905
		- 2. * Z_l.y * tmp1255
		- 4. * Z_l.y * tmp1258
		+ Z_l.y * tmp1250
		- 2. * Z_l.y * tmp1902
		+ Z_l.x * tmp1231
		- 2. * Z_l.x * tmp1235
		- Z_l.x * tmp1238
		- Z_l.y * a_l.x * gamma_uu.xy
		- Z_l.y * a_l.y * gamma_uu.yy
		- Z_l.y * a_l.z * gamma_uu.yz
		- Z_l.y * tmp5025
		- 2. * Z_l.y * tmp5029
		- 2. * Z_l.y * tmp5034
		- Z_l.y * tmp1240
		- 2. * Z_l.y * tmp1896
		+ Z_l.x * tmp1215
		- 2. * Z_l.x * tmp1219
		- Z_l.x * tmp2038
		- 2. * Z_l.x * tmp1222
		- 2. * Z_l.x * tmp1225
		- 2. * Z_l.x * tmp1228
		+ Z_l.x * tmp1200
		- 2. * Z_l.x * tmp2505
		- Z_l.x * tmp3166
		- 2. * Z_l.x * tmp1204
		- 2. * Z_l.x * tmp1207
		- Z_l.x * tmp1210
		- 2. * Z_l.x * tmp1212
		+ 2. * Z_l.x * tmp1195
		- 4. * Z_l.x * tmp1808
		+ Z_l.x * tmp1191
		- 2. * Z_l.x * tmp1805
		+ K_ll.yy * K_ll.zz * tmp279
		- K_ll.yy * K_ll.zz * tmp230
		- K_ll.yy * Theta * gamma_uu.yy
		- 2. * K_ll.yz * Theta * gamma_uu.yz
		- K_ll.zz * Theta * gamma_uu.zz
		- 8. * M_PI * rho
		- Z_l.x * a_l.x * gamma_uu.xx
		- Z_l.x * a_l.y * gamma_uu.xy
		- Z_l.x * a_l.z * gamma_uu.xz
		- Z_l.x * tmp4920
		- 2. * Z_l.x * tmp4924
		- 2. * Z_l.x * tmp4929
		+ 2. * K_ll.xz * K_ll.yz * tmp128
		- 2. * K_ll.xz * Theta * gamma_uu.xz
		+ 2. * K_ll.xz * K_ll.yy * tmp122
		- 2. * K_ll.xz * K_ll.yz * tmp177
		+ 2. * K_ll.xy * K_ll.zz * tmp177
		- 2. * K_ll.xy * K_ll.zz * tmp128
		- 2. * K_ll.xy * Theta * gamma_uu.xy
		- 2. * K_ll.xz * K_ll.yy * tmp116
		+ 2. * K_ll.xy * K_ll.yz * tmp116
		- 2. * K_ll.xy * K_ll.yz * tmp122
		+ 2. * K_ll.xy * K_ll.xz * tmp61
		+ K_ll.xx * K_ll.zz * tmp159
		- K_ll.xx * K_ll.zz * tmp67
		- K_ll.xx * Theta * gamma_uu.xx
		- 2. * K_ll.xy * K_ll.xz * tmp98
		+ 2. * K_ll.xx * K_ll.yz * tmp98
		- 2. * K_ll.xx * K_ll.yz * tmp61
		+ K_ll.xx * K_ll.yy * tmp86
		- K_ll.xx * K_ll.yy * tmp55);
	deriv->Z_l.x += -alpha * (Theta * a_l.x
		+ 8. * M_PI * S_l.x
		+ K_ll.xz * tmp1337
		- K_ll.yy * d_lll.x.xx * tmp55
		- 2. * K_ll.yy * d_lll.x.xy * tmp110
		- 2. * K_ll.yy * d_lll.x.xz * tmp116
		- K_ll.yy * d_lll.x.yy * tmp218
		- 2. * K_ll.yy * tmp1340
		- K_ll.yy * tmp1345
		- 2. * K_ll.yz * d_lll.x.xx * tmp61
		- 2. * K_ll.yz * d_lll.x.xy * tmp116
		- 2. * K_ll.yz * d_lll.x.xy * tmp122
		- 2. * K_ll.yz * d_lll.x.xz * tmp177
		- 2. * K_ll.yz * d_lll.x.xz * tmp128
		- 2. * K_ll.yz * d_lll.x.yy * tmp224
		- 2. * K_ll.yz * d_lll.x.yz * tmp279
		- 2. * K_ll.yz * d_lll.x.yz * tmp230
		- 2. * K_ll.yz * tmp1365
		- K_ll.zz * d_lll.x.xx * tmp67
		- 2. * K_ll.zz * d_lll.x.xy * tmp128
		- 2. * K_ll.zz * d_lll.x.xz * tmp189
		- K_ll.zz * d_lll.x.yy * tmp230
		- 2. * K_ll.zz * d_lll.x.yz * tmp291
		- K_ll.zz * d_lll.x.zz * tmp332
		+ 2. * K_ll.xz * tmp1334
		+ 2. * K_ll.xz * tmp1331
		+ 2. * K_ll.xz * tmp1323
		- K_ll.xz * tmp1326
		+ 4. * K_ll.xz * tmp1320
		+ 2. * K_ll.xz * tmp2079
		- 2. * K_ll.xz * tmp1494
		+ K_ll.xz * tmp1317
		- K_ll.xz * tmp1961
		+ 2. * K_ll.xz * tmp1314
		+ K_ll.xz * tmp1311
		+ 2. * K_ll.xz * tmp2647
		+ 2. * K_ll.xz * tmp1305
		+ 2. * K_ll.xz * tmp1940
		+ 2. * K_ll.xz * tmp7345
		- K_ll.xz * tmp5130
		- 2. * K_ll.xz * d_lll.x.xy * tmp61
		- 2. * K_ll.xz * d_lll.x.xz * tmp67
		- K_ll.xz * tmp1290
		- 2. * K_ll.xz * tmp1295
		- K_ll.xz * tmp1300
		- K_ll.xz * tmp1936
		+ 2. * K_ll.xz * tmp7343
		+ 2. * K_ll.xz * tmp7341
		+ K_ll.xy * tmp1287
		+ 2. * K_ll.xy * tmp1284
		+ K_ll.xy * tmp1281
		+ 2. * K_ll.xy * tmp1278
		+ 2. * K_ll.xy * tmp1436
		+ 2. * K_ll.xy * tmp2060
		+ 2. * K_ll.xy * tmp1272
		- K_ll.xy * tmp1426
		+ 2. * K_ll.xy * tmp1264
		- K_ll.xy * tmp1267
		+ K_ll.xy * tmp1261
		+ 4. * K_ll.xy * tmp1258
		- 2. * K_ll.xy * tmp1405
		+ 2. * K_ll.xy * tmp1255
		+ 2. * K_ll.xy * tmp1905
		+ 2. * K_ll.xy * tmp7262
		- K_ll.xy * tmp5025
		- 2. * K_ll.xy * d_lll.x.xy * tmp55
		- 2. * K_ll.xy * d_lll.x.xz * tmp61
		- K_ll.xy * tmp1240
		- 2. * K_ll.xy * tmp1245
		- K_ll.xy * tmp1250
		- K_ll.xy * tmp2592
		+ 2. * K_ll.xy * tmp7260
		+ 2. * K_ll.xy * tmp7258
		+ K_ll.xx * tmp1238
		+ 2. * K_ll.xx * tmp1235
		+ 2. * K_ll.xx * tmp1228
		- K_ll.xx * tmp1231
		+ 2. * K_ll.xx * tmp1225
		+ 2. * K_ll.xx * tmp1222
		+ K_ll.xx * tmp2038
		+ 2. * K_ll.xx * tmp1219
		+ 2. * K_ll.xx * tmp1212
		- K_ll.xx * tmp1215
		+ K_ll.xx * tmp1210
		+ 2. * K_ll.xx * tmp1207
		+ 2. * K_ll.xx * tmp1204
		+ K_ll.xx * tmp3166
		+ K_ll.xx * tmp2505
		+ 2. * K_ll.xx * tmp1808
		- K_ll.xx * tmp1200
		+ K_ll.xx * tmp1805
		- 2. * K_ll.xx * tmp1195
		+ 2. * K_ll.xx * tmp7196
		- K_ll.xx * tmp1191
		+ 2. * K_ll.xx * tmp7194
		+ 2. * K_ll.xx * tmp7192);
	deriv->Z_l.y += alpha * (K_ll.zz * tmp2010
		- 8. * M_PI * S_l.y
		- Theta * a_l.y
		+ 2. * K_ll.zz * tmp2005
		+ K_ll.zz * tmp2001
		+ 2. * K_ll.zz * tmp1997
		+ 2. * K_ll.zz * tmp1992
		+ K_ll.zz * tmp1990
		+ K_ll.yz * tmp1326
		- 2. * K_ll.yz * tmp1331
		- 2. * K_ll.yz * tmp1334
		- K_ll.yz * tmp1337
		+ 2. * K_ll.yz * tmp1494
		- 4. * K_ll.yz * tmp1320
		- 2. * K_ll.yz * tmp1323
		+ K_ll.yz * tmp1961
		- 2. * K_ll.yz * tmp2079
		+ K_ll.yz * tmp1317
		+ 2. * K_ll.yz * tmp1952
		+ K_ll.yz * tmp1311
		+ 2. * K_ll.yz * tmp1308
		+ 2. * K_ll.yz * tmp1943
		+ K_ll.yz * tmp1936
		+ K_ll.yz * tmp1290
		- 2. * K_ll.yz * tmp5151
		- K_ll.yz * tmp1300
		+ K_ll.yy * tmp1426
		- 2. * K_ll.yy * tmp2060
		- 2. * K_ll.yy * tmp1436
		- 2. * K_ll.yy * tmp1278
		- K_ll.yy * tmp1281
		- 2. * K_ll.yy * tmp1284
		- K_ll.yy * tmp1287
		- 2. * K_ll.yz * tmp7341
		- 2. * K_ll.yz * tmp7343
		- 2. * K_ll.yz * tmp7345
		- K_ll.yz * tmp5130
		- 2. * K_ll.yz * tmp5134
		- 2. * K_ll.yz * tmp5139
		- 2. * K_ll.yz * tmp5144
		+ K_ll.yy * tmp1267
		- K_ll.yy * tmp1272
		+ 2. * K_ll.yy * tmp1405
		+ K_ll.yy * tmp2592
		- K_ll.yy * tmp1905
		- 2. * K_ll.yy * tmp1258
		+ K_ll.yy * tmp1250
		- 2. * K_ll.yy * tmp1902
		+ 2. * K_ll.xz * tmp1875
		- 2. * K_ll.yy * tmp7258
		- 2. * K_ll.yy * tmp7260
		- 2. * K_ll.yy * tmp7262
		- K_ll.yy * tmp5025
		- 2. * K_ll.yy * tmp5029
		- 2. * K_ll.yy * tmp5034
		- K_ll.yy * tmp1240
		- 2. * K_ll.yy * tmp1896
		+ 2. * K_ll.xz * tmp3209
		+ 2. * K_ll.xz * tmp1870
		+ 2. * K_ll.xz * tmp1865
		+ 2. * K_ll.xz * tmp1860
		+ 2. * K_ll.xz * tmp3852
		+ 2. * K_ll.xz * tmp1855
		+ 2. * K_ll.xz * tmp3844
		+ 2. * K_ll.xz * d_lll.y.xx * tmp49
		+ K_ll.xy * tmp1231
		- 2. * K_ll.xy * tmp1235
		- K_ll.xy * tmp1238
		+ K_ll.xy * tmp1215
		- K_ll.xy * tmp2038
		- 2. * K_ll.xy * tmp1222
		- 2. * K_ll.xy * tmp1225
		- 2. * K_ll.xy * tmp1228
		+ 2. * K_ll.xy * tmp1826
		+ K_ll.xy * tmp1210
		+ 2. * K_ll.xy * tmp1816
		+ 2. * K_ll.xy * tmp1811
		+ K_ll.xy * tmp3166
		+ K_ll.xy * tmp1200
		- 2. * K_ll.xy * tmp2505
		+ 2. * K_ll.xy * tmp1195
		- 4. * K_ll.xy * tmp1808
		+ K_ll.xy * tmp1191
		- 2. * K_ll.xy * tmp1805
		+ K_ll.xx * tmp1797
		- 2. * K_ll.xy * tmp7192
		- 2. * K_ll.xy * tmp7194
		- 2. * K_ll.xy * tmp7196
		- K_ll.xy * tmp4920
		- 2. * K_ll.xy * tmp4924
		- 2. * K_ll.xy * tmp4929
		+ 2. * K_ll.xx * tmp1792
		+ K_ll.xx * tmp1788
		+ 2. * K_ll.xx * tmp1784
		+ 2. * K_ll.xx * tmp1779
		+ K_ll.xx * d_lll.y.xx * tmp37);
	deriv->Z_l.z += alpha * (K_ll.zz * tmp1326
		- K_ll.zz * tmp1331
		- 8. * M_PI * S_l.z
		- Theta * a_l.z
		+ 2. * K_ll.zz * tmp1494
		- 2. * K_ll.zz * tmp1320
		+ K_ll.zz * tmp1961
		- K_ll.zz * tmp2079
		+ K_ll.zz * tmp1936
		- 2. * K_ll.zz * tmp1940
		- 2. * K_ll.zz * tmp1305
		- 2. * K_ll.zz * tmp2647
		- K_ll.zz * tmp1311
		- 2. * K_ll.zz * tmp1314
		- K_ll.zz * tmp1317
		+ K_ll.zz * tmp1290
		- 2. * K_ll.zz * tmp5151
		- K_ll.zz * tmp1300
		+ K_ll.yz * tmp1287
		- 2. * K_ll.zz * tmp7341
		- 2. * K_ll.zz * tmp7343
		- 2. * K_ll.zz * tmp7345
		- K_ll.zz * tmp5130
		- 2. * K_ll.zz * tmp5134
		- 2. * K_ll.zz * tmp5139
		- 2. * K_ll.zz * tmp5144
		+ 2. * K_ll.yz * tmp2633
		+ K_ll.yz * tmp1281
		+ 2. * K_ll.yz * tmp2624
		+ 2. * K_ll.yz * tmp1275
		+ K_ll.yz * tmp1426
		+ K_ll.yz * tmp1267
		- 2. * K_ll.yz * tmp1272
		+ 2. * K_ll.yz * tmp1405
		- K_ll.yz * tmp1261
		- 2. * K_ll.yz * tmp1264
		+ K_ll.yz * tmp2592
		- 2. * K_ll.yz * tmp1905
		- 2. * K_ll.yz * tmp1255
		- 4. * K_ll.yz * tmp1258
		+ K_ll.yz * tmp1250
		- 2. * K_ll.yz * tmp1902
		+ K_ll.yy * tmp2583
		- 2. * K_ll.yz * tmp7258
		- 2. * K_ll.yz * tmp7260
		- 2. * K_ll.yz * tmp7262
		- K_ll.yz * tmp5025
		- 2. * K_ll.yz * tmp5029
		- 2. * K_ll.yz * tmp5034
		- K_ll.yz * tmp1240
		- 2. * K_ll.yz * tmp1896
		+ 2. * K_ll.yy * tmp2578
		+ K_ll.yy * tmp2574
		+ 2. * K_ll.yy * tmp2569
		+ 2. * K_ll.yy * tmp2565
		+ K_ll.yy * tmp2563
		+ K_ll.xz * tmp1238
		+ 2. * K_ll.xz * tmp2536
		+ K_ll.xz * tmp1231
		+ 2. * K_ll.xz * tmp2526
		+ 2. * K_ll.xz * tmp1836
		+ K_ll.xz * tmp2038
		+ K_ll.xz * tmp1215
		- 2. * K_ll.xz * tmp1219
		+ K_ll.xz * tmp1200
		- 2. * K_ll.xz * tmp2505
		- K_ll.xz * tmp3166
		- 2. * K_ll.xz * tmp1204
		- 2. * K_ll.xz * tmp1207
		- K_ll.xz * tmp1210
		- 2. * K_ll.xz * tmp1212
		+ 2. * K_ll.xz * tmp1195
		- 4. * K_ll.xz * tmp1808
		+ K_ll.xz * tmp1191
		- 2. * K_ll.xz * tmp1805
		+ 2. * K_ll.xy * tmp2497
		- 2. * K_ll.xz * tmp7192
		- 2. * K_ll.xz * tmp7194
		- 2. * K_ll.xz * tmp7196
		- K_ll.xz * tmp4920
		- 2. * K_ll.xz * tmp4924
		- 2. * K_ll.xz * tmp4929
		+ 2. * K_ll.xy * tmp2492
		+ 2. * K_ll.xy * tmp2700
		+ 2. * K_ll.xy * tmp2487
		+ 2. * K_ll.xy * tmp2482
		+ 2. * K_ll.xy * tmp2692
		+ 2. * K_ll.xy * tmp2477
		+ 2. * K_ll.xy * tmp2685
		+ 2. * K_ll.xy * tmp2683
		+ K_ll.xx * tmp2454
		+ 2. * K_ll.xx * tmp2449
		+ K_ll.xx * tmp2445
		+ 2. * K_ll.xx * tmp2440
		+ 2. * K_ll.xx * tmp2436
		+ K_ll.xx * tmp3890);
	// END CUT

#endif
//decay for 1st deriv hyperbolic state vars constraints:

	//turns out if you "if conv != 0" all these then skipping decay explodes soon in simulation steps, but time grows quickly, so it dies at a high t value 

	// a_x = log(alpha)_,x <=> a_x += eta (log(alpha)_,x - a_x)
	if (solver->a_convCoeff != 0.) 
	{
		<? for i,xi in ipairs(xNames) do ?>{
			<? if i <= solver.dim then ?>
			real di_log_alpha = (
				log(U[solver->stepsize.<?=xi?>].alpha) 
				- log(U[-solver->stepsize.<?=xi?>].alpha)
			) / (2. * solver->grid_dx.<?=xi?>);
			<? else ?>
			real di_log_alpha = 0.;
			<? end ?>
			deriv->a_l.<?=xi?> += solver->a_convCoeff * (di_log_alpha - U->a_l.<?=xi?>);
		}<? end ?>	
	}

	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	if (solver->d_convCoeff != 0.) 
	{
		<? for i,xi in ipairs(xNames) do 
			for jk,xjk in ipairs(symNames) do ?>{
				<? if i <= solver.dim then ?>
			real di_gamma_jk = .5 * (
				U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> 
				- U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>
			) / (2. * solver->grid_dx.<?=xi?>);
				<? else ?>
			real di_gamma_jk = 0;
				<? end ?>
			deriv->d_lll.<?=xi?>.<?=xjk?> += solver->d_convCoeff * (.5 * di_gamma_jk - U->d_lll.<?=xi?>.<?=xjk?>);
		}<? end ?>
		<? end ?>
	}
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

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3x3 d_llu = _3sym3_sym3_mul(U->d_lll, gamma_uu);
	
	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gamma_uu, U->d_lll);

	//e_i = d^j_ji
	real3 e_l = _3sym3_tr12(d_ull);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 conn_ull = {
<? for k,xk in ipairs(xNames) do 
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i],xNames[j]
?>			.<?=xij?> = d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?> - d_ull.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};
	
	//d_l.i = d_i = d_ij^j
	real3 d_l = real3x3x3_tr23(d_llu);

	//partial_d_lll.ij.kl = d_kij,l = d_(k|(ij),|l)
	//so this object's indexes are rearranged compared to the papers 
	sym3sym3 partial_d_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
		local k,l,xk,xl = from6to3x3(kl)
?>	partial_d_llll.<?=xij?>.<?=xkl?> = 0.
<?		if l <= solver.dim then
?>	+ .5 * (	// 1/2 d_kij,l
		U[solver->stepsize.<?=xl?>].d_lll.<?=xk?>.<?=xij?>
		- U[-solver->stepsize.<?=xl?>].d_lll.<?=xk?>.<?=xij?>
	) / (2. * solver->grid_dx.<?=xl?>);
<?
		end
		if k <= solver.dim then
?>	+ .5 * (
		U[solver->stepsize.<?=xk?>].d_lll.<?=xl?>.<?=xij?>
		- U[-solver->stepsize.<?=xk?>].d_lll.<?=xl?>.<?=xij?>
	) / (2. * solver->grid_dx.<?=xk?>);
<?		end
?>	;
<?	end
end
?>

	// R_ll.ij := R_ij
	//	= gamma^kl (-gamma_ij,kl - gamma_kl,ij + gamma_ik,jl + gamma_jl,ik)
	//		+ conn^k_ij (d_k - 2 e_k)
	//		- 2 d^l_ki d^k_lj 
	//		+ 2 d^l_ki d_lj^k 
	//		+ d_il^k d_jk^l
	sym3 R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * (d_l.<?=xk?> - 2. * e_l.<?=xk?>)
<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_ull.<?=xl?>.<?=sym(k,i)?> * (d_llu.<?=xl?>.<?=xj?>.<?=xk?> - d_ull.<?=xk?>.<?=sym(l,j)?>)
			+ d_llu.<?=xi?>.<?=xl?>.<?=xk?> * d_llu.<?=xj?>.<?=xk?>.<?=xl?>
			+ gamma_uu.<?=sym(k,l)?> * (
				- partial_d_llll.<?=xij?>.<?=sym(k,l)?>
				- partial_d_llll.<?=sym(k,l)?>.<?=xij?>
				+ partial_d_llll.<?=sym(i,k)?>.<?=sym(j,l)?>
				+ partial_d_llll.<?=sym(j,l)?>.<?=sym(i,k)?>
			)
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
