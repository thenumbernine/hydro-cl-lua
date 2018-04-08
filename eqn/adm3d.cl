kernel void calcDT(
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma) is only called once
	real f = calc_f(U->alpha);
	real det_gamma = sym3_det(U->gamma);
	real sqrt_f = sqrt(f);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side==0 then ?>
		real gammaUjj = (U->gamma.yy * U->gamma.zz - U->gamma.yz * U->gamma.yz) / det_gamma;
		<? elseif side==1 then ?>
		real gammaUjj = (U->gamma.xx * U->gamma.zz - U->gamma.xz * U->gamma.xz) / det_gamma;
		<? elseif side==2 then ?>
		real gammaUjj = (U->gamma.xx * U->gamma.yy - U->gamma.xy * U->gamma.xy) / det_gamma;
		<? end ?>	
		real lambdaLight = U->alpha * sqrt(gammaUjj);
		
		real lambdaGauge = lambdaLight * sqrt_f;
		real lambda = (real)max(lambdaGauge, lambdaLight);

		<? if eqn.useShift then ?>
		real betaUi = U->beta_u.s<?=side?>;
		<? else ?>
		const real betaUi = 0.;
		<? end ?>
		
		real lambdaMin = (real)min((real)0., -betaUi - lambda);
		real lambdaMax = (real)max((real)0., -betaUi + lambda);

		dt = (real)min((real)dt, (real)(grid_dx<?=side?> / (fabs(lambdaMax - lambdaMin) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt; 
}

//used by PLM
<? 
for side=0,solver.dim-1 do 
	for _,addrAndSuffix in ipairs{
		{[''] = 'local'}, 	-- local addr has 'local' suffix
		{global = ''},		-- global addr has no suffix
	} do
	local addr, suffix = next(addrAndSuffix)
?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?><?=suffix?>(
	<?=addr?> const <?=eqn.cons_t?>* U,
	real3 x 
) {
	<?=eqn.eigen_t?> eig;
	eig.alpha = U->alpha;
	eig.sqrt_f = sqrt(calc_f(U->alpha));
	real det_gamma = sym3_det(U->gamma);
	eig.gammaU = sym3_inv(U->gamma, det_gamma);
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gammaU.xx), sqrt(eig.gammaU.yy), sqrt(eig.gammaU.zz));
	
	<? if eqn.useShift then ?>
	eig.beta_u = U->beta_u;
	<? end ?>

	return eig;
}
<? 
	end
end 
?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	real det_gamma = sym3_det(U->gamma);
	
	<? if side==0 then ?>
	real gammaUjj = (U->gamma.yy * U->gamma.zz - U->gamma.yz * U->gamma.yz) / det_gamma;
	<? elseif side==1 then ?>
	real gammaUjj = (U->gamma.xx * U->gamma.zz - U->gamma.xz * U->gamma.xz) / det_gamma;
	<? elseif side==2 then ?>
	real gammaUjj = (U->gamma.xx * U->gamma.yy - U->gamma.xy * U->gamma.xy) / det_gamma;
	<? end ?>
	
	real lambdaLight = U->alpha * sqrt(gammaUjj);
	
	real f = calc_f(U->alpha);
	real lambdaGauge = lambdaLight * sqrt(f);

	real lambdaMax = max(lambdaGauge, lambdaLight);
	//= lambdaLight * max(sqrt(f), 1)
	real lambdaMin = -lambdaMin;

	<? if eqn.useShift then ?>
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
<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	real alpha = .5 * (UL->alpha + UR->alpha);
	sym3 avg_gamma = (sym3){
		.xx = .5 * (UL->gamma.xx + UR->gamma.xx),
		.xy = .5 * (UL->gamma.xy + UR->gamma.xy),
		.xz = .5 * (UL->gamma.xz + UR->gamma.xz),
		.yy = .5 * (UL->gamma.yy + UR->gamma.yy),
		.yz = .5 * (UL->gamma.yz + UR->gamma.yz),
		.zz = .5 * (UL->gamma.zz + UR->gamma.zz),
	};
	real det_avg_gamma = sym3_det(avg_gamma);

	<?=eqn.eigen_t?> eig;
	eig.alpha = alpha;
	eig.sqrt_f = sqrt(calc_f(alpha));
	eig.gammaU = sym3_inv(avg_gamma, det_avg_gamma);
	eig.sqrt_gammaUjj.x = sqrt(eig.gammaU.xx);
	eig.sqrt_gammaUjj.y = sqrt(eig.gammaU.yy);
	eig.sqrt_gammaUjj.z = sqrt(eig.gammaU.zz);
	
	<? if eqn.useShift then ?>
	eig.beta_u = real3_scale(real3_add(UL->beta_u, UR->beta_u), .5);
	<? end ?>

	return eig;
}
<? end ?>

kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 xR = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize.s<?=side?>;
		
		<?= solver.getULRCode ?>	
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		int indexInt = side + dim * index;	
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		*eig = eigen_forSide_<?=side?>(UL, UR, xInt);
	}<? end ?>
}

<?
local unpack = unpack or table.unpack
for _,addrs in ipairs{
	{'', '', ''},	-- only used by fluxFromCons below
	{'', 'global', ''},
	{'global', 'global', ''},
} do
	local addr0, addr1, addr2 = unpack(addrs)
	for side=0,solver.dim-1 do ?>
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 unused
) {
	<?=addr2?> const <?=eqn.cons_t?>* inputU = (<?=addr2?> const <?=eqn.cons_t?>*)input;

<? if not eqn.noZeroRowsInFlux then ?>

	<? if side == 0 then ?>

	//a_y, a_z
	results[6] = inputU->a.y;
	results[7] = inputU->a.z;

	//d_yij
	results[8] = inputU->d.y.xx;
	results[9] = inputU->d.y.xy;
	results[10] = inputU->d.y.xz;
	results[11] = inputU->d.y.yy;
	results[12] = inputU->d.y.yz;
	results[13] = inputU->d.y.zz;

	//d_zij
	results[14] = inputU->d.z.xx;
	results[15] = inputU->d.z.xy;
	results[16] = inputU->d.z.xz;
	results[17] = inputU->d.z.yy;
	results[18] = inputU->d.z.yz;
	results[19] = inputU->d.z.zz;

	//V_j
	results[20] = inputU->V.x;
	results[21] = inputU->V.y;
	results[22] = inputU->V.z;
	
	sym3 K_sqrt_gammaUxx = sym3_scale(inputU->K, eig->sqrt_gammaUjj.x);

	//a^x - f d^xj_j

	real f = eig->sqrt_f * eig->sqrt_f;
	real d_x_input = sym3_dot(eig->gammaU, inputU->d.x);
	results[23] = inputU->a.x - f * d_x_input;

	//gauge:
	//sqrt(f gamma^xx) K +- (a^x + 2 V^x)

	real ev0a = eig->sqrt_f * sym3_dot(eig->gammaU, K_sqrt_gammaUxx);
	real ev0b = eig->gammaU.xx * (inputU->a.x + 2. * inputU->V.x) 
				+ eig->gammaU.xy * (inputU->a.y + 2. * inputU->V.y)
				+ eig->gammaU.xz * (inputU->a.z + 2. * inputU->V.z);
	results[0] = ev0a - ev0b;
	results[29] = ev0a + ev0b;

	//light:
	//sqrt(gamma^xx) K_xy +- (d^x_xy + .5 (a_y - d_yj^j) + V_y)

	real d_y_input = sym3_dot(eig->gammaU, inputU->d.y);
	real dUx_xy_input = eig->gammaU.xx * inputU->d.x.xy + eig->gammaU.xy * inputU->d.y.xy + eig->gammaU.xz * inputU->d.z.xy;
	real ev1b = .5 * (inputU->a.y - d_y_input) + inputU->V.y + dUx_xy_input;
	results[1] = K_sqrt_gammaUxx.xy - ev1b;
	results[24] = K_sqrt_gammaUxx.xy + ev1b;

	//light:
	//sqrt(gamma^xx) K_xz +- (d^x_xz + .5 (a_z - d_zj^j) + V_z)

	real d_z_input = sym3_dot(eig->gammaU, inputU->d.z);
	real dUx_xz_input = eig->gammaU.xx * inputU->d.x.xz + eig->gammaU.xy * inputU->d.y.xz + eig->gammaU.xz * inputU->d.z.xz;
	real ev2b = .5 * (inputU->a.z - d_z_input) + inputU->V.z + dUx_xz_input;
	results[2] = K_sqrt_gammaUxx.xz - ev2b;
	results[25] = K_sqrt_gammaUxx.xz + ev2b;

	//light:
	//sqrt(gamma^xx) K_yy +- d^x_yy

	real dUx_yy_input = eig->gammaU.xx * inputU->d.x.yy + eig->gammaU.xy * inputU->d.y.yy + eig->gammaU.xz * inputU->d.z.yy;
	results[3] = K_sqrt_gammaUxx.yy - dUx_yy_input;
	results[26] = K_sqrt_gammaUxx.yy + dUx_yy_input;

	//light:
	//sqrt(gamma^xx) K_yz +- d^x_yz

	real dUx_yz_input = eig->gammaU.xx * inputU->d.x.yz + eig->gammaU.xy * inputU->d.y.yz + eig->gammaU.xz * inputU->d.z.yz;
	results[4] = K_sqrt_gammaUxx.yz - dUx_yz_input; 
	results[27] = K_sqrt_gammaUxx.yz + dUx_yz_input;

	//light:
	//sqrt(gamma^xx) K_zz +- d^x_zz

	real dUx_zz_input = eig->gammaU.xx * inputU->d.x.zz + eig->gammaU.xy * inputU->d.y.zz + eig->gammaU.xz * inputU->d.z.zz;
	results[5] = K_sqrt_gammaUxx.zz - dUx_zz_input;
	results[28] = K_sqrt_gammaUxx.zz + dUx_zz_input;

	<? elseif side == 1 then ?>

	//a_x, a_z
	results[6] = inputU->a.x;
	results[7] = inputU->a.z;

	//d_xij
	results[8] = inputU->d.x.xx;
	results[9] = inputU->d.x.xy;
	results[10] = inputU->d.x.xz;
	results[11] = inputU->d.x.yy;
	results[12] = inputU->d.x.yz;
	results[13] = inputU->d.x.zz;
	
	//d_zij
	results[14] = inputU->d.z.xx;
	results[15] = inputU->d.z.xy;
	results[16] = inputU->d.z.xz;
	results[17] = inputU->d.z.yy;
	results[18] = inputU->d.z.yz;
	results[19] = inputU->d.z.zz;
	
	//V_j
	results[20] = inputU->V.x;
	results[21] = inputU->V.y;
	results[22] = inputU->V.z;
	
	sym3 K_sqrt_gammaUyy = sym3_scale(inputU->K, eig->sqrt_gammaUjj.y);

	//a^y - f d^yj_j

	real f = eig->sqrt_f * eig->sqrt_f;
	real d_y_input = sym3_dot(eig->gammaU, inputU->d.y);
	results[23] = inputU->a.y - f * d_y_input;
	
	//gauge:
	//sqrt(f gamma^yy) K +- (a^y + 2 V^y)

	real ev0a = eig->sqrt_f * sym3_dot(eig->gammaU, K_sqrt_gammaUyy);
	real ev0b = eig->gammaU.xy * (inputU->a.x + 2. * inputU->V.x)
				+ eig->gammaU.yy * (inputU->a.y + 2. * inputU->V.y)
				+ eig->gammaU.yz * (inputU->a.z + 2. * inputU->V.z);
	results[0] = ev0a - ev0b;
	results[29] = ev0a + ev0b;

	//light:
	//sqrt(gamma^yy) K_xx +- d^y_xx

	real dUy_xx_input = eig->gammaU.xy * inputU->d.x.xx + eig->gammaU.yy * inputU->d.y.xx + eig->gammaU.yz * inputU->d.z.xx;
	results[1] = K_sqrt_gammaUyy.xx - dUy_xx_input;
	results[24] = K_sqrt_gammaUyy.xx + dUy_xx_input;

	//light:
	//sqrt(gamma^yy) K_xy +- (d^y_xy + .5 (a_x - d_xj^j) + V_x)

	real d_x_input = sym3_dot(eig->gammaU, inputU->d.x);
	real dUy_xy_input = eig->gammaU.xy * inputU->d.x.xy + eig->gammaU.yy * inputU->d.y.xy + eig->gammaU.yz * inputU->d.z.xy;
	real ev2b = dUy_xy_input + .5 * (inputU->a.x - d_x_input) + inputU->V.x;
	results[2] = K_sqrt_gammaUyy.xy - ev2b;
	results[25] = K_sqrt_gammaUyy.xy + ev2b;

	//light:
	//sqrt(gamma^yy) K_xz +- d^y_xz

	real dUy_xz_input = eig->gammaU.xy * inputU->d.x.xz + eig->gammaU.yy * inputU->d.y.xz + eig->gammaU.yz * inputU->d.z.xz;
	results[3] = K_sqrt_gammaUyy.xz - dUy_xz_input;
	results[26] = K_sqrt_gammaUyy.xz + dUy_xz_input;

	//light:
	//sqrt(gamma^yy) K_yz +- (d^y_yz + .5 (a_z - d_zj^j) + V_z)

	real dUy_yz_input = eig->gammaU.xy * inputU->d.x.yz + eig->gammaU.yy * inputU->d.y.yz + eig->gammaU.yz * inputU->d.z.yz;
	real d_z_input = sym3_dot(eig->gammaU, inputU->d.z);
	real ev4b = dUy_yz_input + .5 * (inputU->a.z - d_z_input) + inputU->V.z;
	results[4] = K_sqrt_gammaUyy.yz - ev4b;
	results[27] = K_sqrt_gammaUyy.yz + ev4b;

	//light:
	//sqrt(gamma^yy) K_zz +- d^y_zz

	real dUy_zz_input = eig->gammaU.xy * inputU->d.x.zz + eig->gammaU.yy * inputU->d.y.zz + eig->gammaU.yz * inputU->d.z.zz;
	results[5] = K_sqrt_gammaUyy.zz - dUy_zz_input;
	results[28] = K_sqrt_gammaUyy.zz - dUy_zz_input;
	
	<? elseif side == 2 then ?>

	//a_x, a_y
	results[6] = inputU->a.x;
	results[7] = inputU->a.y;
	
	//d_xij
	results[8] =  inputU->d.x.xx;
	results[9] =  inputU->d.x.xy;
	results[10] = inputU->d.x.xz;
	results[11] = inputU->d.x.yy;
	results[12] = inputU->d.x.yz;
	results[13] = inputU->d.x.zz;
	
	//d_yij
	results[14] = inputU->d.y.xx;
	results[15] = inputU->d.y.xy;
	results[16] = inputU->d.y.xz;
	results[17] = inputU->d.y.yy;
	results[18] = inputU->d.y.yz;
	results[19] = inputU->d.y.zz;
	
	//V_j
	results[20] = inputU->V.x;
	results[21] = inputU->V.y;
	results[22] = inputU->V.z;

	sym3 K_sqrt_gammaUzz = sym3_scale(inputU->K, eig->sqrt_gammaUjj.z);

	//a^z - f d^zj_j

	real f = eig->sqrt_f * eig->sqrt_f;
	real d_z_input = sym3_dot(eig->gammaU, inputU->d.z);
	results[23] = inputU->a.z - f * d_z_input;

	//gauge:
	//sqrt(f gamma^zz) K +- (a^z + 2 V^z)

	real ev0a = eig->sqrt_f * sym3_dot(eig->gammaU, K_sqrt_gammaUzz);
	real ev0b = eig->gammaU.xz * (inputU->a.x + 2. * inputU->V.x)
				+ eig->gammaU.yz * (inputU->a.y + 2. * inputU->V.y)
				+ eig->gammaU.zz * (inputU->a.z + 2. * inputU->V.z);
	results[0] = ev0a - ev0b;
	results[29] = ev0a + ev0b;

	//light:
	//sqrt(gamma^zz) K_xx +- d^z_xx
	
	real dUz_xx_input = eig->gammaU.xz * inputU->d.x.xx + eig->gammaU.yz * inputU->d.y.xx + eig->gammaU.zz * inputU->d.z.xx;
	results[1] = K_sqrt_gammaUzz.xx - dUz_xx_input;
	results[24] = K_sqrt_gammaUzz.xx + dUz_xx_input;

	//light:
	//sqrt(gamma^zz) K_xy +- d^z_xy

	real dUz_xy_input = eig->gammaU.xz * inputU->d.x.xy + eig->gammaU.yz * inputU->d.y.xy + eig->gammaU.zz * inputU->d.z.xy;
	results[2] = K_sqrt_gammaUzz.xy - dUz_xy_input;
	results[25] = K_sqrt_gammaUzz.xy + dUz_xy_input;

	//light:
	//sqrt(gamma^zz) K_xz +- (d^z_xz + .5 (a_x - d_xj^j) + V_x)
	
	real d_x_input = sym3_dot(eig->gammaU, inputU->d.x);
	real dUz_xz_input = eig->gammaU.xz * inputU->d.x.xz + eig->gammaU.yz * inputU->d.y.xz + eig->gammaU.zz * inputU->d.z.xz;
	real ev3b = .5 * (inputU->a.x - d_x_input) + inputU->V.x + dUz_xz_input;
	results[3] = K_sqrt_gammaUzz.xz - ev3b;
	results[26] = K_sqrt_gammaUzz.xz + ev3b;

	//light:
	//sqrt(gamma^zz) K_yy +- d^z_yy

	real dUz_yy_input = eig->gammaU.xz * inputU->d.x.yy + eig->gammaU.yz * inputU->d.y.yy + eig->gammaU.zz * inputU->d.z.yy;
	results[4] = K_sqrt_gammaUzz.yy - dUz_yy_input;
	results[27] = K_sqrt_gammaUzz.yy + dUz_yy_input;
	
	//light:
	//sqrt(gamma^zz) K_yz

	real d_y_input = sym3_dot(eig->gammaU, inputU->d.y);
	real dUz_yz_input = eig->gammaU.xz * inputU->d.x.yz + eig->gammaU.yz * inputU->d.y.yz + eig->gammaU.zz * inputU->d.z.yz;
	real ev5b = .5 * (inputU->a.y - d_y_input) + inputU->V.y + dUz_yz_input;
	results[5] = K_sqrt_gammaUzz.yz - ev5b;
	results[28] = K_sqrt_gammaUzz.yz + ev5b;

	<? end -- side ?>

<? else -- eqn.noZeroRowsInFlux ?>

	real _1_sqrt_f = 1. / eig->sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f; 
	real sqrt_gammaUjj = eig->sqrt_gammaUjj.s<?=side?>;
	real _1_gammaUjj = 1. / eig->gammaU.s<?=side?><?=side?>;

	real a = inputU->a.s<?=side?>;

	//now swap x and side on the sym3's
	sym3 d = sym3_swap<?=side?>(inputU->d.v<?=side?>);
	sym3 K = sym3_swap<?=side?>(inputU->K);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);

	real K_dot_eig_gamma = sym3_dot(K, gammaU);
	real dj_dot_eig_gamma = sym3_dot(d, gammaU);

	results[0] = (a * -sqrt_gammaUjj * _1_sqrt_f + K_dot_eig_gamma) * .5 * _1_gammaUjj;

		<? for i=1,5 do ?>
	results[<?=i?>] = .5 * (-sqrt_gammaUjj * d.s[<?=i?>] + K.s[<?=i?>]);
		<? end ?>

	results[6] = (-a * _1_f + dj_dot_eig_gamma) * _1_gammaUjj;

		<? for i=1,5 do ?>
	results[<?=6+i?>] = .5 * (sqrt_gammaUjj * d.s[<?=i?>] + K.s[<?=i?>]);
		<? end ?>

	results[12] = (a * sqrt_gammaUjj * _1_sqrt_f + K_dot_eig_gamma) * .5 * _1_gammaUjj;
	
<? end -- eqn.noZeroRowsInFlux ?>
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 unused
) {
	<?=addr0?> <?=eqn.cons_t?>* resultU = (<?=addr0?> <?=eqn.cons_t?>*)results;

<? if not eqn.noZeroRowsInFlux then ?>
	
	<? if side == 0 then ?>

	//write zeros to the alpha and gammaLL terms
	resultU->alpha = 0;
	resultU->gamma = _sym3(0,0,0,0,0,0);

	resultU->a.y = input[6];
	resultU->a.z = input[7];
	
	resultU->d.y.xx = input[8];
	resultU->d.y.xy = input[9];
	resultU->d.y.xz = input[10];
	resultU->d.y.yy = input[11];
	resultU->d.y.yz = input[12];
	resultU->d.y.zz = input[13];
	
	resultU->d.z.xx = input[14];
	resultU->d.z.xy = input[15];
	resultU->d.z.xz = input[16];
	resultU->d.z.yy = input[17];
	resultU->d.z.yz = input[18];
	resultU->d.z.zz = input[19];
	
	resultU->V.x = input[20];
	resultU->V.y = input[21];
	resultU->V.z = input[22];

	real sqrt_gammaUxx = eig->sqrt_gammaUjj.x;
	real _1_sqrt_gammaUxx = 1. / sqrt_gammaUxx;
	real _1_gammaUxx = _1_sqrt_gammaUxx * _1_sqrt_gammaUxx;
	real invDenom = .5 * _1_gammaUxx;

	real VUx_input = eig->gammaU.xx * input[20]		//V_x
					+ eig->gammaU.xy * input[21]	//V_y
					+ eig->gammaU.xz * input[22];	//V_z

	real gauge_diff = input[0] - input[29];
	real gauge_sum = input[0] + input[29];

	resultU->a.x = -(
			gauge_diff
			+ 2. * (
				eig->gammaU.xy * input[6]	//a_y
				+ eig->gammaU.xz * input[7]	//a_z
			)	
			+ 4. * VUx_input
		) * invDenom;
	
	real K_input_minus = 
		2. * (
			eig->gammaU.xy * input[1]
			+ eig->gammaU.xz * input[2]
			+ eig->gammaU.yz * input[4]
		) + eig->gammaU.yy * input[3]
		+ eig->gammaU.zz * input[5];

	real K_input_plus = 
		2. * (
			eig->gammaU.xy * input[24]
			+ eig->gammaU.xz * input[25]
			+ eig->gammaU.yz * input[27]
		) + eig->gammaU.yy * input[26]
		+ eig->gammaU.zz * input[28];

	real _1_sqrt_f = 1. / eig->sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f;

	resultU->d.x.xx = -(
			- K_input_minus
			+ K_input_plus

			+ 2. * eig->gammaU.xy * (
				eig->gammaU.xx * input[8]	//d_yxx
				- input[6]		//a_y
				- 2. * input[21]	//V_y
			)
			+ 2. * eig->gammaU.xz * (
				eig->gammaU.xx * input[14]	//d_zxx
				- input[7]		//a_z
				- 2. * input[22]	//V_z
			)
			//(a^x + 2 V^x) / f - 2 gamma^xx d^xj_j
			+ (
				+ gauge_diff	//-ev0b = -a^x - 2 V^x
				+ 2. * (
					+ eig->gammaU.xx * input[23]	//a_x - f d^xj_j
					+ eig->gammaU.xy * input[6]	//a_y
					+ eig->gammaU.xz * input[7]	//a_z
					+ 2. * VUx_input
				)
			) * _1_f
		) * invDenom * _1_gammaUxx;

	//d_yj^j
	real d_y_input = 
		eig->gammaU.xx * input[8]
		+ eig->gammaU.yy * input[11]
		+ eig->gammaU.zz * input[13]
		+ 2. * (
			eig->gammaU.xy * input[9]
			+ eig->gammaU.xz * input[10]
			+ eig->gammaU.yz * input[12]
		);
	
	resultU->d.x.xy = -(
			input[1]
			+ input[6]
			- d_y_input
			+ 2. * eig->gammaU.xy * input[9]
			+ 2. * eig->gammaU.xz * input[15]
			+ 2. * input[21]
			- input[24]
		) * invDenom;
	
	//d_zj^j
	real d_z_input = 
		eig->gammaU.xx * input[14]
		+ eig->gammaU.yy * input[17]
		+ eig->gammaU.zz * input[19]
		+ 2. * (
			eig->gammaU.xy * input[15]
			+ eig->gammaU.xz * input[16]
			+ eig->gammaU.yz * input[18]
		);

	resultU->d.x.xz = -(
			input[2]
			+ input[7]
			- d_z_input
			+ 2. * eig->gammaU.xy * input[10]
			+ 2. * eig->gammaU.xz * input[16]
			+ 2. * input[22]
			- input[25]
		) * invDenom;
	resultU->d.x.yy = -(
			input[3]
			+ 2. * eig->gammaU.xy * input[11]
			+ 2. * eig->gammaU.xz * input[17]
			- input[26]
		) * invDenom;
	resultU->d.x.yz = -(
			input[4]
			+ 2. * eig->gammaU.xy * input[12]
			+ 2. * eig->gammaU.xz * input[18]
			- input[27]
		) * invDenom;
	resultU->d.x.zz = -(
			input[5]
			+ 2. * eig->gammaU.xy * input[13]
			+ 2. * eig->gammaU.xz * input[19]
			- input[28]
		) * invDenom;

	resultU->K.xx = (
			- K_input_minus
			- K_input_plus	
			+ gauge_sum * _1_sqrt_f
		) * invDenom * _1_sqrt_gammaUxx;

	real tmp = .5 * _1_sqrt_gammaUxx;
	resultU->K.xy = (input[1] + input[24]) * tmp;
	resultU->K.xz = (input[2] + input[25]) * tmp;
	resultU->K.yy = (input[3] + input[26]) * tmp;
	resultU->K.yz = (input[4] + input[27]) * tmp;
	resultU->K.zz = (input[5] + input[28]) * tmp;

	<? elseif side == 1 then ?>
	
	//write zeros to the alpha and gammaLL terms
	resultU->alpha = 0;
	resultU->gamma = _sym3(0,0,0,0,0,0);
	
	resultU->a.x = input[6];
	resultU->a.z = input[7];
	
	resultU->d.x.xx = input[8];
	resultU->d.x.xy = input[9];
	resultU->d.x.xz = input[10];
	resultU->d.x.yy = input[11];
	resultU->d.x.yz = input[12];
	resultU->d.x.zz = input[13];
	
	resultU->d.z.xx = input[14];
	resultU->d.z.xy = input[15];
	resultU->d.z.xz = input[16];
	resultU->d.z.yy = input[17];
	resultU->d.z.yz = input[18];
	resultU->d.z.zz = input[19];
	
	resultU->V.x = input[20];
	resultU->V.y = input[21];
	resultU->V.z = input[22];

	real sqrt_gammaUyy = eig->sqrt_gammaUjj.y;
	real _1_sqrt_gammaUyy = 1. / sqrt_gammaUyy;
	real inv_gammaUyy = _1_sqrt_gammaUyy * _1_sqrt_gammaUyy;
	real invDenom = .5 * inv_gammaUyy;

	real VUy_input = eig->gammaU.xy * input[20]
					+ eig->gammaU.yy * input[21]
					+ eig->gammaU.yz * input[22];
	
	real gauge_diff = input[0] - input[29];
	real gauge_sum = input[0] + input[29];

	resultU->a.y = -(
			gauge_diff
			+ 2. * (
				eig->gammaU.xy * input[6]
				+ eig->gammaU.yz * input[7]
			)
			+ 4. * VUy_input
		) * invDenom;

	resultU->d.y.xx = -(
			+ input[1]
			- input[24]
			+ 2. * eig->gammaU.xy * input[8]
			+ 2. * eig->gammaU.yz * input[14]
		) * invDenom;

	//d_xj^j
	real d_x_input = 
		eig->gammaU.xx * input[8]
		+ eig->gammaU.yy * input[11] 
		+ eig->gammaU.zz * input[13]
		+ 2. * (
			eig->gammaU.xy * input[9]
			+ eig->gammaU.xz * input[10]
			+ eig->gammaU.yz * input[12]
		);
	
	resultU->d.y.xy = -(
			+ input[2]
			+ input[6]
			- d_x_input			
			+ 2. * eig->gammaU.xy * input[9]
			+ 2. * eig->gammaU.yz * input[15]
			+ 2. * input[20]
			- input[25]
		) * invDenom;
	
	resultU->d.y.xz = -(
			+ input[3]
			+ 2. * eig->gammaU.xy * input[10]
			+ 2. * eig->gammaU.yz * input[16]
			- input[26]
		) * invDenom;
	
	real K_input_minus = 
		2. * (
			eig->gammaU.xy * input[2]
			+ eig->gammaU.xz * input[3]
			+ eig->gammaU.yz * input[4]
		) + eig->gammaU.xx * input[1]
		+ eig->gammaU.zz * input[5];

	real K_input_plus = 
		2. * (
			eig->gammaU.xy * input[25]
			+ eig->gammaU.xz * input[26]
			+ eig->gammaU.yz * input[27]
		) + eig->gammaU.xx * input[24]
		+ eig->gammaU.zz * input[28];

	real _1_sqrt_f = 1. / eig->sqrt_f;
	real _1_f = _1_sqrt_f  * _1_sqrt_f;

	resultU->d.y.yy = -(
			- K_input_minus	
			+ K_input_plus	
			
			+ 2. * eig->gammaU.xy * (
				input[11] * eig->gammaU.yy
				- input[6]
				- 2. * input[20]
			)	
			+ 2. * eig->gammaU.yz * (
				input[17] * eig->gammaU.yy
				- input[7]
				- 2. * input[22]
			)
			+ (
				+ gauge_diff
				+ 2. * (
					+ eig->gammaU.yy * input[23]
					+ eig->gammaU.xy * input[6]
					+ eig->gammaU.yz * input[7]
					+ 2. * VUy_input
				)
			) * _1_f
		) * invDenom * inv_gammaUyy;

	//gamma_zj^j
	real d_z_input = eig->gammaU.xx * input[14]
		+ eig->gammaU.yy * input[17]
		+ eig->gammaU.zz * input[19]
		+ 2. * (
			eig->gammaU.xy * input[15]
			+ eig->gammaU.xz * input[16]
			+ eig->gammaU.yz * input[18]);

	resultU->d.y.yz = -(
			+ input[4]
			+ input[7]
			- d_z_input	
			+ 2. * eig->gammaU.xy * input[12]
			+ 2. * eig->gammaU.yz * input[18]
			+ 2. * input[22]
			- input[27]
		) * invDenom;
	
	resultU->d.y.zz = -(
			+ input[5]
			- input[28]
			+ 2. * eig->gammaU.xy * input[13]
			+ 2. * eig->gammaU.yz * input[19]
		) * invDenom;

	resultU->K.yy = (
			- K_input_minus
			- K_input_plus
			+ gauge_sum * _1_sqrt_f
		) * invDenom * _1_sqrt_gammaUyy;

	real tmp = .5 * _1_sqrt_gammaUyy;
	resultU->K.xx = (input[1] + input[24]) * tmp;
	resultU->K.xy = (input[2] + input[25]) * tmp;		//once this gets enabled, compiling crashes
	resultU->K.xz = (input[3] + input[26]) * tmp;
	resultU->K.yz = (input[4] + input[27]) * tmp;
	resultU->K.zz = (input[5] + input[28]) * tmp;
	
	<? elseif side == 2 then ?>

	//write zeros to the alpha and gammaLL terms
	resultU->alpha = 0;
	resultU->gamma = _sym3(0,0,0,0,0,0);
	
	resultU->a.x = input[6];
	resultU->a.y = input[7];
	
	resultU->d.x.xx = input[8];
	resultU->d.x.xy = input[9];
	resultU->d.x.xz = input[10];
	resultU->d.x.yy = input[11];
	resultU->d.x.yz = input[12];
	resultU->d.x.zz = input[13];
	
	resultU->d.y.xx = input[14];
	resultU->d.y.xy = input[15];
	resultU->d.y.xz = input[16];
	resultU->d.y.yy = input[17];
	resultU->d.y.yz = input[18];
	resultU->d.y.zz = input[19];

	resultU->V.x = input[20];
	resultU->V.y = input[21];
	resultU->V.z = input[22];

	real sqrt_gammaUzz = eig->sqrt_gammaUjj.z;
	real _1_sqrt_gammaUzz = 1. / sqrt_gammaUzz;
	real inv_gammaUzz = _1_sqrt_gammaUzz * _1_sqrt_gammaUzz;
	real invDenom = .5 * inv_gammaUzz;

	real VUz_input = eig->gammaU.xz * input[20]
					+ eig->gammaU.yz * input[21]
					+ eig->gammaU.zz * input[22];

	real gauge_diff = input[0] - input[29];
	real gauge_sum = input[0] + input[29];
	
	resultU->a.z = -(
			gauge_diff
			+ 2. * (
				eig->gammaU.xz * input[6]
				+ eig->gammaU.yz * input[7]
			)
			+ 4. * VUz_input
		) * invDenom;
	
	resultU->d.z.xx = -(
			+ input[1]
			- input[24]
			+ 2. * eig->gammaU.xz * input[8]
			+ 2. * eig->gammaU.yz * input[14]
		) * invDenom;
	
	resultU->d.z.xy = -(
			+ input[2]
			- input[25]
			+ 2. * eig->gammaU.xz * input[9]
			+ 2. * eig->gammaU.yz * input[15]
		) * invDenom;

	//d_xj^j
	real d_x_input = 
		eig->gammaU.xx * input[8]
		+ eig->gammaU.yy * input[11]
		+ eig->gammaU.zz * input[13]
		+ 2. * (
			eig->gammaU.xy * input[9]
			+ eig->gammaU.xz * input[10]
			+ eig->gammaU.yz * input[12]);

	resultU->d.z.xz = -(
			+ input[3]
			+ input[6]
			- d_x_input
			+ 2. * eig->gammaU.xz * input[10]
			+ 2. * eig->gammaU.yz * input[16]
			+ 2. * input[20]
			- input[26]
		) * invDenom;
	
	resultU->d.z.yy = -(
			+ input[4]
			- input[27]
			+ 2. * eig->gammaU.xz * input[11]
			+ 2. * eig->gammaU.yz * input[17]
		) * invDenom;

	//d_yj^j
	real d_y_input = 
		eig->gammaU.xx * input[14]
		+ eig->gammaU.yy * input[17]
		+ eig->gammaU.zz * input[19]
		+ 2. * (
			eig->gammaU.xy * input[15]
			+ eig->gammaU.xz * input[16]
			+ eig->gammaU.yz * input[18]);

	resultU->d.z.yz = -(
			+ input[5]
			+ input[7]
			- d_y_input	
			+ 2. * eig->gammaU.xz * input[12]
			+ 2. * eig->gammaU.yz * input[18]
			+ 2. * input[21]
			- input[28]
		) * invDenom;

	real K_input_minus = 
		2. * (
			eig->gammaU.xy * input[2]
			+ eig->gammaU.xz * input[3]
			+ eig->gammaU.yz * input[5]
		) + eig->gammaU.xx * input[1]
		+ eig->gammaU.yy * input[4];

	real K_input_plus = 
		2. * (eig->gammaU.xy * input[25]
			+ eig->gammaU.xz * input[26]
			+ eig->gammaU.yz * input[28]
		) + eig->gammaU.xx * input[24]
		+ eig->gammaU.yy * input[27];

	real _1_sqrt_f = 1. / eig->sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f;

	resultU->d.z.zz = -(
			- K_input_minus
			+ K_input_plus
			
			+ 2. * eig->gammaU.xz * (
				input[13] * eig->gammaU.zz
				- input[6]
				- 2. * input[20]
			)
			+ 2. * eig->gammaU.yz * (
				eig->gammaU.zz * input[19]
				- input[7]
				- 2. * input[21]
			)	
			+ (	
				+ gauge_diff
				+ 2. * (
					eig->gammaU.zz * input[23]
					+ eig->gammaU.xz * input[6]
					+ eig->gammaU.yz * input[7]
					+ 2. * VUz_input
				)
			) * _1_f
		) * invDenom * inv_gammaUzz;

	resultU->K.zz = (
			+ gauge_sum * _1_sqrt_f
			- K_input_minus
			- K_input_plus
		) * invDenom * _1_sqrt_gammaUzz;
	
	real tmp = .5 * _1_sqrt_gammaUzz;
	resultU->K.xx = (input[1] + input[24]) * tmp;
	resultU->K.xy = (input[2] + input[25]) * tmp;
	resultU->K.xz = (input[3] + input[26]) * tmp;
	resultU->K.yy = (input[4] + input[27]) * tmp;
	resultU->K.yz = (input[5] + input[28]) * tmp;
	
	<? end ?>

<? else -- eqn.noZeroRowsInFlux ?>

	//TODO swap size inside eigen_t structure
	//instead of doing it here
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);

	real input1_dot_gammaU = input[1] * 2. * gammaU.xy
		+ input[2] * 2. * gammaU.xz
		+ input[3] * gammaU.yy
		+ input[4] * 2. * gammaU.yz
		+ input[5] * gammaU.zz;
	real input7_dot_gammaU = input[7] * 2. * gammaU.xy
		+ input[8] * 2. * gammaU.xz
		+ input[9] * gammaU.yy
		+ input[10] * 2. * gammaU.yz
		+ input[11] * gammaU.zz;

	real _1_sqrt_f = 1. / eig->sqrt_f;
	real sqrt_gammaUjj = eig->sqrt_gammaUjj.s<?=side?>;
	real _1_sqrt_gammaUjj = 1. / sqrt_gammaUjj;
	real _1_gammaUjj = _1_sqrt_gammaUjj * _1_sqrt_gammaUjj; 

	resultU->a.s<?=side?> = eig->sqrt_f * sqrt_gammaUjj * (input[12] - input[0]);

	sym3 d, K;
	d.xx = (
		(input[12] - input[0]) * _1_sqrt_f
		+ (input1_dot_gammaU - input7_dot_gammaU) * _1_gammaUjj
	) * _1_sqrt_gammaUjj + input[6];

	<? for i=1,5 do ?>
	d.s[<?=i?>] = (input[<?=i+6?>] - input[<?=i?>]) * _1_sqrt_gammaUjj;
	<? end ?>

	K.xx = input[0] + input[12] - (input1_dot_gammaU + input7_dot_gammaU) * _1_gammaUjj;

	<? for i=1,5 do ?>
	K.s[<?=i?>] = input[<?=i?>] + input[<?=i+6?>];
	<? end ?>

	//now swap x and side on the sym3's
	resultU->d.v<?=side?> = sym3_swap<?=side?>(d);
	resultU->K = sym3_swap<?=side?>(K);

<? end 	-- eqn.noZeroRowsInFlux ?>
}
<?	end
end

for _,addrs in ipairs{
	{'', '', ''},		-- used by fluxFromCons
	{'', 'global', ''},	-- used by the error stuff
} do
	local addr0, addr1, addr2 = unpack(addrs)
	for side=0,solver.dim-1 do
?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 x
) {
	for (int i = 0; i < numStates; ++i) {
		results[i] = 0;
	}

<? if not eqn.noZeroRowsInFlux then ?>

	// TODO make this a default implementation somewhere.
	// If no one has provided one then just fall back on left / wave / right transform.
	// TODO use that static function for the calc waves as well
	

	real tmp[numWaves];
	eigen_leftTransform_<?=side?>__<?=addr1?>_<?=addr2?>(tmp, eig, input, x);

	<?=eqn:eigenWaveCodePrefix(side, 'eig', 'x')?>

	//hmm, needs access to the wave buf ...
<? for j=0,eqn.numWaves-1 do 
?>	tmp[<?=j?>] *= <?=eqn:eigenWaveCode(side, 'eig', 'x', j)?>;
<? end 
?>
	eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_(results, eig, tmp, x);

<? else -- noZeroRowsInFlux ?>

	<?=addr2?> const <?=eqn.cons_t?>* inputU = (<?=addr2?> const <?=eqn.cons_t?>*)input;
	<?=addr0?> <?=eqn.cons_t?>* resultU = (<?=addr0?> <?=eqn.cons_t?>*)results;
	
	real f = eig->sqrt_f * eig->sqrt_f;
	
	//now swap x and side on the sym3's
	sym3 input_d = sym3_swap<?=side?>(inputU->d.v<?=side?>);
	sym3 input_K = sym3_swap<?=side?>(inputU->K);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);

	resultU->a.s<?=side?> = sym3_dot(input_K, gammaU) * eig->alpha * f;
	sym3 result_d = sym3_scale(input_K, eig->alpha);
	sym3 result_K = sym3_scale(input_d, eig->alpha * gammaU.xx);
	result_K.xx += (inputU->a.s<?=side?> - sym3_dot(input_d, gammaU)) * eig->alpha;

	//now swap x and side on the sym3's
	resultU->d.v<?=side?> = sym3_swap<?=side?>(result_d);
	resultU->K = sym3_swap<?=side?>(result_K);

<? end -- noZeroRowsInFlux ?>
}
<?
	end
end
?>

//TODO make this the default implementation of fluxFromCons
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>local(&U, x);

	<?=eqn.cons_t?> F;
	eigen_fluxTransform_<?=side?>___(F.ptr, &eig, U.ptr, x);

	return F;
}
<? end ?>



kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf)
{
	SETBOUNDS_NOGHOST();
	const global <?=eqn.cons_t?>* U = UBuf + index;
	global <?=eqn.cons_t?>* deriv = derivBuf + index;

	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	real f = calc_f(U->alpha);	//could be based on alpha...

#if 0
	real density = 0;//state[STATE_DENSITY];
	real pressure = 0;//state[STATE_PRESSURE];
	real3 vel3 = _real3(0,0,0);//(real4)(state[STATE_VELOCITY_X], state[STATE_VELOCITY_Y], state[STATE_VELOCITY_Z], 0.);
	real vel3Sq = real3_dot(vel3, vel3);	//TODO use gamma
	real LorentzFactor = 1. / sqrt(1. - vel3Sq);
	real4 vel4_ = (real4)(vel3.x * LorentzFactor, vel3.y * LorentzFactor, vel3.z * LorentzFactor, LorentzFactor);
	real3 beta_ = _real3(0,0,0);
#endif

	//S_terms_ll[ij] = 4 pi (gamma_ij (S - rho) - 2 S_ij)
	//...where rho = n^a n^b T_ab
	//...and S_ij = proj T_ij
	//..and S = gamma^ij S_ij
#if 1
	sym3 S_terms_ll = (sym3){.s={0,0,0,0,0,0}};
#else
	sym3 S_terms_ll = (sym3){
	8. * M_PI * ((density + pressure) * vel4_.x * vel4_.x + .5 * (density - pressure) * U->gamma.xx),
	8. * M_PI * ((density + pressure) * vel4_.x * vel4_.y + .5 * (density - pressure) * U->gamma.xy),
	8. * M_PI * ((density + pressure) * vel4_.x * vel4_.z + .5 * (density - pressure) * U->gamma.xz),
	8. * M_PI * ((density + pressure) * vel4_.y * vel4_.y + .5 * (density - pressure) * U->gamma.yy),
	8. * M_PI * ((density + pressure) * vel4_.y * vel4_.z + .5 * (density - pressure) * U->gamma.yz),
	8. * M_PI * ((density + pressure) * vel4_.z * vel4_.z + .5 * (density - pressure) * U->gamma.zz),
	};
#endif
	
	//momentum term: j_i = n^a proj T_ai
#if 1
	real3 j_l = _real3(0,0,0);
#else
	real3 j_l = (real3){
		8. * M_PI * ((density + pressure) * vel4_.w * vel4_.x + pressure * beta_.x),
		8. * M_PI * ((density + pressure) * vel4_.w * vel4_.y + pressure * beta_.y),
		8. * M_PI * ((density + pressure) * vel4_.w * vel4_.z + pressure * beta_.z),
	};
#endif


	// source terms
	
	mat3 K_ul = sym3_sym3_mul(gammaU, U->K);				//K^i_j
	real trK = mat3_trace(K_ul);				//K^k_k
	sym3 KSq_ll = sym3_mat3_to_sym3_mul(U->K, K_ul);	//KSq_ij = K_ik K^k_j

	//d_llu = d_ij^k = d_ijl * gamma^lk
	mat3 d_llu[3] = {
<? for i,xi in ipairs(xNames) do
?>		sym3_sym3_mul(U->d.<?=xi?>, gammaU),
<? end
?>	};

	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 d_ull = sym3_3sym3_mul(gammaU, U->d);

	//d3_l = d^j_ji
	real3 d3_l = (real3){
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
?>			.<?=xij?> = d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end 
?>	};

	real3 a_V_d3_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = U->a.<?=xi?> + U->V.<?=xi?> - d3_l.<?=xi?>,
<? end
?>	};
	real3 a_V_d3_u = sym3_real3_mul(gammaU, a_V_d3_l);

	//srcK_ij = (-a_i a_j 
	//		+ conn^k_ij (a_k + V_k - d^l_lk) 
	//		+ 2 d_ki^l d^k_jl
	//		- 2 d_ki^l d_lj^k
	//		+ 2 d_ki^l d_jl^k
	//		+ 2 d_il^k d_kj^l
	//		- 3 d_il^k d_jk^l
	//		+ K K_ij - 2 K_ik K^k_j
	//		+ Sterms_ij
	sym3 srcK_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>		.<?=xij?> = 
			- U->a.<?=xi?> * U->a.<?=xj?>

<? 	for k,xk in ipairs(xNames) do 
?>			+ conn_ull.<?=xk?>.<?=xij?> * a_V_d3_u.<?=xk?>
<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(j,l)?>
			- 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=l-1?>].<?=xj?>.<?=xk?>
			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_llu[<?=j-1?>].<?=xl?>.<?=xk?>
			+ 2. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=k-1?>].<?=xj?>.<?=xl?>
			- 3. * d_llu[<?=i-1?>].<?=xl?>.<?=xk?> * d_llu[<?=j-1?>].<?=xk?>.<?=xl?>
<? 		end
	end
?>
			+ trK * U->K.<?=xij?>
			- 2. * KSq_ll.<?=xij?>
			
			- S_terms_ll.<?=xij?>,
<? end
?>	};

	//d1_l = d_ij^j
	real3 d1_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = mat3_trace(d_llu[<?=i-1?>]),
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
?>		.<?=xk?> = -trK * U->a.<?=xk?>
			+ 8. * M_PI * j_l.<?=xk?>
<?	for j,xj in ipairs(xNames) do
?>			+ U->a.<?=xj?> * K_ul.<?=xj?>.<?=xk?>
			- 4. * trK * d_ull.<?=xj?>.<?=sym(j,k)?>
			+ 2. * K_ul.<?=xj?>.<?=xk?> * d3_l.<?=xj?>
			+ 2. * K_ul.<?=xj?>.<?=xk?> * d1_l.<?=xj?>
<?		for i,xi in ipairs(xNames) do
?>			- K_ul.<?=xi?>.<?=xj?> * d_llu[<?=k-1?>].<?=xi?>.<?=xj?>
			+ 2. * K_ul.<?=xi?>.<?=xj?> * d_llu[<?=i-1?>].<?=xk?>.<?=xj?>
<? 		end
	end ?>,
<? end
?>	};


	//alpha_,t = first derivs - alpha^2 f gamma^ij K_ij
	deriv->alpha += -U->alpha * U->alpha * f * trK;
	
	//gamma_ij,t = first derivs - 2 alpha K_ij
	sym3_add(deriv->gamma, sym3_scale(U->K, -2. * U->alpha));
	
	//K_ij,t = first derivs + alpha srcK_ij
	sym3_add(deriv->K, sym3_scale(srcK_ll, U->alpha));

	//V_k,t = first derivs + alpha srcV_k
	real3_add(deriv->V, real3_scale(srcV_l, U->alpha));

<? if eqn.useShift then ?>



	<? if eqn.useShift == '2005 Bona / 2008 Yano' then ?>
	
	
	
	<? elseif eqn.useShift == 'HarmonicShiftCondition-FiniteDifference' then ?>

	//TODO make each useShift option into an object
	//and give it a method for producing the partial_beta_ul
	// so the finite difference ones can use this
	// and the state ones can just assign a variable
	//then all this can be moved into the "if eqn.useShift then" block
	//but then again, any FVS-based method will already have the U_,k beta^k terms incorporated into its eigensystem ...
	//this means the useShift will also determine whether the flux has any shift terms
	// for any '-FiniteDifference' shift, the flux won't have any shift terms.


	//partial_beta_ul[j].i := beta^i_,j
<?=makePartial('beta_u', 'real3')?>	

	//= gamma_ik beta^k_,j
	mat3 partial_beta_u_ll = (mat3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_real3_mul(gammaU, partial_beta_ul[<?=i-1?>]),
<? end
?>	};

	//first add the Lie derivative beta terms to the time derivatives

	//alpha_,t = alpha_,i beta^i + ...
	// = alpha a_i beta^i + ..
	deriv->alpha += real3_dot(U->beta_u, U->a) * U->alpha;

	//gamma_ij,t = gamma_ij,k beta^k + gamma_kj beta^k_,i + gamma_ik beta^k_,j
	// = 2 d_kij beta^k + gamma_kj beta^k_,i + gamma_ik beta^k_,j
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	deriv->gamma.<?=xij?> += partial_beta_u_ll.<?=xi?>.<?=xj?>
		+ partial_beta_u_ll.<?=xj?>.<?=xi?>
<?	for k,xk in ipairs(xNames) do
?> 		+ 2. * U->d.<?=xk?>.<?=xij?> * U->beta_u.<?=xk?>
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
	//partial_d_l[k][l].ij = d_kij,l
<?=makePartial('d', '_3sym3')?>

	//d_kij,t = d_kij,l beta^l + d_lij beta^l_,k + d_klj beta^l_,i + d_kil beta^l_,j
<? for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>	deriv->d.<?=xk?>.<?=xij?> += 0.
<?		for l,xl in ipairs(xNames) do
?>			+ partial_d_l[<?=l-1?>].<?=xk?>.<?=xij?>
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
<?=makePartial('V', 'real3')?>

	//V_i,t = V_i,j beta^j + V_j beta^i_,j
<? for i,xi in ipairs(xNames) do
?>	deriv->V.<?=xi?> += 0. <?
	for j,xj in ipairs(xNames) do
?>		+ partial_V_l[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
		+ U->V.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?	end
?>	;
<? end
?>

#endif

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
	
	real3 a_u = sym3_real3_mul(gammaU, U->a);

	//conn^i = conn^i_jk gamma^jk
	real3 conn_u = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(conn_ull.<?=xi?>, gammaU),
<? end
?>	};

	deriv->beta_u = 
	real3_add(
		deriv->beta_u,
		real3_add(
			real3_add(
				dbeta_beta,
				real3_scale(
					U->beta_u,
					U->alpha * trK * (1. - f)
				)
			),
			real3_scale(
				real3_sub(
					conn_u,
					a_u
				),
				U->alpha * U->alpha
			)
		)
	);

	<?
	elseif eqn.useShift == 'LagrangianCoordinates' then
	?>
	//nothing for LagrangianCoordinates shift -- this is handled by the beta_u advection operation
	<? else
		error("I don't have any source terms implemented for this particular useShift")
	end ?>
<? end ?>

	// and now for the first-order constraints

<? if eqn.guiVars.a_convCoeff.value ~= 0 then ?>
	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	<? for i,xi in ipairs(xNames) do ?>{
		<? if i <= solver.dim then ?>
		real di_alpha = (U[stepsize.<?=xi?>].alpha - U[-stepsize.<?=xi?>].alpha) / (2. * grid_dx<?=i-1?>);
		<? else ?>
		real di_alpha = 0.;
		<? end ?>
		deriv->a.<?=xi?> += gui_a_convCoeff * (di_alpha / U->alpha - U->a.<?=xi?>);
	}<? end ?>	
<? end -- eqn.guiVars.a_convCoeff.value  ?>
	
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
<? if eqn.guiVars.d_convCoeff.value ~= 0 then ?>
	<? 
for i,xi in ipairs(xNames) do 
	for jk,xjk in ipairs(symNames) do ?>{
		<? if i <= solver.dim then ?>
		real di_gamma_jk = (U[stepsize.<?=xi?>].gamma.<?=xjk?> - U[-stepsize.<?=xi?>].gamma.<?=xjk?>) / (2. * grid_dx<?=i-1?>);
		<? else ?>
		real di_gamma_jk = 0;
		<? end ?>
		deriv->d.<?=xi?>.<?=xjk?> += gui_d_convCoeff * (.5 * di_gamma_jk - U->d.<?=xi?>.<?=xjk?>);
	}<? 
	end
end ?>
<? end -- eqn.guiVars.d_convCoeff.value  ?>

<? if eqn.guiVars.V_convCoeff.value ~= 0 then ?>
	//V_i = d_ik^k - d^k_ki <=> V_i += eta (d_ik^k - d^k_ki - V_i)
	deriv->V = real3_add(
		deriv->V,
		real3_scale(
			real3_sub(real3_sub(d1_l, d3_l), U->V),
			gui_V_convCoeff));
<? end -- eqn.guiVars.V_convCoeff.value  ?>

	//Kreiss-Oligar diffusion, for stability's sake?
}

kernel void constrainU(
	global <?=eqn.cons_t?>* UBuf
) {

<?
local constrainVGuiVar = eqn.guiVars['constrain V']
local constrainV = constrainVGuiVar.options[constrainVGuiVar.value]  
if constrainV ~= 'none' then 
?>	//gravitational_wave_sim uses this (for 1D), HydroGPU doesn't (for 2D/3D)
	SETBOUNDS(numGhost,numGhost);	
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);

	real3 delta;
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d.<?=xi?>, gammaU);
		real d2 = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + U->d.<?=xj?>.<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
		end
	end ?>;
		delta.<?=xi?> = U->V.<?=xi?> - (d1 - d2);
	}<? end ?>

<?
	if constrainV == 'replace V' then
?>
	//directly assign to V
	U->V = real3_sub(U->V, delta);
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

	U->V = real3_sub(U->V, real3_scale(delta, v_weight));
<? 
		for i,xi in ipairs(xNames) do 
			for jk,xjk in ipairs(symNames) do
				local j,k = from6to3x3(jk)
				local xk = xNames[k]
?>	U->d.<?=xi?>.<?=xjk?> += (
		delta.<?=xi?> * U->gamma.<?=xjk?> 
		- delta.<?=xk?> * U->gamma.<?=sym(i,j)?>
	) * d_weight;
<?	
			end
		end 
	end
?>

//...or linearly project out the [V_i, U->d_ijk] vector
//...or do a single gradient descent step
<?
end
?>
}
