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
	real det_gamma = sym3_det(U->gamma);
	real f = calc_f(U->alpha);
	real sqrt_f = sqrt(f);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side==0 then ?>
		real gammaUxx = (U->gamma.yy * U->gamma.zz - U->gamma.yz * U->gamma.yz) / det_gamma;
		real lambdaLight = U->alpha * sqrt(gammaUxx);
		<? elseif side==1 then ?>
		real gammaUyy = (U->gamma.xx * U->gamma.zz - U->gamma.xz * U->gamma.xz) / det_gamma;
		real lambdaLight = U->alpha * sqrt(gammaUyy);
		<? elseif side==2 then ?>
		real gammaUzz = (U->gamma.xx * U->gamma.yy - U->gamma.xy * U->gamma.xy) / det_gamma;
		real lambdaLight = U->alpha * sqrt(gammaUzz);
		<? end ?>	
		
		real lambdaGauge = lambdaLight * sqrt_f;
		real lambda = (real)max(lambdaGauge, lambdaLight);
		
		real lambdaMin = (real)min((real)0., -lambda);
		real lambdaMax = (real)max((real)0., lambda);
		dt = (real)min((real)dt, (real)(grid_dx<?=side?> / (fabs(lambdaMax - lambdaMin) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt; 
}

//used by PLM
<? 
for side=0,solver.dim-1 do 
	-- TODO add _'s to the default eigen_forCell name
	-- to make it like eigen_calcWaves_...
	for _,addrAndSuffix in ipairs{
		{[''] = 'local'}, 	-- local addr has 'local' suffix
		{global = ''},		-- global addr has no suffix
	} do
	local addr, suffix = next(addrAndSuffix)
?>
void eigen_forCell_<?=side?><?=suffix?>(
	<?=eqn.eigen_t?>* eig,
	<?=addr?> const <?=eqn.cons_t?>* U,
	real3 x 
) {
	eig->alpha = U->alpha;
	eig->sqrt_f = sqrt(calc_f(U->alpha));
	real det_gamma = sym3_det(U->gamma);
	eig->gammaU = sym3_inv(U->gamma, det_gamma);
}
<? 
	end
end 
?>

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr0?> real* wave,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	real3 x
) {
	<? if side==0 then ?>
	real lambdaLight = eig->alpha * eig->sqrt_gammaUjj.x;
	<? elseif side==1 then ?>
	real lambdaLight = eig->alpha * eig->sqrt_gammaUjj.y;
	<? elseif side==2 then ?>
	real lambdaLight = eig->alpha * eig->sqrt_gammaUjj.z;
	<? end ?>
	real lambdaGauge = lambdaLight * eig->sqrt_f;
	
	wave[0] = -lambdaGauge;
	<? for i=1,5 do ?> wave[<?=i?>] = -lambdaLight; <? end ?>
	<? for i=6,23 do ?> wave[<?=i?>] = 0.; <? end ?>
	<? for i=24,28 do ?> wave[<?=i?>] = lambdaLight; <? end ?>
	wave[29] = lambdaGauge;
}
<?		end
	end
end ?>

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

	return (range_t){
		.min = -lambdaMax, 
		.max = lambdaMax,
	};
}
<? end ?>

//used for interface eigen basis
<?=eqn.eigen_t?> eigen_forSide(
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
	return eig;
}

kernel void calcEigenBasis(
	global real* waveBuf,
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
		*eig = eigen_forSide(UL, UR, xInt);
		
		global real* wave = waveBuf + numWaves * indexInt;
		
		eigen_calcWaves_<?=side?>_global_global(wave, eig, xInt);
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
	
	<? if side == 0 then ?>

	//a_y, a_z
	results[6] = inputU->a.y;
	results[7] = inputU->a.z;

	//d_yij
	results[8] = inputU->d[1].xx;
	results[9] = inputU->d[1].xy;
	results[10] = inputU->d[1].xz;
	results[11] = inputU->d[1].yy;
	results[12] = inputU->d[1].yz;
	results[13] = inputU->d[1].zz;

	//d_zij
	results[14] = inputU->d[2].xx;
	results[15] = inputU->d[2].xy;
	results[16] = inputU->d[2].xz;
	results[17] = inputU->d[2].yy;
	results[18] = inputU->d[2].yz;
	results[19] = inputU->d[2].zz;

	//V_j
	results[20] = inputU->V.x;
	results[21] = inputU->V.y;
	results[22] = inputU->V.z;
	
	sym3 K_sqrt_gammaUxx = sym3_scale(inputU->K, eig->sqrt_gammaUjj.x);

	//a^x - f d^xj_j

	real f = eig->sqrt_f * eig->sqrt_f;
	real d_x_input = sym3_dot(eig->gammaU, inputU->d[0]);
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

	real d_y_input = sym3_dot(eig->gammaU, inputU->d[1]);
	real dUx_xy_input = eig->gammaU.xx * inputU->d[0].xy + eig->gammaU.xy * inputU->d[1].xy + eig->gammaU.xz * inputU->d[2].xy;
	real ev1b = .5 * (inputU->a.y - d_y_input) + inputU->V.y + dUx_xy_input;
	results[1] = K_sqrt_gammaUxx.xy - ev1b;
	results[24] = K_sqrt_gammaUxx.xy + ev1b;

	//light:
	//sqrt(gamma^xx) K_xz +- (d^x_xz + .5 (a_z - d_zj^j) + V_z)

	real d_z_input = sym3_dot(eig->gammaU, inputU->d[2]);
	real dUx_xz_input = eig->gammaU.xx * inputU->d[0].xz + eig->gammaU.xy * inputU->d[1].xz + eig->gammaU.xz * inputU->d[2].xz;
	real ev2b = .5 * (inputU->a.z - d_z_input) + inputU->V.z + dUx_xz_input;
	results[2] = K_sqrt_gammaUxx.xz - ev2b;
	results[25] = K_sqrt_gammaUxx.xz + ev2b;

	//light:
	//sqrt(gamma^xx) K_yy +- d^x_yy

	real dUx_yy_input = eig->gammaU.xx * inputU->d[0].yy + eig->gammaU.xy * inputU->d[1].yy + eig->gammaU.xz * inputU->d[2].yy;
	results[3] = K_sqrt_gammaUxx.yy - dUx_yy_input;
	results[26] = K_sqrt_gammaUxx.yy + dUx_yy_input;

	//light:
	//sqrt(gamma^xx) K_yz +- d^x_yz

	real dUx_yz_input = eig->gammaU.xx * inputU->d[0].yz + eig->gammaU.xy * inputU->d[1].yz + eig->gammaU.xz * inputU->d[2].yz;
	results[4] = K_sqrt_gammaUxx.yz - dUx_yz_input; 
	results[27] = K_sqrt_gammaUxx.yz + dUx_yz_input;

	//light:
	//sqrt(gamma^xx) K_zz +- d^x_zz

	real dUx_zz_input = eig->gammaU.xx * inputU->d[0].zz + eig->gammaU.xy * inputU->d[1].zz + eig->gammaU.xz * inputU->d[2].zz;
	results[5] = K_sqrt_gammaUxx.zz - dUx_zz_input;
	results[28] = K_sqrt_gammaUxx.zz + dUx_zz_input;

	<? elseif side == 1 then ?>

	//a_x, a_z
	results[6] = inputU->a.x;
	results[7] = inputU->a.z;

	//d_xij
	results[8] = inputU->d[0].xx;
	results[9] = inputU->d[0].xy;
	results[10] = inputU->d[0].xz;
	results[11] = inputU->d[0].yy;
	results[12] = inputU->d[0].yz;
	results[13] = inputU->d[0].zz;
	
	//d_zij
	results[14] = inputU->d[2].xx;
	results[15] = inputU->d[2].xy;
	results[16] = inputU->d[2].xz;
	results[17] = inputU->d[2].yy;
	results[18] = inputU->d[2].yz;
	results[19] = inputU->d[2].zz;
	
	//V_j
	results[20] = inputU->V.x;
	results[21] = inputU->V.y;
	results[22] = inputU->V.z;
	
	sym3 K_sqrt_gammaUyy = sym3_scale(inputU->K, eig->sqrt_gammaUjj.y);

	//a^y - f d^yj_j

	real f = eig->sqrt_f * eig->sqrt_f;
	real d_y_input = sym3_dot(eig->gammaU, inputU->d[1]);
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

	real dUy_xx_input = eig->gammaU.xy * inputU->d[0].xx + eig->gammaU.yy * inputU->d[1].xx + eig->gammaU.yz * inputU->d[2].xx;
	results[1] = K_sqrt_gammaUyy.xx - dUy_xx_input;
	results[24] = K_sqrt_gammaUyy.xx + dUy_xx_input;

	//light:
	//sqrt(gamma^yy) K_xy +- (d^y_xy + .5 (a_x - d_xj^j) + V_x)

	real d_x_input = sym3_dot(eig->gammaU, inputU->d[0]);
	real dUy_xy_input = eig->gammaU.xy * inputU->d[0].xy + eig->gammaU.yy * inputU->d[1].xy + eig->gammaU.yz * inputU->d[2].xy;
	real ev2b = dUy_xy_input + .5 * (inputU->a.x - d_x_input) + inputU->V.x;
	results[2] = K_sqrt_gammaUyy.xy - ev2b;
	results[25] = K_sqrt_gammaUyy.xy + ev2b;

	//light:
	//sqrt(gamma^yy) K_xz +- d^y_xz

	real dUy_xz_input = eig->gammaU.xy * inputU->d[0].xz + eig->gammaU.yy * inputU->d[1].xz + eig->gammaU.yz * inputU->d[2].xz;
	results[3] = K_sqrt_gammaUyy.xz - dUy_xz_input;
	results[26] = K_sqrt_gammaUyy.xz + dUy_xz_input;

	//light:
	//sqrt(gamma^yy) K_yz +- (d^y_yz + .5 (a_z - d_zj^j) + V_z)

	real dUy_yz_input = eig->gammaU.xy * inputU->d[0].yz + eig->gammaU.yy * inputU->d[1].yz + eig->gammaU.yz * inputU->d[2].yz;
	real d_z_input = sym3_dot(eig->gammaU, inputU->d[2]);
	real ev4b = dUy_yz_input + .5 * (inputU->a.z - d_z_input) + inputU->V.z;
	results[4] = K_sqrt_gammaUyy.yz - ev4b;
	results[27] = K_sqrt_gammaUyy.yz + ev4b;

	//light:
	//sqrt(gamma^yy) K_zz +- d^y_zz

	real dUy_zz_input = eig->gammaU.xy * inputU->d[0].zz + eig->gammaU.yy * inputU->d[1].zz + eig->gammaU.yz * inputU->d[2].zz;
	results[5] = K_sqrt_gammaUyy.zz - dUy_zz_input;
	results[28] = K_sqrt_gammaUyy.zz - dUy_zz_input;
	
	<? elseif side == 2 then ?>

	//a_x, a_y
	results[6] = inputU->a.x;
	results[7] = inputU->a.y;
	
	//d_xij
	results[8] =  inputU->d[0].xx;
	results[9] =  inputU->d[0].xy;
	results[10] = inputU->d[0].xz;
	results[11] = inputU->d[0].yy;
	results[12] = inputU->d[0].yz;
	results[13] = inputU->d[0].zz;
	
	//d_yij
	results[14] = inputU->d[1].xx;
	results[15] = inputU->d[1].xy;
	results[16] = inputU->d[1].xz;
	results[17] = inputU->d[1].yy;
	results[18] = inputU->d[1].yz;
	results[19] = inputU->d[1].zz;
	
	//V_j
	results[20] = inputU->V.x;
	results[21] = inputU->V.y;
	results[22] = inputU->V.z;

	sym3 K_sqrt_gammaUzz = sym3_scale(inputU->K, eig->sqrt_gammaUjj.z);

	//a^z - f d^zj_j

	real f = eig->sqrt_f * eig->sqrt_f;
	real d_z_input = sym3_dot(eig->gammaU, inputU->d[2]);
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
	
	real dUz_xx_input = eig->gammaU.xz * inputU->d[0].xx + eig->gammaU.yz * inputU->d[1].xx + eig->gammaU.zz * inputU->d[2].xx;
	results[1] = K_sqrt_gammaUzz.xx - dUz_xx_input;
	results[24] = K_sqrt_gammaUzz.xx + dUz_xx_input;

	//light:
	//sqrt(gamma^zz) K_xy +- d^z_xy

	real dUz_xy_input = eig->gammaU.xz * inputU->d[0].xy + eig->gammaU.yz * inputU->d[1].xy + eig->gammaU.zz * inputU->d[2].xy;
	results[2] = K_sqrt_gammaUzz.xy - dUz_xy_input;
	results[25] = K_sqrt_gammaUzz.xy + dUz_xy_input;

	//light:
	//sqrt(gamma^zz) K_xz +- (d^z_xz + .5 (a_x - d_xj^j) + V_x)
	
	real d_x_input = sym3_dot(eig->gammaU, inputU->d[0]);
	real dUz_xz_input = eig->gammaU.xz * inputU->d[0].xz + eig->gammaU.yz * inputU->d[1].xz + eig->gammaU.zz * inputU->d[2].xz;
	real ev3b = .5 * (inputU->a.x - d_x_input) + inputU->V.x + dUz_xz_input;
	results[3] = K_sqrt_gammaUzz.xz - ev3b;
	results[26] = K_sqrt_gammaUzz.xz + ev3b;

	//light:
	//sqrt(gamma^zz) K_yy +- d^z_yy

	real dUz_yy_input = eig->gammaU.xz * inputU->d[0].yy + eig->gammaU.yz * inputU->d[1].yy + eig->gammaU.zz * inputU->d[2].yy;
	results[4] = K_sqrt_gammaUzz.yy - dUz_yy_input;
	results[27] = K_sqrt_gammaUzz.yy + dUz_yy_input;
	
	//light:
	//sqrt(gamma^zz) K_yz

	real d_y_input = sym3_dot(eig->gammaU, inputU->d[1]);
	real dUz_yz_input = eig->gammaU.xz * inputU->d[0].yz + eig->gammaU.yz * inputU->d[1].yz + eig->gammaU.zz * inputU->d[2].yz;
	real ev5b = .5 * (inputU->a.y - d_y_input) + inputU->V.y + dUz_yz_input;
	results[5] = K_sqrt_gammaUzz.yz - ev5b;
	results[28] = K_sqrt_gammaUzz.yz + ev5b;

	<? end ?>
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 unused
) {
	<?=addr0?> <?=eqn.cons_t?>* resultU = (<?=addr0?> <?=eqn.cons_t?>*)results;
	
	<? if side == 0 then ?>

	//write zeros to the alpha and gammaLL terms
	resultU->alpha = 0;
	resultU->gamma = _sym3(0,0,0,0,0,0);

	resultU->a.y = input[6];
	resultU->a.z = input[7];
	
	resultU->d[1].xx = input[8];
	resultU->d[1].xy = input[9];
	resultU->d[1].xz = input[10];
	resultU->d[1].yy = input[11];
	resultU->d[1].yz = input[12];
	resultU->d[1].zz = input[13];
	
	resultU->d[2].xx = input[14];
	resultU->d[2].xy = input[15];
	resultU->d[2].xz = input[16];
	resultU->d[2].yy = input[17];
	resultU->d[2].yz = input[18];
	resultU->d[2].zz = input[19];
	
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

	resultU->d[0].xx = -(
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
	
	resultU->d[0].xy = -(
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

	resultU->d[0].xz = -(
			input[2]
			+ input[7]
			- d_z_input
			+ 2. * eig->gammaU.xy * input[10]
			+ 2. * eig->gammaU.xz * input[16]
			+ 2. * input[22]
			- input[25]
		) * invDenom;
	resultU->d[0].yy = -(
			input[3]
			+ 2. * eig->gammaU.xy * input[11]
			+ 2. * eig->gammaU.xz * input[17]
			- input[26]
		) * invDenom;
	resultU->d[0].yz = -(
			input[4]
			+ 2. * eig->gammaU.xy * input[12]
			+ 2. * eig->gammaU.xz * input[18]
			- input[27]
		) * invDenom;
	resultU->d[0].zz = -(
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
	
	resultU->d[0].xx = input[8];
	resultU->d[0].xy = input[9];
	resultU->d[0].xz = input[10];
	resultU->d[0].yy = input[11];
	resultU->d[0].yz = input[12];
	resultU->d[0].zz = input[13];
	
	resultU->d[2].xx = input[14];
	resultU->d[2].xy = input[15];
	resultU->d[2].xz = input[16];
	resultU->d[2].yy = input[17];
	resultU->d[2].yz = input[18];
	resultU->d[2].zz = input[19];
	
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

	resultU->d[1].xx = -(
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
	
	resultU->d[1].xy = -(
			+ input[2]
			+ input[6]
			- d_x_input			
			+ 2. * eig->gammaU.xy * input[9]
			+ 2. * eig->gammaU.yz * input[15]
			+ 2. * input[20]
			- input[25]
		) * invDenom;
	
	resultU->d[1].xz = -(
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

	resultU->d[1].yy = -(
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

	resultU->d[1].yz = -(
			+ input[4]
			+ input[7]
			- d_z_input	
			+ 2. * eig->gammaU.xy * input[12]
			+ 2. * eig->gammaU.yz * input[18]
			+ 2. * input[22]
			- input[27]
		) * invDenom;
	
	resultU->d[1].zz = -(
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
	
	resultU->d[0].xx = input[8];
	resultU->d[0].xy = input[9];
	resultU->d[0].xz = input[10];
	resultU->d[0].yy = input[11];
	resultU->d[0].yz = input[12];
	resultU->d[0].zz = input[13];
	
	resultU->d[1].xx = input[14];
	resultU->d[1].xy = input[15];
	resultU->d[1].xz = input[16];
	resultU->d[1].yy = input[17];
	resultU->d[1].yz = input[18];
	resultU->d[1].zz = input[19];

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
	
	resultU->d[2].xx = -(
			+ input[1]
			- input[24]
			+ 2. * eig->gammaU.xz * input[8]
			+ 2. * eig->gammaU.yz * input[14]
		) * invDenom;
	
	resultU->d[2].xy = -(
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

	resultU->d[2].xz = -(
			+ input[3]
			+ input[6]
			- d_x_input
			+ 2. * eig->gammaU.xz * input[10]
			+ 2. * eig->gammaU.yz * input[16]
			+ 2. * input[20]
			- input[26]
		) * invDenom;
	
	resultU->d[2].yy = -(
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

	resultU->d[2].yz = -(
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

	resultU->d[2].zz = -(
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
}
<?	end
end
if solver.checkFluxError then 
	for _,addrs in ipairs{
		{'', 'global', ''},
	} do
		local addr0, addr1, addr2 = unpack(addrs)
		for side=0,solver.dim-1 do
?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x,
	real3 unused
) {
	for (int i = 0; i < numIntStates; ++i) {
		*y = 0;
		++y;
	}
}
<?
		end
	end
end
?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	//TODO do the math for this
	//until then, use the eigenbasis ...
	//TODO use this as a default for fluxFromCons as well

	<?=eqn.eigen_t?> eig;
	eigen_forCell_<?=side?>local(&eig, &U, x);

	real wave[numWaves];
	eigen_calcWaves_<?=side?>__(wave, &eig, x);

	real charvars[numWaves];
	eigen_leftTransform_<?=side?>___(charvars, &eig, &U, x);
	
	for (int j = 0; j < numWaves; ++j) {
		charvars[j] *= wave[j];
	}
			
	<?=eqn.cons_t?> F;
	eigen_rightTransform_<?=side?>___(F.ptr, &eig, charvars, x);

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

	// source terms
	
	mat3 KUL = sym3_sym3_mul(gammaU, U->K);				//K^i_j
	real trK = KUL.x.x + KUL.y.y + KUL.z.z;				//K^k_k
	sym3 KSqSymLL = sym3_mat3_to_sym3_mul(U->K, KUL);	//KSq_ij = K_ik K^k_j
	
	real DLUL[3][3][3] = {
	{{U->d[0].xx * gammaU.xx + U->d[0].xy * gammaU.xy + U->d[0].xz * gammaU.xz,
	U->d[0].xy * gammaU.xx + U->d[0].yy * gammaU.xy + U->d[0].yz * gammaU.xz,
	U->d[0].xz * gammaU.xx + U->d[0].yz * gammaU.xy + U->d[0].zz * gammaU.xz,
	},{U->d[0].xx * gammaU.xy + U->d[0].xy * gammaU.yy + U->d[0].xz * gammaU.yz,
	U->d[0].xy * gammaU.xy + U->d[0].yy * gammaU.yy + U->d[0].yz * gammaU.yz,
	U->d[0].xz * gammaU.xy + U->d[0].yz * gammaU.yy + U->d[0].zz * gammaU.yz,
	},{U->d[0].xx * gammaU.xz + U->d[0].xy * gammaU.yz + U->d[0].xz * gammaU.zz,
	U->d[0].xy * gammaU.xz + U->d[0].yy * gammaU.yz + U->d[0].yz * gammaU.zz,
	U->d[0].xz * gammaU.xz + U->d[0].yz * gammaU.yz + U->d[0].zz * gammaU.zz,
	},},{{U->d[1].xx * gammaU.xx + U->d[1].xy * gammaU.xy + U->d[1].xz * gammaU.xz,
	U->d[1].xy * gammaU.xx + U->d[1].yy * gammaU.xy + U->d[1].yz * gammaU.xz,
	U->d[1].xz * gammaU.xx + U->d[1].yz * gammaU.xy + U->d[1].zz * gammaU.xz,
	},{U->d[1].xx * gammaU.xy + U->d[1].xy * gammaU.yy + U->d[1].xz * gammaU.yz,
	U->d[1].xy * gammaU.xy + U->d[1].yy * gammaU.yy + U->d[1].yz * gammaU.yz,
	U->d[1].xz * gammaU.xy + U->d[1].yz * gammaU.yy + U->d[1].zz * gammaU.yz,
	},{U->d[1].xx * gammaU.xz + U->d[1].xy * gammaU.yz + U->d[1].xz * gammaU.zz,
	U->d[1].xy * gammaU.xz + U->d[1].yy * gammaU.yz + U->d[1].yz * gammaU.zz,
	U->d[1].xz * gammaU.xz + U->d[1].yz * gammaU.yz + U->d[1].zz * gammaU.zz,
	},},{{U->d[2].xx * gammaU.xx + U->d[2].xy * gammaU.xy + U->d[2].xz * gammaU.xz,
	U->d[2].xy * gammaU.xx + U->d[2].yy * gammaU.xy + U->d[2].yz * gammaU.xz,
	U->d[2].xz * gammaU.xx + U->d[2].yz * gammaU.xy + U->d[2].zz * gammaU.xz,
	},{U->d[2].xx * gammaU.xy + U->d[2].xy * gammaU.yy + U->d[2].xz * gammaU.yz,
	U->d[2].xy * gammaU.xy + U->d[2].yy * gammaU.yy + U->d[2].yz * gammaU.yz,
	U->d[2].xz * gammaU.xy + U->d[2].yz * gammaU.yy + U->d[2].zz * gammaU.yz,
	},{U->d[2].xx * gammaU.xz + U->d[2].xy * gammaU.yz + U->d[2].xz * gammaU.zz,
	U->d[2].xy * gammaU.xz + U->d[2].yy * gammaU.yz + U->d[2].yz * gammaU.zz,
	U->d[2].xz * gammaU.xz + U->d[2].yz * gammaU.yz + U->d[2].zz * gammaU.zz,
	},},};
	real D1L[3] = {
	DLUL[0][0][0] + DLUL[0][1][1] + DLUL[0][2][2],
	DLUL[1][0][0] + DLUL[1][1][1] + DLUL[1][2][2],
	DLUL[2][0][0] + DLUL[2][1][1] + DLUL[2][2][2],
	};
	real D3L[3] = {
	DLUL[0][0][0] + DLUL[1][1][0] + DLUL[2][2][0],
	DLUL[0][0][1] + DLUL[1][1][1] + DLUL[2][2][1],
	DLUL[0][0][2] + DLUL[1][1][2] + DLUL[2][2][2],
	};
	real DUUL[3][3][3] = {
	{{DLUL[0][0][0] * gammaU.xx + DLUL[1][0][0] * gammaU.xy + DLUL[2][0][0] * gammaU.xz,
	DLUL[0][0][1] * gammaU.xx + DLUL[1][0][1] * gammaU.xy + DLUL[2][0][1] * gammaU.xz,
	DLUL[0][0][2] * gammaU.xx + DLUL[1][0][2] * gammaU.xy + DLUL[2][0][2] * gammaU.xz,
	},{DLUL[0][1][0] * gammaU.xx + DLUL[1][1][0] * gammaU.xy + DLUL[2][1][0] * gammaU.xz,
	DLUL[0][1][1] * gammaU.xx + DLUL[1][1][1] * gammaU.xy + DLUL[2][1][1] * gammaU.xz,
	DLUL[0][1][2] * gammaU.xx + DLUL[1][1][2] * gammaU.xy + DLUL[2][1][2] * gammaU.xz,
	},{DLUL[0][2][0] * gammaU.xx + DLUL[1][2][0] * gammaU.xy + DLUL[2][2][0] * gammaU.xz,
	DLUL[0][2][1] * gammaU.xx + DLUL[1][2][1] * gammaU.xy + DLUL[2][2][1] * gammaU.xz,
	DLUL[0][2][2] * gammaU.xx + DLUL[1][2][2] * gammaU.xy + DLUL[2][2][2] * gammaU.xz,
	},},{{DLUL[0][0][0] * gammaU.xy + DLUL[1][0][0] * gammaU.yy + DLUL[2][0][0] * gammaU.yz,
	DLUL[0][0][1] * gammaU.xy + DLUL[1][0][1] * gammaU.yy + DLUL[2][0][1] * gammaU.yz,
	DLUL[0][0][2] * gammaU.xy + DLUL[1][0][2] * gammaU.yy + DLUL[2][0][2] * gammaU.yz,
	},{DLUL[0][1][0] * gammaU.xy + DLUL[1][1][0] * gammaU.yy + DLUL[2][1][0] * gammaU.yz,
	DLUL[0][1][1] * gammaU.xy + DLUL[1][1][1] * gammaU.yy + DLUL[2][1][1] * gammaU.yz,
	DLUL[0][1][2] * gammaU.xy + DLUL[1][1][2] * gammaU.yy + DLUL[2][1][2] * gammaU.yz,
	},{DLUL[0][2][0] * gammaU.xy + DLUL[1][2][0] * gammaU.yy + DLUL[2][2][0] * gammaU.yz,
	DLUL[0][2][1] * gammaU.xy + DLUL[1][2][1] * gammaU.yy + DLUL[2][2][1] * gammaU.yz,
	DLUL[0][2][2] * gammaU.xy + DLUL[1][2][2] * gammaU.yy + DLUL[2][2][2] * gammaU.yz,
	},},{{DLUL[0][0][0] * gammaU.xz + DLUL[1][0][0] * gammaU.yz + DLUL[2][0][0] * gammaU.zz,
	DLUL[0][0][1] * gammaU.xz + DLUL[1][0][1] * gammaU.yz + DLUL[2][0][1] * gammaU.zz,
	DLUL[0][0][2] * gammaU.xz + DLUL[1][0][2] * gammaU.yz + DLUL[2][0][2] * gammaU.zz,
	},{DLUL[0][1][0] * gammaU.xz + DLUL[1][1][0] * gammaU.yz + DLUL[2][1][0] * gammaU.zz,
	DLUL[0][1][1] * gammaU.xz + DLUL[1][1][1] * gammaU.yz + DLUL[2][1][1] * gammaU.zz,
	DLUL[0][1][2] * gammaU.xz + DLUL[1][1][2] * gammaU.yz + DLUL[2][1][2] * gammaU.zz,
	},{DLUL[0][2][0] * gammaU.xz + DLUL[1][2][0] * gammaU.yz + DLUL[2][2][0] * gammaU.zz,
	DLUL[0][2][1] * gammaU.xz + DLUL[1][2][1] * gammaU.yz + DLUL[2][2][1] * gammaU.zz,
	DLUL[0][2][2] * gammaU.xz + DLUL[1][2][2] * gammaU.yz + DLUL[2][2][2] * gammaU.zz,
	},},};
	real D12SymLL[6] = {
	U->d[0].xx * DUUL[0][0][0] + U->d[0].xy * DUUL[0][1][0] + U->d[0].xz * DUUL[0][2][0] + U->d[1].xx * DUUL[1][0][0] + U->d[1].xy * DUUL[1][1][0] + U->d[1].xz * DUUL[1][2][0] + U->d[2].xx * DUUL[2][0][0] + U->d[2].xy * DUUL[2][1][0] + U->d[2].xz * DUUL[2][2][0],
	U->d[0].xy * DUUL[0][0][0] + U->d[0].yy * DUUL[0][1][0] + U->d[0].yz * DUUL[0][2][0] + U->d[1].xy * DUUL[1][0][0] + U->d[1].yy * DUUL[1][1][0] + U->d[1].yz * DUUL[1][2][0] + U->d[2].xy * DUUL[2][0][0] + U->d[2].yy * DUUL[2][1][0] + U->d[2].yz * DUUL[2][2][0],
	U->d[0].xz * DUUL[0][0][0] + U->d[0].yz * DUUL[0][1][0] + U->d[0].zz * DUUL[0][2][0] + U->d[1].xz * DUUL[1][0][0] + U->d[1].yz * DUUL[1][1][0] + U->d[1].zz * DUUL[1][2][0] + U->d[2].xz * DUUL[2][0][0] + U->d[2].yz * DUUL[2][1][0] + U->d[2].zz * DUUL[2][2][0],
	U->d[0].xy * DUUL[0][0][1] + U->d[0].yy * DUUL[0][1][1] + U->d[0].yz * DUUL[0][2][1] + U->d[1].xy * DUUL[1][0][1] + U->d[1].yy * DUUL[1][1][1] + U->d[1].yz * DUUL[1][2][1] + U->d[2].xy * DUUL[2][0][1] + U->d[2].yy * DUUL[2][1][1] + U->d[2].yz * DUUL[2][2][1],
	U->d[0].xz * DUUL[0][0][1] + U->d[0].yz * DUUL[0][1][1] + U->d[0].zz * DUUL[0][2][1] + U->d[1].xz * DUUL[1][0][1] + U->d[1].yz * DUUL[1][1][1] + U->d[1].zz * DUUL[1][2][1] + U->d[2].xz * DUUL[2][0][1] + U->d[2].yz * DUUL[2][1][1] + U->d[2].zz * DUUL[2][2][1],
	U->d[0].xz * DUUL[0][0][2] + U->d[0].yz * DUUL[0][1][2] + U->d[0].zz * DUUL[0][2][2] + U->d[1].xz * DUUL[1][0][2] + U->d[1].yz * DUUL[1][1][2] + U->d[1].zz * DUUL[1][2][2] + U->d[2].xz * DUUL[2][0][2] + U->d[2].yz * DUUL[2][1][2] + U->d[2].zz * DUUL[2][2][2],
	};
	real GammaLSymLL[3][6] = {
	{U->d[0].xx,
	U->d[1].xx,
	U->d[2].xx,
	((2 * U->d[1].xy) - U->d[0].yy),
	(U->d[2].xy + (U->d[1].xz - U->d[0].yz)),
	((2 * U->d[2].xz) - U->d[0].zz),
	},{((2 * U->d[0].xy) - U->d[1].xx),
	U->d[0].yy,
	(U->d[2].xy + (U->d[0].yz - U->d[1].xz)),
	U->d[1].yy,
	U->d[2].yy,
	((2 * U->d[2].yz) - U->d[1].zz),
	},{((2 * U->d[0].xz) - U->d[2].xx),
	(U->d[1].xz + (U->d[0].yz - U->d[2].xy)),
	U->d[0].zz,
	((2 * U->d[1].yz) - U->d[2].yy),
	U->d[1].zz,
	U->d[2].zz,
	},};
	real GammaUSymLL[3][6] = {
	{gammaU.xx * GammaLSymLL[0][0] + gammaU.xy * GammaLSymLL[1][0] + gammaU.xz * GammaLSymLL[2][0],
	gammaU.xx * GammaLSymLL[0][1] + gammaU.xy * GammaLSymLL[1][1] + gammaU.xz * GammaLSymLL[2][1],
	gammaU.xx * GammaLSymLL[0][2] + gammaU.xy * GammaLSymLL[1][2] + gammaU.xz * GammaLSymLL[2][2],
	gammaU.xx * GammaLSymLL[0][3] + gammaU.xy * GammaLSymLL[1][3] + gammaU.xz * GammaLSymLL[2][3],
	gammaU.xx * GammaLSymLL[0][4] + gammaU.xy * GammaLSymLL[1][4] + gammaU.xz * GammaLSymLL[2][4],
	gammaU.xx * GammaLSymLL[0][5] + gammaU.xy * GammaLSymLL[1][5] + gammaU.xz * GammaLSymLL[2][5],
	},{gammaU.xy * GammaLSymLL[0][0] + gammaU.yy * GammaLSymLL[1][0] + gammaU.yz * GammaLSymLL[2][0],
	gammaU.xy * GammaLSymLL[0][1] + gammaU.yy * GammaLSymLL[1][1] + gammaU.yz * GammaLSymLL[2][1],
	gammaU.xy * GammaLSymLL[0][2] + gammaU.yy * GammaLSymLL[1][2] + gammaU.yz * GammaLSymLL[2][2],
	gammaU.xy * GammaLSymLL[0][3] + gammaU.yy * GammaLSymLL[1][3] + gammaU.yz * GammaLSymLL[2][3],
	gammaU.xy * GammaLSymLL[0][4] + gammaU.yy * GammaLSymLL[1][4] + gammaU.yz * GammaLSymLL[2][4],
	gammaU.xy * GammaLSymLL[0][5] + gammaU.yy * GammaLSymLL[1][5] + gammaU.yz * GammaLSymLL[2][5],
	},{gammaU.xz * GammaLSymLL[0][0] + gammaU.yz * GammaLSymLL[1][0] + gammaU.zz * GammaLSymLL[2][0],
	gammaU.xz * GammaLSymLL[0][1] + gammaU.yz * GammaLSymLL[1][1] + gammaU.zz * GammaLSymLL[2][1],
	gammaU.xz * GammaLSymLL[0][2] + gammaU.yz * GammaLSymLL[1][2] + gammaU.zz * GammaLSymLL[2][2],
	gammaU.xz * GammaLSymLL[0][3] + gammaU.yz * GammaLSymLL[1][3] + gammaU.zz * GammaLSymLL[2][3],
	gammaU.xz * GammaLSymLL[0][4] + gammaU.yz * GammaLSymLL[1][4] + gammaU.zz * GammaLSymLL[2][4],
	gammaU.xz * GammaLSymLL[0][5] + gammaU.yz * GammaLSymLL[1][5] + gammaU.zz * GammaLSymLL[2][5],
	},};
	real Gamma3L[3] = {
	GammaUSymLL[0][0] + GammaUSymLL[1][1] + GammaUSymLL[2][2],
	GammaUSymLL[0][1] + GammaUSymLL[1][3] + GammaUSymLL[2][4],
	GammaUSymLL[0][2] + GammaUSymLL[1][4] + GammaUSymLL[2][5],
	};
	real Gamma31SymLL[6] = {
	Gamma3L[0] * GammaUSymLL[0][0] + Gamma3L[1] * GammaUSymLL[1][0] + Gamma3L[2] * GammaUSymLL[2][0],
	Gamma3L[0] * GammaUSymLL[0][1] + Gamma3L[1] * GammaUSymLL[1][1] + Gamma3L[2] * GammaUSymLL[2][1],
	Gamma3L[0] * GammaUSymLL[0][2] + Gamma3L[1] * GammaUSymLL[1][2] + Gamma3L[2] * GammaUSymLL[2][2],
	Gamma3L[0] * GammaUSymLL[0][3] + Gamma3L[1] * GammaUSymLL[1][3] + Gamma3L[2] * GammaUSymLL[2][3],
	Gamma3L[0] * GammaUSymLL[0][4] + Gamma3L[1] * GammaUSymLL[1][4] + Gamma3L[2] * GammaUSymLL[2][4],
	Gamma3L[0] * GammaUSymLL[0][5] + Gamma3L[1] * GammaUSymLL[1][5] + Gamma3L[2] * GammaUSymLL[2][5],
	};
	real GammaLUL[3][3][3] = {
	{{gammaU.xx * GammaLSymLL[0][0] + gammaU.xy * GammaLSymLL[0][1] + gammaU.xz * GammaLSymLL[0][2],
	gammaU.xx * GammaLSymLL[0][1] + gammaU.xy * GammaLSymLL[0][3] + gammaU.xz * GammaLSymLL[0][4],
	gammaU.xx * GammaLSymLL[0][2] + gammaU.xy * GammaLSymLL[0][4] + gammaU.xz * GammaLSymLL[0][5],
	},{gammaU.xy * GammaLSymLL[0][0] + gammaU.yy * GammaLSymLL[0][1] + gammaU.yz * GammaLSymLL[0][2],
	gammaU.xy * GammaLSymLL[0][1] + gammaU.yy * GammaLSymLL[0][3] + gammaU.yz * GammaLSymLL[0][4],
	gammaU.xy * GammaLSymLL[0][2] + gammaU.yy * GammaLSymLL[0][4] + gammaU.yz * GammaLSymLL[0][5],
	},{gammaU.xz * GammaLSymLL[0][0] + gammaU.yz * GammaLSymLL[0][1] + gammaU.zz * GammaLSymLL[0][2],
	gammaU.xz * GammaLSymLL[0][1] + gammaU.yz * GammaLSymLL[0][3] + gammaU.zz * GammaLSymLL[0][4],
	gammaU.xz * GammaLSymLL[0][2] + gammaU.yz * GammaLSymLL[0][4] + gammaU.zz * GammaLSymLL[0][5],
	},},{{gammaU.xx * GammaLSymLL[1][0] + gammaU.xy * GammaLSymLL[1][1] + gammaU.xz * GammaLSymLL[1][2],
	gammaU.xx * GammaLSymLL[1][1] + gammaU.xy * GammaLSymLL[1][3] + gammaU.xz * GammaLSymLL[1][4],
	gammaU.xx * GammaLSymLL[1][2] + gammaU.xy * GammaLSymLL[1][4] + gammaU.xz * GammaLSymLL[1][5],
	},{gammaU.xy * GammaLSymLL[1][0] + gammaU.yy * GammaLSymLL[1][1] + gammaU.yz * GammaLSymLL[1][2],
	gammaU.xy * GammaLSymLL[1][1] + gammaU.yy * GammaLSymLL[1][3] + gammaU.yz * GammaLSymLL[1][4],
	gammaU.xy * GammaLSymLL[1][2] + gammaU.yy * GammaLSymLL[1][4] + gammaU.yz * GammaLSymLL[1][5],
	},{gammaU.xz * GammaLSymLL[1][0] + gammaU.yz * GammaLSymLL[1][1] + gammaU.zz * GammaLSymLL[1][2],
	gammaU.xz * GammaLSymLL[1][1] + gammaU.yz * GammaLSymLL[1][3] + gammaU.zz * GammaLSymLL[1][4],
	gammaU.xz * GammaLSymLL[1][2] + gammaU.yz * GammaLSymLL[1][4] + gammaU.zz * GammaLSymLL[1][5],
	},},{{gammaU.xx * GammaLSymLL[2][0] + gammaU.xy * GammaLSymLL[2][1] + gammaU.xz * GammaLSymLL[2][2],
	gammaU.xx * GammaLSymLL[2][1] + gammaU.xy * GammaLSymLL[2][3] + gammaU.xz * GammaLSymLL[2][4],
	gammaU.xx * GammaLSymLL[2][2] + gammaU.xy * GammaLSymLL[2][4] + gammaU.xz * GammaLSymLL[2][5],
	},{gammaU.xy * GammaLSymLL[2][0] + gammaU.yy * GammaLSymLL[2][1] + gammaU.yz * GammaLSymLL[2][2],
	gammaU.xy * GammaLSymLL[2][1] + gammaU.yy * GammaLSymLL[2][3] + gammaU.yz * GammaLSymLL[2][4],
	gammaU.xy * GammaLSymLL[2][2] + gammaU.yy * GammaLSymLL[2][4] + gammaU.yz * GammaLSymLL[2][5],
	},{gammaU.xz * GammaLSymLL[2][0] + gammaU.yz * GammaLSymLL[2][1] + gammaU.zz * GammaLSymLL[2][2],
	gammaU.xz * GammaLSymLL[2][1] + gammaU.yz * GammaLSymLL[2][3] + gammaU.zz * GammaLSymLL[2][4],
	gammaU.xz * GammaLSymLL[2][2] + gammaU.yz * GammaLSymLL[2][4] + gammaU.zz * GammaLSymLL[2][5],
	},},};
	real GammaLSymUU[3][6] = {
	{gammaU.xx * GammaLUL[0][0][0] + gammaU.xy * GammaLUL[0][0][1] + gammaU.xz * GammaLUL[0][0][2],
	gammaU.xy * GammaLUL[0][0][0] + gammaU.yy * GammaLUL[0][0][1] + gammaU.yz * GammaLUL[0][0][2],
	gammaU.xz * GammaLUL[0][0][0] + gammaU.yz * GammaLUL[0][0][1] + gammaU.zz * GammaLUL[0][0][2],
	gammaU.xy * GammaLUL[0][1][0] + gammaU.yy * GammaLUL[0][1][1] + gammaU.yz * GammaLUL[0][1][2],
	gammaU.xz * GammaLUL[0][1][0] + gammaU.yz * GammaLUL[0][1][1] + gammaU.zz * GammaLUL[0][1][2],
	gammaU.xz * GammaLUL[0][2][0] + gammaU.yz * GammaLUL[0][2][1] + gammaU.zz * GammaLUL[0][2][2],
	},{gammaU.xx * GammaLUL[1][0][0] + gammaU.xy * GammaLUL[1][0][1] + gammaU.xz * GammaLUL[1][0][2],
	gammaU.xy * GammaLUL[1][0][0] + gammaU.yy * GammaLUL[1][0][1] + gammaU.yz * GammaLUL[1][0][2],
	gammaU.xz * GammaLUL[1][0][0] + gammaU.yz * GammaLUL[1][0][1] + gammaU.zz * GammaLUL[1][0][2],
	gammaU.xy * GammaLUL[1][1][0] + gammaU.yy * GammaLUL[1][1][1] + gammaU.yz * GammaLUL[1][1][2],
	gammaU.xz * GammaLUL[1][1][0] + gammaU.yz * GammaLUL[1][1][1] + gammaU.zz * GammaLUL[1][1][2],
	gammaU.xz * GammaLUL[1][2][0] + gammaU.yz * GammaLUL[1][2][1] + gammaU.zz * GammaLUL[1][2][2],
	},{gammaU.xx * GammaLUL[2][0][0] + gammaU.xy * GammaLUL[2][0][1] + gammaU.xz * GammaLUL[2][0][2],
	gammaU.xy * GammaLUL[2][0][0] + gammaU.yy * GammaLUL[2][0][1] + gammaU.yz * GammaLUL[2][0][2],
	gammaU.xz * GammaLUL[2][0][0] + gammaU.yz * GammaLUL[2][0][1] + gammaU.zz * GammaLUL[2][0][2],
	gammaU.xy * GammaLUL[2][1][0] + gammaU.yy * GammaLUL[2][1][1] + gammaU.yz * GammaLUL[2][1][2],
	gammaU.xz * GammaLUL[2][1][0] + gammaU.yz * GammaLUL[2][1][1] + gammaU.zz * GammaLUL[2][1][2],
	gammaU.xz * GammaLUL[2][2][0] + gammaU.yz * GammaLUL[2][2][1] + gammaU.zz * GammaLUL[2][2][2],
	},};
	real Gamma11SymLL[6] = {
	GammaLSymLL[0][0] * GammaLSymUU[0][0] + GammaLSymLL[0][1] * GammaLSymUU[0][1] + GammaLSymLL[0][2] * GammaLSymUU[0][2] + GammaLSymLL[0][1] * GammaLSymUU[0][1] + GammaLSymLL[0][3] * GammaLSymUU[0][3] + GammaLSymLL[0][4] * GammaLSymUU[0][4] + GammaLSymLL[0][2] * GammaLSymUU[0][2] + GammaLSymLL[0][4] * GammaLSymUU[0][4] + GammaLSymLL[0][5] * GammaLSymUU[0][5],
	GammaLSymLL[0][0] * GammaLSymUU[1][0] + GammaLSymLL[0][1] * GammaLSymUU[1][1] + GammaLSymLL[0][2] * GammaLSymUU[1][2] + GammaLSymLL[0][1] * GammaLSymUU[1][1] + GammaLSymLL[0][3] * GammaLSymUU[1][3] + GammaLSymLL[0][4] * GammaLSymUU[1][4] + GammaLSymLL[0][2] * GammaLSymUU[1][2] + GammaLSymLL[0][4] * GammaLSymUU[1][4] + GammaLSymLL[0][5] * GammaLSymUU[1][5],
	GammaLSymLL[0][0] * GammaLSymUU[2][0] + GammaLSymLL[0][1] * GammaLSymUU[2][1] + GammaLSymLL[0][2] * GammaLSymUU[2][2] + GammaLSymLL[0][1] * GammaLSymUU[2][1] + GammaLSymLL[0][3] * GammaLSymUU[2][3] + GammaLSymLL[0][4] * GammaLSymUU[2][4] + GammaLSymLL[0][2] * GammaLSymUU[2][2] + GammaLSymLL[0][4] * GammaLSymUU[2][4] + GammaLSymLL[0][5] * GammaLSymUU[2][5],
	GammaLSymLL[1][0] * GammaLSymUU[1][0] + GammaLSymLL[1][1] * GammaLSymUU[1][1] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][1] * GammaLSymUU[1][1] + GammaLSymLL[1][3] * GammaLSymUU[1][3] + GammaLSymLL[1][4] * GammaLSymUU[1][4] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][4] * GammaLSymUU[1][4] + GammaLSymLL[1][5] * GammaLSymUU[1][5],
	GammaLSymLL[1][0] * GammaLSymUU[2][0] + GammaLSymLL[1][1] * GammaLSymUU[2][1] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][1] * GammaLSymUU[2][1] + GammaLSymLL[1][3] * GammaLSymUU[2][3] + GammaLSymLL[1][4] * GammaLSymUU[2][4] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][4] * GammaLSymUU[2][4] + GammaLSymLL[1][5] * GammaLSymUU[2][5],
	GammaLSymLL[2][0] * GammaLSymUU[2][0] + GammaLSymLL[2][1] * GammaLSymUU[2][1] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][1] * GammaLSymUU[2][1] + GammaLSymLL[2][3] * GammaLSymUU[2][3] + GammaLSymLL[2][4] * GammaLSymUU[2][4] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][4] * GammaLSymUU[2][4] + GammaLSymLL[2][5] * GammaLSymUU[2][5],
	};
	real ADL[3] = {
	U->a.x - 2 * D3L[0],
	U->a.y - 2 * D3L[1],
	U->a.z - 2 * D3L[2],
	};
	real ADU[3] = {
	gammaU.xx * ADL[0] + gammaU.xy * ADL[1] + gammaU.xz * ADL[2],
	gammaU.xy * ADL[0] + gammaU.yy * ADL[1] + gammaU.yz * ADL[2],
	gammaU.xz * ADL[0] + gammaU.yz * ADL[1] + gammaU.zz * ADL[2],
	};
	real ADDSymLL[6] = {
	ADU[0] * (2 * U->d[0].xx) + ADU[1] * (2 * U->d[0].xy) + ADU[2] * (2 * U->d[0].xz),
	ADU[0] * (U->d[0].xy + U->d[1].xx) + ADU[1] * (U->d[0].yy + U->d[1].xy) + ADU[2] * (U->d[0].yz + U->d[1].xz),
	ADU[0] * (U->d[0].xz + U->d[2].xx) + ADU[1] * (U->d[0].yz + U->d[2].xy) + ADU[2] * (U->d[0].zz + U->d[2].xz),
	ADU[0] * (2 * U->d[1].xy) + ADU[1] * (2 * U->d[1].yy) + ADU[2] * (2 * U->d[1].yz),
	ADU[0] * (U->d[1].xz + U->d[2].xy) + ADU[1] * (U->d[1].yz + U->d[2].yy) + ADU[2] * (U->d[1].zz + U->d[2].yz),
	ADU[0] * (2 * U->d[2].xz) + ADU[1] * (2 * U->d[2].yz) + ADU[2] * (2 * U->d[2].zz),
	};
	real R4SymLL[6] = {
	#if 1
	0,0,0,0,0,0
	#else
	8. * M_PI * ((density + pressure) * vel4_.x * vel4_.x + .5 * (density - pressure) * U->gamma.xx),
	8. * M_PI * ((density + pressure) * vel4_.x * vel4_.y + .5 * (density - pressure) * U->gamma.xy),
	8. * M_PI * ((density + pressure) * vel4_.x * vel4_.z + .5 * (density - pressure) * U->gamma.xz),
	8. * M_PI * ((density + pressure) * vel4_.y * vel4_.y + .5 * (density - pressure) * U->gamma.yy),
	8. * M_PI * ((density + pressure) * vel4_.y * vel4_.z + .5 * (density - pressure) * U->gamma.yz),
	8. * M_PI * ((density + pressure) * vel4_.z * vel4_.z + .5 * (density - pressure) * U->gamma.zz),
	#endif
	};
	real SSymLL[6] = {
	-R4SymLL[0] + trK * U->K.xx - 2 * KSqSymLL.s[0] + 4 * D12SymLL[0] + Gamma31SymLL[0] - Gamma11SymLL[0] + ADDSymLL[0] + (U->a.x * ((2 * U->V.x) - D1L[0])),
	-R4SymLL[1] + trK * U->K.xy - 2 * KSqSymLL.s[1] + 4 * D12SymLL[1] + Gamma31SymLL[1] - Gamma11SymLL[1] + ADDSymLL[1] + ((((2 * U->a.y * U->V.x) - (U->a.y * D1L[0])) + ((2 * U->a.x * U->V.y) - (U->a.x * D1L[1]))) / 2),
	-R4SymLL[2] + trK * U->K.xz - 2 * KSqSymLL.s[2] + 4 * D12SymLL[2] + Gamma31SymLL[2] - Gamma11SymLL[2] + ADDSymLL[2] + ((((2 * U->a.z * U->V.x) - (U->a.z * D1L[0])) + ((2 * U->a.x * U->V.z) - (U->a.x * D1L[2]))) / 2),
	-R4SymLL[3] + trK * U->K.yy - 2 * KSqSymLL.s[3] + 4 * D12SymLL[3] + Gamma31SymLL[3] - Gamma11SymLL[3] + ADDSymLL[3] + (U->a.y * ((2 * U->V.y) - D1L[1])),
	-R4SymLL[4] + trK * U->K.yz - 2 * KSqSymLL.s[4] + 4 * D12SymLL[4] + Gamma31SymLL[4] - Gamma11SymLL[4] + ADDSymLL[4] + ((((2 * U->a.z * U->V.y) - (U->a.z * D1L[1])) + ((2 * U->a.y * U->V.z) - (U->a.y * D1L[2]))) / 2),
	-R4SymLL[5] + trK * U->K.zz - 2 * KSqSymLL.s[5] + 4 * D12SymLL[5] + Gamma31SymLL[5] - Gamma11SymLL[5] + ADDSymLL[5] + (U->a.z * ((2 * U->V.z) - D1L[2])),
	};
	real GU0L[3] = {
	#if 1
	0,0,0
	#else
	8. * M_PI * ((density + pressure) * vel4_.w * vel4_.x + pressure * beta_.x),
	8. * M_PI * ((density + pressure) * vel4_.w * vel4_.y + pressure * beta_.y),
	8. * M_PI * ((density + pressure) * vel4_.w * vel4_.z + pressure * beta_.z),
	#endif
	};
	real AKL[3] = {
	U->a.x * KUL.x.x + U->a.y * KUL.y.x + U->a.z * KUL.z.x,
	U->a.x * KUL.x.y + U->a.y * KUL.y.y + U->a.z * KUL.z.y,
	U->a.x * KUL.x.z + U->a.y * KUL.y.z + U->a.z * KUL.z.z,
	};
	real K12D23L[3] = {
	KUL.x.x * DLUL[0][0][0] +KUL.x.y * DLUL[0][1][0] +KUL.x.z * DLUL[0][2][0] + KUL.y.x * DLUL[0][0][1] +KUL.y.y * DLUL[0][1][1] +KUL.y.z * DLUL[0][2][1] + KUL.z.x * DLUL[0][0][2] +KUL.z.y * DLUL[0][1][2] +KUL.z.z * DLUL[0][2][2],
	KUL.x.x * DLUL[1][0][0] +KUL.x.y * DLUL[1][1][0] +KUL.x.z * DLUL[1][2][0] + KUL.y.x * DLUL[1][0][1] +KUL.y.y * DLUL[1][1][1] +KUL.y.z * DLUL[1][2][1] + KUL.z.x * DLUL[1][0][2] +KUL.z.y * DLUL[1][1][2] +KUL.z.z * DLUL[1][2][2],
	KUL.x.x * DLUL[2][0][0] +KUL.x.y * DLUL[2][1][0] +KUL.x.z * DLUL[2][2][0] + KUL.y.x * DLUL[2][0][1] +KUL.y.y * DLUL[2][1][1] +KUL.y.z * DLUL[2][2][1] + KUL.z.x * DLUL[2][0][2] +KUL.z.y * DLUL[2][1][2] +KUL.z.z * DLUL[2][2][2],
	};
	real KD23L[3] = {
	KUL.x.x * D1L[0] + KUL.y.x * D1L[1] + KUL.z.x * D1L[2],
	KUL.x.y * D1L[0] + KUL.y.y * D1L[1] + KUL.z.y * D1L[2],
	KUL.x.z * D1L[0] + KUL.y.z * D1L[1] + KUL.z.z * D1L[2],
	};
	real K12D12L[3] = {
	KUL.x.x * DLUL[0][0][0] + KUL.x.y * DLUL[0][1][0] + KUL.x.z * DLUL[0][2][0] + KUL.y.x * DLUL[1][0][0] + KUL.y.y * DLUL[1][1][0] + KUL.y.z * DLUL[1][2][0] + KUL.z.x * DLUL[2][0][0] + KUL.z.y * DLUL[2][1][0] + KUL.z.z * DLUL[2][2][0],
	KUL.x.x * DLUL[0][0][1] + KUL.x.y * DLUL[0][1][1] + KUL.x.z * DLUL[0][2][1] + KUL.y.x * DLUL[1][0][1] + KUL.y.y * DLUL[1][1][1] + KUL.y.z * DLUL[1][2][1] + KUL.z.x * DLUL[2][0][1] + KUL.z.y * DLUL[2][1][1] + KUL.z.z * DLUL[2][2][1],
	KUL.x.x * DLUL[0][0][2] + KUL.x.y * DLUL[0][1][2] + KUL.x.z * DLUL[0][2][2] + KUL.y.x * DLUL[1][0][2] + KUL.y.y * DLUL[1][1][2] + KUL.y.z * DLUL[1][2][2] + KUL.z.x * DLUL[2][0][2] + KUL.z.y * DLUL[2][1][2] + KUL.z.z * DLUL[2][2][2],
	};
	real KD12L[3] = {
	KUL.x.x * D3L[0] + KUL.y.x * D3L[1] + KUL.z.x * D3L[2],
	KUL.x.y * D3L[0] + KUL.y.y * D3L[1] + KUL.z.y * D3L[2],
	KUL.x.z * D3L[0] + KUL.y.z * D3L[1] + KUL.z.z * D3L[2],
	};
	real PL[3] = {
	GU0L[0] + AKL[0] - U->a.x * trK + K12D23L[0] + KD23L[0] - 2 * K12D12L[0] + 2 * KD12L[0],
	GU0L[1] + AKL[1] - U->a.y * trK + K12D23L[1] + KD23L[1] - 2 * K12D12L[1] + 2 * KD12L[1],
	GU0L[2] + AKL[2] - U->a.z * trK + K12D23L[2] + KD23L[2] - 2 * K12D12L[2] + 2 * KD12L[2],
	};

	/*
	TODO now that alpha and gamma are moved from the flux eigensystem
	how would we still incorporate the shift terms with them?
	the easy solution is to re-add them back in (even though shift computation is supposed to go on apart from the eigensystem)
	typically, when they are included, they get a flux between cells equal to the shift at that cell
	maybe that has to be done, even if they are not a part of the eigensystem?
	then again, maybe I should be keeping alpha and gamma in the eigensystem,
	and maybe aux vars like the density should be in the eigensystem as well,
	all for no reason more than to be influenced by the shift
	
	Then again, maybe this is an argument for the solver to specify the flux vector size 
	-- especially if it is allowed a custom RoeFluxDeriv function.
	*/
	
	deriv->alpha += -U->alpha * U->alpha * f * trK;
	deriv->gamma.xx += -2. * U->alpha * U->K.xx;
	deriv->gamma.xy += -2. * U->alpha * U->K.xy;
	deriv->gamma.xz += -2. * U->alpha * U->K.xz;
	deriv->gamma.yy += -2. * U->alpha * U->K.yy;
	deriv->gamma.yz += -2. * U->alpha * U->K.yz;
	deriv->gamma.zz += -2. * U->alpha * U->K.zz;

	deriv->K.xx += U->alpha * SSymLL[0];
	deriv->K.xy += U->alpha * SSymLL[1];
	deriv->K.xz += U->alpha * SSymLL[2];
	deriv->K.yy += U->alpha * SSymLL[3];
	deriv->K.yz += U->alpha * SSymLL[4];
	deriv->K.zz += U->alpha * SSymLL[5];
	deriv->V.x += U->alpha * PL[0];
	deriv->V.y += U->alpha * PL[1];
	deriv->V.z += U->alpha * PL[2];
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
		real d1 = sym3_dot(U->d[<?=i-1?>], gammaU);
		real d2 = 0.<?
	for j=1,3 do
		for k,xk in ipairs(xNames) do
?> + U->d[<?=j-1?>].<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
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
?>	U->d[<?=i-1?>].<?=xjk?> += (
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
