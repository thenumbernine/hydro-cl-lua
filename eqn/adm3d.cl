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
	eig->sqrt_gammaUjj = _real3(sqrt(eig->gammaU.xx), sqrt(eig->gammaU.yy), sqrt(eig->gammaU.zz));
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

<? if not eqn.noZeroRowsInFlux then ?>
	wave[0] = -lambdaGauge;
	<? for i=1,5 do ?> wave[<?=i?>] = -lambdaLight; <? end ?>
	<? for i=6,23 do ?> wave[<?=i?>] = 0.; <? end ?>
	<? for i=24,28 do ?> wave[<?=i?>] = lambdaLight; <? end ?>
	wave[29] = lambdaGauge;
<? else ?>
	wave[0] = -lambdaGauge;
	<? for i=1,5 do ?> wave[<?=i?>] = -lambdaLight; <? end ?>
	wave[6] = 0.;
	<? for i=7,11 do ?> wave[<?=i?>] = lambdaLight; <? end ?>
	wave[12] = lambdaGauge;
<? end ?>
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

<? if not eqn.noZeroRowsInFlux then ?>

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

	<? end -- side ?>

<? else -- eqn.noZeroRowsInFlux ?>

	real _1_sqrt_f = 1. / eig->sqrt_f;
	real _1_f = _1_sqrt_f * _1_sqrt_f; 
	real sqrt_gammaUjj = eig->sqrt_gammaUjj.s<?=side?>;
	real _1_gammaUjj = 1. / eig->gammaU.s<?=side?><?=side?>;

	real a = inputU->a.s<?=side?>;

	//now swap x and side on the sym3's
	sym3 d = sym3_swap<?=side?>(inputU->d[<?=side?>]);
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
	resultU->d[<?=side?>] = sym3_swap<?=side?>(d);
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
<? else -- noZeroRowsInFlux ?>

	<?=addr2?> const <?=eqn.cons_t?>* inputU = (<?=addr2?> const <?=eqn.cons_t?>*)input;
	<?=addr0?> <?=eqn.cons_t?>* resultU = (<?=addr0?> <?=eqn.cons_t?>*)results;
	
	real f = eig->sqrt_f * eig->sqrt_f;
	
	//now swap x and side on the sym3's
	sym3 input_d = sym3_swap<?=side?>(inputU->d[<?=side?>]);
	sym3 input_K = sym3_swap<?=side?>(inputU->K);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);

	resultU->a.s<?=side?> = sym3_dot(input_K, gammaU) * eig->alpha * f;
	sym3 result_d = sym3_scale(input_K, eig->alpha);
	sym3 result_K = sym3_scale(input_d, eig->alpha * gammaU.xx);
	result_K.xx += (inputU->a.s<?=side?> - sym3_dot(input_d, gammaU)) * eig->alpha;

	//now swap x and side on the sym3's
	resultU->d[<?=side?>] = sym3_swap<?=side?>(result_d);
	resultU->K = sym3_swap<?=side?>(result_K);

<? end -- noZeroRowsInFlux ?>
}
<?
	end
end
?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
<? if not eqn.noZeroRowsInFlux then ?>
	// !noZeroRowsInFlux doesn't have fluxTransform implemented

	<?=eqn.eigen_t?> eig;
	eigen_forCell_<?=side?>local(&eig, &U, x);

	real wave[numWaves];
	eigen_calcWaves_<?=side?>__(wave, &eig, x);

	real charvars[numWaves];
	eigen_leftTransform_<?=side?>___(charvars, &eig, U.ptr, x);
	
	for (int j = 0; j < numWaves; ++j) {
		charvars[j] *= wave[j];
	}
			
	<?=eqn.cons_t?> F;
	eigen_rightTransform_<?=side?>___(F.ptr, &eig, charvars, x);

<? else	-- noZeroRowsInFlux ?>
	// noZeroRowsInFlux has fluxTransform implemented

	<?=eqn.eigen_t?> eig;
	eigen_forCell_<?=side?>local(&eig, &U, x);

	<?=eqn.cons_t?> F;
	eigen_fluxTransform_<?=side?>___(F.ptr, &eig, U.ptr, x);

<? end -- noZeroRowsInFlux ?>

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
?>		sym3_sym3_mul(U->d[<?=i-1?>], gammaU),
<? end
?>	};

	//d_ull = d^i_jk = gamma^il d_ljk
	sym3 d_ull[3] = {
<? for i,xi in ipairs(xNames) do
?>		(sym3){
<?	for jk,xjk in ipairs(symNames) do
?>			.<?=xjk?> = 0.<?
		for l,xl in ipairs(xNames) do
?> + gammaU.<?=sym(i,l)?> * U->d[<?=l-1?>].<?=xjk?><?
		end ?>,
<? end
?>		},
<? end
?>	};

	//d3_l = d^j_ji
	real3 d3_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		?> + d_ull[<?=j-1?>].<?=sym(j,i)?><?
	end	?>,
<? end
?>	};

	real3 a_V_d3_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = U->a.<?=xi?> + U->V.<?=xi?> - d3_l.<?=xi?>,
<? end
?>	};
	real3 a_V_d3_u = sym3_real3_mul(gammaU, a_V_d3_l);

	//srcK_ij = (-a_i a_j 
	//		+ (d_ij^k + d_ji^k - d^k_ij) (a_k + V_k - d^l_lk) 
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
?>			+ (d_llu[<?=i-1?>].<?=xj?>.<?=xk?> - d_llu[<?=j-1?>].<?=xi?>.<?=xk?> - U->d[<?=k-1?>].<?=xij?>) * a_V_d3_u.<?=xk?>
<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_llu[<?=k-1?>].<?=xi?>.<?=xl?> * d_ull[<?=k-1?>].<?=sym(j,l)?>
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
			- 4. * trK * d_ull[<?=j-1?>].<?=sym(j,k)?>
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

<? if eqn.guiVars.linearConstraintCoeff.value ~= 0 then ?>
	// and now for the first-order constraints

	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	<? for i,xi in ipairs(xNames) do ?>{
		<? if i <= solver.dim then ?>
		real di_alpha = (U[stepsize.<?=xi?>].alpha - U[-stepsize.<?=xi?>].alpha) / (2. * grid_dx<?=i-1?>);
		<? else ?>
		real di_alpha = 0.;
		<? end ?>
		deriv->a.<?=xi?> += gui_linearConstraintCoeff * (di_alpha / U->alpha - U->a.<?=xi?>);
	}<? end ?>	
	
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	<? 
for i,xi in ipairs(xNames) do 
	for jk,xjk in ipairs(symNames) do ?>{
		<? if i <= solver.dim then ?>
		real di_gamma_jk = (U[stepsize.<?=xi?>].gamma.<?=xjk?> - U[-stepsize.<?=xi?>].gamma.<?=xjk?>) / (2. * grid_dx<?=i-1?>);
		<? else ?>
		real di_gamma_jk = 0;
		<? end ?>
		deriv->d[<?=i-1?>].<?=xjk?> += gui_linearConstraintCoeff * (.5 * di_gamma_jk - U->d[<?=i-1?>].<?=xjk?>);
	}<? 
	end
end ?>

	//V_i = d_ik^k - d^k_ki <=> V_i += eta (d_ik^k - d^k_ki - V_i)
	deriv->V = real3_add(
		deriv->V,
		real3_scale(
			real3_sub(real3_sub(d1_l, d3_l), U->V),
			gui_linearConstraintCoeff));

	//Kreiss-Oligar diffusion, for stability's sake?
<? end -- eqn.guiVars.linearConstraintCoeff.value  ?>
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
	for j,xj in ipairs(xNames) do
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
