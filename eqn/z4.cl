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
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) is only called once
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
	eig->gamma = U->gamma;
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
	
	wave[0] = -lambdaGauge;
	<? for i=1,6 do ?> wave[<?=i?>] = -lambdaLight; <? end ?>
	<? for i=7,23 do ?> wave[<?=i?>] = 0.; <? end ?>
	<? for i=24,29 do ?> wave[<?=i?>] = lambdaLight; <? end ?>
	wave[30] = lambdaGauge;
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
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gammaU.xx), sqrt(eig.gammaU.yy), sqrt(eig.gammaU.zz));
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

	//input
	real3 a = real3_swap<?=side?>(inputU->a);							//0-2
	sym3 dx = sym3_swap<?=side?>(inputU->d[<?=side?>]);					//3-8
	sym3 dy = sym3_swap<?=side?>(inputU->d[<?=side==1 and 0 or 1?>]);	//9-14
	sym3 dz = sym3_swap<?=side?>(inputU->d[<?=side==2 and 0 or 2?>]);	//15-20
	sym3 K = sym3_swap<?=side?>(inputU->K);								//21-26
	real Theta = inputU->Theta;											//27
	real3 Z = real3_swap<?=side?>(inputU->Z);							//28-30

	//eig
	real sqrt_f = eig->sqrt_f;
	real f = sqrt_f * sqrt_f;	
	sym3 gamma = sym3_swap<?=side?>(eig->gamma);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);
	
	real lambda_1 = (2. - gui_lambda) / (f - 1.);
	real lambda_2 = (2. * f - gui_lambda) / (f - 1.);

	real tr_K = sym3_dot(gammaU, K);

	real aUx = gammaU.xx * a.x + gammaU.xy * a.y + gammaU.xz * a.z;
	real ZUx = gammaU.xx * Z.x + gammaU.xy * Z.y + gammaU.xz * Z.z;

	real result0_a = sqrt_f * (tr_K + lambda_1 * Theta);
	real result0_b = aUx - lambda_2 * (
		  (gammaU.xy * gammaU.xy - gammaU.xx * gammaU.yy) * dx.yy
		+ (gammaU.xy * gammaU.xz - gammaU.xx * gammaU.yz) * 2. * dx.yz
		+ (gammaU.xz * gammaU.xz - gammaU.xx * gammaU.zz) * dx.zz
		
		+ (gammaU.xx * gammaU.yy - gammaU.xy * gammaU.xy) * dy.xy
		+ (gammaU.xx * gammaU.yz - gammaU.xy * gammaU.xz) * dy.xz
		+ (gammaU.xz * gammaU.yy - gammaU.xy * gammaU.yz) * dy.yz
		+ (gammaU.xz * gammaU.yz - gammaU.xy * gammaU.zz) * dy.zz
		
		+ (gammaU.xx * gammaU.yz - gammaU.xy * gammaU.xz) * dz.xy
		+ (gammaU.xx * gammaU.zz - gammaU.xz * gammaU.xz) * dz.xz
		+ (gammaU.xy * gammaU.yz - gammaU.xz * gammaU.yy) * dz.yy
		+ (gammaU.xy * gammaU.zz - gammaU.xz * gammaU.yz) * dz.yz
	
		+ ZUx
	);

	results[0] = result0_a - result0_b;
	results[30] = result0_a + result0_b;

	real results1b = 
		+ (gammaU.xx * gammaU.yy - gammaU.xy * gammaU.xy) * dx.yy
		+ (gammaU.xx * gammaU.yz - gammaU.xy * gammaU.xz) * 2. * dx.yz
		+ (gammaU.xx * gammaU.zz - gammaU.xz * gammaU.xz) * dx.zz
		
		+ (gammaU.xy * gammaU.xy - gammaU.xx * gammaU.yy) * dy.xy
		+ (gammaU.xy * gammaU.xz - gammaU.xx * gammaU.yz) * dy.xz
		+ (gammaU.xy * gammaU.yz - gammaU.xz * gammaU.yy) * dy.yz
		+ (gammaU.xy * gammaU.zz - gammaU.xz * gammaU.yz) * dy.zz
		
		+ (gammaU.xy * gammaU.xz - gammaU.xx * gammaU.yz) * dz.xy
		+ (gammaU.xz * gammaU.xz - gammaU.xx * gammaU.zz) * dz.xz
		+ (gammaU.xz * gammaU.yy - gammaU.xy * gammaU.yz) * dz.yy
		+ (gammaU.xz * gammaU.yz - gammaU.xy * gammaU.zz) * dz.yz
		
		- ZUx;

	results[1] = Theta - results1b;
	results[29] = Theta + results1b;

	real results2b = .5 * (
		a.y
		
		- 2. * gammaU.xy * dx.yy
		- 2. * gammaU.xz * dx.yz
		
		+ gammaU.xx * dy.xx
		+ 2. * gammaU.xy * dy.xy
		+ 2. * gammaU.xz * dy.xz
		- gammaU.yy * dy.yy
		+ gammaU.zz * dy.zz
		
		- 2. * gammaU.yz * dz.yy
		- 2. * gammaU.zz * dz.yz
		
		- 2. * Z.y	
	);

	results[2] = K.xy - results2b;
	results[24] = K.xy + results2b;

	real results3b = .5 * (
		a.z
		
		- 2. * gammaU.xy * dx.yz
		- 2. * gammaU.xz * dx.zz
		
		- 2. * gammaU.yy * dy.yz
		- 2. * gammaU.yz * dy.zz
		
		+ gammaU.xx * dz.xx
		+ 2. * gammaU.xy * dz.xy
		+ 2. * gammaU.xz * dz.xz
		+ gammaU.yy * dz.yy
		- gammaU.zz * dz.zz
		
		- 2. * Z.z
	);
	
	results[3] = K.xz - results3b;
	results[25] = K.xz + results3b;

	real dUx_yy = 
		 gammaU.xx * dx.yy
		+ gammaU.xy * dy.yy
		+ gammaU.xz * dz.yy;
	
	results[4] = K.yy - dUx_yy;
	results[26] = K.yy + dUx_yy;

	real dUx_yz = 
		  gammaU.xx * dx.yz
		+ gammaU.xy * dy.yz
		+ gammaU.xz * dz.yz;

	results[5] = K.yz - dUx_yz;
	results[27] = K.yz + dUx_yz;

	real dUx_zz = 
		 gammaU.xx * dx.zz
		+ gammaU.xy * dy.zz
		+ gammaU.xz * dz.zz;

	results[6] = K.zz - dUx_zz;
	results[28] = K.zz + dUx_zz;
	
	results[7] = a.y;
	results[8] = a.z;
	results[9] = dy.xx;
	results[10] = dy.xy;
	results[11] = dy.xz;
	results[12] = dy.yy;
	results[13] = dy.yz;
	results[14] = dy.zz;
	results[15] = dz.xx;
	results[16] = dz.xy;
	results[17] = dz.xz;
	results[18] = dz.yy;
	results[19] = dz.yz;
	results[20] = dz.zz;

	real fPlus1 = f + 1.;
	real3 d_1 = _real3(
		sym3_dot(dx, gammaU),
		sym3_dot(dy, gammaU),
		sym3_dot(dz, gammaU));

	results[21] = a.x - fPlus1 * d_1.x

		//g^ij d_ijx
		+ gammaU.xx * dx.xx
		+ gammaU.xy * dx.xy
		+ gammaU.xz * dx.xz
		+ gammaU.xy * dy.xx
		+ gammaU.yy * dy.xy
		+ gammaU.yz * dy.xz
		+ gammaU.xz * dz.xx
		+ gammaU.yz * dz.xy
		+ gammaU.zz * dz.xz
		
		+ Z.x
	;
	results[22] = a.y - fPlus1 * d_1.y
	
		//g^ij d_ijy
		+ gammaU.xx * dx.xy
		+ gammaU.xy * dx.yy
		+ gammaU.xz * dx.yz
		+ gammaU.xy * dy.xy
		+ gammaU.yy * dy.yy
		+ gammaU.yz * dy.yz
		+ gammaU.xz * dz.xy
		+ gammaU.yz * dz.yy
		+ gammaU.zz * dz.yz
		
		+ Z.y
	;
	results[23] = a.z - fPlus1 * d_1.z
		
		//g^ij d_ijz
		+ gammaU.xx * dx.xz
		+ gammaU.xy * dx.yz
		+ gammaU.xz * dx.zz
		+ gammaU.xy * dy.xz
		+ gammaU.yy * dy.yz
		+ gammaU.yz * dy.zz
		+ gammaU.xz * dz.xz
		+ gammaU.yz * dz.yz
		+ gammaU.zz * dz.zz

		+ Z.z
	;
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 unused
) {
	<?=addr0?> <?=eqn.cons_t?>* resultU = (<?=addr0?> <?=eqn.cons_t?>*)results;

	real sqrt_f = eig->sqrt_f;
	real f = sqrt_f * sqrt_f;
	real fSq = f * f;
	sym3 gamma = sym3_swap<?=side?>(eig->gamma);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);

	real fMinus1 = f - 1.;
	real lambda_1 = (2. - gui_lambda) / fMinus1;
	real lambda_2 = (2. * f - gui_lambda) / fMinus1;

	real3 a;
	a.x = -.5 * (
			input[0] 
			- input[30] 
			- lambda_2 * input[1] 
			+ lambda_2 * input[29] 
			+ 2. * gammaU.xy * input[7] 
			+ 2. * gammaU.xz * input[8]
		) / gammaU.xx;
	a.y = input[7];
	a.z = input[8];
	
	sym3 d[3];

	<?=addr2?> const sym3 *sym3_1 = (<?=addr2?> const sym3*)(input+1);
	<?=addr2?> const sym3 *sym3_9 = (<?=addr2?> const sym3*)(input+9);
	<?=addr2?> const sym3 *sym3_15 = (<?=addr2?> const sym3*)(input+15);
	<?=addr2?> const sym3 *sym3_23 = (<?=addr2?> const sym3*)(input+23);

	real tr_1 = sym3_dot(gammaU, *sym3_1);
	real tr_9 = sym3_dot(gammaU, *sym3_9);
	real tr_15 = sym3_dot(gammaU, *sym3_15);
	real tr_23 = sym3_dot(gammaU, *sym3_23);

	real _2fPlus1 = 2. * f + 1.;

	d[0] = sym3_scale(*sym3_9, -gammaU.xy);
	d[0] = sym3_add(d[0], sym3_scale(*sym3_15, -gammaU.xz));
	d[0] = sym3_add(d[0], sym3_scale(sym3_add(*sym3_1, *sym3_23), -.5));

	d[0].xx += -.5 * (
		(
			(
				input[0] 
				
				- input[1] 
				+ input[29] 
				
				- input[30] 
				
				- input[1] * lambda_2
				+ input[29] * lambda_2
				
				+ 2. * gammaU.xx * input[21]
				+ 2. * gammaU.xy * input[22] 
				+ 2. * gammaU.xz * input[23] 
						
				+ fMinus1 * (
					tr_23
					- gammaU.xx * input[23]
				)	
					
			) / f
			
			- 6. * gammaU.xy * input[7] 
			- 6. * gammaU.xz * input[8] 
			
			+ 2. * gammaU.xy * tr_9 * _2fPlus1
			+ 2. * gammaU.xz * tr_15 * _2fPlus1
			
			+ 4. * gammaU.xy * input[22] 
			+ 4. * gammaU.xz * input[23]
		
		) / gammaU.xx
	
		- input[1]
		+ input[23]
	);
	
	d[0].xy += -.5 * (
		3. * input[7] 
		- tr_9 * _2fPlus1
		- 2. * input[22] 
	);
	
	d[0].xz += -.5 * (
		3. * input[8] 
		- tr_15 * _2fPlus1
		- 2. * input[23] 
	);
	

	d[0] = sym3_scale(d[0], 1. / gammaU.xx);
	
	d[1] = *(<?=addr2?> const sym3*)(input+9);
	
	d[2] = *(<?=addr2?> const sym3*)(input+15);
	
	sym3 K = sym3_add( *sym3_1, sym3_scale( *sym3_23, .5));
	
	K.xx += .5 * (
			(
				input[0] 
				+ input[30]
			) / sqrt_f
			
			- input[1] * lambda_1
			+ gammaU.xx * input[1]
			+ gammaU.xx * input[23]

			- tr_1
			- tr_23
			
			- input[29] * lambda_1
			
		) / gammaU.xx
		- input[1]
		- .5 * input[23]
	;

	real Theta = .5 * (
		input[1] 
		+ input[29]
	);
	
	real3 Z;
	
	Z.x = .5 * (
		input[1] 
		- input[29] 
		
		- gammaU.xy * input[2] 
		- gammaU.xz * input[3] 
		- gammaU.yy * input[4] 
		- 2. * gammaU.yz * input[5] 
		- gammaU.zz * input[6] 
		
		- gammaU.xy * input[7] 
		- gammaU.xz * input[8] 
		
		- gammaU.xy * gammaU.xx * input[9] 
		- 2. * gammaU.yy * gammaU.xx * input[10] 
		- 2. * gammaU.yz * gammaU.xx * input[11] 
		- gammaU.xy * gammaU.yy * input[12] 
		- 2. * gammaU.xy * gammaU.yz * input[13] 
		- gammaU.xy * gammaU.zz * input[14] 
		
		- gammaU.xz * gammaU.xx * input[15] 
		- 2. * gammaU.yz * gammaU.xx * input[16] 
		- 2. * gammaU.zz * gammaU.xx * input[17] 
		- gammaU.xz * gammaU.yy * input[18] 
		- 2. * gammaU.xz * gammaU.yz * input[19] 
		- gammaU.xz * gammaU.zz * input[20]
		
		+ gammaU.xy * input[24] 
		+ gammaU.xz * input[25] 
		+ gammaU.yy * input[26] 
		+ 2. * gammaU.yz * input[27] 
		+ gammaU.zz * input[28] 
	) / gammaU.xx;
	
	Z.y = .5 * (
		gammaU.xx * input[2]
		+ gammaU.xy * input[4]
		+ gammaU.xz * input[5]
		+ gammaU.xx * input[7]
		
		+ gammaU.xx * gammaU.xx * input[9]
		+ 2. * gammaU.xy * gammaU.xx * input[10]
		+ 2. * gammaU.xz * gammaU.xx * input[11] 
		- gammaU.xx * gammaU.yy * input[12]
		+ 2. * gammaU.xy * gammaU.xy * input[12]
		+ 2. * gammaU.xy * gammaU.xz * input[13]
		+ gammaU.zz * gammaU.xx * input[14] 
		
		- 2. * gammaU.xx * gammaU.yz * input[18]
		+ 2. * gammaU.xz * gammaU.xy * input[18] 
		+ 2. * gammaU.xz * gammaU.xz * input[19] 
		- 2. * gammaU.zz * gammaU.xx * input[19]
		
		- gammaU.xx * input[24]
		- gammaU.xy * input[26]
		- gammaU.xz * input[27] 
	) / gammaU.xx;
	
	Z.z = .5 * (
		gammaU.xx * input[3]
		+ gammaU.xy * input[5] 
		+ gammaU.xz * input[6] 
		+ gammaU.xx * input[8]
		
		+ 2. * gammaU.xy * gammaU.xy * input[13] 
		- 2. * gammaU.xx * gammaU.yy * input[13] 
		- 2. * gammaU.xx * gammaU.yz * input[14] 
		+ 2. * gammaU.xz * gammaU.xy * input[14] 
		
		+ gammaU.xx * gammaU.xx * input[15] 
		+ 2. * gammaU.xy * gammaU.xx * input[16] 
		+ 2. * gammaU.xz * gammaU.xx * input[17] 
		+ gammaU.yy * gammaU.xx * input[18]
		+ 2. * gammaU.xy * gammaU.xz * input[19] 
		- gammaU.xx * gammaU.zz * input[20] 
		+ 2. * gammaU.xz * gammaU.xz * input[20] 
		
		- gammaU.xx * input[25]
		- gammaU.xy * input[27] 
		- gammaU.xz * input[28] 
	) / gammaU.xx;

	<?=addr0?> <?=eqn.cons_t?>* resultsU = (<?=addr0?> <?=eqn.cons_t?>*)results;
	resultsU->alpha = 0;
	resultsU->gamma = _sym3(0,0,0,0,0,0);
	resultsU->a = real3_swap<?=side?>(a);
	resultsU->d[0] = sym3_swap<?=side?>(d[<?=side?>]);
	resultsU->d[1] = sym3_swap<?=side?>(d[<?=side==1 and 0 or 1?>]);
	resultsU->d[2] = sym3_swap<?=side?>(d[<?=side==2 and 0 or 2?>]);
	resultsU->K = sym3_swap<?=side?>(K);
	resultsU->Theta = Theta;
	resultsU->Z = real3_swap<?=side?>(Z);

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
	eigen_leftTransform_<?=side?>___(charvars, &eig, U.ptr, x);
	
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

	real tr_K = sym3_dot(gammaU, U->K);

	//TODO correct source terms

	deriv->alpha += -U->alpha * U->alpha * f * tr_K;
	deriv->gamma = sym3_add(deriv->gamma, sym3_scale(U->K, -2. * U->alpha));
}
