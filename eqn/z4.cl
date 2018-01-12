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
<?=eqn.eigen_t?> eigen_forCell_<?=side?><?=suffix?>(
	<?=addr?> const <?=eqn.cons_t?>* U,
	real3 x 
) {
	<?=eqn.eigen_t?> eig;
	eig.alpha = U->alpha;
	eig.sqrt_f = sqrt(calc_f(U->alpha));
	eig.gamma = U->gamma;
	real det_gamma = sym3_det(U->gamma);
	eig.gammaU = sym3_inv(U->gamma, det_gamma);
	eig.sqrt_gammaUjj = _real3(sqrt(eig.gammaU.xx), sqrt(eig.gammaU.yy), sqrt(eig.gammaU.zz));
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
	sym3 dx = sym3_swap<?=side?>(inputU->d.v<?=side?>);					//3-8
	sym3 dy = sym3_swap<?=side?>(inputU->d.v<?=side==1 and 0 or 1?>);	//9-14
	sym3 dz = sym3_swap<?=side?>(inputU->d.v<?=side==2 and 0 or 2?>);	//15-20
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
	
	_3sym3 d;

	<?=addr2?> const sym3 *sym3_1 = (<?=addr2?> const sym3*)(input+1);
	<?=addr2?> const sym3 *sym3_9 = (<?=addr2?> const sym3*)(input+9);
	<?=addr2?> const sym3 *sym3_15 = (<?=addr2?> const sym3*)(input+15);
	<?=addr2?> const sym3 *sym3_23 = (<?=addr2?> const sym3*)(input+23);

	real tr_1 = sym3_dot(gammaU, *sym3_1);
	real tr_9 = sym3_dot(gammaU, *sym3_9);
	real tr_15 = sym3_dot(gammaU, *sym3_15);
	real tr_23 = sym3_dot(gammaU, *sym3_23);

	real _2fPlus1 = 2. * f + 1.;

	d.x = sym3_scale(*sym3_9, -gammaU.xy);
	d.x = sym3_add(d.x, sym3_scale(*sym3_15, -gammaU.xz));
	d.x = sym3_add(d.x, sym3_scale(sym3_add(*sym3_1, *sym3_23), -.5));

	d.x.xx += -.5 * (
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
	
	d.x.xy += -.5 * (
		3. * input[7] 
		- tr_9 * _2fPlus1
		- 2. * input[22] 
	);
	
	d.x.xz += -.5 * (
		3. * input[8] 
		- tr_15 * _2fPlus1
		- 2. * input[23] 
	);
	

	d.x = sym3_scale(d.x, 1. / gammaU.xx);
	
	d.y = *(<?=addr2?> const sym3*)(input+9);
	
	d.z = *(<?=addr2?> const sym3*)(input+15);
	
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
	resultsU->d.x = sym3_swap<?=side?>(d.v<?=side?>);
	resultsU->d.y = sym3_swap<?=side?>(d.v<?=side==1 and 0 or 1?>);
	resultsU->d.z = sym3_swap<?=side?>(d.v<?=side==2 and 0 or 2?>);
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

	<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>local(&U, x);

	real charvars[numWaves];
	eigen_leftTransform_<?=side?>___(charvars, &eig, U.ptr, x);

	<?=eqn:eigenWaveCodePrefix(side, '&eig', 'x')?>

<? for j=0,eqn.numWaves-1 do 
?>	charvars[<?=j?>] *= <?=eqn:eigenWaveCode(side, '&eig', 'x', j)?>;
<? end
?>
	<?=eqn.cons_t?> F;
	eigen_rightTransform_<?=side?>___(F.ptr, &eig, charvars, x);

	return F;
}
<? end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf)
{
#if 0
	SETBOUNDS_NOGHOST();
	const global <?=eqn.cons_t?>* U = UBuf + index;
	global <?=eqn.cons_t?>* deriv = derivBuf + index;

	const real xi = 1.;	//which is which parameter?  I alwasy forget .. m .. lambda ... xi ... always changing names

	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	real f = calc_f(U->alpha);	//could be based on alpha...

	//K^i_j = gamma^ik K_kj
	mat3 K_ul = (mat3){
<? 
for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<? 	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.<?
		for k,xk in ipairs(xNames) do
		?> + gammaU.<?=sym(i,k)?> * U->K.<?=sym(k,j)?><? 
		end ?>,
<? 	end
?>		},
<? 
end
?>	};

	real tr_K = mat3_trace(K_ul);

	//TODO correct source terms
	//I'm taking this from the 2008 Yano et al "Flux-Vector Splitting..."
	//which itself references (for the source terms) the 2005 Bona et al "Geometrically Motivated ..." 
	//but notice, 2005 Bona paper shows the flux as densitized, so the eigenvalues without any gamma influence

	//alpha,t + ... = alpha beta^k a_k - alpha^2 f (K - 2 Theta)
	deriv->alpha += -U->alpha * U->alpha * f * (tr_K - 2. * U->Theta);

	//a_i,t + ... = b_i^k a_k - b_k^k a_i

	//beta^i_,t + ... = beta^k b_k^i - alpha Q^i
	//Q^i = alpha (a^i - d^i + 2 V^i)
	//V^i = d^i - e^i - Z^i
	
	//gamma_ij,t + ... = 2 beta^k d_kij + b_ji + b_ij - 2 alpha K_ij
	deriv->gamma = sym3_add(deriv->gamma, sym3_scale(U->K, -2. * U->alpha));

	//d_kij,t + ... = b_k^l d_lij - b_l^l d_kij

	//conn_ijk = .5 * (g_ij,k + g_ik,j - g_jk,i) 
	//= d_kij + d_jik - d_ijk
	_3sym3 conn_lll = {
<? 
for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<? 	for jk,xjk in ipairs(symNames) do 
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
?>			.<?=xjk?> = U->d.<?=xk?>.<?=sym(i,j)?> + U->d.<?=xj?>.<?=sym(i,k)?> - U->d.<?=xi?>.<?=sym(j,k)?>,
<? 	end 
?>		},
<? 
end
?>	};

	_3sym3 conn_ull = sym3_3sym3_mul(gammaU, conn_lll);
	_3sym3 d_ull = sym3_3sym3_mul(gammaU, U->d);

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
				?> + gammaU.<?=sym(l,j)?> * conn_ull.<?=xi?>.<?=sym(l,k)?><?
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
?> + U->d.<?=xi?>.<?=sym(j,k)?> * conn_uul.<?=xi?>.<?=sym(j,k)?> * gammaU.<?=sym(k,l)?><?
			end
		end
	end
end
?>;

	//d_i = d_ik^k
	real3 d = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(U->d.<?=xi?>, gammaU),
<? end
?>	};
	
	//e_i = d^k_ki
	real3 e = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0. <?
	for j,xj in ipairs(xNames) do
?> + d_ull.<?=xj?>.<?=sym(j,i)?><?
	end
?>,
<? end
?>	};

	real3 d_u = sym3_real3_mul(gammaU, d);
	real3 e_u = sym3_real3_mul(gammaU, e);
	real3 a_u = sym3_real3_mul(gammaU, U->a);
	real3 Z_u = sym3_real3_mul(gammaU, U->Z);

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
	sym3 KSq_ll = sym3_mat3_to_sym3_mul(U->K, K_ul);

	real tr_KSq = sym3_dot(KSq_ll, gammaU);

	//dsq_ij = d_ikl gamma^lm d^k_mj
	mat3 dsq = {
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
					?> + U->d.<?=xi?>.<?=sym(j,l)?> * gammaU.<?=sym(l,m)?> * d_ull.<?=xk?>.<?=sym(m,j)?>
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
			+ 1/2 (1 - xi) (a_k d^k_ij - 1/2 (a_j (2 e_i - d_i) + a_i (2 e_j - d_j))
				+ 2 (d_ir^m d^r_mj + d_jr^m d^r_mi) - 2 e_k (d_ij^k + d_ji^k)
			)
			+ (d_k + a_k - 2 Z_k) Gamma^k_ij - Gamma^k_mj Gamma^m_ki - (a_i Z_j + a_j Z_i)
			- 2 K^k_i K_kj + (K - 2 Theta) K_ij
		) - 8 pi alpha (S_ij - 1/2 (S - tau) gamma_ij)
	*/
<? for ij,xij in ipairs(symNames) do 
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i], xNames[j]
?>
	deriv->K.<?=xij?> += 
		//TODO shift terms
	U->alpha * (
		.5 * (1. + xi) * (
		<? for k,xk in ipairs(xNames) do ?> 
			-U->a.<?=xk?> * conn_ull.<?=xk?>.<?=xij?>
		<? end ?>
			+ .5 * (U->a.<?=xi?> * d.<?=xj?> + U->a.<?=xj?> * d.<?=xi?>)
		)
		+ .5 * (1. - xi) * (0.
		<? for k,xk in ipairs(xNames) do ?> 
			+ a_u.<?=xk?> * U->d.<?=xk?>.<?=xij?>
		<? end ?>
			- .5 * ( 
				U->a.<?=xj?> * (2. * e.<?=xi?> - d.<?=xi?>) 
				+ U->a.<?=xi?> * (2. * e.<?=xj?> - d.<?=xj?>)
			)
			+ 2. * (dsq.v[<?=i-1?>].s[<?=j-1?>] + dsq.v[<?=j-1?>].s[<?=i-1?>])
		<? for k,xk in ipairs(xNames) do ?> 
			- 2. * e_u.<?=xk?> * (U->d.<?=xi?>.<?=sym(j,k)?> + U->d.<?=xj?>.<?=sym(i,k)?>)
		<? end ?>
		)
		<? for k,xk in ipairs(xNames) do ?> 
		+ (d.<?=xk?> + U->a.<?=xk?> - 2. * U->Z.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		<? end ?>
		- connSq_ll[<?=j-1?>][<?=i-1?>]
		- U->a.<?=xi?> * U->Z.<?=xj?>
		- U->a.<?=xj?> * U->Z.<?=xi?>
		- 2. * KSq_ll.<?=xij?>
		+ U->K.<?=xij?> * (tr_K - 2. * U->Theta)
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
		+ 2. * U->a.<?=xk?> * (d_u.<?=xk?> - e_u.<?=xk?> - 2. * Z_u.<?=xk?>)
		- d_u.<?=xk?> * (d.<?=xk?> - 2. * U->Z.<?=xk?>)
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
			a_i (K - 2 Theta) 
			- a_k K^k_i 
			- K^k_r Gamma^r_ki 
			+ K^k_i (d_k - 2 Z_k)
		)
	*/
	<? for i,xi in ipairs(xNames) do ?> 
	deriv->Z.<?=xi?> += U->alpha * (
		U->a.<?=xi?> * (tr_K - 2. * U->Theta)
		<? for k,xk in ipairs(xNames) do ?> 
		- a_u.<?=xk?> * U->K.<?=sym(i,k)?>
		+ U->K.<?=sym(i,k)?> * (d_u.<?=xk?> - 2. * Z_u.<?=xk?>)
		<? end ?>
		- K_times_conn.<?=xi?>
	);
	<? end ?>

#endif

}
