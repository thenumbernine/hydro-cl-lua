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
	real fSq = f * f;
	real f_3_2 = f * sqrt_f;
	real f_5_2 = fSq * sqrt_f;

	real f_minus_1 = f - 1.;
	real f_minus_1_sq = f_minus_1 * f_minus_1;

	sym3 gamma = sym3_swap<?=side?>(eig->gamma);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);
	
	real tr_K = sym3_dot(gammaU, K);

	real sqrt_gUxx = eig->sqrt_gammaUjj.s<?=side?>;

	real aUx = a.x * gammaU.xx + a.y * gammaU.xy + a.z * gammaU.xz;
	real ZUx =  Z.x * gammaU.xx + Z.y * gammaU.xy + Z.z * gammaU.xz;

	result[0] = -(
		
		+ f_minus_1_sq * aUx / sqrt_gUxx
		
		+ sqrt_f * f_minus_1 * (
			- tr_K * f_minus_1
			+ Theta * (m * f + 2.)
		)

		+ f * (
	
			+ f_minus_1 * (m - 2.) * (
				
				- gammaU.xx * gammaU.yz * dx.yz
				+ gammaU.xy * gammaU.xz * dx.yz


				- gammaU.xx * gammaU.yy * dx.yy
				- gammaU.xx * gammaU.yz * dx.yz
				- gammaU.xx * gammaU.zz * dx.zz
				
				+ gammaU.xy * gammaU.xy * dx.yy
				+ gammaU.xy * gammaU.xz * dx.yz
				+ gammaU.xz * gammaU.xz * dx.zz

				+ gammaU.xx * gammaU.yy * dy.xy
				- gammaU.xy * gammaU.xy * dy.xy
				
				- gammaU.xy * gammaU.xz * dy.xz
				+ gammaU.xx * gammaU.yz * dy.xz
				
				+ gammaU.xx * gammaU.zz * dz.xz
				- gammaU.xz * gammaU.xz * dz.xz

				- gammaU.xz * gammaU.yy * dy.yz
				- gammaU.xy * gammaU.yz * dy.yz
				
				- gammaU.xy * gammaU.zz * dy.zz
				+ gammaU.xz * gammaU.yz * dy.zz
			
				+ gammaU.xx * gammaU.yz * dz.xy
				- gammaU.xy * gammaU.xz * dz.xy
			
				+ gammaU.xy * gammaU.yz * dz.yy
				- gammaU.xz * gammaU.yy * dz.yy
				
				- gammaU.xz * gammaU.yz * dz.yz
				+ gammaU.xy * gammaU.zz * dz.yz
			)
	
			+ f_minus_1 * (m - 2.) * ZUx	
		) / sqrt_gUxx
		
	) / (
		//this is gonna be singular for f=1... which was the harmonic solution, right?
		2. * gammaU.xx * sqrt_f * f_minus_1_sq 
	);
	result[1] = ((((2. * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * K.xz) - (gammaU.xy * gammaU.xz * dy.zz)) + ((gammaU.xy * gammaU.xz * dz.yz) - ((real)pow((real)gammaU.xy, (real)2.) * dy.yz)) + ((real)pow((real)gammaU.xy, (real)2.) * dz.yy) + (2. * gammaU.xy * K.yz * sqrt(gammaU.xx)) + (2. * gammaU.xz * K.zz * sqrt(gammaU.xx)) + ((gammaU.yy * gammaU.xx * dy.yz) - (gammaU.yy * gammaU.xx * dz.yy)) + (((gammaU.yz * gammaU.xx * dy.zz) - (gammaU.yz * gammaU.xx * dz.yz)) - (a.z * gammaU.xx)) + (2. * Z.z * gammaU.xx)) / (4. * gammaU.xx))
	result[2] = (((2. * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * K.xy) + ((gammaU.xy * gammaU.xz * dy.yz) - (gammaU.xy * gammaU.xz * dz.yy)) + (2. * gammaU.xy * K.yy * sqrt(gammaU.xx)) + (((real)pow((real)gammaU.xz, (real)2.) * dy.zz) - ((real)pow((real)gammaU.xz, (real)2.) * dz.yz)) + ((2. * gammaU.xz * K.yz * sqrt(gammaU.xx)) - (gammaU.yz * gammaU.xx * dy.yz)) + ((gammaU.yz * gammaU.xx * dz.yy) - (gammaU.zz * gammaU.xx * dy.zz)) + ((gammaU.zz * gammaU.xx * dz.yz) - (a.y * gammaU.xx)) + (2. * Z.y * gammaU.xx)) / (4. * gammaU.xx))
	result[3] = ((((gamma * gamma.xy * dz.yz) - (gamma * gamma.xz * dz.yy)) + (gammaU.xy * gammaU.yz * dy.yz) + (gammaU.xy * gammaU.zz * dy.zz) + (((((gammaU.xy * a.y) - (2. * gammaU.xy * K.xy * sqrt(gammaU.xx))) - (gammaU.xz * gammaU.yy * dy.yz)) - (gammaU.xz * gammaU.yz * dy.zz)) - (2. * gammaU.xz * K.xz * sqrt(gammaU.xx))) + ((((gammaU.xz * a.z) - (2. * gammaU.yy * K.yy * sqrt(gammaU.xx))) - (4. * gammaU.yz * K.yz * sqrt(gammaU.xx))) - (2. * gammaU.zz * K.zz * sqrt(gammaU.xx))) + (2. * Theta * sqrt(gammaU.xx)) + (2. * Z.x * gammaU.xx)) / (4. * gammaU.xx))
	result[4] = ((((gamma * gamma.xy * dy.zz) - (gamma * gamma.xz * dy.yz)) + ((gamma * gamma.yy * dz.xz) - (gamma * gamma.yy * dx.zz)) + (2. * gamma * gamma.yz * dx.yz) + ((gamma * gamma.zz * dy.xy) - (gamma * gamma.zz * dx.yy)) + (gammaU.xx * gammaU.yz * dy.xz) + (gammaU.xx * gammaU.yz * dz.xy) + (((gammaU.xx * Z.x) - (gammaU.xy * gammaU.xz * dy.xz)) - (gammaU.xy * gammaU.xz * dz.xy)) + (gammaU.xy * gammaU.yz * dz.yy) + (gammaU.xy * gammaU.zz * dz.yz) + (((gammaU.xy * Z.y) - (gammaU.xz * gammaU.yy * dz.yy)) - (gammaU.xz * gammaU.yz * dz.yz)) + (gammaU.xz * Z.z) + (Theta * sqrt(gammaU.xx))) / (2. * sqrt(gammaU.xx)))
	result[5] = (((((gammaU.xx * dz.xz) - (gammaU.xx * dx.zz)) - (gammaU.xy * dy.zz)) + (gammaU.xy * dz.yz) + (K.zz * sqrt(gammaU.xx))) / (2. * sqrt(gammaU.xx)))
	result[6] = (((gammaU.xx * dy.xz) + (((gammaU.xx * dz.xy) - (2. * gammaU.xx * dx.yz)) - (gammaU.xy * dy.yz)) + (gammaU.xy * dz.yy) + ((gammaU.xz * dy.zz) - (gammaU.xz * dz.yz)) + (2. * K.yz * sqrt(gammaU.xx))) / (4. * sqrt(gammaU.xx)))
	result[7] = ((((gammaU.xy * gammaU.xz * dy.zz) - (gammaU.xy * gammaU.xz * dz.yz)) + ((((real)pow((real)gammaU.xy, (real)2.) * dy.yz) - ((real)pow((real)gammaU.xy, (real)2.) * dz.yy)) - (gammaU.yy * gammaU.xx * dy.yz)) + ((gammaU.yy * gammaU.xx * dz.yy) - (gammaU.yz * gammaU.xx * dy.zz)) + (gammaU.yz * gammaU.xx * dz.yz) + (a.z * gammaU.xx)) / (2. * gammaU.xx))
	result[8] = ((-(((gammaU.xy * gammaU.xz * dy.yz) - (gammaU.xy * gammaU.xz * dz.yy)) + ((((real)pow((real)gammaU.xz, (real)2.) * dy.zz) - ((real)pow((real)gammaU.xz, (real)2.) * dz.yz)) - (gammaU.yz * gammaU.xx * dy.yz)) + ((gammaU.yz * gammaU.xx * dz.yy) - (gammaU.zz * gammaU.xx * dy.zz)) + ((gammaU.zz * gammaU.xx * dz.yz) - (a.y * gammaU.xx)))) / (2. * gammaU.xx))
	result[9] = dz.zz
	result[10] = dz.yz
	result[11] = dz.yy
	result[12] = dz.xz
	result[13] = dz.xy
	result[14] = dz.xx
	result[15] = dy.zz
	result[16] = dy.yz
	result[17] = dy.yy
	result[18] = dy.xz
	result[19] = dy.xy
	result[20] = dy.xx
	result[21] = ((-((gammaU.xy * dy.xz) + ((gammaU.xy * dz.xy) - (2. * gammaU.xy * dx.yz)) + (((2. * gammaU.xz * dz.xz) - (2. * gammaU.xz * dx.zz)) - (gammaU.yy * dy.yz)) + ((gammaU.yy * dz.yy) - (gammaU.yz * dy.zz)) + (gammaU.yz * dz.yz) + ((a.z - (2. * Z.z)) - (2. * dx.xz * gammaU.xx)))) / (2. * gammaU.xx))
	result[22] = ((-(((2. * gammaU.xy * dy.xy) - (2. * gammaU.xy * dx.yy)) + (gammaU.xz * dy.xz) + ((gammaU.xz * dz.xy) - (2. * gammaU.xz * dx.yz)) + ((gammaU.yz * dy.yz) - (gammaU.yz * dz.yy)) + ((gammaU.zz * dy.zz) - (gammaU.zz * dz.yz)) + ((a.y - (2. * Z.y)) - (2. * dx.xy * gammaU.xx)))) / (2. * gammaU.xx))
	result[23] = ((-(((2. * gamma * gamma.xy * dy.zz * f) - (2. * gamma * gamma.xy * dz.yz * f)) + ((gamma * gamma.xy * m * dz.yz * f) - (2. * gamma * gamma.xz * dy.yz * f)) + ((2. * gamma * gamma.xz * dz.yy * f) - (gamma * gamma.xz * m * dz.yy * f)) + (((gammaU.xx * gammaU.yy * dy.xy * f) - (gammaU.xx * gammaU.yy * dx.yy * f)) - (gammaU.xx * gammaU.yy * m * dy.xy * f)) + (gammaU.xx * gammaU.yy * m * dx.yy * f) + (gammaU.xx * gammaU.yz * dy.xz * f) + ((((gammaU.xx * gammaU.yz * dz.xy * f) - (2. * gammaU.xx * gammaU.yz * dx.yz * f)) - (gammaU.xx * gammaU.yz * m * dy.xz * f)) - (gammaU.xx * gammaU.yz * m * dz.xy * f)) + (2. * gammaU.xx * gammaU.yz * m * dx.yz * f) + (((gammaU.xx * gammaU.zz * dz.xz * f) - (gammaU.xx * gammaU.zz * dx.zz * f)) - (gammaU.xx * gammaU.zz * m * dz.xz * f)) + (((gammaU.xx * gammaU.zz * m * dx.zz * f) - (2. * gammaU.xy * gammaU.xz * dy.xz * f)) - (2. * gammaU.xy * gammaU.xz * dz.xy * f)) + (4. * gammaU.xy * gammaU.xz * dx.yz * f) + (gammaU.xy * gammaU.xz * m * dy.xz * f) + (((gammaU.xy * gammaU.xz * m * dz.xy * f) - (2. * gammaU.xy * gammaU.xz * m * dx.yz * f)) - (2. * (real)pow((real)gammaU.xy, (real)2.) * dy.xy * f)) + (gammaU.xy * a.y) + (2. * gammaU.xy * Z.y * f) + (2. * (real)pow((real)gammaU.xy, (real)2.) * dx.yy * f) + (((((real)pow((real)gammaU.xy, (real)2.) * m * dy.xy * f) - (gammaU.xy * m * Z.y * f)) - ((real)pow((real)gammaU.xy, (real)2.) * m * dx.yy * f)) - (2. * (real)pow((real)gammaU.xz, (real)2.) * dz.xz * f)) + (gammaU.xz * a.z) + (2. * gammaU.xz * Z.z * f) + (2. * (real)pow((real)gammaU.xz, (real)2.) * dx.zz * f) + ((((((real)pow((real)gammaU.xz, (real)2.) * m * dz.xz * f) - (gammaU.xz * m * Z.z * f)) - ((real)pow((real)gammaU.xz, (real)2.) * m * dx.zz * f)) - (f * gammaU.xy * a.y)) - (f * gammaU.xz * a.z)) + (((a.x * gammaU.xx) - (dx.xx * (real)pow((real)gammaU.xx, (real)2.) * f)) - (m * gamma * gamma.xy * dy.zz * f)) + ((m * gamma * gamma.xz * dy.yz * f) - (m * Z.x * gammaU.xx * f)))) / ((real)pow((real)gammaU.xx, (real)2.) * f))
	result[24] = ((-((2. * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * K.xz) + ((gammaU.xy * gammaU.xz * dy.zz) - (gammaU.xy * gammaU.xz * dz.yz)) + (((real)pow((real)gammaU.xy, (real)2.) * dy.yz) - ((real)pow((real)gammaU.xy, (real)2.) * dz.yy)) + (2. * gammaU.xy * K.yz * sqrt(gammaU.xx)) + ((2. * gammaU.xz * K.zz * sqrt(gammaU.xx)) - (gammaU.yy * gammaU.xx * dy.yz)) + ((gammaU.yy * gammaU.xx * dz.yy) - (gammaU.yz * gammaU.xx * dy.zz)) + (gammaU.yz * gammaU.xx * dz.yz) + ((a.z * gammaU.xx) - (2. * Z.z * gammaU.xx)))) / (4. * gammaU.xx))
	result[25] = ((-(((2. * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * K.xy) - (gammaU.xy * gammaU.xz * dy.yz)) + (gammaU.xy * gammaU.xz * dz.yy) + ((2. * gammaU.xy * K.yy * sqrt(gammaU.xx)) - ((real)pow((real)gammaU.xz, (real)2.) * dy.zz)) + ((real)pow((real)gammaU.xz, (real)2.) * dz.yz) + (2. * gammaU.xz * K.yz * sqrt(gammaU.xx)) + ((gammaU.yz * gammaU.xx * dy.yz) - (gammaU.yz * gammaU.xx * dz.yy)) + ((gammaU.zz * gammaU.xx * dy.zz) - (gammaU.zz * gammaU.xx * dz.yz)) + ((a.y * gammaU.xx) - (2. * Z.y * gammaU.xx)))) / (4. * gammaU.xx))
	result[26] = ((((gamma * gamma.xy * dz.yz) - (gamma * gamma.xz * dz.yy)) + (gammaU.xy * gammaU.yz * dy.yz) + (gammaU.xy * gammaU.zz * dy.zz) + (gammaU.xy * a.y) + (((2. * gammaU.xy * K.xy * sqrt(gammaU.xx)) - (gammaU.xz * gammaU.yy * dy.yz)) - (gammaU.xz * gammaU.yz * dy.zz)) + (2. * gammaU.xz * K.xz * sqrt(gammaU.xx)) + (gammaU.xz * a.z) + (2. * gammaU.yy * K.yy * sqrt(gammaU.xx)) + (4. * gammaU.yz * K.yz * sqrt(gammaU.xx)) + ((2. * gammaU.zz * K.zz * sqrt(gammaU.xx)) - (2. * Theta * sqrt(gammaU.xx))) + (2. * Z.x * gammaU.xx)) / (4. * gammaU.xx))
	result[27] = (((((gamma * gamma.xy * dz.yz) - (gamma * gamma.xz * dz.yy)) - (gamma * gamma.yy * dz.xz)) + (gamma * gamma.yy * dx.zz) + (gamma * gamma.yz * dy.xz) + ((gamma * gamma.yz * dz.xy) - (gamma * gamma.zz * dy.xy)) + (gamma * gamma.zz * dx.yy) + (((2. * gammaU.xx * gammaU.yz * dx.yz) - (gammaU.xx * Z.x)) - (2. * gammaU.xy * gammaU.xz * dx.yz)) + (gammaU.xy * gammaU.yz * dy.yz) + (((((gammaU.xy * gammaU.zz * dy.zz) - (gammaU.xy * Z.y)) - (gammaU.xz * gammaU.yy * dy.yz)) - (gammaU.xz * gammaU.yz * dy.zz)) - (gammaU.xz * Z.z)) + (Theta * sqrt(gammaU.xx))) / (2. * sqrt(gammaU.xx)))
	result[28] = ((-((((gammaU.xx * dz.xz) - (gammaU.xx * dx.zz)) - (gammaU.xy * dy.zz)) + ((gammaU.xy * dz.yz) - (K.zz * sqrt(gammaU.xx))))) / (2. * sqrt(gammaU.xx)))
	result[29] = ((-((gammaU.xx * dy.xz) + (((gammaU.xx * dz.xy) - (2. * gammaU.xx * dx.yz)) - (gammaU.xy * dy.yz)) + (gammaU.xy * dz.yy) + (((gammaU.xz * dy.zz) - (gammaU.xz * dz.yz)) - (2. * K.yz * sqrt(gammaU.xx))))) / (4. * sqrt(gammaU.xx)))
	result[30] = (((2. * gammaU.xy * a.y * (1. / sqrt(gammaU.xx))) +
		 ((2. * gammaU.xy * a.y * fSq * (1. / sqrt(gammaU.xx))) -
		 (4. * gammaU.xy * a.y * f * (1. / sqrt(gammaU.xx)))) +
		 (4. * gammaU.xy * K.xy * f_5_2) +
		 ((4. * gammaU.xy * K.xy * sqrt(f)) -
		 (8. * gammaU.xy * K.xy * f_3_2)) +
		 (4. * gammaU.xz * K.xz * f_5_2) +
		 ((4. * gammaU.xz * K.xz * sqrt(f)) -
		 (8. * gammaU.xz * K.xz * f_3_2)) +
		 (2. * gammaU.xz * a.z * (1. / sqrt(gammaU.xx))) +
		 ((2. * gammaU.xz * a.z * fSq * (1. / sqrt(gammaU.xx))) -
		 (4. * gammaU.xz * a.z * f * (1. / sqrt(gammaU.xx)))) +
		 ((2. * gammaU.yy * K.yy * f_5_2) -
		 (4. * gammaU.yy * K.yy * f_3_2)) +
		 (2. * gammaU.yy * K.yy * sqrt(f)) +
		 (4. * gammaU.yz * K.yz * sqrt(f)) +
		 ((4. * gammaU.yz * K.yz * f_5_2) -
		 (8. * gammaU.yz * K.yz * f_3_2)) +
		 (2. * gammaU.zz * K.zz * sqrt(f)) +
		 ((2. * gammaU.zz * K.zz * f_5_2) -
		 (4. * gammaU.zz * K.zz * f_3_2)) +
		 ((4. * f * gamma * gamma.xy * dy.zz * (1. / sqrt(gammaU.xx))) -
		 (4. * fSq * gamma * gamma.xy * dy.zz * (1. / sqrt(gammaU.xx)))) +
		 (((2. * fSq * gamma * gamma.xy * m * dy.zz * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gamma * gamma.xy * m * dy.zz * (1. / sqrt(gammaU.xx)))) -
		 (4. * f * gamma * gamma.xz * dy.yz * (1. / sqrt(gammaU.xx)))) +
		 ((4. * fSq * gamma * gamma.xz * dy.yz * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gamma * gamma.xz * m * dy.yz * (1. / sqrt(gammaU.xx)))) +
		 (2. * f * gamma * gamma.xz * m * dy.yz * (1. / sqrt(gammaU.xx))) +
		 ((4. * f * gamma * gamma.yy * dz.xz * (1. / sqrt(gammaU.xx))) -
		 (4. * fSq * gamma * gamma.yy * dz.xz * (1. / sqrt(gammaU.xx)))) +
		 ((4. * fSq * gamma * gamma.yy * dx.zz * (1. / sqrt(gammaU.xx))) -
		 (4. * f * gamma * gamma.yy * dx.zz * (1. / sqrt(gammaU.xx)))) +
		 ((2. * fSq * gamma * gamma.yy * m * dz.xz * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gamma * gamma.yy * m * dz.xz * (1. / sqrt(gammaU.xx)))) +
		 ((2. * f * gamma * gamma.yy * m * dx.zz * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gamma * gamma.yy * m * dx.zz * (1. / sqrt(gammaU.xx)))) +
		 (((8. * f * gamma * gamma.yz * dx.yz * (1. / sqrt(gammaU.xx))) -
		 (8. * fSq * gamma * gamma.yz * dx.yz * (1. / sqrt(gammaU.xx)))) -
		 (4. * f * gamma * gamma.yz * m * dx.yz * (1. / sqrt(gammaU.xx)))) +
		 (4. * fSq * gamma * gamma.yz * m * dx.yz * (1. / sqrt(gammaU.xx))) +
		 ((4. * f * gamma * gamma.zz * dy.xy * (1. / sqrt(gammaU.xx))) -
		 (4. * fSq * gamma * gamma.zz * dy.xy * (1. / sqrt(gammaU.xx)))) +
		 ((4. * fSq * gamma * gamma.zz * dx.yy * (1. / sqrt(gammaU.xx))) -
		 (4. * f * gamma * gamma.zz * dx.yy * (1. / sqrt(gammaU.xx)))) +
		 ((2. * fSq * gamma * gamma.zz * m * dy.xy * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gamma * gamma.zz * m * dy.xy * (1. / sqrt(gammaU.xx)))) +
		 ((2. * f * gamma * gamma.zz * m * dx.yy * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gamma * gamma.zz * m * dx.yy * (1. / sqrt(gammaU.xx)))) +
		 ((4. * f * sqrt(gammaU.xx) * gammaU.yz * dy.xz) -
		 (4. * fSq * sqrt(gammaU.xx) * gammaU.yz * dy.xz)) +
		 (((4. * f * sqrt(gammaU.xx) * gammaU.yz * dz.xy) -
		 (4. * fSq * sqrt(gammaU.xx) * gammaU.yz * dz.xy)) -
		 (2. * f * sqrt(gammaU.xx) * gammaU.yz * m * dy.xz)) +
		 (2. * fSq * sqrt(gammaU.xx) * gammaU.yz * m * dy.xz) +
		 ((2. * fSq * sqrt(gammaU.xx) * gammaU.yz * m * dz.xy) -
		 (2. * f * sqrt(gammaU.xx) * gammaU.yz * m * dz.xy)) +
		 (((4. * fSq * gammaU.xy * gammaU.xz * dy.xz * (1. / sqrt(gammaU.xx))) -
		 (4. * f * gammaU.xy * gammaU.xz * dy.xz * (1. / sqrt(gammaU.xx)))) -
		 (4. * f * gammaU.xy * gammaU.xz * dz.xy * (1. / sqrt(gammaU.xx)))) +
		 (4. * fSq * gammaU.xy * gammaU.xz * dz.xy * (1. / sqrt(gammaU.xx))) +
		 ((2. * f * gammaU.xy * gammaU.xz * m * dy.xz * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gammaU.xy * gammaU.xz * m * dy.xz * (1. / sqrt(gammaU.xx)))) +
		 (((2. * f * gammaU.xy * gammaU.xz * m * dz.xy * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gammaU.xy * gammaU.xz * m * dz.xy * (1. / sqrt(gammaU.xx)))) -
		 (4. * fSq * gammaU.xy * gammaU.yz * dz.yy * (1. / sqrt(gammaU.xx)))) +
		 ((4. * f * gammaU.xy * gammaU.yz * dz.yy * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gammaU.xy * gammaU.yz * m * dz.yy * (1. / sqrt(gammaU.xx)))) +
		 ((2. * fSq * gammaU.xy * gammaU.yz * m * dz.yy * (1. / sqrt(gammaU.xx))) -
		 (4. * fSq * gammaU.xy * gammaU.zz * dz.yz * (1. / sqrt(gammaU.xx)))) +
		 (4. * f * gammaU.xy * gammaU.zz * dz.yz * (1. / sqrt(gammaU.xx))) +
		 ((2. * fSq * gammaU.xy * gammaU.zz * m * dz.yz * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gammaU.xy * gammaU.zz * m * dz.yz * (1. / sqrt(gammaU.xx)))) +
		 ((4. * f * gammaU.xy * Z.y * (1. / sqrt(gammaU.xx))) -
		 (4. * fSq * gammaU.xy * Z.y * (1. / sqrt(gammaU.xx)))) +
		 (((2. * fSq * gammaU.xy * m * Z.y * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gammaU.xy * m * Z.y * (1. / sqrt(gammaU.xx)))) -
		 (4. * f * gammaU.xz * gammaU.yy * dz.yy * (1. / sqrt(gammaU.xx)))) +
		 (4. * fSq * gammaU.xz * gammaU.yy * dz.yy * (1. / sqrt(gammaU.xx))) +
		 ((2. * f * gammaU.xz * gammaU.yy * m * dz.yy * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gammaU.xz * gammaU.yy * m * dz.yy * (1. / sqrt(gammaU.xx)))) +
		 ((4. * fSq * gammaU.xz * gammaU.yz * dz.yz * (1. / sqrt(gammaU.xx))) -
		 (4. * f * gammaU.xz * gammaU.yz * dz.yz * (1. / sqrt(gammaU.xx)))) +
		 ((2. * f * gammaU.xz * gammaU.yz * m * dz.yz * (1. / sqrt(gammaU.xx))) -
		 (2. * fSq * gammaU.xz * gammaU.yz * m * dz.yz * (1. / sqrt(gammaU.xx)))) +
		 ((4. * f * gammaU.xz * Z.z * (1. / sqrt(gammaU.xx))) -
		 (4. * fSq * gammaU.xz * Z.z * (1. / sqrt(gammaU.xx)))) +
		 (((2. * fSq * gammaU.xz * m * Z.z * (1. / sqrt(gammaU.xx))) -
		 (2. * f * gammaU.xz * m * Z.z * (1. / sqrt(gammaU.xx)))) -
		 (4. * fSq * Z.x * sqrt(gammaU.xx))) +
		 (4. * f * Z.x * sqrt(gammaU.xx)) +
		 ((2. * f_3_2 * m * Theta) -
		 (2. * f_5_2 * m * Theta)) +
		 ((2. * fSq * m * Z.x * sqrt(gammaU.xx)) -
		 (2. * f * m * Z.x * sqrt(gammaU.xx))) +
		 (2. * a.x * sqrt(gammaU.xx)) +
		 (((2. * a.x * sqrt(gammaU.xx) * fSq) -
		 (4. * a.x * sqrt(gammaU.xx) * f)) -
		 (4. * K.xx * gammaU.xx * f_3_2)) +
		 (2. * K.xx * gammaU.xx * sqrt(f)) +
		 ((2. * K.xx * f_5_2 * gammaU.xx) -
		 (4. * Theta * sqrt(f))) +
		 (4. * Theta * f_3_2)) / (gammaU.xx * sqrt(f) * ((4. -
		 (8. * f)) +
		 (4. * fSq))))


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
	real f_3_2 = f * sqrt_f;
	sym3 gamma = sym3_swap<?=side?>(eig->gamma);
	sym3 gammaU = sym3_swap<?=side?>(eig->gammaU);

	result[0] = (((((gamma * gamma.xy * input[10] * (1. / sqrt(gammaU.xx))) - (gamma * gamma.xy * input[10] * f * (1. / sqrt(gammaU.xx)))) - (gamma * gamma.xz * input[11] * (1. / sqrt(gammaU.xx)))) + (gamma * gamma.xz * input[11] * f * (1. / sqrt(gammaU.xx))) + ((gammaU.xy * gammaU.yz * input[16] * (1. / sqrt(gammaU.xx))) - (gammaU.xy * gammaU.yz * input[16] * f * (1. / sqrt(gammaU.xx)))) + (((gammaU.xy * gammaU.zz * input[15] * (1. / sqrt(gammaU.xx))) - (gammaU.xy * gammaU.zz * input[15] * f * (1. / sqrt(gammaU.xx)))) - (2. * gammaU.xy * input[8] * (1. / sqrt(gammaU.xx)))) + ((2. * gammaU.xy * input[8] * f * (1. / sqrt(gammaU.xx))) - (gammaU.xz * gammaU.yy * input[16] * (1. / sqrt(gammaU.xx)))) + ((gammaU.xz * gammaU.yy * input[16] * f * (1. / sqrt(gammaU.xx))) - (gammaU.xz * gammaU.yz * input[15] * (1. / sqrt(gammaU.xx)))) + ((gammaU.xz * gammaU.yz * input[15] * f * (1. / sqrt(gammaU.xx))) - (2. * gammaU.xz * input[7] * (1. / sqrt(gammaU.xx)))) + ((2. * gammaU.xz * input[7] * f * (1. / sqrt(gammaU.xx))) - (sqrt(f) * gammaU.xx * input[0])) + (f_3_2 * gammaU.xx * input[0]) + ((sqrt(f) * gammaU.xx * input[30]) - (f_3_2 * gammaU.xx * input[30])) + (((2. * f * input[27]) - (2. * f * input[4])) - (f * m * input[27])) + (f * m * input[4])) / (sqrt(gammaU.xx) * (1. - f)))
	result[1] = ((-(((gammaU.xy * gammaU.xz * input[11]) - (gammaU.xy * gammaU.xz * input[16])) + ((((real)pow((real)gammaU.xz, (real)2.) * input[10]) - ((real)pow((real)gammaU.xz, (real)2.) * input[15])) - (gammaU.yz * gammaU.xx * input[11])) + ((gammaU.yz * gammaU.xx * input[16]) - (gammaU.zz * gammaU.xx * input[10])) + ((gammaU.zz * gammaU.xx * input[15]) - (2. * input[8] * gammaU.xx)))) / gammaU.xx)
	result[2] = ((((gammaU.xy * gammaU.xz * input[10]) - (gammaU.xy * gammaU.xz * input[15])) + ((((real)pow((real)gammaU.xy, (real)2.) * input[11]) - ((real)pow((real)gammaU.xy, (real)2.) * input[16])) - (gammaU.yy * gammaU.xx * input[11])) + ((gammaU.yy * gammaU.xx * input[16]) - (gammaU.yz * gammaU.xx * input[10])) + (gammaU.yz * gammaU.xx * input[15]) + (2. * input[7] * gammaU.xx)) / gammaU.xx)
	result[3] = (((((gamma * gamma.zz * input[26] * gammaU.xx * gammaU.yy * sqrt(f)) - (gamma * gamma.zz * input[26] * gammaU.xx * f_3_2 * gammaU.yy)) - (gamma * gamma.zz * input[26] * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f))) + (gamma * gamma.zz * input[26] * f_3_2 * (real)pow((real)gammaU.xy, (real)2.)) + ((gamma * gamma.zz * input[27] * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f)) - (gamma * gamma.zz * input[27] * sqrt(gammaU.xx) * gammaU.yy * sqrt(f))) + ((gamma * gamma.zz * input[3] * gammaU.xx * gammaU.yy * sqrt(f)) - (gamma * gamma.zz * input[3] * gammaU.xx * gammaU.yy * f_3_2)) + ((gamma * gamma.zz * input[3] * (real)pow((real)gammaU.xy, (real)2.) * f_3_2) - (gamma * gamma.zz * input[3] * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f))) + ((((real)pow((real)gammaU.xx, (real)(3. / 2.)) * (real)pow((real)gammaU.yy, (real)2.) * f_3_2 * input[4]) - ((real)pow((real)gammaU.xx, (real)(3. / 2.)) * (real)pow((real)gammaU.yy, (real)2.) * f_3_2 * m * input[4])) - (sqrt(gammaU.xx) * gammaU.yy * input[4] * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f))) + ((real)pow((real)gammaU.xx, (real)(3. / 2.)) * (real)pow((real)gammaU.yy, (real)2.) * input[4] * sqrt(f)) + (((2. * gammaU.xy * gammaU.xz * (real)pow((real)gammaU.yy, (real)2.) * input[6] * sqrt(gammaU.xx) * sqrt(f)) - (2. * (real)pow((real)gammaU.xy, (real)3.) * gammaU.xz * gammaU.yy * input[6] * (1. / sqrt(gammaU.xx)) * sqrt(f))) - (2. * gammaU.xy * gammaU.xz * (real)pow((real)gammaU.yy, (real)2.) * input[6] * f_3_2 * sqrt(gammaU.xx))) + ((2. * (real)pow((real)gammaU.xy, (real)3.) * gammaU.xz * gammaU.yy * input[6] * f_3_2 * (1. / sqrt(gammaU.xx))) - (gammaU.xy * gammaU.yy * input[25] * gamma * gamma.zz * sqrt(f))) + (gammaU.xy * gammaU.yy * input[25] * f_3_2 * gamma * gamma.zz) + (((gammaU.xy * gammaU.yy * input[2] * f_3_2 * gamma * gamma.zz) - (gammaU.xy * gammaU.yy * input[2] * sqrt(f) * gamma * gamma.zz)) - (2. * (real)pow((real)gammaU.xy, (real)2.) * gammaU.yz * input[6] * sqrt(gammaU.xx) * gammaU.yy * sqrt(f))) + ((2. * (real)pow((real)gammaU.xy, (real)4.) * gammaU.yz * input[6] * (1. / sqrt(gammaU.xx)) * sqrt(f)) - (2. * (real)pow((real)gammaU.xy, (real)4.) * gammaU.yz * input[6] * f_3_2 * (1. / sqrt(gammaU.xx)))) + (2. * (real)pow((real)gammaU.xy, (real)2.) * gammaU.yz * input[6] * f_3_2 * sqrt(gammaU.xx) * gammaU.yy) + (((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[28] * sqrt(gammaU.xx) * gammaU.yy * sqrt(f)) - ((real)pow((real)gammaU.xy, (real)4.) * gammaU.zz * input[28] * (1. / sqrt(gammaU.xx)) * sqrt(f))) + ((((real)pow((real)gammaU.xy, (real)4.) * gammaU.zz * input[28] * f_3_2 * (1. / sqrt(gammaU.xx))) - ((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[28] * f_3_2 * sqrt(gammaU.xx) * gammaU.yy)) - ((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[5] * sqrt(gammaU.xx) * gammaU.yy * sqrt(f))) + (((real)pow((real)gammaU.xy, (real)4.) * gammaU.zz * input[5] * (1. / sqrt(gammaU.xx)) * sqrt(f)) - ((real)pow((real)gammaU.xy, (real)4.) * gammaU.zz * input[5] * f_3_2 * (1. / sqrt(gammaU.xx)))) + ((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[5] * f_3_2 * sqrt(gammaU.xx) * gammaU.yy) + ((2. * gammaU.xy * gamma.xz * input[29] * gamma * sqrt(gammaU.xx) * gammaU.yy * sqrt(f)) - (2. * (real)pow((real)gammaU.xy, (real)3.) * gamma.xz * input[29] * gamma * (1. / sqrt(gammaU.xx)) * sqrt(f))) + (((2. * (real)pow((real)gammaU.xy, (real)3.) * gamma.xz * input[29] * f_3_2 * gamma * (1. / sqrt(gammaU.xx))) - (2. * gammaU.xy * gamma.xz * input[29] * f_3_2 * gamma * sqrt(gammaU.xx) * gammaU.yy)) - ((real)pow((real)gammaU.xy, (real)4.) * f_3_2 * input[27] * (1. / sqrt(gammaU.xx)))) + (((real)pow((real)gammaU.xy, (real)2.) * f_3_2 * input[27] * sqrt(gammaU.xx) * gammaU.yy) - ((real)pow((real)gammaU.xy, (real)2.) * input[26] * gammaU.xx * gammaU.yy * sqrt(f))) + ((real)pow((real)gammaU.xy, (real)2.) * input[26] * gammaU.xx * f_3_2 * gammaU.yy) + ((((real)pow((real)gammaU.xy, (real)4.) * input[26] * sqrt(f)) - ((real)pow((real)gammaU.xy, (real)4.) * input[26] * f_3_2)) - ((real)pow((real)gammaU.xy, (real)2.) * input[27] * sqrt(gammaU.xx) * gammaU.yy * sqrt(f))) + ((real)pow((real)gammaU.xy, (real)4.) * input[27] * (1. / sqrt(gammaU.xx)) * sqrt(f)) + (((real)pow((real)gammaU.xy, (real)2.) * input[3] * gammaU.xx * gammaU.yy * f_3_2) - ((real)pow((real)gammaU.xy, (real)2.) * input[3] * gammaU.xx * gammaU.yy * sqrt(f))) + (((real)pow((real)gammaU.xy, (real)4.) * input[3] * sqrt(f)) - ((real)pow((real)gammaU.xy, (real)4.) * input[3] * f_3_2)) + (((gammaU.xz * gammaU.yy * input[1] * f_3_2 * gamma * gamma.zz) - (gammaU.xz * gammaU.yy * input[1] * sqrt(f) * gamma * gamma.zz)) - (gammaU.xz * gammaU.yy * input[24] * gamma * gamma.zz * sqrt(f))) + (gammaU.xz * gammaU.yy * input[24] * f_3_2 * gamma * gamma.zz) + (((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[28] * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f)) - ((real)pow((real)gammaU.xz, (real)2.) * (real)pow((real)gammaU.yy, (real)2.) * input[28] * sqrt(gammaU.xx) * sqrt(f))) + ((((real)pow((real)gammaU.xz, (real)2.) * (real)pow((real)gammaU.yy, (real)2.) * input[28] * f_3_2 * sqrt(gammaU.xx)) - ((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[28] * f_3_2 * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.))) - ((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[5] * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f))) + (((real)pow((real)gammaU.xz, (real)2.) * (real)pow((real)gammaU.yy, (real)2.) * input[5] * sqrt(gammaU.xx) * sqrt(f)) - ((real)pow((real)gammaU.xz, (real)2.) * (real)pow((real)gammaU.yy, (real)2.) * input[5] * f_3_2 * sqrt(gammaU.xx))) + ((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[5] * f_3_2 * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.)) + ((f_3_2 * gamma * gamma.zz * input[27] * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.)) - (f_3_2 * gamma * gamma.zz * input[27] * sqrt(gammaU.xx) * gammaU.yy)) + (((2. * f_3_2 * (real)pow((real)gammaU.xy, (real)4.) * input[4] * (1. / sqrt(gammaU.xx))) - (3. * f_3_2 * (real)pow((real)gammaU.xy, (real)2.) * input[4] * sqrt(gammaU.xx) * gammaU.yy)) - (f_3_2 * (real)pow((real)gammaU.xy, (real)4.) * m * input[4] * (1. / sqrt(gammaU.xx)))) + ((2. * f_3_2 * (real)pow((real)gammaU.xy, (real)2.) * m * input[4] * sqrt(gammaU.xx) * gammaU.yy) - (f_3_2 * m * gamma * gamma.zz * input[27] * (1. / sqrt(gammaU.xx)) * (real)pow((real)gammaU.xy, (real)2.))) + ((f_3_2 * m * gamma * gamma.zz * input[27] * sqrt(gammaU.xx) * gammaU.yy) - (input[0] * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.) * gamma * gamma.zz)) + (input[0] * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.) * gamma * gamma.zz * f) + ((input[0] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy * gamma * gamma.zz) - (input[0] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy * gamma * gamma.zz * f)) + (((input[23] * gammaU.xx * gamma * gamma.zz * (real)pow((real)gammaU.xy, (real)2.) * sqrt(f)) - (input[23] * (real)pow((real)gammaU.xx, (real)2.) * gamma * gamma.zz * gammaU.yy * sqrt(f))) - (input[23] * gammaU.xx * f_3_2 * gamma * gamma.zz * (real)pow((real)gammaU.xy, (real)2.))) + (input[23] * (real)pow((real)gammaU.xx, (real)2.) * f_3_2 * gamma * gamma.zz * gammaU.yy) + (((input[30] * sqrt(gammaU.xx) * gamma * gamma.zz * (real)pow((real)gammaU.xy, (real)2.)) - (input[30] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gamma * gamma.zz * gammaU.yy)) - (input[30] * sqrt(gammaU.xx) * f * gamma * gamma.zz * (real)pow((real)gammaU.xy, (real)2.))) + (input[30] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * f * gamma * gamma.zz * gammaU.yy)) / (gammaU.xx * gamma * gamma.zz * sqrt(f) * ((((real)pow((real)gammaU.xy, (real)2.) - (f * (real)pow((real)gammaU.xy, (real)2.))) - (gammaU.xx * gammaU.yy)) + (gammaU.xx * f * gammaU.yy))))
	result[4] = ((((((2. * sqrt(gammaU.xx) * gammaU.xy * gammaU.yz * input[29]) - (2. * sqrt(gammaU.xx) * gammaU.xy * gammaU.yz * input[6])) - (gammaU.xy * gammaU.xz * input[1])) - (gammaU.xy * gammaU.xz * input[24])) + (((((gammaU.xy * gamma.yy * input[28] * (1. / sqrt(gammaU.xx)) * gamma) - (gammaU.xy * gamma.yy * input[5] * (1. / sqrt(gammaU.xx)) * gamma)) - (gammaU.xy * input[26] * gammaU.xx)) - (gammaU.xy * input[27] * sqrt(gammaU.xx))) - (gammaU.xy * input[3] * gammaU.xx)) + ((gammaU.xy * input[4] * sqrt(gammaU.xx)) - (gammaU.xz * sqrt(gammaU.xx) * gammaU.yy * input[29])) + ((gammaU.xz * sqrt(gammaU.xx) * gammaU.yy * input[6]) - (gammaU.xz * (real)pow((real)gammaU.xy, (real)2.) * input[29] * (1. / sqrt(gammaU.xx)))) + (((gammaU.xz * (real)pow((real)gammaU.xy, (real)2.) * input[6] * (1. / sqrt(gammaU.xx))) - (gammaU.yy * input[25] * gammaU.xx)) - (gammaU.yy * input[2] * gammaU.xx)) + (input[22] * gammaU.xx * gamma * gamma.zz)) / (gammaU.xx * gamma * gamma.zz))
	result[5] = ((-(((gammaU.xy * input[29] * (1. / sqrt(gammaU.xx))) - (gammaU.xy * input[6] * (1. / sqrt(gammaU.xx)))) + ((gammaU.xz * input[28] * (1. / sqrt(gammaU.xx))) - (gammaU.xz * input[5] * (1. / sqrt(gammaU.xx)))) + (input[1] - (input[21] * gammaU.xx)) + input[24])) / gammaU.xx)
	result[6] = ((((real)pow((real)gammaU.xx, (real)(3. / 2.)) * input[26]) + (gammaU.xx * input[27]) + (((real)pow((real)gammaU.xx, (real)(3. / 2.)) * input[3]) - (gammaU.xx * input[4])) + (gammaU.xy * input[25] * sqrt(gammaU.xx)) + ((gammaU.xy * input[2] * sqrt(gammaU.xx)) - (gammaU.xz * input[11] * (1. / sqrt(gammaU.xx)) * gamma * gamma.zz)) + (gammaU.xz * input[16] * (1. / sqrt(gammaU.xx)) * gamma * gamma.zz) + (gammaU.xz * input[1] * sqrt(gammaU.xx)) + ((gammaU.xz * input[24] * sqrt(gammaU.xx)) - (gamma.yy * input[28] * gamma)) + (gamma.yy * input[5] * gamma) + ((2. * gamma.yz * input[29] * gamma) - (2. * gamma.yz * input[6] * gamma)) + (input[19] * sqrt(gammaU.xx) * gamma * gamma.zz)) / (sqrt(gammaU.xx) * gamma * gamma.zz))
	result[7] = (((((gammaU.xy * input[11] * (1. / sqrt(gammaU.xx))) - (gammaU.xy * input[16] * (1. / sqrt(gammaU.xx)))) - (gammaU.xz * input[10] * (1. / sqrt(gammaU.xx)))) + (gammaU.xz * input[15] * (1. / sqrt(gammaU.xx))) + (input[13] * sqrt(gammaU.xx)) + (input[18] * sqrt(gammaU.xx)) + ((2. * input[29]) - (2. * input[6]))) / (2. * sqrt(gammaU.xx)))
	result[8] = ((((gammaU.xy * input[10] * (1. / sqrt(gammaU.xx))) - (gammaU.xy * input[15] * (1. / sqrt(gammaU.xx)))) + (input[12] * sqrt(gammaU.xx)) + (input[28] - input[5])) / sqrt(gammaU.xx))
	result[9] = input[20]
	result[10] = input[19]
	result[11] = input[18]
	result[12] = input[17]
	result[13] = input[16]
	result[14] = input[15]
	result[15] = input[14]
	result[16] = input[13]
	result[17] = input[12]
	result[18] = input[11]
	result[19] = input[10]
	result[20] = input[9]
	result[21] = ((-((((gamma * gamma.zz * input[27] * gammaU.xx * gammaU.yy) - (gamma * gamma.zz * input[27] * (real)pow((real)gammaU.xy, (real)2.))) - (gamma * gamma.zz * input[3] * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.))) + (gamma * gamma.zz * input[3] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy) + ((gamma * gamma.zz * input[3] * f * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.)) - (gamma * gamma.zz * input[3] * f * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy)) + (((gamma * gamma.zz * input[4] * gammaU.xx * gammaU.yy) - (gamma * gamma.zz * input[4] * (real)pow((real)gammaU.xy, (real)2.))) - ((real)pow((real)gammaU.xx, (real)(5. / 2.)) * (real)pow((real)gammaU.yy, (real)2.) * input[26])) + ((real)pow((real)gammaU.xx, (real)(5. / 2.)) * (real)pow((real)gammaU.yy, (real)2.) * input[26] * f) + ((gammaU.xy * (real)pow((real)gammaU.yy, (real)2.) * input[25] * (real)pow((real)gammaU.xx, (real)(3. / 2.))) - ((real)pow((real)gammaU.xy, (real)3.) * gammaU.yy * input[25] * sqrt(gammaU.xx))) + (((real)pow((real)gammaU.xy, (real)3.) * gammaU.yy * input[25] * sqrt(gammaU.xx) * f) - (gammaU.xy * (real)pow((real)gammaU.yy, (real)2.) * input[25] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * f)) + ((((real)pow((real)gammaU.xy, (real)3.) * gammaU.yy * input[2] * sqrt(gammaU.xx)) - (gammaU.xy * (real)pow((real)gammaU.yy, (real)2.) * input[2] * (real)pow((real)gammaU.xx, (real)(3. / 2.)))) - ((real)pow((real)gammaU.xy, (real)3.) * gammaU.yy * input[2] * f * sqrt(gammaU.xx))) + ((gammaU.xy * (real)pow((real)gammaU.yy, (real)2.) * input[2] * f * (real)pow((real)gammaU.xx, (real)(3. / 2.))) - ((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[28] * gamma * gamma.zz)) + (((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[28] * f * gamma * gamma.zz) - ((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[5] * gamma * gamma.zz)) + ((real)pow((real)gammaU.xy, (real)2.) * gammaU.zz * input[5] * f * gamma * gamma.zz) + ((2. * (real)pow((real)gammaU.xy, (real)3.) * gamma.xz * input[29] * gamma) - (2. * gammaU.xy * gamma.xz * input[29] * gammaU.xx * gammaU.yy * gamma)) + ((2. * gammaU.xy * gamma.xz * input[29] * gammaU.xx * gammaU.yy * f * gamma) - (2. * (real)pow((real)gammaU.xy, (real)3.) * gamma.xz * input[29] * f * gamma)) + ((2. * (real)pow((real)gammaU.xy, (real)3.) * gamma.xz * input[6] * gamma) - (2. * gammaU.xy * gamma.xz * input[6] * gammaU.xx * gammaU.yy * gamma)) + ((2. * gammaU.xy * gamma.xz * input[6] * gammaU.xx * gammaU.yy * f * gamma) - (2. * (real)pow((real)gammaU.xy, (real)3.) * gamma.xz * input[6] * f * gamma)) + (((real)pow((real)gammaU.xy, (real)4.) * f * input[27]) - ((real)pow((real)gammaU.xy, (real)2.) * f * input[27] * gammaU.xx * gammaU.yy)) + ((((real)pow((real)gammaU.xy, (real)4.) * f * input[4]) - ((real)pow((real)gammaU.xy, (real)2.) * f * input[4] * gammaU.xx * gammaU.yy)) - (2. * (real)pow((real)gammaU.xy, (real)4.) * input[26] * sqrt(gammaU.xx))) + ((3. * (real)pow((real)gammaU.xy, (real)2.) * input[26] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy) - (3. * (real)pow((real)gammaU.xy, (real)2.) * input[26] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy * f)) + ((2. * (real)pow((real)gammaU.xy, (real)4.) * input[26] * sqrt(gammaU.xx) * f) - ((real)pow((real)gammaU.xy, (real)4.) * input[27])) + ((real)pow((real)gammaU.xy, (real)2.) * input[27] * gammaU.xx * gammaU.yy) + ((((real)pow((real)gammaU.xy, (real)4.) * input[3] * sqrt(gammaU.xx)) - ((real)pow((real)gammaU.xy, (real)2.) * input[3] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy)) - ((real)pow((real)gammaU.xy, (real)4.) * input[3] * f * sqrt(gammaU.xx))) + (((real)pow((real)gammaU.xy, (real)2.) * input[3] * f * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * gammaU.yy) - ((real)pow((real)gammaU.xy, (real)4.) * input[4])) + (((real)pow((real)gammaU.xy, (real)2.) * input[4] * gammaU.xx * gammaU.yy) - (gammaU.xz * (real)pow((real)gammaU.yy, (real)2.) * input[1] * (real)pow((real)gammaU.xx, (real)(3. / 2.)))) + (gammaU.xz * gammaU.yy * input[1] * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.)) + ((gammaU.xz * (real)pow((real)gammaU.yy, (real)2.) * input[1] * f * (real)pow((real)gammaU.xx, (real)(3. / 2.))) - (gammaU.xz * gammaU.yy * input[1] * f * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.))) + ((gammaU.xz * (real)pow((real)gammaU.yy, (real)2.) * input[24] * (real)pow((real)gammaU.xx, (real)(3. / 2.))) - (gammaU.xz * gammaU.yy * input[24] * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.))) + ((gammaU.xz * gammaU.yy * input[24] * sqrt(gammaU.xx) * (real)pow((real)gammaU.xy, (real)2.) * f) - (gammaU.xz * (real)pow((real)gammaU.yy, (real)2.) * input[24] * (real)pow((real)gammaU.xx, (real)(3. / 2.)) * f)) + (((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[28] * gamma * gamma.zz) - ((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[28] * f * gamma * gamma.zz)) + (((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[5] * gamma * gamma.zz) - ((real)pow((real)gammaU.xz, (real)2.) * gammaU.yy * input[5] * f * gamma * gamma.zz)) + ((f * gamma * gamma.zz * input[27] * gammaU.xx * gammaU.yy) - (f * gamma * gamma.zz * input[27] * (real)pow((real)gammaU.xy, (real)2.))) + (((f * gamma * gamma.zz * input[4] * gammaU.xx * gammaU.yy) - (f * gamma * gamma.zz * input[4] * (real)pow((real)gammaU.xy, (real)2.))) - (f * m * gamma * gamma.zz * input[27] * gammaU.xx * gammaU.yy)) + ((f * m * gamma * gamma.zz * input[27] * (real)pow((real)gammaU.xy, (real)2.)) - (f * m * gamma * gamma.zz * input[4] * gammaU.xx * gammaU.yy)) + ((f * m * gamma * gamma.zz * input[4] * (real)pow((real)gammaU.xy, (real)2.)) - (input[0] * gammaU.xx * gamma * gamma.zz * (real)pow((real)gammaU.xy, (real)2.))) + (input[0] * (real)pow((real)gammaU.xx, (real)2.) * gamma * gamma.zz * gammaU.yy) + (((input[0] * gammaU.xx * gamma * gamma.zz * f * (real)pow((real)gammaU.xy, (real)2.)) - (input[0] * (real)pow((real)gammaU.xx, (real)2.) * gamma * gamma.zz * f * gammaU.yy)) - (input[30] * gammaU.xx * (real)pow((real)gammaU.xy, (real)2.) * gamma * gamma.zz)) + (input[30] * gammaU.xx * (real)pow((real)gammaU.xy, (real)2.) * f * gamma * gamma.zz) + ((input[30] * (real)pow((real)gammaU.xx, (real)2.) * gammaU.yy * gamma * gamma.zz) - (input[30] * (real)pow((real)gammaU.xx, (real)2.) * gammaU.yy * f * gamma * gamma.zz)))) / (gammaU.xx * gamma * gamma.zz * ((((real)pow((real)gammaU.xy, (real)2.) - ((real)pow((real)gammaU.xy, (real)2.) * f)) - (gammaU.xx * gammaU.yy)) + (gammaU.xx * gammaU.yy * f))))
	result[22] = (((2. * sqrt(gammaU.xx) * gammaU.xy * gammaU.yz * input[29]) + ((2. * sqrt(gammaU.xx) * gammaU.xy * gammaU.yz * input[6]) - (gammaU.xx * gammaU.xy * input[26])) + ((gammaU.xx * gammaU.xy * input[3]) - (gammaU.xx * gammaU.yy * input[25])) + (gammaU.xx * gammaU.yy * input[2]) + ((gammaU.xy * gammaU.xz * input[1]) - (gammaU.xy * gammaU.xz * input[24])) + (gammaU.xy * gamma.yy * input[28] * (1. / sqrt(gammaU.xx)) * gamma) + (((((((gammaU.xy * gamma.yy * input[5] * (1. / sqrt(gammaU.xx)) * gamma) - (gammaU.xy * input[27] * sqrt(gammaU.xx))) - (gammaU.xy * input[4] * sqrt(gammaU.xx))) - (gammaU.xz * sqrt(gammaU.xx) * gammaU.yy * input[29])) - (gammaU.xz * sqrt(gammaU.xx) * gammaU.yy * input[6])) - (gammaU.xz * (real)pow((real)gammaU.xy, (real)2.) * input[29] * (1. / sqrt(gammaU.xx)))) - (gammaU.xz * (real)pow((real)gammaU.xy, (real)2.) * input[6] * (1. / sqrt(gammaU.xx))))) / (sqrt(gammaU.xx) * gamma * gamma.zz))
	result[23] = ((-((gammaU.xy * input[29] * (1. / sqrt(gammaU.xx))) + (gammaU.xy * input[6] * (1. / sqrt(gammaU.xx))) + (gammaU.xz * input[28] * (1. / sqrt(gammaU.xx))) + ((gammaU.xz * input[5] * (1. / sqrt(gammaU.xx))) - input[1]) + input[24])) / sqrt(gammaU.xx))
	result[24] = (((((sqrt(gammaU.xx) * gammaU.xy * input[25]) - (sqrt(gammaU.xx) * gammaU.xy * input[2])) - (sqrt(gammaU.xx) * gammaU.xz * input[1])) + (((sqrt(gammaU.xx) * gammaU.xz * input[24]) - (2. * gammaU.xx * gammaU.yz * input[29])) - (2. * gammaU.xx * gammaU.yz * input[6])) + ((real)pow((real)gammaU.xx, (real)(3. / 2.)) * input[26]) + ((gammaU.xx * input[27]) - ((real)pow((real)gammaU.xx, (real)(3. / 2.)) * input[3])) + (gammaU.xx * input[4]) + (2. * gammaU.xy * gammaU.xz * input[29]) + (((2. * gammaU.xy * gammaU.xz * input[6]) - (gamma.yy * input[28] * gamma)) - (gamma.yy * input[5] * gamma))) / (gamma * gamma.zz))
	result[25] = (input[29] + input[6])
	result[26] = (input[28] + input[5])
	result[27] = (input[27] + input[4])
	result[28] = ((-((gammaU.xy * input[8]) + (((gammaU.xz * input[7]) - (input[26] * gammaU.xx)) - (input[3] * gammaU.xx)))) / gammaU.xx)
	result[29] = (input[25] + input[2] + input[8])
	result[30] = (input[1] + input[24] + input[7])

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

	result[0] = (f * ((gammaU.xx * input[21]) + (2. * gammaU.xy * input[22]) + (2. * gammaU.xz * input[23]) + (gammaU.yy * input[24]) + (2. * gammaU.yz * input[25]) + ((gammaU.zz * input[26]) - (m * input[27]))))
	result[1] = 0.
	result[2] = 0.
	result[3] = input[21]
	result[4] = input[22]
	result[5] = input[23]
	result[6] = input[24]
	result[7] = input[25]
	result[8] = input[26]
	result[9] = 0.
	result[10] = 0.
	result[11] = 0.
	result[12] = 0.
	result[13] = 0.
	result[14] = 0.
	result[15] = 0.
	result[16] = 0.
	result[17] = 0.
	result[18] = 0.
	result[19] = 0.
	result[20] = 0.
	result[21] = (-(((gammaU.yy * input[10]) - (gammaU.yy * input[6])) + (gammaU.yz * input[11]) + ((gammaU.yz * input[16]) - (2. * gammaU.yz * input[7])) + (((gammaU.zz * input[17]) - (gammaU.zz * input[8])) - input[0]) + (2. * input[28])))
	result[22] = ((((2. * gammaU.xy * input[10]) - (2. * gammaU.xy * input[6])) + (gammaU.xz * input[11]) + ((gammaU.xz * input[16]) - (2. * gammaU.xz * input[7])) + ((gammaU.yz * input[13]) - (gammaU.yz * input[18])) + ((gammaU.zz * input[14]) - (gammaU.zz * input[19])) + (input[1] - (2. * input[29]))) / 2.)
	result[23] = (((gammaU.xy * input[11]) + ((gammaU.xy * input[16]) - (2. * gammaU.xy * input[7])) + (((2. * gammaU.xz * input[17]) - (2. * gammaU.xz * input[8])) - (gammaU.yy * input[13])) + ((gammaU.yy * input[18]) - (gammaU.yz * input[14])) + (gammaU.yz * input[19]) + (input[2] - (2. * input[30]))) / 2.)
	result[24] = (-(((gammaU.xx * input[10]) - (gammaU.xx * input[6])) + ((gammaU.xz * input[13]) - (gammaU.xz * input[18]))))
	result[25] = ((-((gammaU.xx * input[11]) + (((gammaU.xx * input[16]) - (2. * gammaU.xx * input[7])) - (gammaU.xy * input[13])) + (gammaU.xy * input[18]) + ((gammaU.xz * input[14]) - (gammaU.xz * input[19])))) / 2.)
	result[26] = (-((((gammaU.xx * input[17]) - (gammaU.xx * input[8])) - (gammaU.xy * input[14])) + (gammaU.xy * input[19])))
	result[27] = ((((gamma * gamma.xy * input[19]) - (gamma * gamma.xz * input[18])) - (gamma * gamma.yy * input[17])) + (gamma * gamma.yy * input[8]) + (gamma * gamma.yz * input[11]) + ((gamma * gamma.yz * input[16]) - (gamma * gamma.zz * input[10])) + (gamma * gamma.zz * input[6]) + (((2. * gammaU.xx * gammaU.yz * input[7]) - (gammaU.xx * input[28])) - (2. * gammaU.xy * gammaU.xz * input[7])) + (gammaU.xy * gammaU.yz * input[13]) + (((((gammaU.xy * gammaU.zz * input[14]) - (gammaU.xy * input[29])) - (gammaU.xz * gammaU.yy * input[13])) - (gammaU.xz * gammaU.yz * input[14])) - (gammaU.xz * input[30])))
	result[28] = ((gammaU.xy * input[22]) + (gammaU.xz * input[23]) + (gammaU.yy * input[24]) + (2. * gammaU.yz * input[25]) + ((gammaU.zz * input[26]) - input[27]))
	result[29] = (-((gammaU.xx * input[22]) + (gammaU.xy * input[24]) + (gammaU.xz * input[25])))
	result[30] = (-((gammaU.xx * input[23]) + (gammaU.xy * input[25]) + (gammaU.xz * input[26])))
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
