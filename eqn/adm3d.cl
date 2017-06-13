inline real calcMaxEigenvalue(real alpha, real gammaUii) {
	return alpha * sqrt(calc_f(alpha) * gammaUii);
}

kernel void calcDT(
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(2,2)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const global <?=eqn.cons_t?>* U = UBuf + index;
	real gamma = sym3_det(U->gamma);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		<? if side==0 then ?>
		real gammaUxx = (U->gamma.yy * U->gamma.zz - U->gamma.yz * U->gamma.yz) / gamma;
		real lambda = calcMaxEigenvalue(U->alpha, gammaUxx);
		<? elseif side==1 then ?>
		real gammaUyy = (U->gamma.xx * U->gamma.zz - U->gamma.xz * U->gamma.xz) / gamma;
		real lambda = calcMaxEigenvalue(U->alpha, gammaUyy);
		<? elseif side==2 then ?>
		real gammaUzz = (U->gamma.xx * U->gamma.yy - U->gamma.xy * U->gamma.xy) / gamma;
		lambda = calcMaxEigenvalue(U->alpha, gammaUzz);
		<? end ?>	
	
		real lambdaMin = (real)min((real)0., -lambda);
		real lambdaMax = (real)max((real)0., lambda);
		dt = (real)min(dt, (real)(dx<?=side?>_at(i) / (fabs(lambdaMax - lambdaMin) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt; 
}

//used by PLM
<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	<?=eqn.eigen_t?>* eig,
	const global <?=eqn.cons_t?>* U,
	real3 x 
) {
	eig->alpha = U->alpha;
	eig->f = calc_f(U->alpha);
	real gamma_det = sym3_det(U->gamma);
	eig->gammaU = sym3_inv(gamma_det, U->gamma);
}
<? end ?>

//used by PLM
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
	real lambdaLight = eig->alpha * sqrt(eig->gammaU.xx);
	<? elseif side==1 then ?>
	real lambdaLight = eig->alpha * sqrt(eig->gammaU.yy);
	<? elseif side==2 then ?>
	real lambdaLight = eig->alpha * sqrt(eig->gammaU.zz);
	<? end ?>
	real lambdaGauge = lambdaLight * sqrt(eig->f);
			
	wave[0] = -lambdaGauge;
	<? for i=1,5 do ?> wave[<?=i?>] = -lambdaLight; <? end ?>
	<? for i=6,23 do ?> wave[<?=i?>] = 0.; <? end ?>
	<? for i=24,28 do ?> wave[<?=i?>] = lambdaLight; <? end ?>
	wave[29] = lambdaGauge;
}
<?		end
	end
end ?>

//used for interface eigen basis
void eigen_forSide(
	global <?=eqn.eigen_t?>* eig,
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR
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
	real avg_gamma_det = sym3_det(avg_gamma);
	
	eig->alpha = alpha;
	eig->f = calc_f(alpha);
	eig->gammaU = sym3_inv(avg_gamma_det, avg_gamma);
}

kernel void calcEigenBasis(
	global real* waveBuf,
	global <?=eqn.eigen_t?>* eigenBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(2,1);
	real3 x = cell_x(i);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		<?= solver.getULRCode ?>	
		int indexInt = side + dim * index;	
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		eigen_forSide(eig, UL, UR);
		global real* wave = waveBuf + numWaves * indexInt;
		eigen_calcWaves_<?=side?>_global_global(wave, eig, x);
	}<? end ?>
}

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,solver.dim-1 do ?>
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 unused
) {
	real f = eig->f;
	sym3 gammaU = eig->gammaU;
	real sqrt_f = sqrt(f);
	
	//input of left eigenvectors is the state
	//so skip the first 7 of input
	input += 7;

	<? if side == 0 then ?>

	real sqrt_gammaUxx = sqrt(gammaU.xx);
	real gammaUxx_toThe_3_2 = sqrt_gammaUxx * gammaU.xx;
	
	results[0] = ((sqrt_f * gammaUxx_toThe_3_2 * input[21]) + (2. * sqrt_f * sqrt_gammaUxx * gammaU.xy * input[22]) + (2. * sqrt_f * sqrt_gammaUxx * gammaU.xz * input[23]) + (sqrt_f * sqrt_gammaUxx * gammaU.yy * input[24]) + (2. * sqrt_f * sqrt_gammaUxx * gammaU.yz * input[25]) + (((((((sqrt_f * sqrt_gammaUxx * gammaU.zz * input[26]) - (gammaU.xx * input[0])) - (2. * gammaU.xx * input[27])) - (gammaU.xy * input[1])) - (2. * gammaU.xy * input[28])) - (gammaU.xz * input[2])) - (2. * gammaU.xz * input[29])));
	results[1] = ((-(input[1] + (2. * input[28]) + ((((2. * gammaU.xx * input[4]) - (gammaU.xx * input[9])) - (2. * sqrt_gammaUxx * input[22])) - (2. * gammaU.xz * input[11])) + ((((2. * gammaU.xz * input[16]) - (gammaU.yy * input[12])) - (2. * gammaU.yz * input[13])) - (gammaU.zz * input[14])))) / 2.);
	results[2] = ((-(input[2] + (2. * input[29]) + (((2. * gammaU.xx * input[5]) - (gammaU.xx * input[15])) - (2. * sqrt_gammaUxx * input[23])) + (((((2. * gammaU.xy * input[11]) - (2. * gammaU.xy * input[16])) - (gammaU.yy * input[18])) - (2. * gammaU.yz * input[19])) - (gammaU.zz * input[20])))) / 2.);
	results[3] = (-(((gammaU.xx * input[6]) - (sqrt_gammaUxx * input[24])) + (gammaU.xy * input[12]) + (gammaU.xz * input[18])));
	results[4] = (-(((gammaU.xx * input[7]) - (sqrt_gammaUxx * input[25])) + (gammaU.xy * input[13]) + (gammaU.xz * input[19])));
	results[5] = (-(((gammaU.xx * input[8]) - (sqrt_gammaUxx * input[26])) + (gammaU.xy * input[14]) + (gammaU.xz * input[20])));
	results[6] = input[1];
	results[7] = input[2];
	results[8] = input[9];
	results[9] = input[10];
	results[10] = input[11];
	results[11] = input[12];
	results[12] = input[13];
	results[13] = input[14];
	results[14] = input[15];
	results[15] = input[16];
	results[16] = input[17];
	results[17] = input[18];
	results[18] = input[19];
	results[19] = input[20];
	results[20] = input[27];
	results[21] = input[28];
	results[22] = input[29];
	results[23] = ((((((input[0] - (f * gammaU.xx * input[3])) - (2. * f * gammaU.xy * input[4])) - (2. * f * gammaU.xz * input[5])) - (f * gammaU.yy * input[6])) - (2. * f * gammaU.yz * input[7])) - (f * gammaU.zz * input[8]));
	results[24] = ((input[1] + (2. * input[28]) + ((2. * gammaU.xx * input[4]) - (gammaU.xx * input[9])) + ((2. * sqrt_gammaUxx * input[22]) - (2. * gammaU.xz * input[11])) + ((((2. * gammaU.xz * input[16]) - (gammaU.yy * input[12])) - (2. * gammaU.yz * input[13])) - (gammaU.zz * input[14]))) / 2.);
	results[25] = ((input[2] + (2. * input[29]) + ((2. * gammaU.xx * input[5]) - (gammaU.xx * input[15])) + (2. * sqrt_gammaUxx * input[23]) + (((((2. * gammaU.xy * input[11]) - (2. * gammaU.xy * input[16])) - (gammaU.yy * input[18])) - (2. * gammaU.yz * input[19])) - (gammaU.zz * input[20]))) / 2.);
	results[26] = ((gammaU.xx * input[6]) + (sqrt_gammaUxx * input[24]) + (gammaU.xy * input[12]) + (gammaU.xz * input[18]));
	results[27] = ((gammaU.xx * input[7]) + (sqrt_gammaUxx * input[25]) + (gammaU.xy * input[13]) + (gammaU.xz * input[19]));
	results[28] = ((gammaU.xx * input[8]) + (sqrt_gammaUxx * input[26]) + (gammaU.xy * input[14]) + (gammaU.xz * input[20]));
	results[29] = ((sqrt_f * gammaUxx_toThe_3_2 * input[21]) + (2. * sqrt_f * sqrt_gammaUxx * gammaU.xy * input[22]) + (2. * sqrt_f * sqrt_gammaUxx * gammaU.xz * input[23]) + (sqrt_f * sqrt_gammaUxx * gammaU.yy * input[24]) + (2. * sqrt_f * sqrt_gammaUxx * gammaU.yz * input[25]) + (sqrt_f * sqrt_gammaUxx * gammaU.zz * input[26]) + (gammaU.xx * input[0]) + (2. * gammaU.xx * input[27]) + (gammaU.xy * input[1]) + (2. * gammaU.xy * input[28]) + (gammaU.xz * input[2]) + (2. * gammaU.xz * input[29]));
	
	<? elseif side == 1 then ?>
	
	real sqrt_gammaUyy = sqrt(gammaU.yy);
	real gammaUyy_toThe_3_2 = sqrt_gammaUyy * gammaU.yy;
	
	results[0] = ((sqrt_f * gammaUyy_toThe_3_2 * input[24]) + (sqrt_f * sqrt_gammaUyy * gammaU.xx * input[21]) + (2. * sqrt_f * sqrt_gammaUyy * gammaU.xy * input[22]) + (2. * sqrt_f * sqrt_gammaUyy * gammaU.xz * input[23]) + (2. * sqrt_f * sqrt_gammaUyy * gammaU.yz * input[25]) + (((((((sqrt_f * sqrt_gammaUyy * gammaU.zz * input[26]) - (gammaU.xy * input[0])) - (2. * gammaU.xy * input[27])) - (gammaU.yy * input[1])) - (2. * gammaU.yy * input[28])) - (gammaU.yz * input[2])) - (2. * gammaU.yz * input[29])));
	results[1] = (-((gammaU.xy * input[3]) + ((gammaU.yy * input[9]) - (sqrt_gammaUyy * input[21])) + (gammaU.yz * input[15])));
	results[2] = ((-(input[0] + ((((2. * input[27]) - (gammaU.xx * input[3])) - (2. * gammaU.xz * input[5])) - (gammaU.yy * input[6])) + (((2. * gammaU.yy * input[10]) - (2. * sqrt_gammaUyy * input[22])) - (2. * gammaU.yz * input[7])) + ((2. * gammaU.yz * input[16]) - (gammaU.zz * input[8])))) / 2.);
	results[3] = (-((gammaU.xy * input[5]) + ((gammaU.yy * input[11]) - (sqrt_gammaUyy * input[23])) + (gammaU.yz * input[17])));
	results[4] = ((-(input[2] + ((2. * input[29]) - (gammaU.xx * input[15])) + (((2. * gammaU.xy * input[7]) - (2. * gammaU.xy * input[16])) - (2. * gammaU.xz * input[17])) + ((((2. * gammaU.yy * input[13]) - (gammaU.yy * input[18])) - (2. * sqrt_gammaUyy * input[25])) - (gammaU.zz * input[20])))) / 2.);
	results[5] = (-((gammaU.xy * input[8]) + ((gammaU.yy * input[14]) - (sqrt_gammaUyy * input[26])) + (gammaU.yz * input[20])));
	results[6] = input[0];
	results[7] = input[2];
	results[8] = input[3];
	results[9] = input[4];
	results[10] = input[5];
	results[11] = input[6];
	results[12] = input[7];
	results[13] = input[8];
	results[14] = input[15];
	results[15] = input[16];
	results[16] = input[17];
	results[17] = input[18];
	results[18] = input[19];
	results[19] = input[20];
	results[20] = input[27];
	results[21] = input[28];
	results[22] = input[29];
	results[23] = ((((((input[1] - (f * gammaU.xx * input[9])) - (2. * f * gammaU.xy * input[10])) - (2. * f * gammaU.xz * input[11])) - (f * gammaU.yy * input[12])) - (2. * f * gammaU.yz * input[13])) - (f * gammaU.zz * input[14]));
	results[24] = ((gammaU.xy * input[3]) + (gammaU.yy * input[9]) + (sqrt_gammaUyy * input[21]) + (gammaU.yz * input[15]));
	results[25] = ((input[0] + ((((2. * input[27]) - (gammaU.xx * input[3])) - (2. * gammaU.xz * input[5])) - (gammaU.yy * input[6])) + (2. * gammaU.yy * input[10]) + ((2. * sqrt_gammaUyy * input[22]) - (2. * gammaU.yz * input[7])) + ((2. * gammaU.yz * input[16]) - (gammaU.zz * input[8]))) / 2.);
	results[26] = ((gammaU.xy * input[5]) + (gammaU.yy * input[11]) + (sqrt_gammaUyy * input[23]) + (gammaU.yz * input[17]));
	results[27] = ((input[2] + ((2. * input[29]) - (gammaU.xx * input[15])) + (((2. * gammaU.xy * input[7]) - (2. * gammaU.xy * input[16])) - (2. * gammaU.xz * input[17])) + ((2. * gammaU.yy * input[13]) - (gammaU.yy * input[18])) + ((2. * sqrt_gammaUyy * input[25]) - (gammaU.zz * input[20]))) / 2.);
	results[28] = ((gammaU.xy * input[8]) + (gammaU.yy * input[14]) + (sqrt_gammaUyy * input[26]) + (gammaU.yz * input[20]));
	results[29] = ((sqrt_f * gammaUyy_toThe_3_2 * input[24]) + (sqrt_f * sqrt_gammaUyy * gammaU.xx * input[21]) + (2. * sqrt_f * sqrt_gammaUyy * gammaU.xy * input[22]) + (2. * sqrt_f * sqrt_gammaUyy * gammaU.xz * input[23]) + (2. * sqrt_f * sqrt_gammaUyy * gammaU.yz * input[25]) + (sqrt_f * sqrt_gammaUyy * gammaU.zz * input[26]) + (gammaU.xy * input[0]) + (2. * gammaU.xy * input[27]) + (gammaU.yy * input[1]) + (2. * gammaU.yy * input[28]) + (gammaU.yz * input[2]) + (2. * gammaU.yz * input[29]));
	
	<? elseif side == 2 then ?>
	
	real sqrt_gammaUzz = sqrt(gammaU.zz);
	real gammaUzz_toThe_3_2 = sqrt_gammaUzz * gammaU.zz;
	
	results[0] = ((sqrt_f * gammaUzz_toThe_3_2 * input[26]) + (sqrt_f * sqrt_gammaUzz * gammaU.xx * input[21]) + (2. * sqrt_f * sqrt_gammaUzz * gammaU.xy * input[22]) + (2. * sqrt_f * sqrt_gammaUzz * gammaU.xz * input[23]) + (sqrt_f * sqrt_gammaUzz * gammaU.yy * input[24]) + (((((((2. * sqrt_f * sqrt_gammaUzz * gammaU.yz * input[25]) - (gammaU.xz * input[0])) - (2. * gammaU.xz * input[27])) - (gammaU.yz * input[1])) - (2. * gammaU.yz * input[28])) - (gammaU.zz * input[2])) - (2. * gammaU.zz * input[29])));
	results[1] = (-((gammaU.xz * input[3]) + (gammaU.yz * input[9]) + ((gammaU.zz * input[15]) - (sqrt_gammaUzz * input[21]))));
	results[2] = (-((gammaU.xz * input[4]) + (gammaU.yz * input[10]) + ((gammaU.zz * input[16]) - (sqrt_gammaUzz * input[22]))));
	results[3] = ((-(input[0] + (((((2. * input[27]) - (gammaU.xx * input[3])) - (2. * gammaU.xy * input[4])) - (gammaU.yy * input[6])) - (2. * gammaU.yz * input[7])) + ((2. * gammaU.yz * input[11]) - (gammaU.zz * input[8])) + ((2. * gammaU.zz * input[17]) - (2. * sqrt_gammaUzz * input[23])))) / 2.);
	results[4] = (-((gammaU.xz * input[6]) + (gammaU.yz * input[12]) + ((gammaU.zz * input[18]) - (sqrt_gammaUzz * input[24]))));
	results[5] = ((-(input[1] + (((2. * input[28]) - (gammaU.xx * input[9])) - (2. * gammaU.xy * input[10])) + ((((2. * gammaU.xz * input[7]) - (2. * gammaU.xz * input[11])) - (gammaU.yy * input[12])) - (gammaU.zz * input[14])) + ((2. * gammaU.zz * input[19]) - (2. * sqrt_gammaUzz * input[25])))) / 2.);
	results[6] = input[0];
	results[7] = input[1];
	results[8] = input[3];
	results[9] = input[4];
	results[10] = input[5];
	results[11] = input[6];
	results[12] = input[7];
	results[13] = input[8];
	results[14] = input[9];
	results[15] = input[10];
	results[16] = input[11];
	results[17] = input[12];
	results[18] = input[13];
	results[19] = input[14];
	results[20] = input[27];
	results[21] = input[28];
	results[22] = input[29];
	results[23] = ((((((input[2] - (f * gammaU.xx * input[15])) - (2. * f * gammaU.xy * input[16])) - (2. * f * gammaU.xz * input[17])) - (f * gammaU.yy * input[18])) - (2. * f * gammaU.yz * input[19])) - (f * gammaU.zz * input[20]));
	results[24] = ((gammaU.xz * input[3]) + (gammaU.yz * input[9]) + (gammaU.zz * input[15]) + (sqrt_gammaUzz * input[21]));
	results[25] = ((gammaU.xz * input[4]) + (gammaU.yz * input[10]) + (gammaU.zz * input[16]) + (sqrt_gammaUzz * input[22]));
	results[26] = ((input[0] + (((((2. * input[27]) - (gammaU.xx * input[3])) - (2. * gammaU.xy * input[4])) - (gammaU.yy * input[6])) - (2. * gammaU.yz * input[7])) + ((2. * gammaU.yz * input[11]) - (gammaU.zz * input[8])) + (2. * gammaU.zz * input[17]) + (2. * sqrt_gammaUzz * input[23])) / 2.);
	results[27] = ((gammaU.xz * input[6]) + (gammaU.yz * input[12]) + (gammaU.zz * input[18]) + (sqrt_gammaUzz * input[24]));
	results[28] = ((input[1] + (((2. * input[28]) - (gammaU.xx * input[9])) - (2. * gammaU.xy * input[10])) + ((((2. * gammaU.xz * input[7]) - (2. * gammaU.xz * input[11])) - (gammaU.yy * input[12])) - (gammaU.zz * input[14])) + (2. * gammaU.zz * input[19]) + (2. * sqrt_gammaUzz * input[25])) / 2.);
	results[29] = ((sqrt_f * gammaUzz_toThe_3_2 * input[26]) + (sqrt_f * sqrt_gammaUzz * gammaU.xx * input[21]) + (2. * sqrt_f * sqrt_gammaUzz * gammaU.xy * input[22]) + (2. * sqrt_f * sqrt_gammaUzz * gammaU.xz * input[23]) + (sqrt_f * sqrt_gammaUzz * gammaU.yy * input[24]) + (2. * sqrt_f * sqrt_gammaUzz * gammaU.yz * input[25]) + (gammaU.xz * input[0]) + (2. * gammaU.xz * input[27]) + (gammaU.yz * input[1]) + (2. * gammaU.yz * input[28]) + (gammaU.zz * input[2]) + (2. * gammaU.zz * input[29]));

	<? end ?>
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* results,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 unused
) {
	real f = eig->f;
	sym3 gammaU = eig->gammaU;
	real sqrt_f = sqrt(f);

	//write zeros to the alpha and gammaLL terms
	for (int i = 0; i < 7; ++i) {
		*results = 0;
		++results;
	}

	<? if side == 0 then ?>

	real sqrt_gammaUxx = sqrt(gammaU.xx);
	real gammaUxx_toThe_3_2 = sqrt_gammaUxx * gammaU.xx;
	real gammaUxxSq = gammaU.xx * gammaU.xx;
	
	results[0] = ((input[0] + ((4. * input[20] * gammaU.xx) - input[29]) + (2. * gammaU.xy * input[6]) + (4. * gammaU.xy * input[21]) + (2. * gammaU.xz * input[7]) + (4. * gammaU.xz * input[22])) / (-(2. * gammaU.xx)));
	results[1] = input[6];
	results[2] = input[7];
	results[3] = ((-(input[0] + (4. * input[20] * gammaU.xx) + (((2. * input[23] * gammaU.xx) - input[29]) - (2. * gammaU.xy * input[1] * f)) + (2. * gammaU.xy * input[6]) + (2. * gammaU.xy * input[8] * gammaU.xx * f) + (4. * gammaU.xy * input[21]) + ((((2. * gammaU.xy * input[24] * f) - (2. * gammaU.xy * f * input[6])) - (4. * gammaU.xy * f * input[21])) - (2. * gammaU.xz * input[2] * f)) + (2. * gammaU.xz * input[7]) + (2. * gammaU.xz * input[14] * gammaU.xx * f) + (4. * gammaU.xz * input[22]) + ((((2. * gammaU.xz * input[25] * f) - (2. * gammaU.xz * f * input[7])) - (4. * gammaU.xz * f * input[22])) - (gammaU.yy * input[3] * f)) + ((gammaU.yy * input[26] * f) - (2. * gammaU.yz * input[4] * f)) + ((2. * gammaU.yz * input[27] * f) - (gammaU.zz * input[5] * f)) + (gammaU.zz * input[28] * f))) / (2. * gammaUxxSq * f));
	results[4] = ((input[1] + (input[6] - (input[8] * gammaU.xx)) + (((2. * input[21]) - input[24]) - (2. * gammaU.xz * input[10])) + ((((2. * gammaU.xz * input[15]) - (gammaU.yy * input[11])) - (2. * gammaU.yz * input[12])) - (gammaU.zz * input[13]))) / (-(2. * gammaU.xx)));
	results[5] = ((input[2] + (input[7] - (input[14] * gammaU.xx)) + ((2. * input[22]) - input[25]) + (((((2. * gammaU.xy * input[10]) - (2. * gammaU.xy * input[15])) - (gammaU.yy * input[17])) - (2. * gammaU.yz * input[18])) - (gammaU.zz * input[19]))) / (-(2. * gammaU.xx)));
	results[6] = (((input[3] - input[26]) + (2. * gammaU.xy * input[11]) + (2. * gammaU.xz * input[17])) / (-(2. * gammaU.xx)));
	results[7] = (((input[4] - input[27]) + (2. * gammaU.xy * input[12]) + (2. * gammaU.xz * input[18])) / (-(2. * gammaU.xx)));
	results[8] = (((input[5] - input[28]) + (2. * gammaU.xy * input[13]) + (2. * gammaU.xz * input[19])) / (-(2. * gammaU.xx)));
	results[9] = input[8];
	results[10] = input[9];
	results[11] = input[10];
	results[12] = input[11];
	results[13] = input[12];
	results[14] = input[13];
	results[15] = input[14];
	results[16] = input[15];
	results[17] = input[16];
	results[18] = input[17];
	results[19] = input[18];
	results[20] = input[19];
	results[21] = ((input[0] + ((((((((((input[29] - (2. * gammaU.xy * input[1] * sqrt_f)) - (2. * gammaU.xy * input[24] * sqrt_f)) - (2. * gammaU.xz * input[2] * sqrt_f)) - (2. * gammaU.xz * input[25] * sqrt_f)) - (gammaU.yy * input[3] * sqrt_f)) - (gammaU.yy * input[26] * sqrt_f)) - (2. * gammaU.yz * input[4] * sqrt_f)) - (2. * gammaU.yz * input[27] * sqrt_f)) - (gammaU.zz * input[5] * sqrt_f)) - (gammaU.zz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUxx_toThe_3_2));
	results[22] = ((input[1] + input[24]) / (2. * sqrt_gammaUxx));
	results[23] = ((input[2] + input[25]) / (2. * sqrt_gammaUxx));
	results[24] = ((input[3] + input[26]) / (2. * sqrt_gammaUxx));
	results[25] = ((input[4] + input[27]) / (2. * sqrt_gammaUxx));
	results[26] = ((input[5] + input[28]) / (2. * sqrt_gammaUxx));
	results[27] = input[20];
	results[28] = input[21];
	results[29] = input[22];

	<? elseif side == 1 then ?>

	real sqrt_gammaUyy = sqrt(gammaU.yy);
	real gammaUyy_toThe_3_2 = sqrt_gammaUyy * gammaU.yy;
	real gammaUyySq = gammaU.yy * gammaU.yy;
	
	results[0] = input[6];
	results[1] = ((-(input[0] + ((4. * input[21] * gammaU.yy) - input[29]) + (2. * gammaU.xy * input[6]) + (4. * gammaU.xy * input[20]) + (2. * gammaU.yz * input[7]) + (4. * gammaU.yz * input[22]))) / (2. * gammaU.yy));
	results[2] = input[7];
	results[3] = input[8];
	results[4] = input[9];
	results[5] = input[10];
	results[6] = input[11];
	results[7] = input[12];
	results[8] = input[13];
	results[9] = (((input[1] - input[24]) + (2. * gammaU.xy * input[8]) + (2. * gammaU.yz * input[14])) / (-(2. * gammaU.yy)));
	results[10] = ((-(input[2] + (input[6] - (input[11] * gammaU.yy)) + (((((2. * input[20]) - input[25]) - (gammaU.xx * input[8])) - (2. * gammaU.xz * input[10])) - (2. * gammaU.yz * input[12])) + ((2. * gammaU.yz * input[15]) - (gammaU.zz * input[13])))) / (2. * gammaU.yy));
	results[11] = (((input[3] - input[26]) + (2. * gammaU.xy * input[10]) + (2. * gammaU.yz * input[16])) / (-(2. * gammaU.yy)));
	results[12] = ((-(input[0] + (4. * input[21] * gammaU.yy) + (((2. * input[23] * gammaU.yy) - input[29]) - (gammaU.xx * input[1] * f)) + ((gammaU.xx * input[24] * f) - (2. * gammaU.xy * input[2] * f)) + (2. * gammaU.xy * input[6]) + (2. * gammaU.xy * input[11] * gammaU.yy * f) + (4. * gammaU.xy * input[20]) + ((((2. * gammaU.xy * input[25] * f) - (2. * gammaU.xy * f * input[6])) - (4. * gammaU.xy * f * input[20])) - (2. * gammaU.xz * input[3] * f)) + ((2. * gammaU.xz * input[26] * f) - (2. * gammaU.yz * input[4] * f)) + (2. * gammaU.yz * input[7]) + (2. * gammaU.yz * input[17] * gammaU.yy * f) + (4. * gammaU.yz * input[22]) + ((((2. * gammaU.yz * input[27] * f) - (2. * gammaU.yz * f * input[7])) - (4. * gammaU.yz * f * input[22])) - (gammaU.zz * input[5] * f)) + (gammaU.zz * input[28] * f))) / (2. * gammaUyySq * f));
	results[13] = ((input[4] + (input[7] - (input[17] * gammaU.yy)) + (((2. * input[22]) - input[27]) - (gammaU.xx * input[14])) + ((((2. * gammaU.xy * input[12]) - (2. * gammaU.xy * input[15])) - (2. * gammaU.xz * input[16])) - (gammaU.zz * input[19]))) / (-(2. * gammaU.yy)));
	results[14] = (((input[5] - input[28]) + (2. * gammaU.xy * input[13]) + (2. * gammaU.yz * input[19])) / (-(2. * gammaU.yy)));
	results[15] = input[14];
	results[16] = input[15];
	results[17] = input[16];
	results[18] = input[17];
	results[19] = input[18];
	results[20] = input[19];
	results[21] = ((input[1] + input[24]) / (2. * sqrt_gammaUyy));
	results[22] = ((input[2] + input[25]) / (2. * sqrt_gammaUyy));
	results[23] = ((input[3] + input[26]) / (2. * sqrt_gammaUyy));
	results[24] = ((input[0] + ((((((((((input[29] - (gammaU.xx * input[1] * sqrt_f)) - (gammaU.xx * input[24] * sqrt_f)) - (2. * gammaU.xy * input[2] * sqrt_f)) - (2. * gammaU.xy * input[25] * sqrt_f)) - (2. * gammaU.xz * input[3] * sqrt_f)) - (2. * gammaU.xz * input[26] * sqrt_f)) - (2. * gammaU.yz * input[4] * sqrt_f)) - (2. * gammaU.yz * input[27] * sqrt_f)) - (gammaU.zz * input[5] * sqrt_f)) - (gammaU.zz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUyy_toThe_3_2));
	results[25] = ((input[4] + input[27]) / (2. * sqrt_gammaUyy));
	results[26] = ((input[5] + input[28]) / (2. * sqrt_gammaUyy));
	results[27] = input[20];
	results[28] = input[21];
	results[29] = input[22];
	
	<? elseif side == 2 then ?>

	real sqrt_gammaUzz = sqrt(gammaU.zz);
	real gammaUzz_toThe_3_2 = sqrt_gammaUzz * gammaU.zz;
	real gammaUzzSq = gammaU.zz * gammaU.zz;
	
	results[0] = input[6];
	results[1] = input[7];
	results[2] = ((-(input[0] + ((4. * input[22] * gammaU.zz) - input[29]) + (2. * gammaU.xz * input[6]) + (4. * gammaU.xz * input[20]) + (2. * gammaU.yz * input[7]) + (4. * gammaU.yz * input[21]))) / (2. * gammaU.zz));
	results[3] = input[8];
	results[4] = input[9];
	results[5] = input[10];
	results[6] = input[11];
	results[7] = input[12];
	results[8] = input[13];
	results[9] = input[14];
	results[10] = input[15];
	results[11] = input[16];
	results[12] = input[17];
	results[13] = input[18];
	results[14] = input[19];
	results[15] = (((input[1] - input[24]) + (2. * gammaU.xz * input[8]) + (2. * gammaU.yz * input[14])) / (-(2. * gammaU.zz)));
	results[16] = (((input[2] - input[25]) + (2. * gammaU.xz * input[9]) + (2. * gammaU.yz * input[15])) / (-(2. * gammaU.zz)));
	results[17] = ((-(input[3] + (input[6] - (input[13] * gammaU.zz)) + ((((((2. * input[20]) - input[26]) - (gammaU.xx * input[8])) - (2. * gammaU.xy * input[9])) - (gammaU.yy * input[11])) - (2. * gammaU.yz * input[12])) + (2. * gammaU.yz * input[16]))) / (2. * gammaU.zz));
	results[18] = (((input[4] - input[27]) + (2. * gammaU.xz * input[11]) + (2. * gammaU.yz * input[17])) / (-(2. * gammaU.zz)));
	results[19] = ((-(input[5] + (input[7] - (input[19] * gammaU.zz)) + ((((2. * input[21]) - input[28]) - (gammaU.xx * input[14])) - (2. * gammaU.xy * input[15])) + (((2. * gammaU.xz * input[12]) - (2. * gammaU.xz * input[16])) - (gammaU.yy * input[17])))) / (2. * gammaU.zz));
	results[20] = ((-(input[0] + (4. * input[22] * gammaU.zz) + (((2. * input[23] * gammaU.zz) - input[29]) - (gammaU.xx * input[1] * f)) + ((gammaU.xx * input[24] * f) - (2. * gammaU.xy * input[2] * f)) + ((2. * gammaU.xy * input[25] * f) - (2. * gammaU.xz * input[3] * f)) + (2. * gammaU.xz * input[6]) + (2. * gammaU.xz * input[13] * gammaU.zz * f) + (4. * gammaU.xz * input[20]) + ((((2. * gammaU.xz * input[26] * f) - (2. * gammaU.xz * f * input[6])) - (4. * gammaU.xz * f * input[20])) - (gammaU.yy * input[4] * f)) + ((gammaU.yy * input[27] * f) - (2. * gammaU.yz * input[5] * f)) + (2. * gammaU.yz * input[7]) + (2. * gammaU.yz * input[19] * gammaU.zz * f) + (4. * gammaU.yz * input[21]) + (((2. * gammaU.yz * input[28] * f) - (2. * gammaU.yz * f * input[7])) - (4. * gammaU.yz * f * input[21])))) / (2. * gammaUzzSq * f));
	results[21] = ((input[1] + input[24]) / (2. * sqrt_gammaUzz));
	results[22] = ((input[2] + input[25]) / (2. * sqrt_gammaUzz));
	results[23] = ((input[3] + input[26]) / (2. * sqrt_gammaUzz));
	results[24] = ((input[4] + input[27]) / (2. * sqrt_gammaUzz));
	results[25] = ((input[5] + input[28]) / (2. * sqrt_gammaUzz));
	results[26] = ((input[0] + ((((((((((input[29] - (gammaU.xx * input[1] * sqrt_f)) - (gammaU.xx * input[24] * sqrt_f)) - (2. * gammaU.xy * input[2] * sqrt_f)) - (2. * gammaU.xy * input[25] * sqrt_f)) - (2. * gammaU.xz * input[3] * sqrt_f)) - (2. * gammaU.xz * input[26] * sqrt_f)) - (gammaU.yy * input[4] * sqrt_f)) - (gammaU.yy * input[27] * sqrt_f)) - (2. * gammaU.yz * input[5] * sqrt_f)) - (2. * gammaU.yz * input[28] * sqrt_f))) / (2. * sqrt_f * gammaUzz_toThe_3_2));
	results[27] = input[20];
	results[28] = input[21];
	results[29] = input[22];

	<? end ?>
}
<?				if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x,
	real3 unused
) {
	for (int i = 0; i < numStates; ++i) {
		*y = 0;
		++y;
	}
}
<?				end
			end
		end
	end
end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf)
{
	SETBOUNDS(2,2);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	global <?=eqn.cons_t?>* deriv = derivBuf + index;

	real gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(gamma, U->gamma);
	real f = calc_f(U->alpha);	//could be based on alpha...

real density = 0;//state[STATE_DENSITY];
real pressure = 0;//state[STATE_PRESSURE];
real3 vel3 = _real3(0,0,0);//(real4)(state[STATE_VELOCITY_X], state[STATE_VELOCITY_Y], state[STATE_VELOCITY_Z], 0.);
real vel3Sq = real3_dot(vel3, vel3);	//TODO use gamma
real LorentzFactor = 1. / sqrt(1. - vel3Sq);
real4 vel4_ = (real4)(vel3.x * LorentzFactor, vel3.y * LorentzFactor, vel3.z * LorentzFactor, LorentzFactor);

real3 beta_ = _real3(0,0,0);

// source terms
real KUL[3][3] = {
{gammaU.xx * U->K.xx + gammaU.xy * U->K.xy + gammaU.xz * U->K.xz,
gammaU.xx * U->K.xy + gammaU.xy * U->K.yy + gammaU.xz * U->K.yz,
gammaU.xx * U->K.xz + gammaU.xy * U->K.yz + gammaU.xz * U->K.zz,
},{gammaU.xy * U->K.xx + gammaU.yy * U->K.xy + gammaU.yz * U->K.xz,
gammaU.xy * U->K.xy + gammaU.yy * U->K.yy + gammaU.yz * U->K.yz,
gammaU.xy * U->K.xz + gammaU.yy * U->K.yz + gammaU.yz * U->K.zz,
},{gammaU.xz * U->K.xx + gammaU.yz * U->K.xy + gammaU.zz * U->K.xz,
gammaU.xz * U->K.xy + gammaU.yz * U->K.yy + gammaU.zz * U->K.yz,
gammaU.xz * U->K.xz + gammaU.yz * U->K.yz + gammaU.zz * U->K.zz,
},};
real trK = KUL[0][0] + KUL[1][1] + KUL[2][2];
real KSqSymLL[6] = {
U->K.xx * KUL[0][0] + U->K.xy * KUL[1][0] + U->K.xz * KUL[2][0],
U->K.xx * KUL[0][1] + U->K.xy * KUL[1][1] + U->K.xz * KUL[2][1],
U->K.xx * KUL[0][2] + U->K.xy * KUL[1][2] + U->K.xz * KUL[2][2],
U->K.xy * KUL[0][1] + U->K.yy * KUL[1][1] + U->K.yz * KUL[2][1],
U->K.xy * KUL[0][2] + U->K.yy * KUL[1][2] + U->K.yz * KUL[2][2],
U->K.xz * KUL[0][2] + U->K.yz * KUL[1][2] + U->K.zz * KUL[2][2],
};
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
8. * M_PI * ((density + pressure) * vel4_.x * vel4_.x + .5 * (density - pressure) * U->gamma.xx),
8. * M_PI * ((density + pressure) * vel4_.x * vel4_.y + .5 * (density - pressure) * U->gamma.xy),
8. * M_PI * ((density + pressure) * vel4_.x * vel4_.z + .5 * (density - pressure) * U->gamma.xz),
8. * M_PI * ((density + pressure) * vel4_.y * vel4_.y + .5 * (density - pressure) * U->gamma.yy),
8. * M_PI * ((density + pressure) * vel4_.y * vel4_.z + .5 * (density - pressure) * U->gamma.yz),
8. * M_PI * ((density + pressure) * vel4_.z * vel4_.z + .5 * (density - pressure) * U->gamma.zz),
};
real SSymLL[6] = {
-R4SymLL[0] + trK * U->K.xx - 2 * KSqSymLL[0] + 4 * D12SymLL[0] + Gamma31SymLL[0] - Gamma11SymLL[0] + ADDSymLL[0] + (U->a.x * ((2 * U->V.x) - D1L[0])),
-R4SymLL[1] + trK * U->K.xy - 2 * KSqSymLL[1] + 4 * D12SymLL[1] + Gamma31SymLL[1] - Gamma11SymLL[1] + ADDSymLL[1] + ((((2 * U->a.y * U->V.x) - (U->a.y * D1L[0])) + ((2 * U->a.x * U->V.y) - (U->a.x * D1L[1]))) / 2),
-R4SymLL[2] + trK * U->K.xz - 2 * KSqSymLL[2] + 4 * D12SymLL[2] + Gamma31SymLL[2] - Gamma11SymLL[2] + ADDSymLL[2] + ((((2 * U->a.z * U->V.x) - (U->a.z * D1L[0])) + ((2 * U->a.x * U->V.z) - (U->a.x * D1L[2]))) / 2),
-R4SymLL[3] + trK * U->K.yy - 2 * KSqSymLL[3] + 4 * D12SymLL[3] + Gamma31SymLL[3] - Gamma11SymLL[3] + ADDSymLL[3] + (U->a.y * ((2 * U->V.y) - D1L[1])),
-R4SymLL[4] + trK * U->K.yz - 2 * KSqSymLL[4] + 4 * D12SymLL[4] + Gamma31SymLL[4] - Gamma11SymLL[4] + ADDSymLL[4] + ((((2 * U->a.z * U->V.y) - (U->a.z * D1L[1])) + ((2 * U->a.y * U->V.z) - (U->a.y * D1L[2]))) / 2),
-R4SymLL[5] + trK * U->K.zz - 2 * KSqSymLL[5] + 4 * D12SymLL[5] + Gamma31SymLL[5] - Gamma11SymLL[5] + ADDSymLL[5] + (U->a.z * ((2 * U->V.z) - D1L[2])),
};
real GU0L[3] = {
8. * M_PI * ((density + pressure) * vel4_.w * vel4_.x + pressure * beta_.x),
8. * M_PI * ((density + pressure) * vel4_.w * vel4_.y + pressure * beta_.y),
8. * M_PI * ((density + pressure) * vel4_.w * vel4_.z + pressure * beta_.z),
};
real AKL[3] = {
U->a.x * KUL[0][0] + U->a.y * KUL[1][0] + U->a.z * KUL[2][0],
U->a.x * KUL[0][1] + U->a.y * KUL[1][1] + U->a.z * KUL[2][1],
U->a.x * KUL[0][2] + U->a.y * KUL[1][2] + U->a.z * KUL[2][2],
};
real K12D23L[3] = {
KUL[0][0] * DLUL[0][0][0] +KUL[0][1] * DLUL[0][1][0] +KUL[0][2] * DLUL[0][2][0] + KUL[1][0] * DLUL[0][0][1] +KUL[1][1] * DLUL[0][1][1] +KUL[1][2] * DLUL[0][2][1] + KUL[2][0] * DLUL[0][0][2] +KUL[2][1] * DLUL[0][1][2] +KUL[2][2] * DLUL[0][2][2],
KUL[0][0] * DLUL[1][0][0] +KUL[0][1] * DLUL[1][1][0] +KUL[0][2] * DLUL[1][2][0] + KUL[1][0] * DLUL[1][0][1] +KUL[1][1] * DLUL[1][1][1] +KUL[1][2] * DLUL[1][2][1] + KUL[2][0] * DLUL[1][0][2] +KUL[2][1] * DLUL[1][1][2] +KUL[2][2] * DLUL[1][2][2],
KUL[0][0] * DLUL[2][0][0] +KUL[0][1] * DLUL[2][1][0] +KUL[0][2] * DLUL[2][2][0] + KUL[1][0] * DLUL[2][0][1] +KUL[1][1] * DLUL[2][1][1] +KUL[1][2] * DLUL[2][2][1] + KUL[2][0] * DLUL[2][0][2] +KUL[2][1] * DLUL[2][1][2] +KUL[2][2] * DLUL[2][2][2],
};
real KD23L[3] = {
KUL[0][0] * D1L[0] + KUL[1][0] * D1L[1] + KUL[2][0] * D1L[2],
KUL[0][1] * D1L[0] + KUL[1][1] * D1L[1] + KUL[2][1] * D1L[2],
KUL[0][2] * D1L[0] + KUL[1][2] * D1L[1] + KUL[2][2] * D1L[2],
};
real K12D12L[3] = {
KUL[0][0] * DLUL[0][0][0] + KUL[0][1] * DLUL[0][1][0] + KUL[0][2] * DLUL[0][2][0] + KUL[1][0] * DLUL[1][0][0] + KUL[1][1] * DLUL[1][1][0] + KUL[1][2] * DLUL[1][2][0] + KUL[2][0] * DLUL[2][0][0] + KUL[2][1] * DLUL[2][1][0] + KUL[2][2] * DLUL[2][2][0],
KUL[0][0] * DLUL[0][0][1] + KUL[0][1] * DLUL[0][1][1] + KUL[0][2] * DLUL[0][2][1] + KUL[1][0] * DLUL[1][0][1] + KUL[1][1] * DLUL[1][1][1] + KUL[1][2] * DLUL[1][2][1] + KUL[2][0] * DLUL[2][0][1] + KUL[2][1] * DLUL[2][1][1] + KUL[2][2] * DLUL[2][2][1],
KUL[0][0] * DLUL[0][0][2] + KUL[0][1] * DLUL[0][1][2] + KUL[0][2] * DLUL[0][2][2] + KUL[1][0] * DLUL[1][0][2] + KUL[1][1] * DLUL[1][1][2] + KUL[1][2] * DLUL[1][2][2] + KUL[2][0] * DLUL[2][0][2] + KUL[2][1] * DLUL[2][1][2] + KUL[2][2] * DLUL[2][2][2],
};
real KD12L[3] = {
KUL[0][0] * D3L[0] + KUL[1][0] * D3L[1] + KUL[2][0] * D3L[2],
KUL[0][1] * D3L[0] + KUL[1][1] * D3L[1] + KUL[2][1] * D3L[2],
KUL[0][2] * D3L[0] + KUL[1][2] * D3L[1] + KUL[2][2] * D3L[2],
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

// the 1D version has no problems, but at 2D we get instabilities ... 
kernel void constrain(
	global <?=eqn.cons_t?>* UBuf
) {
<? if false then ?>
#if 0	//use constraints at all?
	SETBOUNDS(2,2);	
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(gamma, U->gamma);

	real D3_D1_x = 
		(gammaU.xy * U->d[0].xy)
		+ (gammaU.xz * U->d[0].xz)
		+ (gammaU.yy * U->d[0].yy)
		+ (2. * gammaU.yz * U->d[0].yz)
		+ (gammaU.zz * U->d[0].zz)
		- (gammaU.xy * U->d[1].xx)
		- (gammaU.xz * U->d[2].xx)
		- (gammaU.yy * U->d[1].xy)
		- (gammaU.yz * U->d[2].xy)
		- (gammaU.yz * U->d[1].xz)
		- (gammaU.zz * U->d[2].xz);
	real D3_D1_y = 
		(gammaU.xx * U->d[1].xx)
		+ (gammaU.xy * U->d[1].xy)
		+ (2. * gammaU.xz * U->d[1].xz)
		+ (gammaU.yz * U->d[1].yz)
		+ (gammaU.zz * U->d[1].zz)
		- (gammaU.xx * U->d[0].xy)
		- (gammaU.xz * U->d[2].xy)
		- (gammaU.xy * U->d[0].yy)
		- (gammaU.yz * U->d[2].yy)
		- (gammaU.xz * U->d[0].yz)
		- (gammaU.zz * U->d[2].yz);
	real D3_D1_z = 
		(gammaU.xx * U->d[2].xx)
		+ (2. * gammaU.xy * U->d[2].xy)
		+ (gammaU.xz * U->d[2].xz)
		+ (gammaU.yy * U->d[2].yy)
		+ (gammaU.yz * U->d[2].yz)
		- (gammaU.xx * U->d[0].xz)
		- (gammaU.xy * U->d[1].xz)
		- (gammaU.xy * U->d[0].yz)
		- (gammaU.yy * U->d[1].yz)
		- (gammaU.xz * U->d[0].zz)
		- (gammaU.yz * U->d[1].zz);

#if 0	//directly assign V_i's
	U->V.x = D3_D1_x;
	U->V.y = D3_D1_y;
	U->V.z = D3_D1_z;
#endif
#if 0	//linearly project out the [V_i, U->d_ijk] vector
#endif
#if 0	//do a single gradient descent step
/*
V_i = d_ik^k - d^k_ki
f_i = V_i - (d_ijk - d_jki) gamma^jk
*/
	U->V.x -= epsilon;
	U->V.y -= epsilon;
#endif

#endif
<? end ?>
}
