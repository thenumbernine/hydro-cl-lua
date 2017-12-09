#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x,<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	
	real3 B_u = U.B;
	real3 epsE_u = U.epsE;

	//this only works for fluxFromCons
	real3 B_l = sym3_real3_mul(gamma, B_u);
	real3 epsE_l = sym3_real3_mul(gamma, epsE_u);
	
	real mu = U.mu;
	real eps = U.eps;
	
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.epsE = _real3(0., B_l.z / mu, -B_l.y / mu),
		.B = _real3(0., -epsE_l.z / eps, epsE_l.y / eps),
	<? elseif side == 1 then ?>
		.epsE = _real3(-B_l.z / mu, 0., B_l.x / mu),
		.B = _real3(epsE_l.z / eps, 0., -epsE_l.x / eps),
	<? elseif side == 2 then ?>
		.epsE = _real3(B_l.y / mu, -B_l.x / mu, 0.),
		.B = _real3(-epsE_l.y / eps, epsE_l.x / eps, 0.),
	<? end ?>
		.BPot = 0.,
		.sigma = 0.,
		.eps = 0.,
		.mu = 0.,
	};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	
	real det_gamma = sym3_det(gamma);	
	real det_gamma2 = det_gamma * det_gamma;
	real det_gamma3 = det_gamma * det_gamma2;
	
	<? if side == 0 then ?>
	real detg_gUjj = gamma.yy * gamma.zz - gamma.yz * gamma.yz;
	<? elseif side == 1 then ?>
	real detg_gUjj = gamma.xx * gamma.zz - gamma.xz * gamma.xz;
	<? elseif side == 2 then ?>
	real detg_gUjj = gamma.xx * gamma.yy - gamma.xy * gamma.xy;
	<? end ?>

	real lambda = alpha * sqrt(detg_gUjj / (det_gamma3 * U->eps * U->mu));
	return (range_t){-lambda, lambda};
}
<? end ?>

//TODO HLL needs eigen_forSide, 
//but it would have to pass the extra ADM args into it
/*
<?=eqn.eigen_t?> eigen_forSide(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
}
*/

//TODO this is exactly like calcCellMinMax, so just have this call it
// ... but then the addr space of U ...
<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr1?> real* wave,
	const <?=addr0?> <?=eqn.eigen_t?>* eig,
	real3 x
) {
	real det_gamma = eig->det_gamma;
	real det_gamma2 = det_gamma * det_gamma;
	real det_gamma3 = det_gamma * det_gamma2;
	
	real lambda = eig->alpha / sqrt(eig->detg_gUjj / (det_gamma3 * eig->eps * eig->mu));
	
	wave[0] = -lambda;
	wave[1] = -lambda;
	wave[2] = 0;
	wave[3] = 0;
	wave[4] = lambda;
	wave[5] = lambda;
}
<?
		end
	end
end
?>

//same as in eqn/euler.cl
kernel void calcEigenBasis(
	global real* waveBuf,
	global <?=eqn.eigen_t?>* eigenBuf,
	<?= solver.getULRArg ?><?=
	solver:getADMArgs()?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;
	
	<?=solver:getADMVarCode{suffix='R'} --[[ produce alphaR, betaR, gammaR at indexR ]] ?>
	real det_gammaR = sym3_det(gammaR);

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
		
		<?= solver.getULRCode ?>
		
		<?=solver:getADMVarCode{suffix='L'} --[[ produce alphaL, betaL, gammaL at indexL ]] ?>
		real det_gammaL = sym3_det(gammaL);
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		//*eig = eigen_forSide(UL, UR, xInt);
		eig->eps = .5 * (UL->eps + UR->eps),
		eig->mu = .5 * (UL->mu + UR->mu),
		eig->alpha = .5 * (alphaL + alphaR);
		eig->det_gamma = .5 * (det_gammaL + det_gammaR);

		<? if side == 0 then ?>
		eig->detg_gUjj = .5 * (
			gammaL.yy * gammaL.zz - gammaL.yz * gammaL.yz
			+ gammaR.yy * gammaR.zz - gammaR.yz * gammaR.yz
		);
		<? elseif side == 1 then ?>
		eig->detg_gUjj = .5 * (
			gammaL.xx * gammaL.zz - gammaL.xz * gammaL.xz
			+ gammaR.xx * gammaR.zz - gammaR.xz * gammaR.xz
		);
		<? elseif side == 2 then ?>
		eig->detg_gUjj = .5 * (
			gammaL.xx * gammaL.yy - gammaL.xy * gammaL.xy
			+ gammaR.xx * gammaR.yy - gammaR.xy * gammaR.xy
		);
		<? end ?>

		global real* wave = waveBuf + numWaves * indexInt;
		eigen_calcWaves_<?=side?>_global_global(wave, eig, xInt);
	}<? end ?>
}
		
<? 
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,solver.dim-1 do 
?>

void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) {
	const real ise = sqrt_1_2 / sqrt(eig->eps);
	const real isu = sqrt_1_2 / sqrt(eig->mu);

	<? if side == 0 then ?>
	
	Y[0] = X[2] *  ise + X[4] * isu;
	Y[1] = X[1] * -ise + X[5] * isu;
	Y[2] = X[0] * -ise + X[3] * isu;
	Y[3] = X[0] *  ise + X[3] * isu;
	Y[4] = X[1] *  ise + X[5] * isu;
	Y[5] = X[2] * -ise + X[4] * isu;
	
	<? elseif side == 1 then ?>
	
	Y[0] = X[0] *  ise + X[5] * isu;
	Y[1] = X[2] * -ise + X[3] * isu;
	Y[2] = X[1] * -ise + X[4] * isu;
	Y[3] = X[1] *  ise + X[4] * isu;
	Y[4] = X[2] *  ise + X[3] * isu;
	Y[5] = X[0] * -ise + X[5] * isu;
	
	<? elseif side == 2 then ?>
	
	Y[0] = X[1] *  ise + X[3] * isu;
	Y[1] = X[0] * -ise + X[4] * isu;
	Y[2] = X[2] * -ise + X[5] * isu;
	Y[3] = X[2] *  ise + X[5] * isu;
	Y[4] = X[0] *  ise + X[4] * isu;
	Y[5] = X[1] * -ise + X[3] * isu;
	
	<? end ?>
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) {
	const real se = sqrt_1_2 * sqrt(eig->eps);
	const real su = sqrt_1_2 * sqrt(eig->mu);

	<? if side==0 then ?>
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
	Y[0] = se * (-X[2] + X[3]);
	Y[1] = se * (-X[1] + X[4]);
	Y[2] = se * (X[0] + -X[5]);
	Y[3] = su * (X[2] + X[3]);
	Y[4] = su * (X[0] + X[5]);
	Y[5] = su * (X[1] + X[4]);
	
	<? elseif side==1 then ?>

/*
x, -z, -y, y, z, -x
z,  x,  y, y, x,  z

1  0  0 0 0 -1
0  0 -1 1 0  0
0 -1  0 0 1  0
0  1  0 0 1  0
0  0  1 1 0  0
1  0  0 0 0  1
*/
	Y[0] = se * (X[0] - X[5]);
	Y[1] = se * (-X[2] + X[3]);
	Y[2] = se * (-X[1] + X[4]);
	Y[3] = su * (X[1] + X[4]);
	Y[4] = su * (X[2] + X[3]);
	Y[5] = su * (X[0] + X[5]);
	
	<? elseif side==2 then ?>

/*
y, -x, -z, z, x, -y
x,  y,  z, z,  y,  x

0 -1  0 0 1  0
1  0  0 0 0 -1
0  0 -1 1 0  0
1  0  0 0 0  1
0  1  0 0 1  0
0  0  1 1 0  0
*/
	Y[0] = se * (-X[1] + X[4]);
	Y[1] = se * (X[0] - X[5]);
	Y[2] = se * (-X[2] + X[3]);
	Y[3] = su * (X[0] + X[5]);
	Y[4] = su * (X[1] + X[4]);
	Y[5] = su * (X[2] + X[3]);
	
	<? end ?>
	
	Y[6] = 0;	//BPot
}

<? 
				if solver.checkFluxError then 
?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X_,
	real3 x
) {
	//swap input dim x<->side
	<?=addr2?> const <?=eqn.cons_t?>* X = (<?=addr2?> const <?=eqn.cons_t?>*)X_;
	real3 epsE = X->epsE;
	real3 B = X->B;
	real eps = eig->eps;
	real mu = eig->mu;

	<? if side==0 then ?>
	
	Y[0] = 0;
	Y[1] = B.z / mu;
	Y[2] = -B.y / mu;
	Y[3] = 0;
	Y[4] = -epsE.z / eps;
	Y[5] = epsE.y / eps;

	<? elseif side==1 then ?>
		
	Y[0] = -B.z / mu;
	Y[1] = 0;
	Y[2] = B.x / mu;
	Y[3] = epsE.z / eps;
	Y[4] = 0;
	Y[5] = -epsE.x / eps;
		
	<? elseif side==2 then ?>
		
	Y[0] = B.y / mu;
	Y[1] = -B.x / mu;
	Y[2] = 0;
	Y[3] = -epsE.y / eps;
	Y[4] = epsE.x / eps;
	Y[5] = 0;
		
	<? end ?>
}
<?
				end
			end
		end
	end
end
?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(U->epsE, 1. / U->eps * U->sigma));
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	<?=eqn.eigen_t?>* eig,
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	eig->eps = U->eps;
	eig->mu = U->mu;
}
<? end ?>

void apply_dU_dW(
	<?=eqn.cons_t?>* U, 
	const <?=eqn.prim_t?>* WA, 
	const <?=eqn.prim_t?>* W, 
	real3 x
) { *U = *W; }

void apply_dW_dU(
	<?=eqn.prim_t?>* W,
	const <?=eqn.prim_t?>* WA,
	const <?=eqn.cons_t?>* U,
	real3 x
) { *W = *U; }
