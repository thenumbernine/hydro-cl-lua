#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 B = U.B;
	real3 epsE = U.epsE;
	real mu = U.mu;
	real eps = U.eps;
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.epsE = _real3(0., B.z / mu, -B.y / mu),
		.B = _real3(0., -epsE.z / eps, epsE.y / eps),
	<? elseif side == 1 then ?>
		.epsE = _real3(-B.z / mu, 0., B.x / mu),
		.B = _real3(epsE.z / eps, 0., -epsE.x / eps),
	<? elseif side == 2 then ?>
		.epsE = _real3(B.y / mu, -B.x / mu, 0.),
		.B = _real3(-epsE.y / eps, epsE.x / eps, 0.),
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
	real3 x
) {
	real lambda = 1. / sqrt(U->eps * U->mu);
	return (range_t){-lambda, lambda};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_eps = sqrt(.5 * (UL->eps + UR->eps)),
		.sqrt_mu = sqrt(.5 * (UL->mu + UR->mu)),
	};
}
<? end ?>

//same as in eqn/euler.cl
kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
		
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		*eig = eigen_forSide_<?=side?>(UL, UR, xInt);
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
	const real ise = sqrt_1_2 / eig->sqrt_eps;
	const real isu = sqrt_1_2 / eig->sqrt_mu;

	<? if side==0 then ?>
	
	Y[0] = X[2] *  ise + X[4] * isu;
	Y[1] = X[1] * -ise + X[5] * isu;
	Y[2] = X[0] * -ise + X[3] * isu;
	Y[3] = X[0] *  ise + X[3] * isu;
	Y[4] = X[1] *  ise + X[5] * isu;
	Y[5] = X[2] * -ise + X[4] * isu;
	
	<? elseif side==1 then ?>
	
	Y[0] = X[0] *  ise + X[5] * isu;
	Y[1] = X[2] * -ise + X[3] * isu;
	Y[2] = X[1] * -ise + X[4] * isu;
	Y[3] = X[1] *  ise + X[4] * isu;
	Y[4] = X[2] *  ise + X[3] * isu;
	Y[5] = X[0] * -ise + X[5] * isu;
	
	<? elseif side==2 then ?>
	
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
	const real se = sqrt_1_2 * eig->sqrt_eps;
	const real su = sqrt_1_2 * eig->sqrt_mu;

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
	
	for (int i = numWaves; i < numStates; ++i) {
		Y[i] = 0;
	}
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
	real eps = eig->sqrt_eps * eig->sqrt_eps;
	real mu = eig->sqrt_mu * eig->sqrt_mu;

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
	
	for (int i = numWaves; i < numStates; ++i) {
		Y[i] = 0;
	}
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
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_eps = sqrt(U->eps),
		.sqrt_mu = sqrt(U->mu),
	};
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
