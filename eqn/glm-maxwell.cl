#define chi 			1.
#define gamma			1.
#define sqrt_mu0 		1.
#define sqrt_epsilon0 	1.
#define mu_0			(sqrt_mu0 * sqrt_mu0)
#define epsilon_0		(sqrt_epsilon0 * sqrt_epsilon0)
#define speedOfLight	(1./(sqrt_mu0 * sqrt_epsilon0))
#define speedOfLightSq (speedOfLight * speedOfLight)

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 B = U.B;
	real3 E = U.E;
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.E = _real3(speedOfLightSq * chi * U.phi, speedOfLightSq * B.z, -speedOfLightSq * B.y),
		.B = _real3(gamma * U.psi, -E.z, E.y),
	<? elseif side == 1 then ?>
		.E = _real3(-speedOfLightSq * B.z, speedOfLightSq * chi * U.phi, speedOfLightSq * B.x),
		.B = _real3(E.z, gamma * U.psi, -E.x),
	<? elseif side == 2 then ?>
		.E = _real3(speedOfLightSq * B.y, -speedOfLightSq * B.x, speedOfLightSq * chi * U.phi),
		.B = _real3(-E.y, E.x, gamma * U.psi),
	<? end ?>
		.phi = chi * E.s<?=side?>,
		.psi = speedOfLightSq * gamma * B.s<?=side?>,
	
		.conductivity = 0.,
		.charge = 0.,
	};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	real lambda = max(max(gamma, chi), 1.) * speedOfLight;
	return (range_t){-lambda, lambda};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	return (<?=eqn.eigen_t?>){};
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
	<? if side==0 then ?>

	Y[0] = (X[0] - X[6] * speedOfLight) * .5;
	Y[1] = (X[3] - X[7] / speedOfLight) * .5;
	Y[2] = (X[1] - X[5] * speedOfLight) * .5;
	Y[3] = (X[4] * speedOfLight + X[2]) * .5;
	Y[4] = (X[5] * speedOfLight + X[1]) * .5;
	Y[5] = (X[2] - X[4] * speedOfLight) * .5;
	Y[6] = (X[7] / speedOfLight + X[3]) * .5;
	Y[7] = (X[6] * speedOfLight + X[0]) * .5;
   
   <? elseif side==1 then ?>
   
	Y[0] = (X[1] - X[6] * speedOfLight) * .5;
	Y[1] = (X[4] - X[7] / speedOfLight) * .5;
	Y[2] = (X[5] * speedOfLight + X[0]) * .5;
	Y[3] = (X[2] - X[3] * speedOfLight) * .5;
	Y[4] = (X[0] - X[5] * speedOfLight) * .5;
	Y[5] = (X[3] * speedOfLight + X[2]) * .5;
	Y[6] = (X[7] / speedOfLight + X[4]) * .5;
	Y[7] = (X[6] * speedOfLight + X[1]) * .5;
   
   <? elseif side==2 then ?>
   
	Y[0] = (X[2] - X[6] * speedOfLight) * .5;
	Y[1] = (X[5] - X[7] / speedOfLight) * .5;
	Y[2] = (X[0] - X[4] * speedOfLight) * .5;
	Y[3] = (X[3] * speedOfLight + X[1]) * .5;
	Y[4] = (X[4] * speedOfLight + X[0]) * .5;
	Y[5] = (X[1] - X[3] * speedOfLight) * .5;
	Y[6] = (X[7] / speedOfLight + X[5]) * .5;
	Y[7] = (X[6] * speedOfLight + X[2]) * .5;
	
	<? end ?>
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) {

	<? if side==0 then ?>
	
	Y[0] = X[7] + X[0];
	Y[1] = X[4] + X[2];
	Y[2] = X[5] + X[3];
	Y[3] = X[6] + X[1];
	Y[4] = (X[3] - X[5]) / speedOfLight;
	Y[5] = (X[4] - X[2]) / speedOfLight;
	Y[6] = (X[7] - X[0]) / speedOfLight;
	Y[7] = (X[6] - X[1]) * speedOfLight;
   
   <? elseif side==1 then ?>

	Y[0] = X[4] + X[2];
	Y[1] = X[7] + X[0];
	Y[2] = X[5] + X[3];
	Y[3] = (X[5] - X[3]) / speedOfLight;
	Y[4] = X[6] + X[1];
	Y[5] = (X[2] - X[4]) / speedOfLight;
	Y[6] = (X[7] - X[0]) / speedOfLight;
	Y[7] = (X[6] - X[1]) * speedOfLight;
   
   <? elseif side==2 then ?>

	Y[0] = X[4] + X[2];
	Y[1] = X[5] + X[3];
	Y[2] = X[7] + X[0];
	Y[3] = (X[3] - X[5]) / speedOfLight;
	Y[4] = (X[4] - X[2]) / speedOfLight;
	Y[5] = X[6] + X[1];
	Y[6] = (X[7] - X[0]) / speedOfLight;
	Y[7] = (X[6] - X[1]) * speedOfLight;
	
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
	<?=addr2?> const <?=eqn.cons_t?>* X = (<?=addr2?> const <?=eqn.cons_t?>*)X_;
	*(<?=addr0?> <?=eqn.cons_t?>*)Y = fluxFromCons_<?=side?>(*X, x);
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
	real3 J = real3_scale(U->E, mu_0 / U->conductivity);
	deriv->E = real3_sub(deriv->E, J);
	deriv->phi += U->charge * chi / epsilon_0;
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){};
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
