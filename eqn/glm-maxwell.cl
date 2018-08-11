<?
local solver = eqn.solver
local common = require 'common'()
local xNames = common.xNames
local sym = common.sym
?>

#define divPhiWavespeed 	10.
#define divPsiWavespeed		10.

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>
<? local sqrt2 = math.sqrt(.5) ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 E = calc_E(U);
	real3 H = calc_H(U);
	real v_p = U.sqrt_1_eps * U.sqrt_1_eps * U.sqrt_1_mu * U.sqrt_1_mu;
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.D = _real3(v_p * divPhiWavespeed * U.phi,  H.z, -H.y),
		.B = _real3(v_p * divPsiWavespeed * U.psi, -E.z,  E.y),
	<? elseif side == 1 then ?>
		.D = _real3(-H.z, v_p * divPhiWavespeed * U.phi,  H.x),
		.B = _real3( E.z, v_p * divPsiWavespeed * U.psi, -E.x),
	<? elseif side == 2 then ?>
		.D = _real3( H.y, -H.x, v_p * divPhiWavespeed * U.phi),
		.B = _real3(-E.y,  E.x, v_p * divPsiWavespeed * U.psi),
	<? end ?>
		.phi = divPhiWavespeed * U.D.s<?=side?>,
		.psi = divPsiWavespeed * U.B.s<?=side?>,
	
		.sigma = 0.,
		.rhoCharge = 0.,
		.sqrt_1_eps = 0.,
		.sqrt_1_mu = 0.,
	};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	real v_p = U->sqrt_1_eps * U->sqrt_1_eps * U->sqrt_1_mu * U->sqrt_1_mu;
	real lambda = max(max(divPsiWavespeed, divPhiWavespeed), 1.) * v_p;
	return (range_t){-lambda, lambda};
}
<? end ?>

<?=eqn.eigen_t?> eigen_forInterface(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	real3 n
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_1_eps = .5 * (UL.sqrt_1_eps + UR.sqrt_1_eps),
		.sqrt_1_mu = .5 * (UL.sqrt_1_mu + UR.sqrt_1_mu),
	};
}

<? for side=0,solver.dim-1 do ?>
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=eqn.waves_t?> Y;
	real* Yp = Y.ptr;

	real sqrt_1_eps = eig.sqrt_1_eps;
	real sqrt_1_mu = eig.sqrt_1_mu;
	real sqrt_eps = 1. / sqrt_1_eps;
	real sqrt_mu = 1. / sqrt_1_mu;
	const real sqrt2 = <?=sqrt2?>;

	<? if side==0 then ?>

	Yp[0] = X.ptr[0] * sqrt_1_2 * sqrt_mu - X.ptr[6] * sqrt_1_2 * sqrt_1_eps;
	Yp[1] = X.ptr[3] * sqrt_1_2 * sqrt_mu - X.ptr[7] * sqrt_1_2 * sqrt_1_eps;
	Yp[2] = X.ptr[1] * sqrt_1_2 * sqrt_mu - X.ptr[5] * sqrt_eps * sqrt_1_2;
	Yp[3] = X.ptr[4] * sqrt_eps * sqrt_1_2 + X.ptr[2] * sqrt_1_2 * sqrt_mu;
	Yp[4] = X.ptr[5] * sqrt_eps * sqrt_1_2 + X.ptr[1] * sqrt_1_2 * sqrt_mu;
	Yp[5] = X.ptr[2] * sqrt_1_2 * sqrt_mu - X.ptr[4] * sqrt_eps * sqrt_1_2;
	Yp[6] = X.ptr[0] * sqrt_1_2 * sqrt_mu + X.ptr[6] * sqrt_1_2 * sqrt_1_eps;
	Yp[7] = X.ptr[3] * sqrt_1_2 * sqrt_mu + X.ptr[7] * sqrt_1_2 * sqrt_1_eps;
	
	<? elseif side==1 then ?>
	
	Yp[0] = X.ptr[1] * sqrt_1_2 * sqrt_mu - X.ptr[6] * sqrt_1_2 * sqrt_1_eps;
	Yp[1] = X.ptr[4] * sqrt_1_2 * sqrt_mu - X.ptr[7] * sqrt_1_2 * sqrt_1_eps;
	Yp[2] = X.ptr[5] * sqrt_eps * sqrt_1_2 + X.ptr[0] * sqrt_1_2 * sqrt_mu;
	Yp[3] = X.ptr[2] * sqrt_1_2 * sqrt_mu - X.ptr[3] * sqrt_eps * sqrt_1_2;
	Yp[4] = X.ptr[0] * sqrt_1_2 * sqrt_mu - X.ptr[5] * sqrt_eps * sqrt_1_2;
	Yp[5] = X.ptr[3] * sqrt_eps * sqrt_1_2 + X.ptr[2] * sqrt_1_2 * sqrt_mu;
	Yp[6] = X.ptr[1] * sqrt_1_2 * sqrt_mu + X.ptr[6] * sqrt_1_2 * sqrt_1_eps;
	Yp[7] = X.ptr[4] * sqrt_1_2 * sqrt_mu + X.ptr[7] * sqrt_1_2 * sqrt_1_eps;
	
	<? elseif side==2 then ?>
	
	Yp[0] = X.ptr[2] * sqrt_1_2 * sqrt_mu - X.ptr[6] * sqrt_1_2 * sqrt_1_eps;
	Yp[1] = X.ptr[5] * sqrt_1_2 * sqrt_mu - X.ptr[7] * sqrt_1_2 * sqrt_1_eps;
	Yp[2] = X.ptr[0] * sqrt_1_2 * sqrt_mu - X.ptr[4] * sqrt_eps * sqrt_1_2;
	Yp[3] = X.ptr[3] * sqrt_eps * sqrt_1_2 + X.ptr[1] * sqrt_1_2 * sqrt_mu;
	Yp[4] = X.ptr[4] * sqrt_eps * sqrt_1_2 + X.ptr[0] * sqrt_1_2 * sqrt_mu;
	Yp[5] = X.ptr[1] * sqrt_1_2 * sqrt_mu - X.ptr[3] * sqrt_eps * sqrt_1_2;
	Yp[6] = X.ptr[2] * sqrt_1_2 * sqrt_mu + X.ptr[6] * sqrt_1_2 * sqrt_1_eps;
	Yp[7] = X.ptr[5] * sqrt_1_2 * sqrt_mu + X.ptr[7] * sqrt_1_2 * sqrt_1_eps;

	<? end ?>

	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;
	real* Yp = Y.ptr;

	real sqrt_1_eps = eig.sqrt_1_eps;
	real sqrt_1_mu = eig.sqrt_1_mu;
	real sqrt_eps = 1. / eig.sqrt_1_eps;
	real sqrt_mu = 1. / eig.sqrt_1_mu;
	const real sqrt2 = <?=sqrt2?>;

	<? if side==0 then ?>

	Yp[0] = X.ptr[6] * sqrt_1_mu * sqrt_1_2 + X.ptr[0] * sqrt_1_mu * sqrt_1_2;
	Yp[1] = X.ptr[4] * sqrt_1_mu * sqrt_1_2 + X.ptr[2] * sqrt_1_mu * sqrt_1_2;
	Yp[2] = X.ptr[5] * sqrt_1_mu * sqrt_1_2 + X.ptr[3] * sqrt_1_mu * sqrt_1_2;
	Yp[3] = X.ptr[7] * sqrt_1_mu * sqrt_1_2 + X.ptr[1] * sqrt_1_mu * sqrt_1_2;
	Yp[4] = X.ptr[3] * sqrt_1_2 * sqrt_1_eps - X.ptr[5] * sqrt_1_2 * sqrt_1_eps;
	Yp[5] = X.ptr[4] * sqrt_1_2 * sqrt_1_eps - X.ptr[2] * sqrt_1_2 * sqrt_1_eps;
	Yp[6] = X.ptr[6] * sqrt_eps * sqrt_1_2 - X.ptr[0] * sqrt_eps * sqrt_1_2;
	Yp[7] = X.ptr[7] * sqrt_eps * sqrt_1_2 - X.ptr[1] * sqrt_eps * sqrt_1_2;
	
	<? elseif side==1 then ?>
	
	Yp[0] = X.ptr[4] * sqrt_1_mu * sqrt_1_2 + X.ptr[2] * sqrt_1_mu * sqrt_1_2;
	Yp[1] = X.ptr[6] * sqrt_1_mu * sqrt_1_2 + X.ptr[0] * sqrt_1_mu * sqrt_1_2;
	Yp[2] = X.ptr[5] * sqrt_1_mu * sqrt_1_2 + X.ptr[3] * sqrt_1_mu * sqrt_1_2;
	Yp[3] = X.ptr[5] * sqrt_1_2 * sqrt_1_eps - X.ptr[3] * sqrt_1_2 * sqrt_1_eps;
	Yp[4] = X.ptr[7] * sqrt_1_mu * sqrt_1_2 + X.ptr[1] * sqrt_1_mu * sqrt_1_2;
	Yp[5] = X.ptr[2] * sqrt_1_2 * sqrt_1_eps - X.ptr[4] * sqrt_1_2 * sqrt_1_eps;
	Yp[6] = X.ptr[6] * sqrt_eps * sqrt_1_2 - X.ptr[0] * sqrt_eps * sqrt_1_2;
	Yp[7] = X.ptr[7] * sqrt_eps * sqrt_1_2 - X.ptr[1] * sqrt_eps * sqrt_1_2;
	
	<? elseif side==2 then ?>
	
	Yp[0] = X.ptr[4] * sqrt_1_mu * sqrt_1_2 + X.ptr[2] * sqrt_1_mu * sqrt_1_2;
	Yp[1] = X.ptr[5] * sqrt_1_mu * sqrt_1_2 + X.ptr[3] * sqrt_1_mu * sqrt_1_2;
	Yp[2] = X.ptr[6] * sqrt_1_mu * sqrt_1_2 + X.ptr[0] * sqrt_1_mu * sqrt_1_2;
	Yp[3] = X.ptr[3] * sqrt_1_2 * sqrt_1_eps - X.ptr[5] * sqrt_1_2 * sqrt_1_eps;
	Yp[4] = X.ptr[4] * sqrt_1_2 * sqrt_1_eps - X.ptr[2] * sqrt_1_2 * sqrt_1_eps;
	Yp[5] = X.ptr[7] * sqrt_1_mu * sqrt_1_2 + X.ptr[1] * sqrt_1_mu * sqrt_1_2;
	Yp[6] = X.ptr[6] * sqrt_eps * sqrt_1_2 - X.ptr[0] * sqrt_eps * sqrt_1_2;
	Yp[7] = X.ptr[7] * sqrt_eps * sqrt_1_2 - X.ptr[1] * sqrt_eps * sqrt_1_2;

	<? end ?>

	for (int i = 8; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	return fluxFromCons_<?=side?>(X, x);
}

<? end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);

	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

	real3 J = real3_real_mul(U->D, U->sqrt_1_eps * U->sqrt_1_eps * U->sigma);
	deriv->D = real3_sub(deriv->D, J);

	real3 grad_1_mu = real3_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = real_mul(
		real_sub(
			U[stepsize.<?=xj?>].sqrt_1_mu * U[stepsize.<?=xj?>].sqrt_1_mu,
			U[-stepsize.<?=xj?>].sqrt_1_mu * U[-stepsize.<?=xj?>].sqrt_1_mu
		), 1. / grid_dx<?=j?>);
	<? end ?>
	
	real3 grad_1_eps = real3_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = real_mul(
		real_sub(
			U[stepsize.<?=xj?>].sqrt_1_eps * U[stepsize.<?=xj?>].sqrt_1_eps,
			U[-stepsize.<?=xj?>].sqrt_1_eps * U[-stepsize.<?=xj?>].sqrt_1_eps
		), 1. / grid_dx<?=j?>);
	<? end ?>
	
	real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>{
		<?=eqn.cons_t?> flux = fluxFromCons_<?=j?>(*U, x);
		flux.D = real3_real_mul(coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = real3_real_mul(coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = real_sub(deriv->D.<?=xj?>, real3_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = real_sub(deriv->B.<?=xj?>, real3_dot(flux.B, grad_1_eps));
	}<? end ?>
	
	deriv->phi += U->rhoCharge * divPhiWavespeed;
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_1_eps = U.sqrt_1_eps,
		.sqrt_1_mu = U.sqrt_1_mu,
	};
}
<? end ?>
