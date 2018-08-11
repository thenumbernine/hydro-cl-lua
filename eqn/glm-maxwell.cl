<?
local solver = eqn.solver
local common = require 'common'()
local xNames = common.xNames
local sym = common.sym
?>

#define divPhiWavespeed 	1.
#define divPsiWavespeed		1.

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=vec3?> E = calc_E(U);
	<?=vec3?> H = calc_H(U);
	<?=scalar?> _1_eps = <?=mul?>(U.sqrt_1_eps, U.sqrt_1_eps);
	<?=scalar?> _1_mu = <?=mul?>(U.sqrt_1_mu, U.sqrt_1_mu);
	<?=scalar?> v_pSq = <?=mul?>(_1_eps, _1_mu);
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.D = _real3(v_pSq * divPhiWavespeed * U.phi,  H.z, -H.y),
		.B = _real3(v_pSq * divPsiWavespeed * U.psi, -E.z,  E.y),
	<? elseif side == 1 then ?>
		.D = _real3(-H.z, v_pSq * divPhiWavespeed * U.phi,  H.x),
		.B = _real3( E.z, v_pSq * divPsiWavespeed * U.psi, -E.x),
	<? elseif side == 2 then ?>
		.D = _real3( H.y, -H.x, v_pSq * divPhiWavespeed * U.phi),
		.B = _real3(-E.y,  E.x, v_pSq * divPsiWavespeed * U.psi),
	<? end ?>
		.phi = divPhiWavespeed * U.D.s<?=side?>,
		.psi = divPsiWavespeed * U.B.s<?=side?>,
	
		.sigma = <?=zero?>,
		.rhoCharge = <?=zero?>,
		.sqrt_1_eps = <?=zero?>,
		.sqrt_1_mu = <?=zero?>,
	};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
<? if scalar == 'real' then ?>
	real v_p = U->sqrt_1_eps * U->sqrt_1_mu;
<? else ?>
	real v_p = <?=abs?>(<?=mul?>(U->sqrt_1_eps, U->sqrt_1_mu));
<? end ?>
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
		.sqrt_1_eps = <?=real_mul?>(<?=add?>(UL.sqrt_1_eps, UR.sqrt_1_eps), .5),
		.sqrt_1_mu = <?=real_mul?>(<?=add?>(UL.sqrt_1_mu, UR.sqrt_1_mu), .5),
	};
}

<? for side=0,solver.dim-1 do ?>
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=eqn.waves_t?> Y;
	<?=scalar?>* Yp = (<?=scalar?>*)Y.ptr;

	<?=scalar?> sqrt_1_eps = eig.sqrt_1_eps;
	<?=scalar?> sqrt_1_mu = eig.sqrt_1_mu;
	<?=scalar?> sqrt_eps = <?=inv?>(sqrt_1_eps);
	<?=scalar?> sqrt_mu = <?=inv?>(sqrt_1_mu);

	<? if side==0 then ?>

	Yp[0] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[0], sqrt_mu ), <?=mul?>(X.ptr[6], sqrt_1_eps)));
	Yp[1] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[3], sqrt_mu ), <?=mul?>(X.ptr[7], sqrt_1_eps)));
	Yp[2] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[1], sqrt_mu ), <?=mul?>(X.ptr[5], sqrt_eps)));
	Yp[3] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[4], sqrt_eps), <?=mul?>(X.ptr[2], sqrt_mu)));
	Yp[4] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[5], sqrt_eps), <?=mul?>(X.ptr[1], sqrt_mu)));
	Yp[5] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[2], sqrt_mu ), <?=mul?>(X.ptr[4], sqrt_eps)));
	Yp[6] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[0], sqrt_mu ), <?=mul?>(X.ptr[6], sqrt_1_eps)));
	Yp[7] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[3], sqrt_mu ), <?=mul?>(X.ptr[7], sqrt_1_eps)));
	
	<? elseif side==1 then ?>
	
	Yp[0] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[1], sqrt_mu ), <?=mul?>(X.ptr[6], sqrt_1_eps)));
	Yp[1] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[4], sqrt_mu ), <?=mul?>(X.ptr[7], sqrt_1_eps)));
	Yp[2] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[5], sqrt_eps), <?=mul?>(X.ptr[0], sqrt_mu)));
	Yp[3] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[2], sqrt_mu ), <?=mul?>(X.ptr[3], sqrt_eps)));
	Yp[4] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[0], sqrt_mu ), <?=mul?>(X.ptr[5], sqrt_eps)));
	Yp[5] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[3], sqrt_eps), <?=mul?>(X.ptr[2], sqrt_mu)));
	Yp[6] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[1], sqrt_mu ), <?=mul?>(X.ptr[6], sqrt_1_eps)));
	Yp[7] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[4], sqrt_mu ), <?=mul?>(X.ptr[7], sqrt_1_eps)));
	
	<? elseif side==2 then ?>
	
	Yp[0] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[2], sqrt_mu ), <?=mul?>(X.ptr[6], sqrt_1_eps)));
	Yp[1] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[5], sqrt_mu ), <?=mul?>(X.ptr[7], sqrt_1_eps)));
	Yp[2] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[0], sqrt_mu ), <?=mul?>(X.ptr[4], sqrt_eps)));
	Yp[3] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[3], sqrt_eps), <?=mul?>(X.ptr[1], sqrt_mu)));
	Yp[4] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[4], sqrt_eps), <?=mul?>(X.ptr[0], sqrt_mu)));
	Yp[5] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[1], sqrt_mu ), <?=mul?>(X.ptr[3], sqrt_eps)));
	Yp[6] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[2], sqrt_mu ), <?=mul?>(X.ptr[6], sqrt_1_eps)));
	Yp[7] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[5], sqrt_mu ), <?=mul?>(X.ptr[7], sqrt_1_eps)));

	<? end ?>

	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;
	<?=scalar?>* Yp = (<?=scalar?>*)Y.ptr;

	<?=scalar?> sqrt_1_eps = eig.sqrt_1_eps;
	<?=scalar?> sqrt_1_mu = eig.sqrt_1_mu;
	<?=scalar?> sqrt_eps = <?=inv?>(eig.sqrt_1_eps);
	<?=scalar?> sqrt_mu = <?=inv?>(eig.sqrt_1_mu);

	<? if side==0 then ?>

	Yp[0] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[6], sqrt_1_mu ), <?=mul?>(X.ptr[0], sqrt_1_mu)));
	Yp[1] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[4], sqrt_1_mu ), <?=mul?>(X.ptr[2], sqrt_1_mu)));
	Yp[2] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[5], sqrt_1_mu ), <?=mul?>(X.ptr[3], sqrt_1_mu)));
	Yp[3] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[7], sqrt_1_mu ), <?=mul?>(X.ptr[1], sqrt_1_mu)));
	Yp[4] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[3], sqrt_1_eps), <?=mul?>(X.ptr[5], sqrt_1_eps)));
	Yp[5] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[4], sqrt_1_eps), <?=mul?>(X.ptr[2], sqrt_1_eps)));
	Yp[6] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[6], sqrt_eps  ), <?=mul?>(X.ptr[0], sqrt_eps)));
	Yp[7] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[7], sqrt_eps  ), <?=mul?>(X.ptr[1], sqrt_eps)));
	
	<? elseif side==1 then ?>
	
	Yp[0] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[4], sqrt_1_mu ), <?=mul?>(X.ptr[2], sqrt_1_mu)));
	Yp[1] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[6], sqrt_1_mu ), <?=mul?>(X.ptr[0], sqrt_1_mu)));
	Yp[2] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[5], sqrt_1_mu ), <?=mul?>(X.ptr[3], sqrt_1_mu)));
	Yp[3] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[5], sqrt_1_eps), <?=mul?>(X.ptr[3], sqrt_1_eps)));
	Yp[4] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[7], sqrt_1_mu ), <?=mul?>(X.ptr[1], sqrt_1_mu)));
	Yp[5] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[2], sqrt_1_eps), <?=mul?>(X.ptr[4], sqrt_1_eps)));
	Yp[6] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[6], sqrt_eps  ), <?=mul?>(X.ptr[0], sqrt_eps)));
	Yp[7] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[7], sqrt_eps  ), <?=mul?>(X.ptr[1], sqrt_eps)));
	
	<? elseif side==2 then ?>
	
	Yp[0] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[4], sqrt_1_mu ), <?=mul?>(X.ptr[2], sqrt_1_mu)));
	Yp[1] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[5], sqrt_1_mu ), <?=mul?>(X.ptr[3], sqrt_1_mu)));
	Yp[2] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[6], sqrt_1_mu ), <?=mul?>(X.ptr[0], sqrt_1_mu)));
	Yp[3] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[3], sqrt_1_eps), <?=mul?>(X.ptr[5], sqrt_1_eps)));
	Yp[4] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[4], sqrt_1_eps), <?=mul?>(X.ptr[2], sqrt_1_eps)));
	Yp[5] = <?=mul?>(sqrt_1_2, <?=add?>(<?=mul?>(X.ptr[7], sqrt_1_mu ), <?=mul?>(X.ptr[1], sqrt_1_mu)));
	Yp[6] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[6], sqrt_eps  ), <?=mul?>(X.ptr[0], sqrt_eps)));
	Yp[7] = <?=mul?>(sqrt_1_2, <?=sub?>(<?=mul?>(X.ptr[7], sqrt_eps  ), <?=mul?>(X.ptr[1], sqrt_eps)));

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

	//TODO J = J_f + J_b = J_f + J_P + J_M = J_f + dP/dt + curl M
	deriv->D = <?=vec3?>_sub(
		deriv->D, 
		<?=vec3?>_<?=scalar?>_mul(
			U->D, 
			<?=mul?>(<?=mul?>(U->sqrt_1_eps, U->sqrt_1_eps), U->sigma)
		)
	);


	//for non-time-varying susceptibilities, here's the source term:
	//D_i,t ... = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_,j B_k - J_i
	//B_i,t ... = 1/sqrt(g) g_il epsBar^ljk (1/eps)_,j B_k

	<?=vec3?> grad_1_mu = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[stepsize.<?=xj?>].sqrt_1_mu, U[stepsize.<?=xj?>].sqrt_1_mu),
			<?=mul?>(U[-stepsize.<?=xj?>].sqrt_1_mu, U[-stepsize.<?=xj?>].sqrt_1_mu)
		), 1. / grid_dx<?=j?>);
	<? end ?>
	
	<?=vec3?> grad_1_eps = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[stepsize.<?=xj?>].sqrt_1_eps, U[stepsize.<?=xj?>].sqrt_1_eps),
			<?=mul?>(U[-stepsize.<?=xj?>].sqrt_1_eps, U[-stepsize.<?=xj?>].sqrt_1_eps)
		), 1. / grid_dx<?=j?>);
	<? end ?>

	real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>{
		<?=eqn.cons_t?> flux = fluxFromCons_<?=j?>(*U, x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = <?=sub?>(deriv->D.<?=xj?>, <?=vec3?>_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = <?=sub?>(deriv->B.<?=xj?>, <?=vec3?>_dot(flux.B, grad_1_eps));
	}<? end ?>
	
	deriv->phi = <?=add?>(deriv->phi, <?=real_mul?>(U->rhoCharge, divPhiWavespeed));
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
