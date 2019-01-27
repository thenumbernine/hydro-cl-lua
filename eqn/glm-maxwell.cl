<?
local solver = eqn.solver
local common = require 'common'()
local xNames = common.xNames
local sym = common.sym
?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
cons_t fluxFromCons_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	<?=vec3?> E = calc_E(U);
	<?=vec3?> H = calc_H(U);
	<?=scalar?> _1_eps = <?=mul?>(U.sqrt_1_eps, U.sqrt_1_eps);
	<?=scalar?> _1_mu = <?=mul?>(U.sqrt_1_mu, U.sqrt_1_mu);
	return (cons_t){
	<? if side == 0 then ?>
		.D = _<?=vec3?>(<?=real_mul?>(U.phi, solver->divPhiWavespeed),  H.z, <?=neg?>(H.y)),
		.B = _<?=vec3?>(<?=real_mul?>(U.psi, solver->divPsiWavespeed), <?=neg?>(E.z),  E.y),
	<? elseif side == 1 then ?>
		.D = _<?=vec3?>(<?=neg?>(H.z), <?=real_mul?>(U.phi, solver->divPhiWavespeed),  H.x),
		.B = _<?=vec3?>( E.z,          <?=real_mul?>(U.psi, solver->divPsiWavespeed), <?=neg?>(E.x)),
	<? elseif side == 2 then ?>
		.D = _<?=vec3?>( H.y, <?=neg?>(H.x), <?=real_mul?>(U.phi, solver->divPhiWavespeed)),
		.B = _<?=vec3?>(<?=neg?>(E.y),  E.x, <?=real_mul?>(U.psi, solver->divPsiWavespeed)),
	<? end ?>
		.phi = <?=real_mul?>(U.D.s<?=side?>, solver->divPhiWavespeed),
		.psi = <?=real_mul?>(U.B.s<?=side?>, solver->divPsiWavespeed),
	
		.sigma = <?=zero?>,
		.rhoCharge = <?=zero?>,
		.sqrt_1_eps = <?=zero?>,
		.sqrt_1_mu = <?=zero?>,
	};
}
<? end ?>

eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3 n
) {
	return (eigen_t){
		.sqrt_1_eps = <?=real_mul?>(<?=add?>(UL.sqrt_1_eps, UR.sqrt_1_eps), .5),
		.sqrt_1_mu = <?=real_mul?>(<?=add?>(UL.sqrt_1_mu, UR.sqrt_1_mu), .5),
	};
}

<? for side=0,solver.dim-1 do ?>
waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	waves_t Y;
	<?=scalar?>* Xp = (<?=scalar?>*)X.ptr;
	<?=scalar?>* Yp = (<?=scalar?>*)Y.ptr;

	<?=scalar?> sqrt_1_eps = eig.sqrt_1_eps;
	<?=scalar?> sqrt_1_mu = eig.sqrt_1_mu;
	<?=scalar?> sqrt_eps = <?=inv?>(sqrt_1_eps);
	<?=scalar?> sqrt_mu = <?=inv?>(sqrt_1_mu);
	<?=scalar?> sqrt_2 = <?=inv?>(sqrt_1_2);

	<? if side==0 then ?>

	Yp[0] = ((-(sqrt_eps * (Xp[0] - Xp[6]))) / sqrt_2);
	Yp[1] = ((-(sqrt_eps * (Xp[3] - Xp[7]))) / sqrt_2);
	Yp[2] = (((Xp[2] * sqrt_mu) + (Xp[4] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));
	Yp[3] = (((Xp[1] * sqrt_mu) - (Xp[5] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
	Yp[4] = ((-((Xp[2] * sqrt_mu) - (Xp[4] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));
	Yp[5] = (((Xp[1] * sqrt_mu) + (Xp[5] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));
	Yp[6] = ((sqrt_eps * (Xp[0] + Xp[6])) / sqrt_2);
	Yp[7] = ((sqrt_eps * (Xp[3] + Xp[7])) / sqrt_2);

	<? elseif side==1 then ?>

	Yp[0] = ((-(sqrt_eps * (Xp[1] - Xp[6]))) / sqrt_2);
	Yp[1] = ((-(sqrt_eps * (Xp[4] - Xp[7]))) / sqrt_2);
	Yp[2] = (((Xp[2] * sqrt_mu) - (Xp[3] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
	Yp[3] = (((Xp[0] * sqrt_mu) + (Xp[5] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));
	Yp[4] = (((Xp[2] * sqrt_mu) + (Xp[3] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));
	Yp[5] = ((-((Xp[0] * sqrt_mu) - (Xp[5] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));
	Yp[6] = ((sqrt_eps * (Xp[1] + Xp[6])) / sqrt_2);
	Yp[7] = ((sqrt_eps * (Xp[4] + Xp[7])) / sqrt_2);

	<? elseif side==2 then ?>
	
	Yp[0] = ((-(sqrt_eps * (Xp[2] - Xp[6]))) / sqrt_2);
	Yp[1] = ((-(sqrt_eps * (Xp[5] - Xp[7]))) / sqrt_2);
	Yp[2] = (((Xp[1] * sqrt_mu) + (Xp[3] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));
	Yp[3] = (((Xp[0] * sqrt_mu) - (Xp[4] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
	Yp[4] = (((Xp[1] * sqrt_mu) - (Xp[3] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
	Yp[5] = (((Xp[0] * sqrt_mu) + (Xp[4] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));
	Yp[6] = ((sqrt_eps * (Xp[2] + Xp[6])) / sqrt_2);
	Yp[7] = ((sqrt_eps * (Xp[5] + Xp[7])) / sqrt_2);

	<? end ?>

	return Y;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x
) {
	cons_t Y;
	<?=scalar?>* Xp = (<?=scalar?>*)X.ptr;
	<?=scalar?>* Yp = (<?=scalar?>*)Y.ptr;

	<?=scalar?> sqrt_1_eps = eig.sqrt_1_eps;
	<?=scalar?> sqrt_1_mu = eig.sqrt_1_mu;
	<?=scalar?> sqrt_eps = <?=inv?>(eig.sqrt_1_eps);
	<?=scalar?> sqrt_mu = <?=inv?>(eig.sqrt_1_mu);
	<?=scalar?> sqrt_2 = <?=inv?>(sqrt_1_2);

	<? if side==0 then ?>

	Yp[0] = ((-(Xp[0] - Xp[6])) / (sqrt_2 * sqrt_eps));
	Yp[1] = ((-(sqrt_eps * (Xp[3] - Xp[5]))) / sqrt_2);
	Yp[2] = ((sqrt_eps * (Xp[2] - Xp[4])) / sqrt_2);
	Yp[3] = ((-(Xp[1] - Xp[7])) / (sqrt_2 * sqrt_eps));
	Yp[4] = ((sqrt_mu * (Xp[2] + Xp[4])) / sqrt_2);
	Yp[5] = ((sqrt_mu * (Xp[3] + Xp[5])) / sqrt_2);
	Yp[6] = ((Xp[0] + Xp[6]) / (sqrt_2 * sqrt_eps));
	Yp[7] = ((Xp[1] + Xp[7]) / (sqrt_2 * sqrt_eps));
	
	<? elseif side==1 then ?>

	Yp[0] = ((sqrt_eps * (Xp[3] - Xp[5])) / sqrt_2);
	Yp[1] = ((-(Xp[0] - Xp[6])) / (sqrt_2 * sqrt_eps));
	Yp[2] = ((-(sqrt_eps * (Xp[2] - Xp[4]))) / sqrt_2);
	Yp[3] = ((sqrt_mu * (Xp[2] + Xp[4])) / sqrt_2);
	Yp[4] = ((-(Xp[1] - Xp[7])) / (sqrt_2 * sqrt_eps));
	Yp[5] = ((sqrt_mu * (Xp[3] + Xp[5])) / sqrt_2);
	Yp[6] = ((Xp[0] + Xp[6]) / (sqrt_2 * sqrt_eps));
	Yp[7] = ((Xp[1] + Xp[7]) / (sqrt_2 * sqrt_eps));
	
	<? elseif side==2 then ?>

	Yp[0] = ((-(sqrt_eps * (Xp[3] - Xp[5]))) / sqrt_2);
	Yp[1] = ((sqrt_eps * (Xp[2] - Xp[4])) / sqrt_2);
	Yp[2] = ((-(Xp[0] - Xp[6])) / (sqrt_2 * sqrt_eps));
	Yp[3] = ((sqrt_mu * (Xp[2] + Xp[4])) / sqrt_2);
	Yp[4] = ((sqrt_mu * (Xp[3] + Xp[5])) / sqrt_2);
	Yp[5] = ((-(Xp[1] - Xp[7])) / (sqrt_2 * sqrt_eps));
	Yp[6] = ((Xp[0] + Xp[6]) / (sqrt_2 * sqrt_eps));
	Yp[7] = ((Xp[1] + Xp[7]) / (sqrt_2 * sqrt_eps));

	<? end ?>

	for (int i = <?=eqn.numWaves?>; i < <?=eqn.numStates?>; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	return fluxFromCons_<?=side?>(solver, X, x);
}

<? end ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);

	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

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
			<?=mul?>(U[solver->stepsize.<?=xj?>].sqrt_1_mu, U[solver->stepsize.<?=xj?>].sqrt_1_mu),
			<?=mul?>(U[-solver->stepsize.<?=xj?>].sqrt_1_mu, U[-solver->stepsize.<?=xj?>].sqrt_1_mu)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>
	
	<?=vec3?> grad_1_eps = <?=vec3?>_zero;
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = <?=real_mul?>(
		<?=sub?>(
			<?=mul?>(U[solver->stepsize.<?=xj?>].sqrt_1_eps, U[solver->stepsize.<?=xj?>].sqrt_1_eps),
			<?=mul?>(U[-solver->stepsize.<?=xj?>].sqrt_1_eps, U[-solver->stepsize.<?=xj?>].sqrt_1_eps)
		), 1. / solver->grid_dx.s<?=j?>);
	<? end ?>

	real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>{
		cons_t flux = fluxFromCons_<?=j?>(solver, *U, x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> = <?=sub?>(deriv->D.<?=xj?>, <?=vec3?>_dot(flux.D, grad_1_mu));
		deriv->B.<?=xj?> = <?=sub?>(deriv->B.<?=xj?>, <?=vec3?>_dot(flux.B, grad_1_eps));
	}<? end ?>
	
	deriv->phi = <?=add?>(deriv->phi, <?=real_mul?>(U->rhoCharge, solver->divPhiWavespeed));
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
eigen_t eigen_forCell_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	return (eigen_t){
		.sqrt_1_eps = U.sqrt_1_eps,
		.sqrt_1_mu = U.sqrt_1_mu,
	};
}
<? end ?>
