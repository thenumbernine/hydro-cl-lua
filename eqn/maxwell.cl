<?
local solver = eqn.solver

local common = require 'common'()
local xNames = common.xNames
local sym = common.sym
?>
#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 B = U.B;
	real3 D = U.D;
	real _1_mu = U._1_mu;
	real _1_eps = U._1_eps;
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.D = _real3(0., B.z * _1_mu, -B.y * _1_mu),
		.B = _real3(0., -D.z * _1_eps, D.y * _1_eps),
	<? elseif side == 1 then ?>
		.D = _real3(-B.z * _1_mu, 0., B.x * _1_mu),
		.B = _real3(D.z * _1_eps, 0., -D.x * _1_eps),
	<? elseif side == 2 then ?>
		.D = _real3(B.y * _1_mu, -B.x * _1_mu, 0.),
		.B = _real3(-D.y * _1_eps, D.x * _1_eps, 0.),
	<? end ?>
		.BPot = 0.,
		.sigma = 0.,
		._1_eps = 0.,
		._1_mu = 0.,
	};
}
<? end ?>

<?=eqn.eigen_t?> eigen_forInterface(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	real3 n
) {
	return (<?=eqn.eigen_t?>){
		.sqrt_1_eps = sqrt(.5 * (UL._1_eps + UR._1_eps)),
		.sqrt_1_mu = sqrt(.5 * (UL._1_mu + UR._1_mu)),
	};
}

/*
TODO update this for Einstein-Maxwell (take the metric into consideration
*/
<? for side=0,solver.dim-1 do ?>

<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=eqn.waves_t?> Y;

	const real ise = sqrt_1_2 * eig.sqrt_1_eps;
	const real isu = sqrt_1_2 * eig.sqrt_1_mu;

	<? if side==0 then ?>
	
	Y.ptr[0] = X.ptr[2] *  ise + X.ptr[4] * isu;
	Y.ptr[1] = X.ptr[1] * -ise + X.ptr[5] * isu;
	Y.ptr[2] = X.ptr[0] * -ise + X.ptr[3] * isu;
	Y.ptr[3] = X.ptr[0] *  ise + X.ptr[3] * isu;
	Y.ptr[4] = X.ptr[1] *  ise + X.ptr[5] * isu;
	Y.ptr[5] = X.ptr[2] * -ise + X.ptr[4] * isu;
	
	<? elseif side==1 then ?>
	
	Y.ptr[0] = X.ptr[0] *  ise + X.ptr[5] * isu;
	Y.ptr[1] = X.ptr[2] * -ise + X.ptr[3] * isu;
	Y.ptr[2] = X.ptr[1] * -ise + X.ptr[4] * isu;
	Y.ptr[3] = X.ptr[1] *  ise + X.ptr[4] * isu;
	Y.ptr[4] = X.ptr[2] *  ise + X.ptr[3] * isu;
	Y.ptr[5] = X.ptr[0] * -ise + X.ptr[5] * isu;
	
	<? elseif side==2 then ?>
	
	Y.ptr[0] = X.ptr[1] *  ise + X.ptr[3] * isu;
	Y.ptr[1] = X.ptr[0] * -ise + X.ptr[4] * isu;
	Y.ptr[2] = X.ptr[2] * -ise + X.ptr[5] * isu;
	Y.ptr[3] = X.ptr[2] *  ise + X.ptr[5] * isu;
	Y.ptr[4] = X.ptr[0] *  ise + X.ptr[4] * isu;
	Y.ptr[5] = X.ptr[1] * -ise + X.ptr[3] * isu;
	
	<? end ?>

	return Y;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;
	
	const real se = sqrt_1_2 / eig.sqrt_1_eps;
	const real su = sqrt_1_2 / eig.sqrt_1_mu;

	<? if side==0 then ?>
/*
z, -y, -x, x, y, -z
y,  z,  x, x, z, y
*/
	Y.ptr[0] = se * (-X.ptr[2] + X.ptr[3]);
	Y.ptr[1] = se * (-X.ptr[1] + X.ptr[4]);
	Y.ptr[2] = se * (X.ptr[0] + -X.ptr[5]);
	Y.ptr[3] = su * (X.ptr[2] + X.ptr[3]);
	Y.ptr[4] = su * (X.ptr[0] + X.ptr[5]);
	Y.ptr[5] = su * (X.ptr[1] + X.ptr[4]);
	
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
	Y.ptr[0] = se * (X.ptr[0] - X.ptr[5]);
	Y.ptr[1] = se * (-X.ptr[2] + X.ptr[3]);
	Y.ptr[2] = se * (-X.ptr[1] + X.ptr[4]);
	Y.ptr[3] = su * (X.ptr[1] + X.ptr[4]);
	Y.ptr[4] = su * (X.ptr[2] + X.ptr[3]);
	Y.ptr[5] = su * (X.ptr[0] + X.ptr[5]);
	
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
	Y.ptr[0] = se * (-X.ptr[1] + X.ptr[4]);
	Y.ptr[1] = se * (X.ptr[0] - X.ptr[5]);
	Y.ptr[2] = se * (-X.ptr[2] + X.ptr[3]);
	Y.ptr[3] = su * (X.ptr[0] + X.ptr[5]);
	Y.ptr[4] = su * (X.ptr[1] + X.ptr[4]);
	Y.ptr[5] = su * (X.ptr[2] + X.ptr[3]);
	
	<? end ?>
	
	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=eqn.cons_t?> Y;

	real3 D = X.D;
	real3 B = X.B;
	real _1_eps = eig.sqrt_1_eps * eig.sqrt_1_eps;
	real _1_mu = eig.sqrt_1_mu * eig.sqrt_1_mu;

	<? if side==0 then ?>
	
	Y.ptr[0] = 0;
	Y.ptr[1] = B.z * _1_mu;
	Y.ptr[2] = -B.y * _1_mu;
	Y.ptr[3] = 0;
	Y.ptr[4] = -D.z * _1_eps;
	Y.ptr[5] = D.y * _1_eps;

	<? elseif side==1 then ?>
		
	Y.ptr[0] = -B.z * _1_mu;
	Y.ptr[1] = 0;
	Y.ptr[2] = B.x * _1_mu;
	Y.ptr[3] = D.z * _1_eps;
	Y.ptr[4] = 0;
	Y.ptr[5] = -D.x * _1_eps;
		
	<? elseif side==2 then ?>
		
	Y.ptr[0] = B.y * _1_mu;
	Y.ptr[1] = -B.x * _1_mu;
	Y.ptr[2] = 0;
	Y.ptr[3] = -D.y * _1_eps;
	Y.ptr[4] = D.x * _1_eps;
	Y.ptr[5] = 0;
		
	<? end ?>
	
	for (int i = numWaves; i < numStates; ++i) {
		Y.ptr[i] = 0;
	}

	return Y;
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
	deriv->D = real3_sub(deriv->D, real3_scale(U->D, U->_1_eps * U->sigma));


	//for non-time-varying susceptibilities, here's the source term:
	//D_i,t ... = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_,j B_k - J_i
	//B_i,t ... = 1/sqrt(g) g_il epsBar^ljk (1/eps)_,j B_k

	real3 grad_1_mu = _real3(0,0,0);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = (
		U[stepsize.<?=xj?>]._1_mu
		- U[-stepsize.<?=xj?>]._1_mu
	) / grid_dx<?=j?>;
	<? end ?>
	
	real3 grad_1_eps = _real3(0,0,0);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>
	grad_1_mu.<?=xj?> = (
		U[stepsize.<?=xj?>]._1_eps
		- U[-stepsize.<?=xj?>]._1_eps
	) / grid_dx<?=j?>;
	<? end ?>

	real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
	<? for j=0,solver.dim-1 do 
		local xj = xNames[j+1] ?>{
		<?=eqn.cons_t?> flux = fluxFromCons_<?=j?>(*U, x);
		flux.D = real3_scale(coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = real3_scale(coord_lower(flux.B, x), _1_sqrt_det_g);
		deriv->D.<?=xj?> -= real3_dot(flux.D, grad_1_mu);
		deriv->B.<?=xj?> -= real3_dot(flux.B, grad_1_eps);
	}<? end ?>
}


//used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.eigen_t?> eig;
	eig.sqrt_1_eps = sqrt(U._1_eps);
	eig.sqrt_1_mu = sqrt(U._1_mu);
	return eig;
}
<? end ?>
