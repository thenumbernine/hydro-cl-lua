<? 
local common = require 'common'
local sym = common.sym
local xNames = common.xNames

local solver = eqn.solver
local scalar = eqn.scalar
local vec3 = eqn.vec3
?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;
//typedef <?=scalar?> scalar;
//typedef <?=vec3?> vec3;

<? for side=0,2 do ?>
real3 coord_g_uu<?=side?>(real3 x) {
	return _real3(<?
	for j=0,2 do
		?>coord_g_uu<?=side <= j and side..j or j..side?>(x)<?
		if j < 2 then ?>, <? end
	end ?>);
}
<? end ?>

<? for side=0,solver.dim-1 do 
	local xside = xNames[side+1]
?>
cons_t fluxFromCons_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	cons_t F;
	real alpha = metric_alpha(x);
	real3 beta_u = metric_beta_u(x);
	F.Pi = <?=scalar?>_real_mul(
		real_<?=scalar?>_add(
			<?=scalar?>_real_mul(U.Pi, beta_u.<?=xside?>), 
			<?=scalar?>_real_mul(<?=vec3?>_real3_dot(U.Psi_l, coord_g_uu<?=side?>(x)), alpha)
		), 
		-solver->wavespeed
	);
	F.Psi_l.x = <?=scalar?>_from_real(-solver->wavespeed * beta_u.<?=xside?>);
	F.Psi_l.y = <?=scalar?>_from_real(-solver->wavespeed * beta_u.<?=xside?>);
	F.Psi_l.z = <?=scalar?>_from_real(-solver->wavespeed * beta_u.<?=xside?>);
	F.Psi_l.s<?=side?> = <?=scalar?>_sub(F.Psi_l.s<?=side?>, <?=scalar?>_real_mul(U.Pi, solver->wavespeed * alpha));
	return F;
}
<? end ?>

<? for side=0,solver.dim-1 do 
	local xside = xNames[side+1]
?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x
) {
	real nLen = coord_sqrt_g_uu<?=side..side?>(x);
	real alpha_nLen = metric_alpha(x) * nLen;
	real3 beta_u = metric_beta_u(x);
	return (range_t){
		.min = solver->wavespeed * (-beta_u.<?=xside?> - alpha_nLen),
		.max = solver->wavespeed * (-beta_u.<?=xside?> + alpha_nLen),
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
	return (eigen_t){};
}

<? for side=0,solver.dim-1 do 
	local side1 = (side+1) % 3	-- solver.dim
	local side2 = (side+2) % 3	-- solver.dim
?>
waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) { 
	waves_t Y;

	sym3 g_uu = coord_g_uu(x);

	real nLen = coord_sqrt_g_uu<?=side..side?>(x);
	real nLenSq = nLen * nLen;
	real invDenom = 1. / nLenSq;

	<? if side == 0 then ?>
		//n = [1,0,0], n2 = [0,1,0], n3 = [0,0,1]
		//nU = [g^xx, g^xy, g^xz], n2U = [g^yx, g^yy, g^yz], n3U = [g^zx, g^zy, g^zz]
		Y.ptr[0] = .5 * invDenom * (
			X.ptr[0]
			+ X.ptr[1] * g_uu.xx
			+ X.ptr[2] * g_uu.xy
			+ X.ptr[3] * g_uu.xz
		);
		Y.ptr[1] = invDenom * (
			X.ptr[1] * g_uu.xy
			+ X.ptr[2] * g_uu.yy
			+ X.ptr[3] * g_uu.yz
		);
		Y.ptr[2] = invDenom * (
			X.ptr[1] * g_uu.xz
			+ X.ptr[2] * g_uu.yz
			+ X.ptr[3] * g_uu.zz
		);
		Y.ptr[3] = .5 * invDenom * (
			-X.ptr[0]
			+ X.ptr[1] * g_uu.xx
			+ X.ptr[2] * g_uu.xy
			+ X.ptr[3] * g_uu.xz
		);
	<? elseif side == 1 then ?>
		//n = [0,1,0], n2 = [0,0,1], n3 = [1,0,0]
		//nU = [g^yx, g^yy, g^yz], n2U = [g^zx, g^zy, g^zz], n3U = [g^xx, g^xy, g^xz]
		Y.ptr[0] = .5 * invDenom * (
			X.ptr[0]
			+ X.ptr[1] * g_uu.xy
			+ X.ptr[2] * g_uu.yy
			+ X.ptr[3] * g_uu.yz
		);
		Y.ptr[1] = invDenom * (
			X.ptr[1] * g_uu.xz
			+ X.ptr[2] * g_uu.yz
			+ X.ptr[3] * g_uu.zz
		);
		Y.ptr[2] = invDenom * (
			X.ptr[1] * g_uu.xx
			+ X.ptr[2] * g_uu.xy
			+ X.ptr[3] * g_uu.xz
		);
		Y.ptr[3] = .5 * invDenom * (
			-X.ptr[0]
			+ X.ptr[1] * g_uu.xy
			+ X.ptr[2] * g_uu.yy
			+ X.ptr[3] * g_uu.yz
		);
	<? elseif side == 2 then ?>
		//n = [0,0,1], n2 = [1,0,0], n3 = [0,1,0]
		//nU = [g^zx, g^zy, g^zz], n2U = [g^xx, g^xy, g^xz], n3U = [g^yx, g^yy, g^yz]
		Y.ptr[0] = .5 * invDenom * (
			X.ptr[0]
			+ X.ptr[1] * g_uu.xz
			+ X.ptr[2] * g_uu.yz
			+ X.ptr[3] * g_uu.zz
		);
		Y.ptr[1] = invDenom * (
			X.ptr[1] * g_uu.xx
			+ X.ptr[2] * g_uu.xy
			+ X.ptr[3] * g_uu.xz
		);
		Y.ptr[2] = invDenom * (
			X.ptr[1] * g_uu.xy
			+ X.ptr[2] * g_uu.yy
			+ X.ptr[3] * g_uu.yz
		);
		Y.ptr[3] = .5 * invDenom * (
			-X.ptr[0]
			+ X.ptr[1] * g_uu.xz
			+ X.ptr[2] * g_uu.yz
			+ X.ptr[3] * g_uu.yz
		);
	<? end ?>

	return Y;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x
) {
	<?=scalar?>* X = (<?=scalar?>*)X_.ptr;
	cons_t Y;
	
	real nLen = coord_sqrt_g_uu<?=side..side?>(x);
	real nLenSq = nLen * nLen;

	<? if side == 0 then ?>
		//n = [1,0,0], n2 = [0,1,0], n3 = [0,0,1]
		Y.ptr[0] = X[0] - X[3];
		Y.ptr[1] = X[0] + X[3];
		Y.ptr[2] = X[1];
		Y.ptr[3] = X[2];
	<? elseif side == 1 then ?>
		//n = [0,1,0], n2 = [0,0,1], n3 = [1,0,0]
		Y.ptr[0] = X[0] - X[3];
		Y.ptr[1] = X[2];
		Y.ptr[2] = X[0] + X[3];
		Y.ptr[3] = X[1];
	<? elseif side == 2 then ?>
		//n = [0,0,1], n2 = [1,0,0], n3 = [0,1,0]
		Y.ptr[0] = X[0] - X[3];
		Y.ptr[1] = X[1];
		Y.ptr[2] = X[2];
		Y.ptr[3] = X[0] + X[3];
	<? end ?>
	
	return Y;
}

// What's the difference between eigen_fluxTransform and fluxFromCons?
// The difference is that the flux matrix of this is based on 'eig', which is derived from U's ... especially UL & UR in the case of the Roe solver
// whereas that of fluxFromCons is based purely on 'U'.
// Since eqn/wave has no eigen_t info derived from U, the two functions are identical.
cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	return fluxFromCons_<?=side?>(solver, X, x);
}

<? end ?>

<? for side=0,solver.dim-1 do ?>
eigen_t eigen_forCell_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	return (eigen_t){};
}
<? end ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
#if 0
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;
	//TODO make use of this
	//real c = solver->wavespeed / unit_m_per_s;

	real alpha = metric_alpha(x);
	real K = metric_K(x);
	real3x3 partial_beta_ul = metric_partial_beta_ul(x);
	real3 partial_alpha_l = metric_partial_alpha_l(x);
	real3 conn23 = coord_conn_trace23(x);
	real f = metric_f(x);

	real3 Psi_u = coord_raise(U->Psi_l, x);

	deriv->Pi += 
		real3_dot(partial_alpha_l, Psi_u)
		+ alpha * K * U->Pi
		- alpha * real3_dot(U->Psi_l, conn23)
		- alpha * f 						//... for □Φ=f
	;

	deriv->Psi_l = real3_add3(
		deriv->Psi_l,
		real3_real3x3_mul(
			deriv->Psi_l,
			partial_beta_ul
		),
		real3_real_mul(partial_alpha_l, U->Pi)
	);
#endif
}
