<? 
local common = require 'common'
local sym = common.sym
local xNames = common.xNames

local solver = eqn.solver
?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

<? for side=0,2 do ?>
static inline real3 coord_g_uu<?=side?>(real3 x) {
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
	F.phi_t = -solver->wavespeed * (beta_u.<?=xside?> + real3_dot(coord_g_uu<?=side?>(x), U.phi_i)),
	F.phi_i = _real3(
		-solver->wavespeed * beta_u.<?=xside?>,
		-solver->wavespeed * beta_u.<?=xside?>,
		-solver->wavespeed * beta_u.<?=xside?>
	);
	F.phi_i.s<?=side?> -= solver->wavespeed * alpha * U.phi_t;
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
	real alpha_sqrt_gUii = metric_alpha(x) * coord_sqrt_g_uu<?=side..side?>(x);
	real3 beta_u = metric_beta_u(x);
	return (range_t){
		.min = solver->wavespeed * (-beta_u.<?=xside?> - alpha_sqrt_gUii),
		.max = solver->wavespeed * (-beta_u.<?=xside?> + alpha_sqrt_gUii),
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
	
	real sqrt_gUii = coord_sqrt_g_uu<?=side..side?>(x);
	real gUii = sqrt_gUii * sqrt_gUii;
	real phiUi = real3_dot(X.phi_i, coord_g_uu<?=side?>(x));
	Y.ptr[0] = .5 * (phiUi / gUii + X.phi_t / sqrt_gUii);
	Y.ptr[1] = X.phi_i.s<?=side1?> / gUii;
	Y.ptr[2] = X.phi_i.s<?=side2?> / gUii;
	Y.ptr[3] = .5 * (phiUi / gUii - X.phi_t / sqrt_gUii);
	return Y;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x
) {
	cons_t Y;
	
	real sqrt_gUii = coord_sqrt_g_uu<?=side..side?>(x);
	real gUii = sqrt_gUii * sqrt_gUii;
	Y.phi_t = sqrt_gUii * (X.ptr[0] - X.ptr[3]);
	Y.phi_i.s<?=side?> = X.ptr[0] + X.ptr[3] - (
		coord_g_uu<?=side <= side1 and side..side1 or side1..side?>(x) * X.ptr[1] 
		+ coord_g_uu<?=side <= side2 and side..side2 or side2..side?>(x) * X.ptr[2]
	);
	Y.phi_i.s<?=side1?> = gUii * X.ptr[1];
	Y.phi_i.s<?=side2?> = gUii * X.ptr[2];
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
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;
	real c = solver->wavespeed / unit_m_per_s;

<? if not solver.coord.anholonomic then ?>
<? if not eqn.weightFluxByGridVolume then ?>

	real alpha = metric_alpha(x);
	real dalpha_t = metric_dalpha_t(x);
	real K = metric_K(x);
	real3x3 dbeta_ul = metric_dbeta_ul(x);
	real3 dalpha_l = metric_dalpha_l(x);
	
	/*	
	\Pi (
		\frac{1}{\alpha} \alpha_{,t} (1 - \frac{1}{\alpha})
		+ K \alpha 
	)
	
	+ \Psi_i (
		\alpha_{,j} \gamma^{ij}
		- \alpha \cdot {}^{(3)} \Gamma^i 
	)

	- \alpha \frac{dV}{d|\Phi|^2} \Phi
	*/
	deriv->phi_t += c * (
		U->phi_t * (
			dalpha_t * (1. - 1. / alpha) / alpha
			+ alpha * K
		)
		+ real3_dot(coord_raise(dalpha_l, x), U->phi_i)
		- alpha * (
			real3_dot(coord_conn_trace23(x), U->phi_i)
			+ eqn_source(x)
		)
	);

	//\alpha_{,i} \Pi + {\beta^k}_{,i} \Psi_k
	deriv->phi_i = real3_add(deriv->phi_i, 
		real3_real_mul(
			real3_add(
				real3_real_mul(dalpha_l, U->phi_t),
				real3x3_real3_mul(dbeta_ul, U->phi_i)
			),
			c
		)
	);

<? else ?>
	real3 conn12 = coord_conn_trace12(x);
	deriv->phi_i.x -= c * conn12.x * U->phi_t;
	deriv->phi_i.y -= c * conn12.y * U->phi_t;
	deriv->phi_i.z -= c * conn12.z * U->phi_t;
<? end ?>
<? end -- anholonomic ?>
}
