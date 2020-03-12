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
	F.Pi = <?=scalar?>_real_mul(real_<?=scalar?>_add(beta_u.<?=xside?>, <?=vec3?>_real3_dot(U.Psi_l, coord_g_uu<?=side?>(x))), -solver->wavespeed),
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
	waves_t Y_;
	<?=scalar?>* Y = (<?=scalar?>*)Y_.ptr;
	
	real sqrt_gUii = coord_sqrt_g_uu<?=side..side?>(x);
	real gUii = sqrt_gUii * sqrt_gUii;
	<?=scalar?> Psi_u = <?=vec3?>_real3_dot(X.Psi_l, coord_g_uu<?=side?>(x));
	Y[0] = <?=scalar?>_add(
		<?=scalar?>_real_mul(Psi_u, .5 / gUii), 
		<?=scalar?>_real_mul(X.Pi, .5 / sqrt_gUii)
	);
	Y[1] = <?=scalar?>_real_mul(X.Psi_l.s<?=side1?>, 1. / gUii);
	Y[2] = <?=scalar?>_real_mul(X.Psi_l.s<?=side2?>, 1. / gUii);
	Y[3] = <?=scalar?>_sub(
		<?=scalar?>_real_mul(Psi_u, .5 / gUii),
		<?=scalar?>_real_mul(X.Pi, .5 / sqrt_gUii)
	);
	return Y_;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x
) {
	<?=scalar?>* X = (<?=scalar?>*)X_.ptr;
	cons_t Y;
	
	real sqrt_gUii = coord_sqrt_g_uu<?=side..side?>(x);
	real gUii = sqrt_gUii * sqrt_gUii;
	
	Y.Pi = <?=scalar?>_real_mul(
		<?=scalar?>_sub(X[0], X[3]),
		sqrt_gUii
	);
	
	Y.Psi_l.s<?=side?> = <?=scalar?>_sub(
		<?=scalar?>_add(X[0], X[3]),
		<?=scalar?>_add(
			<?=scalar?>_real_mul(X[1], coord_g_uu<?=side <= side1 and side..side1 or side1..side?>(x)), 
			<?=scalar?>_real_mul(X[2], coord_g_uu<?=side <= side2 and side..side2 or side2..side?>(x))
		)
	);
	
	Y.Psi_l.s<?=side1?> = <?=scalar?>_real_mul(X[1], gUii);
	
	Y.Psi_l.s<?=side2?> = <?=scalar?>_real_mul(X[2], gUii);
	
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
<? 
if solver.coord.vectorComponent ~= 'cartesian' 
and not require 'coord.cartesian'.is(solver.coord)
then ?>
	
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;
	real c = solver->wavespeed / unit_m_per_s;

<? 
	if not solver.coord.vectorComponent == 'anholonomic' 
	and not eqn.weightFluxByGridVolume 
	then 
?>

	real alpha = metric_alpha(x);
	real K = metric_K(x);
	real3x3 partial_beta_ul = metric_partial_beta_ul(x);
	real3 partial_alpha_l = metric_partial_alpha_l(x);
	real3 conn23 = coord_conn_trace23(x);

	real3 Psi_u = coord_raise(U->Psi_l, x);

	deriv->Pi += 
		real3_dot(partial_alpha_l, Psi_u)
		+ alpha * K * U->Pi
		- alpha * real3_dot(U->Psi_l, conn23)
		//- alpha * f 						//... for □ Φ = f
	;

	deriv->Psi_l = real3_add3(
		deriv->Psi_l,
		real3_real3x3_mul(
			deriv->Psi_l,
			partial_beta_ul
		),
		real3_real_mul(partial_alpha_l, U->Pi)
	);


<? 
	elseif not solver.coord.vectorComponent == 'anholonomic' 
	and eqn.weightFluxByGridVolume 
	then
?>
	
	real3 conn12 = coord_conn_trace12(x);
	deriv->Psi_l.x = <?=scalar?>_sub(deriv->Psi_l.x, <?=scalar?>_real_mul(U->Pi, c * conn12.x));
	deriv->Psi_l.y = <?=scalar?>_sub(deriv->Psi_l.y, <?=scalar?>_real_mul(U->Pi, c * conn12.y));
	deriv->Psi_l.z = <?=scalar?>_sub(deriv->Psi_l.z, <?=scalar?>_real_mul(U->Pi, c * conn12.z));

<? 
	elseif solver.coord.vectorComponent == 'anholonomic' 
	and not eqn.weightFluxByGridVolume 
	then 
?>

#error I haven't calculated this.  Feel free to comment it out and run it anyways. 

<? 	end ?>

<? 
end	-- solver.coord.vectorComponent ~= 'cartesian' 
?>
}
