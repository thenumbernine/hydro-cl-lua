<? 
local common = require 'common'()
local sym = common.sym

local solver = eqn.solver
?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

<? for side=0,2 do ?>
static inline real3 coord_gU<?=side?>(real3 x) {
	return _real3(<?
	for j=0,2 do
		?>coord_gU<?=side <= j and side..j or j..side?>(x)<?
		if j < 2 then ?>, <? end
	end ?>);
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
cons_t fluxFromCons_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	cons_t F;
	F.phi_t = -solver->wavespeed * solver->wavespeed * real3_dot(coord_gU<?=side?>(x), U.phi_i),
	F.phi_i = _real3(0,0,0);
	F.phi_i.s<?=side?> = -U.phi_t;
	return F;
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x
) {
	real c_sqrt_gUii = solver->wavespeed * coord_sqrt_gU<?=side..side?>(x);
	return (range_t){
		.min = -c_sqrt_gUii,
		.max = c_sqrt_gUii,
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
	real c = solver->wavespeed;
	real cSq = c * c;
	
	real sqrt_gUii = coord_sqrt_gU<?=side..side?>(x);
	real gUii = sqrt_gUii * sqrt_gUii;
	real phiUi = real3_dot(X.phi_i, coord_gU<?=side?>(x));
	Y.ptr[0] = .5 * (phiUi / gUii + X.phi_t / (c * sqrt_gUii));
	Y.ptr[1] = X.phi_i.s<?=side1?> / (cSq * gUii);
	Y.ptr[2] = X.phi_i.s<?=side2?> / (cSq * gUii);
	Y.ptr[3] = .5 * (phiUi / gUii - X.phi_t / (c * sqrt_gUii));
	
	return Y;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x
) {
	cons_t Y;
	real c = solver->wavespeed;
	real cSq = c * c;
	
	real sqrt_gUii = coord_sqrt_gU<?=side..side?>(x);
	real gUii = sqrt_gUii * sqrt_gUii;
	Y.phi_t = c * sqrt_gUii * (X.ptr[0] - X.ptr[3]);
	Y.phi_i.s<?=side?> = X.ptr[0] + X.ptr[3] - cSq * (
		coord_gU<?=side <= side1 and side..side1 or side1..side?>(x) * X.ptr[1] 
		+ coord_gU<?=side <= side2 and side..side2 or side2..side?>(x) * X.ptr[2]
	);
	Y.phi_i.s<?=side1?> = cSq * gUii * X.ptr[1];
	Y.phi_i.s<?=side2?> = cSq * gUii * X.ptr[2];
	
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

