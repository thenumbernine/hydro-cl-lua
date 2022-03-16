<?
local useConstrainU = true -- constrains alpha to alphamin and calcs H and M^i
local useAddSource = true
local useKreissOligarDissipation = true	-- depends on useAddSource
?>
//// MODULE_NAME: <?=calc_gamma_ll?>

#define <?=calc_gamma_ll?>(U, x)	((U)->gamma_ll)

//// MODULE_NAME: <?=calc_gamma_uu?>
//// MODULE_DEPENDS: <?=cons_t?>

static inline sym3 <?=calc_gamma_uu?>(
	global <?=cons_t?> const * const U,
	real3 const pt
) {
	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	return gamma_uu;
}

//// MODULE_NAME: <?=calc_d_lll?>
//// MODULE_DEPENDS: _3sym3

static inline _3sym3 <?=calc_d_lll?>(
	global <?=cons_t?> const * const U,
	real3 const pt
) {
	_3sym3 const dHat_lll = _3sym3_zero;	//TODO from the metric partial
	return _3sym3_add(U->dDelta_lll, dHat_lll);
}

//// MODULE_NAME: <?=calcFromGrad_a_l?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?>

//NOTICE THIS DOES NOT BOUNDS CHECK
real3 <?=calcFromGrad_a_l?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U
) {
	real3 a_l;
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>	{
		global <?=cons_t?> const * const UL = U - solver->stepsize.<?=xi?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.<?=xi?>;
		a_l.<?=xi?> = (log(UR->alpha) - log(UL->alpha)) / (2. * solver->grid_dx.s<?=i-1?>);
	}
<?
end
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>	a_l.<?=xi?> = 0;
<?
end
?>
	return a_l;
}

//// MODULE_NAME: <?=calcFromGrad_d_lll?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?>

//NOTICE THIS DOES NOT BOUNDS CHECK
_3sym3 <?=calcFromGrad_d_lll?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	_3sym3 d_lll;
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>	{
		global <?=cons_t?> const * const UL = U - solver->stepsize.<?=xi?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.<?=xi?>;
		// needed for dHat_lll
		//global <?=cell_t?> const * const cellL = cell - solver->stepsize.<?=xi?>;
		//global <?=cell_t?> const * const cellR = cell + solver->stepsize.<?=xi?>;
<? 	for jk,xjk in ipairs(symNames) do
?>		d_lll.<?=xi?>.<?=xjk?> = .5 * (UR->gamma_ll.<?=xjk?> - UL->gamma_ll.<?=xjk?>) / (2. * solver->grid_dx.s<?=i-1?>);
<? 	end
?>	}
<?
end
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>	d_lll.<?=xi?> = sym3_zero;
<?
end
?>
	return d_lll;
}

//// MODULE_NAME: <?=calcFromGrad_b_ul?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?>

//NOTICE THIS DOES NOT BOUNDS CHECK
real3x3 <?=calcFromGrad_b_ul?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U
) {
	real3x3 b_ul;
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>	{
		global <?=cons_t?> const * const UL = U - solver->stepsize.<?=xi?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.<?=xi?>;
<?	for j,xj in ipairs(xNames) do
?>		b_ul.<?=xj?>.<?=xi?> = (UR->beta_u.<?=xj?> - UL->beta_u.<?=xj?>) / (2. * solver->grid_dx.s<?=i-1?>);
<?	end
?>	}
<?
end
for i=solver.dim+1,3 do
	local xi = xNames[i]
	for j,xj in ipairs(xNames) do
?>	b_ul.<?=xj?>.<?=xi?> = 0;
<?	end
end
?>
	return b_ul;

}

//// MODULE_NAME: <?=initDeriv_numeric_and_useBSSNVars?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;

//// MODULE_DEPENDS: <?=calcFromGrad_a_l?>
	U->a_l = <?=calcFromGrad_a_l?>(solver, U);
//// MODULE_DEPENDS: <?=calcFromGrad_d_lll?>
	_3sym3 const dHat_lll = _3sym3_zero;
	_3sym3 const d_lll = <?=calcFromGrad_d_lll?>(solver, U, cell);
	U->dDelta_lll = _3sym3_sub(d_lll, dHat_lll);
<? if eqn.consStruct.vars:find(nil, function(var) return var.name == "b_ul" end) then ?>
//// MODULE_DEPENDS: <?=calcFromGrad_b_ul?>
	U->b_ul = <?=calcFromGrad_b_ul?>(solver, U);
<? end ?>
}

<?
if eqn.initCond.initAnalytical then
	error("TODO - can't handle analytical initial conditions yet")
end
?>

<? if eqn.initCond.useBSSNVars then ?>

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=coordMap?> <?=coord_gHol_ll?> <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real W = 1.;
	real K = 0.;		// I call this 'tr_K' elsewhere in this file to avoid ambiguity, but in useBSSNVars init it is just 'K'
	real3 LambdaBar_U = real3_zero;
	real3 beta_U = real3_zero;
	real3 B_U = real3_zero;
	sym3 epsilon_LL = sym3_zero;
	sym3 ABar_LL = sym3_zero;

	real rho = 0.;

	<?=initCode()?>

	//for (int i = 0; i < numStates; ++i) {
	//	U->ptr[i] = 0.;
	//}
	*U = (<?=cons_t?>){.ptr={ 0. / 0. }};

	U->alpha = alpha;

//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?>
	// gammaHat_IJ = delta_IJ
	// gamma_ij = e_i^I e_j^J (epsilon_IJ + gammaHat_IJ) / W^2
	sym3 const gammaBar_LL = sym3_add(epsilon_LL, sym3_ident);
	sym3 const gamma_LL = sym3_real_mul(gammaBar_LL, 1. / (W*W));
	U->gamma_ll = sym3_rescaleToCoord_LL(gamma_LL, x);
	
	// K_ij = e_i^I e_j^J (ABar_IJ + gammaBar_IJ K/3) / W^2
	U->K_ll = sym3_rescaleToCoord_LL(
		sym3_add(
			sym3_real_mul(ABar_LL, 1. / (W*W)),
			sym3_real_mul(gamma_LL, K / 3.)
		), x);

	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if eqn.consStruct.vars:find(nil, function(var) return var.name == "beta_u" end) then 
?>	U->beta_u = real3_rescaleToCoord_U(beta_U, x);
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "betaLap_u" end) then 
?>	U->betaLap_u = real3_zero;
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end) then 
?>	U->B_u = real3_zero;
<? end
?>
/*<?--[[
b_ul is initialized in initDerivs, beta_u is initialized in the bssn vars
NOTICE LambdaBar_U isn't used
it is defined as _Λ^i = ΔΓ^i = ΔΓ^i_jk _γ^jk = (_Γ^i_jk - ^Γ^i_jk) _γ^jk
where _γ_ij is the conformal metric, _Γ^i_jk is the conformal connection, ^Γ^i_jk is the background (grid) metric connection
--]]?>*/

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	
	U->H = 0;
	U->M_u = real3_zero;
}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=initDeriv_numeric_and_useBSSNVars?>

<? else	-- not eqn.initCond.useBSSNVars ?>

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=coordMap?> <?=coord_gHol_ll?> <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real3 beta_u = real3_zero;

	sym3 gamma_ll = coord_gHol_ll(x);

	sym3 K_ll = sym3_zero;

	//TODO more stress-energy vars
	real rho = 0.;

	<?=initCode()?>

	*U = (<?=cons_t?>){.ptr={ 0. / 0. }};
	
	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;

	//Z_u n^u = 0
	//Theta = alpha n_u Z^u = alpha Z^u
	//for n_a = (-alpha, 0)
	//n^a_l = (1/alpha, -beta^i/alpha)
	//(Z_t - Z_i beta^i) / alpha = Theta ... = ?
	//Z^t n_t + Z^i n_i = -alpha Z^t = Theta
	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if eqn.consStruct.vars:find(nil, function(var) return var.name == "beta_u" end) then 
?>	U->beta_u = beta_u;
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "betaLap_u" end) then 
?>	U->betaLap_u = real3_zero;
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end) then 
?>	U->B_u = real3_zero;
<? end
?>

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	
	U->H = 0;
	U->M_u = real3_zero;
}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=initDeriv_numeric_and_useBSSNVars?>

<? end	-- eqn.initCond.useBSSNVars ?>

//// MODULE_NAME: <?=setFlatSpace?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?>

static inline void <?=setFlatSpace?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	(U)->alpha = 1.;
	(U)->gamma_ll = sym3_ident;
	(U)->a_l = real3_zero;
	(U)->dDelta_lll = _3sym3_zero;
	(U)->K_ll = sym3_zero;
	(U)->Theta = 0.;
	(U)->Z_l = real3_zero;

<? if eqn.consStruct.vars:find(nil, function(var) return var.name == "beta_u" end) then 
?>	(U)->beta_u = beta_u;
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "betaLap_u" end) then 
?>	(U)->betaLap_u = real3_zero;
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end) then 
?>	(U)->B_u = real3_zero;
<? end
if eqn.consStruct.vars:find(nil, function(var) return var.name == "b_ul" end) then 
?>	(U)->b_ul = real3x3_zero;
<? end
?>

<? if eqn.useStressEnergyTerms then ?>
	//what to do with the constraint vars and the source vars?
	(U)->rho = 0;
	(U)->S_u = real3_zero;
	(U)->S_ll = sym3_zero;
<? end ?>
	
	(U)->H = 0;
	(U)->M_u = real3_zero;
}

//// MODULE_NAME: <?=calcDTCell?>
//// MODULE_DEPENDS: <?=SETBOUNDS?> <?=cons_t?> <?=initCond_codeprefix?>

#define <?=calcDTCell?>(\
	/*global real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell\
) {\
	/* the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma_ij) is only called once */\
	real const f_alphaSq = calc_f_alphaSq(U->alpha);\
	real const det_gamma = sym3_det(U->gamma_ll);\
	real const alpha_sqrt_f = sqrt(f_alphaSq);\
\
	<? for side=0,solver.dim-1 do ?>{\
\
		<? if side == 0 then ?>\
		real const gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;\
		<? elseif side == 1 then ?>\
		real const gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;\
		<? elseif side == 2 then ?>\
		real const gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;\
		<? end ?>	\
		real const sqrt_gammaUjj = sqrt(gammaUjj);\
		real const lambdaLight = sqrt_gammaUjj * U->alpha;\
		real const lambdaGauge = sqrt_gammaUjj * alpha_sqrt_f;\
		real const lambda = (real)max(lambdaGauge, lambdaLight);\
\
		<? if eqn.useShift ~= "none" then ?>\
		real const betaUi = U->beta_u.s<?=side?>;\
		<? else ?>\
		real const betaUi = 0.;\
		<? end ?>\
\
		real const lambdaMin = (real)min((real)0., -betaUi - lambda);\
		real const lambdaMax = (real)max((real)0., -betaUi + lambda);\
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
		absLambdaMax = max((real)1e-9, absLambdaMax);\
		*(dt) = (real)min(*(dt), solver->grid_dx.s<?=side?> / absLambdaMax);\
	}<? end ?>\
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=cons_t?> <?=solver_t?> <?=normal_t?> rotate sym3_rotate real3x3_rotate _3sym3_rotate

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const f_alpha = calc_f_alpha((U)->alpha);\
\
	real const det_gamma = sym3_det((U)->gamma_ll);\
\
	real alpha = (U)->alpha;\
	real Theta = (U)->Theta;\
\
	real3 Z_l = (U)->Z_l;\
	real3 a_l = (U)->a_l;\
	sym3 gamma_ll = (U)->gamma_ll;\
	sym3 K_ll = (U)->K_ll;\
	_3sym3 dDelta_lll = (U)->dDelta_lll;\
	_3sym3 dHat_lll = _3sym3_zero;\
<? if eqn.useShift ~= "none" then --\
?>	real3 beta_u = (U)->beta_u;\
<? end --\
if eqn.useShift == "2005 Bona / 2008 Yano" then --\
?>	real3x3 b_ul = (U)->b_ul;\
<? end --\
?>\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
\
/*  I'm using .side for holonomic(coordinate) and anholonomic(orthonormal) */\
/* but for cartesian vector componets there is no .side, just .n, which is covariant iirc */\
/* and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes) */\
/* so here I'm going to just wing it */\
	real3 const n_l = normal_l1(n);\
\
	Z_l = real3_rotateFrom(Z_l, n_l);\
	a_l = real3_rotateFrom(a_l, n_l);\
	gamma_ll = sym3_rotateFrom(gamma_ll, n_l);\
	K_ll = sym3_rotateFrom(K_ll, n_l);\
	dDelta_lll = _3sym3_rotateFrom(dDelta_lll, n_l);\
	dHat_lll = _3sym3_rotateFrom(dHat_lll, n_l);\
<? if eqn.useShift ~= "none" then --\
?>	beta_u = real3_rotateFrom(beta_u);\
<? end --\
if eqn.useShift == "2005 Bona / 2008 Yano" then --\
?>	b_ul = real3x3_rotateFrom(b_ul);\
<? end --\
?>\
\
<? else ?>\
\
	if (false) {}\
	<? for side=0,solver.dim-1 do ?>\
	else if (n.side == <?=side?>) {\
		Z_l = real3_swap<?=side?>(Z_l);\
		a_l = real3_swap<?=side?>(a_l);\
		gamma_ll = sym3_swap<?=side?>(gamma_ll);\
		K_ll = sym3_swap<?=side?>(K_ll);\
		dDelta_lll = _3sym3_swap<?=side?>(dDelta_lll);\
		dHat_lll = _3sym3_swap<?=side?>(dHat_lll);\
	}\
	<? end ?>\
	else {\
		alpha = 0./0.;\
		Theta = 0./0.;\
<? for k,xk in ipairs(xNames) do --\
?>		a_l.<?=xk?> = 0./0.;\
		Z_l.<?=xk?> = 0./0.;\
<? end --\
?>\
<? for ij,xij in ipairs(symNames) do --\
?>		K_ll.<?=xij?> = 0./0.;\
		gamma_ll.<?=xij?> = 0./0.;\
<? end --\
?>\
<? for k,xk in ipairs(xNames) do --\
	for ij,xij in ipairs(symNames) do --\
?>		dDelta_lll.<?=xk?>.<?=xij?> = 0./0.;\
		dHat_lll.<?=xk?>.<?=xij?> = 0./0.;\
<?	end --\
end --\
?>	}\
\
<? end ?>\
\
	sym3 const gamma_uu = sym3_inv(gamma_ll, det_gamma);\
\
	for (int i = 0; i < numStates; ++i) {\
		(resultFlux)->ptr[i] = 0./0.;\
	}\
\
	_3sym3 d_lll = _3sym3_add(dDelta_lll, dHat_lll);\
<? if false then ?>\
	/* BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/z4_noZeroRows.html */\
	real const tmp1 = K_ll.xy * gamma_uu.xy;\
	real const tmp2 = K_ll.xz * gamma_uu.xz;\
	real const tmp3 = K_ll.yz * gamma_uu.yz;\
	real const tmp4 = 2. * tmp3;\
	real const tmp6 = K_ll.yy * gamma_uu.yy;\
	real const tmp7 = K_ll.zz * gamma_uu.zz;\
	real const tmp103 = gamma_uu.xx * gamma_uu.yy;\
	real const tmp106 = gamma_uu.xy * gamma_uu.xy;\
	real const tmp107 = gamma_uu.xx * gamma_uu.yz;\
	real const tmp111 = gamma_uu.xy * gamma_uu.xz;\
	real const tmp113 = gamma_uu.xx * gamma_uu.zz;\
	real const tmp116 = gamma_uu.xz * gamma_uu.xz;\
	real const tmp125 = gamma_uu.xy * gamma_uu.yz;\
	real const tmp128 = gamma_uu.xz * gamma_uu.yy;\
	real const tmp129 = gamma_uu.xy * gamma_uu.zz;\
	real const tmp132 = gamma_uu.xz * gamma_uu.yz;\
	(resultFlux)->alpha = 0.;\
	(resultFlux)->gamma_ll.xx = 0.;\
	(resultFlux)->gamma_ll.xy = 0.;\
	(resultFlux)->gamma_ll.xz = 0.;\
	(resultFlux)->gamma_ll.yy = 0.;\
	(resultFlux)->gamma_ll.yz = 0.;\
	(resultFlux)->gamma_ll.zz = 0.;\
	(resultFlux)->a_l.x = f_alpha * (tmp6 + tmp7 + K_ll.xx * gamma_uu.xx + tmp4 - 2. * Theta + 2. * tmp2 + 2. * tmp1);\
	(resultFlux)->a_l.y = 0.;\
	(resultFlux)->a_l.z = 0.;\
	(resultFlux)->dDelta_lll.x.xx = K_ll.xx * alpha;\
	(resultFlux)->dDelta_lll.x.xy = K_ll.xy * alpha;\
	(resultFlux)->dDelta_lll.x.xz = K_ll.xz * alpha;\
	(resultFlux)->dDelta_lll.x.yy = K_ll.yy * alpha;\
	(resultFlux)->dDelta_lll.x.yz = K_ll.yz * alpha;\
	(resultFlux)->dDelta_lll.x.zz = K_ll.zz * alpha;\
	(resultFlux)->dDelta_lll.y.xx = 0.;\
	(resultFlux)->dDelta_lll.y.xy = 0.;\
	(resultFlux)->dDelta_lll.y.xz = 0.;\
	(resultFlux)->dDelta_lll.y.yy = 0.;\
	(resultFlux)->dDelta_lll.y.yz = 0.;\
	(resultFlux)->dDelta_lll.y.zz = 0.;\
	(resultFlux)->dDelta_lll.z.xx = 0.;\
	(resultFlux)->dDelta_lll.z.xy = 0.;\
	(resultFlux)->dDelta_lll.z.xz = 0.;\
	(resultFlux)->dDelta_lll.z.yy = 0.;\
	(resultFlux)->dDelta_lll.z.yz = 0.;\
	(resultFlux)->dDelta_lll.z.zz = 0.;\
	(resultFlux)->K_ll.xx = -alpha * (2. * d_lll.z.xz * gamma_uu.zz + 2. * d_lll.z.xy * gamma_uu.yz + d_lll.z.xx * gamma_uu.xz + 2. * d_lll.y.xz * gamma_uu.yz + 2. * d_lll.y.xy * gamma_uu.yy + d_lll.y.xx * gamma_uu.xy + 2. * Z_l.x - a_l.x - d_lll.x.yy * gamma_uu.yy - 2. * d_lll.x.yz * gamma_uu.yz - d_lll.x.zz * gamma_uu.zz);\
	(resultFlux)->K_ll.xy = (-alpha * (2. * d_lll.z.yz * gamma_uu.zz + 2. * d_lll.z.yy * gamma_uu.yz + 2. * d_lll.y.yz * gamma_uu.yz + 2. * d_lll.y.yy * gamma_uu.yy + 2. * d_lll.x.yz * gamma_uu.xz + 2. * d_lll.x.yy * gamma_uu.xy + 2. * Z_l.y - a_l.y)) / 2.;\
	(resultFlux)->K_ll.xz = (-alpha * (2. * d_lll.z.zz * gamma_uu.zz + 2. * d_lll.z.yz * gamma_uu.yz + 2. * d_lll.y.zz * gamma_uu.yz + 2. * d_lll.y.yz * gamma_uu.yy + 2. * d_lll.x.zz * gamma_uu.xz + 2. * d_lll.x.yz * gamma_uu.xy + 2. * Z_l.z - a_l.z)) / 2.;\
	(resultFlux)->K_ll.yy = alpha * (d_lll.z.yy * gamma_uu.xz + d_lll.y.yy * gamma_uu.xy + d_lll.x.yy * gamma_uu.xx);\
	(resultFlux)->K_ll.yz = alpha * (d_lll.z.yz * gamma_uu.xz + d_lll.y.yz * gamma_uu.xy + d_lll.x.yz * gamma_uu.xx);\
	(resultFlux)->K_ll.zz = alpha * (d_lll.z.zz * gamma_uu.xz + d_lll.y.zz * gamma_uu.xy + d_lll.x.zz * gamma_uu.xx);\
	(resultFlux)->Theta = -alpha * (d_lll.z.yz * tmp129 - d_lll.z.yz * tmp132 + d_lll.z.yy * tmp125 - d_lll.z.yy * tmp128 + d_lll.z.xz * tmp113 - d_lll.z.xz * tmp116 + d_lll.z.xy * tmp107 - d_lll.z.xy * tmp111 + d_lll.y.zz * tmp132 + d_lll.y.yz * tmp128 - d_lll.y.zz * tmp129 + d_lll.y.xz * tmp107 - d_lll.y.xz * tmp111 - d_lll.y.yz * tmp125 + d_lll.y.xy * tmp103 - d_lll.y.xy * tmp106 + d_lll.x.zz * tmp116 + 2. * d_lll.x.yz * tmp111 - d_lll.x.zz * tmp113 + d_lll.x.yy * tmp106 - 2. * d_lll.x.yz * tmp107 + Z_l.z * gamma_uu.xz - d_lll.x.yy * tmp103 + Z_l.y * gamma_uu.xy + Z_l.x * gamma_uu.xx);\
	(resultFlux)->Z_l.x = alpha * (tmp1 + tmp2 + tmp6 + tmp4 + tmp7 - Theta);\
	(resultFlux)->Z_l.y = -alpha * (K_ll.yz * gamma_uu.xz + K_ll.yy * gamma_uu.xy + K_ll.xy * gamma_uu.xx);\
	(resultFlux)->Z_l.z = -alpha * (K_ll.zz * gamma_uu.xz + K_ll.yz * gamma_uu.xy + K_ll.xz * gamma_uu.xx);	\
	/* END CUT */\
<? end ?>\
<? if true then ?>\
	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);\
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);\
	real3 const e_l = _3sym3_tr12(d_ull);\
	real3 const d_l = real3x3x3_tr23(d_llu);\
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);\
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);\
	real3 const Z_u = sym3_real3_mul(gamma_uu, Z_l);\
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, K_ll);\
	real const tr_K = real3x3_trace(K_ul);\
	sym3 const dHat_t_ll = sym3_zero;\
<? 	if eqn.useShift == "2005 Bona / 2008 Yano" then ?>\
	real3 const a_u = sym3_real3_mul(gamma_uu, a_l);\
	sym3 const b_ll = sym3_real3x3_to_sym3_mul(gamma_ll, b_ul);\
<? 	end ?>\
	/* BEGIN CUT from symmath/tests/output/Z4.html */\
	{\
		(resultFlux)->alpha = 0.;\
		(resultFlux)->gamma_ll.xx = 0.;\
		(resultFlux)->gamma_ll.xy = 0.;\
		(resultFlux)->gamma_ll.xz = 0.;\
		(resultFlux)->gamma_ll.yy = 0.;\
		(resultFlux)->gamma_ll.yz = 0.;\
		(resultFlux)->gamma_ll.zz = 0.;\
		(resultFlux)->a_l.x = f_alpha * (tr_K + -2. * Theta);\
		(resultFlux)->a_l.y = 0.;\
		(resultFlux)->a_l.z = 0.;\
		(resultFlux)->dDelta_lll.x.xx = dHat_t_ll.xx + K_ll.xx * alpha;\
		(resultFlux)->dDelta_lll.x.xy = dHat_t_ll.xy + K_ll.xy * alpha;\
		(resultFlux)->dDelta_lll.x.xz = dHat_t_ll.xz + K_ll.xz * alpha;\
		(resultFlux)->dDelta_lll.x.yy = dHat_t_ll.yy + K_ll.yy * alpha;\
		(resultFlux)->dDelta_lll.x.yz = dHat_t_ll.yz + K_ll.yz * alpha;\
		(resultFlux)->dDelta_lll.x.zz = dHat_t_ll.zz + K_ll.zz * alpha;\
		(resultFlux)->dDelta_lll.y.xx = 0.;\
		(resultFlux)->dDelta_lll.y.xy = 0.;\
		(resultFlux)->dDelta_lll.y.xz = 0.;\
		(resultFlux)->dDelta_lll.y.yy = 0.;\
		(resultFlux)->dDelta_lll.y.yz = 0.;\
		(resultFlux)->dDelta_lll.y.zz = 0.;\
		(resultFlux)->dDelta_lll.z.xx = 0.;\
		(resultFlux)->dDelta_lll.z.xy = 0.;\
		(resultFlux)->dDelta_lll.z.xz = 0.;\
		(resultFlux)->dDelta_lll.z.yy = 0.;\
		(resultFlux)->dDelta_lll.z.yz = 0.;\
		(resultFlux)->dDelta_lll.z.zz = 0.;\
		(resultFlux)->K_ll.xx = alpha * (a_l.x + d_l.x + d_ull.x.xx + -e_l.x + -2. * Z_l.x + -d_llu.x.x.x);\
		(resultFlux)->K_ll.xy = (alpha * (a_l.y + d_l.y + -e_l.y + -d_llu.x.x.y + -d_llu.y.x.x + 2. * d_ull.x.xy + -2. * Z_l.y)) / 2.;\
		(resultFlux)->K_ll.xz = (alpha * (a_l.z + d_l.z + -e_l.z + -d_llu.x.x.z + -d_llu.z.x.x + 2. * d_ull.x.xz + -2. * Z_l.z)) / 2.;\
		(resultFlux)->K_ll.yy = alpha * (d_ull.x.yy + -d_llu.y.x.y);\
		(resultFlux)->K_ll.yz = (alpha * (-d_llu.y.x.z + 2. * d_ull.x.yz + -d_llu.z.x.y)) / 2.;\
		(resultFlux)->K_ll.zz = alpha * (d_ull.x.zz + -d_llu.z.x.z);\
		(resultFlux)->Theta = alpha * (d_u.x + -Z_u.x + -e_u.x);\
		(resultFlux)->Z_l.x = alpha * (tr_K + -K_ul.x.x + -Theta);\
		(resultFlux)->Z_l.y = -K_ul.x.y * alpha;\
		(resultFlux)->Z_l.z = -K_ul.x.z * alpha;\
	}\
	<? if eqn.useShift == "2005 Bona / 2008 Yano" then ?>\
	{\
		real const tmp1 = alpha * alpha;\
		(resultFlux)->a_l.x += -beta_u.x * a_l.x;\
		(resultFlux)->a_l.y += -beta_u.x * a_l.y;\
		(resultFlux)->a_l.z += -beta_u.x * a_l.z;\
		(resultFlux)->dDelta_lll.x.xx += -b_ll.xx + -beta_u.x * d_lll.x.xx + beta_u.y * dDelta_lll.y.xx + -beta_u.y * d_lll.y.xx + -beta_u.z * d_lll.z.xx + beta_u.z * dDelta_lll.z.xx;\
		(resultFlux)->dDelta_lll.x.xy += (-b_ll.xy + -b_ll.xy + -2. * beta_u.x * d_lll.x.xy + 2. * beta_u.y * dDelta_lll.y.xy + -2. * beta_u.y * d_lll.y.xy + -2. * beta_u.z * d_lll.z.xy + 2. * beta_u.z * dDelta_lll.z.xy) / 2.;\
		(resultFlux)->dDelta_lll.x.xz += (-b_ll.xz + -b_ll.xz + -2. * beta_u.x * d_lll.x.xz + 2. * beta_u.y * dDelta_lll.y.xz + -2. * beta_u.y * d_lll.y.xz + -2. * beta_u.z * d_lll.z.xz + 2. * beta_u.z * dDelta_lll.z.xz) / 2.;\
		(resultFlux)->dDelta_lll.x.yy += -b_ll.yy + -beta_u.x * d_lll.x.yy + beta_u.y * dDelta_lll.y.yy + -beta_u.y * d_lll.y.yy + -beta_u.z * d_lll.z.yy + beta_u.z * dDelta_lll.z.yy;\
		(resultFlux)->dDelta_lll.x.yz += (-b_ll.yz + -b_ll.yz + -2. * beta_u.x * d_lll.x.yz + 2. * beta_u.y * dDelta_lll.y.yz + -2. * beta_u.y * d_lll.y.yz + -2. * beta_u.z * d_lll.z.yz + 2. * beta_u.z * dDelta_lll.z.yz) / 2.;\
		(resultFlux)->dDelta_lll.x.zz += -b_ll.zz + -beta_u.x * d_lll.x.zz + beta_u.y * dDelta_lll.y.zz + -beta_u.y * d_lll.y.zz + -beta_u.z * d_lll.z.zz + beta_u.z * dDelta_lll.z.zz;\
		(resultFlux)->dDelta_lll.y.xx += -beta_u.x * dDelta_lll.y.xx;\
		(resultFlux)->dDelta_lll.y.xy += -beta_u.x * dDelta_lll.y.xy;\
		(resultFlux)->dDelta_lll.y.xz += -beta_u.x * dDelta_lll.y.xz;\
		(resultFlux)->dDelta_lll.y.yy += -beta_u.x * dDelta_lll.y.yy;\
		(resultFlux)->dDelta_lll.y.yz += -beta_u.x * dDelta_lll.y.yz;\
		(resultFlux)->dDelta_lll.y.zz += -beta_u.x * dDelta_lll.y.zz;\
		(resultFlux)->dDelta_lll.z.xx += -beta_u.x * dDelta_lll.z.xx;\
		(resultFlux)->dDelta_lll.z.xy += -beta_u.x * dDelta_lll.z.xy;\
		(resultFlux)->dDelta_lll.z.xz += -beta_u.x * dDelta_lll.z.xz;\
		(resultFlux)->dDelta_lll.z.yy += -beta_u.x * dDelta_lll.z.yy;\
		(resultFlux)->dDelta_lll.z.yz += -beta_u.x * dDelta_lll.z.yz;\
		(resultFlux)->dDelta_lll.z.zz += -beta_u.x * dDelta_lll.z.zz;\
		(resultFlux)->K_ll.xx += -beta_u.x * K_ll.xx;\
		(resultFlux)->K_ll.xy += -beta_u.x * K_ll.xy;\
		(resultFlux)->K_ll.xz += -beta_u.x * K_ll.xz;\
		(resultFlux)->K_ll.yy += -beta_u.x * K_ll.yy;\
		(resultFlux)->K_ll.yz += -beta_u.x * K_ll.yz;\
		(resultFlux)->K_ll.zz += -beta_u.x * K_ll.zz;\
		(resultFlux)->Theta += -beta_u.x * Theta;\
		(resultFlux)->Z_l.x += -Z_l.x * beta_u.x;\
		(resultFlux)->Z_l.y += -Z_l.y * beta_u.x;\
		(resultFlux)->Z_l.z += -Z_l.z * beta_u.x;\
		(resultFlux)->beta_u.x = 0.;\
		(resultFlux)->beta_u.y = 0.;\
		(resultFlux)->beta_u.z = 0.;\
		(resultFlux)->b_ul.x.x = -beta_u.x * b_ul.x.x + -beta_u.y * b_ul.x.y + -beta_u.z * b_ul.x.z + a_u.x * tmp1 + -2. * e_u.x * tmp1 + d_u.x * tmp1;\
		(resultFlux)->b_ul.x.y = 0.;\
		(resultFlux)->b_ul.x.z = 0.;\
		(resultFlux)->b_ul.y.x = -beta_u.x * b_ul.y.x + -beta_u.y * b_ul.y.y + -beta_u.z * b_ul.y.z + a_u.y * tmp1 + -2. * e_u.y * tmp1 + d_u.y * tmp1;\
		(resultFlux)->b_ul.y.y = 0.;\
		(resultFlux)->b_ul.y.z = 0.;\
		(resultFlux)->b_ul.z.x = -beta_u.x * b_ul.z.x + -beta_u.y * b_ul.z.y + -beta_u.z * b_ul.z.z + a_u.z * tmp1 + -2. * e_u.z * tmp1 + d_u.z * tmp1;\
		(resultFlux)->b_ul.z.y = 0.;\
		(resultFlux)->b_ul.z.z = 0.;\
	}\
	<? end ?>/* eqn.useShift == "2005 Bona / 2008 Yano" */\
	/* END CUT */\
<? end ?>\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
\
	(resultFlux)->Z_l = real3_rotateFrom((resultFlux)->Z_l, n_l);\
	(resultFlux)->a_l = real3_rotateFrom((resultFlux)->a_l, n_l);\
	(resultFlux)->gamma_ll = sym3_rotateFrom((resultFlux)->gamma_ll, n_l);\
	(resultFlux)->K_ll = sym3_rotateFrom((resultFlux)->K_ll, n_l);\
	(resultFlux)->dDelta_lll = _3sym3_rotateFrom((resultFlux)->dDelta_lll, n_l);\
<? if eqn.useShift ~= "none" then --\
?>	(resultFlux)->beta_u = real3_rotateFrom((resultFlux)->beta_u, n_l);\
<? end --\
?>\
<?	if eqn.useShift == "2005 Bona / 2008 Yano" then --\
?>	(resultFlux)->b_ul = real3x3_rotateFrom((resultFlux)->b_ul, n_l);\
<? end --\
?>\
<? else ?>\
\
	if (false) {}\
	<? for side=0,solver.dim-1 do ?>\
	else if (n.side == <?=side?>) {\
		(resultFlux)->Z_l = real3_swap<?=side?>((resultFlux)->Z_l);\
		(resultFlux)->a_l = real3_swap<?=side?>((resultFlux)->a_l);\
		(resultFlux)->gamma_ll = sym3_swap<?=side?>((resultFlux)->gamma_ll);\
		(resultFlux)->K_ll = sym3_swap<?=side?>((resultFlux)->K_ll);\
		(resultFlux)->dDelta_lll = _3sym3_swap<?=side?>((resultFlux)->dDelta_lll);\
<? if eqn.useShift ~= "none" then --\
?>		(resultFlux)->beta_u = real3_swap<?=side?>((resultFlux)->beta_u);\
<? end --\
?>\
<?	if eqn.useShift == "2005 Bona / 2008 Yano" then --\
?>		(resultFlux)->b_ul = real3x3_swap<?=side?>((resultFlux)->b_ul);\
<? end --\
?>	}\
	<? end ?>\
	else {\
		alpha = 0./0.;\
		Theta = 0./0.;\
<? for k,xk in ipairs(xNames) do --\
?>		a_l.<?=xk?> = 0./0.;\
		Z_l.<?=xk?> = 0./0.;\
<? end --\
?>\
<? for ij,xij in ipairs(symNames) do --\
?>		K_ll.<?=xij?> = 0./0.;\
		gamma_ll.<?=xij?> = 0./0.;\
<? end --\
?>\
<? for k,xk in ipairs(xNames) do --\
	for ij,xij in ipairs(symNames) do --\
?>		dDelta_lll.<?=xk?>.<?=xij?> = 0./0.;\
<?	end --\
end --\
?>	}\
\
<? end ?>\
\
<? --\
if eqn.useShift ~= "none" then --\
?>	/* beta^i_,t = 0 + source terms */\
	(resultFlux)->beta_u = real3_zero;\
<?	if eqn.useShift == "MinimalDistortionElliptic" --\
	or eqn.useShift == "MinimalDistortionEllipticEvolve" --\
	then --\
?>	(resultFlux)->betaLap_u = real3_zero;\
<?	elseif eqn.useShift == "2005 Bona / 2008 Yano" then --\
?>	(resultFlux)->b_ul = real3x3_zero;\
<?	end --\
end --\
?>\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=solver_t?> <?=eigen_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>
// used by hll, roe, weno, plm ... anything that uses eigenvalues or eigenvector transforms

//used for interface eigen basis
#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	(resultEig)->alpha = .5 * ((UL)->alpha + (UR)->alpha);\
	(resultEig)->alpha_sqrt_f = sqrt(calc_f_alphaSq((resultEig)->alpha));\
\
	sym3 const avg_gamma = sym3_real_mul(sym3_add((UL)->gamma_ll, (UR)->gamma_ll), .5);\
	real const det_avg_gamma = sym3_det(avg_gamma);\
	(resultEig)->gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
/*  I'm using .side for holonomic(coordinate) and anholonomic(orthonormal) */\
/* but for cartesian vector componets there is no .side, just .n, which is covariant iirc */\
/* and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes) */\
/* so here I'm going to just wing it */\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, (resultEig)->gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = (resultEig)->gamma_uu.xx;\
	} else if (n.side == 1) {\
		gammaUnn = (resultEig)->gamma_uu.yy;\
	} else if (n.side == 2) {\
		gammaUnn = (resultEig)->gamma_uu.zz;\
	}\
<? end ?>\
\
	(resultEig)->sqrt_gammaUnn = sqrt(gammaUnn);\
\
<? if eqn.useShift ~= "none" then ?>\
	(resultEig)->beta_u = real3_real_mul(real3_add((UL)->beta_u, (UR)->beta_u), .5);\
<? end ?>\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: <?=range_t?> <?=normal_t?> cons_t <?=initCond_codeprefix?>
// not used anymore, replaced in calcDT by eqn:consMinWaveCode/eqn:consMaxWaveCode eigenvalue inlining

#define <?=calcCellMinMaxEigenvalues?>(\
	/*<?=range_t?> * const */result,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const det_gamma = sym3_det((U)->gamma_ll);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
	sym3 const gamma_uu = sym3_inv((U)->gamma_ll, det_gamma);\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = ((U)->gamma_ll.yy * (U)->gamma_ll.zz - (U)->gamma_ll.yz * (U)->gamma_ll.yz) / det_gamma;\
	} else if (n.side == 1) {\
		gammaUnn = ((U)->gamma_ll.xx * (U)->gamma_ll.zz - (U)->gamma_ll.xz * (U)->gamma_ll.xz) / det_gamma;\
	} else if (n.side == 2) {\
		gammaUnn = ((U)->gamma_ll.xx * (U)->gamma_ll.yy - (U)->gamma_ll.xy * (U)->gamma_ll.xy) / det_gamma;\
	}\
<? end ?>\
	real const sqrt_gammaUnn = sqrt(gammaUnn);\
	real const lambdaLight = (U)->alpha * sqrt_gammaUnn;\
\
	real const f_alphaSq = calc_f_alphaSq((U)->alpha);\
	real const lambdaGauge = sqrt(f_alphaSq) * sqrt_gammaUnn;\
\
	real lambdaMax = max(lambdaGauge, lambdaLight);\
	real lambdaMin = -lambdaMin;\
\
	<? if eqn.useShift ~= "none" then ?>\
	lambdaMin -= normal_vecDotN1(n, (U)->beta_u);\
	lambdaMax -= normal_vecDotN1(n, (U)->beta_u);\
	<? end ?>\
\
	(result)->min = lambdaMin;\
	(result)->max = lambdaMax;\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>

//used by PLM, and by the default <?=fluxFromCons?> (used by hll, or roe when roeUseFluxFromCons is set)
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	(resultEig)->alpha = (U)->alpha;\
	(resultEig)->alpha_sqrt_f = sqrt(calc_f_alphaSq((U)->alpha));\
\
	real const det_gamma = sym3_det((U)->gamma_ll);\
	(resultEig)->gamma_uu = sym3_inv((U)->gamma_ll, det_gamma);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
/*  I'm using .side for holonomic(coordinate) and anholonomic(orthonormal) */\
/* but for cartesian vector componets there is no .side, just .n, which is covariant iirc */\
/* and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes) */\
/* so here I'm going to just wing it */\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, (resultEig)->gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = (resultEig)->gamma_uu.xx;\
	} else if (n.side == 1) {\
		gammaUnn = (resultEig)->gamma_uu.yy;\
	} else if (n.side == 2) {\
		gammaUnn = (resultEig)->gamma_uu.zz;\
	}\
<? end ?>\
\
	(resultEig)->sqrt_gammaUnn = sqrt(gammaUnn);\
\
	<? if eqn.useShift ~= "none" then ?>\
	(resultEig)->beta_u = (U)->beta_u;\
	<? end ?>\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> rotate <?=initCond_codeprefix?>
// used by roe, weno, some plm

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */inputU,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	for (int i = 0; i < numWaves; ++i) {\
		(result)->ptr[i] = 0./0.;\
	}\
\
#error "this is broke for now. use a solver like hll, not anything that requires eigensystems like roe"\
<? if false then ?>	/* don't enable this.  it's made for > waves than I'm using, so it will cause buffer corruption */\
\
	real3 const a_l = real3_swap((inputU)->a_l, n.side);							/* 0-2 */\
	_3sym3 const d_lll = _3sym3_swap((inputU)->d_lll, n.side);						/* 3-20 ... .x = 3-8, .y = 9-14, .z = 15-20 */\
	sym3 const K_ll = sym3_swap((inputU)->K_ll, n.side);							/* 21-26 */\
	real const Theta = (inputU)->Theta;												/* 27 */\
	real3 const Z_l = real3_swap((inputU)->Z_l, n.side);							/* 28-30 */\
\
	sym3 const gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 const gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
\
	real const sqrt_gammaUUxx = sqrt(gamma_uu.xx);\
	real const gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;\
	real const gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;\
\
	real const sqrt_f = (eig)->alpha_sqrt_f / (eig)->alpha;\
	real const f = sqrt_f * sqrt_f;\
	/*  from 2004 Bona et al, "A symmetry breaking..." eqn A.20 */\
	/*  mind you the 'm' in that form is for alpha_,t = -alpha^2 (f K - m Theta) */\
	/*  while in more modern Z4 papers it is alpha_,t = -alpha^2 f (K - m Theta) */\
	/*  the difference means that this eigendecomposition is probably incorrect. */\
	real const lambda_1 = (2. - solver->m) / (f - 1);\
	real const lambda_2 = (2. * f - solver->m) / (f - 1);\
\
	(result)->ptr[0] = (((lambda_2 * gamma_uu.xx * Z_l.x) - (lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (lambda_2 * gamma_uu.xy * Z_l.y) + ((((2. * lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (lambda_2 * gamma_uu.xz * Z_l.z) + ((lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)) + (sqrt_f * lambda_1 * Theta) + (sqrt_f * gamma_uu.xx * K_ll.xx) + (2. * sqrt_f * gamma_uu.xy * K_ll.xy) + (2. * sqrt_f * gamma_uu.xz * K_ll.xz) + (sqrt_f * gamma_uu.yy * K_ll.yy) + (2. * sqrt_f * gamma_uu.yz * K_ll.yz) + ((((sqrt_f * gamma_uu.zz * K_ll.zz) - (gamma_uu.xx * a_l.x)) - (gamma_uu.xy * a_l.y)) - (gamma_uu.xz * a_l.z)));\
	(result)->ptr[1] = (Theta + ((gamma_uu.xx * Z_l.x) - (gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (gamma_uu.xy * Z_l.y) + ((((2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (gamma_uu.xz * Z_l.z) + ((gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)));\
	(result)->ptr[2] = ((-(((a_l.y - (2. * K_ll.xy)) - (2. * Z_l.y)) + ((gamma_uu.xx * d_lll.y.xx) - (2. * gamma_uu.xy * d_lll.x.yy)) + ((2. * gamma_uu.xy * d_lll.y.xy) - (2. * gamma_uu.xz * d_lll.x.yz)) + (((2. * gamma_uu.xz * d_lll.y.xz) - (gamma_uu.yy * d_lll.y.yy)) - (2. * gamma_uu.yz * d_lll.z.yy)) + ((gamma_uu.zz * d_lll.y.zz) - (2. * gamma_uu.zz * d_lll.z.yz)))) / 2.);\
	(result)->ptr[3] = ((-(((a_l.z - (2. * K_ll.xz)) - (2. * Z_l.z)) + ((gamma_uu.xx * d_lll.z.xx) - (2. * gamma_uu.xy * d_lll.x.yz)) + ((2. * gamma_uu.xy * d_lll.z.xy) - (2. * gamma_uu.xz * d_lll.x.zz)) + ((2. * gamma_uu.xz * d_lll.z.xz) - (2. * gamma_uu.yy * d_lll.y.yz)) + (((gamma_uu.yy * d_lll.z.yy) - (2. * gamma_uu.yz * d_lll.y.zz)) - (gamma_uu.zz * d_lll.z.zz)))) / 2.);\
	(result)->ptr[4] = (((K_ll.yy - (gamma_uu.xx * d_lll.x.yy)) - (gamma_uu.xy * d_lll.y.yy)) - (gamma_uu.xz * d_lll.z.yy));\
	(result)->ptr[5] = (((K_ll.yz - (gamma_uu.xx * d_lll.x.yz)) - (gamma_uu.xy * d_lll.y.yz)) - (gamma_uu.xz * d_lll.z.yz));\
	(result)->ptr[6] = (((K_ll.zz - (gamma_uu.xx * d_lll.x.zz)) - (gamma_uu.xy * d_lll.y.zz)) - (gamma_uu.xz * d_lll.z.zz));\
	(result)->ptr[7] = a_l.y;\
	(result)->ptr[8] = a_l.z;\
	(result)->ptr[9] = d_lll.y.xx;\
	(result)->ptr[10] = d_lll.y.xy;\
	(result)->ptr[11] = d_lll.y.xz;\
	(result)->ptr[12] = d_lll.y.yy;\
	(result)->ptr[13] = d_lll.y.yz;\
	(result)->ptr[14] = d_lll.y.zz;\
	(result)->ptr[15] = d_lll.z.xx;\
	(result)->ptr[16] = d_lll.z.xy;\
	(result)->ptr[17] = d_lll.z.xz;\
	(result)->ptr[18] = d_lll.z.yy;\
	(result)->ptr[19] = d_lll.z.yz;\
	(result)->ptr[20] = d_lll.z.zz;\
	(result)->ptr[21] = (a_l.x + ((Z_l.x - (f * gamma_uu.xx * d_lll.x.xx)) - (gamma_uu.xy * d_lll.x.xy)) + (((gamma_uu.xy * d_lll.y.xx) - (2. * gamma_uu.xy * f * d_lll.x.xy)) - (gamma_uu.xz * d_lll.x.xz)) + (((gamma_uu.xz * d_lll.z.xx) - (2. * gamma_uu.xz * f * d_lll.x.xz)) - (gamma_uu.yy * d_lll.x.yy)) + (((gamma_uu.yy * d_lll.y.xy) - (gamma_uu.yy * f * d_lll.x.yy)) - (2. * gamma_uu.yz * d_lll.x.yz)) + (gamma_uu.yz * d_lll.y.xz) + (((gamma_uu.yz * d_lll.z.xy) - (2. * gamma_uu.yz * f * d_lll.x.yz)) - (gamma_uu.zz * d_lll.x.zz)) + ((gamma_uu.zz * d_lll.z.xz) - (gamma_uu.zz * f * d_lll.x.zz)));\
	(result)->ptr[22] = (a_l.y + (Z_l.y - (f * gamma_uu.yy * d_lll.y.yy)) + (((gamma_uu.xx * d_lll.x.xy) - (gamma_uu.xx * d_lll.y.xx)) - (gamma_uu.xx * f * d_lll.y.xx)) + (((gamma_uu.xy * d_lll.x.yy) - (gamma_uu.xy * d_lll.y.xy)) - (2. * gamma_uu.xy * f * d_lll.y.xy)) + ((gamma_uu.xz * d_lll.x.yz) - (2. * gamma_uu.xz * d_lll.y.xz)) + (((gamma_uu.xz * d_lll.z.xy) - (2. * gamma_uu.xz * f * d_lll.y.xz)) - (gamma_uu.yz * d_lll.y.yz)) + (((gamma_uu.yz * d_lll.z.yy) - (2. * gamma_uu.yz * f * d_lll.y.yz)) - (gamma_uu.zz * d_lll.y.zz)) + ((gamma_uu.zz * d_lll.z.yz) - (gamma_uu.zz * f * d_lll.y.zz)));\
	(result)->ptr[23] = (a_l.z + (Z_l.z - (f * gamma_uu.zz * d_lll.z.zz)) + (((gamma_uu.xx * d_lll.x.xz) - (gamma_uu.xx * d_lll.z.xx)) - (gamma_uu.xx * f * d_lll.z.xx)) + (gamma_uu.xy * d_lll.x.yz) + (((gamma_uu.xy * d_lll.y.xz) - (2. * gamma_uu.xy * d_lll.z.xy)) - (2. * gamma_uu.xy * f * d_lll.z.xy)) + (((gamma_uu.xz * d_lll.x.zz) - (gamma_uu.xz * d_lll.z.xz)) - (2. * gamma_uu.xz * f * d_lll.z.xz)) + (((gamma_uu.yy * d_lll.y.yz) - (gamma_uu.yy * d_lll.z.yy)) - (gamma_uu.yy * f * d_lll.z.yy)) + (((gamma_uu.yz * d_lll.y.zz) - (gamma_uu.yz * d_lll.z.yz)) - (2. * gamma_uu.yz * f * d_lll.z.yz)));\
	(result)->ptr[24] = ((a_l.y + ((2. * K_ll.xy) - (2. * Z_l.y)) + ((gamma_uu.xx * d_lll.y.xx) - (2. * gamma_uu.xy * d_lll.x.yy)) + ((2. * gamma_uu.xy * d_lll.y.xy) - (2. * gamma_uu.xz * d_lll.x.yz)) + (((2. * gamma_uu.xz * d_lll.y.xz) - (gamma_uu.yy * d_lll.y.yy)) - (2. * gamma_uu.yz * d_lll.z.yy)) + ((gamma_uu.zz * d_lll.y.zz) - (2. * gamma_uu.zz * d_lll.z.yz))) / 2.);\
	(result)->ptr[25] = ((a_l.z + ((2. * K_ll.xz) - (2. * Z_l.z)) + ((gamma_uu.xx * d_lll.z.xx) - (2. * gamma_uu.xy * d_lll.x.yz)) + ((2. * gamma_uu.xy * d_lll.z.xy) - (2. * gamma_uu.xz * d_lll.x.zz)) + ((2. * gamma_uu.xz * d_lll.z.xz) - (2. * gamma_uu.yy * d_lll.y.yz)) + (((gamma_uu.yy * d_lll.z.yy) - (2. * gamma_uu.yz * d_lll.y.zz)) - (gamma_uu.zz * d_lll.z.zz))) / 2.);\
	(result)->ptr[26] = (K_ll.yy + (gamma_uu.xx * d_lll.x.yy) + (gamma_uu.xy * d_lll.y.yy) + (gamma_uu.xz * d_lll.z.yy));\
	(result)->ptr[27] = (K_ll.yz + (gamma_uu.xx * d_lll.x.yz) + (gamma_uu.xy * d_lll.y.yz) + (gamma_uu.xz * d_lll.z.yz));\
	(result)->ptr[28] = (K_ll.zz + (gamma_uu.xx * d_lll.x.zz) + (gamma_uu.xy * d_lll.y.zz) + (gamma_uu.xz * d_lll.z.zz));\
	(result)->ptr[29] = ((Theta - (gamma_uu.xx * Z_l.x)) + ((gamma_uu.xx * gamma_uu.yy * d_lll.x.yy) - (gamma_uu.xx * gamma_uu.yy * d_lll.y.xy)) + (((2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz) - (gamma_uu.xx * gamma_uu.yz * d_lll.y.xz)) - (gamma_uu.xx * gamma_uu.yz * d_lll.z.xy)) + (((gamma_uu.xx * gamma_uu.zz * d_lll.x.zz) - (gamma_uu.xx * gamma_uu.zz * d_lll.z.xz)) - (gamma_uu.xy * gamma_uu.xy * d_lll.x.yy)) + (((gamma_uu.xy * gamma_uu.xy * d_lll.y.xy) - (gamma_uu.xy * Z_l.y)) - (2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz)) + (gamma_uu.xy * gamma_uu.xz * d_lll.y.xz) + (gamma_uu.xy * gamma_uu.xz * d_lll.z.xy) + ((gamma_uu.xy * gamma_uu.yz * d_lll.y.yz) - (gamma_uu.xy * gamma_uu.yz * d_lll.z.yy)) + (((gamma_uu.xy * gamma_uu.zz * d_lll.y.zz) - (gamma_uu.xy * gamma_uu.zz * d_lll.z.yz)) - (gamma_uu.xz * gamma_uu.xz * d_lll.x.zz)) + (((gamma_uu.xz * gamma_uu.xz * d_lll.z.xz) - (gamma_uu.xz * Z_l.z)) - (gamma_uu.xz * gamma_uu.yy * d_lll.y.yz)) + ((gamma_uu.xz * gamma_uu.yy * d_lll.z.yy) - (gamma_uu.xz * gamma_uu.yz * d_lll.y.zz)) + (gamma_uu.xz * gamma_uu.yz * d_lll.z.yz));\
	(result)->ptr[30] = (-(((lambda_2 * gamma_uu.xx * Z_l.x) - (lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (lambda_2 * gamma_uu.xy * Z_l.y) + ((((2. * lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (lambda_2 * gamma_uu.xz * Z_l.z) + ((lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((((((((((((lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)) - (sqrt_f * lambda_1 * Theta)) - (sqrt_f * gamma_uu.xx * K_ll.xx)) - (2. * sqrt_f * gamma_uu.xy * K_ll.xy)) - (2. * sqrt_f * gamma_uu.xz * K_ll.xz)) - (sqrt_f * gamma_uu.yy * K_ll.yy)) - (2. * sqrt_f * gamma_uu.yz * K_ll.yz)) - (sqrt_f * gamma_uu.zz * K_ll.zz)) - (gamma_uu.xx * a_l.x)) - (gamma_uu.xy * a_l.y)) - (gamma_uu.xz * a_l.z))));\
\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> rotate <?=initCond_codeprefix?>

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */input,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	for (int j = 0; j < numStates; ++j) {\
		(result)->ptr[j] = 0./0.;\
	}\
\
#error until you regenerate these left/right transforms, don't use roe with z4\
<? if false then ?>\
	sym3 const gamma_ll = sym3_swap((eig)->gamma_ll, n.side);\
	sym3 const gamma_uu = sym3_swap((eig)->gamma_uu, n.side);\
\
	real const sqrt_gammaUUxx = sqrt(gamma_uu.xx);\
	real const gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;\
	real const gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;\
\
	real const sqrt_f = (eig)->alpha_sqrt_f / (eig)->alpha;\
	real const f = sqrt_f * sqrt_f;\
	real const fSq = f * f;\
	/*  from 2004 Bona et al, "A symmetry breaking..." eqn A.20 */\
	real const lambda_1 = (2. - solver->m) / (f - 1);\
	real const lambda_2 = (2. * f - solver->m) / (f - 1);\
\
	(result)->ptr[0] = ((((input.ptr[0] - input.ptr[30]) - (lambda_2 * input.ptr[1])) + (lambda_2 * input.ptr[29]) + (2. * gamma_uu.xy * input.ptr[7]) + (2. * gamma_uu.xz * input.ptr[8])) / (-(2. * gamma_uu.xx)));\
	(result)->ptr[1] = input.ptr[7];\
	(result)->ptr[2] = input.ptr[8];\
	(result)->ptr[3] = ((((2. * input.ptr[0]) - (2. * input.ptr[1])) + (4. * input.ptr[21] * gamma_uu.xx) + (((2. * input.ptr[29]) - (2. * input.ptr[30])) - (2. * lambda_2 * input.ptr[1])) + (((2. * lambda_2 * input.ptr[29]) - (4. * gamma_uu.xy * input.ptr[2] * f)) - (12. * gamma_uu.xy * input.ptr[7] * f)) + (8. * gamma_uu.xy * input.ptr[9] * gamma_uu.xx * f) + (8. * gamma_uu.xy * gamma_uu.xy * input.ptr[10] * f) + (4. * gamma_uu.xy * input.ptr[22]) + (4. * gamma_uu.xy * input.ptr[24] * f) + (8. * gamma_uu.xy * fSq * input.ptr[9] * gamma_uu.xx) + (16. * gamma_uu.xy * gamma_uu.xy * fSq * input.ptr[10]) + (8. * gamma_uu.xy * f * input.ptr[22]) + (8. * gamma_uu.xy * gamma_uu.xz * input.ptr[11] * f) + (8. * gamma_uu.xy * gamma_uu.xz * input.ptr[16] * f) + (16. * gamma_uu.xy * gamma_uu.xz * fSq * input.ptr[11]) + (16. * gamma_uu.xy * gamma_uu.xz * fSq * input.ptr[16]) + (8. * gamma_uu.xy * gamma_uu.yz * input.ptr[13] * f) + (((16. * gamma_uu.xy * gamma_uu.yz * fSq * input.ptr[13]) - (4. * gamma_uu.xz * input.ptr[3] * f)) - (12. * gamma_uu.xz * input.ptr[8] * f)) + (8. * gamma_uu.xz * input.ptr[15] * gamma_uu.xx * f) + (8. * gamma_uu.xz * gamma_uu.xz * input.ptr[17] * f) + (4. * gamma_uu.xz * input.ptr[23]) + (4. * gamma_uu.xz * input.ptr[25] * f) + (8. * gamma_uu.xz * fSq * input.ptr[15] * gamma_uu.xx) + (16. * gamma_uu.xz * gamma_uu.xz * fSq * input.ptr[17]) + (8. * gamma_uu.xz * f * input.ptr[23]) + (8. * gamma_uu.xz * gamma_uu.yz * input.ptr[19] * f) + ((16. * gamma_uu.xz * gamma_uu.yz * fSq * input.ptr[19]) - (2. * gamma_uu.yy * input.ptr[4] * f)) + (2. * gamma_uu.yy * input.ptr[26] * f) + (4. * gamma_uu.yy * gamma_uu.xy * input.ptr[12] * f) + (8. * gamma_uu.yy * gamma_uu.xy * fSq * input.ptr[12]) + (4. * gamma_uu.yy * gamma_uu.xz * input.ptr[18] * f) + ((8. * gamma_uu.yy * gamma_uu.xz * fSq * input.ptr[18]) - (4. * gamma_uu.yz * input.ptr[5] * f)) + ((4. * gamma_uu.yz * input.ptr[27] * f) - (2. * gamma_uu.zz * input.ptr[6] * f)) + (2. * gamma_uu.zz * input.ptr[28] * f) + (4. * gamma_uu.zz * gamma_uu.xy * input.ptr[14] * f) + (8. * gamma_uu.zz * gamma_uu.xy * fSq * input.ptr[14]) + (4. * gamma_uu.zz * gamma_uu.xz * input.ptr[20] * f) + (8. * gamma_uu.zz * gamma_uu.xz * fSq * input.ptr[20])) / (-(4. * gammaUUxxSq * f)));\
	(result)->ptr[4] = ((-(input.ptr[2] + (((((((3. * input.ptr[7]) - (input.ptr[9] * gamma_uu.xx)) - (2. * input.ptr[22])) - input.ptr[24]) - (2. * f * input.ptr[9] * gamma_uu.xx)) - (4. * gamma_uu.xy * f * input.ptr[10])) - (2. * gamma_uu.xz * input.ptr[11])) + ((((((((2. * gamma_uu.xz * input.ptr[16]) - (4. * gamma_uu.xz * f * input.ptr[11])) - (gamma_uu.yy * input.ptr[12])) - (2. * gamma_uu.yy * f * input.ptr[12])) - (2. * gamma_uu.yz * input.ptr[13])) - (4. * gamma_uu.yz * f * input.ptr[13])) - (gamma_uu.zz * input.ptr[14])) - (2. * gamma_uu.zz * f * input.ptr[14])))) / (2. * gamma_uu.xx));\
	(result)->ptr[5] = ((-(input.ptr[3] + (((((3. * input.ptr[8]) - (input.ptr[15] * gamma_uu.xx)) - (2. * input.ptr[23])) - input.ptr[25]) - (2. * f * input.ptr[15] * gamma_uu.xx)) + ((((((((((2. * gamma_uu.xy * input.ptr[11]) - (2. * gamma_uu.xy * input.ptr[16])) - (4. * gamma_uu.xy * f * input.ptr[16])) - (4. * gamma_uu.xz * f * input.ptr[17])) - (gamma_uu.yy * input.ptr[18])) - (2. * gamma_uu.yy * f * input.ptr[18])) - (2. * gamma_uu.yz * input.ptr[19])) - (4. * gamma_uu.yz * f * input.ptr[19])) - (gamma_uu.zz * input.ptr[20])) - (2. * gamma_uu.zz * f * input.ptr[20])))) / (2. * gamma_uu.xx));\
	(result)->ptr[6] = ((-((input.ptr[4] - input.ptr[26]) + (2. * gamma_uu.xy * input.ptr[12]) + (2. * gamma_uu.xz * input.ptr[18]))) / (2. * gamma_uu.xx));\
	(result)->ptr[7] = ((-((input.ptr[5] - input.ptr[27]) + (2. * gamma_uu.xy * input.ptr[13]) + (2. * gamma_uu.xz * input.ptr[19]))) / (2. * gamma_uu.xx));\
	(result)->ptr[8] = ((-((input.ptr[6] - input.ptr[28]) + (2. * gamma_uu.xy * input.ptr[14]) + (2. * gamma_uu.xz * input.ptr[20]))) / (2. * gamma_uu.xx));\
	(result)->ptr[9] = input.ptr[9];\
	(result)->ptr[10] = input.ptr[10];\
	(result)->ptr[11] = input.ptr[11];\
	(result)->ptr[12] = input.ptr[12];\
	(result)->ptr[13] = input.ptr[13];\
	(result)->ptr[14] = input.ptr[14];\
	(result)->ptr[15] = input.ptr[15];\
	(result)->ptr[16] = input.ptr[16];\
	(result)->ptr[17] = input.ptr[17];\
	(result)->ptr[18] = input.ptr[18];\
	(result)->ptr[19] = input.ptr[19];\
	(result)->ptr[20] = input.ptr[20];\
	(result)->ptr[21] = ((input.ptr[0] + ((((((((((((input.ptr[30] - (lambda_1 * input.ptr[1] * sqrt_f)) - (lambda_1 * input.ptr[29] * sqrt_f)) - (2. * gamma_uu.xy * input.ptr[2] * sqrt_f)) - (2. * gamma_uu.xy * input.ptr[24] * sqrt_f)) - (2. * gamma_uu.xz * input.ptr[3] * sqrt_f)) - (2. * gamma_uu.xz * input.ptr[25] * sqrt_f)) - (gamma_uu.yy * input.ptr[4] * sqrt_f)) - (gamma_uu.yy * input.ptr[26] * sqrt_f)) - (2. * gamma_uu.yz * input.ptr[5] * sqrt_f)) - (2. * gamma_uu.yz * input.ptr[27] * sqrt_f)) - (gamma_uu.zz * input.ptr[6] * sqrt_f)) - (gamma_uu.zz * input.ptr[28] * sqrt_f))) / (2. * sqrt_f * gamma_uu.xx));\
	(result)->ptr[22] = ((input.ptr[2] + input.ptr[24]) / 2.);\
	(result)->ptr[23] = ((input.ptr[3] + input.ptr[25]) / 2.);\
	(result)->ptr[24] = ((input.ptr[4] + input.ptr[26]) / 2.);\
	(result)->ptr[25] = ((input.ptr[5] + input.ptr[27]) / 2.);\
	(result)->ptr[26] = ((input.ptr[6] + input.ptr[28]) / 2.);\
	(result)->ptr[27] = ((input.ptr[1] + input.ptr[29]) / 2.);\
	(result)->ptr[28] = ((((((input.ptr[1] - input.ptr[29]) - (gamma_uu.xy * input.ptr[2])) - (gamma_uu.xy * input.ptr[7])) - (gamma_uu.xy * input.ptr[9] * gamma_uu.xx)) + (((((gamma_uu.xy * input.ptr[24]) - (2. * gamma_uu.xy * gamma_uu.yz * input.ptr[13])) - (gamma_uu.xz * input.ptr[3])) - (gamma_uu.xz * input.ptr[8])) - (gamma_uu.xz * input.ptr[15] * gamma_uu.xx)) + (((gamma_uu.xz * input.ptr[25]) - (gamma_uu.yy * input.ptr[4])) - (2. * gamma_uu.yy * input.ptr[10] * gamma_uu.xx)) + ((((((gamma_uu.yy * input.ptr[26]) - (gamma_uu.yy * gamma_uu.xy * input.ptr[12])) - (gamma_uu.yy * gamma_uu.xz * input.ptr[18])) - (2. * gamma_uu.yz * input.ptr[5])) - (2. * gamma_uu.yz * input.ptr[11] * gamma_uu.xx)) - (2. * gamma_uu.yz * input.ptr[16] * gamma_uu.xx)) + ((((2. * gamma_uu.yz * input.ptr[27]) - (2. * gamma_uu.yz * gamma_uu.xz * input.ptr[19])) - (gamma_uu.zz * input.ptr[6])) - (2. * gamma_uu.zz * input.ptr[17] * gamma_uu.xx)) + (((gamma_uu.zz * input.ptr[28]) - (gamma_uu.zz * gamma_uu.xy * input.ptr[14])) - (gamma_uu.zz * gamma_uu.xz * input.ptr[20]))) / (2. * gamma_uu.xx));\
	(result)->ptr[29] = (((input.ptr[2] * gamma_uu.xx) + ((input.ptr[7] * gamma_uu.xx) - (input.ptr[24] * gamma_uu.xx)) + (((gammaUUxxSq * input.ptr[9]) - (gamma_uu.xx * gamma_uu.yy * input.ptr[12])) - (2. * gamma_uu.xx * gamma_uu.yz * input.ptr[18])) + (gamma_uu.xy * input.ptr[4]) + (2. * gamma_uu.xy * input.ptr[10] * gamma_uu.xx) + ((2. * gamma_uu.xy * gamma_uu.xy * input.ptr[12]) - (gamma_uu.xy * input.ptr[26])) + (2. * gamma_uu.xy * gamma_uu.xz * input.ptr[13]) + (gamma_uu.xz * input.ptr[5]) + (2. * gamma_uu.xz * input.ptr[11] * gamma_uu.xx) + ((2. * gamma_uu.xz * gamma_uu.xz * input.ptr[19]) - (gamma_uu.xz * input.ptr[27])) + (2. * gamma_uu.xz * gamma_uu.xy * input.ptr[18]) + ((gamma_uu.zz * input.ptr[14] * gamma_uu.xx) - (2. * gamma_uu.zz * gamma_uu.xx * input.ptr[19]))) / (2. * gamma_uu.xx));\
	(result)->ptr[30] = (((input.ptr[3] * gamma_uu.xx) + ((input.ptr[8] * gamma_uu.xx) - (input.ptr[25] * gamma_uu.xx)) + ((((gammaUUxxSq * input.ptr[15]) - (2. * gamma_uu.xx * gamma_uu.yy * input.ptr[13])) - (2. * gamma_uu.xx * gamma_uu.yz * input.ptr[14])) - (gamma_uu.xx * gamma_uu.zz * input.ptr[20])) + (gamma_uu.xy * input.ptr[5]) + (2. * gamma_uu.xy * gamma_uu.xy * input.ptr[13]) + ((2. * gamma_uu.xy * input.ptr[16] * gamma_uu.xx) - (gamma_uu.xy * input.ptr[27])) + (2. * gamma_uu.xy * gamma_uu.xz * input.ptr[19]) + (gamma_uu.xz * input.ptr[6]) + (2. * gamma_uu.xz * input.ptr[17] * gamma_uu.xx) + ((2. * gamma_uu.xz * gamma_uu.xz * input.ptr[20]) - (gamma_uu.xz * input.ptr[28])) + (2. * gamma_uu.xz * gamma_uu.xy * input.ptr[14]) + (gamma_uu.yy * input.ptr[18] * gamma_uu.xx)) / (2. * gamma_uu.xx));\
\
	(result)->a_l = real3_swap((result)->a_l, n.side);			/* 0-2 */\
	(result)->d_lll = _3sym3_swap((result)->d_lll, n.side);		/* 3-20 */\
	(result)->K_ll = sym3_swap((result)->K_ll, n.side);			/* 21-26 */\
	(result)->Theta = (result)->Theta;							/* 27 */\
	(result)->Z_l = real3_swap((result)->Z_l, n.side);			/* 28-30 */\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?>
// used by roe, some plm
//so long as roeUseFluxFromCons isn't set for the roe solver,
// and fluxFromCons is provided/unused,
// eigen_fluxTransform isn't needed.
// but some solvers do use a boilerplate right(lambda(left(U)))
//however if you want to use the HLL solver then fluxFromCons is needed
//...however fluxFromCons is not provided by this eqn.

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */inputU,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
<? if false then ?>\
	/* default */\
	<?=waves_t?> waves;\
	<?=eigen_leftTransform?>(&waves, solver, eig, inputU, x, n);\
	<?=eqn:eigenWaveCodePrefix{n="n", eig="eig", pt="(cell)->pos"}?>\
<? 	for j=0,eqn.numWaves-1 do --\
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode{n="n", eig="eig", pt="(cell)->pos", waveIndex=j}?>;\
<? 	end --\
?>	<?=eigen_rightTransform?>(resultFlux, solver, eig, waves, (cell)->pos, n);\
<? else ?>\
	for (int j = 0; j < numStates; ++j) {\
		(resultFlux)->ptr[j] = 0./0.;\
	}\
<? end ?>\
}

//// MODULE_NAME: applyKreissOligar
//// MODULE_DEPENDS: <?=coordMapR?>

// TODO this is finite-difference KO
// but the 2009 Alic paper says to use finite-volume KO
static void applyKreissOligar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	global <?=cons_t?> * const deriv,
	int const * const fields,
	int const numFields
) {
	//Kreiss-Oligar dissipation
#if 0	// 2013 Baumgarte et al, section IIIB
	real const coeff = solver->dissipationCoeff / dt;
#elif 1	// 2017 Ruchlin et al eqn 44
	real const r = coordMapR(cell->pos);
	// all params > 0
	real const rKO = 2.;		//rKO/M = 2
	real const wKO = .17;		//wKO/M = .17
	//epsKO0 = .99 doubles as my 'dissipationCoeff' var
	real const coeff = .5 * (erf((r - rKO) / wKO) + 1.) * solver->dissipationCoeff;
#endif

	real3 dy = solver->grid_dx;
<? if require "hydro.coord.sphere_sinh_radial":isa(coord) then ?>
	real3 const yR = _real3(cell->r, cell->pos.y, cell->pos.z);
<? for i=1,solver.dim do
	local xi = xNames[i]
?>{
	global <?=cell_t?> const * const cellL = cell - solver->stepsize.<?=xi?>;
	real3 const yL = _real3(cellL->r, cellL->pos.y, cellL->pos.z);
	dy.<?=xi?> = real3_len(real3_sub(yR, yL));
}<? end ?>
<? end ?>

	real3 const _1_dy = _real3(1./dy.x, 1./dy.y, 1./dy.z);

	//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
	//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
	int const * const endOfFields = fields + numFields;
	for (int const * ip = fields; ip < endOfFields; ++ip) {
		int const i = *ip;
		deriv->ptr[i] += coeff * (0.
<? for j=1,solver.dim do
	local xj = xNames[j]
?>			+ (
				-20. * U->ptr[i]
				+ 15. * (U[1 * solver->stepsize.<?=xj?>].ptr[i] + U[-1 * solver->stepsize.<?=xj?>].ptr[i])
				+ -6. * (U[2 * solver->stepsize.<?=xj?>].ptr[i] + U[-2 * solver->stepsize.<?=xj?>].ptr[i])
				+ U[3 * solver->stepsize.<?=xj?>].ptr[i] + U[-3 * solver->stepsize.<?=xj?>].ptr[i]
			) * _1_dy.<?=xj?>
<? end
?>		) * (1. / 64.);
	}
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=cell_t?> <?=initCond_codeprefix?>


kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useAddSource then ?>
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;

	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real3 const S_l = real3_zero;
	sym3 const S_ll = sym3_zero;
	real const S = 0.;
	real const rho = 0.;

<? if true then ?>
//hand-rolled
	// source terms
	
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real const trK = real3x3_trace(K_ul);							//K^k_k
//	sym3 KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);			//KSq_ij = K_ik K^k_j

//// MODULE_DEPENDS: <?=calc_d_lll?>
	_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);

	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);
	
	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);
	
	//e_l = e_i = d^j_ji
	real3 const e_l = _3sym3_tr12(d_ull);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 const conn_ull = {
<? for k,xk in ipairs(xNames) do
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>			.<?=xij?> = d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?> - d_ull.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end
?>	};

	/* d_l = d_i = d_ij^j */
	real3 const d_l = real3x3x3_tr23(d_llu);
	
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 const Z_u = sym3_real3_mul(gamma_uu, U->Z_l);

	/* d_luu = d_i^jk = gamma^jl d_il^k */
	_3sym3 const d_luu = (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_real3x3_to_sym3_mul(gamma_uu, d_llu.<?=xi?>),
<? end
?>	};

	/* alpha_,t = shift terms - alpha^2 f (gamma^ij K_ij - m Theta) */
	real const f_alphaSq = calc_f_alphaSq(U->alpha);
	deriv->alpha += -f_alphaSq * (trK - solver->m * U->Theta);
	
	/* gamma_ij,t = shift terms - 2 alpha K_ij */
	deriv->gamma_ll = sym3_add(deriv->gamma_ll, sym3_real_mul(U->K_ll, -2. * U->alpha));

	/* 2005 Bona et al A.1 */
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	deriv->K_ll.<?=xij?> += U->alpha * (
		.5 * (U->a_l.<?=xi?> * d_l.<?=xj?> + U->a_l.<?=xj?> * d_l.<?=xi?>)
		- .25 * (U->a_l.<?=xj?> * (2. * e_l.<?=xi?> - d_l.<?=xi?>) + U->a_l.<?=xi?> * (2. * e_l.<?=xj?> - d_l.<?=xj?>))
		- U->a_l.<?=xi?> * U->Z_l.<?=xj?> - U->a_l.<?=xj?> * U->Z_l.<?=xi?>
		+ (trK - 2. * U->Theta) * U->K_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- .5 * U->a_l.<?=xk?> * conn_ull.<?=xk?>.<?=xij?>
		+ .5 * U->a_l.<?=xk?> * d_ull.<?=xk?>.<?=xij?>
		- 2. * e_l.<?=xk?> * (d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?>)
		+ (d_l.<?=xk?> + U->a_l.<?=xk?> - 2. * U->Z_l.<?=xk?>) * conn_ull.<?=xk?>.<?=xij?>
		- 2. * K_ul.<?=xk?>.<?=xi?> * U->K_ll.<?=sym(k,j)?>
<?		for l,xl in ipairs(xNames) do
?>		+ 2. * (d_llu.<?=xi?>.<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,j)?> + d_llu.<?=xj?>.<?=xk?>.<?=xl?> * d_ull.<?=xk?>.<?=sym(l,i)?>)
		- conn_ull.<?=xk?>.<?=sym(l,j)?> * conn_ull.<?=xl?>.<?=sym(k,i)?>
<?		end
	end
?> - 8. * M_PI * (S_ll.<?=xij?> - .5 * U->gamma_ll.<?=xij?> * (S - rho)));
<? end
?>

	/* 2005 Bona et al A.2 */
<? for i,xi in ipairs(xNames) do
?>	deriv->Z_l.<?=xi?> += U->alpha * (
		U->a_l.<?=xi?> * (trK - 2. * U->Theta)
<?	for k,xk in ipairs(xNames) do
?>		- U->a_l.<?=xk?> * K_ul.<?=xk?>.<?=xi?>
		+ K_ul.<?=xk?>.<?=xi?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?		for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * conn_ull.<?=xr?>.<?=sym(k,i)?>
<?		end
	end
?>
	) - 8 * M_PI * U->alpha * S_l.<?=xi?>;
<? end
?>
	/* 2005 Bona et al A.3 */
	deriv->Theta += U->alpha * .5 * ( trK * (trK - 2. * U->Theta)
<?
for k,xk in ipairs(xNames) do
?>		+ 2. * U->a_l.<?=xk?> * (d_u.<?=xk?> - e_u.<?=xk?> - 2. * Z_u.<?=xk?>)
		- d_u.<?=xk?> * (d_l.<?=xk?> - 2. * U->Z_l.<?=xk?>)
<?	for r,xr in ipairs(xNames) do
?>		- K_ul.<?=xk?>.<?=xr?> * K_ul.<?=xr?>.<?=xk?>
<?		for s,xs in ipairs(xNames) do
?>		+ d_luu.<?=xk?>.<?=sym(r,k)?> * conn_ull.<?=xk?>.<?=sym(r,s)?>
<?		end
	end
end?>
	) - 8. * M_PI * U->alpha * rho;


<? if eqn.useShift == "2005 Bona / 2008 Yano" then ?>
	
	real const tr_b = real3x3_trace(U->b_ul);
	real3 const a_u = sym3_real3_mul(gamma_uu, U->a_l);

	/* alpha_,t += beta^i alpha a_i */
<? for k,xk in ipairs(xNames) do
?>	deriv->alpha += U->alpha * real3_dot(U->a_l, U->beta_u);
<? end
?>

	sym3 const dHat_t_ll = sym3_zero;\
	real3x3 const b_ll = sym3_real3x3_mul(U->gamma_ll, U->b_ul);

	/* gamma_ij += 2 beta^k d_kij + b_ij + b_ji */
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	deriv->gamma_ll.<?=xij?> += 0.
		+ b_ll.<?=xi?>.<?=xj?>
		+ b_ll.<?=xj?>.<?=xi?>
<?	for k,xk in ipairs(xNames) do
?>		+ 2. * U->beta_u.<?=xk?> * d_lll.<?=xk?>.<?=xij?>
<?	end
?>		- 2. * dHat_t_ll.<?=xij?>;
<? end
?>

	/* 2005 Bona et al eqn 28 */
	/* a_i,t += b^l_i a_l - b^l_l a_i */
<? for k,xk in ipairs(xNames) do
?>	deriv->a_l.<?=xk?> += 0. 
<? 	for l,xl in ipairs(xNames) do
?>		+ U->b_ul.<?=xl?>.<?=xk?> * U->a_l.<?=xl?>
<?	end
?>		- tr_b * U->a_l.<?=xk?>;
<? end
?>

	/* 2005 Bona et al eqn 30 */
	/* d_kij,t += b^l_k d_lijk - b^l_l d_kij */
<? for k,xk in ipairs(xNames) do	
	for ij,xij in ipairs(symNames) do
?>	deriv->dDelta_lll.<?=xk?>.<?=xij?> += 0.
<?		for l,xl in ipairs(xNames) do
?>		+ U->b_ul.<?=xl?>.<?=xk?> * d_lll.<?=xl?>.<?=xij?>
<?		end
?>		- tr_b * d_lll.<?=xk?>.<?=xij?>;
<? 	end 
end
?>

	/* 2005 Bona et al eqn A.1 */
	/* K_ij,t += K_ik b^k_j + K_jk b^k_i - K_ij b^k_k */
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	deriv->K_ll.<?=xij?> += 0.
<?	for k,xk in ipairs(xNames) do
?>		+ U->K_ll.<?=sym(i,k)?> * U->b_ul.<?=xk?>.<?=xj?>
		+ U->K_ll.<?=sym(j,k)?> * U->b_ul.<?=xk?>.<?=xi?>
<?	end
?>		- tr_b * U->K_ll.<?=xij?>;
<? end
?>

	deriv->Theta += -tr_b * U->Theta;

<? for k,xk in ipairs(xNames) do
?>	deriv->Z_l.<?=xk?> += 0.
		- tr_b * U->Z_l.<?=xk?>
<? 	for l,xl in ipairs(xNames) do
?>		+ U->Z_l.<?=xl?> * U->b_ul.<?=xl?>.<?=xk?>
<? 	end
?>	;
<? end
?>

<? for l,xl in ipairs(xNames) do
?>	deriv->beta_u.<?=xl?> += U->alpha * U->alpha * (2. * e_u.<?=xl?> - d_u.<?=xl?> - a_u.<?=xl?>)
<?	for k,xk in ipairs(xNames) do
?>		+ U->beta_u.<?=xk?> * U->b_ul.<?=xl?>.<?=xk?>
<?	end
?>	;
<? end
?>

<? end -- eqn.useShift == "2005 Bona / 2008 Yano" ?>


<? end ?>
<? if false then ?>
//code-generated

	real const alpha = U->alpha;
	real const Theta = U->Theta;
	real3 const Z_l = U->Z_l;
	real3 const a_l = U->a_l;
	sym3 const gamma_ll = U->gamma_ll;
	sym3 const K_ll = U->K_ll;

	//for 1 + log(alpha) slicing:
	//f = 2 / alpha
	//f alpha = 2
	//df/dalpha = -2 / alpha^2
	//df/dalpha alpha = -2 / alpha
	//df/dalpha alpha^2 = -2
	real const f_alpha = calc_f_alpha(alpha);
	real const f_alphaSq = calc_f_alphaSq(alpha);
	real const alphaSq_dalpha_f = calc_alphaSq_dalpha_f(alpha);

	//TODO:
	//alpha_,t = f alpha^2 ( ... )
	//a_k,t = f alpha ( ... ) + ...
	//
	//  BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/z4_noZeroRows.html
	real const tmp1 = K_ll.xy * gamma_uu.xy;
	real const tmp2 = K_ll.xz * gamma_uu.xz;
	real const tmp3 = K_ll.yz * gamma_uu.yz;
	real const tmp4 = K_ll.zz * gamma_uu.zz;
	real const tmp8 = K_ll.yy * gamma_uu.yy;
	real const tmp37 = gamma_uu.xx * gamma_uu.xx;
	real const tmp43 = gamma_uu.xx * gamma_uu.xy;
	real const tmp49 = gamma_uu.xx * gamma_uu.xz;
	real const tmp55 = gamma_uu.xy * gamma_uu.xy;
	real const tmp61 = gamma_uu.xy * gamma_uu.xz;
	real const tmp67 = gamma_uu.xz * gamma_uu.xz;
	real const tmp86 = gamma_uu.xx * gamma_uu.yy;
	real const tmp98 = gamma_uu.xx * gamma_uu.yz;
	real const tmp110 = gamma_uu.xy * gamma_uu.yy;
	real const tmp116 = gamma_uu.xy * gamma_uu.yz;
	real const tmp122 = gamma_uu.xz * gamma_uu.yy;
	real const tmp128 = gamma_uu.xz * gamma_uu.yz;
	real const tmp159 = gamma_uu.xx * gamma_uu.zz;
	real const tmp177 = gamma_uu.xy * gamma_uu.zz;
	real const tmp189 = gamma_uu.xz * gamma_uu.zz;
	real const tmp218 = gamma_uu.yy * gamma_uu.yy;
	real const tmp224 = gamma_uu.yy * gamma_uu.yz;
	real const tmp230 = gamma_uu.yz * gamma_uu.yz;
	real const tmp279 = gamma_uu.yy * gamma_uu.zz;
	real const tmp291 = gamma_uu.yz * gamma_uu.zz;
	real const tmp332 = gamma_uu.zz * gamma_uu.zz;
	real const tmp1049 = a_l.x * alpha;
	real const tmp1061 = a_l.y * alpha;
	real const tmp1073 = a_l.z * alpha;
	real const tmp1109 = d_lll.x.xx * gamma_uu.xx;
	real const tmp1111 = d_lll.x.xy * gamma_uu.xy;
	real const tmp1113 = d_lll.x.xz * gamma_uu.xz;
	real const tmp1115 = d_lll.y.xx * gamma_uu.xy;
	real const tmp1119 = d_lll.z.xx * gamma_uu.xz;
	real const tmp1123 = d_lll.x.xx * gamma_uu.xy;
	real const tmp1125 = d_lll.x.xy * gamma_uu.yy;
	real const tmp1127 = d_lll.x.xz * gamma_uu.yz;
	real const tmp1129 = d_lll.y.xx * gamma_uu.yy;
	real const tmp1133 = d_lll.z.xx * gamma_uu.yz;
	real const tmp1137 = d_lll.x.xx * gamma_uu.xz;
	real const tmp1139 = d_lll.x.xy * gamma_uu.yz;
	real const tmp1141 = d_lll.x.xz * gamma_uu.zz;
	real const tmp1143 = d_lll.y.xx * gamma_uu.yz;
	real const tmp1147 = d_lll.z.xx * gamma_uu.zz;
	real const tmp1191 = d_lll.x.yy * tmp86;
	real const tmp1195 = d_lll.x.yz * tmp98;
	real const tmp1200 = d_lll.x.zz * tmp159;
	real const tmp1204 = d_lll.y.xy * tmp86;
	real const tmp1207 = d_lll.y.xz * tmp98;
	real const tmp1210 = d_lll.y.yy * tmp110;
	real const tmp1212 = d_lll.y.yz * tmp122;
	real const tmp1215 = d_lll.y.zz * tmp177;
	real const tmp1219 = d_lll.y.zz * tmp128;
	real const tmp1222 = d_lll.z.xy * tmp98;
	real const tmp1225 = d_lll.z.xz * tmp159;
	real const tmp1228 = d_lll.z.yy * tmp116;
	real const tmp1231 = d_lll.z.yy * tmp122;
	real const tmp1235 = d_lll.z.yz * tmp177;
	real const tmp1238 = d_lll.z.zz * tmp189;
	real const tmp1240 = d_lll.x.yy * tmp110;
	real const tmp1245 = d_lll.x.yz * tmp116;
	real const tmp1250 = d_lll.x.zz * tmp177;
	real const tmp1255 = d_lll.y.xy * tmp110;
	real const tmp1258 = d_lll.y.xz * tmp116;
	real const tmp1261 = d_lll.y.yy * tmp218;
	real const tmp1264 = d_lll.y.yz * tmp224;
	real const tmp1267 = d_lll.y.zz * tmp279;
	real const tmp1272 = d_lll.y.zz * tmp230;
	real const tmp1275 = d_lll.z.xy * tmp116;
	real const tmp1278 = d_lll.z.xz * tmp177;
	real const tmp1281 = d_lll.z.yy * tmp224;
	real const tmp1284 = d_lll.z.yz * tmp279;
	real const tmp1287 = d_lll.z.zz * tmp291;
	real const tmp1290 = d_lll.x.yy * tmp122;
	real const tmp1295 = d_lll.x.yz * tmp128;
	real const tmp1300 = d_lll.x.zz * tmp189;
	real const tmp1305 = d_lll.y.xy * tmp122;
	real const tmp1308 = d_lll.y.xz * tmp128;
	real const tmp1311 = d_lll.y.yy * tmp224;
	real const tmp1314 = d_lll.y.yz * tmp279;
	real const tmp1317 = d_lll.y.zz * tmp291;
	real const tmp1320 = d_lll.z.xy * tmp128;
	real const tmp1323 = d_lll.z.xz * tmp189;
	real const tmp1326 = d_lll.z.yy * tmp279;
	real const tmp1331 = d_lll.z.yy * tmp230;
	real const tmp1334 = d_lll.z.yz * tmp291;
	real const tmp1337 = d_lll.z.zz * tmp332;
	real const tmp1340 = d_lll.x.yz * tmp224;
	real const tmp1345 = d_lll.x.zz * tmp230;
	real const tmp1359 = d_lll.z.xx * tmp122;
	real const tmp1360 = d_lll.x.yy * d_lll.x.yy;
	real const tmp1363 = tmp1360 * tmp218;
	real const tmp1365 = d_lll.x.zz * tmp291;
	real const tmp1375 = d_lll.z.xx * tmp177;
	real const tmp1379 = d_lll.x.yz * d_lll.x.yz;
	real const tmp1381 = tmp1379 * tmp230;
	real const tmp1395 = d_lll.x.zz * d_lll.x.zz;
	real const tmp1398 = tmp1395 * tmp332;
	real const tmp1405 = d_lll.y.xz * tmp122;
	real const tmp1426 = d_lll.z.xx * tmp98;
	real const tmp1436 = d_lll.z.xy * tmp122;
	real const tmp1465 = d_lll.z.xx * tmp128;
	real const tmp1470 = d_lll.z.xy * tmp279;
	real const tmp1473 = d_lll.z.xy * tmp230;
	real const tmp1477 = d_lll.y.xz * d_lll.y.xz;
	real const tmp1479 = tmp1477 * tmp230;
	real const tmp1494 = d_lll.z.xy * tmp177;
	real const tmp1519 = d_lll.z.xy * d_lll.z.xy;
	real const tmp1521 = tmp1519 * tmp230;
	real const tmp1523 = d_lll.y.xx * d_lll.y.xx;
	real const tmp1527 = d_lll.z.xx * d_lll.z.xx;
	real const tmp1531 = K_ll.xy * K_ll.xy;
	real const tmp1532 = gamma_uu.yy * tmp1531;
	real const tmp1534 = gamma_uu.zz * tmp1379;
	real const tmp1535 = gamma_uu.yy * tmp1534;
	real const tmp1539 = gamma_uu.zz * tmp1477;
	real const tmp1540 = gamma_uu.yy * tmp1539;
	real const tmp1544 = gamma_uu.zz * tmp1519;
	real const tmp1545 = gamma_uu.yy * tmp1544;
	real const tmp1548 = K_ll.xz * K_ll.xz;
	real const tmp1549 = gamma_uu.zz * tmp1548;
	real const tmp1679 = K_ll.yy * gamma_uu.xy;
	real const tmp1681 = K_ll.yz * gamma_uu.xz;
	real const tmp1699 = d_lll.x.yy * gamma_uu.xy;
	real const tmp1701 = d_lll.x.yz * gamma_uu.xz;
	real const tmp1703 = d_lll.y.xx * gamma_uu.xx;
	real const tmp1705 = d_lll.y.xz * gamma_uu.xz;
	real const tmp1707 = d_lll.z.xy * gamma_uu.xz;
	real const tmp1711 = d_lll.x.yy * gamma_uu.yy;
	real const tmp1713 = d_lll.x.yz * gamma_uu.yz;
	real const tmp1717 = d_lll.y.xz * gamma_uu.yz;
	real const tmp1719 = d_lll.z.xy * gamma_uu.yz;
	real const tmp1723 = d_lll.x.yy * gamma_uu.yz;
	real const tmp1725 = d_lll.x.yz * gamma_uu.zz;
	real const tmp1727 = d_lll.y.xx * gamma_uu.xz;
	real const tmp1729 = d_lll.y.xz * gamma_uu.zz;
	real const tmp1731 = d_lll.z.xy * gamma_uu.zz;
	real const tmp1775 = d_lll.x.yy * tmp43;
	real const tmp1777 = d_lll.x.yz * tmp49;
	real const tmp1779 = d_lll.y.xy * tmp43;
	real const tmp1784 = d_lll.y.xz * tmp49;
	real const tmp1788 = d_lll.y.yy * tmp55;
	real const tmp1792 = d_lll.y.yz * tmp61;
	real const tmp1797 = d_lll.y.zz * tmp67;
	real const tmp1801 = d_lll.z.xy * tmp49;
	real const tmp1805 = d_lll.x.yy * tmp55;
	real const tmp1808 = d_lll.x.yz * tmp61;
	real const tmp1811 = d_lll.y.xy * tmp55;
	real const tmp1816 = d_lll.y.xz * tmp61;
	real const tmp1826 = d_lll.y.yz * tmp116;
	real const tmp1836 = d_lll.z.xy * tmp61;
	real const tmp1855 = d_lll.y.xy * tmp61;
	real const tmp1860 = d_lll.y.xz * tmp67;
	real const tmp1865 = d_lll.y.yy * tmp116;
	real const tmp1870 = d_lll.y.yz * tmp177;
	real const tmp1875 = d_lll.y.zz * tmp189;
	real const tmp1880 = d_lll.z.xy * tmp67;
	real const tmp1885 = d_lll.z.yy * tmp177;
	real const tmp1888 = d_lll.z.yy * tmp128;
	real const tmp1896 = d_lll.x.yz * tmp122;
	real const tmp1902 = d_lll.x.zz * tmp128;
	real const tmp1905 = d_lll.y.xx * tmp55;
	real const tmp1936 = d_lll.y.xx * tmp98;
	real const tmp1940 = d_lll.y.xx * tmp61;
	real const tmp1943 = d_lll.y.xy * tmp116;
	real const tmp1952 = d_lll.y.yz * tmp230;
	real const tmp1961 = d_lll.z.xx * tmp159;
	real const tmp1990 = d_lll.y.xx * tmp67;
	real const tmp1992 = d_lll.y.xy * tmp128;
	real const tmp1997 = d_lll.y.xz * tmp189;
	real const tmp2001 = d_lll.y.yy * tmp230;
	real const tmp2005 = d_lll.y.yz * tmp291;
	real const tmp2010 = d_lll.y.zz * tmp332;
	real const tmp2038 = d_lll.z.xx * tmp49;
	real const tmp2060 = d_lll.z.xx * tmp61;
	real const tmp2079 = d_lll.z.xx * tmp67;
	real const tmp2100 = d_lll.z.xy * tmp224;
	real const tmp2121 = d_lll.z.xy * tmp159;
	real const tmp2153 = gamma_uu.yy * tmp1360;
	real const tmp2161 = gamma_uu.yz * tmp1519;
	real const tmp2162 = gamma_uu.xz * tmp2161;
	real const tmp2331 = K_ll.xz * gamma_uu.xx;
	real const tmp2332 = K_ll.yz * gamma_uu.xy;
	real const tmp2334 = K_ll.zz * gamma_uu.xz;
	real const tmp2336 = K_ll.yz * gamma_uu.yy;
	real const tmp2338 = K_ll.zz * gamma_uu.yz;
	real const tmp2352 = d_lll.x.yz * gamma_uu.xy;
	real const tmp2354 = d_lll.x.zz * gamma_uu.xz;
	real const tmp2356 = d_lll.y.xz * gamma_uu.xy;
	real const tmp2360 = d_lll.z.xx * gamma_uu.xx;
	real const tmp2362 = d_lll.z.xy * gamma_uu.xy;
	real const tmp2364 = d_lll.x.yz * gamma_uu.yy;
	real const tmp2366 = d_lll.x.zz * gamma_uu.yz;
	real const tmp2368 = d_lll.y.xz * gamma_uu.yy;
	real const tmp2372 = d_lll.z.xx * gamma_uu.xy;
	real const tmp2374 = d_lll.z.xy * gamma_uu.yy;
	real const tmp2378 = d_lll.x.zz * gamma_uu.zz;
	real const tmp2428 = d_lll.x.yz * tmp43;
	real const tmp2430 = d_lll.x.zz * tmp49;
	real const tmp2432 = d_lll.y.xz * tmp43;
	real const tmp2436 = d_lll.z.xy * tmp43;
	real const tmp2440 = d_lll.z.xz * tmp49;
	real const tmp2445 = d_lll.z.yy * tmp55;
	real const tmp2449 = d_lll.z.yz * tmp61;
	real const tmp2454 = d_lll.z.zz * tmp67;
	real const tmp2461 = d_lll.x.zz * tmp61;
	real const tmp2464 = d_lll.y.xz * tmp55;
	real const tmp2469 = d_lll.y.zz * tmp116;
	real const tmp2474 = d_lll.y.zz * tmp122;
	real const tmp2477 = d_lll.z.xy * tmp55;
	real const tmp2482 = d_lll.z.xz * tmp61;
	real const tmp2487 = d_lll.z.yy * tmp110;
	real const tmp2492 = d_lll.z.yz * tmp122;
	real const tmp2497 = d_lll.z.zz * tmp128;
	real const tmp2505 = d_lll.x.zz * tmp67;
	real const tmp2526 = d_lll.z.xz * tmp67;
	real const tmp2536 = d_lll.z.yz * tmp128;
	real const tmp2555 = d_lll.y.xz * tmp110;
	real const tmp2559 = d_lll.z.xx * tmp86;
	real const tmp2563 = d_lll.z.xx * tmp55;
	real const tmp2565 = d_lll.z.xy * tmp110;
	real const tmp2569 = d_lll.z.xz * tmp116;
	real const tmp2574 = d_lll.z.yy * tmp218;
	real const tmp2578 = d_lll.z.yz * tmp224;
	real const tmp2583 = d_lll.z.zz * tmp230;
	real const tmp2592 = d_lll.y.xx * tmp86;
	real const tmp2624 = d_lll.z.xz * tmp128;
	real const tmp2633 = d_lll.z.yz * tmp230;
	real const tmp2647 = d_lll.y.xz * tmp177;
	real const tmp2671 = d_lll.y.xz * tmp86;
	real const tmp2683 = d_lll.z.xx * tmp43;
	real const tmp2685 = d_lll.z.xy * tmp86;
	real const tmp2692 = d_lll.z.xz * tmp98;
	real const tmp2700 = d_lll.z.yz * tmp116;
	real const tmp2806 = gamma_uu.yz * tmp1477;
	real const tmp2811 = gamma_uu.yy * tmp1379;
	real const tmp2814 = gamma_uu.yy * tmp1519;
	real const tmp2817 = gamma_uu.zz * tmp1395;
	real const tmp3012 = d_lll.x.yy * gamma_uu.xx;
	real const tmp3014 = d_lll.y.xy * gamma_uu.xx;
	real const tmp3018 = d_lll.y.yy * gamma_uu.xy;
	real const tmp3022 = d_lll.y.yz * gamma_uu.xz;
	real const tmp3026 = d_lll.z.yy * gamma_uu.xz;
	real const tmp3030 = d_lll.y.xy * gamma_uu.xy;
	real const tmp3034 = d_lll.y.yy * gamma_uu.yy;
	real const tmp3038 = d_lll.y.yz * gamma_uu.yz;
	real const tmp3042 = d_lll.z.yy * gamma_uu.yz;
	real const tmp3044 = d_lll.x.yy * gamma_uu.xz;
	real const tmp3046 = d_lll.y.xy * gamma_uu.xz;
	real const tmp3050 = d_lll.y.yy * gamma_uu.yz;
	real const tmp3054 = d_lll.y.yz * gamma_uu.zz;
	real const tmp3058 = d_lll.z.yy * gamma_uu.zz;
	real const tmp3122 = d_lll.y.yy * tmp86;
	real const tmp3127 = d_lll.y.yz * tmp98;
	real const tmp3132 = d_lll.z.yy * tmp98;
	real const tmp3153 = d_lll.z.yy * tmp159;
	real const tmp3166 = d_lll.y.xx * tmp43;
	real const tmp3204 = d_lll.y.yy * tmp122;
	real const tmp3209 = d_lll.y.yz * tmp128;
	real const tmp3231 = tmp1379 * tmp67;
	real const tmp3255 = d_lll.z.yy * tmp189;
	real const tmp3275 = d_lll.z.yy * tmp61;
	real const tmp3324 = tmp1477 * tmp67;
	real const tmp3386 = d_lll.y.zz * d_lll.y.zz;
	real const tmp3390 = tmp3386 * tmp332;
	real const tmp3400 = tmp1519 * tmp67;
	real const tmp3429 = d_lll.z.yy * d_lll.z.yy;
	real const tmp3431 = K_ll.yz * K_ll.yz;
	real const tmp3432 = gamma_uu.zz * tmp3431;
	real const tmp3479 = tmp1523 * tmp37;
	real const tmp3611 = d_lll.x.yz * gamma_uu.xx;
	real const tmp3613 = d_lll.y.xz * gamma_uu.xx;
	real const tmp3617 = d_lll.y.zz * gamma_uu.xz;
	real const tmp3621 = d_lll.z.xy * gamma_uu.xx;
	real const tmp3625 = d_lll.z.yy * gamma_uu.xy;
	real const tmp3635 = d_lll.y.zz * gamma_uu.yz;
	real const tmp3643 = d_lll.z.yy * gamma_uu.yy;
	real const tmp3653 = d_lll.y.zz * gamma_uu.zz;
	real const tmp3715 = d_lll.y.zz * tmp98;
	real const tmp3725 = d_lll.z.yy * tmp86;
	real const tmp3755 = d_lll.x.zz * tmp98;
	real const tmp3844 = d_lll.y.xy * tmp98;
	real const tmp3852 = d_lll.y.xz * tmp159;
	real const tmp3890 = d_lll.z.xx * tmp37;
	real const tmp3968 = d_lll.y.zz * tmp224;
	real const tmp4086 = gamma_uu.zz * tmp3386;
	real const tmp4256 = d_lll.x.zz * gamma_uu.xx;
	real const tmp4258 = d_lll.y.zz * gamma_uu.xy;
	real const tmp4260 = d_lll.z.xz * gamma_uu.xx;
	real const tmp4264 = d_lll.z.yz * gamma_uu.xy;
	real const tmp4268 = d_lll.z.zz * gamma_uu.xz;
	real const tmp4272 = d_lll.x.zz * gamma_uu.xy;
	real const tmp4274 = d_lll.y.zz * gamma_uu.yy;
	real const tmp4276 = d_lll.z.xz * gamma_uu.xy;
	real const tmp4280 = d_lll.z.yz * gamma_uu.yy;
	real const tmp4284 = d_lll.z.zz * gamma_uu.yz;
	real const tmp4292 = d_lll.z.xz * gamma_uu.xz;
	real const tmp4296 = d_lll.z.yz * gamma_uu.yz;
	real const tmp4300 = d_lll.z.zz * gamma_uu.zz;
	real const tmp4361 = d_lll.y.zz * tmp86;
	real const tmp4364 = d_lll.z.xz * tmp43;
	real const tmp4369 = d_lll.z.yz * tmp86;
	real const tmp4374 = d_lll.z.zz * tmp98;
	real const tmp4390 = d_lll.z.yz * tmp98;
	real const tmp4395 = d_lll.z.zz * tmp159;
	real const tmp4407 = d_lll.y.zz * tmp110;
	real const tmp4409 = d_lll.z.xz * tmp86;
	real const tmp4417 = d_lll.z.yz * tmp110;
	real const tmp4427 = d_lll.z.zz * tmp122;
	real const tmp4458 = d_lll.z.zz * tmp177;
	real const tmp4464 = tmp1379 * tmp55;
	real const tmp4565 = tmp1477 * tmp55;
	real const tmp4598 = d_lll.z.zz * tmp279;
	real const tmp4644 = tmp1519 * tmp55;
	real const tmp4657 = tmp3429 * tmp218;
	real const tmp4705 = tmp1527 * tmp37;
	real const tmp4920 = d_lll.x.xx * tmp37;
	real const tmp4924 = d_lll.x.xy * tmp43;
	real const tmp4929 = d_lll.x.xz * tmp49;
	real const tmp5025 = d_lll.x.xx * tmp43;
	real const tmp5029 = d_lll.x.xy * tmp86;
	real const tmp5034 = d_lll.x.xz * tmp98;
	real const tmp5130 = d_lll.x.xx * tmp49;
	real const tmp5134 = d_lll.x.xy * tmp98;
	real const tmp5139 = d_lll.x.xz * tmp159;
	real const tmp5144 = d_lll.x.yy * tmp116;
	real const tmp5151 = d_lll.x.yz * tmp177;
	real const tmp5232 = gamma_uu.xx * tmp55;
	real const tmp5237 = gamma_uu.yy * tmp37;
	real const tmp5240 = gamma_uu.xx * tmp61;
	real const tmp5246 = gamma_uu.yz * tmp37;
	real const tmp5250 = gamma_uu.xx * tmp67;
	real const tmp5255 = gamma_uu.zz * tmp37;
	real const tmp5278 = gamma_uu.xx * tmp110;
	real const tmp5283 = gamma_uu.xy * tmp55;
	real const tmp5286 = gamma_uu.xx * tmp122;
	real const tmp5292 = gamma_uu.xz * tmp55;
	real const tmp5296 = gamma_uu.xx * tmp177;
	real const tmp5299 = gamma_uu.xx * tmp128;
	real const tmp5305 = gamma_uu.xy * tmp67;
	real const tmp5328 = gamma_uu.xx * tmp116;
	real const tmp5350 = gamma_uu.xx * tmp189;
	real const tmp5355 = gamma_uu.xz * tmp67;
	real const tmp5389 = d_lll.y.xy * tmp5278;
	real const tmp5395 = d_lll.y.xy * tmp5283;
	real const tmp5399 = d_lll.y.xz * tmp5328;
	real const tmp5405 = d_lll.y.xz * tmp5292;
	real const tmp5408 = gamma_uu.xx * tmp218;
	real const tmp5409 = d_lll.y.yy * tmp5408;
	real const tmp5414 = gamma_uu.yy * tmp55;
	real const tmp5415 = d_lll.y.yy * tmp5414;
	real const tmp5418 = gamma_uu.xx * tmp224;
	real const tmp5419 = d_lll.y.yz * tmp5418;
	real const tmp5424 = gamma_uu.yz * tmp55;
	real const tmp5425 = d_lll.y.yz * tmp5424;
	real const tmp5428 = gamma_uu.xx * tmp279;
	real const tmp5429 = d_lll.y.zz * tmp5428;
	real const tmp5432 = gamma_uu.xx * tmp230;
	real const tmp5433 = d_lll.y.zz * tmp5432;
	real const tmp5438 = gamma_uu.xy * tmp128;
	real const tmp5439 = d_lll.y.zz * tmp5438;
	real const tmp5442 = gamma_uu.yy * tmp67;
	real const tmp5443 = d_lll.y.zz * tmp5442;
	real const tmp5449 = d_lll.z.xy * tmp5328;
	real const tmp5455 = d_lll.z.xy * tmp5292;
	real const tmp5459 = d_lll.z.xz * tmp5296;
	real const tmp5465 = d_lll.z.xz * tmp5305;
	real const tmp5469 = d_lll.z.yy * tmp5418;
	real const tmp5474 = gamma_uu.xy * tmp122;
	real const tmp5479 = d_lll.z.yy * tmp5424;
	real const tmp5485 = d_lll.z.yz * tmp5428;
	real const tmp5491 = d_lll.z.yz * tmp5442;
	real const tmp5494 = gamma_uu.xx * tmp291;
	real const tmp5495 = d_lll.z.zz * tmp5494;
	real const tmp5500 = gamma_uu.yz * tmp67;
	real const tmp5501 = d_lll.z.zz * tmp5500;
	real const tmp5570 = gamma_uu.zz * tmp55;
	real const tmp5580 = gamma_uu.xy * tmp189;
	real const tmp5591 = d_lll.z.xy * tmp5299;
	real const tmp5597 = d_lll.z.xy * tmp5305;
	real const tmp5601 = d_lll.z.xz * tmp5350;
	real const tmp5607 = d_lll.z.xz * tmp5355;
	real const tmp5611 = d_lll.z.yy * tmp5428;
	real const tmp5615 = d_lll.z.yy * tmp5432;
	real const tmp5621 = d_lll.z.yy * tmp5438;
	real const tmp5625 = d_lll.z.yy * tmp5570;
	real const tmp5631 = d_lll.z.yz * tmp5494;
	real const tmp5637 = d_lll.z.yz * tmp5500;
	real const tmp5640 = gamma_uu.xx * tmp332;
	real const tmp5641 = d_lll.z.zz * tmp5640;
	real const tmp5646 = gamma_uu.zz * tmp67;
	real const tmp5647 = d_lll.z.zz * tmp5646;
	real const tmp5698 = d_lll.y.xz * tmp5474;
	real const tmp5702 = d_lll.y.xz * tmp5424;
	real const tmp5707 = gamma_uu.xy * tmp279;
	real const tmp5708 = d_lll.y.zz * tmp5707;
	real const tmp5710 = gamma_uu.xy * tmp230;
	real const tmp5711 = d_lll.y.zz * tmp5710;
	real const tmp5716 = d_lll.z.xx * tmp5286;
	real const tmp5719 = d_lll.z.xx * tmp5292;
	real const tmp5724 = d_lll.z.xy * tmp5474;
	real const tmp5728 = d_lll.z.xy * tmp5424;
	real const tmp5734 = d_lll.z.xz * tmp5428;
	real const tmp5738 = d_lll.z.xz * tmp5432;
	real const tmp5753 = gamma_uu.xy * tmp224;
	real const tmp5759 = gamma_uu.xz * tmp218;
	real const tmp5764 = d_lll.z.yz * tmp5707;
	real const tmp5770 = d_lll.z.yz * tmp5710;
	real const tmp5775 = gamma_uu.xz * tmp224;
	real const tmp5779 = gamma_uu.xy * tmp291;
	real const tmp5785 = gamma_uu.xz * tmp279;
	real const tmp5786 = d_lll.z.zz * tmp5785;
	real const tmp5788 = gamma_uu.xz * tmp230;
	real const tmp5789 = d_lll.z.zz * tmp5788;
	real const tmp5852 = d_lll.y.yy * tmp5753;
	real const tmp5856 = d_lll.y.yy * tmp5759;
	real const tmp5862 = d_lll.y.yz * tmp5710;
	real const tmp5866 = d_lll.y.yz * tmp5775;
	real const tmp5872 = d_lll.y.zz * tmp5779;
	real const tmp5876 = d_lll.y.zz * tmp5788;
	real const tmp5882 = d_lll.z.xx * tmp5299;
	real const tmp5886 = d_lll.z.xx * tmp5305;
	real const tmp5892 = d_lll.z.xy * tmp5428;
	real const tmp5898 = d_lll.z.xy * tmp5432;
	real const tmp5902 = d_lll.z.xy * tmp5438;
	real const tmp5908 = d_lll.z.xy * tmp5442;
	real const tmp5912 = d_lll.z.xy * tmp5570;
	real const tmp5916 = d_lll.z.xz * tmp5580;
	real const tmp5922 = d_lll.z.xz * tmp5500;
	real const tmp5926 = d_lll.z.yy * tmp5710;
	real const tmp5932 = d_lll.z.yy * tmp5775;
	real const tmp5936 = d_lll.z.yz * tmp5779;
	real const tmp5942 = d_lll.z.yz * tmp5788;
	real const tmp5945 = gamma_uu.xy * tmp332;
	real const tmp5946 = d_lll.z.zz * tmp5945;
	real const tmp5951 = gamma_uu.xz * tmp291;
	real const tmp5952 = d_lll.z.zz * tmp5951;
	real const tmp6050 = d_lll.z.yy * tmp5785;
	real const tmp6053 = d_lll.z.yy * tmp5788;
	real const tmp6386 = gamma_uu.yy * tmp230;
	real const tmp6391 = gamma_uu.zz * tmp218;
	real const tmp6436 = gamma_uu.yy * tmp291;
	real const tmp6441 = gamma_uu.yz * tmp230;
	real const tmp6495 = d_lll.z.yz * tmp6436;
	real const tmp6501 = d_lll.z.yz * tmp6441;
	real const tmp6504 = gamma_uu.yy * tmp332;
	real const tmp6505 = d_lll.z.zz * tmp6504;
	real const tmp6510 = gamma_uu.zz * tmp230;
	real const tmp6511 = d_lll.z.zz * tmp6510;
	real const tmp7192 = Z_l.x * gamma_uu.xx;
	real const tmp7194 = Z_l.y * gamma_uu.xy;
	real const tmp7196 = Z_l.z * gamma_uu.xz;
	real const tmp7258 = Z_l.x * gamma_uu.xy;
	real const tmp7260 = Z_l.y * gamma_uu.yy;
	real const tmp7262 = Z_l.z * gamma_uu.yz;
	real const tmp7341 = Z_l.x * gamma_uu.xz;
	real const tmp7343 = Z_l.y * gamma_uu.yz;
	real const tmp7345 = Z_l.z * gamma_uu.zz;
	deriv->alpha += -f_alphaSq * (tmp8 + tmp4 - 2. * Theta + 2. * tmp3 + 2. * tmp2 + 2. * tmp1 + K_ll.xx * gamma_uu.xx);
	deriv->gamma_ll.xx += -2. * K_ll.xx * alpha;
	deriv->gamma_ll.xy += -2. * K_ll.xy * alpha;
	deriv->gamma_ll.xz += -2. * K_ll.xz * alpha;
	deriv->gamma_ll.yy += -2. * K_ll.yy * alpha;
	deriv->gamma_ll.yz += -2. * K_ll.yz * alpha;
	deriv->gamma_ll.zz += -2. * K_ll.zz * alpha;
	deriv->a_l.x += -alphaSq_dalpha_f * ( K_ll.zz * a_l.x * gamma_uu.zz - 2. * Theta * a_l.x + 2. * K_ll.yz * a_l.x * gamma_uu.yz + K_ll.yy * a_l.x * gamma_uu.yy + 2. * K_ll.xz * a_l.x * gamma_uu.xz + 2. * K_ll.xy * a_l.x * gamma_uu.xy + K_ll.xx * a_l.x * gamma_uu.xx) - f_alpha * ( - 2. * K_ll.zz * d_lll.x.xx * tmp67 - 4. * K_ll.zz * d_lll.x.xy * tmp128 - 4. * K_ll.zz * d_lll.x.xz * tmp189 - 2. * K_ll.zz * d_lll.x.yy * tmp230 - 4. * K_ll.zz * d_lll.x.yz * tmp291 - 2. * K_ll.zz * d_lll.x.zz * tmp332 - 2. * Theta * a_l.x + K_ll.zz * a_l.x * gamma_uu.zz - 4. * K_ll.yz * d_lll.x.xx * tmp61 - 4. * K_ll.yz * d_lll.x.xy * tmp116 - 4. * K_ll.yz * d_lll.x.xy * tmp122 - 4. * K_ll.yz * d_lll.x.xz * tmp177 - 4. * K_ll.yz * d_lll.x.xz * tmp128 - 4. * K_ll.yz * d_lll.x.yy * tmp224 - 4. * K_ll.yz * d_lll.x.yz * tmp279 - 4. * K_ll.yz * d_lll.x.yz * tmp230 - 4. * K_ll.yz * d_lll.x.zz * tmp291 + 2. * K_ll.yz * a_l.x * gamma_uu.yz - 2. * K_ll.yy * d_lll.x.xx * tmp55 - 4. * K_ll.yy * d_lll.x.xy * tmp110 - 4. * K_ll.yy * d_lll.x.xz * tmp116 - 2. * K_ll.yy * d_lll.x.yy * tmp218 - 4. * K_ll.yy * d_lll.x.yz * tmp224 - 2. * K_ll.yy * d_lll.x.zz * tmp230 + K_ll.yy * a_l.x * gamma_uu.yy - 4. * K_ll.xz * d_lll.x.xx * tmp49 - 4. * K_ll.xz * d_lll.x.xy * tmp98 - 4. * K_ll.xz * d_lll.x.xy * tmp61 - 4. * K_ll.xz * d_lll.x.xz * tmp159 - 4. * K_ll.xz * d_lll.x.xz * tmp67 - 4. * K_ll.xz * d_lll.x.yy * tmp116 - 4. * K_ll.xz * d_lll.x.yz * tmp177 - 4. * K_ll.xz * d_lll.x.yz * tmp128 - 4. * K_ll.xz * d_lll.x.zz * tmp189 + 2. * K_ll.xz * a_l.x * gamma_uu.xz - 4. * K_ll.xy * d_lll.x.xx * tmp43 - 4. * K_ll.xy * d_lll.x.xy * tmp86 - 4. * K_ll.xy * d_lll.x.xy * tmp55 - 4. * K_ll.xy * d_lll.x.xz * tmp98 - 4. * K_ll.xy * d_lll.x.xz * tmp61 - 4. * K_ll.xy * d_lll.x.yy * tmp110 - 4. * K_ll.xy * d_lll.x.yz * tmp116 - 4. * K_ll.xy * d_lll.x.yz * tmp122 - 4. * K_ll.xy * d_lll.x.zz * tmp128 + 2. * K_ll.xy * a_l.x * gamma_uu.xy - 2. * K_ll.xx * d_lll.x.xx * tmp37 - 4. * K_ll.xx * d_lll.x.xy * tmp43 - 4. * K_ll.xx * d_lll.x.xz * tmp49 - 2. * K_ll.xx * d_lll.x.yy * tmp55 - 4. * K_ll.xx * d_lll.x.yz * tmp61 - 2. * K_ll.xx * d_lll.x.zz * tmp67 + K_ll.xx * a_l.x * gamma_uu.xx);
	deriv->a_l.y += -alphaSq_dalpha_f * ( K_ll.zz * a_l.y * gamma_uu.zz - 2. * Theta * a_l.y + 2. * K_ll.yz * a_l.y * gamma_uu.yz + K_ll.yy * a_l.y * gamma_uu.yy + 2. * K_ll.xz * a_l.y * gamma_uu.xz + 2. * K_ll.xy * a_l.y * gamma_uu.xy + K_ll.xx * a_l.y * gamma_uu.xx) - f_alpha * ( - 2. * K_ll.zz * d_lll.y.xx * tmp67 - 4. * K_ll.zz * d_lll.y.xy * tmp128 - 4. * K_ll.zz * d_lll.y.xz * tmp189 - 2. * K_ll.zz * d_lll.y.yy * tmp230 - 4. * K_ll.zz * d_lll.y.yz * tmp291 - 2. * K_ll.zz * d_lll.y.zz * tmp332 - 2. * Theta * a_l.y + K_ll.zz * a_l.y * gamma_uu.zz - 4. * K_ll.yz * d_lll.y.xx * tmp61 - 4. * K_ll.yz * d_lll.y.xy * tmp116 - 4. * K_ll.yz * d_lll.y.xy * tmp122 - 4. * K_ll.yz * d_lll.y.xz * tmp177 - 4. * K_ll.yz * d_lll.y.xz * tmp128 - 4. * K_ll.yz * d_lll.y.yy * tmp224 - 4. * K_ll.yz * d_lll.y.yz * tmp279 - 4. * K_ll.yz * d_lll.y.yz * tmp230 - 4. * K_ll.yz * d_lll.y.zz * tmp291 + 2. * K_ll.yz * a_l.y * gamma_uu.yz - 2. * K_ll.yy * d_lll.y.xx * tmp55 - 4. * K_ll.yy * d_lll.y.xy * tmp110 - 4. * K_ll.yy * d_lll.y.xz * tmp116 - 2. * K_ll.yy * d_lll.y.yy * tmp218 - 4. * K_ll.yy * d_lll.y.yz * tmp224 - 2. * K_ll.yy * d_lll.y.zz * tmp230 + K_ll.yy * a_l.y * gamma_uu.yy - 4. * K_ll.xz * d_lll.y.xx * tmp49 - 4. * K_ll.xz * d_lll.y.xy * tmp98 - 4. * K_ll.xz * d_lll.y.xy * tmp61 - 4. * K_ll.xz * d_lll.y.xz * tmp159 - 4. * K_ll.xz * d_lll.y.xz * tmp67 - 4. * K_ll.xz * d_lll.y.yy * tmp116 - 4. * K_ll.xz * d_lll.y.yz * tmp177 - 4. * K_ll.xz * d_lll.y.yz * tmp128 - 4. * K_ll.xz * d_lll.y.zz * tmp189 + 2. * K_ll.xz * a_l.y * gamma_uu.xz - 4. * K_ll.xy * d_lll.y.xx * tmp43 - 4. * K_ll.xy * d_lll.y.xy * tmp86 - 4. * K_ll.xy * d_lll.y.xy * tmp55 - 4. * K_ll.xy * d_lll.y.xz * tmp98 - 4. * K_ll.xy * d_lll.y.xz * tmp61 - 4. * K_ll.xy * d_lll.y.yy * tmp110 - 4. * K_ll.xy * d_lll.y.yz * tmp116 - 4. * K_ll.xy * d_lll.y.yz * tmp122 - 4. * K_ll.xy * d_lll.y.zz * tmp128 + 2. * K_ll.xy * a_l.y * gamma_uu.xy - 2. * K_ll.xx * d_lll.y.xx * tmp37 - 4. * K_ll.xx * d_lll.y.xy * tmp43 - 4. * K_ll.xx * d_lll.y.xz * tmp49 - 2. * K_ll.xx * d_lll.y.yy * tmp55 - 4. * K_ll.xx * d_lll.y.yz * tmp61 - 2. * K_ll.xx * d_lll.y.zz * tmp67 + K_ll.xx * a_l.y * gamma_uu.xx);
	deriv->a_l.z += -alphaSq_dalpha_f * ( K_ll.zz * a_l.z * gamma_uu.zz - 2. * Theta * a_l.z + 2. * K_ll.yz * a_l.z * gamma_uu.yz + K_ll.yy * a_l.z * gamma_uu.yy + 2. * K_ll.xz * a_l.z * gamma_uu.xz + 2. * K_ll.xy * a_l.z * gamma_uu.xy + K_ll.xx * a_l.z * gamma_uu.xx) - f_alpha * ( - 2. * K_ll.zz * d_lll.z.xx * tmp67 - 4. * K_ll.zz * d_lll.z.xy * tmp128 - 4. * K_ll.zz * d_lll.z.xz * tmp189 - 2. * K_ll.zz * d_lll.z.yy * tmp230 - 4. * K_ll.zz * d_lll.z.yz * tmp291 - 2. * K_ll.zz * d_lll.z.zz * tmp332 - 2. * Theta * a_l.z + K_ll.zz * a_l.z * gamma_uu.zz - 4. * K_ll.yz * d_lll.z.xx * tmp61 - 4. * K_ll.yz * d_lll.z.xy * tmp116 - 4. * K_ll.yz * d_lll.z.xy * tmp122 - 4. * K_ll.yz * d_lll.z.xz * tmp177 - 4. * K_ll.yz * d_lll.z.xz * tmp128 - 4. * K_ll.yz * d_lll.z.yy * tmp224 - 4. * K_ll.yz * d_lll.z.yz * tmp279 - 4. * K_ll.yz * d_lll.z.yz * tmp230 - 4. * K_ll.yz * d_lll.z.zz * tmp291 + 2. * K_ll.yz * a_l.z * gamma_uu.yz - 2. * K_ll.yy * d_lll.z.xx * tmp55 - 4. * K_ll.yy * d_lll.z.xy * tmp110 - 4. * K_ll.yy * d_lll.z.xz * tmp116 - 2. * K_ll.yy * d_lll.z.yy * tmp218 - 4. * K_ll.yy * d_lll.z.yz * tmp224 - 2. * K_ll.yy * d_lll.z.zz * tmp230 + K_ll.yy * a_l.z * gamma_uu.yy - 4. * K_ll.xz * d_lll.z.xx * tmp49 - 4. * K_ll.xz * d_lll.z.xy * tmp98 - 4. * K_ll.xz * d_lll.z.xy * tmp61 - 4. * K_ll.xz * d_lll.z.xz * tmp159 - 4. * K_ll.xz * d_lll.z.xz * tmp67 - 4. * K_ll.xz * d_lll.z.yy * tmp116 - 4. * K_ll.xz * d_lll.z.yz * tmp177 - 4. * K_ll.xz * d_lll.z.yz * tmp128 - 4. * K_ll.xz * d_lll.z.zz * tmp189 + 2. * K_ll.xz * a_l.z * gamma_uu.xz - 4. * K_ll.xy * d_lll.z.xx * tmp43 - 4. * K_ll.xy * d_lll.z.xy * tmp86 - 4. * K_ll.xy * d_lll.z.xy * tmp55 - 4. * K_ll.xy * d_lll.z.xz * tmp98 - 4. * K_ll.xy * d_lll.z.xz * tmp61 - 4. * K_ll.xy * d_lll.z.yy * tmp110 - 4. * K_ll.xy * d_lll.z.yz * tmp116 - 4. * K_ll.xy * d_lll.z.yz * tmp122 - 4. * K_ll.xy * d_lll.z.zz * tmp128 + 2. * K_ll.xy * a_l.z * gamma_uu.xy - 2. * K_ll.xx * d_lll.z.xx * tmp37 - 4. * K_ll.xx * d_lll.z.xy * tmp43 - 4. * K_ll.xx * d_lll.z.xz * tmp49 - 2. * K_ll.xx * d_lll.z.yy * tmp55 - 4. * K_ll.xx * d_lll.z.yz * tmp61 - 2. * K_ll.xx * d_lll.z.zz * tmp67 + K_ll.xx * a_l.z * gamma_uu.xx);
	deriv->dDelta_lll.x.xx += -K_ll.xx * tmp1049;
	deriv->dDelta_lll.x.xy += -K_ll.xy * tmp1049;
	deriv->dDelta_lll.x.xz += -K_ll.xz * tmp1049;
	deriv->dDelta_lll.x.yy += -K_ll.yy * tmp1049;
	deriv->dDelta_lll.x.yz += -K_ll.yz * tmp1049;
	deriv->dDelta_lll.x.zz += -K_ll.zz * tmp1049;
	deriv->dDelta_lll.y.xx += -K_ll.xx * tmp1061;
	deriv->dDelta_lll.y.xy += -K_ll.xy * tmp1061;
	deriv->dDelta_lll.y.xz += -K_ll.xz * tmp1061;
	deriv->dDelta_lll.y.yy += -K_ll.yy * tmp1061;
	deriv->dDelta_lll.y.yz += -K_ll.yz * tmp1061;
	deriv->dDelta_lll.y.zz += -K_ll.zz * tmp1061;
	deriv->dDelta_lll.z.xx += -K_ll.xx * tmp1073;
	deriv->dDelta_lll.z.xy += -K_ll.xy * tmp1073;
	deriv->dDelta_lll.z.xz += -K_ll.xz * tmp1073;
	deriv->dDelta_lll.z.yy += -K_ll.yy * tmp1073;
	deriv->dDelta_lll.z.yz += -K_ll.yz * tmp1073;
	deriv->dDelta_lll.z.zz += -K_ll.zz * tmp1073;
	deriv->K_ll.xx += -alpha * (2. * tmp1549 + 2. * tmp1532 - 2. * tmp1535 - 2. * tmp1540 - 2. * tmp1545 + gamma_uu.xx * K_ll.xx * K_ll.xx - gamma_uu.xx * gamma_uu.yy * tmp1523 - gamma_uu.xx * gamma_uu.zz * tmp1527 + 2. * tmp1521 + d_lll.z.xx * tmp1326 - 2. * d_lll.z.xx * tmp1331 - 2. * d_lll.z.xx * tmp1334 - d_lll.z.xx * tmp1337 + 2. * tmp1479 - d_lll.y.yy * d_lll.z.xx * tmp224 - 2. * d_lll.y.yz * d_lll.z.xx * tmp279 - d_lll.y.zz * d_lll.z.xx * tmp291 - 2. * d_lll.z.xx * tmp1494 - 2. * d_lll.z.xx * tmp1323 + 4. * d_lll.y.xz * tmp1470 - 4. * d_lll.y.xz * tmp1473 + 2. * d_lll.y.xz * tmp1375 - 4. * d_lll.y.xz * tmp1465 + 2. * d_lll.y.xx * tmp1436 - 2. * d_lll.y.xx * tmp1278 - d_lll.y.xx * tmp1281 - 2. * d_lll.y.xx * tmp1284 - d_lll.y.xx * tmp1287 - 2. * d_lll.y.xy * tmp1359 + d_lll.y.xx * tmp1267 - 2. * d_lll.y.xx * tmp1272 - 2. * d_lll.y.xx * tmp1426 - 4. * d_lll.y.xx * tmp1275 + d_lll.x.zz * d_lll.y.xx * tmp177 - 2. * d_lll.x.zz * d_lll.y.xx * tmp128 - d_lll.x.zz * d_lll.z.xx * tmp189 - tmp1398 - 2. * d_lll.y.xx * tmp1255 - 2. * d_lll.y.xx * tmp1405 - d_lll.y.xx * tmp1261 - 2. * d_lll.y.xx * tmp1264 + d_lll.x.yy * tmp1359 - tmp1363 - 4. * d_lll.x.yz * tmp1365 - 2. * d_lll.x.yz * d_lll.y.xx * tmp122 - 2. * d_lll.x.yz * tmp1375 - 2. * tmp1381 + 2. * d_lll.x.xz * tmp1337 - 4. * d_lll.x.yy * tmp1340 - 2. * d_lll.x.yy * tmp1345 - d_lll.x.yy * d_lll.y.xx * tmp110 - 2. * d_lll.x.yy * d_lll.z.xx * tmp116 + 4. * d_lll.x.xz * tmp1334 + 4. * d_lll.x.xz * tmp1331 + 4. * d_lll.x.xz * tmp1323 - 2. * d_lll.x.xz * tmp1326 + 4. * d_lll.x.xz * tmp1320 + 2. * d_lll.x.xz * tmp1317 + 4. * d_lll.x.xz * tmp1314 + 2. * d_lll.x.xz * tmp1311 + 4. * d_lll.x.xz * tmp1308 + 4. * d_lll.x.xz * tmp1305 + 2. * d_lll.x.xy * tmp1287 - 2. * d_lll.x.xz * tmp1290 - 4. * d_lll.x.xz * tmp1295 - 2. * d_lll.x.xz * tmp1300 + 4. * d_lll.x.xy * tmp1284 + 2. * d_lll.x.xy * tmp1281 + 4. * d_lll.x.xy * tmp1278 + 4. * d_lll.x.xy * tmp1275 + 4. * d_lll.x.xy * tmp1272 + 4. * d_lll.x.xy * tmp1264 - 2. * d_lll.x.xy * tmp1267 + 2. * d_lll.x.xy * tmp1261 + 4. * d_lll.x.xy * tmp1258 + 4. * d_lll.x.xy * tmp1255 + d_lll.x.xx * tmp1238 - 2. * d_lll.x.xy * tmp1240 - 4. * d_lll.x.xy * tmp1245 - 2. * d_lll.x.xy * tmp1250 + 2. * d_lll.x.xx * tmp1235 + 2. * d_lll.x.xx * tmp1228 - d_lll.x.xx * tmp1231 + 2. * d_lll.x.xx * tmp1225 + 2. * d_lll.x.xx * tmp1222 + 2. * d_lll.x.xx * tmp1219 + 2. * d_lll.x.xx * tmp1212 - d_lll.x.xx * tmp1215 + d_lll.x.xx * tmp1210 + 2. * d_lll.x.xx * tmp1207 + 2. * d_lll.x.xx * tmp1204 + a_l.x * a_l.x - d_lll.x.xx * tmp1191 - 2. * d_lll.x.xx * tmp1195 - d_lll.x.xx * tmp1200 + a_l.z * tmp1147 + a_l.z * tmp1143 + a_l.y * tmp1133 - a_l.z * tmp1137 - 2. * a_l.z * tmp1139 - 2. * a_l.z * tmp1141 + a_l.y * tmp1129 + a_l.x * tmp1119 - a_l.y * tmp1123 - 2. * a_l.y * tmp1125 - 2. * a_l.y * tmp1127 + a_l.x * tmp1115 + 4. * Z_l.z * tmp1141 - 2. * Z_l.z * tmp1143 - 2. * Z_l.z * tmp1147 - a_l.x * tmp1109 - 2. * a_l.x * tmp1111 - 2. * a_l.x * tmp1113 + 4. * Z_l.z * tmp1139 + 2. * Z_l.z * tmp1137 + 4. * Z_l.y * tmp1127 - 2. * Z_l.y * tmp1129 - 2. * Z_l.y * tmp1133 + 4. * Z_l.y * tmp1125 + 2. * Z_l.y * tmp1123 + 4. * Z_l.x * tmp1113 - 2. * Z_l.x * tmp1115 - 2. * Z_l.x * tmp1119 + 4. * Z_l.x * tmp1111 + 2. * Z_l.x * tmp1109 + 4. * M_PI * gamma_ll.xx * rho + 8. * M_PI * S_ll.xx + 4. * K_ll.xy * K_ll.xz * gamma_uu.yz - 4. * M_PI * S * gamma_ll.xx + 2. * K_ll.xx * Theta + 2. * K_ll.xx * tmp2 - K_ll.xx * tmp8 - 2. * K_ll.xx * tmp3 - K_ll.xx * tmp4 + 2. * K_ll.xx * tmp1);
	deriv->K_ll.xy += -alpha * (2. * gamma_uu.xy * tmp1539 - 2. * tmp2162 + 2. * gamma_uu.xy * tmp1534 + gamma_uu.xy * tmp2153 + gamma_uu.xx * gamma_uu.xy * tmp1523 + 2. * d_lll.z.xx * tmp1888 - 2. * d_lll.z.xy * tmp1323 - d_lll.z.xy * tmp1326 - 2. * d_lll.z.xy * tmp1334 - d_lll.z.xy * tmp1337 + 2. * d_lll.y.yz * tmp1375 - 2. * d_lll.y.yz * tmp1465 - 2. * d_lll.y.yz * tmp1473 - d_lll.y.zz * d_lll.z.xy * tmp291 - d_lll.z.xx * tmp2121 - 2. * d_lll.z.xx * tmp1885 + d_lll.y.xz * tmp1337 - d_lll.y.yy * tmp2100 + 2. * d_lll.y.xz * tmp1334 + d_lll.y.xz * tmp1326 + 2. * d_lll.y.xz * tmp1323 + 2. * d_lll.y.xz * tmp1320 + 2. * d_lll.y.xz * tmp2079 - 2. * d_lll.y.xz * tmp1494 + d_lll.y.xz * tmp1317 - d_lll.y.xz * tmp1961 + 2. * d_lll.y.xz * tmp1952 + d_lll.y.xz * tmp1311 + 2. * d_lll.y.xy * tmp2060 - 2. * d_lll.y.xy * tmp1275 + 2. * d_lll.y.xy * tmp1258 - 2. * d_lll.y.xy * tmp1426 + d_lll.y.xx * tmp1238 + 2. * d_lll.y.xx * tmp1235 + d_lll.y.xx * tmp1231 + 2. * d_lll.y.xx * tmp1225 + d_lll.y.xx * tmp1222 + d_lll.y.xx * tmp2038 + 2. * d_lll.y.xx * tmp1219 + 2. * d_lll.y.xx * tmp1826 - d_lll.y.xx * tmp1215 + d_lll.y.xx * tmp1210 + 2. * d_lll.y.xx * tmp1816 + d_lll.y.xx * tmp1207 + 2. * d_lll.y.xx * tmp1811 + d_lll.x.zz * tmp1990 - 2. * d_lll.x.zz * tmp1992 - d_lll.x.zz * tmp1997 - d_lll.x.zz * tmp2001 - 2. * d_lll.x.zz * tmp2005 - d_lll.x.zz * tmp2010 - d_lll.x.zz * d_lll.z.xy * tmp189 + d_lll.x.yz * tmp1337 - d_lll.x.zz * d_lll.y.xx * tmp159 + 2. * d_lll.x.yz * tmp1334 + 2. * d_lll.x.yz * tmp1331 + 2. * d_lll.x.yz * tmp1323 - d_lll.x.yz * tmp1326 + 2. * d_lll.x.yz * tmp1320 + d_lll.x.yz * tmp1961 - 2. * d_lll.x.yz * tmp1494 + 2. * d_lll.x.yz * tmp1940 - 2. * d_lll.x.yz * tmp1943 - d_lll.x.yz * tmp1311 - 2. * d_lll.x.yz * tmp1952 - d_lll.x.yz * tmp1317 + d_lll.x.yz * tmp1300 - d_lll.x.yz * tmp1936 + d_lll.x.yy * tmp1287 + 2. * d_lll.x.yy * tmp1284 + d_lll.x.yy * tmp1281 + 2. * d_lll.x.yy * tmp1278 + d_lll.x.yy * tmp1436 + d_lll.x.yy * tmp1426 + d_lll.x.yy * tmp1272 + 2. * d_lll.x.yy * tmp1258 - d_lll.x.yy * tmp1405 - d_lll.x.yy * tmp1267 + d_lll.x.yy * tmp1905 + 2. * d_lll.x.yy * tmp1902 + d_lll.x.yy * tmp1896 - d_lll.x.yy * tmp1250 + 2. * d_lll.x.yy * tmp1245 + 2. * d_lll.x.xz * tmp1885 - 2. * d_lll.x.xz * tmp1888 + 2. * d_lll.x.xz * d_lll.x.yz * tmp67 - 4. * d_lll.x.xz * tmp1855 - 2. * d_lll.x.xz * tmp1860 - 2. * d_lll.x.xz * tmp1865 - 4. * d_lll.x.xz * tmp1870 - 2. * d_lll.x.xz * tmp1875 - 2. * d_lll.x.xz * tmp1880 + 2. * d_lll.x.xz * d_lll.x.yy * tmp61 + 2. * d_lll.x.xy * tmp1228 - 2. * d_lll.x.xy * tmp1231 + 2. * d_lll.x.xy * tmp1808 - 4. * d_lll.x.xy * tmp1811 - 2. * d_lll.x.xy * tmp1816 - 2. * d_lll.x.xy * tmp1210 - 4. * d_lll.x.xy * tmp1826 - 2. * d_lll.x.xy * tmp1219 - 2. * d_lll.x.xy * tmp1836 + 2. * d_lll.x.xy * tmp1805 + d_lll.x.xx * tmp1777 - 2. * d_lll.x.xx * tmp1779 - d_lll.x.xx * tmp1784 - d_lll.x.xx * tmp1788 - 2. * d_lll.x.xx * tmp1792 - d_lll.x.xx * tmp1797 - d_lll.x.xx * tmp1801 + d_lll.x.xx * tmp1775 + a_l.z * tmp1731 + a_l.y * tmp1719 - a_l.z * tmp1723 - a_l.z * tmp1725 - a_l.z * tmp1727 - a_l.z * tmp1729 + a_l.x * tmp1707 - a_l.y * tmp1711 - a_l.y * tmp1713 - a_l.y * tmp1115 - a_l.y * tmp1717 + a_l.x * a_l.y - a_l.x * tmp1699 - a_l.x * tmp1701 - a_l.x * tmp1703 - a_l.x * tmp1705 + 2. * Z_l.z * tmp1729 - 2. * Z_l.z * tmp1731 + 2. * Z_l.z * tmp1727 + 2. * Z_l.z * tmp1725 + 2. * Z_l.z * tmp1723 + 2. * Z_l.y * tmp1717 - 2. * Z_l.y * tmp1719 + 2. * Z_l.y * tmp1115 + 2. * Z_l.y * tmp1713 + 2. * Z_l.y * tmp1711 + 2. * Z_l.x * tmp1705 - 2. * Z_l.x * tmp1707 + 2. * Z_l.x * tmp1703 + 2. * Z_l.x * tmp1701 + 2. * Z_l.x * tmp1699 + 4. * M_PI * gamma_ll.xy * rho + 8. * M_PI * S_ll.xy + 2. * K_ll.xz * K_ll.yz * gamma_uu.zz - 4. * M_PI * S * gamma_ll.xy + 2. * K_ll.xz * K_ll.yy * gamma_uu.yz + 2. * K_ll.xy * Theta + K_ll.xy * tmp8 - K_ll.xy * tmp4 + 2. * K_ll.xx * tmp1681 + 2. * K_ll.xx * tmp1679 + K_ll.xx * K_ll.xy * gamma_uu.xx);
	deriv->K_ll.xz += -alpha * (gamma_uu.xz * tmp2817 + 2. * gamma_uu.xz * tmp2814 + 2. * gamma_uu.xz * tmp2811 + gamma_uu.xx * gamma_uu.xz * tmp1527 - 2. * gamma_uu.xy * tmp2806 + d_lll.z.xy * tmp1287 + 2. * d_lll.z.xy * tmp2633 + d_lll.z.xy * tmp1281 + 2. * d_lll.z.xy * tmp2624 + d_lll.z.xx * tmp1238 + 2. * d_lll.z.xx * tmp2536 + 2. * d_lll.z.xx * tmp1228 - d_lll.z.xx * tmp1231 + 2. * d_lll.z.xx * tmp2526 + 2. * d_lll.z.xx * tmp1836 + d_lll.z.xx * tmp1222 + d_lll.y.zz * tmp1470 + d_lll.y.zz * tmp1375 + 2. * d_lll.y.yz * tmp2100 + 2. * d_lll.y.yz * tmp1359 + d_lll.y.yy * d_lll.z.xy * tmp218 + d_lll.y.yy * d_lll.z.xx * tmp110 + 2. * d_lll.y.xz * tmp1275 - 2. * d_lll.y.xz * tmp1436 - 2. * d_lll.y.xz * tmp2624 - d_lll.y.xz * tmp1281 - 2. * d_lll.y.xz * tmp2633 - d_lll.y.xz * tmp1287 + d_lll.y.xz * tmp1426 + 2. * d_lll.y.xy * tmp2565 - d_lll.y.xz * tmp1261 - 2. * d_lll.y.xz * tmp1264 - d_lll.y.xz * tmp1267 + 2. * d_lll.y.xy * tmp2559 + 2. * d_lll.y.xx * tmp2492 - 2. * d_lll.y.xy * tmp2555 + 2. * d_lll.y.xx * tmp2482 - 2. * d_lll.y.xx * tmp2700 + 2. * d_lll.y.xx * tmp2477 - 2. * d_lll.y.xx * tmp2692 + d_lll.y.xx * tmp2683 - d_lll.y.xx * tmp2685 + 2. * d_lll.y.xx * tmp2469 - 2. * d_lll.y.xx * tmp2474 + d_lll.x.zz * tmp1331 - d_lll.y.xx * tmp2671 + 2. * d_lll.x.zz * tmp1320 - d_lll.x.zz * tmp1326 + d_lll.x.zz * tmp2079 - d_lll.x.zz * tmp1494 + d_lll.x.zz * tmp1317 + 2. * d_lll.x.zz * tmp1314 + d_lll.x.zz * tmp1311 + d_lll.x.zz * tmp2647 + 2. * d_lll.x.zz * tmp1305 + d_lll.x.zz * tmp1936 + 2. * d_lll.x.yz * tmp2060 - 2. * d_lll.x.yz * tmp2624 - d_lll.x.yz * tmp1281 - 2. * d_lll.x.yz * tmp2633 - d_lll.x.yz * tmp1287 + 2. * d_lll.x.yz * tmp1272 - d_lll.x.yz * tmp1426 + 2. * d_lll.x.yz * tmp1264 - d_lll.x.yz * tmp1267 + d_lll.x.yz * tmp1261 + 2. * d_lll.x.yz * tmp1258 - 2. * d_lll.x.yz * tmp1405 + 2. * d_lll.x.yz * tmp1255 + d_lll.x.yz * tmp2592 + 2. * d_lll.x.yz * tmp1902 + d_lll.x.yz * tmp1250 + d_lll.x.yy * tmp2563 - d_lll.x.yy * tmp2565 - 2. * d_lll.x.yy * tmp2569 - d_lll.x.yy * tmp2574 - 2. * d_lll.x.yy * tmp2578 - d_lll.x.yy * tmp2583 + 2. * d_lll.x.yy * d_lll.x.zz * tmp116 - d_lll.x.yy * d_lll.x.zz * tmp122 - d_lll.x.yy * tmp2555 - d_lll.x.yy * tmp2559 + d_lll.x.yy * d_lll.x.yz * tmp110 + 2. * d_lll.x.xz * tmp1219 - 2. * d_lll.x.xz * tmp1836 - 4. * d_lll.x.xz * tmp2526 - 2. * d_lll.x.xz * tmp1228 - 4. * d_lll.x.xz * tmp2536 - 2. * d_lll.x.xz * tmp1238 + 2. * d_lll.x.xz * tmp2505 - 2. * d_lll.x.xz * tmp1816 - 2. * d_lll.x.xz * tmp1215 + 2. * d_lll.x.xz * tmp1808 + 2. * d_lll.x.xy * tmp2474 - 2. * d_lll.x.xy * tmp2477 - 4. * d_lll.x.xy * tmp2482 - 2. * d_lll.x.xy * tmp2487 - 4. * d_lll.x.xy * tmp2492 - 2. * d_lll.x.xy * tmp2497 + 2. * d_lll.x.xy * tmp2461 - 2. * d_lll.x.xy * tmp2464 - 2. * d_lll.x.xy * tmp2469 + 2. * d_lll.x.xy * d_lll.x.yz * tmp55 + d_lll.x.xx * tmp2430 - d_lll.x.xx * tmp2432 - d_lll.x.xx * tmp2436 - 2. * d_lll.x.xx * tmp2440 - d_lll.x.xx * tmp2445 - 2. * d_lll.x.xx * tmp2449 - d_lll.x.xx * tmp2454 + d_lll.x.xx * tmp2428 + a_l.z * tmp1717 - a_l.z * tmp1119 - a_l.z * tmp1719 + a_l.y * tmp2368 - a_l.y * tmp2372 - a_l.y * tmp2374 - a_l.z * tmp1713 - a_l.z * tmp2378 + a_l.x * tmp2356 - a_l.x * tmp2360 - a_l.x * tmp2362 - a_l.y * tmp2364 - a_l.y * tmp2366 + a_l.x * a_l.z - a_l.x * tmp2352 - a_l.x * tmp2354 + 2. * Z_l.z * tmp1719 + 2. * Z_l.z * tmp1119 + 2. * Z_l.z * tmp2378 - 2. * Z_l.z * tmp1717 + 2. * Z_l.z * tmp1713 + 2. * Z_l.y * tmp2374 + 2. * Z_l.y * tmp2372 + 2. * Z_l.y * tmp2366 - 2. * Z_l.y * tmp2368 + 2. * Z_l.y * tmp2364 + 2. * Z_l.x * tmp2362 + 2. * Z_l.x * tmp2360 + 2. * Z_l.x * tmp2354 - 2. * Z_l.x * tmp2356 + 2. * Z_l.x * tmp2352 + 4. * M_PI * gamma_ll.xz * rho + 8. * M_PI * S_ll.xz + 2. * K_ll.xz * Theta - 4. * M_PI * S * gamma_ll.xz + K_ll.xz * tmp4 + 2. * K_ll.xy * tmp2338 - K_ll.xz * tmp8 + 2. * K_ll.xy * tmp2336 + 2. * K_ll.xx * tmp2334 + 2. * K_ll.xx * tmp2332 + K_ll.xx * tmp2331);
	deriv->K_ll.yy += alpha * (tmp3479 + gamma_uu.yy * gamma_uu.zz * tmp3429 - 2. * tmp3432 + 2. * gamma_uu.xx * tmp1544 - gamma_uu.yy * K_ll.yy * K_ll.yy + 2. * gamma_uu.xx * tmp1539 + 2. * gamma_uu.xx * tmp1534 + gamma_uu.xx * tmp2153 + d_lll.z.yy * tmp1337 - 2. * gamma_uu.xx * tmp1531 + 2. * d_lll.z.yy * tmp1334 + 2. * d_lll.z.xz * tmp3255 + 2. * d_lll.z.xy * tmp1885 - 2. * tmp3400 + 2. * d_lll.z.xx * d_lll.z.yy * tmp67 + tmp3390 - d_lll.z.xx * tmp3153 + d_lll.y.zz * d_lll.z.yy * tmp291 + 2. * d_lll.y.yz * tmp1961 - 4. * d_lll.y.yz * tmp2079 - 4. * d_lll.y.yz * tmp1320 - 4. * d_lll.y.yz * tmp1323 - 4. * d_lll.y.yz * tmp1334 - 2. * d_lll.y.yz * tmp1337 + 2. * d_lll.y.yz * tmp1317 + d_lll.y.yy * tmp1426 - 2. * d_lll.y.yy * tmp2060 - 2. * d_lll.y.yy * tmp1436 - 2. * d_lll.y.yy * tmp1278 - 2. * d_lll.y.yy * tmp1284 - d_lll.y.yy * tmp1287 + d_lll.y.yy * tmp1267 + 2. * tmp3324 + 2. * d_lll.y.xz * tmp1885 + 4. * d_lll.y.xz * tmp1875 + 4. * d_lll.y.xz * tmp3209 + 2. * d_lll.y.xz * tmp3204 + 2. * d_lll.y.xy * tmp1215 - 2. * d_lll.y.xy * tmp2038 - 4. * d_lll.y.xy * tmp1836 - 4. * d_lll.y.xy * tmp1225 - 4. * d_lll.y.xy * tmp1235 - 2. * d_lll.y.xy * tmp1238 + 4. * d_lll.y.xy * tmp1816 + 2. * d_lll.y.xx * tmp3275 + 2. * d_lll.y.xx * tmp1797 - d_lll.y.xx * tmp3132 + 2. * d_lll.y.xx * tmp3127 + d_lll.y.xx * tmp3122 + 4. * d_lll.y.xx * tmp1784 + 2. * d_lll.y.xx * tmp1779 + d_lll.x.zz * tmp3255 + d_lll.x.zz * d_lll.y.yy * tmp177 - 2. * d_lll.x.zz * d_lll.y.yy * tmp128 - 2. * d_lll.x.zz * d_lll.y.yz * tmp189 + 2. * d_lll.x.zz * d_lll.y.xy * tmp159 - 4. * d_lll.x.zz * d_lll.y.xy * tmp67 + 4. * d_lll.x.yz * tmp1888 - 2. * tmp3231 + 4. * d_lll.x.yz * tmp1880 - 2. * d_lll.x.yz * tmp1885 + d_lll.x.yy * tmp1238 - 4. * d_lll.x.yz * tmp1855 - 2. * d_lll.x.yz * tmp3204 - 4. * d_lll.x.yz * tmp3209 - 4. * d_lll.x.yz * tmp2121 + 2. * d_lll.x.yy * tmp1235 + 2. * d_lll.x.yy * tmp1231 + 2. * d_lll.x.yy * tmp1225 + 4. * d_lll.x.yy * tmp1836 + d_lll.x.yy * tmp2038 - 2. * d_lll.x.yy * tmp1222 + 2. * d_lll.x.yy * tmp1219 + 2. * d_lll.x.yy * tmp1207 - d_lll.x.yy * tmp1215 + d_lll.x.yy * tmp3166 + 2. * d_lll.x.yy * tmp2505 + 2. * d_lll.x.yy * tmp1195 - d_lll.x.yy * tmp1200 + 2. * d_lll.x.xz * tmp3153 + 2. * d_lll.x.xz * d_lll.x.yy * tmp49 - 4. * d_lll.x.xz * d_lll.y.xy * tmp49 - 2. * d_lll.x.xz * d_lll.y.yy * tmp98 - 4. * d_lll.x.xz * d_lll.y.yz * tmp159 + 2. * d_lll.x.xy * tmp3132 + 2. * d_lll.x.xy * tmp1775 - 4. * d_lll.x.xy * tmp1779 - 2. * d_lll.x.xy * tmp3122 - 4. * d_lll.x.xy * tmp3127 + d_lll.x.xx * d_lll.z.yy * tmp49 + d_lll.x.xx * d_lll.x.yy * tmp37 - 2. * d_lll.x.xx * d_lll.y.xy * tmp37 - d_lll.x.xx * d_lll.y.yy * tmp43 - 2. * d_lll.x.xx * d_lll.y.yz * tmp49 + 2. * a_l.z * tmp3054 - a_l.z * tmp3058 - a_l.y * a_l.y + a_l.z * tmp3050 + 2. * a_l.z * tmp3046 + 2. * a_l.y * tmp3038 - a_l.y * tmp3042 - a_l.z * tmp3044 + a_l.y * tmp3034 + 2. * a_l.y * tmp3030 + 2. * a_l.x * tmp3022 - a_l.x * tmp3026 - a_l.y * tmp1699 + a_l.x * tmp3018 + 2. * a_l.x * tmp3014 + 2. * Z_l.z * tmp3058 - a_l.x * tmp3012 + 2. * Z_l.z * tmp3044 - 4. * Z_l.z * tmp3046 - 2. * Z_l.z * tmp3050 - 4. * Z_l.z * tmp3054 + 2. * Z_l.y * tmp3042 + 2. * Z_l.y * tmp1699 - 4. * Z_l.y * tmp3030 - 2. * Z_l.y * tmp3034 - 4. * Z_l.y * tmp3038 + 2. * Z_l.x * tmp3026 + 2. * Z_l.x * tmp3012 - 4. * Z_l.x * tmp3014 - 2. * Z_l.x * tmp3018 - 4. * Z_l.x * tmp3022 + 4. * M_PI * S * gamma_ll.yy - 8. * M_PI * S_ll.yy - 4. * M_PI * gamma_ll.yy * rho + K_ll.yy * tmp4 - 2. * K_ll.yy * Theta + 2. * K_ll.xz * K_ll.yy * gamma_uu.xz - 2. * K_ll.yy * tmp3 + K_ll.xx * K_ll.yy * gamma_uu.xx - 2. * K_ll.xy * tmp1679 - 4. * K_ll.xy * tmp1681);
	deriv->K_ll.yz += alpha * (2. * gamma_uu.xy * gamma_uu.xz * tmp1379 - gamma_uu.yy * gamma_uu.yz * tmp3429 - gamma_uu.yz * tmp4086 + d_lll.z.xx * tmp3132 - 2. * d_lll.z.xx * tmp3275 - 2. * d_lll.z.xy * tmp2526 - 2. * d_lll.z.xy * tmp1228 - d_lll.z.xy * tmp1231 - 2. * d_lll.z.xy * tmp2536 - d_lll.z.xy * tmp1238 - 2. * d_lll.z.xz * tmp1888 - 2. * d_lll.z.yy * tmp2633 - d_lll.z.yy * tmp1287 - 2. * gamma_uu.xx * tmp2806 - 2. * gamma_uu.xx * tmp2161 + d_lll.y.zz * tmp1494 - 2. * d_lll.y.zz * tmp1320 - d_lll.y.zz * tmp1331 - d_lll.z.xx * tmp1801 + d_lll.y.zz * tmp1961 - d_lll.y.zz * tmp2079 + 2. * d_lll.y.yz * tmp1287 + 4. * d_lll.y.yz * tmp2633 + 4. * d_lll.y.yz * tmp2624 + 2. * d_lll.y.yz * tmp1275 + 2. * d_lll.y.yz * tmp2060 + d_lll.y.yy * tmp2583 - 2. * d_lll.y.yz * tmp1272 + 2. * d_lll.y.yy * tmp2578 + 2. * d_lll.y.yy * tmp2569 + d_lll.y.yy * tmp2565 + d_lll.y.yy * tmp2563 + d_lll.y.xz * tmp1238 - d_lll.y.yy * tmp3968 + 2. * d_lll.y.xz * tmp2536 + d_lll.y.xz * tmp1231 + 2. * d_lll.y.xz * tmp2526 - 2. * d_lll.y.xz * tmp1228 + d_lll.y.xz * tmp2038 + 2. * d_lll.y.xy * tmp2497 - d_lll.y.xz * tmp1210 - 2. * d_lll.y.xz * tmp1826 - d_lll.y.xz * tmp1215 - 2. * d_lll.y.xz * tmp1219 + 4. * d_lll.y.xy * tmp2700 + 4. * d_lll.y.xy * tmp2692 + 2. * d_lll.y.xy * tmp2477 + 2. * d_lll.y.xy * tmp2683 + d_lll.y.xx * tmp2454 - 2. * d_lll.y.xy * tmp2464 - 2. * d_lll.y.xy * tmp2469 + 2. * d_lll.y.xx * tmp2449 + d_lll.y.xx * tmp3725 - d_lll.y.xx * tmp2445 + 2. * d_lll.y.xx * tmp2440 + d_lll.y.xx * tmp2436 + d_lll.y.xx * tmp3890 + d_lll.y.xx * tmp3715 - 2. * d_lll.y.xx * d_lll.y.zz * tmp61 + 2. * d_lll.x.zz * tmp1870 - 2. * d_lll.x.zz * tmp3209 - d_lll.x.zz * tmp1875 - d_lll.x.zz * tmp2121 - d_lll.x.zz * tmp1885 - d_lll.y.xx * tmp2432 + d_lll.x.zz * tmp3852 - 2. * d_lll.x.zz * tmp1860 + 2. * d_lll.x.zz * tmp1855 + d_lll.x.yz * tmp1238 - 2. * d_lll.x.zz * tmp3844 + 2. * d_lll.x.yz * tmp2536 + 2. * d_lll.x.yz * tmp2526 - d_lll.x.yz * tmp1231 + 2. * d_lll.x.yz * tmp1222 - 2. * d_lll.x.yz * tmp1836 + d_lll.x.yz * tmp2038 + 2. * d_lll.x.yz * tmp1826 - d_lll.x.yz * tmp1215 + d_lll.x.yz * tmp1210 + 2. * d_lll.x.yz * tmp1207 - 2. * d_lll.x.yz * tmp1816 + 2. * d_lll.x.yz * tmp1811 + d_lll.x.yz * tmp3166 + d_lll.x.yz * tmp1200 + 2. * d_lll.x.yy * tmp2492 + 2. * d_lll.x.yy * tmp2482 - d_lll.x.yy * tmp2487 - 2. * d_lll.x.yy * tmp2700 + d_lll.x.yy * tmp2685 - 2. * d_lll.x.yy * tmp2477 - 2. * d_lll.x.yy * tmp2692 + 2. * d_lll.x.yy * tmp3755 - 2. * d_lll.x.yy * tmp2461 - d_lll.x.yy * tmp2671 - d_lll.x.yy * tmp2474 + d_lll.x.yy * d_lll.x.yz * tmp86 + 2. * d_lll.x.xz * tmp1777 - 2. * d_lll.x.xz * tmp1784 - 2. * d_lll.x.xz * d_lll.y.zz * tmp159 - 2. * d_lll.x.xz * tmp1801 - 2. * d_lll.x.xz * tmp3132 + 2. * d_lll.x.xy * tmp2428 - 2. * d_lll.x.xy * tmp2432 - 2. * d_lll.x.xy * tmp3715 - 2. * d_lll.x.xy * tmp2436 - 2. * d_lll.x.xy * tmp3725 + d_lll.x.xx * d_lll.x.yz * tmp37 - d_lll.x.xx * d_lll.y.xz * tmp37 - d_lll.x.xx * d_lll.y.zz * tmp49 - d_lll.x.xx * d_lll.z.xy * tmp37 - d_lll.x.xx * d_lll.z.yy * tmp43 + a_l.z * tmp3042 + a_l.z * tmp1707 + a_l.z * tmp3653 + a_l.z * tmp1705 + a_l.y * tmp3643 - a_l.z * tmp1701 + a_l.y * tmp2362 + a_l.y * tmp3635 + a_l.y * tmp2356 + a_l.x * tmp3625 - a_l.y * a_l.z - a_l.y * tmp2352 + a_l.x * tmp3621 + a_l.x * tmp3617 + a_l.x * tmp3613 + 2. * Z_l.z * tmp1701 - 2. * Z_l.z * tmp1705 - 2. * Z_l.z * tmp3653 - 2. * Z_l.z * tmp1707 - 2. * Z_l.z * tmp3042 - a_l.x * tmp3611 + 2. * Z_l.y * tmp2352 - 2. * Z_l.y * tmp2356 - 2. * Z_l.y * tmp3635 - 2. * Z_l.y * tmp2362 - 2. * Z_l.y * tmp3643 + 2. * Z_l.x * tmp3611 - 2. * Z_l.x * tmp3613 - 2. * Z_l.x * tmp3617 - 2. * Z_l.x * tmp3621 - 2. * Z_l.x * tmp3625 + 4. * M_PI * S * gamma_ll.yz - 8. * M_PI * S_ll.yz - 4. * M_PI * gamma_ll.yz * rho + K_ll.xx * K_ll.yz * gamma_uu.xx - 2. * K_ll.xy * tmp2331 - 2. * K_ll.xy * tmp2334 - 2. * K_ll.xz * tmp1679 - K_ll.yy * tmp2336 - 2. * K_ll.yy * tmp2338 - K_ll.yz * tmp4 - 2. * K_ll.yz * Theta);
	deriv->K_ll.zz += alpha * (tmp4705 + gamma_uu.yy * tmp4086 - gamma_uu.zz * K_ll.zz * K_ll.zz + gamma_uu.xx * tmp2817 - 2. * gamma_uu.yy * tmp3431 + 2. * gamma_uu.xx * tmp2814 + 2. * gamma_uu.xx * gamma_uu.yy * tmp1477 + 2. * gamma_uu.xx * tmp2811 + tmp4657 - 2. * gamma_uu.xx * tmp1548 + d_lll.z.yy * tmp4598 + 2. * d_lll.z.yy * tmp2578 + 2. * d_lll.z.xz * tmp1231 + 2. * tmp4644 + 2. * d_lll.z.xy * tmp4458 + 4. * d_lll.z.xy * tmp2700 + 4. * d_lll.z.xy * tmp2487 + 4. * d_lll.z.xy * tmp2482 + d_lll.z.xx * tmp4395 + 2. * d_lll.z.xx * tmp4390 + 2. * d_lll.z.xx * tmp2445 + 2. * d_lll.z.xx * tmp2440 + 4. * d_lll.z.xx * tmp2436 + d_lll.y.zz * tmp1281 + 2. * d_lll.y.zz * tmp1436 + 2. * d_lll.y.zz * tmp2060 + 2. * d_lll.y.yz * tmp3968 - 4. * d_lll.y.yz * d_lll.z.xz * tmp122 - 4. * d_lll.y.yz * tmp2578 - 2. * d_lll.y.yz * tmp4598 - d_lll.y.zz * tmp1426 + d_lll.y.yy * d_lll.y.zz * tmp218 - 2. * d_lll.y.yy * d_lll.z.xz * tmp110 - 2. * d_lll.y.yy * d_lll.z.yz * tmp218 - d_lll.y.yy * d_lll.z.zz * tmp224 + 2. * d_lll.y.xz * tmp2474 - 4. * d_lll.y.xz * tmp2482 - 4. * d_lll.y.xz * tmp2700 - 2. * d_lll.y.xz * tmp4458 - 2. * tmp4565 + 2. * d_lll.y.xy * tmp4407 - 4. * d_lll.y.xy * tmp4409 - 4. * d_lll.y.xy * tmp4417 - 2. * d_lll.y.xy * tmp4427 + d_lll.y.xx * tmp4374 - 2. * d_lll.y.xx * d_lll.z.zz * tmp61 + 2. * d_lll.y.xx * tmp4369 - 4. * d_lll.y.xx * d_lll.z.yz * tmp55 + 2. * d_lll.y.xx * d_lll.y.zz * tmp55 - 2. * d_lll.y.xx * tmp4364 + 2. * d_lll.x.zz * tmp1228 - d_lll.x.zz * tmp1231 - d_lll.y.xx * tmp4361 + 2. * d_lll.x.zz * tmp1222 + d_lll.x.zz * tmp2038 + 2. * d_lll.x.zz * tmp1215 + 2. * d_lll.x.zz * tmp1212 + d_lll.x.zz * tmp1210 + 4. * d_lll.x.zz * tmp1816 + 2. * d_lll.x.zz * tmp1204 - 2. * d_lll.x.zz * tmp1207 + d_lll.x.zz * tmp3166 + 4. * d_lll.x.yz * tmp2469 - 2. * d_lll.x.yz * tmp2474 - 4. * d_lll.x.yz * tmp2482 - 4. * d_lll.x.yz * tmp2700 - 2. * d_lll.x.yz * tmp4458 - 2. * tmp4464 + 4. * d_lll.x.yz * tmp2464 + 2. * d_lll.x.yz * tmp3755 - 4. * d_lll.x.yz * tmp2671 + d_lll.x.yy * tmp4427 + 2. * d_lll.x.yy * tmp4409 - 4. * d_lll.x.yy * d_lll.z.xz * tmp55 - 2. * d_lll.x.yy * tmp4417 - 2. * d_lll.x.yy * d_lll.z.zz * tmp116 + d_lll.x.yy * tmp4407 + 2. * d_lll.x.yy * d_lll.x.zz * tmp55 + 2. * d_lll.x.xz * tmp3715 - 4. * d_lll.x.xz * tmp2440 - 4. * d_lll.x.xz * tmp4390 - 2. * d_lll.x.xz * tmp4395 - d_lll.x.yy * d_lll.x.zz * tmp86 + 2. * d_lll.x.xz * tmp2430 + 2. * d_lll.x.xy * tmp4361 - 4. * d_lll.x.xy * tmp4364 - 4. * d_lll.x.xy * tmp4369 - 2. * d_lll.x.xy * tmp4374 + 2. * d_lll.x.xy * d_lll.x.zz * tmp43 + d_lll.x.xx * d_lll.y.zz * tmp43 - 2. * d_lll.x.xx * d_lll.z.xz * tmp37 - 2. * d_lll.x.xx * d_lll.z.yz * tmp43 - d_lll.x.xx * d_lll.z.zz * tmp49 + d_lll.x.xx * d_lll.x.zz * tmp37 + a_l.z * tmp4300 - a_l.z * a_l.z + 2. * a_l.z * tmp4296 + 2. * a_l.z * tmp4292 + a_l.y * tmp4284 - a_l.z * tmp2354 - a_l.z * tmp3635 + 2. * a_l.y * tmp4280 + 2. * a_l.y * tmp4276 + a_l.x * tmp4268 - a_l.y * tmp4272 - a_l.y * tmp4274 + 2. * a_l.x * tmp4264 + 2. * a_l.x * tmp4260 + 2. * Z_l.z * tmp3635 - 4. * Z_l.z * tmp4292 - 4. * Z_l.z * tmp4296 - 2. * Z_l.z * tmp4300 - a_l.x * tmp4256 - a_l.x * tmp4258 + 2. * Z_l.z * tmp2354 + 2. * Z_l.y * tmp4274 - 4. * Z_l.y * tmp4276 - 4. * Z_l.y * tmp4280 - 2. * Z_l.y * tmp4284 + 2. * Z_l.y * tmp4272 + 2. * Z_l.x * tmp4258 - 4. * Z_l.x * tmp4260 - 4. * Z_l.x * tmp4264 - 2. * Z_l.x * tmp4268 + 2. * Z_l.x * tmp4256 + 4. * M_PI * S * gamma_ll.zz - 8. * M_PI * S_ll.zz - 4. * M_PI * gamma_ll.zz * rho + K_ll.yy * K_ll.zz * gamma_uu.yy - 2. * K_ll.yz * tmp2338 - 2. * K_ll.zz * Theta + 2. * K_ll.xy * K_ll.zz * gamma_uu.xy - 4. * K_ll.xz * tmp2332 - 2. * K_ll.xz * tmp2334 + K_ll.xx * K_ll.zz * gamma_uu.xx);
	deriv->Theta += alpha * (gamma_uu.zz * tmp4657 + gamma_uu.zz * tmp4644 + gamma_uu.zz * tmp4705 + gamma_uu.yy * tmp3390 - 3. * gamma_uu.yy * tmp3400 - gamma_uu.yy * tmp3429 * tmp230 - gamma_uu.yy * tmp3432 - 3. * gamma_uu.zz * tmp4464 - gamma_uu.zz * tmp1395 * tmp67 - 3. * gamma_uu.zz * tmp4565 - gamma_uu.zz * tmp3386 * tmp230 + gamma_uu.yy * tmp3324 + gamma_uu.yy * tmp3479 + 2. * gamma_uu.xy * tmp2162 - gamma_uu.yy * tmp1360 * tmp55 - 3. * gamma_uu.yy * tmp3231 + 2. * gamma_uu.xy * gamma_uu.xz * tmp2806 + 2. * gamma_uu.xy * gamma_uu.xz * gamma_uu.yz * tmp1379 + 3. * gamma_uu.xx * tmp1545 - gamma_uu.xx * tmp1549 + 3. * gamma_uu.xx * tmp1540 + 3. * gamma_uu.xx * tmp1535 + gamma_uu.xx * tmp1398 - gamma_uu.xx * tmp1523 * tmp55 - 3. * gamma_uu.xx * tmp1479 - gamma_uu.xx * tmp1527 * tmp67 - 3. * gamma_uu.xx * tmp1521 - gamma_uu.xx * tmp1532 + gamma_uu.xx * tmp1381 + gamma_uu.xx * tmp1363 + d_lll.z.yy * tmp6505 - d_lll.z.yy * tmp6511 + 2. * d_lll.z.yy * tmp6495 - 2. * d_lll.z.yy * tmp6501 + 2. * d_lll.z.xz * tmp6050 - 2. * d_lll.z.xz * tmp6053 + 2. * d_lll.z.xy * tmp5946 - 2. * d_lll.z.xy * tmp5952 + 4. * d_lll.z.xy * tmp5936 - 4. * d_lll.z.xy * tmp5942 + 4. * d_lll.z.xy * d_lll.z.yy * tmp5707 - 2. * d_lll.z.xy * tmp5926 - 2. * d_lll.z.xy * tmp5932 + 4. * d_lll.z.xy * tmp5916 - 4. * d_lll.z.xy * tmp5922 + d_lll.z.xx * tmp5641 - d_lll.z.xx * tmp5647 + 2. * d_lll.z.xx * tmp5631 - 2. * d_lll.z.xx * tmp5637 + 3. * d_lll.z.xx * tmp5625 + 2. * d_lll.z.xx * d_lll.z.yy * tmp5442 + 2. * d_lll.z.xx * tmp5615 - 6. * d_lll.z.xx * tmp5621 + 2. * d_lll.z.xx * tmp5601 - 2. * d_lll.z.xx * tmp5607 - d_lll.z.xx * tmp5611 + 4. * d_lll.z.xx * d_lll.z.xy * tmp5296 - 2. * d_lll.z.xx * tmp5591 - 2. * d_lll.z.xx * tmp5597 + d_lll.y.zz * d_lll.z.yy * tmp6436 - d_lll.y.zz * d_lll.z.yy * tmp6441 + 2. * d_lll.y.zz * d_lll.z.xy * tmp5779 - 2. * d_lll.y.zz * d_lll.z.xy * tmp5788 + d_lll.y.zz * d_lll.z.xx * tmp5494 - d_lll.y.zz * d_lll.z.xx * tmp5500 + 2. * d_lll.y.yz * tmp6511 + 4. * d_lll.y.yz * tmp6501 - 2. * d_lll.y.yz * tmp6505 + 4. * d_lll.y.yz * d_lll.z.xz * tmp5788 - 4. * d_lll.y.yz * tmp6495 + 4. * d_lll.y.yz * d_lll.z.xy * tmp5710 - 4. * d_lll.y.yz * d_lll.z.xy * tmp5775 - 4. * d_lll.y.yz * d_lll.z.xz * tmp5785 + 4. * d_lll.y.yz * d_lll.z.xx * tmp5438 - 4. * d_lll.y.yz * d_lll.z.xx * tmp5442 - 2. * d_lll.y.yz * d_lll.z.xx * tmp5570 + 2. * d_lll.y.yz * d_lll.z.xx * tmp5428 + 2. * d_lll.y.yz * d_lll.y.zz * tmp6436 - 2. * d_lll.y.yz * d_lll.y.zz * tmp6441 + d_lll.y.yy * d_lll.z.zz * tmp6441 + 2. * d_lll.y.yy * d_lll.z.yz * tmp6386 - 2. * d_lll.y.yy * d_lll.z.yz * tmp6391 - d_lll.y.yy * d_lll.z.zz * tmp6436 + 2. * d_lll.y.yy * d_lll.z.xz * tmp5710 + 2. * d_lll.y.yy * d_lll.z.xy * tmp5753 - 2. * d_lll.y.yy * d_lll.z.xy * tmp5759 - 2. * d_lll.y.yy * d_lll.z.xz * tmp5707 + d_lll.y.yy * d_lll.z.xx * tmp5424 + d_lll.y.yy * d_lll.z.xx * tmp5418 - 2. * d_lll.y.yy * d_lll.z.xx * tmp5474 + d_lll.y.yy * d_lll.y.zz * tmp6391 + 2. * d_lll.y.xz * tmp5952 - d_lll.y.yy * d_lll.y.zz * tmp6386 + 4. * d_lll.y.xz * tmp5942 - 2. * d_lll.y.xz * tmp5946 + 2. * d_lll.y.xz * tmp5932 - 4. * d_lll.y.xz * tmp5936 + 4. * d_lll.y.xz * tmp5922 - 2. * d_lll.y.xz * tmp5926 + 2. * d_lll.y.xz * tmp5912 - 4. * d_lll.y.xz * tmp5916 + 2. * d_lll.y.xz * tmp5908 + 2. * d_lll.y.xz * tmp5898 - 4. * d_lll.y.xz * tmp5902 + 2. * d_lll.y.xz * tmp5882 - 2. * d_lll.y.xz * tmp5886 - 2. * d_lll.y.xz * tmp5892 + 4. * d_lll.y.xz * d_lll.y.zz * tmp5785 - 2. * d_lll.y.xz * tmp5876 + 4. * d_lll.y.xz * tmp5866 - 2. * d_lll.y.xz * tmp5872 + 2. * d_lll.y.xz * tmp5856 - 4. * d_lll.y.xz * tmp5862 + 2. * d_lll.y.xy * tmp5789 - 2. * d_lll.y.xz * tmp5852 + 4. * d_lll.y.xy * tmp5770 - 2. * d_lll.y.xy * tmp5786 + 4. * d_lll.y.xy * tmp5738 - 4. * d_lll.y.xy * tmp5764 + 4. * d_lll.y.xy * tmp5728 - 4. * d_lll.y.xy * tmp5734 + 4. * d_lll.y.xy * d_lll.z.xx * tmp5328 - 2. * d_lll.y.xy * tmp5716 - 2. * d_lll.y.xy * tmp5719 - 4. * d_lll.y.xy * tmp5724 + 2. * d_lll.y.xy * tmp5708 - 2. * d_lll.y.xy * tmp5711 + 4. * d_lll.y.xy * tmp5698 - 4. * d_lll.y.xy * tmp5702 + d_lll.y.xx * tmp5501 + d_lll.y.xx * tmp5495 - 2. * d_lll.y.xx * d_lll.z.zz * tmp5580 + 4. * d_lll.y.xx * d_lll.z.yz * tmp5438 - 2. * d_lll.y.xx * tmp5491 - 4. * d_lll.y.xx * d_lll.z.yz * tmp5570 + 2. * d_lll.y.xx * tmp5485 + d_lll.y.xx * tmp5469 - d_lll.y.xx * tmp5479 + 4. * d_lll.y.xx * d_lll.z.xz * tmp5299 - 2. * d_lll.y.xx * tmp5465 + 2. * d_lll.y.xx * tmp5449 - 2. * d_lll.y.xx * tmp5455 - 2. * d_lll.y.xx * tmp5459 + 2. * d_lll.y.xx * d_lll.z.xx * tmp5246 + 2. * d_lll.y.xx * d_lll.y.zz * tmp5570 - 2. * d_lll.y.xx * d_lll.z.xx * tmp5240 + 3. * d_lll.y.xx * tmp5443 + 2. * d_lll.y.xx * tmp5433 - 6. * d_lll.y.xx * tmp5439 + 2. * d_lll.y.xx * tmp5419 - 2. * d_lll.y.xx * tmp5425 - d_lll.y.xx * tmp5429 + d_lll.y.xx * tmp5409 - d_lll.y.xx * tmp5415 + 4. * d_lll.y.xx * d_lll.y.xz * tmp5286 - 2. * d_lll.y.xx * tmp5405 + 2. * d_lll.y.xx * tmp5389 - 2. * d_lll.y.xx * tmp5395 - 2. * d_lll.y.xx * tmp5399 + d_lll.x.zz * tmp6050 - d_lll.x.zz * tmp6053 + 2. * d_lll.x.zz * d_lll.z.xy * tmp5580 - 2. * d_lll.x.zz * d_lll.z.xy * tmp5500 + d_lll.x.zz * d_lll.z.xx * tmp5350 - d_lll.x.zz * d_lll.z.xx * tmp5355 + 2. * d_lll.x.zz * d_lll.y.zz * tmp5945 - 2. * d_lll.x.zz * d_lll.y.zz * tmp5951 + 4. * d_lll.x.zz * d_lll.y.yz * tmp5779 - 2. * d_lll.x.zz * d_lll.y.yz * tmp5785 - 2. * d_lll.x.zz * d_lll.y.yz * tmp5788 + d_lll.x.zz * d_lll.y.yy * tmp5710 - 2. * d_lll.x.zz * d_lll.y.yy * tmp5775 + d_lll.x.zz * d_lll.y.yy * tmp5707 + 2. * d_lll.x.zz * d_lll.y.xz * tmp5580 - 2. * d_lll.x.zz * d_lll.y.xz * tmp5500 + 4. * d_lll.x.zz * d_lll.y.xy * tmp5438 - 4. * d_lll.x.zz * d_lll.y.xy * tmp5442 + 2. * d_lll.x.zz * d_lll.y.xy * tmp5428 - 2. * d_lll.x.zz * d_lll.y.xy * tmp5432 + d_lll.x.zz * d_lll.y.xx * tmp5296 - d_lll.x.zz * d_lll.y.xx * tmp5305 + 2. * d_lll.x.yz * tmp5952 + 4. * d_lll.x.yz * tmp5942 - 2. * d_lll.x.yz * tmp5946 + 2. * d_lll.x.yz * tmp5932 - 4. * d_lll.x.yz * tmp5936 + 4. * d_lll.x.yz * tmp5922 - 2. * d_lll.x.yz * tmp5926 + 2. * d_lll.x.yz * tmp5912 - 4. * d_lll.x.yz * tmp5916 + 2. * d_lll.x.yz * tmp5908 + 2. * d_lll.x.yz * tmp5898 - 4. * d_lll.x.yz * tmp5902 + 2. * d_lll.x.yz * tmp5882 - 2. * d_lll.x.yz * tmp5886 - 2. * d_lll.x.yz * tmp5892 + 2. * d_lll.x.yz * tmp5872 - 2. * d_lll.x.yz * tmp5876 + 4. * d_lll.x.yz * tmp5862 - 4. * d_lll.x.yz * tmp5866 + 2. * d_lll.x.yz * tmp5852 - 2. * d_lll.x.yz * tmp5856 + 2. * d_lll.x.yz * d_lll.y.xz * tmp5570 + 2. * d_lll.x.yz * d_lll.y.xz * tmp5442 + 2. * d_lll.x.yz * d_lll.y.xz * tmp5432 - 4. * d_lll.x.yz * d_lll.y.xz * tmp5438 + 4. * d_lll.x.yz * d_lll.y.xy * tmp5424 - 2. * d_lll.x.yz * d_lll.y.xz * tmp5428 + 2. * d_lll.x.yz * d_lll.y.xx * tmp5328 - 2. * d_lll.x.yz * d_lll.y.xx * tmp5292 - 4. * d_lll.x.yz * d_lll.y.xy * tmp5474 + 4. * d_lll.x.yz * d_lll.x.zz * tmp5494 - 2. * d_lll.x.yz * d_lll.x.zz * tmp5580 - 2. * d_lll.x.yz * d_lll.x.zz * tmp5500 + d_lll.x.yy * tmp5789 + d_lll.x.yy * tmp5786 + 4. * d_lll.x.yy * d_lll.z.yz * tmp5775 - 2. * d_lll.x.yy * d_lll.z.zz * tmp5779 + 2. * d_lll.x.yy * d_lll.z.yy * tmp5759 - 2. * d_lll.x.yy * tmp5764 - 2. * d_lll.x.yy * tmp5770 + 4. * d_lll.x.yy * d_lll.z.xz * tmp5438 - 4. * d_lll.x.yy * d_lll.z.xz * tmp5570 - 2. * d_lll.x.yy * d_lll.z.yy * tmp5753 + 2. * d_lll.x.yy * tmp5734 - 2. * d_lll.x.yy * tmp5738 + 2. * d_lll.x.yy * tmp5724 - 2. * d_lll.x.yy * tmp5728 + d_lll.x.yy * tmp5716 - d_lll.x.yy * tmp5719 + d_lll.x.yy * tmp5708 - d_lll.x.yy * tmp5711 + 2. * d_lll.x.yy * tmp5698 - 2. * d_lll.x.yy * tmp5702 + d_lll.x.yy * d_lll.y.xx * tmp5278 - d_lll.x.yy * d_lll.y.xx * tmp5283 + 2. * d_lll.x.yy * d_lll.x.zz * tmp5570 + 2. * d_lll.x.yy * d_lll.x.zz * tmp5442 + 3. * d_lll.x.yy * d_lll.x.zz * tmp5432 - 6. * d_lll.x.yy * d_lll.x.zz * tmp5438 + 4. * d_lll.x.yy * d_lll.x.yz * tmp5418 - 2. * d_lll.x.yy * d_lll.x.yz * tmp5474 - 2. * d_lll.x.yy * d_lll.x.yz * tmp5424 - d_lll.x.yy * d_lll.x.zz * tmp5428 + 2. * d_lll.x.xz * tmp5647 + 4. * d_lll.x.xz * tmp5637 - 2. * d_lll.x.xz * tmp5641 + 4. * d_lll.x.xz * tmp5621 - 2. * d_lll.x.xz * tmp5625 - 4. * d_lll.x.xz * tmp5631 + 2. * d_lll.x.xz * tmp5611 - 4. * d_lll.x.xz * tmp5615 + 4. * d_lll.x.xz * tmp5607 + 4. * d_lll.x.xz * tmp5597 - 4. * d_lll.x.xz * tmp5601 + 4. * d_lll.x.xz * d_lll.y.zz * tmp5580 - 2. * d_lll.x.xz * d_lll.y.zz * tmp5500 - 4. * d_lll.x.xz * tmp5591 + 4. * d_lll.x.xz * d_lll.y.yz * tmp5570 - 2. * d_lll.x.xz * d_lll.y.zz * tmp5494 + 2. * d_lll.x.xz * d_lll.y.yy * tmp5424 - 4. * d_lll.x.xz * d_lll.y.yz * tmp5428 + 4. * d_lll.x.xz * d_lll.y.xz * tmp5305 - 2. * d_lll.x.xz * d_lll.y.yy * tmp5418 + 4. * d_lll.x.xz * d_lll.y.xy * tmp5292 - 4. * d_lll.x.xz * d_lll.y.xz * tmp5299 + 2. * d_lll.x.xz * d_lll.x.zz * tmp5350 - 2. * d_lll.x.xz * d_lll.x.zz * tmp5355 - 4. * d_lll.x.xz * d_lll.y.xy * tmp5286 + 4. * d_lll.x.xz * d_lll.x.yz * tmp5299 - 4. * d_lll.x.xz * d_lll.x.yz * tmp5305 + 2. * d_lll.x.xz * d_lll.x.yy * tmp5286 - 2. * d_lll.x.xz * d_lll.x.yy * tmp5292 + 2. * d_lll.x.xy * tmp5501 + 4. * d_lll.x.xy * tmp5491 - 2. * d_lll.x.xy * tmp5495 + 4. * d_lll.x.xy * d_lll.z.yy * tmp5474 - 2. * d_lll.x.xy * tmp5479 - 4. * d_lll.x.xy * tmp5485 + 4. * d_lll.x.xy * tmp5465 - 2. * d_lll.x.xy * tmp5469 + 4. * d_lll.x.xy * tmp5455 - 4. * d_lll.x.xy * tmp5459 + 4. * d_lll.x.xy * tmp5439 - 2. * d_lll.x.xy * tmp5443 - 4. * d_lll.x.xy * tmp5449 + 2. * d_lll.x.xy * tmp5429 - 4. * d_lll.x.xy * tmp5433 + 4. * d_lll.x.xy * tmp5425 + 2. * d_lll.x.xy * tmp5415 - 4. * d_lll.x.xy * tmp5419 + 4. * d_lll.x.xy * tmp5405 - 2. * d_lll.x.xy * tmp5409 + 4. * d_lll.x.xy * tmp5395 - 4. * d_lll.x.xy * tmp5399 + 2. * d_lll.x.xy * d_lll.x.zz * tmp5296 - 2. * d_lll.x.xy * d_lll.x.zz * tmp5305 - 4. * d_lll.x.xy * tmp5389 + 4. * d_lll.x.xy * d_lll.x.yz * tmp5328 - 4. * d_lll.x.xy * d_lll.x.yz * tmp5292 + 2. * d_lll.x.xy * d_lll.x.yy * tmp5278 - 2. * d_lll.x.xy * d_lll.x.yy * tmp5283 + d_lll.x.xx * d_lll.z.zz * tmp5355 + 2. * d_lll.x.xx * d_lll.z.yz * tmp5305 - d_lll.x.xx * d_lll.z.zz * tmp5350 + d_lll.x.xx * d_lll.z.yy * tmp5292 - 2. * d_lll.x.xx * d_lll.z.yz * tmp5296 + d_lll.x.xx * d_lll.z.yy * tmp5286 + 2. * d_lll.x.xx * d_lll.z.xz * tmp5250 - 2. * d_lll.x.xx * d_lll.z.xz * tmp5255 - 2. * d_lll.x.xx * d_lll.z.yy * tmp5328 + 2. * d_lll.x.xx * d_lll.z.xy * tmp5240 - 2. * d_lll.x.xx * d_lll.z.xy * tmp5246 + d_lll.x.xx * d_lll.y.zz * tmp5305 + d_lll.x.xx * d_lll.y.zz * tmp5296 - 2. * d_lll.x.xx * d_lll.y.zz * tmp5299 + 2. * d_lll.x.xx * d_lll.y.yz * tmp5292 + d_lll.x.xx * d_lll.y.yy * tmp5283 - 2. * d_lll.x.xx * d_lll.y.yz * tmp5286 + 2. * d_lll.x.xx * d_lll.y.xz * tmp5240 - 2. * d_lll.x.xx * d_lll.y.xz * tmp5246 - d_lll.x.xx * d_lll.y.yy * tmp5278 + 2. * d_lll.x.xx * d_lll.y.xy * tmp5232 - 2. * d_lll.x.xx * d_lll.y.xy * tmp5237 + d_lll.x.xx * d_lll.x.zz * tmp5255 + 2. * d_lll.x.xx * d_lll.x.yz * tmp5246 - d_lll.x.xx * d_lll.x.zz * tmp5250 + d_lll.x.xx * d_lll.x.yy * tmp5237 - 2. * d_lll.x.xx * d_lll.x.yz * tmp5240 + tmp3431 * tmp230 - d_lll.x.xx * d_lll.x.yy * tmp5232 + tmp1548 * tmp67 + tmp1531 * tmp55 + Z_l.z * tmp1326 - 2. * Z_l.z * tmp1331 - 2. * Z_l.z * tmp1334 - Z_l.z * tmp1337 + 2. * Z_l.z * tmp1494 - 4. * Z_l.z * tmp1320 - 2. * Z_l.z * tmp1323 + Z_l.z * tmp1961 - 2. * Z_l.z * tmp2079 + Z_l.z * tmp1936 - 2. * Z_l.z * tmp1940 - 2. * Z_l.z * tmp1305 - 2. * Z_l.z * tmp2647 - Z_l.z * tmp1311 - 2. * Z_l.z * tmp1314 - Z_l.z * tmp1317 + Z_l.z * tmp1290 - 2. * Z_l.z * tmp5151 - Z_l.z * tmp1300 + Z_l.y * tmp1426 - 2. * Z_l.y * tmp2060 - 2. * Z_l.y * tmp1436 - 2. * Z_l.y * tmp1278 - Z_l.y * tmp1281 - 2. * Z_l.y * tmp1284 - Z_l.y * tmp1287 - Z_l.z * a_l.x * gamma_uu.xz - Z_l.z * a_l.y * gamma_uu.yz - Z_l.z * a_l.z * gamma_uu.zz - Z_l.z * tmp5130 - 2. * Z_l.z * tmp5134 - 2. * Z_l.z * tmp5139 - 2. * Z_l.z * tmp5144 + Z_l.y * tmp1267 - 2. * Z_l.y * tmp1272 + 2. * Z_l.y * tmp1405 - Z_l.y * tmp1261 - 2. * Z_l.y * tmp1264 + Z_l.y * tmp2592 - 2. * Z_l.y * tmp1905 - 2. * Z_l.y * tmp1255 - 4. * Z_l.y * tmp1258 + Z_l.y * tmp1250 - 2. * Z_l.y * tmp1902 + Z_l.x * tmp1231 - 2. * Z_l.x * tmp1235 - Z_l.x * tmp1238 - Z_l.y * a_l.x * gamma_uu.xy - Z_l.y * a_l.y * gamma_uu.yy - Z_l.y * a_l.z * gamma_uu.yz - Z_l.y * tmp5025 - 2. * Z_l.y * tmp5029 - 2. * Z_l.y * tmp5034 - Z_l.y * tmp1240 - 2. * Z_l.y * tmp1896 + Z_l.x * tmp1215 - 2. * Z_l.x * tmp1219 - Z_l.x * tmp2038 - 2. * Z_l.x * tmp1222 - 2. * Z_l.x * tmp1225 - 2. * Z_l.x * tmp1228 + Z_l.x * tmp1200 - 2. * Z_l.x * tmp2505 - Z_l.x * tmp3166 - 2. * Z_l.x * tmp1204 - 2. * Z_l.x * tmp1207 - Z_l.x * tmp1210 - 2. * Z_l.x * tmp1212 + 2. * Z_l.x * tmp1195 - 4. * Z_l.x * tmp1808 + Z_l.x * tmp1191 - 2. * Z_l.x * tmp1805 + K_ll.yy * K_ll.zz * tmp279 - K_ll.yy * K_ll.zz * tmp230 - K_ll.yy * Theta * gamma_uu.yy - 2. * K_ll.yz * Theta * gamma_uu.yz - K_ll.zz * Theta * gamma_uu.zz - 8. * M_PI * rho - Z_l.x * a_l.x * gamma_uu.xx - Z_l.x * a_l.y * gamma_uu.xy - Z_l.x * a_l.z * gamma_uu.xz - Z_l.x * tmp4920 - 2. * Z_l.x * tmp4924 - 2. * Z_l.x * tmp4929 + 2. * K_ll.xz * K_ll.yz * tmp128 - 2. * K_ll.xz * Theta * gamma_uu.xz + 2. * K_ll.xz * K_ll.yy * tmp122 - 2. * K_ll.xz * K_ll.yz * tmp177 + 2. * K_ll.xy * K_ll.zz * tmp177 - 2. * K_ll.xy * K_ll.zz * tmp128 - 2. * K_ll.xy * Theta * gamma_uu.xy - 2. * K_ll.xz * K_ll.yy * tmp116 + 2. * K_ll.xy * K_ll.yz * tmp116 - 2. * K_ll.xy * K_ll.yz * tmp122 + 2. * K_ll.xy * K_ll.xz * tmp61 + K_ll.xx * K_ll.zz * tmp159 - K_ll.xx * K_ll.zz * tmp67 - K_ll.xx * Theta * gamma_uu.xx - 2. * K_ll.xy * K_ll.xz * tmp98 + 2. * K_ll.xx * K_ll.yz * tmp98 - 2. * K_ll.xx * K_ll.yz * tmp61 + K_ll.xx * K_ll.yy * tmp86 - K_ll.xx * K_ll.yy * tmp55);
	deriv->Z_l.x += -alpha * (Theta * a_l.x + 8. * M_PI * S_l.x + K_ll.xz * tmp1337 - K_ll.yy * d_lll.x.xx * tmp55 - 2. * K_ll.yy * d_lll.x.xy * tmp110 - 2. * K_ll.yy * d_lll.x.xz * tmp116 - K_ll.yy * d_lll.x.yy * tmp218 - 2. * K_ll.yy * tmp1340 - K_ll.yy * tmp1345 - 2. * K_ll.yz * d_lll.x.xx * tmp61 - 2. * K_ll.yz * d_lll.x.xy * tmp116 - 2. * K_ll.yz * d_lll.x.xy * tmp122 - 2. * K_ll.yz * d_lll.x.xz * tmp177 - 2. * K_ll.yz * d_lll.x.xz * tmp128 - 2. * K_ll.yz * d_lll.x.yy * tmp224 - 2. * K_ll.yz * d_lll.x.yz * tmp279 - 2. * K_ll.yz * d_lll.x.yz * tmp230 - 2. * K_ll.yz * tmp1365 - K_ll.zz * d_lll.x.xx * tmp67 - 2. * K_ll.zz * d_lll.x.xy * tmp128 - 2. * K_ll.zz * d_lll.x.xz * tmp189 - K_ll.zz * d_lll.x.yy * tmp230 - 2. * K_ll.zz * d_lll.x.yz * tmp291 - K_ll.zz * d_lll.x.zz * tmp332 + 2. * K_ll.xz * tmp1334 + 2. * K_ll.xz * tmp1331 + 2. * K_ll.xz * tmp1323 - K_ll.xz * tmp1326 + 4. * K_ll.xz * tmp1320 + 2. * K_ll.xz * tmp2079 - 2. * K_ll.xz * tmp1494 + K_ll.xz * tmp1317 - K_ll.xz * tmp1961 + 2. * K_ll.xz * tmp1314 + K_ll.xz * tmp1311 + 2. * K_ll.xz * tmp2647 + 2. * K_ll.xz * tmp1305 + 2. * K_ll.xz * tmp1940 + 2. * K_ll.xz * tmp7345 - K_ll.xz * tmp5130 - 2. * K_ll.xz * d_lll.x.xy * tmp61 - 2. * K_ll.xz * d_lll.x.xz * tmp67 - K_ll.xz * tmp1290 - 2. * K_ll.xz * tmp1295 - K_ll.xz * tmp1300 - K_ll.xz * tmp1936 + 2. * K_ll.xz * tmp7343 + 2. * K_ll.xz * tmp7341 + K_ll.xy * tmp1287 + 2. * K_ll.xy * tmp1284 + K_ll.xy * tmp1281 + 2. * K_ll.xy * tmp1278 + 2. * K_ll.xy * tmp1436 + 2. * K_ll.xy * tmp2060 + 2. * K_ll.xy * tmp1272 - K_ll.xy * tmp1426 + 2. * K_ll.xy * tmp1264 - K_ll.xy * tmp1267 + K_ll.xy * tmp1261 + 4. * K_ll.xy * tmp1258 - 2. * K_ll.xy * tmp1405 + 2. * K_ll.xy * tmp1255 + 2. * K_ll.xy * tmp1905 + 2. * K_ll.xy * tmp7262 - K_ll.xy * tmp5025 - 2. * K_ll.xy * d_lll.x.xy * tmp55 - 2. * K_ll.xy * d_lll.x.xz * tmp61 - K_ll.xy * tmp1240 - 2. * K_ll.xy * tmp1245 - K_ll.xy * tmp1250 - K_ll.xy * tmp2592 + 2. * K_ll.xy * tmp7260 + 2. * K_ll.xy * tmp7258 + K_ll.xx * tmp1238 + 2. * K_ll.xx * tmp1235 + 2. * K_ll.xx * tmp1228 - K_ll.xx * tmp1231 + 2. * K_ll.xx * tmp1225 + 2. * K_ll.xx * tmp1222 + K_ll.xx * tmp2038 + 2. * K_ll.xx * tmp1219 + 2. * K_ll.xx * tmp1212 - K_ll.xx * tmp1215 + K_ll.xx * tmp1210 + 2. * K_ll.xx * tmp1207 + 2. * K_ll.xx * tmp1204 + K_ll.xx * tmp3166 + K_ll.xx * tmp2505 + 2. * K_ll.xx * tmp1808 - K_ll.xx * tmp1200 + K_ll.xx * tmp1805 - 2. * K_ll.xx * tmp1195 + 2. * K_ll.xx * tmp7196 - K_ll.xx * tmp1191 + 2. * K_ll.xx * tmp7194 + 2. * K_ll.xx * tmp7192);
	deriv->Z_l.y += alpha * (K_ll.zz * tmp2010 - 8. * M_PI * S_l.y - Theta * a_l.y + 2. * K_ll.zz * tmp2005 + K_ll.zz * tmp2001 + 2. * K_ll.zz * tmp1997 + 2. * K_ll.zz * tmp1992 + K_ll.zz * tmp1990 + K_ll.yz * tmp1326 - 2. * K_ll.yz * tmp1331 - 2. * K_ll.yz * tmp1334 - K_ll.yz * tmp1337 + 2. * K_ll.yz * tmp1494 - 4. * K_ll.yz * tmp1320 - 2. * K_ll.yz * tmp1323 + K_ll.yz * tmp1961 - 2. * K_ll.yz * tmp2079 + K_ll.yz * tmp1317 + 2. * K_ll.yz * tmp1952 + K_ll.yz * tmp1311 + 2. * K_ll.yz * tmp1308 + 2. * K_ll.yz * tmp1943 + K_ll.yz * tmp1936 + K_ll.yz * tmp1290 - 2. * K_ll.yz * tmp5151 - K_ll.yz * tmp1300 + K_ll.yy * tmp1426 - 2. * K_ll.yy * tmp2060 - 2. * K_ll.yy * tmp1436 - 2. * K_ll.yy * tmp1278 - K_ll.yy * tmp1281 - 2. * K_ll.yy * tmp1284 - K_ll.yy * tmp1287 - 2. * K_ll.yz * tmp7341 - 2. * K_ll.yz * tmp7343 - 2. * K_ll.yz * tmp7345 - K_ll.yz * tmp5130 - 2. * K_ll.yz * tmp5134 - 2. * K_ll.yz * tmp5139 - 2. * K_ll.yz * tmp5144 + K_ll.yy * tmp1267 - K_ll.yy * tmp1272 + 2. * K_ll.yy * tmp1405 + K_ll.yy * tmp2592 - K_ll.yy * tmp1905 - 2. * K_ll.yy * tmp1258 + K_ll.yy * tmp1250 - 2. * K_ll.yy * tmp1902 + 2. * K_ll.xz * tmp1875 - 2. * K_ll.yy * tmp7258 - 2. * K_ll.yy * tmp7260 - 2. * K_ll.yy * tmp7262 - K_ll.yy * tmp5025 - 2. * K_ll.yy * tmp5029 - 2. * K_ll.yy * tmp5034 - K_ll.yy * tmp1240 - 2. * K_ll.yy * tmp1896 + 2. * K_ll.xz * tmp3209 + 2. * K_ll.xz * tmp1870 + 2. * K_ll.xz * tmp1865 + 2. * K_ll.xz * tmp1860 + 2. * K_ll.xz * tmp3852 + 2. * K_ll.xz * tmp1855 + 2. * K_ll.xz * tmp3844 + 2. * K_ll.xz * d_lll.y.xx * tmp49 + K_ll.xy * tmp1231 - 2. * K_ll.xy * tmp1235 - K_ll.xy * tmp1238 + K_ll.xy * tmp1215 - K_ll.xy * tmp2038 - 2. * K_ll.xy * tmp1222 - 2. * K_ll.xy * tmp1225 - 2. * K_ll.xy * tmp1228 + 2. * K_ll.xy * tmp1826 + K_ll.xy * tmp1210 + 2. * K_ll.xy * tmp1816 + 2. * K_ll.xy * tmp1811 + K_ll.xy * tmp3166 + K_ll.xy * tmp1200 - 2. * K_ll.xy * tmp2505 + 2. * K_ll.xy * tmp1195 - 4. * K_ll.xy * tmp1808 + K_ll.xy * tmp1191 - 2. * K_ll.xy * tmp1805 + K_ll.xx * tmp1797 - 2. * K_ll.xy * tmp7192 - 2. * K_ll.xy * tmp7194 - 2. * K_ll.xy * tmp7196 - K_ll.xy * tmp4920 - 2. * K_ll.xy * tmp4924 - 2. * K_ll.xy * tmp4929 + 2. * K_ll.xx * tmp1792 + K_ll.xx * tmp1788 + 2. * K_ll.xx * tmp1784 + 2. * K_ll.xx * tmp1779 + K_ll.xx * d_lll.y.xx * tmp37);
	deriv->Z_l.z += alpha * (K_ll.zz * tmp1326 - K_ll.zz * tmp1331 - 8. * M_PI * S_l.z - Theta * a_l.z + 2. * K_ll.zz * tmp1494 - 2. * K_ll.zz * tmp1320 + K_ll.zz * tmp1961 - K_ll.zz * tmp2079 + K_ll.zz * tmp1936 - 2. * K_ll.zz * tmp1940 - 2. * K_ll.zz * tmp1305 - 2. * K_ll.zz * tmp2647 - K_ll.zz * tmp1311 - 2. * K_ll.zz * tmp1314 - K_ll.zz * tmp1317 + K_ll.zz * tmp1290 - 2. * K_ll.zz * tmp5151 - K_ll.zz * tmp1300 + K_ll.yz * tmp1287 - 2. * K_ll.zz * tmp7341 - 2. * K_ll.zz * tmp7343 - 2. * K_ll.zz * tmp7345 - K_ll.zz * tmp5130 - 2. * K_ll.zz * tmp5134 - 2. * K_ll.zz * tmp5139 - 2. * K_ll.zz * tmp5144 + 2. * K_ll.yz * tmp2633 + K_ll.yz * tmp1281 + 2. * K_ll.yz * tmp2624 + 2. * K_ll.yz * tmp1275 + K_ll.yz * tmp1426 + K_ll.yz * tmp1267 - 2. * K_ll.yz * tmp1272 + 2. * K_ll.yz * tmp1405 - K_ll.yz * tmp1261 - 2. * K_ll.yz * tmp1264 + K_ll.yz * tmp2592 - 2. * K_ll.yz * tmp1905 - 2. * K_ll.yz * tmp1255 - 4. * K_ll.yz * tmp1258 + K_ll.yz * tmp1250 - 2. * K_ll.yz * tmp1902 + K_ll.yy * tmp2583 - 2. * K_ll.yz * tmp7258 - 2. * K_ll.yz * tmp7260 - 2. * K_ll.yz * tmp7262 - K_ll.yz * tmp5025 - 2. * K_ll.yz * tmp5029 - 2. * K_ll.yz * tmp5034 - K_ll.yz * tmp1240 - 2. * K_ll.yz * tmp1896 + 2. * K_ll.yy * tmp2578 + K_ll.yy * tmp2574 + 2. * K_ll.yy * tmp2569 + 2. * K_ll.yy * tmp2565 + K_ll.yy * tmp2563 + K_ll.xz * tmp1238 + 2. * K_ll.xz * tmp2536 + K_ll.xz * tmp1231 + 2. * K_ll.xz * tmp2526 + 2. * K_ll.xz * tmp1836 + K_ll.xz * tmp2038 + K_ll.xz * tmp1215 - 2. * K_ll.xz * tmp1219 + K_ll.xz * tmp1200 - 2. * K_ll.xz * tmp2505 - K_ll.xz * tmp3166 - 2. * K_ll.xz * tmp1204 - 2. * K_ll.xz * tmp1207 - K_ll.xz * tmp1210 - 2. * K_ll.xz * tmp1212 + 2. * K_ll.xz * tmp1195 - 4. * K_ll.xz * tmp1808 + K_ll.xz * tmp1191 - 2. * K_ll.xz * tmp1805 + 2. * K_ll.xy * tmp2497 - 2. * K_ll.xz * tmp7192 - 2. * K_ll.xz * tmp7194 - 2. * K_ll.xz * tmp7196 - K_ll.xz * tmp4920 - 2. * K_ll.xz * tmp4924 - 2. * K_ll.xz * tmp4929 + 2. * K_ll.xy * tmp2492 + 2. * K_ll.xy * tmp2700 + 2. * K_ll.xy * tmp2487 + 2. * K_ll.xy * tmp2482 + 2. * K_ll.xy * tmp2692 + 2. * K_ll.xy * tmp2477 + 2. * K_ll.xy * tmp2685 + 2. * K_ll.xy * tmp2683 + K_ll.xx * tmp2454 + 2. * K_ll.xx * tmp2449 + K_ll.xx * tmp2445 + 2. * K_ll.xx * tmp2440 + 2. * K_ll.xx * tmp2436 + K_ll.xx * tmp3890);
	// END CUT
<? end ?>
<? if false then ?>
	// BEGIN CUT from symmath/tests/output/Z4.html
	{
		real const tmp1 = S * M_PI;
		real const tmp2 = rho * M_PI;
		(deriv)->K_ll.xx += alpha * (-d_ull.x.xx * d_ull.x.xx + -d_ull.y.xy * d_ull.y.xy + -d_ull.z.xz * d_ull.z.xz + -d_llu.x.y.y * d_llu.x.y.y + -d_llu.x.z.z * d_llu.x.z.z + -d_llu.y.x.y * d_llu.y.x.y + -d_llu.z.x.z * d_llu.z.x.z + -2. * d_llu.x.x.x * d_llu.x.x.x + a_u.x * d_lll.x.xx + a_u.y * d_lll.x.xy + a_u.z * d_lll.x.xz + a_l.x * d_l.x + -a_l.x * e_l.x + d_u.x * d_lll.x.xx + -d_u.y * d_lll.y.xx + -d_u.z * d_lll.z.xx + K_ll.xx * tr_K + -2. * Z_u.x * d_lll.x.xx + 2. * Z_u.y * d_lll.y.xx + 2. * Z_u.z * d_lll.z.xx + -2. * Z_l.x * a_l.x + 2. * d_u.y * d_lll.x.xy + 2. * d_u.z * d_lll.x.xz + -2. * e_u.x * d_lll.x.xx + -2. * e_u.y * d_lll.x.xy + -2. * e_u.z * d_lll.x.xz + -2. * K_ul.x.x * K_ll.xx + -2. * K_ul.y.x * K_ll.xy + -2. * K_ul.z.x * K_ll.xz + -2. * K_ll.xx * Theta + 2. * d_ull.x.xx * d_llu.x.x.x + -2. * d_ull.x.xy * d_ull.y.xx + 2. * d_ull.x.xy * d_llu.x.x.y + -2. * d_ull.x.xz * d_ull.z.xx + 2. * d_ull.x.xz * d_llu.x.x.z + 2. * d_ull.y.xx * d_llu.y.x.x + 2. * d_ull.y.xy * d_llu.y.x.y + -2. * d_ull.y.xz * d_ull.z.xy + 2. * d_ull.y.xz * d_llu.y.x.z + 2. * d_ull.z.xx * d_llu.z.x.x + 2. * d_ull.z.xy * d_llu.z.x.y + 2. * d_ull.z.xz * d_llu.z.x.z + 2. * d_luu.x.xx * d_lll.x.xx + 2. * d_luu.x.xy * d_lll.x.xy + 2. * d_luu.x.xy * d_lll.y.xx + 2. * d_luu.x.xz * d_lll.x.xz + 2. * d_luu.x.xz * d_lll.z.xx + 2. * d_luu.x.yy * d_lll.y.xy + 2. * d_luu.x.yz * d_lll.y.xz + 2. * d_luu.x.yz * d_lll.z.xy + 2. * d_luu.x.zz * d_lll.z.xz + -2. * d_llu.x.x.y * d_llu.x.y.x + -2. * d_llu.x.x.y * d_llu.y.x.x + -2. * d_llu.x.x.z * d_llu.x.z.x + -2. * d_llu.x.x.z * d_llu.z.x.x + -2. * d_llu.x.y.z * d_llu.x.z.y + -2. * d_llu.y.x.z * d_llu.z.x.y + -4. * Z_u.y * d_lll.x.xy + -4. * Z_u.z * d_lll.x.xz + -8. * S_ll.xx * M_PI + -4. * gamma_ll.xx * tmp2 + 4. * gamma_ll.xx * tmp1);
		(deriv)->K_ll.xy += (alpha * (a_u.x * d_lll.x.xy + a_u.x * d_lll.y.xx + a_u.y * d_lll.x.yy + a_u.y * d_lll.y.xy + a_u.z * d_lll.x.yz + a_u.z * d_lll.y.xz + a_l.x * d_l.y + -a_l.x * e_l.y + a_l.y * d_l.x + -a_l.y * e_l.x + -2. * Z_l.x * a_l.y + -2. * Z_l.y * a_l.x + 2. * d_u.x * d_lll.y.xx + 2. * d_u.y * d_lll.x.yy + 2. * d_u.z * d_lll.x.yz + 2. * d_u.z * d_lll.y.xz + -2. * d_u.z * d_lll.z.xy + -2. * e_u.x * d_lll.x.xy + -2. * e_u.x * d_lll.y.xx + -2. * e_u.y * d_lll.x.yy + -2. * e_u.y * d_lll.y.xy + -2. * e_u.z * d_lll.x.yz + -2. * e_u.z * d_lll.y.xz + 2. * K_ll.xy * tr_K + -2. * d_ull.x.xx * d_ull.x.xy + 2. * d_ull.x.xx * d_llu.x.y.x + -2. * d_ull.x.xy * d_ull.y.xy + 2. * d_ull.x.xy * d_llu.x.x.x + 2. * d_ull.x.xy * d_llu.x.y.y + -2. * d_ull.x.xz * d_ull.z.xy + 2. * d_ull.x.xz * d_llu.x.y.z + -2. * d_ull.x.yy * d_ull.y.xx + 2. * d_ull.x.yy * d_llu.x.x.y + -2. * d_ull.x.yz * d_ull.z.xx + 2. * d_ull.x.yz * d_llu.x.x.z + 2. * d_ull.y.xx * d_llu.y.y.x + -2. * d_ull.y.xy * d_ull.y.yy + 2. * d_ull.y.xy * d_llu.y.x.x + 2. * d_ull.y.xy * d_llu.y.y.y + -2. * d_ull.y.xz * d_ull.z.yy + 2. * d_ull.y.xz * d_llu.y.y.z + 2. * d_ull.y.yy * d_llu.y.x.y + -2. * d_ull.y.yz * d_ull.z.xy + 2. * d_ull.y.yz * d_llu.y.x.z + 2. * d_ull.z.xx * d_llu.z.y.x + 2. * d_ull.z.xy * d_llu.z.x.x + 2. * d_ull.z.xy * d_llu.z.y.y + -2. * d_ull.z.xz * d_ull.z.yz + 2. * d_ull.z.xz * d_llu.z.y.z + 2. * d_ull.z.yy * d_llu.z.x.y + 2. * d_ull.z.yz * d_llu.z.x.z + 2. * d_luu.x.xx * d_lll.x.xy + 2. * d_luu.x.xy * d_lll.x.yy + 2. * d_luu.x.xy * d_lll.y.xy + 2. * d_luu.x.xz * d_lll.x.yz + 2. * d_luu.x.xz * d_lll.z.xy + 2. * d_luu.x.yy * d_lll.y.yy + 2. * d_luu.x.yz * d_lll.y.yz + 2. * d_luu.x.yz * d_lll.z.yy + 2. * d_luu.x.zz * d_lll.z.yz + -2. * d_llu.x.x.x * d_llu.x.y.x + -2. * d_llu.x.x.x * d_llu.y.x.x + -2. * d_llu.x.x.z * d_llu.y.z.x + -2. * d_llu.x.x.z * d_llu.z.y.x + 2. * d_lll.x.xx * d_luu.y.xx + 2. * d_lll.x.xy * d_luu.y.xy + 2. * d_lll.x.xz * d_luu.y.xz + -2. * d_llu.x.y.x * d_llu.y.x.y + -2. * d_llu.x.y.y * d_llu.y.x.x + -2. * d_llu.x.y.y * d_llu.y.y.y + -2. * d_llu.x.y.z * d_llu.y.z.y + -2. * d_llu.x.y.z * d_llu.z.x.x + -2. * d_llu.x.z.x * d_llu.y.x.z + -2. * d_llu.x.z.y * d_llu.y.y.z + -2. * d_llu.x.z.z * d_llu.y.z.z + 2. * d_luu.y.xy * d_lll.y.xx + 2. * d_luu.y.xz * d_lll.z.xx + 2. * d_luu.y.yy * d_lll.y.xy + 2. * d_luu.y.yz * d_lll.y.xz + 2. * d_luu.y.yz * d_lll.z.xy + 2. * d_luu.y.zz * d_lll.z.xz + -2. * d_llu.y.x.y * d_llu.y.y.y + -2. * d_llu.y.x.z * d_llu.z.y.y + -2. * d_llu.y.y.z * d_llu.z.x.y + -2. * d_llu.z.x.z * d_llu.z.y.z + -4. * Z_u.x * d_lll.y.xx + -4. * Z_u.y * d_lll.x.yy + -4. * Z_u.z * d_lll.x.yz + -4. * Z_u.z * d_lll.y.xz + 4. * Z_u.z * d_lll.z.xy + -4. * K_ul.x.x * K_ll.xy + -4. * K_ul.y.x * K_ll.yy + -4. * K_ul.z.x * K_ll.yz + -4. * K_ll.xy * Theta + -4. * d_llu.x.x.y * d_llu.y.y.x + -16. * S_ll.xy * M_PI + -8. * gamma_ll.xy * tmp2 + 8. * gamma_ll.xy * tmp1)) / 2.;
		(deriv)->K_ll.xz += (alpha * (a_u.x * d_lll.x.xz + a_u.x * d_lll.z.xx + a_u.y * d_lll.x.yz + a_u.y * d_lll.z.xy + a_u.z * d_lll.x.zz + a_u.z * d_lll.z.xz + a_l.x * d_l.z + -a_l.x * e_l.z + a_l.z * d_l.x + -a_l.z * e_l.x + -2. * Z_l.x * a_l.z + -2. * Z_l.z * a_l.x + 2. * d_u.x * d_lll.z.xx + 2. * d_u.y * d_lll.x.yz + -2. * d_u.y * d_lll.y.xz + 2. * d_u.y * d_lll.z.xy + 2. * d_u.z * d_lll.x.zz + -2. * e_u.x * d_lll.x.xz + -2. * e_u.x * d_lll.z.xx + -2. * e_u.y * d_lll.x.yz + -2. * e_u.y * d_lll.z.xy + -2. * e_u.z * d_lll.x.zz + -2. * e_u.z * d_lll.z.xz + 2. * K_ll.xz * tr_K + -2. * d_ull.x.xx * d_ull.x.xz + 2. * d_ull.x.xx * d_llu.x.z.x + -2. * d_ull.x.xy * d_ull.y.xz + 2. * d_ull.x.xy * d_llu.x.z.y + -2. * d_ull.x.xz * d_ull.z.xz + 2. * d_ull.x.xz * d_llu.x.x.x + 2. * d_ull.x.xz * d_llu.x.z.z + -2. * d_ull.x.yz * d_ull.y.xx + 2. * d_ull.x.yz * d_llu.x.x.y + -2. * d_ull.x.zz * d_ull.z.xx + 2. * d_ull.x.zz * d_llu.x.x.z + 2. * d_ull.y.xx * d_llu.y.z.x + -2. * d_ull.y.xy * d_ull.y.yz + 2. * d_ull.y.xy * d_llu.y.z.y + -2. * d_ull.y.xz * d_ull.z.yz + 2. * d_ull.y.xz * d_llu.y.x.x + 2. * d_ull.y.xz * d_llu.y.z.z + 2. * d_ull.y.yz * d_llu.y.x.y + -2. * d_ull.y.zz * d_ull.z.xy + 2. * d_ull.y.zz * d_llu.y.x.z + 2. * d_ull.z.xx * d_llu.z.z.x + 2. * d_ull.z.xy * d_llu.z.z.y + -2. * d_ull.z.xz * d_ull.z.zz + 2. * d_ull.z.xz * d_llu.z.x.x + 2. * d_ull.z.xz * d_llu.z.z.z + 2. * d_ull.z.yz * d_llu.z.x.y + 2. * d_ull.z.zz * d_llu.z.x.z + 2. * d_luu.x.xx * d_lll.x.xz + 2. * d_luu.x.xy * d_lll.x.yz + 2. * d_luu.x.xy * d_lll.y.xz + 2. * d_luu.x.xz * d_lll.x.zz + 2. * d_luu.x.xz * d_lll.z.xz + 2. * d_luu.x.yy * d_lll.y.yz + 2. * d_luu.x.yz * d_lll.y.zz + 2. * d_luu.x.yz * d_lll.z.yz + 2. * d_luu.x.zz * d_lll.z.zz + -2. * d_llu.x.x.x * d_llu.x.z.x + -2. * d_llu.x.x.x * d_llu.z.x.x + -2. * d_llu.x.x.y * d_llu.y.z.x + -2. * d_llu.x.x.y * d_llu.z.y.x + 2. * d_lll.x.xx * d_luu.z.xx + 2. * d_lll.x.xy * d_luu.z.xy + 2. * d_lll.x.xz * d_luu.z.xz + -2. * d_llu.x.y.x * d_llu.z.x.y + -2. * d_llu.x.y.y * d_llu.z.y.y + -2. * d_llu.x.y.z * d_llu.z.z.y + -2. * d_llu.x.z.x * d_llu.z.x.z + -2. * d_llu.x.z.y * d_llu.y.x.x + -2. * d_llu.x.z.y * d_llu.z.y.z + -2. * d_llu.x.z.z * d_llu.z.x.x + -2. * d_llu.x.z.z * d_llu.z.z.z + -2. * d_llu.y.x.y * d_llu.y.z.y + -2. * d_llu.y.x.z * d_llu.z.z.y + 2. * d_lll.y.xx * d_luu.z.xy + 2. * d_lll.y.xy * d_luu.z.yy + 2. * d_lll.y.xz * d_luu.z.yz + -2. * d_llu.y.z.z * d_llu.z.x.y + 2. * d_luu.z.xz * d_lll.z.xx + 2. * d_luu.z.yz * d_lll.z.xy + 2. * d_luu.z.zz * d_lll.z.xz + -2. * d_llu.z.x.z * d_llu.z.z.z + -4. * Z_u.x * d_lll.z.xx + -4. * Z_u.y * d_lll.x.yz + 4. * Z_u.y * d_lll.y.xz + -4. * Z_u.y * d_lll.z.xy + -4. * Z_u.z * d_lll.x.zz + -4. * K_ul.x.x * K_ll.xz + -4. * K_ul.y.x * K_ll.yz + -4. * K_ul.z.x * K_ll.zz + -4. * K_ll.xz * Theta + -4. * d_llu.x.x.z * d_llu.z.z.x + -16. * S_ll.xz * M_PI + -8. * gamma_ll.xz * tmp2 + 8. * gamma_ll.xz * tmp1)) / 2.;
		(deriv)->K_ll.yy += alpha * (-d_ull.x.xy * d_ull.x.xy + -d_ull.y.yy * d_ull.y.yy + -d_ull.z.yz * d_ull.z.yz + -d_llu.x.y.x * d_llu.x.y.x + -d_llu.y.x.x * d_llu.y.x.x + -d_llu.y.z.z * d_llu.y.z.z + -d_llu.z.y.z * d_llu.z.y.z + -2. * d_llu.y.y.y * d_llu.y.y.y + a_u.x * d_lll.y.xy + a_u.y * d_lll.y.yy + a_u.z * d_lll.y.yz + a_l.y * d_l.y + -a_l.y * e_l.y + -d_u.x * d_lll.x.yy + d_u.y * d_lll.y.yy + -d_u.z * d_lll.z.yy + K_ll.yy * tr_K + 2. * Z_u.x * d_lll.x.yy + -2. * Z_u.y * d_lll.y.yy + 2. * Z_u.z * d_lll.z.yy + -2. * Z_l.y * a_l.y + 2. * d_u.x * d_lll.y.xy + 2. * d_u.z * d_lll.y.yz + -2. * e_u.x * d_lll.y.xy + -2. * e_u.y * d_lll.y.yy + -2. * e_u.z * d_lll.y.yz + -2. * K_ll.xy * K_ul.x.y + -2. * K_ul.y.y * K_ll.yy + -2. * K_ul.z.y * K_ll.yz + -2. * K_ll.yy * Theta + 2. * d_ull.x.xy * d_llu.x.y.x + -2. * d_ull.x.yy * d_ull.y.xy + 2. * d_ull.x.yy * d_llu.x.y.y + -2. * d_ull.x.yz * d_ull.z.xy + 2. * d_ull.x.yz * d_llu.x.y.z + 2. * d_ull.y.xy * d_llu.y.y.x + 2. * d_ull.y.yy * d_llu.y.y.y + -2. * d_ull.y.yz * d_ull.z.yy + 2. * d_ull.y.yz * d_llu.y.y.z + 2. * d_ull.z.xy * d_llu.z.y.x + 2. * d_ull.z.yy * d_llu.z.y.y + 2. * d_ull.z.yz * d_llu.z.y.z + 2. * d_lll.x.xy * d_luu.y.xx + -2. * d_llu.x.y.y * d_llu.y.y.x + -2. * d_llu.x.y.z * d_llu.z.y.x + 2. * d_lll.x.yy * d_luu.y.xy + 2. * d_lll.x.yz * d_luu.y.xz + 2. * d_luu.y.xy * d_lll.y.xy + 2. * d_luu.y.xz * d_lll.z.xy + 2. * d_luu.y.yy * d_lll.y.yy + 2. * d_luu.y.yz * d_lll.y.yz + 2. * d_luu.y.yz * d_lll.z.yy + 2. * d_luu.y.zz * d_lll.z.yz + -2. * d_llu.y.x.y * d_llu.y.y.x + -2. * d_llu.y.x.z * d_llu.y.z.x + -2. * d_llu.y.y.z * d_llu.y.z.y + -2. * d_llu.y.y.z * d_llu.z.y.y + -4. * Z_u.x * d_lll.y.xy + -4. * Z_u.z * d_lll.y.yz + -8. * S_ll.yy * M_PI + -4. * gamma_ll.yy * tmp2 + 4. * gamma_ll.yy * tmp1);
		(deriv)->K_ll.yz += (alpha * (a_u.x * d_lll.y.xz + a_u.x * d_lll.z.xy + a_u.y * d_lll.y.yz + a_u.y * d_lll.z.yy + a_u.z * d_lll.y.zz + a_u.z * d_lll.z.yz + a_l.y * d_l.z + -a_l.y * e_l.z + a_l.z * d_l.y + -a_l.z * e_l.y + -2. * Z_l.y * a_l.z + -2. * Z_l.z * a_l.y + -2. * d_u.x * d_lll.x.yz + 2. * d_u.x * d_lll.y.xz + 2. * d_u.x * d_lll.z.xy + 2. * d_u.y * d_lll.z.yy + 2. * d_u.z * d_lll.y.zz + -2. * e_u.x * d_lll.y.xz + -2. * e_u.x * d_lll.z.xy + -2. * e_u.y * d_lll.y.yz + -2. * e_u.y * d_lll.z.yy + -2. * e_u.z * d_lll.y.zz + -2. * e_u.z * d_lll.z.yz + 2. * K_ll.yz * tr_K + -2. * d_ull.x.xy * d_ull.x.xz + 2. * d_ull.x.xy * d_llu.x.z.x + 2. * d_ull.x.xz * d_llu.x.y.x + -2. * d_ull.x.yy * d_ull.y.xz + 2. * d_ull.x.yy * d_llu.x.z.y + -2. * d_ull.x.yz * d_ull.y.xy + -2. * d_ull.x.yz * d_ull.z.xz + 2. * d_ull.x.yz * d_llu.x.y.y + 2. * d_ull.x.yz * d_llu.x.z.z + -2. * d_ull.x.zz * d_ull.z.xy + 2. * d_ull.x.zz * d_llu.x.y.z + 2. * d_ull.y.xy * d_llu.y.z.x + 2. * d_ull.y.xz * d_llu.y.y.x + -2. * d_ull.y.yy * d_ull.y.yz + 2. * d_ull.y.yy * d_llu.y.z.y + -2. * d_ull.y.yz * d_ull.z.yz + 2. * d_ull.y.yz * d_llu.y.y.y + 2. * d_ull.y.yz * d_llu.y.z.z + -2. * d_ull.y.zz * d_ull.z.yy + 2. * d_ull.y.zz * d_llu.y.y.z + 2. * d_ull.z.xy * d_llu.z.z.x + 2. * d_ull.z.xz * d_llu.z.y.x + 2. * d_ull.z.yy * d_llu.z.z.y + -2. * d_ull.z.yz * d_ull.z.zz + 2. * d_ull.z.yz * d_llu.z.y.y + 2. * d_ull.z.yz * d_llu.z.z.z + 2. * d_ull.z.zz * d_llu.z.y.z + 2. * d_lll.x.xy * d_luu.z.xx + 2. * d_lll.x.xz * d_luu.y.xx + -2. * d_llu.x.y.x * d_llu.x.z.x + -2. * d_llu.x.y.y * d_llu.y.z.x + -2. * d_llu.x.y.z * d_llu.z.z.x + 2. * d_lll.x.yy * d_luu.z.xy + 2. * d_lll.x.yz * d_luu.y.xy + 2. * d_lll.x.yz * d_luu.z.xz + -2. * d_llu.x.z.y * d_llu.y.y.x + -2. * d_llu.x.z.z * d_llu.z.y.x + 2. * d_lll.x.zz * d_luu.y.xz + 2. * d_luu.y.xy * d_lll.y.xz + 2. * d_luu.y.xz * d_lll.z.xz + 2. * d_luu.y.yy * d_lll.y.yz + 2. * d_luu.y.yz * d_lll.y.zz + 2. * d_luu.y.yz * d_lll.z.yz + 2. * d_luu.y.zz * d_lll.z.zz + -2. * d_llu.y.x.x * d_llu.z.x.x + -2. * d_llu.y.x.y * d_llu.z.y.x + -2. * d_llu.y.x.z * d_llu.z.z.x + 2. * d_lll.y.xy * d_luu.z.xy + -2. * d_llu.y.y.x * d_llu.z.x.y + -2. * d_llu.y.y.y * d_llu.y.z.y + -2. * d_llu.y.y.y * d_llu.z.y.y + 2. * d_lll.y.yy * d_luu.z.yy + 2. * d_lll.y.yz * d_luu.z.yz + -2. * d_llu.y.z.x * d_llu.z.x.z + -2. * d_llu.y.z.y * d_llu.z.y.z + -2. * d_llu.y.z.z * d_llu.z.y.y + -2. * d_llu.y.z.z * d_llu.z.z.z + 2. * d_luu.z.xz * d_lll.z.xy + 2. * d_luu.z.yz * d_lll.z.yy + 2. * d_luu.z.zz * d_lll.z.yz + -2. * d_llu.z.y.z * d_llu.z.z.z + 4. * Z_u.x * d_lll.x.yz + -4. * Z_u.x * d_lll.y.xz + -4. * Z_u.x * d_lll.z.xy + -4. * Z_u.y * d_lll.z.yy + -4. * Z_u.z * d_lll.y.zz + -4. * K_ll.xz * K_ul.x.y + -4. * K_ul.y.y * K_ll.yz + -4. * K_ul.z.y * K_ll.zz + -4. * K_ll.yz * Theta + -4. * d_llu.y.y.z * d_llu.z.z.y + -16. * S_ll.yz * M_PI + -8. * gamma_ll.yz * tmp2 + 8. * gamma_ll.yz * tmp1)) / 2.;
		(deriv)->K_ll.zz += alpha * (-d_ull.x.xz * d_ull.x.xz + -d_ull.y.yz * d_ull.y.yz + -d_ull.z.zz * d_ull.z.zz + -d_llu.x.z.x * d_llu.x.z.x + -d_llu.y.z.y * d_llu.y.z.y + -d_llu.z.x.x * d_llu.z.x.x + -d_llu.z.y.y * d_llu.z.y.y + -2. * d_llu.z.z.z * d_llu.z.z.z + a_u.x * d_lll.z.xz + a_u.y * d_lll.z.yz + a_u.z * d_lll.z.zz + a_l.z * d_l.z + -a_l.z * e_l.z + -d_u.x * d_lll.x.zz + -d_u.y * d_lll.y.zz + d_u.z * d_lll.z.zz + K_ll.zz * tr_K + 2. * Z_u.x * d_lll.x.zz + 2. * Z_u.y * d_lll.y.zz + -2. * Z_u.z * d_lll.z.zz + -2. * Z_l.z * a_l.z + 2. * d_u.x * d_lll.z.xz + 2. * d_u.y * d_lll.z.yz + -2. * e_u.x * d_lll.z.xz + -2. * e_u.y * d_lll.z.yz + -2. * e_u.z * d_lll.z.zz + -2. * K_ll.xz * K_ul.x.z + -2. * K_ll.yz * K_ul.y.z + -2. * K_ul.z.z * K_ll.zz + -2. * K_ll.zz * Theta + 2. * d_ull.x.xz * d_llu.x.z.x + -2. * d_ull.x.yz * d_ull.y.xz + 2. * d_ull.x.yz * d_llu.x.z.y + -2. * d_ull.x.zz * d_ull.z.xz + 2. * d_ull.x.zz * d_llu.x.z.z + 2. * d_ull.y.xz * d_llu.y.z.x + 2. * d_ull.y.yz * d_llu.y.z.y + -2. * d_ull.y.zz * d_ull.z.yz + 2. * d_ull.y.zz * d_llu.y.z.z + 2. * d_ull.z.xz * d_llu.z.z.x + 2. * d_ull.z.yz * d_llu.z.z.y + 2. * d_ull.z.zz * d_llu.z.z.z + 2. * d_lll.x.xz * d_luu.z.xx + 2. * d_lll.x.yz * d_luu.z.xy + -2. * d_llu.x.z.y * d_llu.y.z.x + -2. * d_llu.x.z.z * d_llu.z.z.x + 2. * d_lll.x.zz * d_luu.z.xz + 2. * d_lll.y.xz * d_luu.z.xy + 2. * d_lll.y.yz * d_luu.z.yy + -2. * d_llu.y.z.z * d_llu.z.z.y + 2. * d_lll.y.zz * d_luu.z.yz + 2. * d_luu.z.xz * d_lll.z.xz + 2. * d_luu.z.yz * d_lll.z.yz + 2. * d_luu.z.zz * d_lll.z.zz + -2. * d_llu.z.x.y * d_llu.z.y.x + -2. * d_llu.z.x.z * d_llu.z.z.x + -2. * d_llu.z.y.z * d_llu.z.z.y + -4. * Z_u.x * d_lll.z.xz + -4. * Z_u.y * d_lll.z.yz + -8. * S_ll.zz * M_PI + -4. * gamma_ll.zz * tmp2 + 4. * gamma_ll.zz * tmp1);
		(deriv)->Theta += (alpha * (2. * tr_K * tr_K + a_l.x * d_u.x + -a_l.x * e_u.x + a_l.y * d_u.y + -a_l.y * e_u.y + a_l.z * d_u.z + -a_l.z * e_u.z + -2. * d_u.x * d_l.x + -2. * d_u.y * d_l.y + -2. * d_u.z * d_l.z + -2. * K_uu.xx * K_ll.xx + -2. * K_uu.yy * K_ll.yy + -2. * K_uu.zz * K_ll.zz + 2. * d_uuu.x.xx * d_lll.x.xx + -2. * d_uuu.x.xy * d_lll.x.xy + -2. * d_uuu.x.xz * d_lll.x.xz + -2. * d_uul.x.xx * d_llu.x.x.x + -2. * d_uul.x.xy * d_llu.x.x.y + -2. * d_uul.x.xz * d_llu.x.x.z + -2. * d_uul.x.yx * d_llu.x.y.x + -2. * d_uul.x.yy * d_llu.x.y.y + -2. * d_uul.x.yz * d_llu.x.y.z + -2. * d_uul.x.zx * d_llu.x.z.x + -2. * d_uul.x.zy * d_llu.x.z.y + -2. * d_uul.x.zz * d_llu.x.z.z + -2. * d_ulu.x.xx * d_ull.x.xx + 2. * d_ulu.x.xx * d_llu.x.x.x + -2. * d_ulu.x.xy * d_ull.x.xy + 2. * d_ulu.x.xy * d_llu.x.y.x + -2. * d_ulu.x.xz * d_ull.x.xz + 2. * d_ulu.x.xz * d_llu.x.z.x + 2. * d_ull.x.xx * d_luu.x.xx + -2. * d_ull.x.xy * d_ulu.y.xx + -2. * d_ull.x.xz * d_ulu.z.xx + -2. * d_ulu.x.yx * d_ull.y.xx + 2. * d_ulu.x.yx * d_llu.x.x.y + -2. * d_ulu.x.yy * d_ull.y.xy + 2. * d_ulu.x.yy * d_llu.x.y.y + -2. * d_ulu.x.yz * d_ull.y.xz + 2. * d_ulu.x.yz * d_llu.x.z.y + -2. * d_ull.x.yy * d_ulu.y.xy + 2. * d_ull.x.yy * d_luu.x.yy + -2. * d_ull.x.yz * d_ulu.y.xz + -2. * d_ull.x.yz * d_ulu.z.xy + -2. * d_ulu.x.zx * d_ull.z.xx + 2. * d_ulu.x.zx * d_llu.x.x.z + -2. * d_ulu.x.zy * d_ull.z.xy + 2. * d_ulu.x.zy * d_llu.x.y.z + -2. * d_ulu.x.zz * d_ull.z.xz + 2. * d_ulu.x.zz * d_llu.x.z.z + -2. * d_ull.x.zz * d_ulu.z.xz + 2. * d_ull.x.zz * d_luu.x.zz + -2. * d_uuu.y.xy * d_lll.y.xy + -2. * d_uul.y.xx * d_llu.y.x.x + -2. * d_uul.y.xy * d_llu.y.x.y + -2. * d_uul.y.xz * d_llu.y.x.z + 2. * d_uuu.y.yy * d_lll.y.yy + -2. * d_uuu.y.yz * d_lll.y.yz + -2. * d_uul.y.yx * d_llu.y.y.x + -2. * d_uul.y.yy * d_llu.y.y.y + -2. * d_uul.y.yz * d_llu.y.y.z + -2. * d_uul.y.zx * d_llu.y.z.x + -2. * d_uul.y.zy * d_llu.y.z.y + -2. * d_uul.y.zz * d_llu.y.z.z + 2. * d_ulu.y.xx * d_llu.y.x.x + 2. * d_ulu.y.xy * d_llu.y.y.x + 2. * d_ulu.y.xz * d_llu.y.z.x + 2. * d_ull.y.xx * d_luu.y.xx + -2. * d_ull.y.xy * d_ulu.y.yx + -2. * d_ull.y.xz * d_ulu.z.yx + 2. * d_ulu.y.yx * d_llu.y.x.y + -2. * d_ulu.y.yy * d_ull.y.yy + 2. * d_ulu.y.yy * d_llu.y.y.y + -2. * d_ulu.y.yz * d_ull.y.yz + 2. * d_ulu.y.yz * d_llu.y.z.y + 2. * d_ull.y.yy * d_luu.y.yy + -2. * d_ull.y.yz * d_ulu.z.yy + -2. * d_ulu.y.zx * d_ull.z.xy + 2. * d_ulu.y.zx * d_llu.y.x.z + -2. * d_ulu.y.zy * d_ull.z.yy + 2. * d_ulu.y.zy * d_llu.y.y.z + -2. * d_ulu.y.zz * d_ull.z.yz + 2. * d_ulu.y.zz * d_llu.y.z.z + -2. * d_ull.y.zz * d_ulu.z.yz + 2. * d_ull.y.zz * d_luu.y.zz + -2. * d_uuu.z.xz * d_lll.z.xz + -2. * d_uul.z.xx * d_llu.z.x.x + -2. * d_uul.z.xy * d_llu.z.x.y + -2. * d_uul.z.xz * d_llu.z.x.z + -2. * d_uuu.z.yz * d_lll.z.yz + -2. * d_uul.z.yx * d_llu.z.y.x + -2. * d_uul.z.yy * d_llu.z.y.y + -2. * d_uul.z.yz * d_llu.z.y.z + 2. * d_uuu.z.zz * d_lll.z.zz + -2. * d_uul.z.zx * d_llu.z.z.x + -2. * d_uul.z.zy * d_llu.z.z.y + -2. * d_uul.z.zz * d_llu.z.z.z + 2. * d_ulu.z.xx * d_llu.z.x.x + 2. * d_ulu.z.xy * d_llu.z.y.x + 2. * d_ulu.z.xz * d_llu.z.z.x + 2. * d_ull.z.xx * d_luu.z.xx + -2. * d_ull.z.xz * d_ulu.z.zx + 2. * d_ulu.z.yx * d_llu.z.x.y + 2. * d_ulu.z.yy * d_llu.z.y.y + 2. * d_ulu.z.yz * d_llu.z.z.y + 2. * d_ull.z.yy * d_luu.z.yy + -2. * d_ull.z.yz * d_ulu.z.zy + 2. * d_ulu.z.zx * d_llu.z.x.z + 2. * d_ulu.z.zy * d_llu.z.y.z + -2. * d_ulu.z.zz * d_ull.z.zz + 2. * d_ulu.z.zz * d_llu.z.z.z + 2. * d_ull.z.zz * d_luu.z.zz + 3. * a_u.x * d_l.x + -3. * a_u.x * e_l.x + 3. * a_u.y * d_l.y + -3. * a_u.y * e_l.y + 3. * a_u.z * d_l.z + -3. * a_u.z * e_l.z + -4. * Z_l.x * conn_uul.x.xx + -4. * Z_l.x * conn_uul.x.yy + -4. * Z_l.x * conn_uul.x.zz + -4. * Z_l.y * conn_uul.y.xx + -4. * Z_l.y * conn_uul.y.yy + -4. * Z_l.y * conn_uul.y.zz + -4. * Z_l.z * conn_uul.z.xx + -4. * Z_l.z * conn_uul.z.yy + -4. * Z_l.z * conn_uul.z.zz + 4. * d_u.x * e_l.x + 4. * d_u.y * e_l.y + 4. * d_u.z * e_l.z + -4. * d_l.x * e_u.x + -4. * d_l.y * e_u.y + -4. * d_l.z * e_u.z + -4. * K_uu.xy * K_ll.xy + -4. * K_uu.xz * K_ll.xz + -4. * K_uu.yz * K_ll.yz + -4. * d_uuu.x.yy * d_lll.x.yy + -4. * d_uuu.x.zz * d_lll.x.zz + 4. * d_ull.x.xy * d_luu.x.xy + 4. * d_ull.x.xz * d_luu.x.xz + 4. * d_ull.x.yz * d_luu.x.yz + -4. * d_uuu.y.xx * d_lll.y.xx + -4. * d_uuu.y.zz * d_lll.y.zz + 4. * d_ull.y.xy * d_luu.y.xy + 4. * d_ull.y.xz * d_luu.y.xz + 4. * d_ull.y.yz * d_luu.y.yz + -4. * d_uuu.z.xx * d_lll.z.xx + -4. * d_uuu.z.yy * d_lll.z.yy + 4. * d_ull.z.xy * d_luu.z.xy + 4. * d_ull.z.xz * d_luu.z.xz + 4. * d_ull.z.yz * d_luu.z.yz + -4. * tr_K * Theta + -8. * Z_l.x * a_u.x + 8. * Z_l.x * e_u.x + -8. * Z_l.y * a_u.y + 8. * Z_l.y * e_u.y + -8. * Z_l.z * a_u.z + 8. * Z_l.z * e_u.z + -8. * d_uuu.x.yz * d_lll.x.yz + -8. * d_uuu.y.xz * d_lll.y.xz + -8. * d_uuu.z.xy * d_lll.z.xy + -32. * tmp2 + 6. * d_uuu.x.xy * d_lll.y.xx + 6. * d_uuu.x.xz * d_lll.z.xx + 6. * d_uuu.x.yy * d_lll.y.xy + 6. * d_uuu.x.yz * d_lll.y.xz + 6. * d_uuu.x.yz * d_lll.z.xy + 6. * d_uuu.x.zz * d_lll.z.xz + 6. * d_uuu.y.xx * d_lll.x.xy + 6. * d_uuu.y.xy * d_lll.x.yy + 6. * d_uuu.y.xz * d_lll.x.yz + 6. * d_uuu.y.xz * d_lll.z.xy + 6. * d_uuu.y.yz * d_lll.z.yy + 6. * d_uuu.y.zz * d_lll.z.yz + 6. * d_uuu.z.xx * d_lll.x.xz + 6. * d_uuu.z.xy * d_lll.x.yz + 6. * d_uuu.z.xy * d_lll.y.xz + 6. * d_uuu.z.xz * d_lll.x.zz + 6. * d_uuu.z.yz * d_lll.y.zz + 6. * d_uuu.z.yy * d_lll.y.yz)) / 4.;
		(deriv)->Z_l.x += alpha * (-a_l.x * K_ul.x.x + a_l.x * tr_K + -a_l.y * K_ul.y.x + -a_l.z * K_ul.z.x + d_l.x * K_ul.x.x + d_l.y * K_ul.y.x + d_l.z * K_ul.z.x + -K_uu.xx * d_lll.x.xx + -K_uu.yy * d_lll.x.yy + -K_uu.zz * d_lll.x.zz + -2. * Z_l.x * K_ul.x.x + -2. * Z_l.y * K_ul.y.x + -2. * Z_l.z * K_ul.z.x + -2. * a_l.x * Theta + -2. * K_uu.xy * d_lll.x.xy + -2. * K_uu.xz * d_lll.x.xz + -8. * S_l.x * M_PI + -2. * K_uu.yz * d_lll.x.yz);
		(deriv)->Z_l.y += alpha * (-a_l.x * K_ul.x.y + -a_l.y * K_ul.y.y + a_l.y * tr_K + -a_l.z * K_ul.z.y + d_l.x * K_ul.x.y + d_l.y * K_ul.y.y + d_l.z * K_ul.z.y + -K_uu.xx * d_lll.y.xx + -K_uu.yy * d_lll.y.yy + -K_uu.zz * d_lll.y.zz + -2. * Z_l.x * K_ul.x.y + -2. * Z_l.y * K_ul.y.y + -2. * Z_l.z * K_ul.z.y + -2. * a_l.y * Theta + -2. * K_uu.xy * d_lll.y.xy + -2. * K_uu.xz * d_lll.y.xz + -8. * S_l.y * M_PI + -2. * K_uu.yz * d_lll.y.yz);
		(deriv)->Z_l.z += alpha * (-a_l.x * K_ul.x.z + -a_l.y * K_ul.y.z + -a_l.z * K_ul.z.z + a_l.z * tr_K + d_l.x * K_ul.x.z + d_l.y * K_ul.y.z + d_l.z * K_ul.z.z + -K_uu.xx * d_lll.z.xx + -K_uu.yy * d_lll.z.yy + -K_uu.zz * d_lll.z.zz + -2. * Z_l.x * K_ul.x.z + -2. * Z_l.y * K_ul.y.z + -2. * Z_l.z * K_ul.z.z + -2. * a_l.z * Theta + -2. * K_uu.xy * d_lll.z.xy + -2. * K_uu.xz * d_lll.z.xz + -8. * S_l.z * M_PI + -2. * K_uu.yz * d_lll.z.yz);
	}
	<? if false then ?>//useShift
	{
		real const tmp1 = alpha * alpha;
		(deriv)->a_l.x += -a_l.x * b_ul.y.y + -a_l.x * b_ul.z.z + a_l.z * b_ul.z.x + a_l.y * b_ul.y.x;
		(deriv)->a_l.y += a_l.x * b_ul.x.y + -a_l.y * b_ul.x.x + a_l.z * b_ul.z.y + -a_l.y * b_ul.z.z;
		(deriv)->a_l.z += a_l.x * b_ul.x.z + a_l.y * b_ul.y.z + -a_l.z * b_ul.y.y + -a_l.z * b_ul.x.x;
		(deriv)->d_lll.x.xx += b_ul.y.x * dDelta_lll.y.xx + -b_ul.y.y * dDelta_lll.x.xx + -b_ul.z.z * dDelta_lll.x.xx + b_ul.z.x * dDelta_lll.z.xx;
		(deriv)->d_lll.x.xy += b_ul.y.x * dDelta_lll.y.xy + -b_ul.y.y * dDelta_lll.x.xy + -b_ul.z.z * dDelta_lll.x.xy + b_ul.z.x * dDelta_lll.z.xy;
		(deriv)->d_lll.x.xz += b_ul.y.x * dDelta_lll.y.xz + -b_ul.y.y * dDelta_lll.x.xz + -b_ul.z.z * dDelta_lll.x.xz + b_ul.z.x * dDelta_lll.z.xz;
		(deriv)->d_lll.x.yy += b_ul.y.x * dDelta_lll.y.yy + -b_ul.y.y * dDelta_lll.x.yy + -b_ul.z.z * dDelta_lll.x.yy + b_ul.z.x * dDelta_lll.z.yy;
		(deriv)->d_lll.x.yz += b_ul.y.x * dDelta_lll.y.yz + -b_ul.y.y * dDelta_lll.x.yz + -b_ul.z.z * dDelta_lll.x.yz + b_ul.z.x * dDelta_lll.z.yz;
		(deriv)->d_lll.x.zz += b_ul.y.x * dDelta_lll.y.zz + -b_ul.y.y * dDelta_lll.x.zz + -b_ul.z.z * dDelta_lll.x.zz + b_ul.z.x * dDelta_lll.z.zz;
		(deriv)->d_lll.y.xx += -b_ul.x.x * dDelta_lll.y.xx + b_ul.x.y * dDelta_lll.x.xx + -b_ul.z.z * dDelta_lll.y.xx + b_ul.z.y * dDelta_lll.z.xx;
		(deriv)->d_lll.y.xy += -b_ul.x.x * dDelta_lll.y.xy + b_ul.x.y * dDelta_lll.x.xy + -b_ul.z.z * dDelta_lll.y.xy + b_ul.z.y * dDelta_lll.z.xy;
		(deriv)->d_lll.y.xz += -b_ul.x.x * dDelta_lll.y.xz + b_ul.x.y * dDelta_lll.x.xz + -b_ul.z.z * dDelta_lll.y.xz + b_ul.z.y * dDelta_lll.z.xz;
		(deriv)->d_lll.y.yy += -b_ul.x.x * dDelta_lll.y.yy + b_ul.x.y * dDelta_lll.x.yy + -b_ul.z.z * dDelta_lll.y.yy + b_ul.z.y * dDelta_lll.z.yy;
		(deriv)->d_lll.y.yz += -b_ul.x.x * dDelta_lll.y.yz + b_ul.x.y * dDelta_lll.x.yz + -b_ul.z.z * dDelta_lll.y.yz + b_ul.z.y * dDelta_lll.z.yz;
		(deriv)->d_lll.y.zz += -b_ul.x.x * dDelta_lll.y.zz + b_ul.x.y * dDelta_lll.x.zz + -b_ul.z.z * dDelta_lll.y.zz + b_ul.z.y * dDelta_lll.z.zz;
		(deriv)->d_lll.z.xx += -b_ul.x.x * dDelta_lll.z.xx + b_ul.x.z * dDelta_lll.x.xx + b_ul.y.z * dDelta_lll.y.xx + -b_ul.y.y * dDelta_lll.z.xx;
		(deriv)->d_lll.z.xy += -b_ul.x.x * dDelta_lll.z.xy + b_ul.x.z * dDelta_lll.x.xy + b_ul.y.z * dDelta_lll.y.xy + -b_ul.y.y * dDelta_lll.z.xy;
		(deriv)->d_lll.z.xz += -b_ul.x.x * dDelta_lll.z.xz + b_ul.x.z * dDelta_lll.x.xz + b_ul.y.z * dDelta_lll.y.xz + -b_ul.y.y * dDelta_lll.z.xz;
		(deriv)->d_lll.z.yy += -b_ul.x.x * dDelta_lll.z.yy + b_ul.x.z * dDelta_lll.x.yy + b_ul.y.z * dDelta_lll.y.yy + -b_ul.y.y * dDelta_lll.z.yy;
		(deriv)->d_lll.z.yz += -b_ul.x.x * dDelta_lll.z.yz + b_ul.x.z * dDelta_lll.x.yz + b_ul.y.z * dDelta_lll.y.yz + -b_ul.y.y * dDelta_lll.z.yz;
		(deriv)->d_lll.z.zz += -b_ul.x.x * dDelta_lll.z.zz + b_ul.x.z * dDelta_lll.x.zz + b_ul.y.z * dDelta_lll.y.zz + -b_ul.y.y * dDelta_lll.z.zz;
		(deriv)->K_ll.xx += beta_u.x * partial_K_lll.xx.x + beta_u.y * partial_K_lll.xx.y + beta_u.z * partial_K_lll.xx.z + 2. * K_ll.xx * b_ul.x.x + 2. * K_ll.xz * b_ul.z.x + 2. * K_ll.xy * b_ul.y.x;
		(deriv)->K_ll.xy += beta_u.x * partial_K_lll.xy.x + beta_u.y * partial_K_lll.xy.y + beta_u.z * partial_K_lll.xy.z + K_ll.xx * b_ul.x.y + K_ll.xy * b_ul.x.x + K_ll.xy * b_ul.y.y + K_ll.xz * b_ul.z.y + K_ll.yz * b_ul.z.x + K_ll.yy * b_ul.y.x;
		(deriv)->K_ll.xz += beta_u.x * partial_K_lll.xz.x + beta_u.y * partial_K_lll.xz.y + beta_u.z * partial_K_lll.xz.z + K_ll.xx * b_ul.x.z + K_ll.xy * b_ul.y.z + K_ll.xz * b_ul.x.x + K_ll.xz * b_ul.z.z + K_ll.zz * b_ul.z.x + K_ll.yz * b_ul.y.x;
		(deriv)->K_ll.yy += beta_u.x * partial_K_lll.yy.x + beta_u.y * partial_K_lll.yy.y + beta_u.z * partial_K_lll.yy.z + 2. * K_ll.xy * b_ul.x.y + 2. * K_ll.yz * b_ul.z.y + 2. * K_ll.yy * b_ul.y.y;
		(deriv)->K_ll.yz += beta_u.x * partial_K_lll.yz.x + beta_u.y * partial_K_lll.yz.y + beta_u.z * partial_K_lll.yz.z + K_ll.xy * b_ul.x.z + K_ll.xz * b_ul.x.y + K_ll.yy * b_ul.y.z + K_ll.yz * b_ul.y.y + K_ll.zz * b_ul.z.y + K_ll.yz * b_ul.z.z;
		(deriv)->K_ll.zz += beta_u.x * partial_K_lll.zz.x + beta_u.y * partial_K_lll.zz.y + beta_u.z * partial_K_lll.zz.z + 2. * K_ll.xz * b_ul.x.z + 2. * K_ll.zz * b_ul.z.z + 2. * K_ll.yz * b_ul.y.z;
		(deriv)->Theta += -Theta * (b_ul.x.x + b_ul.y.y + b_ul.z.z);
		(deriv)->Z_l.x += -Z_l.x * b_ul.y.y + -Z_l.x * b_ul.z.z + Z_l.z * b_ul.z.x + Z_l.y * b_ul.y.x;
		(deriv)->Z_l.y += Z_l.x * b_ul.x.y + -Z_l.y * b_ul.x.x + Z_l.z * b_ul.z.y + -Z_l.y * b_ul.z.z;
		(deriv)->Z_l.z += Z_l.x * b_ul.x.z + Z_l.y * b_ul.y.z + -Z_l.z * b_ul.y.y + -Z_l.z * b_ul.x.x;
		(deriv)->beta_u.x += beta_u.x * b_ul.x.x + beta_u.y * b_ul.x.y + beta_u.z * b_ul.x.z + -a_u.x * tmp1 + 2. * e_u.x * tmp1 + -d_u.x * tmp1;
		(deriv)->beta_u.y += beta_u.x * b_ul.y.x + beta_u.y * b_ul.y.y + beta_u.z * b_ul.y.z + -a_u.y * tmp1 + 2. * e_u.y * tmp1 + -d_u.y * tmp1;
		(deriv)->beta_u.z += beta_u.x * b_ul.z.x + beta_u.y * b_ul.z.y + beta_u.z * b_ul.z.z + -a_u.z * tmp1 + 2. * e_u.z * tmp1 + -d_u.z * tmp1;
	}
	<? end ?>
	// END CUT
<? end ?>
//decay for 1st deriv hyperbolic state vars constraints:

	//turns out if you "if conv != 0" all these then skipping decay explodes soon in simulation steps, but time grows quickly, so it dies at a high t value

//// MODULE_DEPENDS: <?=calcFromGrad_a_l?>
	// a_x = log(alpha)_,x <=> a_x += eta (log(alpha)_,x - a_x)
	if (solver->a_convCoeff != 0.) {
		real3 const target_a_l = <?=calcFromGrad_a_l?>(solver, U);
<?
for i,xi in ipairs(xNames) do
?>		deriv->a_l.<?=xi?> += solver->a_convCoeff * (target_a_l.<?=xi?> - U->a_l.<?=xi?>);
<? end
?>	}

//// MODULE_DEPENDS: <?=calcFromGrad_d_lll?>
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	if (solver->d_convCoeff != 0.) {
		_3sym3 const target_d_lll = <?=calcFromGrad_d_lll?>(solver, U, cell);
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
?>		deriv->dDelta_lll.<?=xi?>.<?=xjk?> += solver->d_convCoeff * (target_d_lll.<?=xi?>.<?=xjk?> - d_lll.<?=xi?>.<?=xjk?>);
<? 	end
end
?>	}

<? if useKreissOligarDissipation then ?>
//// MODULE_DEPENDS: applyKreissOligar numberof
	// Kreiss-Oligar dissipation:
	if (solver->dissipationCoeff != 0.) {
		int fields[numIntStates] = {<?=require "ext.range"(0,eqn.numIntStates-1):concat", "?>};
		applyKreissOligar(solver, U, cell, deriv, fields, numberof(fields));
	}
<? end -- useKreissOligarDissipation ?>
<? end -- useAddSource ?>
}

//// MODULE_NAME: <?=constrainU?>
//// MODULE_DEPENDS: sym3sym3

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useConstrainU then ?>
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cons_t?> * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	
	real const det_gamma = sym3_det(U->gamma_ll);
	sym3 const gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real const tr_K = real3x3_trace(K_ul);							//K^k_k
	sym3 const KSq_ll = sym3_real3x3_to_sym3_mul(U->K_ll, K_ul);	//KSq_ij = K_ik K^k_j
	sym3 const K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);		//K^ij

//// MODULE_DEPENDS: <?=calc_d_lll?>
	_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);
	
	//d_llu = d_ij^k = d_ijl * gamma^lk
	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);
	
	//d_ull = d^i_jk = gamma^il d_ljk
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);

	//e_i = d^j_ji
	real3 const e_l = _3sym3_tr12(d_ull);

	//conn^k_ij = d_ij^k + d_ji^k - d^k_ij
	_3sym3 const conn_ull = {
<? for k,xk in ipairs(xNames) do
?>		.<?=xk?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>			.<?=xij?> = d_llu.<?=xi?>.<?=xj?>.<?=xk?> + d_llu.<?=xj?>.<?=xi?>.<?=xk?> - d_ull.<?=xk?>.<?=xij?>,
<? end
?>		},
<? end
?>	};
	
	//d_l.i = d_i = d_ij^j
	real3 const d_l = real3x3x3_tr23(d_llu);

	//partial_d_lll.ij.kl = d_kij,l = d_(k|(ij),|l)
	//so this object's indexes are rearranged compared to the papers
	sym3sym3 partial_d_llll;
<?
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
		local k,l,xk,xl = from6to3x3(kl)
?>	partial_d_llll.<?=xij?>.<?=xkl?> = 0.
<?		if l <= solver.dim then
?>	+ .5 * (	// 1/2 d_kij,l
		U[solver->stepsize.<?=xl?>].dDelta_lll.<?=xk?>.<?=xij?>
		//+ dHatR_lll.<?=xk?>.<?=xij?>	//TODO
		- U[-solver->stepsize.<?=xl?>].dDelta_lll.<?=xk?>.<?=xij?>
		//- dHatL_lll.<?=xk?>.<?=xij?>	//TODO
	) / (2. * solver->grid_dx.<?=xl?>)
<?
		end
		if k <= solver.dim then
?>	+ .5 * (
		U[solver->stepsize.<?=xk?>].dDelta_lll.<?=xl?>.<?=xij?>
		//+ dHatR_lll.<?=xk?>.<?=xij?>	//TODO
		- U[-solver->stepsize.<?=xk?>].dDelta_lll.<?=xl?>.<?=xij?>
		//- dHatL_lll.<?=xk?>.<?=xij?>	//TODO
	) / (2. * solver->grid_dx.<?=xk?>)
<?		end
?>	;
<?	end
end
?>

	// R_ll.ij := R_ij
	//	= gamma^kl (-gamma_ij,kl - gamma_kl,ij + gamma_ik,jl + gamma_jl,ik)
	//		+ conn^k_ij (d_k - 2 e_k)
	//		- 2 d^l_ki d^k_lj
	//		+ 2 d^l_ki d_lj^k
	//		+ d_il^k d_jk^l
	sym3 const R_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 0.
<? 	for k,xk in ipairs(xNames) do
?>			+ conn_ull.<?=xk?>.<?=xij?> * (d_l.<?=xk?> - 2. * e_l.<?=xk?>)
<?		for l,xl in ipairs(xNames) do
?>			+ 2. * d_ull.<?=xl?>.<?=sym(k,i)?> * (d_llu.<?=xl?>.<?=xj?>.<?=xk?> - d_ull.<?=xk?>.<?=sym(l,j)?>)
			+ d_llu.<?=xi?>.<?=xl?>.<?=xk?> * d_llu.<?=xj?>.<?=xk?>.<?=xl?>
			+ gamma_uu.<?=sym(k,l)?> * (
				- partial_d_llll.<?=xij?>.<?=sym(k,l)?>
				- partial_d_llll.<?=sym(k,l)?>.<?=xij?>
				+ partial_d_llll.<?=sym(i,k)?>.<?=sym(j,l)?>
				+ partial_d_llll.<?=sym(j,l)?>.<?=sym(i,k)?>
			)
<? 		end
	end
?>		,
<? end
?>	};

	//scaled down by 1/8 to match B&S BSSNOK equations ... maybe I'll scale theirs up by 8 ...
	//B&S eqn 2.125 ... divded by two
	//Alcubierre eqn 2.5.9, also Alcubierre 2.4.10 divided by two
	//H = 1/2 (R + K^2 - K_ij K^ij) - 8 pi rho
	real const R = sym3_dot(R_ll, gamma_uu);
	real const tr_KSq = sym3_dot(KSq_ll, gamma_uu);
	U->H = .5 * (R + tr_K * tr_K - tr_KSq) <?
if eqn.useStressEnergyTerms then ?>
	- 8. * M_PI * U->rho <?
end ?>;

<?=eqn:makePartial1"K_ll"?>

	/*
	momentum constraint
	Alcubierre eqn 2.4.11
	M^i = K^ij_;j
		- gamma^ij K_,j
		- 8 pi S^i
	M^i = gamma^im gamma^jn (
			K_mn,j
			- conn^k_mj K_kn
			- conn^k_nj K_mk
		)
		- gamma^ij (gamma^mn_,j K_mn + gamma^mn K_mn,j)
		- 8 pi S^i
	M^i = gamma^ij (
			gamma^mn (
				K_jn,m
				- K_mn,j
				- conn^k_jm K_kn
				- conn^k_nm K_jk
			)
			+ 2 d_jmn K^mn
		)
		- 8 pi S^i
	*/
<? for i,xi in ipairs(xNames) do
?>	U->M_u.<?=xi?> =
<?	for j,xj in ipairs(xNames) do
		for m,xm in ipairs(xNames) do
			for n,xn in ipairs(xNames) do
?>		+ gamma_uu.<?=sym(i,j)?> * (
			gamma_uu.<?=sym(m,n)?> * (0.
				+ partial_K_lll[<?=m-1?>].<?=sym(j,n)?>
				- partial_K_lll[<?=j-1?>].<?=sym(m,n)?>
<?				for k,xk in ipairs(xNames) do
?>				- conn_ull.<?=xk?>.<?=sym(j,m)?> * U->K_ll.<?=sym(k,n)?>
				- conn_ull.<?=xk?>.<?=sym(n,m)?> * U->K_ll.<?=sym(j,k)?>
<?				end
?>			)
			+ 2. * d_lll.<?=xj?>.<?=sym(m,n)?> * K_uu.<?=sym(m,n)?>
		)
<?			end
		end
	end
	if eqn.useStressEnergyTerms then
?>		- 8. * M_PI * U->S_u.<?=xi?>
<?	end
?>	;
<? end
?>

	U->alpha = max(U->alpha, solver->alphaMin);
<? end -- useConstrainU ?>
}
