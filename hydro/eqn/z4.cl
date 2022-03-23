<?
local useConstrainU = true -- constrains alpha to alphamin and calcs H and M^i
local useAddSource = true
local useKreissOligarDissipation = true	-- depends on useAddSource

local has_beta_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "beta_u" end)
local has_betaLap_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "betaLap_u" end)
local has_b_ul = eqn.consStruct.vars:find(nil, function(var) return var.name == "b_ul" end)
local has_B_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end)
?>

<? -- ok so ideally you should be able to subtract out the background metric ...
?>

//// MODULE_NAME: <?=calc_gammaHat_ll?>
//// MODULE_DEPENDS: sym3

static inline sym3 <?=calc_gammaHat_ll?>(real3 const pt) {
<? if false then -- get background grid separation working
?>
//// MODULE_DEPENDS: <?=coord_gHol_ll?>
	return coord_gHol_ll(pt);
<? else ?>
	return sym3_zero;
<? end ?>
}

//// MODULE_NAME: <?=calc_dHat_lll?>
//// MODULE_DEPENDS: _3sym3

static inline _3sym3 <?=calc_dHat_lll?>(real3 const pt) {
<? if false then 	-- get background grid separation working?
?>
//// MODULE_DEPENDS: <?=coord_partial_gHol_lll?>
	return _3sym3_real_mul(coord_partial_gHol_lll(pt), .5);
<? else ?>
	return _3sym3_zero;
<? end ?>
}

//// MODULE_NAME: <?=calc_gamma_ll?>
//// MODULE_DEPENDS: <?=cons_t?> <?=calc_gammaHat_ll?>

#define /*sym3*/ <?=calc_gamma_ll?>(\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt\
) (\
	sym3_add(\
		(U)->gammaDelta_ll,\
		<?=calc_gammaHat_ll?>(pt)\
	)\
)

//// MODULE_NAME: <?=calc_gamma_uu?>
//// MODULE_DEPENDS: <?=cons_t?> <?=calc_gamma_ll?>

#define /*sym3*/ <?=calc_gamma_uu?>(\
	/*global <?=cons_t?> const * const*/ U,\
	/*real3 const*/ pt\
) (sym3_inv_nodet(<?=calc_gamma_ll?>(U, pt)))

//// MODULE_NAME: <?=calc_d_lll?>

static inline _3sym3 <?=calc_d_lll?>(
	global <?=cons_t?> const * const U,
	real3 const pt
) {
//// MODULE_DEPENDS: <?=calc_dHat_lll?>
	_3sym3 const dHat_lll = <?=calc_dHat_lll?>(pt);
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
//// MODULE_DEPENDS: <?=calc_gamma_ll?>
		global <?=cons_t?> const * const UL = U - solver->stepsize.<?=xi?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.<?=xi?>;
		global <?=cell_t?> const * const cellL = cell - solver->stepsize.<?=xi?>;
		global <?=cell_t?> const * const cellR = cell + solver->stepsize.<?=xi?>;
		sym3 const gammaR_ll = <?=calc_gamma_ll?>(UR, cellR->pos);
		sym3 const gammaL_ll = <?=calc_gamma_ll?>(UL, cellL->pos);
<? 	for jk,xjk in ipairs(symNames) do
?>		d_lll.<?=xi?>.<?=xjk?> = .5 * (gammaR_ll.<?=xjk?> - gammaL_ll.<?=xjk?>) / (2. * solver->grid_dx.s<?=i-1?>);
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

	//d_kij = 1/2 γ_ij,k
	//Δd_kij + ^d_kij = 1/2 (Δγ_ij,k + ^γ_ij,k)
	//^d_kij = 1/2 ^γ_ij,k
	// => Δd_kij = 1/2 Δγ_ij,k
//// MODULE_DEPENDS: <?=calc_dHat_lll?>
	_3sym3 const dHat_lll = <?=calc_dHat_lll?>(cell->pos);
//// MODULE_DEPENDS: <?=calcFromGrad_d_lll?>
	_3sym3 const d_lll = <?=calcFromGrad_d_lll?>(solver, U, cell);
	U->dDelta_lll = _3sym3_sub(d_lll, dHat_lll);

<? if has_b_ul then ?>
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
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
//// MODULE_DEPENDS: <?=coordMap?>
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
	// ^γ_IJ = δ_IJ
	// γ_ij = e_i^I e_j^J (ε_IJ + ^γ_IJ) / W^2
	sym3 const gammaBar_LL = sym3_add(epsilon_LL, sym3_ident);
	sym3 const gamma_LL = sym3_real_mul(gammaBar_LL, 1. / (W*W));
	sym3 const gamma_ll = sym3_rescaleToCoord_LL(gamma_LL, x);
//// MODULE_DEPENDS: <?=calc_gammaHat_ll?>
	sym3 const gammaHat_ll = <?=calc_gammaHat_ll?>(x);
	U->gammaDelta_ll = sym3_sub(gamma_ll, gammaHat_ll);

	// K_ij = e_i^I e_j^J (_A_IJ + _γ_IJ K/3) / W^2
	U->K_ll = sym3_rescaleToCoord_LL(
		sym3_add(
			sym3_real_mul(ABar_LL, 1. / (W*W)),
			sym3_real_mul(gamma_LL, K / 3.)
		), x);

	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if has_beta_u then
?>	U->beta_u = real3_rescaleToCoord_U(beta_U, x);
<? end
if has_betaLap_u then
?>	U->betaLap_u = real3_zero;
<? end
if has_B_u then
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
//// MODULE_DEPENDS: <?=solver_t?> <?=initCond_t?> <?=cons_t?> <?=cell_t?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real3 beta_u = real3_zero;
//// MODULE_DEPENDS: <?=coord_gHol_ll?>
	sym3 gamma_ll = coord_gHol_ll(x);

	sym3 K_ll = sym3_zero;

	//TODO more stress-energy vars
	real rho = 0.;

	<?=initCode()?>

	*U = (<?=cons_t?>){.ptr={ 0. / 0. }};
	
	U->alpha = alpha;
	
//// MODULE_DEPENDS: <?=calc_gammaHat_ll?>
	sym3 const gammaHat_ll = <?=calc_gammaHat_ll?>(x);
	U->gammaDelta_ll = sym3_sub(gamma_ll, gammaHat_ll);
	
	U->K_ll = K_ll;

	//Z_u n^u = 0
	//Θ = α n_u Z^u = α Z^u
	//for n_a = (-α, 0)
	//n^a_l = (1/α, -β^i/α)
	//(Z_t - Z_i β^i) / α = Θ ... = ?
	//Z^t n_t + Z^i n_i = -α Z^t = Θ
	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if has_beta_u then
?>	U->beta_u = beta_u;
<? end
if has_betaLap_u then
?>	U->betaLap_u = real3_zero;
<? end
if has_B_u then
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

// irl this should be zero, but for now gammaHat_ll is set to zero so it does have a diff
//// MODULE_DEPENDS: <?=coord_gHol_ll?>
	sym3 const gamma_ll = coord_gHol_ll(x);
//// MODULE_DEPENDS: <?=calc_gammaHat_ll?>
	sym3 const gammaHat_ll = <?=calc_gammaHat_ll?>(x);
	U->gammaDelta_ll = sym3_sub(gamma_ll, gammaHat_ll);
	
	(U)->a_l = real3_zero;
	(U)->dDelta_lll = _3sym3_zero;
	(U)->K_ll = sym3_zero;
	(U)->Theta = 0.;
	(U)->Z_l = real3_zero;

<? if has_beta_u then
?>	(U)->beta_u = beta_u;
<? end
if has_betaLap_u then
?>	(U)->betaLap_u = real3_zero;
<? end
if has_B_u then
?>	(U)->B_u = real3_zero;
<? end
if has_b_ul then
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
//// MODULE_DEPENDS: <?=calc_gamma_ll?>

#define <?=calcDTCell?>(\
	/*global real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell\
) {\
	/* the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma_ij) is only called once */\
	real const f_alphaSq = calc_f_alphaSq(U->alpha);\
	sym3 const gamma_ll = <?=calc_gamma_ll?>(U, (cell)->pos);\
	real const det_gamma = sym3_det(gamma_ll);\
	real const alpha_sqrt_f = sqrt(f_alphaSq);\
\
	<? for side=0,solver.dim-1 do ?>{\
\
		<? if side == 0 then ?>\
		real const gammaUjj = (gamma_ll.yy * gamma_ll.zz - gamma_ll.yz * gamma_ll.yz) / det_gamma;\
		<? elseif side == 1 then ?>\
		real const gammaUjj = (gamma_ll.xx * gamma_ll.zz - gamma_ll.xz * gamma_ll.xz) / det_gamma;\
		<? elseif side == 2 then ?>\
		real const gammaUjj = (gamma_ll.xx * gamma_ll.yy - gamma_ll.xy * gamma_ll.xy) / det_gamma;\
		<? end ?>	\
		real const sqrt_gammaUjj = sqrt(gammaUjj);\
		real const lambdaLight = sqrt_gammaUjj * U->alpha;\
		real const lambdaGauge = sqrt_gammaUjj * alpha_sqrt_f;\
		real const lambda = (real)max(lambdaGauge, lambdaLight);\
\
		<? if has_beta_u then ?>\
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
//// MODULE_DEPENDS: <?=calc_dHat_lll?>
//// MODULE_DEPENDS: <?=calc_gamma_ll?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real const f_alpha = calc_f_alpha((U)->alpha);\
\
	sym3 gamma_ll = <?=calc_gamma_ll?>(U, (cell)->pos);\
	real const det_gamma = sym3_det(gamma_ll);\
\
	real alpha = (U)->alpha;\
	real Theta = (U)->Theta;\
\
	real3 Z_l = (U)->Z_l;\
	real3 a_l = (U)->a_l;\
	sym3 K_ll = (U)->K_ll;\
	_3sym3 dDelta_lll = (U)->dDelta_lll;\
	_3sym3 dHat_lll = <?=calc_dHat_lll?>((cell)->pos);\
<? if has_beta_u then --\
?>	real3 beta_u = (U)->beta_u;\
<? end --\
if has_b_ul then --\
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
<? if has_beta_u then --\
?>	beta_u = real3_rotateFrom(beta_u);\
<? end --\
if has_b_ul then --\
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
	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);\
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);\
	_3sym3 const conn_ull = conn_ull_from_d_llu_d_ull(d_llu, d_ull);\
	real3 const e_l = _3sym3_tr12(d_ull);\
	real3 const d_l = real3x3x3_tr23(d_llu);\
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);\
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);\
	real3 const Z_u = sym3_real3_mul(gamma_uu, Z_l);\
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, K_ll);\
	real const tr_K = real3x3_trace(K_ul);\
	sym3 const dHat_t_ll = sym3_zero;\
<? 	if eqn.useShift == "HarmonicShiftParabolic" then ?>\
	real3 const a_u = sym3_real3_mul(gamma_uu, a_l);\
	sym3 const b_ll = sym3_real3x3_to_sym3_mul(gamma_ll, b_ul);\
<? 	end ?>\
	/* BEGIN CUT from symmath/tests/output/Z4.html */\
	{\
		(resultFlux)->alpha = 0.;\
		(resultFlux)->gammaDelta_ll.xx = 0.;\
		(resultFlux)->gammaDelta_ll.xy = 0.;\
		(resultFlux)->gammaDelta_ll.xz = 0.;\
		(resultFlux)->gammaDelta_ll.yy = 0.;\
		(resultFlux)->gammaDelta_ll.yz = 0.;\
		(resultFlux)->gammaDelta_ll.zz = 0.;\
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
		(resultFlux)->K_ll.xx = alpha * (a_l.x + d_l.x + -2. * Z_l.x + -conn_ull.x.xx);\
		(resultFlux)->K_ll.xy = (alpha * (a_l.y + d_l.y + -2. * conn_ull.x.xy + -2. * Z_l.y)) / 2.;\
		(resultFlux)->K_ll.xz = (alpha * (a_l.z + d_l.z + -2. * conn_ull.x.xz + -2. * Z_l.z)) / 2.;\
		(resultFlux)->K_ll.yy = -conn_ull.x.yy * alpha;\
		(resultFlux)->K_ll.yz = -conn_ull.x.yz * alpha;\
		(resultFlux)->K_ll.zz = -conn_ull.x.zz * alpha;\
		(resultFlux)->Theta = alpha * (d_u.x + -Z_u.x + -e_u.x);\
		(resultFlux)->Z_l.x = alpha * (tr_K + -K_ul.x.x + -Theta);\
		(resultFlux)->Z_l.y = -K_ul.x.y * alpha;\
		(resultFlux)->Z_l.z = -K_ul.x.z * alpha;\
	}\
	<? if eqn.useShift ~= "none" then ?>\
	{\
		real const tmp1 = alpha * alpha;\
		(resultFlux)->a_l.x += -beta_u.x * a_l.x;\
		(resultFlux)->a_l.y += -beta_u.x * a_l.y;\
		(resultFlux)->a_l.z += -beta_u.x * a_l.z;\
		(resultFlux)->dDelta_lll.x.xx += -(b_ll.xx + beta_u.x * dHat_lll.x.xx + beta_u.x * dDelta_lll.x.xx + beta_u.z * dHat_lll.z.xx + beta_u.y * dHat_lll.y.xx);\
		(resultFlux)->dDelta_lll.x.xy += -(b_ll.xy + b_ll.xy + 2. * beta_u.x * dHat_lll.x.xy + 2. * beta_u.x * dDelta_lll.x.xy + 2. * beta_u.z * dHat_lll.z.xy + 2. * beta_u.y * dHat_lll.y.xy) / 2.;\
		(resultFlux)->dDelta_lll.x.xz += -(b_ll.xz + b_ll.xz + 2. * beta_u.x * dHat_lll.x.xz + 2. * beta_u.x * dDelta_lll.x.xz + 2. * beta_u.z * dHat_lll.z.xz + 2. * beta_u.y * dHat_lll.y.xz) / 2.;\
		(resultFlux)->dDelta_lll.x.yy += -(b_ll.yy + beta_u.x * dHat_lll.x.yy + beta_u.x * dDelta_lll.x.yy + beta_u.z * dHat_lll.z.yy + beta_u.y * dHat_lll.y.yy);\
		(resultFlux)->dDelta_lll.x.yz += -(b_ll.yz + b_ll.yz + 2. * beta_u.x * dHat_lll.x.yz + 2. * beta_u.x * dDelta_lll.x.yz + 2. * beta_u.z * dHat_lll.z.yz + 2. * beta_u.y * dHat_lll.y.yz) / 2.;\
		(resultFlux)->dDelta_lll.x.zz += -(b_ll.zz + beta_u.x * dHat_lll.x.zz + beta_u.x * dDelta_lll.x.zz + beta_u.z * dHat_lll.z.zz + beta_u.y * dHat_lll.y.zz);\
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
		<? if eqn.useShift == "HarmonicShiftParabolic" then ?>\
		(resultFlux)->b_ul.x.x += a_u.x * tmp1 + d_u.x * tmp1 + -beta_u.x * b_ul.x.x + -beta_u.y * b_ul.x.y + -2. * e_u.x * tmp1 + -beta_u.z * b_ul.x.z;\
		(resultFlux)->b_ul.y.x += a_u.y * tmp1 + d_u.y * tmp1 + -beta_u.x * b_ul.y.x + -beta_u.y * b_ul.y.y + -2. * e_u.y * tmp1 + -beta_u.z * b_ul.y.z;\
		(resultFlux)->b_ul.z.x += a_u.z * tmp1 + d_u.z * tmp1 + -beta_u.x * b_ul.z.x + -beta_u.y * b_ul.z.y + -2. * e_u.z * tmp1 + -beta_u.z * b_ul.z.z;\
		<? elseif eqn.useShift == "HarmonicShiftHyperbolic" then ?>\
		(resultFlux)->b_ul.x.x += -(B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y);\
		(resultFlux)->b_ul.y.x += -(B_u.y + beta_u.x * b_ul.y.x + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y);\
		(resultFlux)->b_ul.z.x += -(B_u.z + beta_u.x * b_ul.z.x + beta_u.z * b_ul.z.z + beta_u.y * b_ul.z.y);\
		<? end ?>\
	}\
	<? end ?>/* eqn.useShift ~= "none" */\
	/* END CUT */\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
\
	(resultFlux)->Z_l = real3_rotateFrom((resultFlux)->Z_l, n_l);\
	(resultFlux)->a_l = real3_rotateFrom((resultFlux)->a_l, n_l);\
	(resultFlux)->gammaDelta_ll = sym3_rotateFrom((resultFlux)->gammaDelta_ll, n_l);\
	(resultFlux)->K_ll = sym3_rotateFrom((resultFlux)->K_ll, n_l);\
	(resultFlux)->dDelta_lll = _3sym3_rotateFrom((resultFlux)->dDelta_lll, n_l);\
<? if has_beta_u then --\
?>	(resultFlux)->beta_u = real3_rotateFrom((resultFlux)->beta_u, n_l);\
<? end --\
if has_b_ul then --\
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
		(resultFlux)->gammaDelta_ll = sym3_swap<?=side?>((resultFlux)->gammaDelta_ll);\
		(resultFlux)->K_ll = sym3_swap<?=side?>((resultFlux)->K_ll);\
		(resultFlux)->dDelta_lll = _3sym3_swap<?=side?>((resultFlux)->dDelta_lll);\
<? if has_beta_u then --\
?>		(resultFlux)->beta_u = real3_swap<?=side?>((resultFlux)->beta_u);\
<? end --\
if has_B_u then --\
?>		(resultFlux)->B_u = real3_swap<?=side?>((resultFlux)->B_u);\
<? end --\
if has_b_ul then --\
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
<? if has_betaLap_u then --\
?>	(resultFlux)->betaLap_u = real3_zero;\
<? end --\
?>\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=solver_t?> <?=eigen_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>
//// MODULE_DEPENDS: <?=calc_gamma_ll?>
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
	sym3 const gammaR_ll = <?=calc_gamma_ll?>(UR, (cellR)->pos);\
	sym3 const gammaL_ll = <?=calc_gamma_ll?>(UL, (cellL)->pos);\
	sym3 const avg_gamma = sym3_real_mul(sym3_add(gammaL_ll, gammaR_ll), .5);\
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
<? if has_beta_u then ?>\
	(resultEig)->beta_u = real3_real_mul(real3_add((UL)->beta_u, (UR)->beta_u), .5);\
<? end ?>\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: <?=range_t?> <?=normal_t?> cons_t <?=initCond_codeprefix?>
//// MODULE_DEPENDS: <?=calc_gamma_ll?>
// not used anymore, replaced in calcDT by eqn:consMinWaveCode/eqn:consMaxWaveCode eigenvalue inlining

#define <?=calcCellMinMaxEigenvalues?>(\
	/*<?=range_t?> * const */result,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	sym3 const gamma_ll = <?=calc_gamma_ll?>(U, pt);\
	real const det_gamma = sym3_det(gamma_ll);\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
	sym3 const gamma_uu = sym3_inv(gamma_ll, det_gamma);\
	real3 const n_l = normal_l1(n);\
	real const gammaUnn = real3_weightedLenSq(n_l, gamma_uu);\
<? else ?>\
	real gammaUnn = 0./0.;\
	if (n.side == 0) {\
		gammaUnn = (gamma_ll.yy * gamma_ll.zz - gamma_ll.yz * gamma_ll.yz) / det_gamma;\
	} else if (n.side == 1) {\
		gammaUnn = (gamma_ll.xx * gamma_ll.zz - gamma_ll.xz * gamma_ll.xz) / det_gamma;\
	} else if (n.side == 2) {\
		gammaUnn = (gamma_ll.xx * gamma_ll.yy - gamma_ll.xy * gamma_ll.xy) / det_gamma;\
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
	<? if has_beta_u then ?>\
	lambdaMin -= normal_vecDotN1(n, (U)->beta_u);\
	lambdaMax -= normal_vecDotN1(n, (U)->beta_u);\
	<? end ?>\
\
	(result)->min = lambdaMin;\
	(result)->max = lambdaMax;\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> <?=initCond_codeprefix?>
//// MODULE_DEPENDS: <?=calc_gamma_uu?>

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
	(resultEig)->gamma_uu = <?=calc_gamma_uu?>(U, pt);\
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
	<? if has_beta_u then ?>\
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
	/*  mind you the 'm' in that form is for α_,t = -α^2 (f K - m Θ) */\
	/*  while in more modern Z4 papers it is α_,t = -α^2 f (K - m Θ) */\
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

	//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r ρ (D-)^r / 2^(2r)
	//...for r=2... -σ h^3 (D+)^2 ρ (D-)^2 / 16 ... and ρ=1, except ρ=0 at borders maybe.
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

//// MODULE_DEPENDS: <?=calc_gamma_ll?>
	sym3 const gamma_ll = <?=calc_gamma_ll?>(U, cell->pos);
	real const det_gamma = sym3_det(gamma_ll);
	sym3 const gamma_uu = sym3_inv(gamma_ll, det_gamma);

	real3 const S_l = real3_zero;
	sym3 const S_ll = sym3_zero;
	real const S = 0.;
	real const rho = 0.;

<? if false then ?>//hand-rolled
	
	// source terms
	
	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real const trK = real3x3_trace(K_ul);							//K^k_k

//// MODULE_DEPENDS: <?=calc_d_lll?>
	_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);

	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);	//d_llu = d_ij^k = d_ijl * γ^lk
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);		//d_ull = d^i_jk = γ^il d_ljk
	real3 const e_l = _3sym3_tr12(d_ull);						//e_l = e_i = d^j_ji
	_3sym3 const conn_ull = conn_ull_from_d_llu_d_ull(d_llu, d_ull);		//Γ^k_ij = d_ij^k + d_ji^k - d^k_ij	
	real3 const d_l = real3x3x3_tr23(d_llu);			/* d_l = d_i = d_ij^j */	
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 const Z_u = sym3_real3_mul(gamma_uu, U->Z_l);
	_3sym3 const d_luu = sym3_real3x3x3_mul2_to_3sym3(gamma_uu, d_llu);		/* d_luu = d_i^jk = γ^jl d_il^k */
	real const f_alphaSq = calc_f_alphaSq(U->alpha);
	
	/* α_,t = shift terms - α^2 f (γ^ij K_ij - m Θ) */
	deriv->alpha += -f_alphaSq * (trK - solver->m * U->Theta);
	
	sym3 const dt_gammaHat_ll = sym3_zero;
	
	/* γ_ij,t = shift terms - 2 α K_ij */
	/* Δγ_ij,t = shift terms - 2 α K_ij - ^γ_ij,t */
	deriv->gammaDelta_ll = sym3_add(
		deriv->gammaDelta_ll,
		sym3_sub(
			sym3_real_mul(U->K_ll, -2. * U->alpha),
			dt_gammaHat_ll
		)
	);

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
?> - 8. * M_PI * (S_ll.<?=xij?> - .5 * gamma_ll.<?=xij?> * (S - rho)));
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


<? if eqn.useShift == "HarmonicShiftParabolic" then ?>
	
	real const tr_b = real3x3_trace(U->b_ul);
	real3 const a_u = sym3_real3_mul(gamma_uu, U->a_l);

	/* α_,t += β^i α a_i */
<? for k,xk in ipairs(xNames) do
?>	deriv->alpha += U->alpha * real3_dot(U->a_l, U->beta_u);
<? end
?>

	sym3 const dHat_t_ll = sym3_zero;\
	real3x3 const b_ll = sym3_real3x3_mul(gamma_ll, U->b_ul);

	/* γ_ij += 2 β^k d_kij + b_ij + b_ji */
	/* so Δγ_ij += b_ij + b_ji + 2 β^k Δd_kij - 2 ^d_tij */
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	deriv->gammaDelta_ll.<?=xij?> += 0.
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
	/* so Δd_kij,t += b^l_k Δd_lijk - b^l_l Δd_kij */
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

	/* Θ_,t += b^m_m Θ */
	deriv->Theta += -tr_b * U->Theta;

	/* Z_k,t += b^m_k Z_m - b^m_m Z_k */
<? for k,xk in ipairs(xNames) do
?>	deriv->Z_l.<?=xk?> += 0.
		- tr_b * U->Z_l.<?=xk?>
<? 	for m,xm in ipairs(xNames) do
?>		+ U->Z_l.<?=xm?> * U->b_ul.<?=xm?>.<?=xk?>
<? 	end
?>	;
<? end
?>

	/* β^l,t += α^2 (2 e^l - d^l - a^l) + β^m b^l_m */
<? for l,xl in ipairs(xNames) do
?>	deriv->beta_u.<?=xl?> += U->alpha * U->alpha * (2. * e_u.<?=xl?> - d_u.<?=xl?> - a_u.<?=xl?>)
<?	for m,xm in ipairs(xNames) do
?>		+ U->beta_u.<?=xm?> * U->b_ul.<?=xl?>.<?=xm?>
<?	end
?>	;
<? end
?>

<? end -- eqn.useShift == "HarmonicShiftParabolic" ?>


<? end ?>
<? if true then ?>
	//code-generated from symmath/tests/Z4.lua

	real const alpha = U->alpha;
	real3 const a_l = U->a_l;
	sym3 const K_ll = U->K_ll;
	real const Theta = U->Theta;
	real3 const Z_l = U->Z_l;

	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, K_ll);			//K^i_j
	real const tr_K = real3x3_trace(K_ul);							//K^k_k
	sym3 const K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);		//K^ij

//// MODULE_DEPENDS: <?=calc_d_lll?>
	_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);
	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);
	_3sym3 const d_luu = sym3_real3x3x3_mul2_to_3sym3(gamma_uu, d_llu);
	_3sym3 d_uuu = sym3_3sym3_mul(gamma_uu, d_luu);
	real3 const e_l = _3sym3_tr12(d_ull);
	real3 const d_l = real3x3x3_tr23(d_llu);
	real3 const a_u = sym3_real3_mul(gamma_uu, a_l);
	real3 const d_u = sym3_real3_mul(gamma_uu, d_l);
	real3 const e_u = sym3_real3_mul(gamma_uu, e_l);
	real3 const Z_u = sym3_real3_mul(gamma_uu, Z_l);
	_3sym3 const conn_ull = conn_ull_from_d_llu_d_ull(d_llu, d_ull);	//Γ^k_ij = d_ij^k + d_ji^k - d^k_ij
	sym3 const dHat_t_ll = sym3_zero;
	
	real const f_alpha = calc_f_alpha((U)->alpha);

	<? if eqn.useShift ~= "none" then ?>
	<? if has_beta_u then ?>
	real3 const beta_u = U->beta_u;
	<? end ?>
	<? if has_b_ul then ?>
	real3x3 const b_ul = U->b_ul;
	real const tr_b = real3x3_trace(b_ul);
	sym3 b_ll = sym3_real3x3_to_sym3_mul(gamma_ll, b_ul);
	<? end ?>
	<? if has_B_u then ?>
	real3 const B_u = U->B_u;
	<? end ?>
	_3sym3 const dDelta_lll = U->dDelta_lll;
	<? end ?>

	// BEGIN CUT from symmath/tests/output/Z4.html
	{
		real const tmp1 = S * M_PI;
		real const tmp2 = rho * M_PI;
		(deriv)->alpha += f_alpha * alpha * (2. * Theta + -tr_K);
		(deriv)->gammaDelta_ll.xx += -2. * (dHat_t_ll.xx + K_ll.xx * alpha);
		(deriv)->gammaDelta_ll.xy += -2. * (dHat_t_ll.xy + K_ll.xy * alpha);
		(deriv)->gammaDelta_ll.xz += -2. * (dHat_t_ll.xz + K_ll.xz * alpha);
		(deriv)->gammaDelta_ll.yy += -2. * (dHat_t_ll.yy + K_ll.yy * alpha);
		(deriv)->gammaDelta_ll.yz += -2. * (dHat_t_ll.yz + K_ll.yz * alpha);
		(deriv)->gammaDelta_ll.zz += -2. * (dHat_t_ll.zz + K_ll.zz * alpha);
		(deriv)->K_ll.xx += alpha * (-conn_ull.x.xx * conn_ull.x.xx + -conn_ull.y.xy * conn_ull.y.xy + -conn_ull.z.xz * conn_ull.z.xz + a_l.x * d_l.x + d_l.x * conn_ull.x.xx + d_l.y * conn_ull.y.xx + d_l.z * conn_ull.z.xx + K_ll.xx * tr_K + -2. * Z_l.x * a_l.x + -2. * Z_l.x * conn_ull.x.xx + -2. * Z_l.y * conn_ull.y.xx + -2. * Z_l.z * conn_ull.z.xx + -2. * K_ul.x.x * K_ll.xx + -2. * K_ul.x.y * K_ll.xy + -2. * K_ul.x.z * K_ll.xz + -2. * K_ll.xx * Theta + -2. * conn_ull.x.xy * conn_ull.y.xx + -2. * conn_ull.x.xz * conn_ull.z.xx + -2. * conn_ull.y.xz * conn_ull.z.xy + -8. * S_ll.xx * M_PI + -4. * gamma_ll.xx * tmp2 + 4. * gamma_ll.xx * tmp1);
		(deriv)->K_ll.xy += (alpha * (a_l.x * d_l.y + a_l.y * d_l.x + -2. * Z_l.x * a_l.y + -2. * Z_l.y * a_l.x + 2. * d_l.x * conn_ull.x.xy + 2. * d_l.y * conn_ull.y.xy + 2. * d_l.z * conn_ull.z.xy + 2. * K_ll.xy * tr_K + -2. * conn_ull.x.xx * conn_ull.x.xy + -2. * conn_ull.x.xy * conn_ull.y.xy + -2. * conn_ull.x.xz * conn_ull.z.xy + -2. * conn_ull.x.yy * conn_ull.y.xx + -2. * conn_ull.x.yz * conn_ull.z.xx + -2. * conn_ull.y.xy * conn_ull.y.yy + -2. * conn_ull.y.xz * conn_ull.z.yy + -2. * conn_ull.y.yz * conn_ull.z.xy + -2. * conn_ull.z.xz * conn_ull.z.yz + -4. * Z_l.x * conn_ull.x.xy + -4. * Z_l.y * conn_ull.y.xy + -4. * Z_l.z * conn_ull.z.xy + -4. * K_ul.x.x * K_ll.xy + -4. * K_ul.x.y * K_ll.yy + -4. * K_ul.x.z * K_ll.yz + -4. * K_ll.xy * Theta + -16. * S_ll.xy * M_PI + -8. * gamma_ll.xy * tmp2 + 8. * gamma_ll.xy * tmp1)) / 2.;
		(deriv)->K_ll.xz += (alpha * (a_l.x * d_l.z + a_l.z * d_l.x + -2. * Z_l.x * a_l.z + -2. * Z_l.z * a_l.x + 2. * d_l.x * conn_ull.x.xz + 2. * d_l.y * conn_ull.y.xz + 2. * d_l.z * conn_ull.z.xz + 2. * K_ll.xz * tr_K + -2. * conn_ull.x.xx * conn_ull.x.xz + -2. * conn_ull.x.xy * conn_ull.y.xz + -2. * conn_ull.x.xz * conn_ull.z.xz + -2. * conn_ull.x.yz * conn_ull.y.xx + -2. * conn_ull.x.zz * conn_ull.z.xx + -2. * conn_ull.y.xy * conn_ull.y.yz + -2. * conn_ull.y.xz * conn_ull.z.yz + -2. * conn_ull.y.zz * conn_ull.z.xy + -2. * conn_ull.z.xz * conn_ull.z.zz + -4. * Z_l.x * conn_ull.x.xz + -4. * Z_l.y * conn_ull.y.xz + -4. * Z_l.z * conn_ull.z.xz + -4. * K_ul.x.x * K_ll.xz + -4. * K_ul.x.y * K_ll.yz + -4. * K_ul.x.z * K_ll.zz + -4. * K_ll.xz * Theta + -16. * S_ll.xz * M_PI + -8. * gamma_ll.xz * tmp2 + 8. * gamma_ll.xz * tmp1)) / 2.;
		(deriv)->K_ll.yy += alpha * (-conn_ull.x.xy * conn_ull.x.xy + -conn_ull.y.yy * conn_ull.y.yy + -conn_ull.z.yz * conn_ull.z.yz + a_l.y * d_l.y + d_l.x * conn_ull.x.yy + d_l.y * conn_ull.y.yy + d_l.z * conn_ull.z.yy + K_ll.yy * tr_K + -2. * Z_l.x * conn_ull.x.yy + -2. * Z_l.y * a_l.y + -2. * Z_l.y * conn_ull.y.yy + -2. * Z_l.z * conn_ull.z.yy + -2. * K_ul.x.y * K_ll.xy + -2. * K_ul.y.y * K_ll.yy + -2. * K_ul.y.z * K_ll.yz + -2. * K_ll.yy * Theta + -2. * conn_ull.x.yy * conn_ull.y.xy + -2. * conn_ull.x.yz * conn_ull.z.xy + -2. * conn_ull.y.yz * conn_ull.z.yy + -8. * S_ll.yy * M_PI + -4. * gamma_ll.yy * tmp2 + 4. * gamma_ll.yy * tmp1);
		(deriv)->K_ll.yz += (alpha * (a_l.y * d_l.z + a_l.z * d_l.y + -2. * Z_l.y * a_l.z + -2. * Z_l.z * a_l.y + 2. * d_l.x * conn_ull.x.yz + 2. * d_l.y * conn_ull.y.yz + 2. * d_l.z * conn_ull.z.yz + 2. * K_ll.yz * tr_K + -2. * conn_ull.x.xy * conn_ull.x.xz + -2. * conn_ull.x.yy * conn_ull.y.xz + -2. * conn_ull.x.yz * conn_ull.y.xy + -2. * conn_ull.x.yz * conn_ull.z.xz + -2. * conn_ull.x.zz * conn_ull.z.xy + -2. * conn_ull.y.yy * conn_ull.y.yz + -2. * conn_ull.y.yz * conn_ull.z.yz + -2. * conn_ull.y.zz * conn_ull.z.yy + -2. * conn_ull.z.yz * conn_ull.z.zz + -4. * Z_l.x * conn_ull.x.yz + -4. * Z_l.y * conn_ull.y.yz + -4. * Z_l.z * conn_ull.z.yz + -4. * K_ul.x.y * K_ll.xz + -4. * K_ul.y.y * K_ll.yz + -4. * K_ul.y.z * K_ll.zz + -4. * K_ll.yz * Theta + -16. * S_ll.yz * M_PI + -8. * gamma_ll.yz * tmp2 + 8. * gamma_ll.yz * tmp1)) / 2.;
		(deriv)->K_ll.zz += alpha * (-conn_ull.x.xz * conn_ull.x.xz + -conn_ull.y.yz * conn_ull.y.yz + -conn_ull.z.zz * conn_ull.z.zz + a_l.z * d_l.z + d_l.x * conn_ull.x.zz + d_l.y * conn_ull.y.zz + d_l.z * conn_ull.z.zz + K_ll.zz * tr_K + -2. * Z_l.x * conn_ull.x.zz + -2. * Z_l.y * conn_ull.y.zz + -2. * Z_l.z * a_l.z + -2. * Z_l.z * conn_ull.z.zz + -2. * K_ul.x.z * K_ll.xz + -2. * K_ul.y.z * K_ll.yz + -2. * K_ul.z.z * K_ll.zz + -2. * K_ll.zz * Theta + -2. * conn_ull.x.yz * conn_ull.y.xz + -2. * conn_ull.x.zz * conn_ull.z.xz + -2. * conn_ull.y.zz * conn_ull.z.yz + -8. * S_ll.zz * M_PI + -4. * gamma_ll.zz * tmp2 + 4. * gamma_ll.zz * tmp1);
		(deriv)->Theta += (alpha * (-K_ul.x.x * K_ul.x.x + -K_ul.y.y * K_ul.y.y + -K_ul.z.z * K_ul.z.z + tr_K * tr_K + -d_u.x * d_l.x + -d_u.y * d_l.y + -d_u.z * d_l.z + d_uuu.x.xx * d_lll.x.xx + -d_uuu.x.yy * d_lll.x.yy + -d_uuu.x.zz * d_lll.x.zz + -d_uuu.y.xx * d_lll.y.xx + d_uuu.y.yy * d_lll.y.yy + -d_uuu.y.zz * d_lll.y.zz + -d_uuu.z.xx * d_lll.z.xx + -d_uuu.z.yy * d_lll.z.yy + d_uuu.z.zz * d_lll.z.zz + 2. * Z_l.x * d_u.x + 2. * Z_l.y * d_u.y + 2. * Z_l.z * d_u.z + 2. * a_l.x * d_u.x + -2. * a_l.x * e_u.x + 2. * a_l.y * d_u.y + -2. * a_l.y * e_u.y + 2. * a_l.z * d_u.z + -2. * a_l.z * e_u.z + -2. * K_ul.x.y * K_ul.x.y + -2. * K_ul.x.z * K_ul.x.z + -2. * K_ul.y.z * K_ul.y.z + 2. * d_uuu.x.xy * d_lll.y.xx + 2. * d_uuu.x.xz * d_lll.z.xx + 2. * d_uuu.x.yy * d_lll.y.xy + -2. * d_uuu.x.yz * d_lll.x.yz + 2. * d_uuu.x.yz * d_lll.y.xz + 2. * d_uuu.x.yz * d_lll.z.xy + 2. * d_uuu.x.zz * d_lll.z.xz + 2. * d_uuu.y.xx * d_lll.x.xy + 2. * d_uuu.y.xy * d_lll.x.yy + 2. * d_uuu.y.xz * d_lll.x.yz + -2. * d_uuu.y.xz * d_lll.y.xz + 2. * d_uuu.y.xz * d_lll.z.xy + 2. * d_uuu.y.yz * d_lll.z.yy + 2. * d_uuu.y.zz * d_lll.z.yz + 2. * d_uuu.z.xx * d_lll.x.xz + 2. * d_uuu.z.xy * d_lll.x.yz + 2. * d_uuu.z.xy * d_lll.y.xz + -2. * d_uuu.z.xy * d_lll.z.xy + 2. * d_uuu.z.xz * d_lll.x.zz + 2. * d_uuu.z.yy * d_lll.y.yz + 2. * d_uuu.z.yz * d_lll.y.zz + -2. * Theta * tr_K + -4. * Z_l.x * a_u.x + -4. * Z_l.y * a_u.y + -16. * tmp2 + -4. * Z_l.z * a_u.z)) / 2.;
		(deriv)->Z_l.x += alpha * (-a_u.x * K_ll.xx + -a_u.y * K_ll.xy + -a_u.z * K_ll.xz + a_l.x * tr_K + d_u.x * K_ll.xx + d_u.y * K_ll.xy + d_u.z * K_ll.xz + -K_uu.xx * d_lll.x.xx + -K_uu.yy * d_lll.x.yy + -K_uu.zz * d_lll.x.zz + -2. * Z_u.x * K_ll.xx + -2. * Z_u.y * K_ll.xy + -2. * Z_u.z * K_ll.xz + -2. * a_l.x * Theta + -2. * K_uu.xy * d_lll.x.xy + -2. * K_uu.xz * d_lll.x.xz + -8. * S_l.x * M_PI + -2. * K_uu.yz * d_lll.x.yz);
		(deriv)->Z_l.y += alpha * (-a_u.x * K_ll.xy + -a_u.y * K_ll.yy + -a_u.z * K_ll.yz + a_l.y * tr_K + d_u.x * K_ll.xy + d_u.y * K_ll.yy + d_u.z * K_ll.yz + -K_uu.xx * d_lll.y.xx + -K_uu.yy * d_lll.y.yy + -K_uu.zz * d_lll.y.zz + -2. * Z_u.x * K_ll.xy + -2. * Z_u.y * K_ll.yy + -2. * Z_u.z * K_ll.yz + -2. * a_l.y * Theta + -2. * K_uu.xy * d_lll.y.xy + -2. * K_uu.xz * d_lll.y.xz + -8. * S_l.y * M_PI + -2. * K_uu.yz * d_lll.y.yz);
		(deriv)->Z_l.z += alpha * (-a_u.x * K_ll.xz + -a_u.y * K_ll.yz + -a_u.z * K_ll.zz + a_l.z * tr_K + d_u.x * K_ll.xz + d_u.y * K_ll.yz + d_u.z * K_ll.zz + -K_uu.xx * d_lll.z.xx + -K_uu.yy * d_lll.z.yy + -K_uu.zz * d_lll.z.zz + -2. * Z_u.x * K_ll.xz + -2. * Z_u.y * K_ll.yz + -2. * Z_u.z * K_ll.zz + -2. * a_l.z * Theta + -2. * K_uu.xy * d_lll.z.xy + -2. * K_uu.xz * d_lll.z.xz + -8. * S_l.z * M_PI + -2. * K_uu.yz * d_lll.z.yz);
	}
	<? if eqn.useShift ~= "none" then ?>
	{
		real const tmp1 = alpha * alpha;
		real const tmp2 = a_u.x * tmp1;
		real const tmp3 = d_u.x * tmp1;
		real const tmp4 = e_u.x * tmp1;
		real const tmp5 = 2. * tmp4;
		real const tmp6 = -tmp3;
		real const tmp7 = -tmp2;
		real const tmp8 = a_u.y * tmp1;
		real const tmp9 = d_u.y * tmp1;
		real const tmp10 = e_u.y * tmp1;
		real const tmp11 = 2. * tmp10;
		real const tmp12 = -tmp9;
		real const tmp13 = -tmp8;
		real const tmp14 = a_u.z * tmp1;
		real const tmp15 = d_u.z * tmp1;
		real const tmp16 = e_u.z * tmp1;
		real const tmp17 = 2. * tmp16;
		real const tmp18 = -tmp15;
		real const tmp19 = -tmp14;
		(deriv)->alpha += alpha * (beta_u.x * a_l.x + beta_u.z * a_l.z + beta_u.y * a_l.y);
		(deriv)->gammaDelta_ll.xx += 2. * (b_ll.xx + beta_u.x * d_lll.x.xx + beta_u.z * d_lll.z.xx + beta_u.y * d_lll.y.xx);
		(deriv)->gammaDelta_ll.xy += b_ll.xy + b_ll.xy + 2. * beta_u.x * d_lll.x.xy + 2. * beta_u.z * d_lll.z.xy + 2. * beta_u.y * d_lll.y.xy;
		(deriv)->gammaDelta_ll.xz += b_ll.xz + b_ll.xz + 2. * beta_u.x * d_lll.x.xz + 2. * beta_u.z * d_lll.z.xz + 2. * beta_u.y * d_lll.y.xz;
		(deriv)->gammaDelta_ll.yy += 2. * (b_ll.yy + beta_u.x * d_lll.x.yy + beta_u.z * d_lll.z.yy + beta_u.y * d_lll.y.yy);
		(deriv)->gammaDelta_ll.yz += b_ll.yz + b_ll.yz + 2. * beta_u.x * d_lll.x.yz + 2. * beta_u.z * d_lll.z.yz + 2. * beta_u.y * d_lll.y.yz;
		(deriv)->gammaDelta_ll.zz += 2. * (b_ll.zz + beta_u.x * d_lll.x.zz + beta_u.z * d_lll.z.zz + beta_u.y * d_lll.y.zz);
		(deriv)->a_l.x += a_l.x * b_ul.x.x + -a_l.x * tr_b + a_l.z * b_ul.z.x + a_l.y * b_ul.y.x;
		(deriv)->a_l.y += a_l.x * b_ul.x.y + a_l.y * b_ul.y.y + a_l.z * b_ul.z.y + -a_l.y * tr_b;
		(deriv)->a_l.z += a_l.x * b_ul.x.z + a_l.y * b_ul.y.z + -a_l.z * tr_b + a_l.z * b_ul.z.z;
		(deriv)->dDelta_lll.x.xx += b_ul.x.x * dDelta_lll.x.xx + b_ul.y.x * dDelta_lll.y.xx + -dDelta_lll.x.xx * tr_b + b_ul.z.x * dDelta_lll.z.xx;
		(deriv)->dDelta_lll.x.xy += b_ul.x.x * dDelta_lll.x.xy + b_ul.y.x * dDelta_lll.y.xy + -dDelta_lll.x.xy * tr_b + b_ul.z.x * dDelta_lll.z.xy;
		(deriv)->dDelta_lll.x.xz += b_ul.x.x * dDelta_lll.x.xz + b_ul.y.x * dDelta_lll.y.xz + -dDelta_lll.x.xz * tr_b + b_ul.z.x * dDelta_lll.z.xz;
		(deriv)->dDelta_lll.x.yy += b_ul.x.x * dDelta_lll.x.yy + b_ul.y.x * dDelta_lll.y.yy + -dDelta_lll.x.yy * tr_b + b_ul.z.x * dDelta_lll.z.yy;
		(deriv)->dDelta_lll.x.yz += b_ul.x.x * dDelta_lll.x.yz + b_ul.y.x * dDelta_lll.y.yz + -dDelta_lll.x.yz * tr_b + b_ul.z.x * dDelta_lll.z.yz;
		(deriv)->dDelta_lll.x.zz += b_ul.x.x * dDelta_lll.x.zz + b_ul.y.x * dDelta_lll.y.zz + -dDelta_lll.x.zz * tr_b + b_ul.z.x * dDelta_lll.z.zz;
		(deriv)->dDelta_lll.y.xx += b_ul.x.y * dDelta_lll.x.xx + b_ul.y.y * dDelta_lll.y.xx + -dDelta_lll.y.xx * tr_b + b_ul.z.y * dDelta_lll.z.xx;
		(deriv)->dDelta_lll.y.xy += b_ul.x.y * dDelta_lll.x.xy + b_ul.y.y * dDelta_lll.y.xy + -dDelta_lll.y.xy * tr_b + b_ul.z.y * dDelta_lll.z.xy;
		(deriv)->dDelta_lll.y.xz += b_ul.x.y * dDelta_lll.x.xz + b_ul.y.y * dDelta_lll.y.xz + -dDelta_lll.y.xz * tr_b + b_ul.z.y * dDelta_lll.z.xz;
		(deriv)->dDelta_lll.y.yy += b_ul.x.y * dDelta_lll.x.yy + b_ul.y.y * dDelta_lll.y.yy + -dDelta_lll.y.yy * tr_b + b_ul.z.y * dDelta_lll.z.yy;
		(deriv)->dDelta_lll.y.yz += b_ul.x.y * dDelta_lll.x.yz + b_ul.y.y * dDelta_lll.y.yz + -dDelta_lll.y.yz * tr_b + b_ul.z.y * dDelta_lll.z.yz;
		(deriv)->dDelta_lll.y.zz += b_ul.x.y * dDelta_lll.x.zz + b_ul.y.y * dDelta_lll.y.zz + -dDelta_lll.y.zz * tr_b + b_ul.z.y * dDelta_lll.z.zz;
		(deriv)->dDelta_lll.z.xx += b_ul.x.z * dDelta_lll.x.xx + b_ul.y.z * dDelta_lll.y.xx + -dDelta_lll.z.xx * tr_b + b_ul.z.z * dDelta_lll.z.xx;
		(deriv)->dDelta_lll.z.xy += b_ul.x.z * dDelta_lll.x.xy + b_ul.y.z * dDelta_lll.y.xy + -dDelta_lll.z.xy * tr_b + b_ul.z.z * dDelta_lll.z.xy;
		(deriv)->dDelta_lll.z.xz += b_ul.x.z * dDelta_lll.x.xz + b_ul.y.z * dDelta_lll.y.xz + -dDelta_lll.z.xz * tr_b + b_ul.z.z * dDelta_lll.z.xz;
		(deriv)->dDelta_lll.z.yy += b_ul.x.z * dDelta_lll.x.yy + b_ul.y.z * dDelta_lll.y.yy + -dDelta_lll.z.yy * tr_b + b_ul.z.z * dDelta_lll.z.yy;
		(deriv)->dDelta_lll.z.yz += b_ul.x.z * dDelta_lll.x.yz + b_ul.y.z * dDelta_lll.y.yz + -dDelta_lll.z.yz * tr_b + b_ul.z.z * dDelta_lll.z.yz;
		(deriv)->dDelta_lll.z.zz += b_ul.x.z * dDelta_lll.x.zz + b_ul.y.z * dDelta_lll.y.zz + -dDelta_lll.z.zz * tr_b + b_ul.z.z * dDelta_lll.z.zz;
		(deriv)->K_ll.xx += -K_ll.xx * tr_b + 2. * K_ll.xx * b_ul.x.x + 2. * K_ll.xz * b_ul.z.x + 2. * K_ll.xy * b_ul.y.x;
		(deriv)->K_ll.xy += K_ll.xx * b_ul.x.y + K_ll.xy * b_ul.x.x + K_ll.xy * b_ul.y.y + -K_ll.xy * tr_b + K_ll.xz * b_ul.z.y + K_ll.yz * b_ul.z.x + K_ll.yy * b_ul.y.x;
		(deriv)->K_ll.xz += K_ll.xx * b_ul.x.z + K_ll.xy * b_ul.y.z + K_ll.xz * b_ul.x.x + K_ll.xz * b_ul.z.z + -K_ll.xz * tr_b + K_ll.zz * b_ul.z.x + K_ll.yz * b_ul.y.x;
		(deriv)->K_ll.yy += -K_ll.yy * tr_b + 2. * K_ll.xy * b_ul.x.y + 2. * K_ll.yz * b_ul.z.y + 2. * K_ll.yy * b_ul.y.y;
		(deriv)->K_ll.yz += K_ll.xy * b_ul.x.z + K_ll.xz * b_ul.x.y + K_ll.yy * b_ul.y.z + K_ll.yz * b_ul.y.y + K_ll.yz * b_ul.z.z + K_ll.zz * b_ul.z.y + -K_ll.yz * tr_b;
		(deriv)->K_ll.zz += -K_ll.zz * tr_b + 2. * K_ll.xz * b_ul.x.z + 2. * K_ll.zz * b_ul.z.z + 2. * K_ll.yz * b_ul.y.z;
		(deriv)->Theta += -Theta * tr_b;
		(deriv)->Z_l.x += Z_l.x * b_ul.x.x + -Z_l.x * tr_b + Z_l.z * b_ul.z.x + Z_l.y * b_ul.y.x;
		(deriv)->Z_l.y += Z_l.x * b_ul.x.y + Z_l.y * b_ul.y.y + Z_l.z * b_ul.z.y + -Z_l.y * tr_b;
		(deriv)->Z_l.z += Z_l.x * b_ul.x.z + Z_l.y * b_ul.y.z + -Z_l.z * tr_b + Z_l.z * b_ul.z.z;
		<? if eqn.useShift == "HarmonicShiftParabolic" then ?>\
		(deriv)->beta_u.x += beta_u.x * b_ul.x.x + beta_u.y * b_ul.x.y + beta_u.z * b_ul.x.z + tmp5 + tmp6 + tmp7;
		(deriv)->beta_u.y += beta_u.x * b_ul.y.x + beta_u.y * b_ul.y.y + beta_u.z * b_ul.y.z + tmp11 + tmp12 + tmp13;
		(deriv)->beta_u.z += beta_u.x * b_ul.z.x + beta_u.y * b_ul.z.y + beta_u.z * b_ul.z.z + tmp17 + tmp18 + tmp19;
		<? elseif eqn.useShift == "HarmonicShiftHyperbolic" then ?>\
		(deriv)->beta_u.x += B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y;
		(deriv)->beta_u.y += B_u.y + beta_u.x * b_ul.y.x + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y;
		(deriv)->beta_u.z += B_u.z + beta_u.x * b_ul.z.x + beta_u.z * b_ul.z.z + beta_u.y * b_ul.z.y;
		(deriv)->B_u.x += B_u.x * b_ul.x.x + B_u.y * b_ul.x.y + B_u.z * b_ul.x.z + tmp5 + tmp6 + tmp7;
		(deriv)->B_u.y += B_u.x * b_ul.y.x + B_u.y * b_ul.y.y + B_u.z * b_ul.y.z + tmp11 + tmp12 + tmp13;
		(deriv)->B_u.z += B_u.x * b_ul.z.x + B_u.y * b_ul.z.y + B_u.z * b_ul.z.z + tmp17 + tmp18 + tmp19;
		<? end ?>\
	}
	<? end ?>/* eqn.useShift ~= "none" */	
	// END CUT
<? end ?>
//decay for 1st deriv hyperbolic state vars constraints:

	//turns out if you "if conv != 0" all these then skipping decay explodes soon in simulation steps, but time grows quickly, so it dies at a high t value

//// MODULE_DEPENDS: <?=calcFromGrad_a_l?>
	// a_x = log(α)_,x <=> a_x += η (log(α)_,x - a_x)
	if (solver->a_convCoeff != 0.) {
		real3 const target_a_l = <?=calcFromGrad_a_l?>(solver, U);
<?
for i,xi in ipairs(xNames) do
?>		deriv->a_l.<?=xi?> += solver->a_convCoeff * (target_a_l.<?=xi?> - U->a_l.<?=xi?>);
<? end
?>	}

//// MODULE_DEPENDS: <?=calcFromGrad_d_lll?>
	// d_xxx = .5 γ_xx,x <=> d_xxx += η (.5 γ_xx,x - d_xxx)
	if (solver->d_convCoeff != 0.) {
		_3sym3 const target_d_lll = <?=calcFromGrad_d_lll?>(solver, U, cell);
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
?>		deriv->dDelta_lll.<?=xi?>.<?=xjk?> += solver->d_convCoeff * (target_d_lll.<?=xi?>.<?=xjk?> - d_lll.<?=xi?>.<?=xjk?>);
<? 	end
end
?>	}

<? if has_b_ul then ?>
//// MODULE_DEPENDS: <?=calcFromGrad_b_ul?>
	// b_ul.j.i = beta_u.j,i <=> b_ul.j.i += η (β_u.j,i - b_ul.j.i)
	if (solver->b_convCoeff != 0.) {
		real3x3 const target_b_ul = <?=calcFromGrad_b_ul?>(solver, U);
<?
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>		deriv->b_ul.<?=xi?>.<?=xj?> += solver->b_convCoeff * (target_b_ul.<?=xj?>.<?=xi?> - U->b_ul.<?=xj?>.<?=xi?>);
<?	end
end
?>	}
<? end ?>

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

//// MODULE_DEPENDS: <?=calc_gamma_uu?>
	sym3 const gamma_uu = <?=calc_gamma_uu?>(U, cell->pos);

	real3x3 const K_ul = sym3_sym3_mul(gamma_uu, U->K_ll);			//K^i_j
	real const tr_K = real3x3_trace(K_ul);							//K^k_k
	sym3 const K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);		//K^ij

//// MODULE_DEPENDS: <?=calc_d_lll?>
	_3sym3 const d_lll = <?=calc_d_lll?>(U, cell->pos);	
	real3x3x3 const d_llu = _3sym3_sym3_mul(d_lll, gamma_uu);	//d_llu = d_ij^k = d_ijl * γ^lk	
	_3sym3 const d_ull = sym3_3sym3_mul(gamma_uu, d_lll);	//d_ull = d^i_jk = γ^il d_ljk
	_3sym3 const conn_ull = conn_ull_from_d_llu_d_ull(d_llu, d_ull);	//Γ^k_ij = d_ij^k + d_ji^k - d^k_ij
	real3 const e_l = _3sym3_tr12(d_ull);	//e_i = d^j_ji
	real3 const d_l = real3x3x3_tr23(d_llu);	//d_l.i = d_i = d_ij^j

	//partial_d_lll.ij.kl = d_kij,l = d_(k|(ij),|l)
	//so this object's indexes are rearranged compared to the papers
	//sym3sym3 partial_d_llll = sym3sym3_zero; //sym3sym3_zero doesn't exist
	sym3sym3 partial_d_llll = {
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = sym3_zero,
<? end
?>	};

	<?
for k=1,solver.dim do	-- beyond dim and the finite-difference will be zero
	local xk = xNames[k]
	?>{
		global <?=cons_t?> const * const UR = U + solver->stepsize.<?=xk?>;
		global <?=cons_t?> const * const UL = U - solver->stepsize.<?=xk?>;
		global <?=cell_t?> const * const cellR = cell + solver->stepsize.<?=xk?>;
		global <?=cell_t?> const * const cellL = cell - solver->stepsize.<?=xk?>;
		_3sym3 const dR_lll = <?=calc_d_lll?>(UR, cellR->pos);
		_3sym3 const dL_lll = <?=calc_d_lll?>(UL, cellL->pos);
<?	for l=k,3 do	-- since we are writing to xl, only iterate through symmetric terms
		local xl = xNames[l]
		for ij,xij in ipairs(symNames) do
?>		partial_d_llll.<?=xij?>.<?=sym(k,l)?> += (dR_lll.<?=xk?>.<?=xij?> - dL_lll.<?=xk?>.<?=xij?>) / (2. * solver->grid_dx.<?=xk?>);
<?		end
	end
?>	}<?
end ?>

	// R_ll.ij := R_ij
	//	= γ^kl (-γ_ij,kl - γ_kl,ij + γ_ik,jl + γ_jl,ik)
	//		+ Γ^k_ij (d_k - 2 e_k)
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
	//H = 1/2 (R + K^2 - K_ij K^ij) - 8 π ρ
	real const R = sym3_dot(R_ll, gamma_uu);
	real const tr_KSq = sym3_dot(U->K_ll, K_uu);
	U->H = .5 * (R + tr_K * tr_K - tr_KSq) <?
if eqn.useStressEnergyTerms then ?>
	- 8. * M_PI * U->rho <?
end ?>;

<?=eqn:makePartial1"K_ll"?>

	/*
	momentum constraint
	Alcubierre eqn 2.4.11
	M^i = K^ij_;j
		- γ^ij K_,j
		- 8 π S^i
	M^i = γ^im γ^jn (
			K_mn,j
			- Γ^k_mj K_kn
			- Γ^k_nj K_mk
		)
		- γ^ij (γ^mn_,j K_mn + γ^mn K_mn,j)
		- 8 π S^i
	M^i = γ^ij (
			γ^mn (
				K_jn,m
				- K_mn,j
				- Γ^k_jm K_kn
				- Γ^k_nm K_jk
			)
			+ 2 d_jmn K^mn
		)
		- 8 π S^i
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
