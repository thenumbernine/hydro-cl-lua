<?
local useConstrainU = true -- constrains alpha to alphamin and calcs H and M^i
local useAddSource = true
local useKreissOligarDissipation = true	-- depends on useAddSource

local has_beta_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "beta_u" end)
local has_betaLap_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "betaLap_u" end)
local has_b_ul = eqn.consStruct.vars:find(nil, function(var) return var.name == "b_ul" end)
local has_B_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end)
?>

//// MODULE_NAME: mdeShiftEpsilon

<? if eqn.useShift == "MinimalDistortionHyperbolic" then ?>
constant real const mdeShiftEpsilon = 1.;
<? end ?>


<? -- ok so ideally you should be able to subtract out the background metric ...
?>

//// MODULE_NAME: <?=calc_gammaHat_ll?>
//// MODULE_DEPENDS: sym3

static inline sym3 <?=calc_gammaHat_ll?>(real3 const pt) {
<? if true then -- get background grid separation working
?>
//// MODULE_DEPENDS: <?=coord_gHol_ll?>
	return coord_gHol_ll(pt);
<? else ?>
	return sym3_ident;
<? end ?>
}

//// MODULE_NAME: <?=calc_dHat_lll?>
//// MODULE_DEPENDS: _3sym3

static inline _3sym3 <?=calc_dHat_lll?>(real3 const pt) {
<? if true then 	-- get background grid separation working?
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

<? if has_betaLap_u then
?>	U->betaLap_u = real3_zero;
<? end
if has_beta_u then
?>	U->beta_u = real3_rescaleToCoord_U(beta_U, x);
<? end
if has_B_u then
?>	U->B_u = real3_rescaleToCoord_U(B_U, x);
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
		
<? if eqn.useShift == "MinimalDistortionHyperbolic" then ?>
//// MODULE_DEPENDS: mdeShiftEpsilon
<? end ?>

<? if eqn.useShift == "GammaDriverHyperbolic" then ?>
//// MODULE_DEPENDS: <?=calc_gammaHat_ll?>
//// MODULE_DEPENDS: <?=coord_conn_ull?>
<? end ?>

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
if has_B_u then --\
?>	real3 B_u = (U)->B_u;\
<? 	end --\
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
end --\
if has_B_u then --\
?>	B_u = real3_rotateFrom(B_u, n_l);\
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
\
	<? if eqn.useShift ~= "none" then ?>\
	sym3 const b_ll = sym3_real3x3_to_sym3_mul(gamma_ll, b_ul);\
	<? end ?>\
\
	<? if eqn.useShift == "HarmonicParabolic" then ?>\
	real3 const a_u = sym3_real3_mul(gamma_uu, a_l);\
	<? end ?>\
\
	<? if eqn.useShift == "MinimalDistortionHyperbolic" then ?>\
	real3x3 const DBeta_ul = real3x3_add(b_ul, real3_3sym3_dot2(beta_u, conn_ull));\
	sym3 const DBeta_uu = real3x3_sym3_to_sym3_mul(DBeta_ul, gamma_uu);\
	real const tr_DBeta = real3x3_trace(DBeta_ul);\
	sym3 const K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);\
	sym3 const A_uu = sym3_sub(K_uu, sym3_real_mul(gamma_uu, tr_K / 3.));\
	<? end ?>\
\
	<? if eqn.useShift == "GammaDriverHyperbolic" then ?>\
	sym3 const gammaHat_ll = <?=calc_gammaHat_ll?>((cell)->pos);\
	real const det_gammaHat = sym3_det(gammaHat_ll);\
	real const W = pow(det_gammaHat / det_gamma, 1./6.);\
	real const invW = 1. / W;\
	_3sym3 const connHat_ull = coord_conn_ull((cell)->pos);\
	real3x3 const DHatBeta_ul = real3x3_add(b_ul, real3_3sym3_dot2(beta_u, connHat_ull));\
	real const tr_DHatBeta = real3x3_trace(DHatBeta_ul);\
	<? end ?>\
\
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
		{\
			(resultFlux)->a_l.x += -beta_u.x * a_l.x;\
			(resultFlux)->a_l.y += -beta_u.x * a_l.y;\
			(resultFlux)->a_l.z += -beta_u.x * a_l.z;\
			(resultFlux)->dDelta_lll.x.xx += -(b_ll.xx + beta_u.x * dHat_lll.x.xx + beta_u.x * dDelta_lll.x.xx + beta_u.z * dHat_lll.z.xx + beta_u.y * dHat_lll.y.xx);\
			(resultFlux)->dDelta_lll.x.xy += -(b_ll.xy + beta_u.x * dHat_lll.x.xy + beta_u.x * dDelta_lll.x.xy + beta_u.z * dHat_lll.z.xy + beta_u.y * dHat_lll.y.xy);\
			(resultFlux)->dDelta_lll.x.xz += -(b_ll.xz + beta_u.x * dHat_lll.x.xz + beta_u.x * dDelta_lll.x.xz + beta_u.z * dHat_lll.z.xz + beta_u.y * dHat_lll.y.xz);\
			(resultFlux)->dDelta_lll.x.yy += -(b_ll.yy + beta_u.x * dHat_lll.x.yy + beta_u.x * dDelta_lll.x.yy + beta_u.z * dHat_lll.z.yy + beta_u.y * dHat_lll.y.yy);\
			(resultFlux)->dDelta_lll.x.yz += -(b_ll.yz + beta_u.x * dHat_lll.x.yz + beta_u.x * dDelta_lll.x.yz + beta_u.z * dHat_lll.z.yz + beta_u.y * dHat_lll.y.yz);\
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
		}\
		<? if eqn.useShift == "HarmonicParabolic" then ?>\
		{\
			real const tmp1 = alpha * alpha;\
			(resultFlux)->b_ul.x.x += -beta_u.x * b_ul.x.x + -beta_u.y * b_ul.x.y + -beta_u.z * b_ul.x.z + a_u.x * tmp1 + -2. * e_u.x * tmp1 + d_u.x * tmp1;\
			(resultFlux)->b_ul.y.x += -beta_u.x * b_ul.x.y + -beta_u.y * b_ul.y.y + -beta_u.z * b_ul.y.z + a_u.y * tmp1 + -2. * e_u.y * tmp1 + d_u.y * tmp1;\
			(resultFlux)->b_ul.z.x += -beta_u.x * b_ul.x.z + -beta_u.y * b_ul.y.z + -beta_u.z * b_ul.z.z + a_u.z * tmp1 + -2. * e_u.z * tmp1 + d_u.z * tmp1;\
		}\
		<? end ?>/* eqn.useShift == "HarmonicParabolic" */\
		<? if eqn.useShift == "HarmonicHyperbolic" then ?>\
		{\
			(resultFlux)->b_ul.x.x += -(B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y);\
			(resultFlux)->b_ul.y.x += -(B_u.y + beta_u.x * b_ul.x.y + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y);\
			(resultFlux)->b_ul.z.x += -(B_u.z + beta_u.x * b_ul.x.z + beta_u.z * b_ul.z.z + beta_u.y * b_ul.y.z);\
		}\
		<? end ?>/* eqn.useShift == "HarmonicHyperbolic" */\
		<? if eqn.useShift == "MinimalDistortionHyperbolic" then ?>\
		{\
			real const tmp1 = d_l.y * gamma_uu.xy;\
			real const tmp2 = d_l.z * gamma_uu.xz;\
			real const tmp3 = gamma_uu.xy * conn_ull.x.xy;\
			real const tmp4 = gamma_uu.xz * conn_ull.x.xz;\
			real const tmp5 = gamma_uu.yz * conn_ull.x.yz;\
			(resultFlux)->b_ul.x.x += -(B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y);\
			(resultFlux)->b_ul.y.x += -(B_u.y + beta_u.x * b_ul.x.y + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y);\
			(resultFlux)->b_ul.z.x += -(B_u.z + beta_u.x * b_ul.x.z + beta_u.z * b_ul.z.z + beta_u.y * b_ul.y.z);\
			(resultFlux)->B_u.x += (mdeShiftEpsilon * (-6. * DBeta_uu.xx + -2. * gamma_uu.xx * tr_DBeta + 3. * beta_u.x * tmp1 + 3. * beta_u.x * tmp2 + 3. * beta_u.y * d_l.y * gamma_uu.xx + 3. * beta_u.z * d_l.z * gamma_uu.xx + 12. * A_uu.xx * alpha + 6. * beta_u.x * d_l.x * gamma_uu.xx + -6. * beta_u.x * gamma_uu.xx * conn_ull.x.xx + -6. * beta_u.x * tmp3 + -6. * beta_u.x * tmp4 + -6. * beta_u.y * gamma_uu.xx * conn_ull.x.xy + -6. * beta_u.y * gamma_uu.xy * conn_ull.x.yy + -6. * beta_u.y * gamma_uu.xz * conn_ull.x.yz + -6. * beta_u.z * gamma_uu.xx * conn_ull.x.xz + -6. * beta_u.z * gamma_uu.xz * conn_ull.x.zz + -6. * beta_u.z * gamma_uu.xy * conn_ull.x.yz)) / 6.;\
			(resultFlux)->B_u.y += (mdeShiftEpsilon * (-6. * DBeta_uu.xy + -2. * gamma_uu.xy * tr_DBeta + 3. * beta_u.x * d_l.y * gamma_uu.yy + 3. * beta_u.x * d_l.z * gamma_uu.yz + 3. * beta_u.y * tmp1 + 3. * beta_u.z * d_l.z * gamma_uu.xy + 12. * A_uu.xy * alpha + 6. * beta_u.x * d_l.x * gamma_uu.xy + -6. * beta_u.x * gamma_uu.xy * conn_ull.x.xx + -6. * beta_u.x * gamma_uu.yy * conn_ull.x.xy + -6. * beta_u.x * gamma_uu.yz * conn_ull.x.xz + -6. * beta_u.y * tmp3 + -6. * beta_u.y * gamma_uu.yy * conn_ull.x.yy + -6. * beta_u.y * tmp5 + -6. * beta_u.z * gamma_uu.xy * conn_ull.x.xz + -6. * beta_u.z * gamma_uu.yz * conn_ull.x.zz + -6. * beta_u.z * gamma_uu.yy * conn_ull.x.yz)) / 6.;\
			(resultFlux)->B_u.z += (mdeShiftEpsilon * (-6. * DBeta_uu.xz + -2. * gamma_uu.xz * tr_DBeta + 3. * beta_u.x * d_l.y * gamma_uu.yz + 3. * beta_u.x * d_l.z * gamma_uu.zz + 3. * beta_u.y * d_l.y * gamma_uu.xz + 3. * beta_u.z * tmp2 + 12. * A_uu.xz * alpha + 6. * beta_u.x * d_l.x * gamma_uu.xz + -6. * beta_u.x * gamma_uu.xz * conn_ull.x.xx + -6. * beta_u.x * gamma_uu.yz * conn_ull.x.xy + -6. * beta_u.x * gamma_uu.zz * conn_ull.x.xz + -6. * beta_u.y * gamma_uu.xz * conn_ull.x.xy + -6. * beta_u.y * gamma_uu.yz * conn_ull.x.yy + -6. * beta_u.y * gamma_uu.zz * conn_ull.x.yz + -6. * beta_u.z * tmp4 + -6. * beta_u.z * gamma_uu.zz * conn_ull.x.zz + -6. * beta_u.z * tmp5)) / 6.;\
		}\
		<? end ?>/* eqn.useShift == "MinimalDistortionHyperbolic" */\
		<? if eqn.useShift == "GammaDriverHyperbolic" then ?>\
		{\
			real const tmp1 = alpha * tr_K;\
			real const tmp2 = invW * invW;\
			(resultFlux)->b_ul.x.x += -(B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y);\
			(resultFlux)->b_ul.y.x += -(B_u.y + beta_u.x * b_ul.x.y + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y);\
			(resultFlux)->b_ul.z.x += -(B_u.z + beta_u.x * b_ul.x.z + beta_u.z * b_ul.z.z + beta_u.y * b_ul.y.z);\
			(resultFlux)->B_u.x += ((-gamma_uu.xx * tr_DHatBeta + -3. * gamma_uu.xx * DHatBeta_ul.x.x + -3. * gamma_uu.xy * DHatBeta_ul.x.y + 4. * gamma_uu.xx * tmp1 + -3. * gamma_uu.xz * DHatBeta_ul.x.z) * tmp2) / 4.;\
			(resultFlux)->B_u.y += ((-gamma_uu.xy * tr_DHatBeta + -3. * gamma_uu.xx * DHatBeta_ul.y.x + -3. * gamma_uu.xy * DHatBeta_ul.y.y + 4. * gamma_uu.xy * tmp1 + -3. * gamma_uu.xz * DHatBeta_ul.y.z) * tmp2) / 4.;\
			(resultFlux)->B_u.z += ((-gamma_uu.xz * tr_DHatBeta + -3. * gamma_uu.xx * DHatBeta_ul.z.x + -3. * gamma_uu.xy * DHatBeta_ul.z.y + 4. * gamma_uu.xz * tmp1 + -3. * gamma_uu.xz * DHatBeta_ul.z.z) * tmp2) / 4.;\
		}\
		<? end ?>/* eqn.useShift == "GammaDriverHyperbolic" */\
	}\
	<? end ?>/* eqn.useShift ~= "none" */\
	/* END CUT from symmath/tests/output/Z4.html */\
\
<? if solver.coord.vectorComponent == "cartesian" then ?>\
\
	(resultFlux)->Z_l = real3_rotateTo((resultFlux)->Z_l, n_l);\
	(resultFlux)->a_l = real3_rotateTo((resultFlux)->a_l, n_l);\
	(resultFlux)->gammaDelta_ll = sym3_rotateTo((resultFlux)->gammaDelta_ll, n_l);\
	(resultFlux)->K_ll = sym3_rotateTo((resultFlux)->K_ll, n_l);\
	(resultFlux)->dDelta_lll = _3sym3_rotateTo((resultFlux)->dDelta_lll, n_l);\
<? if has_beta_u then --\
?>	(resultFlux)->beta_u = real3_rotateTo((resultFlux)->beta_u, n_l);\
<? end --\
if has_b_ul then --\
?>	(resultFlux)->b_ul = real3x3_rotateTo((resultFlux)->b_ul, n_l);\
<? end --\
if has_B_u then --\
?>	(resultFlux)->B_u = real3_rotateTo((resultFlux)->B_u, n_l);\
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
if has_b_ul then --\
?>		(resultFlux)->b_ul = real3x3_swap<?=side?>((resultFlux)->b_ul);\
<? end --\
if has_B_u then --\
?>		(resultFlux)->B_u = real3_swap<?=side?>((resultFlux)->B_u);\
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
	(resultEig)->gamma_ll = sym3_real_mul(sym3_add(gammaL_ll, gammaR_ll), .5);\
	real const det_avg_gamma = sym3_det((resultEig)->gamma_ll);\
	(resultEig)->gamma_uu = sym3_inv((resultEig)->gamma_ll, det_avg_gamma);\
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
	/* rotate from flux direction to x-axis */\
	/* TODO copy this from fluxFromCons or make its own function */\
	real Theta = (inputU)->Theta;												/* 27 */\
\
	real3 a_l = (inputU)->a_l;\
	_3sym3 dDelta_lll = (inputU)->dDelta_lll;\
	_3sym3 dHat_lll = <?=calc_dHat_lll?>(pt);\
	sym3 K_ll = (inputU)->K_ll;\
	real3 Z_l = (inputU)->Z_l;\
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
	K_ll = sym3_rotateFrom(K_ll, n_l);\
	dDelta_lll = _3sym3_rotateFrom(dDelta_lll, n_l);\
	dHat_lll = _3sym3_rotateFrom(dHat_lll, n_l);\
<? if has_beta_u then --\
?>	beta_u = real3_rotateFrom(beta_u);\
<? end --\
if has_b_ul then --\
?>	b_ul = real3x3_rotateFrom(b_ul);\
end --\
if has_B_u then --\
?>	B_u = real3_rotateFrom(B_u, n_l);\
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
		K_ll = sym3_swap<?=side?>(K_ll);\
		dDelta_lll = _3sym3_swap<?=side?>(dDelta_lll);\
		dHat_lll = _3sym3_swap<?=side?>(dHat_lll);\
	}\
	<? end ?>\
	else {\
		Theta = 0./0.;\
<? for k,xk in ipairs(xNames) do --\
?>		a_l.<?=xk?> = 0./0.;\
		Z_l.<?=xk?> = 0./0.;\
<? end --\
?>\
<? for ij,xij in ipairs(symNames) do --\
?>		K_ll.<?=xij?> = 0./0.;\
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
	_3sym3 d_lll = _3sym3_add(dDelta_lll, dHat_lll);\
\
	/* the eigen_t vars should already be rotated into flux normal frame */\
	sym3 const gamma_ll = (eig)->gamma_ll;\
	sym3 const gamma_uu = (eig)->gamma_uu;\
\
	real const sqrt_gammaUUxx = sqrt(gamma_uu.xx);\
	real const gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;\
	real const gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;\
\
	real const sqrt_f = (eig)->alpha_sqrt_f / (eig)->alpha;\
	real const f = sqrt_f * sqrt_f;\
\
<? if eqn.noZeroRowsInFlux then ?>\
	real const tmp1 = sqrt_f * sqrt_f;\
	real const tmp2 = gamma_uu.xy * d_lll.x.yy;\
	real const tmp3 = gamma_uu.xz * d_lll.x.yz;\
	real const tmp4 = gamma_uu.xy * d_lll.x.yz;\
	real const tmp5 = gamma_uu.xz * d_lll.x.zz;\
	real const tmp6 = 1. / 2.;\
	real const tmp7 = sqrt(gamma_uu.xx);\
	real const tmp8 = K_ll.yy * tmp7;\
	real const tmp9 = gamma_uu.xx * d_lll.x.yy;\
	real const tmp10 = K_ll.yz * tmp7;\
	real const tmp11 = gamma_uu.xx * d_lll.x.yz;\
	real const tmp12 = K_ll.zz * tmp7;\
	real const tmp13 = gamma_uu.xx * d_lll.x.zz;\
	real const tmp14 = -1. * tmp3;\
	real const tmp15 = K_ll.xy * tmp7;\
	real const tmp16 = -1. * tmp2;\
	real const tmp17 = -1. * Z_l.y;\
	real const tmp18 = -1. * tmp5;\
	real const tmp19 = K_ll.xz * tmp7;\
	real const tmp20 = -1. * tmp4;\
	real const tmp21 = -1. * Z_l.z;\
	real const tmp22 = gamma_uu.xx * Z_l.x;\
	real const tmp23 = gamma_uu.xy * Z_l.y;\
	real const tmp24 = gamma_uu.xz * Z_l.z;\
	real const tmp25 = gamma_uu.xy * gamma_uu.xy;\
	real const tmp26 = d_lll.x.yy * tmp25;\
	real const tmp27 = gamma_uu.xz * gamma_uu.xz;\
	real const tmp28 = d_lll.x.zz * tmp27;\
	real const tmp29 = gamma_uu.yy * d_lll.x.yy;\
	real const tmp30 = gamma_uu.zz * d_lll.x.zz;\
	real const tmp31 = gamma_uu.yz * d_lll.x.yz;\
	real const tmp32 = gamma_uu.xx * tmp31;\
	real const tmp33 = gamma_uu.xy * tmp3;\
	real const tmp34 = 2. * tmp32;\
	real const tmp35 = -2. * tmp33;\
	real const tmp36 = gamma_uu.xx * tmp30;\
	real const tmp37 = tmp34 + tmp35;\
	real const tmp38 = gamma_uu.xx * tmp29;\
	real const tmp39 = tmp36 + tmp37;\
	real const tmp40 = -1. * tmp28;\
	real const tmp41 = tmp38 + tmp39;\
	real const tmp42 = -1. * tmp26;\
	real const tmp43 = tmp40 + tmp41;\
	real const tmp44 = Theta * tmp7;\
	real const tmp45 = tmp42 + tmp43;\
	real const tmp46 = -1. * tmp24;\
	real const tmp47 = -1. * tmp23;\
	real const tmp48 = -1. * tmp22;\
	real const tmp49 = sqrt_f * tmp44;\
	real const tmp50 = 3. / 2.;\
	real const tmp51 = tmp7 * tmp7 * tmp7;\
	real const tmp52 = K_ll.xx * tmp51;\
	real const tmp53 = sqrt_f * tmp8;\
	real const tmp54 = sqrt_f * tmp12;\
	real const tmp55 = sqrt_f * tmp15;\
	real const tmp56 = gamma_uu.xy * tmp55;\
	real const tmp57 = sqrt_f * tmp19;\
	real const tmp58 = gamma_uu.xz * tmp57;\
	real const tmp59 = sqrt_f * tmp10;\
	real const tmp60 = gamma_uu.yz * tmp59;\
	real const tmp61 = gamma_uu.zz * tmp54;\
	real const tmp62 = gamma_uu.yy * tmp53;\
	real const tmp63 = sqrt_f * tmp52;\
	real const tmp64 = gamma_uu.xx * a_l.x;\
	(result)->ptr[0] = -1. * a_l.x + 2. * Z_l.x * tmp1 + gamma_uu.xx * d_lll.x.xx * tmp1 + -1. * gamma_uu.yy * d_lll.x.yy * tmp1 + -2. * gamma_uu.yz * d_lll.x.yz * tmp1 + -1. * gamma_uu.zz * d_lll.x.zz * tmp1;\
	(result)->ptr[1] = Z_l.y + gamma_uu.xx * d_lll.x.xy + tmp2 + tmp3;\
	(result)->ptr[2] = Z_l.z + gamma_uu.xx * d_lll.x.xz + tmp4 + tmp5;\
	(result)->ptr[3] = tmp8 + tmp9;\
	(result)->ptr[4] = -1. * tmp8 + tmp9;\
	(result)->ptr[5] = tmp10 + tmp11;\
	(result)->ptr[6] = -1. * tmp10 + tmp11;\
	(result)->ptr[7] = tmp12 + tmp13;\
	(result)->ptr[8] = -1. * tmp12 + tmp13;\
	(result)->ptr[9] = tmp14 + tmp15 + tmp16 + tmp17;\
	(result)->ptr[10] = -1. * tmp15 + tmp14 + tmp16 + tmp17;\
	(result)->ptr[11] = tmp18 + tmp19 + tmp20 + tmp21;\
	(result)->ptr[12] = -1. * tmp19 + tmp18 + tmp20 + tmp21;\
	(result)->ptr[13] = tmp44 + tmp45 + tmp46 + tmp47 + tmp48;\
	(result)->ptr[14] = -1. * tmp44 + tmp45 + tmp46 + tmp47 + tmp48;\
	(result)->ptr[15] = -2. * tmp49 + 2. * tmp56 + 2. * tmp60 + 2. * tmp58 + tmp61 + tmp62 + tmp63 + tmp64;\
	(result)->ptr[16] = 2. * tmp49 + -1. * tmp63 + -1. * tmp62 + -1. * tmp61 + -2. * tmp56 + -2. * tmp60 + -2. * tmp58 + tmp64;\
<? else ?>\
	/* BEGIN CUT from symmath/tests/Z4 - compute flux eigenmodes.symmath */\
	real const tmp1 = 1. / 2.;\
	real const tmp2 = 3. * tmp1;\
	real const tmp6 = sqrt(gamma_uu.xx);\
	real const tmp3 = tmp6 * tmp6 * tmp6;\
	real const tmp4 = K_ll.xx * tmp3;\
	real const tmp5 = sqrt_f * tmp4;\
	real const tmp7 = Theta * tmp6;\
	real const tmp8 = sqrt_f * tmp7;\
	real const tmp9 = K_ll.yy * tmp6;\
	real const tmp10 = sqrt_f * tmp9;\
	real const tmp11 = gamma_uu.yy * tmp10;\
	real const tmp12 = K_ll.zz * tmp6;\
	real const tmp13 = sqrt_f * tmp12;\
	real const tmp14 = gamma_uu.zz * tmp13;\
	real const tmp15 = K_ll.xy * tmp6;\
	real const tmp16 = sqrt_f * tmp15;\
	real const tmp17 = gamma_uu.xy * tmp16;\
	real const tmp18 = K_ll.xz * tmp6;\
	real const tmp19 = sqrt_f * tmp18;\
	real const tmp20 = gamma_uu.xz * tmp19;\
	real const tmp21 = K_ll.yz * tmp6;\
	real const tmp22 = sqrt_f * tmp21;\
	real const tmp23 = gamma_uu.yz * tmp22;\
	real const tmp24 = gamma_uu.xz * a_l.z;\
	real const tmp25 = gamma_uu.xy * a_l.y;\
	real const tmp26 = gamma_uu.xx * a_l.x;\
	real const tmp27 = gamma_uu.xy * d_lll.y.yy;\
	real const tmp28 = gamma_uu.xx * d_lll.y.xy;\
	real const tmp29 = gamma_uu.xz * d_lll.y.yz;\
	real const tmp30 = -2. * tmp28;\
	real const tmp31 = -2. * tmp29;\
	real const tmp32 = tmp30 + tmp31;\
	real const tmp33 = gamma_uu.xz * d_lll.z.yy;\
	real const tmp34 = -1. * tmp27;\
	real const tmp35 = gamma_uu.xx * d_lll.x.yy;\
	real const tmp36 = gamma_uu.xx * d_lll.y.xz;\
	real const tmp37 = gamma_uu.xx * d_lll.z.xy;\
	real const tmp38 = gamma_uu.xy * d_lll.z.yy;\
	real const tmp39 = gamma_uu.xz * d_lll.y.zz;\
	real const tmp40 = -1. * tmp39;\
	real const tmp41 = -1. * tmp38;\
	real const tmp42 = -1. * tmp37;\
	real const tmp43 = -1. * tmp36;\
	real const tmp44 = gamma_uu.xx * d_lll.x.yz;\
	real const tmp45 = gamma_uu.xz * d_lll.z.zz;\
	real const tmp46 = gamma_uu.xx * d_lll.z.xz;\
	real const tmp47 = gamma_uu.xy * d_lll.z.yz;\
	real const tmp48 = -2. * tmp46;\
	real const tmp49 = -2. * tmp47;\
	real const tmp50 = tmp48 + tmp49;\
	real const tmp51 = -1. * tmp45;\
	real const tmp52 = gamma_uu.xy * d_lll.y.zz;\
	real const tmp53 = gamma_uu.xx * d_lll.x.zz;\
	real const tmp54 = gamma_uu.xx * d_lll.y.xx;\
	real const tmp55 = tmp1 * tmp54;\
	real const tmp56 = gamma_uu.yy * d_lll.y.yy;\
	real const tmp57 = gamma_uu.zz * d_lll.y.zz;\
	real const tmp58 = gamma_uu.xy * d_lll.x.yy;\
	real const tmp59 = gamma_uu.xz * d_lll.x.yz;\
	real const tmp60 = gamma_uu.yz * d_lll.y.yz;\
	real const tmp61 = gamma_uu.xz * d_lll.z.xy;\
	real const tmp62 = -1. * tmp59;\
	real const tmp63 = gamma_uu.xy * d_lll.y.xy;\
	real const tmp64 = -1. * tmp58;\
	real const tmp65 = tmp1 * tmp57;\
	real const tmp66 = tmp1 * tmp56;\
	real const tmp67 = -1. * tmp55;\
	real const tmp68 = -1. * Z_l.y;\
	real const tmp69 = a_l.y * tmp1;\
	real const tmp70 = gamma_uu.xx * d_lll.z.xx;\
	real const tmp71 = tmp1 * tmp70;\
	real const tmp72 = gamma_uu.yy * d_lll.z.yy;\
	real const tmp73 = gamma_uu.zz * d_lll.z.zz;\
	real const tmp74 = gamma_uu.xy * d_lll.x.yz;\
	real const tmp75 = gamma_uu.xz * d_lll.x.zz;\
	real const tmp76 = gamma_uu.yz * d_lll.z.yz;\
	real const tmp77 = gamma_uu.xz * d_lll.z.xz;\
	real const tmp78 = -1. * tmp75;\
	real const tmp79 = gamma_uu.xy * d_lll.y.xz;\
	real const tmp80 = -1. * tmp74;\
	real const tmp81 = tmp1 * tmp73;\
	real const tmp82 = tmp1 * tmp72;\
	real const tmp83 = -1. * tmp71;\
	real const tmp84 = -1. * Z_l.z;\
	real const tmp85 = a_l.z * tmp1;\
	real const tmp86 = gamma_uu.xx * Z_l.x;\
	real const tmp87 = gamma_uu.xy * Z_l.y;\
	real const tmp88 = gamma_uu.xz * Z_l.z;\
	real const tmp89 = gamma_uu.xy * gamma_uu.xy;\
	real const tmp90 = d_lll.x.yy * tmp89;\
	real const tmp91 = gamma_uu.xz * gamma_uu.xz;\
	real const tmp92 = d_lll.x.zz * tmp91;\
	real const tmp93 = gamma_uu.yy * d_lll.x.yy;\
	real const tmp94 = gamma_uu.yy * d_lll.y.xy;\
	real const tmp95 = gamma_uu.xx * tmp94;\
	real const tmp96 = gamma_uu.yz * d_lll.y.xz;\
	real const tmp97 = gamma_uu.xx * tmp96;\
	real const tmp98 = gamma_uu.yz * d_lll.z.xy;\
	real const tmp99 = gamma_uu.xx * tmp98;\
	real const tmp100 = gamma_uu.zz * d_lll.x.zz;\
	real const tmp101 = gamma_uu.zz * d_lll.z.xz;\
	real const tmp102 = gamma_uu.xx * tmp101;\
	real const tmp103 = gamma_uu.xz * d_lll.y.xz;\
	real const tmp104 = gamma_uu.yz * d_lll.z.yy;\
	real const tmp105 = gamma_uu.xy * tmp104;\
	real const tmp106 = gamma_uu.zz * d_lll.z.yz;\
	real const tmp107 = gamma_uu.xy * tmp106;\
	real const tmp108 = gamma_uu.yy * d_lll.y.yz;\
	real const tmp109 = gamma_uu.xz * tmp108;\
	real const tmp110 = gamma_uu.yz * d_lll.y.zz;\
	real const tmp111 = gamma_uu.xz * tmp110;\
	real const tmp112 = gamma_uu.yz * d_lll.x.yz;\
	real const tmp113 = gamma_uu.xx * tmp112;\
	real const tmp114 = gamma_uu.xy * tmp59;\
	real const tmp115 = 2. * tmp113;\
	real const tmp116 = -2. * tmp114;\
	real const tmp117 = gamma_uu.xz * tmp76;\
	real const tmp118 = tmp115 + tmp116;\
	real const tmp119 = -1. * tmp111;\
	real const tmp120 = tmp117 + tmp118;\
	real const tmp121 = gamma_uu.xz * tmp72;\
	real const tmp122 = tmp119 + tmp120;\
	real const tmp123 = -1. * tmp109;\
	real const tmp124 = tmp121 + tmp122;\
	real const tmp125 = -1. * tmp107;\
	real const tmp126 = tmp123 + tmp124;\
	real const tmp127 = gamma_uu.xy * tmp57;\
	real const tmp128 = tmp125 + tmp126;\
	real const tmp129 = -1. * tmp105;\
	real const tmp130 = tmp127 + tmp128;\
	real const tmp131 = gamma_uu.xy * tmp60;\
	real const tmp132 = tmp129 + tmp130;\
	real const tmp133 = gamma_uu.xy * tmp61;\
	real const tmp134 = tmp131 + tmp132;\
	real const tmp135 = gamma_uu.xy * tmp103;\
	real const tmp136 = tmp133 + tmp134;\
	real const tmp137 = -1. * tmp102;\
	real const tmp138 = tmp135 + tmp136;\
	real const tmp139 = gamma_uu.xx * tmp100;\
	real const tmp140 = tmp137 + tmp138;\
	real const tmp141 = -1. * tmp99;\
	real const tmp142 = tmp139 + tmp140;\
	real const tmp143 = -1. * tmp97;\
	real const tmp144 = tmp141 + tmp142;\
	real const tmp145 = -1. * tmp95;\
	real const tmp146 = tmp143 + tmp144;\
	real const tmp147 = gamma_uu.xx * tmp93;\
	real const tmp148 = tmp145 + tmp146;\
	real const tmp149 = d_lll.z.xz * tmp91;\
	real const tmp150 = tmp147 + tmp148;\
	real const tmp151 = d_lll.y.xy * tmp89;\
	real const tmp152 = tmp149 + tmp150;\
	real const tmp153 = -1. * tmp92;\
	real const tmp154 = tmp151 + tmp152;\
	real const tmp155 = -1. * tmp90;\
	real const tmp156 = tmp153 + tmp154;\
	real const tmp157 = tmp155 + tmp156;\
	real const tmp158 = -1. * tmp88;\
	real const tmp159 = -1. * tmp87;\
	real const tmp160 = -1. * tmp86;\
	real const tmp161 = sqrt_f * sqrt_f;\
	(result)->ptr[0] = -1. * tmp5 + 2. * tmp8 + -1. * tmp11 + -1. * tmp14 + -2. * tmp17 + -2. * tmp23 + -2. * tmp20 + tmp24 + tmp25 + tmp26;\
	(result)->ptr[1] = -1. * tmp9 + tmp32 + tmp33 + tmp34 + tmp35;\
	(result)->ptr[2] = -1. * tmp21 + tmp40 + tmp41 + tmp42 + tmp43 + tmp44;\
	(result)->ptr[3] = -1. * tmp12 + tmp50 + tmp51 + tmp52 + tmp53;\
	(result)->ptr[4] = -1. * tmp15 + tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65 + tmp66 + tmp67 + tmp68 + tmp69;\
	(result)->ptr[5] = -1. * tmp18 + tmp76 + tmp77 + tmp78 + tmp79 + tmp80 + tmp81 + tmp82 + tmp83 + tmp84 + tmp85;\
	(result)->ptr[6] = -1. * tmp7 + tmp157 + tmp158 + tmp159 + tmp160;\
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
	(result)->ptr[21] = -1. * a_l.x + 2. * Z_l.x * tmp161 + gamma_uu.xx * d_lll.x.xx * tmp161 + -1. * gamma_uu.yy * d_lll.x.yy * tmp161 + -2. * gamma_uu.yz * d_lll.x.yz * tmp161 + -1. * gamma_uu.zz * d_lll.x.zz * tmp161;\
	(result)->ptr[22] = Z_l.y + gamma_uu.xx * d_lll.x.xy + tmp58 + tmp59;\
	(result)->ptr[23] = Z_l.z + gamma_uu.xx * d_lll.x.xz + tmp74 + tmp75;\
	(result)->ptr[24] = tmp9 + tmp32 + tmp33 + tmp34 + tmp35;\
	(result)->ptr[25] = tmp21 + tmp40 + tmp41 + tmp42 + tmp43 + tmp44;\
	(result)->ptr[26] = tmp12 + tmp50 + tmp51 + tmp52 + tmp53;\
	(result)->ptr[27] = tmp15 + tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65 + tmp66 + tmp67 + tmp68 + tmp69;\
	(result)->ptr[28] = tmp18 + tmp76 + tmp77 + tmp78 + tmp79 + tmp80 + tmp81 + tmp82 + tmp83 + tmp84 + tmp85;\
	(result)->ptr[29] = tmp7 + tmp157 + tmp158 + tmp159 + tmp160;\
	(result)->ptr[30] = -2. * tmp8 + 2. * tmp17 + 2. * tmp23 + 2. * tmp20 + tmp5 + tmp11 + tmp14 + tmp24 + tmp25 + tmp26;\
	/* END CUT from symmath/tests/Z4 - compute flux eigenmodes.symmath */\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=normal_t?> rotate <?=initCond_codeprefix?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */inputU,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	for (int j = 0; j < numStates; ++j) {\
		(result)->ptr[j] = 0./0.;\
	}\
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
	real const fSq = f * f;\
\
<? if eqn.noZeroRowsInFlux then ?>\
	real const tmp1 = 1. / gamma_uu.xx;\
	real const tmp2 = 1. / 2.;\
	real const tmp3 = gamma_uu.xx * gamma_uu.xx;\
	real const tmp4 = sqrt_f * sqrt_f;\
	real const tmp5 = 1. / tmp3;\
	real const tmp6 = 1. / tmp4;\
	real const tmp7 = tmp5 * tmp6;\
	real const tmp8 = (inputU)->ptr[13] * tmp1;\
	real const tmp9 = (inputU)->ptr[4] * tmp1;\
	real const tmp10 = (inputU)->ptr[14] * tmp1;\
	real const tmp11 = (inputU)->ptr[5] * tmp1;\
	real const tmp12 = (inputU)->ptr[10] * tmp1;\
	real const tmp13 = (inputU)->ptr[1] * tmp1;\
	real const tmp14 = (inputU)->ptr[11] * tmp1;\
	real const tmp15 = (inputU)->ptr[2] * tmp1;\
	real const tmp16 = (inputU)->ptr[12] * tmp1;\
	real const tmp17 = (inputU)->ptr[3] * tmp1;\
	real const tmp18 = 3. / 2.;\
	real const tmp23 = sqrt(gamma_uu.xx);\
	real const tmp19 = tmp23 * tmp23 * tmp23;\
	real const tmp20 = 1. / sqrt_f;\
	real const tmp21 = 1. / tmp19;\
	real const tmp22 = tmp20 * tmp21;\
	real const tmp24 = 1. / tmp23;\
	(result)->ptr[7] = (inputU)->ptr[16] * tmp1 * tmp2 + (inputU)->ptr[0] * tmp1 * tmp2;\
	(result)->ptr[10] = (inputU)->ptr[0] * tmp2 * tmp7 + (inputU)->ptr[16] * tmp2 * tmp7 + (inputU)->ptr[7] * tmp1 * tmp6 + (inputU)->ptr[15] * tmp5 + (inputU)->ptr[6] * tmp5 + -1. * gamma_uu.yy * (inputU)->ptr[10] * tmp2 * tmp5 + -1. * gamma_uu.yy * (inputU)->ptr[1] * tmp2 * tmp5 + -1. * gamma_uu.zz * (inputU)->ptr[12] * tmp2 * tmp5 + -1. * gamma_uu.zz * (inputU)->ptr[3] * tmp2 * tmp5 + -1. * gamma_uu.xy * (inputU)->ptr[13] * tmp5 + -1. * gamma_uu.xy * (inputU)->ptr[4] * tmp5 + -1. * gamma_uu.xz * (inputU)->ptr[14] * tmp5 + -1. * gamma_uu.xz * (inputU)->ptr[5] * tmp5 + -1. * gamma_uu.yz * (inputU)->ptr[2] * tmp5 + -1. * gamma_uu.yz * (inputU)->ptr[11] * tmp5;\
	(result)->ptr[11] = tmp2 * tmp8 + (inputU)->ptr[8] * tmp1 + tmp2 * tmp9;\
	(result)->ptr[12] = tmp2 * tmp10 + (inputU)->ptr[9] * tmp1 + tmp2 * tmp11;\
	(result)->ptr[13] = tmp2 * tmp13 + tmp2 * tmp12;\
	(result)->ptr[14] = tmp2 * tmp15 + tmp2 * tmp14;\
	(result)->ptr[15] = tmp2 * tmp17 + tmp2 * tmp16;\
	(result)->ptr[28] = -1. * (inputU)->ptr[0] * tmp2 * tmp22 + (inputU)->ptr[16] * tmp2 * tmp22 + (inputU)->ptr[15] * tmp21 + -1. * (inputU)->ptr[6] * tmp21 + -1. * gamma_uu.yy * (inputU)->ptr[10] * tmp2 * tmp21 + gamma_uu.yy * (inputU)->ptr[1] * tmp2 * tmp21 + -1. * gamma_uu.zz * (inputU)->ptr[12] * tmp2 * tmp21 + gamma_uu.zz * (inputU)->ptr[3] * tmp2 * tmp21 + -1. * gamma_uu.xy * (inputU)->ptr[13] * tmp21 + gamma_uu.xy * (inputU)->ptr[4] * tmp21 + -1. * gamma_uu.xz * (inputU)->ptr[14] * tmp21 + gamma_uu.xz * (inputU)->ptr[5] * tmp21 + gamma_uu.yz * (inputU)->ptr[2] * tmp21 + -1. * gamma_uu.yz * (inputU)->ptr[11] * tmp21;\
	(result)->ptr[29] = -1. * (inputU)->ptr[4] * tmp2 * tmp24 + (inputU)->ptr[13] * tmp2 * tmp24;\
	(result)->ptr[30] = -1. * (inputU)->ptr[5] * tmp2 * tmp24 + (inputU)->ptr[14] * tmp2 * tmp24;\
	(result)->ptr[31] = -1. * (inputU)->ptr[1] * tmp2 * tmp24 + (inputU)->ptr[10] * tmp2 * tmp24;\
	(result)->ptr[32] = -1. * (inputU)->ptr[2] * tmp2 * tmp24 + (inputU)->ptr[11] * tmp2 * tmp24;\
	(result)->ptr[33] = -1. * (inputU)->ptr[3] * tmp2 * tmp24 + (inputU)->ptr[12] * tmp2 * tmp24;\
	(result)->ptr[34] = -1. * (inputU)->ptr[6] * tmp2 * tmp24 + (inputU)->ptr[15] * tmp2 * tmp24;\
	(result)->ptr[35] = -1. * (inputU)->ptr[15] * tmp1 * tmp2 + -1. * (inputU)->ptr[6] * tmp1 * tmp2 + gamma_uu.xy * tmp2 * tmp8 + gamma_uu.xy * tmp2 * tmp9 + gamma_uu.xz * tmp2 * tmp10 + gamma_uu.xz * tmp2 * tmp11 + gamma_uu.yy * tmp2 * tmp12 + gamma_uu.yy * tmp2 * tmp13 + gamma_uu.zz * tmp2 * tmp16 + gamma_uu.zz * tmp2 * tmp17 + gamma_uu.yz * tmp15 + gamma_uu.yz * tmp14;\
	(result)->ptr[36] = -1. * (inputU)->ptr[13] * tmp2 + -1. * (inputU)->ptr[4] * tmp2 + -1. * gamma_uu.xy * tmp2 * tmp12 + -1. * gamma_uu.xy * tmp2 * tmp13 + -1. * gamma_uu.xz * tmp2 * tmp15 + -1. * gamma_uu.xz * tmp2 * tmp14;\
	(result)->ptr[37] = -1. * (inputU)->ptr[14] * tmp2 + -1. * (inputU)->ptr[5] * tmp2 + -1. * gamma_uu.xy * tmp2 * tmp14 + -1. * gamma_uu.xy * tmp2 * tmp15 + -1. * gamma_uu.xz * tmp2 * tmp17 + -1. * gamma_uu.xz * tmp2 * tmp16;\
<? else ?>\
	/* BEGIN CUT from symmath/tests/Z4 - compute flux eigenmodes.symmath */\
	real const tmp1 = 1. / gamma_uu.xx;\
	real const tmp2 = 1. / 2.;\
	real const tmp3 = (inputU)->ptr[7] * tmp1;\
	real const tmp4 = gamma_uu.xy * tmp3;\
	real const tmp5 = (inputU)->ptr[8] * tmp1;\
	real const tmp6 = gamma_uu.xz * tmp5;\
	real const tmp7 = gamma_uu.xx * gamma_uu.xx;\
	real const tmp8 = sqrt_f * sqrt_f;\
	real const tmp9 = 1. / tmp7;\
	real const tmp10 = 1. / tmp8;\
	real const tmp11 = tmp9 * tmp10;\
	real const tmp12 = (inputU)->ptr[27] * tmp1;\
	real const tmp13 = (inputU)->ptr[4] * tmp1;\
	real const tmp14 = (inputU)->ptr[12] * tmp1;\
	real const tmp15 = gamma_uu.yy * tmp14;\
	real const tmp16 = (inputU)->ptr[14] * tmp1;\
	real const tmp17 = gamma_uu.zz * tmp16;\
	real const tmp18 = (inputU)->ptr[13] * tmp1;\
	real const tmp19 = (inputU)->ptr[28] * tmp1;\
	real const tmp20 = (inputU)->ptr[5] * tmp1;\
	real const tmp21 = (inputU)->ptr[18] * tmp1;\
	real const tmp22 = gamma_uu.yy * tmp21;\
	real const tmp23 = (inputU)->ptr[20] * tmp1;\
	real const tmp24 = gamma_uu.zz * tmp23;\
	real const tmp25 = (inputU)->ptr[19] * tmp1;\
	real const tmp26 = (inputU)->ptr[1] * tmp1;\
	real const tmp27 = (inputU)->ptr[24] * tmp1;\
	real const tmp28 = gamma_uu.xz * tmp18;\
	real const tmp29 = (inputU)->ptr[25] * tmp1;\
	real const tmp30 = (inputU)->ptr[2] * tmp1;\
	real const tmp31 = (inputU)->ptr[26] * tmp1;\
	real const tmp32 = (inputU)->ptr[3] * tmp1;\
	real const tmp33 = 3. / 2.;\
	real const tmp38 = sqrt(gamma_uu.xx);\
	real const tmp34 = tmp38 * tmp38 * tmp38;\
	real const tmp35 = 1. / sqrt_f;\
	real const tmp36 = 1. / tmp34;\
	real const tmp37 = tmp35 * tmp36;\
	real const tmp39 = 1. / tmp38;\
	real const tmp40 = gamma_uu.xy * gamma_uu.xy;\
	real const tmp41 = tmp1 * tmp40;\
	real const tmp42 = gamma_uu.xz * gamma_uu.xz;\
	real const tmp43 = tmp1 * tmp42;\
	(result)->ptr[7+0] = (inputU)->ptr[0] * tmp1 * tmp2 + (inputU)->ptr[30] * tmp1 * tmp2 + -1. * tmp6 + -1. * tmp4;\
	(result)->ptr[7+1] = (inputU)->ptr[7];\
	(result)->ptr[7+2] = (inputU)->ptr[8];\
	(result)->ptr[7+3] = (inputU)->ptr[0] * tmp2 * tmp11 + (inputU)->ptr[30] * tmp2 * tmp11 + (inputU)->ptr[21] * tmp1 * tmp10 + -1. * gamma_uu.xy * (inputU)->ptr[7] * tmp11 + -1. * gamma_uu.xz * (inputU)->ptr[8] * tmp11 + (inputU)->ptr[29] * tmp9 + (inputU)->ptr[6] * tmp9 + -1. * gamma_uu.yy * (inputU)->ptr[1] * tmp2 * tmp9 + -1. * gamma_uu.yy * (inputU)->ptr[24] * tmp2 * tmp9 + -1. * gamma_uu.zz * (inputU)->ptr[26] * tmp2 * tmp9 + -1. * gamma_uu.zz * (inputU)->ptr[3] * tmp2 * tmp9 + -1. * gamma_uu.xy * (inputU)->ptr[27] * tmp9 + -1. * gamma_uu.xy * (inputU)->ptr[4] * tmp9 + gamma_uu.xy * (inputU)->ptr[7] * tmp9 + -1. * gamma_uu.xz * (inputU)->ptr[28] * tmp9 + -1. * gamma_uu.xz * (inputU)->ptr[5] * tmp9 + gamma_uu.xz * (inputU)->ptr[8] * tmp9 + -1. * gamma_uu.yz * (inputU)->ptr[25] * tmp9 + -1. * gamma_uu.yz * (inputU)->ptr[2] * tmp9 + -1. * gamma_uu.xz * (inputU)->ptr[15] * tmp1 + -1. * gamma_uu.xy * (inputU)->ptr[9] * tmp1;\
	(result)->ptr[7+4] = tmp2 * tmp12 + tmp2 * tmp13 + -1. * tmp2 * tmp3 + (inputU)->ptr[22] * tmp1 + (inputU)->ptr[9] * tmp2 + -1. * tmp2 * tmp15 + -1. * tmp2 * tmp17 + -1. * gamma_uu.xy * (inputU)->ptr[10] * tmp1 + -1. * gamma_uu.yz * tmp18 + -1. * gamma_uu.xz * (inputU)->ptr[16] * tmp1;\
	(result)->ptr[7+5] = tmp2 * tmp19 + tmp2 * tmp20 + -1. * tmp2 * tmp5 + (inputU)->ptr[23] * tmp1 + (inputU)->ptr[15] * tmp2 + -1. * tmp2 * tmp22 + -1. * tmp2 * tmp24 + -1. * gamma_uu.xy * (inputU)->ptr[11] * tmp1 + -1. * gamma_uu.yz * tmp25 + -1. * gamma_uu.xz * (inputU)->ptr[17] * tmp1;\
	(result)->ptr[7+6] = tmp2 * tmp26 + tmp2 * tmp27 + gamma_uu.xy * tmp14 + -1. * gamma_uu.xz * tmp21 + 2. * tmp28 + 2. * (inputU)->ptr[10];\
	(result)->ptr[7+7] = tmp2 * tmp29 + (inputU)->ptr[11] + tmp2 * tmp30 + (inputU)->ptr[16] + gamma_uu.xz * tmp16 + gamma_uu.xy * tmp21;\
	(result)->ptr[7+8] = tmp2 * tmp31 + tmp2 * tmp32 + -1. * gamma_uu.xy * tmp16 + gamma_uu.xz * tmp23 + 2. * gamma_uu.xy * tmp25 + 2. * (inputU)->ptr[17];\
	(result)->ptr[7+9] = (inputU)->ptr[9];\
	(result)->ptr[7+10] = (inputU)->ptr[10];\
	(result)->ptr[7+11] = (inputU)->ptr[11];\
	(result)->ptr[7+12] = (inputU)->ptr[12];\
	(result)->ptr[7+13] = (inputU)->ptr[13];\
	(result)->ptr[7+14] = (inputU)->ptr[14];\
	(result)->ptr[7+15] = (inputU)->ptr[15];\
	(result)->ptr[7+16] = (inputU)->ptr[16];\
	(result)->ptr[7+17] = (inputU)->ptr[17];\
	(result)->ptr[7+18] = (inputU)->ptr[18];\
	(result)->ptr[7+19] = (inputU)->ptr[19];\
	(result)->ptr[7+20] = (inputU)->ptr[20];\
	(result)->ptr[7+21] = -1. * (inputU)->ptr[0] * tmp2 * tmp37 + (inputU)->ptr[30] * tmp2 * tmp37 + (inputU)->ptr[29] * tmp36 + -1. * (inputU)->ptr[6] * tmp36 + gamma_uu.yy * (inputU)->ptr[1] * tmp2 * tmp36 + -1. * gamma_uu.yy * (inputU)->ptr[24] * tmp2 * tmp36 + -1. * gamma_uu.zz * (inputU)->ptr[26] * tmp2 * tmp36 + gamma_uu.zz * (inputU)->ptr[3] * tmp2 * tmp36 + -1. * gamma_uu.xy * (inputU)->ptr[27] * tmp36 + gamma_uu.xy * (inputU)->ptr[4] * tmp36 + -1. * gamma_uu.xz * (inputU)->ptr[28] * tmp36 + gamma_uu.xz * (inputU)->ptr[5] * tmp36 + gamma_uu.yz * (inputU)->ptr[2] * tmp36 + -1. * gamma_uu.yz * (inputU)->ptr[25] * tmp36;\
	(result)->ptr[7+22] = -1. * (inputU)->ptr[4] * tmp2 * tmp39 + (inputU)->ptr[27] * tmp2 * tmp39;\
	(result)->ptr[7+23] = -1. * (inputU)->ptr[5] * tmp2 * tmp39 + (inputU)->ptr[28] * tmp2 * tmp39;\
	(result)->ptr[7+24] = (inputU)->ptr[24] * tmp2 * tmp39 + -1. * (inputU)->ptr[1] * tmp2 * tmp39;\
	(result)->ptr[7+25] = -1. * (inputU)->ptr[2] * tmp2 * tmp39 + (inputU)->ptr[25] * tmp2 * tmp39;\
	(result)->ptr[7+26] = -1. * (inputU)->ptr[3] * tmp2 * tmp39 + (inputU)->ptr[26] * tmp2 * tmp39;\
	(result)->ptr[7+27] = -1. * (inputU)->ptr[6] * tmp2 * tmp39 + (inputU)->ptr[29] * tmp2 * tmp39;\
	(result)->ptr[7+28] = -1. * (inputU)->ptr[29] * tmp1 * tmp2 + -1. * (inputU)->ptr[6] * tmp1 * tmp2 + gamma_uu.xy * tmp2 * tmp12 + gamma_uu.xy * tmp2 * tmp13 + -1. * tmp2 * tmp4 + gamma_uu.xz * tmp2 * tmp19 + gamma_uu.xz * tmp2 * tmp20 + -1. * tmp2 * tmp6 + gamma_uu.yy * tmp2 * tmp26 + gamma_uu.yy * tmp2 * tmp27 + gamma_uu.zz * tmp2 * tmp31 + gamma_uu.zz * tmp2 * tmp32 + gamma_uu.yz * tmp29 + gamma_uu.yz * tmp30 + gamma_uu.xy * (inputU)->ptr[9] * tmp2 + gamma_uu.xz * (inputU)->ptr[15] * tmp2 + gamma_uu.xy * tmp2 * tmp15 + -1. * gamma_uu.xy * tmp2 * tmp17 + -1. * gamma_uu.xz * tmp2 * tmp22 + gamma_uu.xz * tmp2 * tmp24 + gamma_uu.yy * (inputU)->ptr[10] + gamma_uu.yz * (inputU)->ptr[11] + gamma_uu.yz * (inputU)->ptr[16] + gamma_uu.zz * (inputU)->ptr[17] + gamma_uu.xy * gamma_uu.yz * tmp21 + gamma_uu.xy * gamma_uu.zz * tmp25 + gamma_uu.xz * gamma_uu.yz * tmp16 + gamma_uu.xz * gamma_uu.yy * tmp18;\
	(result)->ptr[7+29] = -1. * (inputU)->ptr[27] * tmp2 + -1. * (inputU)->ptr[4] * tmp2 + (inputU)->ptr[7] * tmp2 + -1. * gamma_uu.xy * tmp2 * tmp26 + -1. * gamma_uu.xy * tmp2 * tmp27 + -1. * gamma_uu.xz * tmp2 * tmp29 + -1. * gamma_uu.xz * tmp2 * tmp30 + -1. * gamma_uu.xx * (inputU)->ptr[9] * tmp2 + -1. * (inputU)->ptr[12] * tmp41 + -1. * (inputU)->ptr[14] * tmp43 + gamma_uu.yy * (inputU)->ptr[12] * tmp2 + gamma_uu.zz * (inputU)->ptr[14] * tmp2 + -1. * gamma_uu.xy * (inputU)->ptr[10] + -1. * gamma_uu.xz * (inputU)->ptr[11] + -2. * gamma_uu.xy * tmp28 + gamma_uu.yz * (inputU)->ptr[13];\
	(result)->ptr[7+30] = -1. * (inputU)->ptr[28] * tmp2 + -1. * (inputU)->ptr[5] * tmp2 + (inputU)->ptr[8] * tmp2 + -1. * gamma_uu.xy * tmp2 * tmp29 + -1. * gamma_uu.xy * tmp2 * tmp30 + -1. * gamma_uu.xz * tmp2 * tmp31 + -1. * gamma_uu.xz * tmp2 * tmp32 + -1. * gamma_uu.xx * (inputU)->ptr[15] * tmp2 + -1. * (inputU)->ptr[18] * tmp41 + -1. * (inputU)->ptr[20] * tmp43 + gamma_uu.yy * (inputU)->ptr[18] * tmp2 + gamma_uu.zz * (inputU)->ptr[20] * tmp2 + -1. * gamma_uu.xy * (inputU)->ptr[16] + -1. * gamma_uu.xz * (inputU)->ptr[17] + -2. * gamma_uu.xy * gamma_uu.xz * tmp25 + gamma_uu.yz * (inputU)->ptr[19];\
	/* END CUT from symmath/tests/Z4 - compute flux eigenmodes.symmath */\
<? end ?>\
\
	/* rotate back into flux direction */\
	/* TODO copy the rotation for non-cartesian from fluxFromCons */\
	/* TODO make it into a function */\
	(result)->a_l = real3_swap((result)->a_l, n.side);			/* 0-2 */\
	(result)->dDelta_lll = _3sym3_swap((result)->dDelta_lll, n.side);		/* 3-20 */\
	(result)->K_ll = sym3_swap((result)->K_ll, n.side);			/* 21-26 */\
	(result)->Theta = (result)->Theta;							/* 27 */\
	(result)->Z_l = real3_swap((result)->Z_l, n.side);			/* 28-30 */\
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

#error "does anyone even use this?"
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


<? if eqn.useShift == "HarmonicParabolic" then ?>
	
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

<? end -- eqn.useShift == "HarmonicParabolic" ?>


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
		
		<? if eqn.useShift == "MinimalDistortionHyperbolic" then ?>
//// MODULE_DEPENDS: mdeShiftEpsilon
	real3x3 const DBeta_ul = real3x3_add(b_ul, real3_3sym3_dot2(beta_u, conn_ull));
	sym3 const DBeta_uu = real3x3_sym3_to_sym3_mul(DBeta_ul, gamma_uu);
	real const tr_DBeta = real3x3_trace(DBeta_ul);
	real3x3x3 const conn_uul = _3sym3_sym3_mul(conn_ull, gamma_uu);
	//TODO why not use b_ll as a state var? then we'd be storing 6 instead of 9 vars.
	sym3 const b_uu = real3x3_sym3_to_sym3_mul(b_ul, gamma_uu);
	sym3 const A_uu = sym3_sub(K_uu, sym3_real_mul(gamma_uu, tr_K / 3.));
		<? end ?>
		
		<? if eqn.useShift == "GammaDriverHyperbolic" then ?>
//// MODULE_DEPENDS: <?=calc_gammaHat_ll?>
	sym3 const gammaHat_ll = <?=calc_gammaHat_ll?>(cell->pos);
	// W = (_γ/γ)^(1/6)
	// log(W)_,i = 1/3 (log(√_γ)_,i - log(√γ)_,i) = 1/3 (_Γ^j_ji - Γ^j_ji)
	real const det_gammaHat = sym3_det(gammaHat_ll);	//TODO use coord module for math simplifications?
	real const W = pow(det_gammaHat / det_gamma, 1./6.);
	real const invW = 1. / W;
//// MODULE_DEPENDS: <?=coord_conn_ull?>
	real const gammaDriver_eta = 1.;

	_3sym3 const connHat_ull = coord_conn_ull((cell)->pos);\
	real3x3 const DHatBeta_ul = real3x3_add(b_ul, real3_3sym3_dot2(beta_u, connHat_ull));
	real const tr_DHatBeta = real3x3_trace(DHatBeta_ul);
	
	//ok so here's an error in my analytical/codegen ... ^Γ^i_jk is raised/lowered by ^γ_ij
	//so this "connHat_uul" had its first index raised by γ_ij and second index raised by ^γ_ij
	real3x3x3 const connHat_uul = _3sym3_sym3_mul(connHat_ull, gamma_uu);
	
	sym3 const A_uu = sym3_sub(K_uu, sym3_real_mul(gamma_uu, tr_K / 3.));
	
	_3sym3 const neg2d_luu = _3sym3_real_mul(d_luu, -2.);	// TODO this is partial_uul ... do this in analytics
	real3 const dDelta_l = _3sym3_sym3_dot23(dDelta_lll, gamma_uu);
	real3 const dDelta_u = sym3_real3_mul(gamma_uu, dDelta_l);

	//(_Γ->Γ)^i_jk = Γ^i_jk - _Γ^i_jk = 1/3 (δ^i_j Δd_k + δ^i_k Δd_j - γ_jk Δd^i)
	_3sym3 const connToBar_ull = (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>			.<?=xjk?> = (0 
				+ <? if i == j then ?>dDelta_l.<?=xk?><? else ?>0.<? end ?>
				+ <? if i == k then ?>dDelta_l.<?=xj?><? else ?>0.<? end ?>
				- gamma_ll.<?=xjk?> * dDelta_u.<?=xi?>
			) / 3.,
<?	end
?>		},
<? end
?>	};

	_3sym3 const connBar_ull =_3sym3_sub(conn_ull, connToBar_ull);

	sym3 const gammaBar_uu = sym3_real_mul(gamma_uu, invW * invW);
	_3sym3 const DeltaGamma_ull = _3sym3_sub(connBar_ull, connHat_ull);
	real3 const DeltaGamma_u = _3sym3_sym3_dot23(DeltaGamma_ull, gammaBar_uu);
	real3 const LambdaBar_u = DeltaGamma_u;	 // ... plus C^i 

	//same problem here as connHat_uul
	real3x3x3 const DeltaGamma_uul = _3sym3_sym3_mul(DeltaGamma_ull, gamma_uu);
		<? end ?>
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
		{
			(deriv)->alpha += alpha * (beta_u.x * a_l.x + beta_u.z * a_l.z + beta_u.y * a_l.y);
			(deriv)->gammaDelta_ll.xx += 2. * (b_ll.xx + beta_u.x * d_lll.x.xx + beta_u.z * d_lll.z.xx + beta_u.y * d_lll.y.xx);
			(deriv)->gammaDelta_ll.xy += 2. * (b_ll.xy + beta_u.x * d_lll.x.xy + beta_u.z * d_lll.z.xy + beta_u.y * d_lll.y.xy);
			(deriv)->gammaDelta_ll.xz += 2. * (b_ll.xz + beta_u.x * d_lll.x.xz + beta_u.z * d_lll.z.xz + beta_u.y * d_lll.y.xz);
			(deriv)->gammaDelta_ll.yy += 2. * (b_ll.yy + beta_u.x * d_lll.x.yy + beta_u.z * d_lll.z.yy + beta_u.y * d_lll.y.yy);
			(deriv)->gammaDelta_ll.yz += 2. * (b_ll.yz + beta_u.x * d_lll.x.yz + beta_u.z * d_lll.z.yz + beta_u.y * d_lll.y.yz);
			(deriv)->gammaDelta_ll.zz += 2. * (b_ll.zz + beta_u.x * d_lll.x.zz + beta_u.z * d_lll.z.zz + beta_u.y * d_lll.y.zz);
			(deriv)->a_l.x += a_l.x * b_ul.x.x + -a_l.x * tr_b + a_l.z * b_ul.x.z + a_l.y * b_ul.x.y;
			(deriv)->a_l.y += a_l.x * b_ul.x.y + a_l.y * b_ul.y.y + a_l.z * b_ul.y.z + -a_l.y * tr_b;
			(deriv)->a_l.z += a_l.x * b_ul.x.z + a_l.y * b_ul.y.z + -a_l.z * tr_b + a_l.z * b_ul.z.z;
			(deriv)->dDelta_lll.x.xx += b_ul.x.x * dDelta_lll.x.xx + b_ul.x.y * dDelta_lll.y.xx + -dDelta_lll.x.xx * tr_b + b_ul.x.z * dDelta_lll.z.xx;
			(deriv)->dDelta_lll.x.xy += b_ul.x.x * dDelta_lll.x.xy + b_ul.x.y * dDelta_lll.y.xy + -dDelta_lll.x.xy * tr_b + b_ul.x.z * dDelta_lll.z.xy;
			(deriv)->dDelta_lll.x.xz += b_ul.x.x * dDelta_lll.x.xz + b_ul.x.y * dDelta_lll.y.xz + -dDelta_lll.x.xz * tr_b + b_ul.x.z * dDelta_lll.z.xz;
			(deriv)->dDelta_lll.x.yy += b_ul.x.x * dDelta_lll.x.yy + b_ul.x.y * dDelta_lll.y.yy + -dDelta_lll.x.yy * tr_b + b_ul.x.z * dDelta_lll.z.yy;
			(deriv)->dDelta_lll.x.yz += b_ul.x.x * dDelta_lll.x.yz + b_ul.x.y * dDelta_lll.y.yz + -dDelta_lll.x.yz * tr_b + b_ul.x.z * dDelta_lll.z.yz;
			(deriv)->dDelta_lll.x.zz += b_ul.x.x * dDelta_lll.x.zz + b_ul.x.y * dDelta_lll.y.zz + -dDelta_lll.x.zz * tr_b + b_ul.x.z * dDelta_lll.z.zz;
			(deriv)->dDelta_lll.y.xx += b_ul.x.y * dDelta_lll.x.xx + b_ul.y.y * dDelta_lll.y.xx + -dDelta_lll.y.xx * tr_b + b_ul.y.z * dDelta_lll.z.xx;
			(deriv)->dDelta_lll.y.xy += b_ul.x.y * dDelta_lll.x.xy + b_ul.y.y * dDelta_lll.y.xy + -dDelta_lll.y.xy * tr_b + b_ul.y.z * dDelta_lll.z.xy;
			(deriv)->dDelta_lll.y.xz += b_ul.x.y * dDelta_lll.x.xz + b_ul.y.y * dDelta_lll.y.xz + -dDelta_lll.y.xz * tr_b + b_ul.y.z * dDelta_lll.z.xz;
			(deriv)->dDelta_lll.y.yy += b_ul.x.y * dDelta_lll.x.yy + b_ul.y.y * dDelta_lll.y.yy + -dDelta_lll.y.yy * tr_b + b_ul.y.z * dDelta_lll.z.yy;
			(deriv)->dDelta_lll.y.yz += b_ul.x.y * dDelta_lll.x.yz + b_ul.y.y * dDelta_lll.y.yz + -dDelta_lll.y.yz * tr_b + b_ul.y.z * dDelta_lll.z.yz;
			(deriv)->dDelta_lll.y.zz += b_ul.x.y * dDelta_lll.x.zz + b_ul.y.y * dDelta_lll.y.zz + -dDelta_lll.y.zz * tr_b + b_ul.y.z * dDelta_lll.z.zz;
			(deriv)->dDelta_lll.z.xx += b_ul.x.z * dDelta_lll.x.xx + b_ul.y.z * dDelta_lll.y.xx + -dDelta_lll.z.xx * tr_b + b_ul.z.z * dDelta_lll.z.xx;
			(deriv)->dDelta_lll.z.xy += b_ul.x.z * dDelta_lll.x.xy + b_ul.y.z * dDelta_lll.y.xy + -dDelta_lll.z.xy * tr_b + b_ul.z.z * dDelta_lll.z.xy;
			(deriv)->dDelta_lll.z.xz += b_ul.x.z * dDelta_lll.x.xz + b_ul.y.z * dDelta_lll.y.xz + -dDelta_lll.z.xz * tr_b + b_ul.z.z * dDelta_lll.z.xz;
			(deriv)->dDelta_lll.z.yy += b_ul.x.z * dDelta_lll.x.yy + b_ul.y.z * dDelta_lll.y.yy + -dDelta_lll.z.yy * tr_b + b_ul.z.z * dDelta_lll.z.yy;
			(deriv)->dDelta_lll.z.yz += b_ul.x.z * dDelta_lll.x.yz + b_ul.y.z * dDelta_lll.y.yz + -dDelta_lll.z.yz * tr_b + b_ul.z.z * dDelta_lll.z.yz;
			(deriv)->dDelta_lll.z.zz += b_ul.x.z * dDelta_lll.x.zz + b_ul.y.z * dDelta_lll.y.zz + -dDelta_lll.z.zz * tr_b + b_ul.z.z * dDelta_lll.z.zz;
			(deriv)->K_ll.xx += -K_ll.xx * tr_b + 2. * K_ll.xx * b_ul.x.x + 2. * K_ll.xz * b_ul.x.z + 2. * K_ll.xy * b_ul.x.y;
			(deriv)->K_ll.xy += K_ll.xx * b_ul.x.y + K_ll.xy * b_ul.x.x + K_ll.xy * b_ul.y.y + -K_ll.xy * tr_b + K_ll.xz * b_ul.y.z + K_ll.yz * b_ul.x.z + K_ll.yy * b_ul.x.y;
			(deriv)->K_ll.xz += K_ll.xx * b_ul.x.z + K_ll.xy * b_ul.y.z + K_ll.xz * b_ul.x.x + K_ll.xz * b_ul.z.z + -K_ll.xz * tr_b + K_ll.zz * b_ul.x.z + K_ll.yz * b_ul.x.y;
			(deriv)->K_ll.yy += -K_ll.yy * tr_b + 2. * K_ll.xy * b_ul.x.y + 2. * K_ll.yz * b_ul.y.z + 2. * K_ll.yy * b_ul.y.y;
			(deriv)->K_ll.yz += K_ll.xy * b_ul.x.z + K_ll.xz * b_ul.x.y + K_ll.yy * b_ul.y.z + K_ll.yz * b_ul.y.y + K_ll.yz * b_ul.z.z + K_ll.zz * b_ul.y.z + -K_ll.yz * tr_b;
			(deriv)->K_ll.zz += -K_ll.zz * tr_b + 2. * K_ll.xz * b_ul.x.z + 2. * K_ll.zz * b_ul.z.z + 2. * K_ll.yz * b_ul.y.z;
			(deriv)->Theta += -Theta * tr_b;
			(deriv)->Z_l.x += Z_l.x * b_ul.x.x + -Z_l.x * tr_b + Z_l.z * b_ul.x.z + Z_l.y * b_ul.x.y;
			(deriv)->Z_l.y += Z_l.x * b_ul.x.y + Z_l.y * b_ul.y.y + Z_l.z * b_ul.y.z + -Z_l.y * tr_b;
			(deriv)->Z_l.z += Z_l.x * b_ul.x.z + Z_l.y * b_ul.y.z + -Z_l.z * tr_b + Z_l.z * b_ul.z.z;
		}
		<? if eqn.useShift == "HarmonicParabolic" then ?>
		{
			real const tmp1 = alpha * alpha;
			(deriv)->beta_u.x += beta_u.x * b_ul.x.x + beta_u.y * b_ul.x.y + beta_u.z * b_ul.x.z + -a_u.x * tmp1 + 2. * e_u.x * tmp1 + -d_u.x * tmp1;
			(deriv)->beta_u.y += beta_u.x * b_ul.x.y + beta_u.y * b_ul.y.y + beta_u.z * b_ul.y.z + -a_u.y * tmp1 + 2. * e_u.y * tmp1 + -d_u.y * tmp1;
			(deriv)->beta_u.z += beta_u.x * b_ul.x.z + beta_u.y * b_ul.y.z + beta_u.z * b_ul.z.z + -a_u.z * tmp1 + 2. * e_u.z * tmp1 + -d_u.z * tmp1;
		}
		<? end ?>/* eqn.useShift == "HarmonicParabolic" */
		<? if eqn.useShift == "HarmonicHyperbolic" then ?>
		{
			real const tmp1 = alpha * alpha;
			(deriv)->beta_u.x += B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y;
			(deriv)->beta_u.y += B_u.y + beta_u.x * b_ul.x.y + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y;
			(deriv)->beta_u.z += B_u.z + beta_u.x * b_ul.x.z + beta_u.z * b_ul.z.z + beta_u.y * b_ul.y.z;
			(deriv)->B_u.x += B_u.x * b_ul.x.x + B_u.y * b_ul.x.y + B_u.z * b_ul.x.z + -a_u.x * tmp1 + 2. * e_u.x * tmp1 + -d_u.x * tmp1;
			(deriv)->B_u.y += B_u.x * b_ul.x.y + B_u.y * b_ul.y.y + B_u.z * b_ul.y.z + -a_u.y * tmp1 + 2. * e_u.y * tmp1 + -d_u.y * tmp1;
			(deriv)->B_u.z += B_u.x * b_ul.x.z + B_u.y * b_ul.y.z + B_u.z * b_ul.z.z + -a_u.z * tmp1 + 2. * e_u.z * tmp1 + -d_u.z * tmp1;
		}
		<? end ?>/* eqn.useShift == "HarmonicHyperbolic" */
		<? if eqn.useShift == "MinimalDistortionHyperbolic" then ?>
		{
			real const tmp1 = b_ul.x.x * mdeShiftEpsilon;
			real const tmp2 = b_ul.y.y * mdeShiftEpsilon;
			real const tmp3 = b_ul.z.z * mdeShiftEpsilon;
			real const tmp4 = b_uu.xy * mdeShiftEpsilon;
			real const tmp5 = b_uu.xz * mdeShiftEpsilon;
			real const tmp6 = tr_DBeta * mdeShiftEpsilon;
			real const tmp7 = DBeta_uu.xy * mdeShiftEpsilon;
			real const tmp8 = DBeta_uu.xz * mdeShiftEpsilon;
			real const tmp9 = conn_ull.x.xx * mdeShiftEpsilon;
			real const tmp10 = conn_uul.x.x.x * mdeShiftEpsilon;
			real const tmp11 = conn_uul.y.x.x * mdeShiftEpsilon;
			real const tmp12 = conn_uul.z.x.x * mdeShiftEpsilon;
			real const tmp13 = conn_uul.y.x.y * mdeShiftEpsilon;
			real const tmp14 = conn_uul.z.x.y * mdeShiftEpsilon;
			real const tmp15 = conn_uul.z.x.z * mdeShiftEpsilon;
			real const tmp16 = conn_uul.x.x.y * mdeShiftEpsilon;
			real const tmp17 = conn_uul.x.x.z * mdeShiftEpsilon;
			real const tmp18 = conn_uul.y.x.z * mdeShiftEpsilon;
			real const tmp19 = conn_ull.x.xy * mdeShiftEpsilon;
			real const tmp20 = conn_ull.x.xz * mdeShiftEpsilon;
			real const tmp21 = e_u.x * mdeShiftEpsilon;
			real const tmp22 = d_luu.x.xx * mdeShiftEpsilon;
			real const tmp23 = d_luu.x.xy * mdeShiftEpsilon;
			real const tmp24 = d_luu.x.xz * mdeShiftEpsilon;
			real const tmp25 = conn_ull.y.xx * mdeShiftEpsilon;
			real const tmp26 = conn_ull.z.xx * mdeShiftEpsilon;
			real const tmp27 = conn_ull.y.xy * mdeShiftEpsilon;
			real const tmp28 = conn_ull.z.xy * mdeShiftEpsilon;
			real const tmp29 = conn_ull.z.xz * mdeShiftEpsilon;
			real const tmp30 = d_luu.y.xx * mdeShiftEpsilon;
			real const tmp31 = d_luu.y.xy * mdeShiftEpsilon;
			real const tmp32 = d_luu.y.xz * mdeShiftEpsilon;
			real const tmp33 = conn_ull.y.yy * mdeShiftEpsilon;
			real const tmp34 = conn_ull.z.yy * mdeShiftEpsilon;
			real const tmp35 = conn_ull.z.yz * mdeShiftEpsilon;
			real const tmp36 = d_luu.z.xx * mdeShiftEpsilon;
			real const tmp37 = d_luu.z.xy * mdeShiftEpsilon;
			real const tmp38 = d_luu.z.xz * mdeShiftEpsilon;
			real const tmp39 = conn_ull.y.xz * mdeShiftEpsilon;
			real const tmp40 = conn_ull.y.yz * mdeShiftEpsilon;
			real const tmp41 = conn_ull.z.zz * mdeShiftEpsilon;
			real const tmp42 = conn_ull.x.xy * tmp23;
			real const tmp43 = conn_ull.x.xz * tmp24;
			real const tmp44 = conn_ull.y.xy * tmp31;
			real const tmp45 = conn_ull.y.xz * tmp32;
			real const tmp46 = conn_ull.z.xy * tmp37;
			real const tmp47 = conn_ull.z.xz * tmp38;
			real const tmp48 = alpha * mdeShiftEpsilon;
			real const tmp49 = conn_ull.y.xy * tmp48;
			real const tmp50 = conn_ull.z.xz * tmp48;
			real const tmp51 = conn_ull.y.yy * tmp48;
			real const tmp52 = conn_ull.z.yz * tmp48;
			real const tmp53 = conn_ull.y.yz * tmp48;
			real const tmp54 = conn_ull.z.zz * tmp48;
			real const tmp55 = conn_ull.x.xy * tmp48;
			real const tmp56 = conn_ull.x.xz * tmp48;
			real const tmp57 = conn_ull.x.xx * tmp48;
			real const tmp58 = b_uu.yz * mdeShiftEpsilon;
			real const tmp59 = DBeta_uu.yz * mdeShiftEpsilon;
			real const tmp60 = conn_uul.x.x.y * mdeShiftEpsilon;
			real const tmp61 = conn_uul.y.x.y * mdeShiftEpsilon;
			real const tmp62 = conn_uul.z.x.y * mdeShiftEpsilon;
			real const tmp63 = conn_uul.y.y.y * mdeShiftEpsilon;
			real const tmp64 = conn_uul.z.y.y * mdeShiftEpsilon;
			real const tmp65 = conn_uul.z.y.z * mdeShiftEpsilon;
			real const tmp66 = conn_uul.x.y.y * mdeShiftEpsilon;
			real const tmp67 = conn_uul.x.y.z * mdeShiftEpsilon;
			real const tmp68 = conn_uul.y.y.z * mdeShiftEpsilon;
			real const tmp69 = e_u.y * mdeShiftEpsilon;
			real const tmp70 = d_luu.x.yy * mdeShiftEpsilon;
			real const tmp71 = d_luu.x.yz * mdeShiftEpsilon;
			real const tmp72 = d_luu.y.yy * mdeShiftEpsilon;
			real const tmp73 = d_luu.y.yz * mdeShiftEpsilon;
			real const tmp74 = d_luu.z.yy * mdeShiftEpsilon;
			real const tmp75 = d_luu.z.yz * mdeShiftEpsilon;
			real const tmp76 = conn_ull.x.yz * tmp71;
			real const tmp77 = conn_ull.y.yz * tmp73;
			real const tmp78 = conn_ull.z.yz * tmp75;
			real const tmp79 = conn_uul.x.x.z * mdeShiftEpsilon;
			real const tmp80 = conn_uul.y.x.z * mdeShiftEpsilon;
			real const tmp81 = conn_uul.z.x.z * mdeShiftEpsilon;
			real const tmp82 = conn_uul.y.y.z * mdeShiftEpsilon;
			real const tmp83 = conn_uul.z.y.z * mdeShiftEpsilon;
			real const tmp84 = conn_uul.z.z.z * mdeShiftEpsilon;
			real const tmp85 = conn_uul.x.y.z * mdeShiftEpsilon;
			real const tmp86 = conn_uul.x.z.z * mdeShiftEpsilon;
			real const tmp87 = conn_uul.y.z.z * mdeShiftEpsilon;
			real const tmp88 = e_u.z * mdeShiftEpsilon;
			real const tmp89 = d_luu.x.zz * mdeShiftEpsilon;
			real const tmp90 = d_luu.y.zz * mdeShiftEpsilon;
			real const tmp91 = d_luu.z.zz * mdeShiftEpsilon;
			(deriv)->beta_u.x += B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y;
			(deriv)->beta_u.y += B_u.y + beta_u.x * b_ul.x.y + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y;
			(deriv)->beta_u.z += B_u.z + beta_u.x * b_ul.x.z + beta_u.z * b_ul.z.z + beta_u.y * b_ul.y.z;
			(deriv)->B_u.x += (6. * B_u.x * b_ul.x.x + 6. * B_u.y * b_ul.x.y + 6. * B_u.z * b_ul.x.z + -3. * d_u.x * tmp1 + -3. * d_u.x * tmp2 + -3. * d_u.x * tmp3 + -3. * d_l.x * b_uu.xx * mdeShiftEpsilon + -3. * d_l.y * tmp4 + -3. * d_l.z * tmp5 + 4. * e_u.x * tmp6 + 6. * d_l.x * DBeta_uu.xx * mdeShiftEpsilon + 6. * d_l.y * tmp7 + 6. * d_l.z * tmp8 + 6. * DBeta_uu.xx * tmp9 + 6. * DBeta_uu.yy * conn_ull.x.yy * mdeShiftEpsilon + 6. * DBeta_uu.zz * conn_ull.x.zz * mdeShiftEpsilon + 6. * b_ul.x.x * tmp10 + 6. * b_ul.x.y * tmp11 + 6. * b_ul.x.z * tmp12 + 6. * b_ul.y.y * tmp13 + 6. * b_ul.y.z * tmp14 + 6. * b_ul.z.z * tmp15 + 6. * b_ul.x.y * tmp16 + 6. * b_ul.x.z * tmp17 + 6. * b_ul.y.z * tmp18 + 12. * DBeta_uu.xy * tmp19 + 12. * DBeta_uu.xz * tmp20 + 12. * DBeta_uu.yz * conn_ull.x.yz * mdeShiftEpsilon + 6. * beta_u.x * d_l.x * tmp21 + 6. * beta_u.x * d_l.x * tmp10 + 6. * beta_u.x * d_l.x * tmp22 + 6. * beta_u.x * d_l.y * tmp11 + 6. * beta_u.x * d_l.y * tmp23 + 6. * beta_u.x * d_l.z * tmp12 + 6. * beta_u.x * d_l.z * tmp24 + -6. * beta_u.x * conn_uul.x.x.x * tmp9 + -6. * beta_u.x * conn_uul.x.x.y * tmp25 + -6. * beta_u.x * conn_uul.x.x.z * tmp26 + -6. * beta_u.x * conn_ull.x.xy * tmp11 + -6. * beta_u.x * conn_ull.x.xz * tmp12 + -6. * beta_u.x * conn_uul.y.x.y * tmp27 + -6. * beta_u.x * conn_uul.y.x.z * tmp28 + -6. * beta_u.x * conn_ull.y.xz * tmp14 + -6. * beta_u.x * conn_uul.z.x.z * tmp29 + 6. * beta_u.y * d_l.x * tmp16 + 6. * beta_u.y * d_l.x * tmp30 + 6. * beta_u.y * d_l.y * tmp21 + 6. * beta_u.y * d_l.y * tmp13 + 6. * beta_u.y * d_l.y * tmp31 + 6. * beta_u.y * d_l.z * tmp14 + 6. * beta_u.y * d_l.z * tmp32 + -6. * beta_u.y * conn_uul.x.x.x * tmp19 + -6. * beta_u.y * conn_uul.x.x.y * tmp27 + -6. * beta_u.y * conn_uul.x.x.z * tmp28 + -6. * beta_u.y * conn_ull.x.yy * tmp11 + -6. * beta_u.y * conn_ull.x.yz * tmp12 + -6. * beta_u.y * conn_uul.y.x.y * tmp33 + -6. * beta_u.y * conn_uul.y.x.z * tmp34 + -6. * beta_u.y * conn_ull.y.yz * tmp14 + -6. * beta_u.y * conn_uul.z.x.z * tmp35 + 6. * beta_u.z * d_l.x * tmp17 + 6. * beta_u.z * d_l.x * tmp36 + 6. * beta_u.z * d_l.y * tmp18 + 6. * beta_u.z * d_l.y * tmp37 + 6. * beta_u.z * d_l.z * tmp21 + 6. * beta_u.z * d_l.z * tmp15 + 6. * beta_u.z * d_l.z * tmp38 + -6. * beta_u.z * conn_uul.x.x.x * tmp20 + -6. * beta_u.z * conn_uul.x.x.y * tmp39 + -6. * beta_u.z * conn_uul.x.x.z * tmp29 + -6. * beta_u.z * conn_ull.x.yz * tmp11 + -6. * beta_u.z * conn_ull.x.zz * tmp12 + -6. * beta_u.z * conn_uul.y.x.y * tmp40 + -6. * beta_u.z * conn_uul.y.x.z * tmp35 + -6. * beta_u.z * conn_ull.y.zz * tmp14 + -6. * beta_u.z * conn_uul.z.x.z * tmp41 + -12. * beta_u.x * conn_ull.x.xx * tmp22 + -12. * beta_u.x * tmp42 + -12. * beta_u.x * tmp43 + -12. * beta_u.x * conn_ull.y.xx * tmp30 + -12. * beta_u.x * tmp44 + -12. * beta_u.x * tmp45 + -12. * beta_u.x * conn_ull.z.xx * tmp36 + -12. * beta_u.x * tmp46 + -12. * beta_u.x * tmp47 + -12. * beta_u.y * conn_ull.x.xy * tmp22 + -12. * beta_u.y * conn_ull.x.yy * tmp23 + -12. * beta_u.y * conn_ull.x.yz * tmp24 + -12. * beta_u.y * conn_ull.y.xy * tmp30 + -12. * beta_u.y * conn_ull.y.yy * tmp31 + -12. * beta_u.y * conn_ull.y.yz * tmp32 + -12. * beta_u.y * conn_ull.z.xy * tmp36 + -12. * beta_u.y * conn_ull.z.yy * tmp37 + -12. * beta_u.y * conn_ull.z.yz * tmp38 + -12. * beta_u.z * conn_ull.x.xz * tmp22 + -12. * beta_u.z * conn_ull.x.yz * tmp23 + -12. * beta_u.z * conn_ull.x.zz * tmp24 + -12. * beta_u.z * conn_ull.y.xz * tmp30 + -12. * beta_u.z * conn_ull.y.yz * tmp31 + -12. * beta_u.z * conn_ull.y.zz * tmp32 + -12. * beta_u.z * conn_ull.z.xz * tmp36 + -12. * beta_u.z * conn_ull.z.yz * tmp37 + -12. * beta_u.z * conn_ull.z.zz * tmp38 + -12. * A_uu.xx * tmp49 + -12. * A_uu.xx * tmp50 + -12. * A_uu.xy * tmp51 + -12. * A_uu.xy * tmp52 + -12. * A_uu.xz * tmp53 + -12. * A_uu.xz * tmp54 + -12. * A_uu.yy * conn_ull.x.yy * tmp48 + -12. * A_uu.zz * conn_ull.x.zz * tmp48 + -36. * A_uu.xy * tmp55 + -36. * A_uu.xz * tmp56 + -24. * A_uu.yz * conn_ull.x.yz * tmp48 + -24. * A_uu.xx * tmp57) / 6.;
			(deriv)->B_u.y += (6. * B_u.x * b_ul.x.y + 6. * B_u.y * b_ul.y.y + 6. * B_u.z * b_ul.y.z + -3. * d_u.y * tmp1 + -3. * d_u.y * tmp2 + -3. * d_u.y * tmp3 + -3. * d_l.x * tmp4 + -3. * d_l.y * b_uu.yy * mdeShiftEpsilon + -3. * d_l.z * tmp58 + 4. * e_u.y * tmp6 + 6. * d_l.x * tmp7 + 6. * d_l.y * DBeta_uu.yy * mdeShiftEpsilon + 6. * d_l.z * tmp59 + 6. * DBeta_uu.xx * tmp25 + 6. * DBeta_uu.yy * tmp33 + 6. * DBeta_uu.zz * conn_ull.y.zz * mdeShiftEpsilon + 6. * b_ul.x.x * tmp60 + 6. * b_ul.x.y * tmp61 + 6. * b_ul.x.z * tmp62 + 6. * b_ul.y.y * tmp63 + 6. * b_ul.y.z * tmp64 + 6. * b_ul.z.z * tmp65 + 6. * b_ul.x.y * tmp66 + 6. * b_ul.x.z * tmp67 + 6. * b_ul.y.z * tmp68 + 12. * DBeta_uu.xy * tmp27 + 12. * DBeta_uu.xz * tmp39 + 12. * DBeta_uu.yz * tmp40 + 6. * beta_u.x * d_l.x * tmp69 + 6. * beta_u.x * d_l.x * tmp60 + 6. * beta_u.x * d_l.x * tmp23 + 6. * beta_u.x * d_l.y * tmp61 + 6. * beta_u.x * d_l.y * tmp70 + 6. * beta_u.x * d_l.z * tmp62 + 6. * beta_u.x * d_l.z * tmp71 + -6. * beta_u.x * conn_uul.x.y.y * tmp25 + -6. * beta_u.x * conn_uul.x.y.z * tmp26 + -6. * beta_u.x * conn_uul.x.x.y * tmp9 + -6. * beta_u.x * conn_ull.x.xy * tmp61 + -6. * beta_u.x * conn_ull.x.xz * tmp62 + -6. * beta_u.x * conn_uul.y.y.y * tmp27 + -6. * beta_u.x * conn_uul.y.y.z * tmp28 + -6. * beta_u.x * conn_ull.y.xz * tmp64 + -6. * beta_u.x * conn_uul.z.y.z * tmp29 + 6. * beta_u.y * d_l.x * tmp66 + 6. * beta_u.y * d_l.x * tmp31 + 6. * beta_u.y * d_l.y * tmp69 + 6. * beta_u.y * d_l.y * tmp63 + 6. * beta_u.y * d_l.y * tmp72 + 6. * beta_u.y * d_l.z * tmp64 + 6. * beta_u.y * d_l.z * tmp73 + -6. * beta_u.y * conn_uul.x.y.y * tmp27 + -6. * beta_u.y * conn_uul.x.y.z * tmp28 + -6. * beta_u.y * conn_uul.x.x.y * tmp19 + -6. * beta_u.y * conn_ull.x.yy * tmp61 + -6. * beta_u.y * conn_ull.x.yz * tmp62 + -6. * beta_u.y * conn_uul.y.y.y * tmp33 + -6. * beta_u.y * conn_uul.y.y.z * tmp34 + -6. * beta_u.y * conn_ull.y.yz * tmp64 + -6. * beta_u.y * conn_uul.z.y.z * tmp35 + 6. * beta_u.z * d_l.x * tmp67 + 6. * beta_u.z * d_l.x * tmp37 + 6. * beta_u.z * d_l.y * tmp68 + 6. * beta_u.z * d_l.y * tmp74 + 6. * beta_u.z * d_l.z * tmp69 + 6. * beta_u.z * d_l.z * tmp65 + 6. * beta_u.z * d_l.z * tmp75 + -6. * beta_u.z * conn_uul.x.y.y * tmp39 + -6. * beta_u.z * conn_uul.x.y.z * tmp29 + -6. * beta_u.z * conn_uul.x.x.y * tmp20 + -6. * beta_u.z * conn_ull.x.yz * tmp61 + -6. * beta_u.z * conn_ull.x.zz * tmp62 + -6. * beta_u.z * conn_uul.y.y.y * tmp40 + -6. * beta_u.z * conn_uul.y.y.z * tmp35 + -6. * beta_u.z * conn_ull.y.zz * tmp64 + -6. * beta_u.z * conn_uul.z.y.z * tmp41 + -12. * beta_u.x * conn_ull.x.xx * tmp23 + -12. * beta_u.x * conn_ull.x.xy * tmp70 + -12. * beta_u.x * conn_ull.x.xz * tmp71 + -12. * beta_u.x * conn_ull.y.xx * tmp31 + -12. * beta_u.x * conn_ull.y.xy * tmp72 + -12. * beta_u.x * conn_ull.y.xz * tmp73 + -12. * beta_u.x * conn_ull.z.xx * tmp37 + -12. * beta_u.x * conn_ull.z.xy * tmp74 + -12. * beta_u.x * conn_ull.z.xz * tmp75 + -12. * beta_u.y * tmp42 + -12. * beta_u.y * conn_ull.x.yy * tmp70 + -12. * beta_u.y * tmp76 + -12. * beta_u.y * tmp44 + -12. * beta_u.y * conn_ull.y.yy * tmp72 + -12. * beta_u.y * tmp77 + -12. * beta_u.y * tmp46 + -12. * beta_u.y * conn_ull.z.yy * tmp74 + -12. * beta_u.y * tmp78 + -12. * beta_u.z * conn_ull.x.xz * tmp23 + -12. * beta_u.z * conn_ull.x.yz * tmp70 + -12. * beta_u.z * conn_ull.x.zz * tmp71 + -12. * beta_u.z * conn_ull.y.xz * tmp31 + -12. * beta_u.z * conn_ull.y.yz * tmp72 + -12. * beta_u.z * conn_ull.y.zz * tmp73 + -12. * beta_u.z * conn_ull.z.xz * tmp37 + -12. * beta_u.z * conn_ull.z.yz * tmp74 + -12. * beta_u.z * conn_ull.z.zz * tmp75 + -12. * A_uu.xx * conn_ull.y.xx * tmp48 + -12. * A_uu.xy * tmp57 + -12. * A_uu.xy * tmp50 + -12. * A_uu.yy * tmp55 + -12. * A_uu.yy * tmp52 + -12. * A_uu.yz * tmp56 + -12. * A_uu.yz * tmp54 + -12. * A_uu.zz * conn_ull.y.zz * tmp48 + -36. * A_uu.xy * tmp49 + -36. * A_uu.yz * tmp53 + -24. * A_uu.yy * tmp51 + -24. * A_uu.xz * conn_ull.y.xz * tmp48) / 6.;
			(deriv)->B_u.z += (6. * B_u.x * b_ul.x.z + 6. * B_u.y * b_ul.y.z + 6. * B_u.z * b_ul.z.z + -3. * d_u.z * tmp1 + -3. * d_u.z * tmp2 + -3. * d_u.z * tmp3 + -3. * d_l.x * tmp5 + -3. * d_l.y * tmp58 + -3. * d_l.z * b_uu.zz * mdeShiftEpsilon + 4. * e_u.z * tmp6 + 6. * d_l.x * tmp8 + 6. * d_l.y * tmp59 + 6. * d_l.z * DBeta_uu.zz * mdeShiftEpsilon + 6. * DBeta_uu.xx * tmp26 + 6. * DBeta_uu.yy * tmp34 + 6. * DBeta_uu.zz * tmp41 + 6. * b_ul.x.x * tmp79 + 6. * b_ul.x.y * tmp80 + 6. * b_ul.x.z * tmp81 + 6. * b_ul.y.y * tmp82 + 6. * b_ul.y.z * tmp83 + 6. * b_ul.z.z * tmp84 + 6. * b_ul.x.y * tmp85 + 6. * b_ul.x.z * tmp86 + 6. * b_ul.y.z * tmp87 + 12. * DBeta_uu.xy * tmp28 + 12. * DBeta_uu.xz * tmp29 + 12. * DBeta_uu.yz * tmp35 + 6. * beta_u.x * d_l.x * tmp88 + 6. * beta_u.x * d_l.x * tmp79 + 6. * beta_u.x * d_l.x * tmp24 + 6. * beta_u.x * d_l.y * tmp80 + 6. * beta_u.x * d_l.y * tmp71 + 6. * beta_u.x * d_l.z * tmp81 + 6. * beta_u.x * d_l.z * tmp89 + -6. * beta_u.x * conn_uul.x.z.z * tmp26 + -6. * beta_u.x * conn_uul.x.x.z * tmp9 + -6. * beta_u.x * conn_ull.x.xy * tmp80 + -6. * beta_u.x * conn_ull.x.xz * tmp81 + -6. * beta_u.x * conn_uul.x.y.z * tmp25 + -6. * beta_u.x * conn_uul.y.z.z * tmp28 + -6. * beta_u.x * conn_ull.y.xy * tmp82 + -6. * beta_u.x * conn_ull.y.xz * tmp83 + -6. * beta_u.x * conn_uul.z.z.z * tmp29 + 6. * beta_u.y * d_l.x * tmp85 + 6. * beta_u.y * d_l.x * tmp32 + 6. * beta_u.y * d_l.y * tmp88 + 6. * beta_u.y * d_l.y * tmp82 + 6. * beta_u.y * d_l.y * tmp73 + 6. * beta_u.y * d_l.z * tmp83 + 6. * beta_u.y * d_l.z * tmp90 + -6. * beta_u.y * conn_uul.x.z.z * tmp28 + -6. * beta_u.y * conn_uul.x.x.z * tmp19 + -6. * beta_u.y * conn_uul.x.y.z * tmp27 + -6. * beta_u.y * conn_ull.x.yy * tmp80 + -6. * beta_u.y * conn_ull.x.yz * tmp81 + -6. * beta_u.y * conn_uul.y.z.z * tmp34 + -6. * beta_u.y * conn_uul.y.y.z * tmp33 + -6. * beta_u.y * conn_ull.y.yz * tmp83 + -6. * beta_u.y * conn_uul.z.z.z * tmp35 + 6. * beta_u.z * d_l.x * tmp86 + 6. * beta_u.z * d_l.x * tmp38 + 6. * beta_u.z * d_l.y * tmp87 + 6. * beta_u.z * d_l.y * tmp75 + 6. * beta_u.z * d_l.z * tmp88 + 6. * beta_u.z * d_l.z * tmp84 + 6. * beta_u.z * d_l.z * tmp91 + -6. * beta_u.z * conn_uul.x.z.z * tmp29 + -6. * beta_u.z * conn_uul.x.x.z * tmp20 + -6. * beta_u.z * conn_uul.x.y.z * tmp39 + -6. * beta_u.z * conn_ull.x.yz * tmp80 + -6. * beta_u.z * conn_ull.x.zz * tmp81 + -6. * beta_u.z * conn_uul.y.z.z * tmp35 + -6. * beta_u.z * conn_uul.y.y.z * tmp40 + -6. * beta_u.z * conn_ull.y.zz * tmp83 + -6. * beta_u.z * conn_uul.z.z.z * tmp41 + -12. * beta_u.x * conn_ull.x.xx * tmp24 + -12. * beta_u.x * conn_ull.x.xy * tmp71 + -12. * beta_u.x * conn_ull.x.xz * tmp89 + -12. * beta_u.x * conn_ull.y.xx * tmp32 + -12. * beta_u.x * conn_ull.y.xy * tmp73 + -12. * beta_u.x * conn_ull.y.xz * tmp90 + -12. * beta_u.x * conn_ull.z.xx * tmp38 + -12. * beta_u.x * conn_ull.z.xy * tmp75 + -12. * beta_u.x * conn_ull.z.xz * tmp91 + -12. * beta_u.y * conn_ull.x.xy * tmp24 + -12. * beta_u.y * conn_ull.x.yy * tmp71 + -12. * beta_u.y * conn_ull.x.yz * tmp89 + -12. * beta_u.y * conn_ull.y.xy * tmp32 + -12. * beta_u.y * conn_ull.y.yy * tmp73 + -12. * beta_u.y * conn_ull.y.yz * tmp90 + -12. * beta_u.y * conn_ull.z.xy * tmp38 + -12. * beta_u.y * conn_ull.z.yy * tmp75 + -12. * beta_u.y * conn_ull.z.yz * tmp91 + -12. * beta_u.z * tmp43 + -12. * beta_u.z * tmp76 + -12. * beta_u.z * conn_ull.x.zz * tmp89 + -12. * beta_u.z * tmp45 + -12. * beta_u.z * tmp77 + -12. * beta_u.z * conn_ull.y.zz * tmp90 + -12. * beta_u.z * tmp47 + -12. * beta_u.z * tmp78 + -12. * beta_u.z * conn_ull.z.zz * tmp91 + -12. * A_uu.xx * conn_ull.z.xx * tmp48 + -12. * A_uu.xz * tmp57 + -12. * A_uu.xz * tmp49 + -12. * A_uu.yy * conn_ull.z.yy * tmp48 + -12. * A_uu.yz * tmp55 + -12. * A_uu.yz * tmp51 + -12. * A_uu.zz * tmp56 + -12. * A_uu.zz * tmp53 + -36. * A_uu.xz * tmp50 + -36. * A_uu.yz * tmp52 + -24. * A_uu.zz * tmp54 + -24. * A_uu.xy * conn_ull.z.xy * tmp48) / 6.;
		}
		<? end ?>/* eqn.useShift == "MinimalDistortionHyperbolic" */
		<? if eqn.useShift == "GammaDriverHyperbolic" then ?>
		{
			real const tmp1 = invW * invW;
			real const tmp2 = tr_DHatBeta * tmp1;
			real const tmp3 = connHat_uul.x.y.y * tmp1;
			real const tmp4 = connHat_uul.x.z.z * tmp1;
			real const tmp5 = connHat_uul.y.x.x * tmp1;
			real const tmp6 = connHat_uul.y.y.y * tmp1;
			real const tmp7 = connHat_uul.y.z.z * tmp1;
			real const tmp8 = connHat_uul.z.x.x * tmp1;
			real const tmp9 = connHat_uul.z.y.y * tmp1;
			real const tmp10 = connHat_uul.z.z.z * tmp1;
			real const tmp11 = DHatBeta_ul.x.x * tmp1;
			real const tmp12 = DHatBeta_ul.x.y * tmp1;
			real const tmp13 = DHatBeta_ul.x.z * tmp1;
			real const tmp14 = A_uu.xy * tmp1;
			real const tmp15 = A_uu.xz * tmp1;
			real const tmp16 = connHat_ull.y.xx * tmp1;
			real const tmp17 = connHat_ull.z.xx * tmp1;
			real const tmp18 = connHat_ull.x.xx * tmp1;
			real const tmp19 = connHat_ull.y.xy * tmp1;
			real const tmp20 = connHat_uul.x.y.y * tmp19;
			real const tmp21 = connHat_ull.z.xy * tmp1;
			real const tmp22 = connHat_ull.z.xz * tmp1;
			real const tmp23 = connHat_uul.x.z.z * tmp22;
			real const tmp24 = connHat_ull.x.xy * tmp1;
			real const tmp25 = connHat_ull.x.xz * tmp1;
			real const tmp26 = connHat_ull.x.xy * tmp5;
			real const tmp27 = connHat_ull.x.xz * tmp8;
			real const tmp28 = connHat_ull.y.xz * tmp1;
			real const tmp29 = connHat_ull.y.yy * tmp1;
			real const tmp30 = connHat_ull.z.yy * tmp1;
			real const tmp31 = connHat_ull.z.yz * tmp1;
			real const tmp32 = connHat_ull.x.yz * tmp1;
			real const tmp33 = connHat_ull.y.yz * tmp1;
			real const tmp34 = connHat_ull.z.zz * tmp1;
			real const tmp35 = connHat_ull.y.zz * tmp1;
			real const tmp36 = tr_K * tmp1;
			real const tmp37 = alpha * tmp36;
			real const tmp38 = alpha * tmp1;
			real const tmp39 = A_uu.xy * tmp38;
			real const tmp40 = A_uu.xz * tmp38;
			real const tmp41 = connHat_uul.y.x.y * tmp1;
			real const tmp42 = connHat_uul.y.x.z * tmp1;
			real const tmp43 = connHat_uul.y.y.z * tmp1;
			real const tmp44 = connHat_uul.x.x.x * tmp1;
			real const tmp45 = DHatBeta_ul.y.x * tmp1;
			real const tmp46 = DHatBeta_ul.y.y * tmp1;
			real const tmp47 = DHatBeta_ul.y.z * tmp1;
			real const tmp48 = A_uu.yz * tmp1;
			real const tmp49 = connHat_uul.y.z.z * tmp31;
			real const tmp50 = connHat_ull.y.yz * tmp9;
			real const tmp51 = A_uu.yz * tmp38;
			real const tmp52 = connHat_uul.z.x.y * tmp1;
			real const tmp53 = connHat_uul.z.x.z * tmp1;
			real const tmp54 = connHat_uul.z.y.z * tmp1;
			real const tmp55 = connHat_uul.z.x.y * tmp1;
			real const tmp56 = DHatBeta_ul.z.x * tmp1;
			real const tmp57 = DHatBeta_ul.z.y * tmp1;
			real const tmp58 = DHatBeta_ul.z.z * tmp1;
			(deriv)->beta_u.x += B_u.x + beta_u.x * b_ul.x.x + beta_u.z * b_ul.x.z + beta_u.y * b_ul.x.y;
			(deriv)->beta_u.y += B_u.y + beta_u.x * b_ul.x.y + beta_u.z * b_ul.y.z + beta_u.y * b_ul.y.y;
			(deriv)->beta_u.z += B_u.z + beta_u.x * b_ul.x.z + beta_u.z * b_ul.z.z + beta_u.y * b_ul.y.z;
			(deriv)->B_u.x += (-9. * LambdaBar_u.x * b_ul.x.x + -9. * LambdaBar_u.y * b_ul.x.y + -9. * LambdaBar_u.z * b_ul.x.z + -2. * dDelta_u.x * tmp2 + 12. * B_u.x * b_ul.x.x + -12. * B_u.x * gammaDriver_eta + 12. * B_u.y * b_ul.x.y + 12. * B_u.z * b_ul.x.z + -9. * b_ul.x.x * tmp3 + -9. * b_ul.x.x * tmp4 + 9. * b_ul.x.y * connHat_uul.x.x.y * tmp1 + -9. * b_ul.x.y * tmp5 + -9. * b_ul.x.y * tmp6 + -9. * b_ul.x.y * tmp7 + 9. * b_ul.x.z * connHat_uul.x.x.z * tmp1 + -9. * b_ul.x.z * tmp8 + -9. * b_ul.x.z * tmp9 + -9. * b_ul.x.z * tmp10 + 9. * b_ul.y.y * tmp3 + 9. * b_ul.y.z * connHat_uul.x.y.z * tmp1 + 9. * b_ul.z.z * tmp4 + 9. * b_ul.x.y * connHat_uul.x.x.y * tmp1 + 9. * b_ul.x.z * connHat_uul.x.x.z * tmp1 + 9. * b_ul.y.z * connHat_uul.x.y.z * tmp1 + -6. * dDelta_u.x * tmp11 + -6. * dDelta_u.y * tmp12 + -6. * dDelta_u.z * tmp13 + 6. * e_u.x * tmp2 + 6. * DeltaGamma_uul.x.x.x * tmp2 + 6. * DeltaGamma_uul.x.y.y * tmp2 + 6. * DeltaGamma_uul.x.z.z * tmp2 + 18. * dDelta_l.x * A_uu.xx * tmp1 + 18. * dDelta_l.y * tmp14 + 18. * dDelta_l.z * tmp15 + 18. * e_u.x * tmp11 + 18. * e_u.y * tmp12 + 18. * e_u.z * tmp13 + 18. * A_uu.xx * DeltaGamma_ull.x.xx * tmp1 + 18. * A_uu.yy * DeltaGamma_ull.x.yy * tmp1 + 18. * A_uu.zz * DeltaGamma_ull.x.zz * tmp1 + 36. * A_uu.xy * DeltaGamma_ull.x.xy * tmp1 + 36. * A_uu.xz * DeltaGamma_ull.x.xz * tmp1 + 36. * A_uu.yz * DeltaGamma_ull.x.yz * tmp1 + 9. * beta_u.x * connHat_uul.x.x.y * tmp16 + 9. * beta_u.x * connHat_uul.x.x.z * tmp17 + -9. * beta_u.x * connHat_uul.x.y.y * tmp18 + 9. * beta_u.x * tmp20 + 9. * beta_u.x * connHat_uul.x.y.z * tmp21 + -9. * beta_u.x * connHat_uul.x.z.z * tmp18 + 9. * beta_u.x * tmp23 + 9. * beta_u.x * connHat_uul.x.x.y * tmp24 + 9. * beta_u.x * connHat_uul.x.x.z * tmp25 + -9. * beta_u.x * tmp26 + -9. * beta_u.x * connHat_ull.x.xy * tmp6 + -9. * beta_u.x * connHat_ull.x.xy * tmp7 + -9. * beta_u.x * tmp27 + -9. * beta_u.x * connHat_ull.x.xz * tmp9 + -9. * beta_u.x * connHat_ull.x.xz * tmp10 + 9. * beta_u.x * connHat_uul.x.y.z * tmp28 + 9. * beta_u.y * connHat_uul.x.x.y * tmp19 + 9. * beta_u.y * connHat_uul.x.x.z * tmp21 + -9. * beta_u.y * connHat_uul.x.y.y * tmp24 + 9. * beta_u.y * connHat_uul.x.y.y * tmp29 + 9. * beta_u.y * connHat_uul.x.y.z * tmp30 + -9. * beta_u.y * connHat_uul.x.z.z * tmp24 + 9. * beta_u.y * connHat_uul.x.z.z * tmp31 + 9. * beta_u.y * connHat_uul.x.x.y * connHat_ull.x.yy * tmp1 + 9. * beta_u.y * connHat_uul.x.x.z * tmp32 + 9. * beta_u.y * connHat_uul.x.y.z * tmp33 + -9. * beta_u.y * connHat_ull.x.yy * tmp5 + -9. * beta_u.y * connHat_ull.x.yy * tmp6 + -9. * beta_u.y * connHat_ull.x.yy * tmp7 + -9. * beta_u.y * connHat_ull.x.yz * tmp8 + -9. * beta_u.y * connHat_ull.x.yz * tmp9 + -9. * beta_u.y * connHat_ull.x.yz * tmp10 + 9. * beta_u.z * connHat_uul.x.x.y * tmp28 + 9. * beta_u.z * connHat_uul.x.x.z * tmp22 + -9. * beta_u.z * connHat_uul.x.y.y * tmp25 + 9. * beta_u.z * connHat_uul.x.y.y * tmp33 + 9. * beta_u.z * connHat_uul.x.y.z * tmp31 + -9. * beta_u.z * connHat_uul.x.z.z * tmp25 + 9. * beta_u.z * connHat_uul.x.z.z * tmp34 + 9. * beta_u.z * connHat_uul.x.x.y * tmp32 + 9. * beta_u.z * connHat_uul.x.x.z * connHat_ull.x.zz * tmp1 + 9. * beta_u.z * connHat_uul.x.y.z * tmp35 + -9. * beta_u.z * connHat_ull.x.yz * tmp5 + -9. * beta_u.z * connHat_ull.x.yz * tmp6 + -9. * beta_u.z * connHat_ull.x.yz * tmp7 + -9. * beta_u.z * connHat_ull.x.zz * tmp8 + -9. * beta_u.z * connHat_ull.x.zz * tmp9 + -9. * beta_u.z * connHat_ull.x.zz * tmp10 + 8. * dDelta_u.x * tmp37 + -18. * a_l.x * A_uu.xx * tmp38 + -18. * a_l.y * tmp39 + -18. * a_l.z * tmp40 + -24. * e_u.x * tmp37 + 12. * a_u.x * tmp37) / 12.;
			(deriv)->B_u.y += (-9. * LambdaBar_u.x * b_ul.x.y + -9. * LambdaBar_u.y * b_ul.y.y + -9. * LambdaBar_u.z * b_ul.y.z + -2. * dDelta_u.y * tmp2 + 12. * B_u.x * b_ul.x.y + 12. * B_u.y * b_ul.y.y + -12. * B_u.y * gammaDriver_eta + 12. * B_u.z * b_ul.y.z + 9. * b_ul.x.x * tmp5 + 9. * b_ul.x.y * tmp41 + 9. * b_ul.x.z * tmp42 + -9. * b_ul.y.y * tmp5 + -9. * b_ul.y.y * tmp7 + 9. * b_ul.y.z * tmp43 + -9. * b_ul.y.z * tmp8 + -9. * b_ul.y.z * tmp9 + -9. * b_ul.y.z * tmp10 + 9. * b_ul.z.z * tmp7 + -9. * b_ul.x.y * tmp44 + -9. * b_ul.x.y * tmp3 + -9. * b_ul.x.y * tmp4 + 9. * b_ul.x.y * connHat_uul.y.x.y * tmp1 + 9. * b_ul.x.z * connHat_uul.y.x.z * tmp1 + 9. * b_ul.y.z * connHat_uul.y.y.z * tmp1 + -6. * dDelta_u.x * tmp45 + -6. * dDelta_u.y * tmp46 + -6. * dDelta_u.z * tmp47 + 6. * e_u.y * tmp2 + 6. * DeltaGamma_uul.y.x.x * tmp2 + 6. * DeltaGamma_uul.y.y.y * tmp2 + 6. * DeltaGamma_uul.y.z.z * tmp2 + 18. * dDelta_l.x * tmp14 + 18. * dDelta_l.y * A_uu.yy * tmp1 + 18. * dDelta_l.z * tmp48 + 18. * e_u.x * tmp45 + 18. * e_u.y * tmp46 + 18. * e_u.z * tmp47 + 18. * A_uu.xx * DeltaGamma_ull.y.xx * tmp1 + 18. * A_uu.yy * DeltaGamma_ull.y.yy * tmp1 + 18. * A_uu.zz * DeltaGamma_ull.y.zz * tmp1 + 36. * A_uu.xy * DeltaGamma_ull.y.xy * tmp1 + 36. * A_uu.xz * DeltaGamma_ull.y.xz * tmp1 + 36. * A_uu.yz * DeltaGamma_ull.y.yz * tmp1 + -9. * beta_u.x * connHat_uul.x.x.x * tmp16 + -9. * beta_u.x * connHat_uul.x.y.y * tmp16 + -9. * beta_u.x * connHat_uul.x.z.z * tmp16 + 9. * beta_u.x * connHat_ull.x.xx * tmp5 + 9. * beta_u.x * connHat_ull.x.xy * tmp41 + 9. * beta_u.x * connHat_ull.x.xz * tmp42 + -9. * beta_u.x * connHat_uul.y.x.x * tmp19 + 9. * beta_u.x * connHat_uul.y.x.y * tmp16 + 9. * beta_u.x * connHat_uul.y.x.z * tmp17 + 9. * beta_u.x * connHat_uul.y.y.z * tmp21 + -9. * beta_u.x * connHat_uul.y.z.z * tmp19 + 9. * beta_u.x * connHat_uul.y.z.z * tmp22 + 9. * beta_u.x * connHat_ull.y.xz * tmp43 + -9. * beta_u.x * connHat_ull.y.xz * tmp8 + -9. * beta_u.x * connHat_ull.y.xz * tmp9 + -9. * beta_u.x * connHat_ull.y.xz * tmp10 + -9. * beta_u.y * connHat_uul.x.x.x * tmp19 + -9. * beta_u.y * tmp20 + -9. * beta_u.y * connHat_uul.x.z.z * tmp19 + 9. * beta_u.y * tmp26 + 9. * beta_u.y * connHat_ull.x.yy * tmp41 + 9. * beta_u.y * connHat_ull.x.yz * tmp42 + -9. * beta_u.y * connHat_uul.y.x.x * tmp29 + 9. * beta_u.y * connHat_uul.y.x.y * tmp19 + 9. * beta_u.y * connHat_uul.y.x.z * tmp21 + 9. * beta_u.y * connHat_uul.y.y.z * tmp30 + -9. * beta_u.y * connHat_uul.y.z.z * tmp29 + 9. * beta_u.y * tmp49 + 9. * beta_u.y * connHat_uul.y.y.z * tmp33 + -9. * beta_u.y * connHat_ull.y.yz * tmp8 + -9. * beta_u.y * tmp50 + -9. * beta_u.y * connHat_ull.y.yz * tmp10 + -9. * beta_u.z * connHat_uul.x.x.x * tmp28 + -9. * beta_u.z * connHat_uul.x.y.y * tmp28 + -9. * beta_u.z * connHat_uul.x.z.z * tmp28 + 9. * beta_u.z * connHat_ull.x.xz * tmp5 + 9. * beta_u.z * connHat_ull.x.yz * tmp41 + 9. * beta_u.z * connHat_ull.x.zz * tmp42 + -9. * beta_u.z * connHat_uul.y.x.x * tmp33 + 9. * beta_u.z * connHat_uul.y.x.y * tmp28 + 9. * beta_u.z * connHat_uul.y.x.z * tmp22 + 9. * beta_u.z * connHat_uul.y.y.z * tmp31 + -9. * beta_u.z * connHat_uul.y.z.z * tmp33 + 9. * beta_u.z * connHat_uul.y.z.z * tmp34 + 9. * beta_u.z * connHat_uul.y.y.z * tmp35 + -9. * beta_u.z * connHat_ull.y.zz * tmp8 + -9. * beta_u.z * connHat_ull.y.zz * tmp9 + -9. * beta_u.z * connHat_ull.y.zz * tmp10 + 8. * dDelta_u.y * tmp37 + -18. * a_l.x * tmp39 + -18. * a_l.y * A_uu.yy * tmp38 + -18. * a_l.z * tmp51 + -24. * e_u.y * tmp37 + 12. * a_u.y * tmp37) / 12.;
			(deriv)->B_u.z += (-9. * LambdaBar_u.x * b_ul.x.z + -9. * LambdaBar_u.y * b_ul.y.z + -9. * LambdaBar_u.z * b_ul.z.z + -2. * dDelta_u.z * tmp2 + 12. * B_u.x * b_ul.x.z + 12. * B_u.y * b_ul.y.z + 12. * B_u.z * b_ul.z.z + -12. * B_u.z * gammaDriver_eta + 9. * b_ul.x.x * tmp8 + 9. * b_ul.x.y * tmp52 + 9. * b_ul.x.z * tmp53 + 9. * b_ul.y.y * tmp9 + 9. * b_ul.y.z * tmp54 + -9. * b_ul.z.z * tmp8 + -9. * b_ul.z.z * tmp9 + 9. * b_ul.x.y * tmp55 + -9. * b_ul.x.z * tmp44 + -9. * b_ul.x.z * tmp3 + -9. * b_ul.x.z * tmp4 + 9. * b_ul.x.z * connHat_uul.z.x.z * tmp1 + -9. * b_ul.y.z * tmp5 + -9. * b_ul.y.z * tmp6 + -9. * b_ul.y.z * tmp7 + 9. * b_ul.y.z * connHat_uul.z.y.z * tmp1 + -6. * dDelta_u.x * tmp56 + -6. * dDelta_u.y * tmp57 + -6. * dDelta_u.z * tmp58 + 6. * e_u.z * tmp2 + 6. * DeltaGamma_uul.z.x.x * tmp2 + 6. * DeltaGamma_uul.z.y.y * tmp2 + 6. * DeltaGamma_uul.z.z.z * tmp2 + 18. * dDelta_l.x * tmp15 + 18. * dDelta_l.y * tmp48 + 18. * dDelta_l.z * A_uu.zz * tmp1 + 18. * e_u.x * tmp56 + 18. * e_u.y * tmp57 + 18. * e_u.z * tmp58 + 18. * A_uu.xx * DeltaGamma_ull.z.xx * tmp1 + 18. * A_uu.yy * DeltaGamma_ull.z.yy * tmp1 + 18. * A_uu.zz * DeltaGamma_ull.z.zz * tmp1 + 36. * A_uu.xy * DeltaGamma_ull.z.xy * tmp1 + 36. * A_uu.xz * DeltaGamma_ull.z.xz * tmp1 + 36. * A_uu.yz * DeltaGamma_ull.z.yz * tmp1 + -9. * beta_u.x * connHat_uul.x.x.x * tmp17 + -9. * beta_u.x * connHat_uul.x.y.y * tmp17 + -9. * beta_u.x * connHat_uul.x.z.z * tmp17 + 9. * beta_u.x * connHat_ull.x.xx * tmp8 + 9. * beta_u.x * connHat_ull.x.xy * tmp52 + 9. * beta_u.x * connHat_ull.x.xz * tmp53 + -9. * beta_u.x * connHat_uul.y.x.x * tmp21 + -9. * beta_u.x * connHat_uul.y.y.y * tmp21 + -9. * beta_u.x * connHat_uul.y.z.z * tmp21 + 9. * beta_u.x * connHat_ull.y.xx * tmp55 + 9. * beta_u.x * connHat_ull.y.xy * tmp9 + 9. * beta_u.x * connHat_ull.y.xz * tmp54 + -9. * beta_u.x * connHat_uul.z.x.x * tmp22 + 9. * beta_u.x * connHat_uul.z.x.z * tmp17 + -9. * beta_u.x * connHat_uul.z.y.y * tmp22 + 9. * beta_u.x * connHat_uul.z.y.z * tmp21 + -9. * beta_u.y * connHat_uul.x.x.x * tmp21 + -9. * beta_u.y * connHat_uul.x.y.y * tmp21 + -9. * beta_u.y * connHat_uul.x.z.z * tmp21 + 9. * beta_u.y * connHat_ull.x.xy * tmp8 + 9. * beta_u.y * connHat_ull.x.yy * tmp52 + 9. * beta_u.y * connHat_ull.x.yz * tmp53 + -9. * beta_u.y * connHat_uul.y.x.x * tmp30 + -9. * beta_u.y * connHat_uul.y.y.y * tmp30 + -9. * beta_u.y * connHat_uul.y.z.z * tmp30 + 9. * beta_u.y * connHat_ull.y.xy * tmp55 + 9. * beta_u.y * connHat_ull.y.yy * tmp9 + 9. * beta_u.y * connHat_ull.y.yz * tmp54 + -9. * beta_u.y * connHat_uul.z.x.x * tmp31 + 9. * beta_u.y * connHat_uul.z.x.z * tmp21 + -9. * beta_u.y * connHat_uul.z.y.y * tmp31 + 9. * beta_u.y * connHat_uul.z.y.z * tmp30 + -9. * beta_u.z * connHat_uul.x.x.x * tmp22 + -9. * beta_u.z * connHat_uul.x.y.y * tmp22 + -9. * beta_u.z * tmp23 + 9. * beta_u.z * tmp27 + 9. * beta_u.z * connHat_ull.x.yz * tmp52 + 9. * beta_u.z * connHat_ull.x.zz * tmp53 + -9. * beta_u.z * connHat_uul.y.x.x * tmp31 + -9. * beta_u.z * connHat_uul.y.y.y * tmp31 + -9. * beta_u.z * tmp49 + 9. * beta_u.z * connHat_ull.y.xz * tmp55 + 9. * beta_u.z * tmp50 + 9. * beta_u.z * connHat_ull.y.zz * tmp54 + -9. * beta_u.z * connHat_uul.z.x.x * tmp34 + 9. * beta_u.z * connHat_uul.z.x.z * tmp22 + -9. * beta_u.z * connHat_uul.z.y.y * tmp34 + 9. * beta_u.z * connHat_uul.z.y.z * tmp31 + 8. * dDelta_u.z * tmp37 + -18. * a_l.x * tmp40 + -18. * a_l.y * tmp51 + -18. * a_l.z * A_uu.zz * tmp38 + -24. * e_u.z * tmp37 + 12. * a_u.z * tmp37) / 12.;
		}
		<? end ?>/* eqn.useShift == "GammaDriverHyperbolic" */
	}
	<? end ?>/* eqn.useShift ~= "none" */
	// END CUT from symmath/tests/output/Z4.html
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
	//TODO should Theta or Z be included into these?
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
