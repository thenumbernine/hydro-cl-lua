//// MODULE_NAME: <?=elecChargeMassRatio?>

// r_e = q_e / m_e
// r_e = q_e / (m_i / (m_i / m_e))
// using m = m_i / m_e
// r_e = m q_e / m_i
// using q_i = q_e
// r_e = m q_i / m_i
// using r_i = q_i / m_i
// r_e = m r_i
// https://en.wikipedia.org/wiki/Mass-to-charge_ratio
// q_e / m_e = -1.758820024e+11 C/kg
// notice this hasn't been converted to units yet, so divide by unit_C_per_kg
#define elecChargeMassRatio			(solver->ionElectronMassRatio * solver->ionChargeMassRatio)

//// MODULE_NAME: <?=sqrt_2_and_1_2?>

#define sqrt_1_2 <?=("%.50f"):format(math.sqrt(.5))?>
#define sqrt_2 <?=("%.50f"):format(math.sqrt(2))?>

//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: units <?=coordLenSq?> <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define /*real3*/ calc_EField(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
)\
	(real3_real_mul(\
		(U)->D,\
		1. / (/*eps = */solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3)\
	))
 
#define /*real3*/ calc_HField(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
)\
	(real3_real_mul(\
		(U)->B,\
		1. / (/*mu = */solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2)\
	))

#define /*real3*/ calc_SField(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
) \
	(real3_cross(\
		calc_EField(solver, U),\
		calc_HField(solver, U)))

static inline real calc_H(
	constant <?=solver_t?> const * const solver,
	real const P
) {
	return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.));
}

static inline real calc_h(
	constant <?=solver_t?> const * const solver,
	real const rho,
	real const P
) {
	return calc_H(solver, P) / rho;
}

static inline real calc_HTotal(
	real const P,
	real const ETotal
) {
	return P + ETotal;
}

static inline real calc_hTotal(
	constant <?=solver_t?> const * const solver,
	real const rho,
	real const P,
	real const ETotal
) {
	return calc_HTotal(P, ETotal) / rho;
}

static inline real calc_rho_from_U(
	global <?=cons_t?> const * const U
) {
	return 0.<? 
for _,fluid in ipairs(fluids) do 
?> + (U)-><?=fluid?>_rho<? 
end 
?>;
}

static inline real calc_rho_from_W(
	<?=prim_t?> const * const W
) {
	return 0.<?
for _,fluid in ipairs(fluids) do 
?> + (W)-><?=fluid?>_rho<?
end 
?>;
}

static inline real calc_EPot(
	global <?=cons_t?> const * const U
) {
	return calc_rho_from_U(U) * (U)->ePot;
}

static inline real calc_EPot_from_W(
	<?=prim_t?> const * const W
) {
	return calc_rho_from_W(W) * (W)->ePot;
}

<? for _,fluid in ipairs(fluids) do ?>

static inline real calc_<?=fluid?>_eKin(
	<?=prim_t?> const * const W,
	real3 const pt
) {
	return .5 * coordLenSq((W)-><?=fluid?>_v, pt);
}

static inline real calc_<?=fluid?>_EKin(
	<?=prim_t?> const * const W,
	real3 const pt
) {
	return (W)-><?=fluid?>_rho * calc_<?=fluid?>_eKin(W, pt);
}

static inline real calc_<?=fluid?>_EInt(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const * const W
) {
	return (W)-><?=fluid?>_P / (solver->heatCapacityRatio - 1.);
}

static inline real calc_<?=fluid?>_eInt(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const * const W
) {
	return calc_<?=fluid?>_EInt(solver, W) / (W)-><?=fluid?>_rho;
}

#define /*real*/ calc_<?=fluid?>_EKin_fromCons(\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
)	(.5 * coordLenSq((U)-><?=fluid?>_m, pt) / (U)-><?=fluid?>_rho)

static inline real calc_<?=fluid?>_ETotal(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const * const W,
	real3 const pt
) {
	return calc_<?=fluid?>_EKin(W, pt) + calc_<?=fluid?>_EInt(solver, W);
}

static inline real calc_<?=fluid?>_Cs(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const * const W
) {
	return sqrt(solver->heatCapacityRatio * (W)-><?=fluid?>_P / (W)-><?=fluid?>_rho);
}

static inline real calc_<?=fluid?>_P(
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const * const U,
	real3 const pt
) {
	return (solver->heatCapacityRatio - 1.) * (/*EInt=*/(U)-><?=fluid?>_ETotal  - /*EKin=*/calc_<?=fluid?>_EKin_fromCons(U, pt));
}

static inline real calc_<?=fluid?>_Cs_fromCons(
	constant <?=solver_t?> const * const solver,
	<?=cons_t?> const * const U,
	real3 const pt
) {
	return sqrt(solver->heatCapacityRatio 
		* calc_<?=fluid?>_P(solver, U, pt)
		/ (U)-><?=fluid?>_rho);
}

<? end ?>

static inline real calc_EM_energy(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const pt
) {
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return .5 * (coordLenSq(U->D, pt) / eps + coordLenSq(U->B, pt) / mu);
}

//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=eqn_common?>

#define <?=primFromCons?>(\
	/*<?=prim_t?> * const */W,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(fluids) do ?>\
	real const <?=fluid?>_EKin = calc_<?=fluid?>_EKin_fromCons(U, pt);\
	real const <?=fluid?>_EInt = (U)-><?=fluid?>_ETotal - <?=fluid?>_EKin;\
<? end ?>\
<? for _,fluid in ipairs(fluids) do ?>\
	(W)-><?=fluid?>_rho = (U)-><?=fluid?>_rho;\
	(W)-><?=fluid?>_v = real3_real_mul((U)-><?=fluid?>_m, 1./(U)-><?=fluid?>_rho);\
	(W)-><?=fluid?>_P = (solver->heatCapacityRatio - 1.) * <?=fluid?>_EInt;\
<? end ?>\
	(W)->D = (U)->D;\
	(W)->B = (U)->B;\
	(W)->psi = (U)->psi;\
	(W)->phi = (U)->phi;\
	(W)->ePot = (U)->ePot;\
}

//// MODULE_NAME: <?=consFromPrim?>
//// MODULE_DEPENDS: <?=solver_t?> <?=eqn_common?>

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */U,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(fluids) do ?>\
	(U)-><?=fluid?>_rho = (W)-><?=fluid?>_rho;\
	(U)-><?=fluid?>_m = real3_real_mul((W)-><?=fluid?>_v, (W)-><?=fluid?>_rho);\
	(U)-><?=fluid?>_ETotal = calc_<?=fluid?>_ETotal(solver, W, pt);\
<? end ?>\
	(U)->D = (W)->D;\
	(U)->B = (W)->B;\
	(U)->psi = (W)->psi;\
	(U)->phi = (W)->phi;\
	(U)->ePot = (W)->ePot;\
}

//// MODULE_NAME: <?=apply_dU_dW?>
//// MODULE_DEPENDS: <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>
// only used by PLM

#define <?=apply_dU_dW?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const WA_<?=fluid?>_vL = coord_lower((WA)-><?=fluid?>_v, pt);\
<? end ?>\
<? for _,fluid in ipairs(fluids) do ?>\
	(result)-><?=fluid?>_rho = (W)-><?=fluid?>_rho;\
	(result)-><?=fluid?>_m = real3_add(\
		real3_real_mul((WA)-><?=fluid?>_v, (W)-><?=fluid?>_rho), \
		real3_real_mul((W)-><?=fluid?>_v, (WA)-><?=fluid?>_rho));\
	(result)-><?=fluid?>_ETotal = (W)-><?=fluid?>_rho * .5 * real3_dot((WA)-><?=fluid?>_v, WA_<?=fluid?>_vL) \
		+ (WA)-><?=fluid?>_rho * real3_dot((W)-><?=fluid?>_v, WA_<?=fluid?>_vL)\
		+ (W)-><?=fluid?>_P / (solver->heatCapacityRatio - 1.);\
<? end ?>\
	(result)->B = (W)->B;\
	(result)->D = (W)->D;\
	(result)->phi = (W)->phi;\
	(result)->psi = (W)->psi;\
	(result)->ePot = (W)->ePot;\
}

//// MODULE_NAME: <?=apply_dW_dU?>
//// MODULE_DEPENDS: <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>

#define <?=apply_dW_dU?>(\
	/*<?=prim_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const WA_<?=fluid?>_vL = coord_lower((WA)-><?=fluid?>_v, pt);\
<? end ?>\
<? for _,fluid in ipairs(fluids) do ?>\
	(result)-><?=fluid?>_rho = (U)-><?=fluid?>_rho;\
	(result)-><?=fluid?>_v = real3_sub(\
		real3_real_mul((U)-><?=fluid?>_m, 1. / (WA)-><?=fluid?>_rho),\
		real3_real_mul((WA)-><?=fluid?>_v, (U)-><?=fluid?>_rho / (WA)-><?=fluid?>_rho));\
	(result)-><?=fluid?>_P = (solver->heatCapacityRatio - 1.) * (\
		.5 * real3_dot((WA)-><?=fluid?>_v, WA_<?=fluid?>_vL) * (U)-><?=fluid?>_rho \
		- real3_dot((U)-><?=fluid?>_m, WA_<?=fluid?>_vL)\
		+ (U)-><?=fluid?>_ETotal);\
<? end ?>\
	(result)->B = (U)->B;\
	(result)->D = (U)->D;\
	(result)->phi = (U)->phi;\
	(result)->psi = (U)->psi;\
	(result)->ePot = (U)->ePot;\
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=cartesianToCoord?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool const lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
<? 
if eqn.useEulerInitState then 
?>
	real rho = 0.;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
<?
else
	 for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_rho = 0;
	real3 <?=fluid?>_v = real3_zero;
	real <?=fluid?>_P = 0;
	real <?=fluid?>_ePot = 0;
<? 
	end 
end
?>	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=fromreal?>(1.);
	<?=scalar?> permittivity = <?=fromreal?>(1. / (4. * M_PI));
	<?=scalar?> permeability = <?=fromreal?>(4. * M_PI);

	<?=initCode()?>

	<?=prim_t?> W;
<? 
if eqn.useEulerInitState then 
?>
	W.ion_rho = rho;
	W.elec_rho = rho / solver->ionElectronMassRatio;

	/* "the electron pressure is taken to be elec_P = 5 ion_rho" */
	/* is that arbitrary? */
	W.elec_P = 5. * rho;
	
	/* "the ion pressure is 1/100th the electron pressure" */
	/* is that from the mass ratio of ion/electron? */
	W.ion_P = P / solver->ionElectronMassRatio;

	W.ion_v = cartesianToCoord(v, x);
	W.elec_v = cartesianToCoord(v, x);

<?	
else	-- expect the initCond to explicitly provide the ion_ and elec_ Euler fluid variables
	for _,fluid in ipairs(fluids) do ?>
	W.<?=fluid?>_rho = <?=fluid?>_rho;
	W.<?=fluid?>_v = cartesianToCoord(<?=fluid?>_v, x);
	W.<?=fluid?>_P = <?=fluid?>_P;
<?
	end
end
?>
	W.D = cartesianToCoord(D, x);
	W.B = cartesianToCoord(B, x);
	W.psi = 0;
	W.phi = 0;
	
	W.ePot = 0;
	<?=consFromPrim?>(U, solver, &W, x);
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: units <?=normal_t?> <?=primFromCons?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, (cell)->pos);\
\
<? --\
for _,fluid in ipairs(fluids) do --\
?>	real <?=fluid?>_vj = normal_vecDotN1(n, W.<?=fluid?>_v);\
	real <?=fluid?>_HTotal = (U)-><?=fluid?>_ETotal + W.<?=fluid?>_P;\
\
	(resultFlux)-><?=fluid?>_rho = normal_vecDotN1(n, (U)-><?=fluid?>_m);\
	(resultFlux)-><?=fluid?>_m = real3_real_mul((U)-><?=fluid?>_m, <?=fluid?>_vj);\
<? 	for i,xi in ipairs(xNames) do --\
?>	(resultFlux)-><?=fluid?>_m.<?=xi?> += normal_u1<?=xi?>(n) * W.<?=fluid?>_P;\
<? 	end --\
?>	(resultFlux)-><?=fluid?>_ETotal = <?=fluid?>_HTotal * <?=fluid?>_vj;\
<? --\
end --\
?>	(resultFlux)->ePot = 0.;\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	/* taken from glm-maxwell instead of the 2014 Abgrall, Kumar */\
	real3 const E = real3_real_mul((U)->D, 1. / eps);\
	real3 const H = real3_real_mul((U)->B, 1. / mu);\
	if (n.side == 0) {\
		(resultFlux)->D = real3((U)->phi * solver->divPhiWavespeed / unit_m_per_s, H.z, -H.y);\
		(resultFlux)->B = real3((U)->psi * solver->divPsiWavespeed / unit_m_per_s, -E.z, E.y);\
	} else if (n.side == 1) {\
		(resultFlux)->D = real3(-H.z, (U)->phi * solver->divPhiWavespeed / unit_m_per_s, H.x);\
		(resultFlux)->B = real3(E.z, (U)->psi * solver->divPsiWavespeed / unit_m_per_s, -E.x);\
	} else if (n.side == 2) {\
		(resultFlux)->D = real3(H.y, -H.x, (U)->phi * solver->divPhiWavespeed / unit_m_per_s);\
		(resultFlux)->B = real3(-E.y, E.x, (U)->psi * solver->divPsiWavespeed / unit_m_per_s);\
	}\
	(resultFlux)->phi = normal_vecDotN1(n, (U)->D) * solver->divPhiWavespeed / unit_m_per_s;\
	(resultFlux)->psi = normal_vecDotN1(n, (U)->B) * solver->divPsiWavespeed / unit_m_per_s;\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=primFromCons?>

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
	<?=prim_t?> WL;\
	<?=primFromCons?>(&WL, solver, UL, (cellL)->pos);\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, (cellR)->pos);\
\
<? for _,fluid in ipairs(fluids) do ?>\
\
	real const <?=fluid?>_sqrtRhoL = sqrt(WL.<?=fluid?>_rho);\
	real3 const <?=fluid?>_vL = WL.<?=fluid?>_v;\
	real const <?=fluid?>_hTotalL = calc_hTotal(solver, WL.<?=fluid?>_rho, WL.<?=fluid?>_P, (UL)-><?=fluid?>_ETotal);\
\
	real const <?=fluid?>_sqrtRhoR = sqrt(WR.<?=fluid?>_rho);\
	real3 const <?=fluid?>_vR = WR.<?=fluid?>_v;\
	real const <?=fluid?>_hTotalR = calc_hTotal(solver, WR.<?=fluid?>_rho, WR.<?=fluid?>_P, (UR)-><?=fluid?>_ETotal);\
\
	real const <?=fluid?>_invDenom = 1./(<?=fluid?>_sqrtRhoL + <?=fluid?>_sqrtRhoR);\
\
	/* Roe-averaged */\
	(resultEig)-><?=fluid?>_rho = <?=fluid?>_sqrtRhoL * <?=fluid?>_sqrtRhoR;\
	(resultEig)-><?=fluid?>_v = real3_add(\
			real3_real_mul(<?=fluid?>_vL, <?=fluid?>_sqrtRhoL * <?=fluid?>_invDenom),\
			real3_real_mul(<?=fluid?>_vR, <?=fluid?>_sqrtRhoR * <?=fluid?>_invDenom));\
	(resultEig)-><?=fluid?>_hTotal = <?=fluid?>_invDenom * (<?=fluid?>_sqrtRhoL * <?=fluid?>_hTotalL + <?=fluid?>_sqrtRhoR * <?=fluid?>_hTotalR);\
\
	/* derived: */\
	(resultEig)-><?=fluid?>_vSq = coordLenSq((resultEig)-><?=fluid?>_v, pt);\
	real const <?=fluid?>_eKin = .5 * (resultEig)-><?=fluid?>_vSq;\
	real const <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * ((resultEig)-><?=fluid?>_hTotal - <?=fluid?>_eKin);\
	(resultEig)-><?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);\
\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=primFromCons?>

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, (cell)->pos);\
<? for _,fluid in ipairs(fluids) do ?>\
	(resultEig)-><?=fluid?>_rho = W.<?=fluid?>_rho;\
	(resultEig)-><?=fluid?>_v = W.<?=fluid?>_v;\
	(resultEig)-><?=fluid?>_vSq = coordLenSq(W.<?=fluid?>_v, (cell)->pos);\
	real const <?=fluid?>_eKin = .5 * (resultEig)-><?=fluid?>_vSq;\
	(resultEig)-><?=fluid?>_hTotal = calc_hTotal(solver, W.<?=fluid?>_rho, W.<?=fluid?>_P, (U)-><?=fluid?>_ETotal);\
	real const <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * ((resultEig)-><?=fluid?>_hTotal - <?=fluid?>_eKin);\
	(resultEig)-><?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: units <?=eigen_t?> <?=waves_t?> <?=coord_lower?> <?=sqrt_2_and_1_2?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */UY,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */UX,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const nLen = normal_len(n);\
	real const nLenSq = nLen * nLen;\
	real const inv_nLen = 1. / nLen;\
	real const inv_nLenSq = 1. / nLenSq;\
\
	/* g^ij for fixed j=side */\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const <?=fluid?>_vL = coord_lower((eig)-><?=fluid?>_v, pt);\
	real const <?=fluid?>_denom = 2. * (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_Cs;\
	real const <?=fluid?>_invDenom = 1. / <?=fluid?>_denom;\
<? end ?>\
\
	real const heatRatioMinusOne = solver->heatCapacityRatio - 1.;\
\
	real3 const nU = normal_u1(n);\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	real const sqrt_eps = sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_mu = sqrt(mu);\
\
	if (n.side == 0) {\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=k+0?>] = (\
			  (UX)->ptr[<?=k+0?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x - (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=k+2?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=k+1?>] = (\
			  (UX)->ptr[<?=k+0?>] * (<?=fluid?>_denom - heatRatioMinusOne * (eig)-><?=fluid?>_vSq)\
			+ (UX)->ptr[<?=k+1?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=k+2?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.x * nU.y * inv_nLenSq - (eig)-><?=fluid?>_v.y)\
			+ (UX)->ptr[<?=k+1?>] * -nU.y * inv_nLenSq\
			+ (UX)->ptr[<?=k+2?>];\
		(UY)->ptr[<?=k+3?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.x * nU.z * inv_nLenSq - (eig)-><?=fluid?>_v.z)\
			+ (UX)->ptr[<?=k+1?>] * -nU.z * inv_nLenSq\
			+ (UX)->ptr[<?=k+3?>];\
		(UY)->ptr[<?=k+4?>] = (\
			  (UX)->ptr[<?=k+0?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x + (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=k+2?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
	k = k + 5 --\
end --\
?>\
		/* EM */\
		(UY)->ptr[<?=k+0?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+0?>] - (UX)->ptr[<?=k+6?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+1?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+3?>] - (UX)->ptr[<?=k+7?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+2?>] = ((((UX)->ptr[<?=k+2?>] * sqrt_mu) + ((UX)->ptr[<?=k+4?>] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+3?>] = ((((UX)->ptr[<?=k+1?>] * sqrt_mu) - ((UX)->ptr[<?=k+5?>] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[<?=k+4?>] = ((-(((UX)->ptr[<?=k+2?>] * sqrt_mu) - ((UX)->ptr[<?=k+4?>] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[<?=k+5?>] = ((((UX)->ptr[<?=k+1?>] * sqrt_mu) + ((UX)->ptr[<?=k+5?>] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[<?=k+6?>] = ((sqrt_eps * ((UX)->ptr[<?=k+0?>] + (UX)->ptr[<?=k+6?>])) / sqrt_2);\
		(UY)->ptr[<?=k+7?>] = ((sqrt_eps * ((UX)->ptr[<?=k+3?>] + (UX)->ptr[<?=k+7?>])) / sqrt_2);\
\
	} else if (n.side == 1) {\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=k+0?>] = (\
			  (UX)->ptr[<?=k+0?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y - (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=k+3?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=k+1?>] = \
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.y * nU.x * inv_nLenSq - (eig)-><?=fluid?>_v.x)\
			+ (UX)->ptr[<?=k+1?>]\
			+ (UX)->ptr[<?=k+2?>] * -nU.x * inv_nLenSq;\
		(UY)->ptr[<?=k+2?>] = (\
			  (UX)->ptr[<?=k+0?>] * (<?=fluid?>_denom - heatRatioMinusOne * (eig)-><?=fluid?>_vSq)\
			+ (UX)->ptr[<?=k+1?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=k+3?>] = \
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.y * nU.z * inv_nLenSq - (eig)-><?=fluid?>_v.z)\
			+ (UX)->ptr[<?=k+2?>] * -nU.z * inv_nLenSq\
			+ (UX)->ptr[<?=k+3?>];\
		(UY)->ptr[<?=k+4?>] = (\
			  (UX)->ptr[<?=k+0?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y + (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=k+3?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
	k = k + 5 --\
end --\
?>\
		/* EM */\
		(UY)->ptr[<?=k+0?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+1?>] - (UX)->ptr[<?=k+6?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+1?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+4?>] - (UX)->ptr[<?=k+7?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+2?>] = ((((UX)->ptr[<?=k+2?>] * sqrt_mu) - ((UX)->ptr[<?=k+3?>] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[<?=k+3?>] = ((((UX)->ptr[<?=k+0?>] * sqrt_mu) + ((UX)->ptr[<?=k+5?>] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+4?>] = ((((UX)->ptr[<?=k+2?>] * sqrt_mu) + ((UX)->ptr[<?=k+3?>] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[<?=k+5?>] = ((-(((UX)->ptr[<?=k+0?>] * sqrt_mu) - ((UX)->ptr[<?=k+5?>] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[<?=k+6?>] = ((sqrt_eps * ((UX)->ptr[<?=k+1?>] + (UX)->ptr[<?=k+6?>])) / sqrt_2);\
		(UY)->ptr[<?=k+7?>] = ((sqrt_eps * ((UX)->ptr[<?=k+4?>] + (UX)->ptr[<?=k+7?>])) / sqrt_2);\
\
	} else if (n.side == 2) {\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=k+0?>] = (\
			  (UX)->ptr[<?=k+0?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z - (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=k+4?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=k+1?>] = \
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.z * nU.x * inv_nLenSq - (eig)-><?=fluid?>_v.x)\
			+ (UX)->ptr[<?=k+1?>]\
			+ (UX)->ptr[<?=k+3?>] * -nU.x * inv_nLenSq;\
		(UY)->ptr[<?=k+2?>] = \
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.z * nU.y * inv_nLenSq - (eig)-><?=fluid?>_v.y)\
			+ (UX)->ptr[<?=k+2?>]\
			+ (UX)->ptr[<?=k+3?>] * -nU.y * inv_nLenSq;\
		(UY)->ptr[<?=k+3?>] = (\
			  (UX)->ptr[<?=k+0?>] * (<?=fluid?>_denom - heatRatioMinusOne * (eig)-><?=fluid?>_vSq)\
			+ (UX)->ptr[<?=k+1?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=k+4?>] = (\
			  (UX)->ptr[<?=k+0?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z + (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=k+4?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
	k = k + 5 --\
end --\
?>\
		/* EM */\
		(UY)->ptr[<?=k+0?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+2?>] - (UX)->ptr[<?=k+6?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+1?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+5?>] - (UX)->ptr[<?=k+7?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+2?>] = ((((UX)->ptr[<?=k+1?>] * sqrt_mu) + ((UX)->ptr[<?=k+3?>] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+3?>] = ((((UX)->ptr[<?=k+0?>] * sqrt_mu) - ((UX)->ptr[<?=k+4?>] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[<?=k+4?>] = ((((UX)->ptr[<?=k+1?>] * sqrt_mu) - ((UX)->ptr[<?=k+3?>] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[<?=k+5?>] = ((((UX)->ptr[<?=k+0?>] * sqrt_mu) + ((UX)->ptr[<?=k+4?>] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[<?=k+6?>] = ((sqrt_eps * ((UX)->ptr[<?=k+2?>] + (UX)->ptr[<?=k+6?>])) / sqrt_2);\
		(UY)->ptr[<?=k+7?>] = ((sqrt_eps * ((UX)->ptr[<?=k+5?>] + (UX)->ptr[<?=k+7?>])) / sqrt_2);\
\
	}\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: units <?=eigen_t?> <?=waves_t?> <?=coord_lower?> <?=sqrt_2_and_1_2?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */UY,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */UX,	/* numWaves = 16 */\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
\
	/* g^ij for fixed j=side */\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const <?=fluid?>_vL = coord_lower((eig)-><?=fluid?>_v, pt);\
<? end ?>\
\
	real3 const nU = normal_u1(n);\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
	real const sqrt_eps = sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_mu = sqrt(mu);\
\
	if (n.side == 0) {\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=k+0?>] =\
			  (UX)->ptr[<?=k+0?>]\
			+ (UX)->ptr[<?=k+1?>]\
			+ (UX)->ptr[<?=k+4?>];\
		(UY)->ptr[<?=k+1?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=k+1?>] * (eig)-><?=fluid?>_v.x\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=k+2?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * nU.y * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * (eig)-><?=fluid?>_v.y\
			+ (UX)->ptr[<?=k+2?>]\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * nU.y * inv_nLen);\
		(UY)->ptr[<?=k+3?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * nU.z * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * (eig)-><?=fluid?>_v.z\
			+ (UX)->ptr[<?=k+3?>]\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * nU.z * inv_nLen);\
		(UY)->ptr[<?=k+4?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * (eig)-><?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=k+2?>] * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen);\
<? --\
	k = k + 5 --\
end --\
?>\
		/* EM */\
		(UY)->ptr[<?=k+0?>] = ((-((UX)->ptr[<?=k+0?>] - (UX)->ptr[<?=k+6?>])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+1?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+3?>] - (UX)->ptr[<?=k+5?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+2?>] = ((sqrt_eps * ((UX)->ptr[<?=k+2?>] - (UX)->ptr[<?=k+4?>])) / sqrt_2);\
		(UY)->ptr[<?=k+3?>] = ((-((UX)->ptr[<?=k+1?>] - (UX)->ptr[<?=k+7?>])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+4?>] = ((sqrt_mu * ((UX)->ptr[<?=k+2?>] + (UX)->ptr[<?=k+4?>])) / sqrt_2);\
		(UY)->ptr[<?=k+5?>] = ((sqrt_mu * ((UX)->ptr[<?=k+3?>] + (UX)->ptr[<?=k+5?>])) / sqrt_2);\
		(UY)->ptr[<?=k+6?>] = (((UX)->ptr[<?=k+0?>] + (UX)->ptr[<?=k+6?>]) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+7?>] = (((UX)->ptr[<?=k+1?>] + (UX)->ptr[<?=k+7?>]) / (sqrt_2 * sqrt_eps));\
\
	} else if (n.side == 1) {\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=k+0?>] =\
			  (UX)->ptr[<?=k+0?>]\
			+ (UX)->ptr[<?=k+2?>]\
			+ (UX)->ptr[<?=k+4?>];\
		(UY)->ptr[<?=k+1?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * nU.x * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>]\
			+ (UX)->ptr[<?=k+2?>] * (eig)-><?=fluid?>_v.x\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * nU.x * inv_nLen);\
		(UY)->ptr[<?=k+2?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=k+2?>] * (eig)-><?=fluid?>_v.y\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=k+3?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * nU.z * inv_nLen)\
			+ (UX)->ptr[<?=k+2?>] * (eig)-><?=fluid?>_v.z\
			+ (UX)->ptr[<?=k+3?>]\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * nU.z * inv_nLen);\
		(UY)->ptr[<?=k+4?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * (eig)-><?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=k+3?>] * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen);\
<? --\
	k = k + 5 --\
end --\
?>\
		/* EM */\
		(UY)->ptr[<?=k+0?>] = ((sqrt_eps * ((UX)->ptr[<?=k+3?>] - (UX)->ptr[<?=k+5?>])) / sqrt_2);\
		(UY)->ptr[<?=k+1?>] = ((-((UX)->ptr[<?=k+0?>] - (UX)->ptr[<?=k+6?>])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+2?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+2?>] - (UX)->ptr[<?=k+4?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+3?>] = ((sqrt_mu * ((UX)->ptr[<?=k+2?>] + (UX)->ptr[<?=k+4?>])) / sqrt_2);\
		(UY)->ptr[<?=k+4?>] = ((-((UX)->ptr[<?=k+1?>] - (UX)->ptr[<?=k+7?>])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+5?>] = ((sqrt_mu * ((UX)->ptr[<?=k+3?>] + (UX)->ptr[<?=k+5?>])) / sqrt_2);\
		(UY)->ptr[<?=k+6?>] = (((UX)->ptr[<?=k+0?>] + (UX)->ptr[<?=k+6?>]) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+7?>] = (((UX)->ptr[<?=k+1?>] + (UX)->ptr[<?=k+7?>]) / (sqrt_2 * sqrt_eps));\
\
	} else if (n.side == 2) {\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=k+0?>] =\
			  (UX)->ptr[<?=k+0?>]\
			+ (UX)->ptr[<?=k+3?>]\
			+ (UX)->ptr[<?=k+4?>];\
		(UY)->ptr[<?=k+1?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * nU.x * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>]\
			+ (UX)->ptr[<?=k+3?>] * (eig)-><?=fluid?>_v.x\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * nU.x * inv_nLen);\
		(UY)->ptr[<?=k+2?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * nU.y * inv_nLen)\
			+ (UX)->ptr[<?=k+2?>]\
			+ (UX)->ptr[<?=k+3?>] * (eig)-><?=fluid?>_v.y\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * nU.y * inv_nLen);\
		(UY)->ptr[<?=k+3?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=k+3?>] * (eig)-><?=fluid?>_v.z\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=k+4?>] =\
			  (UX)->ptr[<?=k+0?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen)\
			+ (UX)->ptr[<?=k+1?>] * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=k+2?>] * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=k+3?>] * (eig)-><?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=k+4?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen);\
<? --\
	k = k + 5 --\
end --\
?>\
		/* EM */\
		(UY)->ptr[<?=k+0?>] = ((-(sqrt_eps * ((UX)->ptr[<?=k+3?>] - (UX)->ptr[<?=k+5?>]))) / sqrt_2);\
		(UY)->ptr[<?=k+1?>] = ((sqrt_eps * ((UX)->ptr[<?=k+2?>] - (UX)->ptr[<?=k+4?>])) / sqrt_2);\
		(UY)->ptr[<?=k+2?>] = ((-((UX)->ptr[<?=k+0?>] - (UX)->ptr[<?=k+6?>])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+3?>] = ((sqrt_mu * ((UX)->ptr[<?=k+2?>] + (UX)->ptr[<?=k+4?>])) / sqrt_2);\
		(UY)->ptr[<?=k+4?>] = ((sqrt_mu * ((UX)->ptr[<?=k+3?>] + (UX)->ptr[<?=k+5?>])) / sqrt_2);\
		(UY)->ptr[<?=k+5?>] = ((-((UX)->ptr[<?=k+1?>] - (UX)->ptr[<?=k+7?>])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+6?>] = (((UX)->ptr[<?=k+0?>] + (UX)->ptr[<?=k+6?>]) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[<?=k+7?>] = (((UX)->ptr[<?=k+1?>] + (UX)->ptr[<?=k+7?>]) / (sqrt_2 * sqrt_eps));\
\
	}\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: units

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */UY,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */UX,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	/* g^ij for fixed j=side */\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const <?=fluid?>_vL = coord_lower((eig)-><?=fluid?>_v, (cell)->pos);\
	real const <?=fluid?>_v_n = normal_vecDotN1(n, (eig)-><?=fluid?>_v);\
<? end ?>\
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
	real3 const nU = normal_u1(n);\
\
<? --\
local k = 0 --\
for i,fluid	in ipairs(fluids) do --\
?>\
	(UY)-><?=fluid?>_rho =\
		  (UX)->ptr[<?=k+1?>] * nx\
		+ (UX)->ptr[<?=k+2?>] * ny\
		+ (UX)->ptr[<?=k+3?>] * nz;\
	(UY)-><?=fluid?>_m.x =\
		  (UX)->ptr[<?=k+0?>] * (-<?=fluid?>_v_n * (eig)-><?=fluid?>_v.x + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.x)\
		+ (UX)->ptr[<?=k+1?>] * ((eig)-><?=fluid?>_v.x * nx - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)\
		+ (UX)->ptr[<?=k+2?>] * ((eig)-><?=fluid?>_v.x * ny - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.y)\
		+ (UX)->ptr[<?=k+3?>] * ((eig)-><?=fluid?>_v.x * nz - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.z)\
		+ (UX)->ptr[<?=k+4?>] * (solver->heatCapacityRatio - 1.) * nx;\
	(UY)-><?=fluid?>_m.y =\
		  (UX)->ptr[<?=k+0?>] * (-<?=fluid?>_v_n * (eig)-><?=fluid?>_v.y + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.y)\
		+ (UX)->ptr[<?=k+1?>] * ((eig)-><?=fluid?>_v.y * nx - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.x)\
		+ (UX)->ptr[<?=k+2?>] * ((eig)-><?=fluid?>_v.y * ny - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)\
		+ (UX)->ptr[<?=k+3?>] * ((eig)-><?=fluid?>_v.y * nz - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.z)\
		+ (UX)->ptr[<?=k+4?>] * (solver->heatCapacityRatio - 1.) * ny;\
	(UY)-><?=fluid?>_m.z =\
		  (UX)->ptr[<?=k+0?>] * (-<?=fluid?>_v_n * (eig)-><?=fluid?>_v.z + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.z)\
		+ (UX)->ptr[<?=k+1?>] * ((eig)-><?=fluid?>_v.z * nx - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.x)\
		+ (UX)->ptr[<?=k+2?>] * ((eig)-><?=fluid?>_v.z * ny - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.y)\
		+ (UX)->ptr[<?=k+3?>] * ((eig)-><?=fluid?>_v.z * nz - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)\
		+ (UX)->ptr[<?=k+4?>] * (solver->heatCapacityRatio - 1.) * nz;\
	(UY)-><?=fluid?>_ETotal =\
		  (UX)->ptr[<?=k+0?>] * <?=fluid?>_v_n * ((solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=k+1?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.x + nx * (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=k+2?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.y + ny * (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=k+3?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.z + nz * (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=k+4?>] * solver->heatCapacityRatio * <?=fluid?>_v_n;\
<? --\
	k = k + 5 --\
end --\
?>\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	real3 const E = real3_real_mul((UX)->D, 1. / eps);\
	real3 const H = real3_real_mul((UX)->B, 1. / mu);\
	if (n.side == 0) {\
		(UY)->D = real3(solver->divPhiWavespeed / unit_m_per_s * (UX)->phi, H.z, -H.y);\
		(UY)->B = real3(solver->divPsiWavespeed / unit_m_per_s * (UX)->psi, -E.z, E.y);\
	} else if (n.side == 1) {\
		(UY)->D = real3(-H.z, solver->divPhiWavespeed / unit_m_per_s * (UX)->phi, H.x);\
		(UY)->B = real3(E.z, solver->divPsiWavespeed / unit_m_per_s * (UX)->psi, -E.x);\
	} else if (n.side == 2) {\
		(UY)->D = real3(H.y, -H.x, solver->divPhiWavespeed / unit_m_per_s * (UX)->phi);\
		(UY)->B = real3(-E.y, E.x, solver->divPsiWavespeed / unit_m_per_s * (UX)->psi);\
	}\
	(UY)->phi = solver->divPhiWavespeed / unit_m_per_s * normal_vecDotN1(n, (UX)->D);\
	(UY)->psi = solver->divPsiWavespeed / unit_m_per_s * normal_vecDotN1(n, (UX)->B);\
	(UY)->ePot = 0;\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: units <?=eqn_common?> <?=elecChargeMassRatio?> <?=primFromCons?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;

	/*
	electric force:
	kg/(m^2 s^2) = C/kg * (kg/m^3 * C/m^2 / ((C^2 s^2)/(kg m^3))
	kg/(m^2 s^2) = C/kg * (kg^2/(C m^2 s^2)
	kg/(m^2 s^2) = kg/(m^2 s^2) 

	magnetic force
	kg/(m^2 s^2) = C/kg * (kg/(m^2 s) kg/(C s))
	kg/(m^2 s^2) = C/kg * (kg^2/(C m^2 s^2)
	kg/(m^2 s^2) = kg/(m^2 s^2) 
	*/
	deriv->ion_m.x += solver->ionChargeMassRatio / unit_C_per_kg * (U->ion_rho * U->D.x / eps + U->ion_m.y * U->B.z - U->ion_m.z * U->B.y);
	deriv->ion_m.y += solver->ionChargeMassRatio / unit_C_per_kg * (U->ion_rho * U->D.y / eps + U->ion_m.z * U->B.x - U->ion_m.x * U->B.z);
	deriv->ion_m.z += solver->ionChargeMassRatio / unit_C_per_kg * (U->ion_rho * U->D.z / eps + U->ion_m.x * U->B.y - U->ion_m.y * U->B.x);
	
	/*
	kg/(m s^3) = C/kg * C/m^2 * kg/(m^2 s) / ((C^2 s^2)/(kg m^3))
	kg/(m s^3) = C/kg * (C kg)/(m^4 s) * (kg m^3)/(C^2 s^2)
	kg/(m s^3) = C/kg * (kg^2)/(C m s^3)
	kg/(m s^3) = kg/(m s^3)
	*/
	deriv->ion_ETotal += solver->ionChargeMassRatio / unit_C_per_kg * real3_dot(U->D, U->ion_m) / eps;

	deriv->elec_m.x -= elecChargeMassRatio / unit_C_per_kg * (U->elec_rho * U->D.x / eps + U->elec_m.y * U->B.z - U->elec_m.z * U->B.y);
	deriv->elec_m.y -= elecChargeMassRatio / unit_C_per_kg * (U->elec_rho * U->D.y / eps + U->elec_m.z * U->B.x - U->elec_m.x * U->B.z);
	deriv->elec_m.z -= elecChargeMassRatio / unit_C_per_kg * (U->elec_rho * U->D.z / eps + U->elec_m.x * U->B.y - U->elec_m.y * U->B.x);
	
	deriv->elec_ETotal -= elecChargeMassRatio / unit_C_per_kg * real3_dot(U->D, U->elec_m) / eps;

	/*
	C/(m^2 s) = kg/(m^2*s) * C/kg
	C/(m^2 s) = C/(m^2*s)
	*/
	real3 J;
	J.x = (U->ion_m.x * solver->ionChargeMassRatio + U->elec_m.x * elecChargeMassRatio) / unit_C_per_kg;
	J.y = (U->ion_m.y * solver->ionChargeMassRatio + U->elec_m.y * elecChargeMassRatio) / unit_C_per_kg;
	J.z = (U->ion_m.z * solver->ionChargeMassRatio + U->elec_m.z * elecChargeMassRatio) / unit_C_per_kg;

	/* source of D is the current */
	deriv->D.x -= J.x;
	deriv->D.y -= J.y;
	deriv->D.z -= J.z;
	
	/* source of phi is the charge */
	deriv->phi += eps * (U->ion_rho * solver->ionChargeMassRatio + U->elec_rho * elecChargeMassRatio) / unit_C_per_kg * solver->divPhiWavespeed / unit_m_per_s;

<? if not require "hydro.coord.cartesian":isa(solver.coord) then ?>
	/* connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system */
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	real3 const conn1_u = coord_conn_trace23(x);
	<? for _,fluid in ipairs(fluids) do ?>{
		real3 const m_conn_vv = coord_conn_apply23(W.<?=fluid?>_v, U-><?=fluid?>_m, x);
		deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, m_conn_vv);	/* -Conn^i_jk rho v^j v^k  */
		deriv-><?=fluid?>_m = real3_add(deriv-><?=fluid?>_m, real3_real_mul(coord_raise(coord_conn_trace13(x), x), W.<?=fluid?>_P));		/* +Conn^j_kj g^ki P */
		deriv-><?=fluid?>_m = real3_add(deriv-><?=fluid?>_m, real3_real_mul(coord_conn_apply13(W.<?=fluid?>_v, U-><?=fluid?>_m, x), (solver->heatCapacityRatio - 1.) ));	/* + (gamma-1) rho v^k v^l Gamma_kjl g^ij	 */
		deriv-><?=fluid?>_ETotal -= (solver->heatCapacityRatio - 1.) * coord_conn_apply123(W.<?=fluid?>_v, W.<?=fluid?>_v, U-><?=fluid?>_m, x);		/* - (gamma-1) rho v^j v^k v^l Gamma_jkl */
	}<? end ?>
<? end ?>
}

//// MODULE_NAME: <?=constrainU?>
//// MODULE_DEPENDS: <?=primFromCons?> <?=consFromPrim?>

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?> * const U = UBuf + index;
	real3 const x = cellBuf[index].pos;
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);

<? for _,fluid in ipairs(fluids) do
?>	W.<?=fluid?>_rho = max((real)W.<?=fluid?>_rho, (real)solver->min_<?=fluid?>_rho);
	W.<?=fluid?>_P = max((real)W.<?=fluid?>_P, (real)solver->min_<?=fluid?>_P);
<? end
?>
	<?=consFromPrim?>(U, solver, &W, x);
}
