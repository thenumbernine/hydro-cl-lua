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

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>
#define sqrt_2 <?=('%.50f'):format(math.sqrt(2))?>

//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: units <?=coordLenSq?>

static inline real3 calc_EField(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U
) {
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	return real3_real_mul((U)->D, 1. / eps);
}
 
static inline real3 calc_HField(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U
) { 
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return real3_real_mul((U)->B, 1. / mu);
}

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
for _,fluid in ipairs(eqn.fluids) do 
?> + (U)-><?=fluid?>_rho<? 
end 
?>;
}

static inline real calc_rho_from_W(
	<?=prim_t?> const * const W
) {
	return 0.<?
for _,fluid in ipairs(eqn.fluids) do 
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

<? for _,fluid in ipairs(eqn.fluids) do ?>

static inline real calc_<?=fluid?>_eKin(
	<?=prim_t?> const * const W,
	real3 const x
) {
	return .5 * coordLenSq((W)-><?=fluid?>_v, x);
}

static inline real calc_<?=fluid?>_EKin(
	<?=prim_t?> const * const W,
	real3 const x
) {
	return (W)-><?=fluid?>_rho * calc_<?=fluid?>_eKin(W, x);
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

static inline real calc_<?=fluid?>_EKin_fromCons(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	return .5 * coordLenSq((U)-><?=fluid?>_m, x) / (U)-><?=fluid?>_rho;
}

static inline real calc_<?=fluid?>_ETotal(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const * const W,
	real3 const x
) {
	return calc_<?=fluid?>_EKin(W, x) + calc_<?=fluid?>_EInt(solver, W);
}

static inline real calc_<?=fluid?>_Cs(
	constant <?=solver_t?> const * const solver,
	<?=prim_t?> const * const W
) {
	return sqrt(solver->heatCapacityRatio * (W)-><?=fluid?>_P / (W)-><?=fluid?>_rho);
}

<? end ?>

static inline real calc_EM_energy(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const x
) {
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return .5 * (coordLenSq(U->D, x) / eps + coordLenSq(U->B, x) / mu);
}

//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=eqn_common?>

#define <?=primFromCons?>(\
	/*<?=prim_t?> * const */W,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
	real const <?=fluid?>_EKin = calc_<?=fluid?>_EKin_fromCons(U, pt);\
	real const <?=fluid?>_EInt = (U)-><?=fluid?>_ETotal - <?=fluid?>_EKin;\
<? end ?>\
	<? for _,fluid in ipairs(eqn.fluids) do ?>\
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
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=eqn_common?>

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */U,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
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
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>
// only used by PLM

#define <?=apply_dU_dW?>(\
	/*<?=cons_t?> * const */result\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA, \
	/*<?=prim_t?> const * const */W, \
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
	real3 WA_<?=fluid?>_vL = coord_lower((WA)-><?=fluid?>_v, pt);\
<? end ?>\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
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
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>

#define <?=apply_dW_dU?>(\
	/*<?=prim_t?> * const */result\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
	real3 WA_<?=fluid?>_vL = coord_lower((WA)-><?=fluid?>_v, pt);\
<? end ?>\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
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
	real ePot = 0;
<? 
if eqn.useEulerInitState then 
?>
	real rho = 0.;
	real3 v = real3_zero;
	real P = 0;
<?
else
	 for _,fluid in ipairs(eqn.fluids) do
?>	real <?=fluid?>_rho = 0;
	real3 <?=fluid?>_v = real3_zero;
	real <?=fluid?>_P = 0;
<? 
	end 
end
?>	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=fromreal?>(1.);
	<?=scalar?> permittivity = <?=fromreal?>(1. / (4. * M_PI));
	<?=scalar?> permeability = <?=fromreal?>(4. * M_PI);

	<?=initCode()?>

	<?=prim_t?> W = {
<? 
if eqn.useEulerInitState then 
?>
		.ion_rho = rho,
		.elec_rho = rho / solver->ionElectronMassRatio, 

		/*  "the electron pressure is taken to be elec_P = 5 ion_rho" */
		/*  is that arbitrary? */
		.elec_P = 5. * rho,
		
		/*  "the ion pressure is 1/100th the electron pressure" */
		/*  is that from the mass ratio of ion/electron? */
		.ion_P = P / solver->ionElectronMassRatio, 

		.ion_v = cartesianToCoord(v, x),
		.elec_v = cartesianToCoord(v, x),
	
<?	
else	-- expect the initCond to explicitly provide the ion_ and elec_ Euler fluid variables
	for _,fluid in ipairs(eqn.fluids) do ?>
		.<?=fluid?>_rho = <?=fluid?>_rho,
		.<?=fluid?>_v = cartesianToCoord(<?=fluid?>_v, x),
		.<?=fluid?>_P = <?=fluid?>_P,
<?
	end
end
?>
		.D = cartesianToCoord(D, x), 
		.B = cartesianToCoord(B, x),
		.psi = 0,
		.phi = 0,
	
		.ePot = 0,
	};
	<?=consFromPrim?>(U, solver, &W, x);
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: units <?=normal_t?> <?=primFromCons?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */F,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
\
<? --\
for _,fluid in ipairs(eqn.fluids) do --\
?>	real <?=fluid?>_vj = normal_vecDotN1(n, W.<?=fluid?>_v);\
	real <?=fluid?>_HTotal = (U)-><?=fluid?>_ETotal + W.<?=fluid?>_P;\
\
	(F)-><?=fluid?>_rho = normal_vecDotN1(n, (U)-><?=fluid?>_m);\
	(F)-><?=fluid?>_m = real3_real_mul((U)-><?=fluid?>_m, <?=fluid?>_vj);\
<? 	for i,xi in ipairs(xNames) do --\
?>	(F)-><?=fluid?>_m.<?=xi?> += normal_u1<?=xi?>(n) * W.<?=fluid?>_P;\
<? 	end --\
?>	(F)-><?=fluid?>_ETotal = <?=fluid?>_HTotal * <?=fluid?>_vj;\
<? --\
end --\
?>	(F)->ePot = 0.;\
	\
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	/* taken from glm-maxwell instead of the 2014 Abgrall, Kumar */\
	real3 E = real3_real_mul((U)->D, 1. / eps);\
	real3 H = real3_real_mul((U)->B, 1. / mu);\
	if (n.side == 0) {\
		(F)->D = _real3((U)->phi * solver->divPhiWavespeed / unit_m_per_s, H.z, -H.y);\
		(F)->B = _real3((U)->psi * solver->divPsiWavespeed / unit_m_per_s, -E.z, E.y);\
	} else if (n.side == 1) {\
		(F)->D = _real3(-H.z, (U)->phi * solver->divPhiWavespeed / unit_m_per_s, H.x);\
		(F)->B = _real3(E.z, (U)->psi * solver->divPsiWavespeed / unit_m_per_s, -E.x);\
	} else if (n.side == 2) {\
		(F)->D = _real3(H.y, -H.x, (U)->phi * solver->divPhiWavespeed / unit_m_per_s);\
		(F)->B = _real3(-E.y, E.x, (U)->psi * solver->divPsiWavespeed / unit_m_per_s);\
	}\
	(F)->phi = normal_vecDotN1(n, (U)->D) * solver->divPhiWavespeed / unit_m_per_s;\
	(F)->psi = normal_vecDotN1(n, (U)->B) * solver->divPsiWavespeed / unit_m_per_s;\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=primFromCons?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> WL;\
	<?=primFromCons?>(&WL, solver, UL, pt);\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, pt);\
\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
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
	(eig)-><?=fluid?>_rho = <?=fluid?>_sqrtRhoL * <?=fluid?>_sqrtRhoR;\
	(eig)-><?=fluid?>_v = real3_add(\
			real3_real_mul(<?=fluid?>_vL, <?=fluid?>_sqrtRhoL * <?=fluid?>_invDenom),\
			real3_real_mul(<?=fluid?>_vR, <?=fluid?>_sqrtRhoR * <?=fluid?>_invDenom));\
	(eig)-><?=fluid?>_hTotal = <?=fluid?>_invDenom * (<?=fluid?>_sqrtRhoL * <?=fluid?>_hTotalL + <?=fluid?>_sqrtRhoR * <?=fluid?>_hTotalR);\
	/* derived: */\
	(eig)-><?=fluid?>_vSq = coordLenSq((eig)-><?=fluid?>_v, pt);\
	real const <?=fluid?>_eKin = .5 * (eig)-><?=fluid?>_vSq;\
	real const <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * ((eig)-><?=fluid?>_hTotal - <?=fluid?>_eKin);\
	(eig)-><?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);\
\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=primFromCons?>

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
	real const <?=fluid?>_vSq = coordLenSq(W.<?=fluid?>_v, pt);\
	real const <?=fluid?>_eKin = .5 * <?=fluid?>_vSq;\
	real const <?=fluid?>_hTotal = calc_hTotal(solver, W.<?=fluid?>_rho, W.<?=fluid?>_P, (U)-><?=fluid?>_ETotal);\
	real const <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * (<?=fluid?>_hTotal - <?=fluid?>_eKin);\
	real const <?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);\
<? end ?>\
<? for _,fluid in ipairs(eqn.fluids) do ?>\
	(eig)-><?=fluid?>_rho = W.<?=fluid?>_rho;\
	(eig)-><?=fluid?>_v = W.<?=fluid?>_v;\
	(eig)-><?=fluid?>_hTotal = <?=fluid?>_hTotal;\
	(eig)-><?=fluid?>_vSq = <?=fluid?>_vSq;\
	(eig)-><?=fluid?>_Cs = <?=fluid?>_Cs;\
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
) { \
	real const nLen = normal_len(n);\
	real const nLenSq = nLen * nLen;\
\
	/* g^ij for fixed j=side */\
	real3 const ion_v = (eig)->ion_v;\
	real3 const ion_vL = coord_lower(ion_v, pt);\
	real const ion_hTotal = (eig)->ion_hTotal;\
	real const ion_vSq = real3_dot(ion_v, ion_vL);\
	real const ion_Cs = (eig)->ion_Cs;\
	real const ion_Cs_over_nLen = ion_Cs / nLen; \
	\
	real3 const elec_v = (eig)->elec_v;\
	real3 const elec_vL = coord_lower(elec_v, pt);\
	real const elec_hTotal = (eig)->elec_hTotal;\
	real const elec_vSq = real3_dot(elec_v, elec_vL);\
	real const elec_Cs = (eig)->elec_Cs;\
	real const elec_Cs_over_nLen = elec_Cs / nLen; \
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
	real const n1x = normal_l2x(n);\
	real const n1y = normal_l2y(n);\
	real const n1z = normal_l2z(n);\
	real const n2x = normal_l3x(n);\
	real const n2y = normal_l3y(n);\
	real const n2z = normal_l3z(n);\
	real3 const nU = normal_u1(n);\
\
	real3 const ion_v_ns = normal_vecDotNs(n, ion_v);\
	real const ion_v_n = ion_v_ns.x, ion_v_n1 = ion_v_ns.y, ion_v_n2 = ion_v_ns.z;\
	real3 const elec_v_ns = normal_vecDotNs(n, elec_v);\
	real const elec_v_n = elec_v_ns.x, elec_v_n1 = elec_v_ns.y, elec_v_n2 = elec_v_ns.z;\
\
\
	real const ion_denom = 2. * ion_Cs * ion_Cs;\
	real const ion_invDenom = 1. / ion_denom;\
	\
	real const elec_denom = 2. * elec_Cs * elec_Cs;\
	real const elec_invDenom = 1. / elec_denom;\
\
	real const heatRatioMinusOne = solver->heatCapacityRatio - 1.;\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
	real const sqrt_eps = sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_mu = sqrt(mu);\
\
	if (n.side == 0) {\
<? --\
					for i,fluid in ipairs(eqn.fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.x / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x - <?=fluid?>_Cs / nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-4?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)\
			+ (UX)->ptr[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-3?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.x * nU.y / nLenSq - <?=fluid?>_v.y)\
			+ (UX)->ptr[<?=5*i-4?>] * -nU.y / nLenSq\
			+ (UX)->ptr[<?=5*i-3?>];\
		(UY)->ptr[<?=5*i-2?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.x * nU.z / nLenSq - <?=fluid?>_v.z)\
			+ (UX)->ptr[<?=5*i-4?>] * -nU.z / nLenSq\
			+ (UX)->ptr[<?=5*i-2?>];\
		(UY)->ptr[<?=5*i-1?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.x / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x + <?=fluid?>_Cs / nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
					end --\
?>\
		/* EM */\
		(UY)->ptr[10] = ((-(sqrt_eps * ((UX)->ptr[10] - (UX)->ptr[16]))) / sqrt_2);\
		(UY)->ptr[11] = ((-(sqrt_eps * ((UX)->ptr[13] - (UX)->ptr[17]))) / sqrt_2);\
		(UY)->ptr[12] = ((((UX)->ptr[12] * sqrt_mu) + ((UX)->ptr[14] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(UY)->ptr[13] = ((((UX)->ptr[11] * sqrt_mu) - ((UX)->ptr[15] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[14] = ((-(((UX)->ptr[12] * sqrt_mu) - ((UX)->ptr[14] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[15] = ((((UX)->ptr[11] * sqrt_mu) + ((UX)->ptr[15] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[16] = ((sqrt_eps * ((UX)->ptr[10] + (UX)->ptr[16])) / sqrt_2);\
		(UY)->ptr[17] = ((sqrt_eps * ((UX)->ptr[13] + (UX)->ptr[17])) / sqrt_2);\
\
	} else if (n.side == 1) {\
	\
<? --\
					for i,fluid in ipairs(eqn.fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.y / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y - <?=fluid?>_Cs / nLen)\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-4?>] = \
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.y * nU.x / nLenSq - <?=fluid?>_v.x)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-3?>] * -nU.x / nLenSq;\
		(UY)->ptr[<?=5*i-3?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)\
			+ (UX)->ptr[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-2?>] = \
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.y * nU.z / nLenSq - <?=fluid?>_v.z)\
			+ (UX)->ptr[<?=5*i-3?>] * -nU.z / nLenSq\
			+ (UX)->ptr[<?=5*i-2?>];\
		(UY)->ptr[<?=5*i-1?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.y / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y + <?=fluid?>_Cs / nLen)\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
					end --\
?>\
		/* EM */\
		(UY)->ptr[10] = ((-(sqrt_eps * ((UX)->ptr[11] - (UX)->ptr[16]))) / sqrt_2);\
		(UY)->ptr[11] = ((-(sqrt_eps * ((UX)->ptr[14] - (UX)->ptr[17]))) / sqrt_2);\
		(UY)->ptr[12] = ((((UX)->ptr[12] * sqrt_mu) - ((UX)->ptr[13] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[13] = ((((UX)->ptr[10] * sqrt_mu) + ((UX)->ptr[15] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(UY)->ptr[14] = ((((UX)->ptr[12] * sqrt_mu) + ((UX)->ptr[13] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[15] = ((-(((UX)->ptr[10] * sqrt_mu) - ((UX)->ptr[15] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[16] = ((sqrt_eps * ((UX)->ptr[11] + (UX)->ptr[16])) / sqrt_2);\
		(UY)->ptr[17] = ((sqrt_eps * ((UX)->ptr[14] + (UX)->ptr[17])) / sqrt_2);\
\
	} else if (n.side == 2) {\
<? --\
					for i,fluid in ipairs(eqn.fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.z / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z - <?=fluid?>_Cs / nLen)\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-4?>] = \
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.z * nU.x / nLenSq - <?=fluid?>_v.x)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-2?>] * -nU.x / nLenSq;\
		(UY)->ptr[<?=5*i-3?>] = \
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.z * nU.y / nLenSq - <?=fluid?>_v.y)\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-2?>] * -nU.y / nLenSq;\
		(UY)->ptr[<?=5*i-2?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)\
			+ (UX)->ptr[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-1?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.z / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z + <?=fluid?>_Cs / nLen)\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
					end --\
?>\
		/* EM */\
		(UY)->ptr[10] = ((-(sqrt_eps * ((UX)->ptr[12] - (UX)->ptr[16]))) / sqrt_2);\
		(UY)->ptr[11] = ((-(sqrt_eps * ((UX)->ptr[15] - (UX)->ptr[17]))) / sqrt_2);\
		(UY)->ptr[12] = ((((UX)->ptr[11] * sqrt_mu) + ((UX)->ptr[13] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));\
		(UY)->ptr[13] = ((((UX)->ptr[10] * sqrt_mu) - ((UX)->ptr[14] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[14] = ((((UX)->ptr[11] * sqrt_mu) - ((UX)->ptr[13] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));\
		(UY)->ptr[15] = ((((UX)->ptr[10] * sqrt_mu) + ((UX)->ptr[14] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));\
		(UY)->ptr[16] = ((sqrt_eps * ((UX)->ptr[12] + (UX)->ptr[16])) / sqrt_2);\
		(UY)->ptr[17] = ((sqrt_eps * ((UX)->ptr[15] + (UX)->ptr[17])) / sqrt_2);\
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
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
	real const sqrt_eps = sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_mu = sqrt(mu);\
\
	real const nLen = normal_len(n);\
	real const nLenSq = nLen * nLen;\
\
	/* g^ij for fixed j=side */\
	real3 const ion_v = (eig)->ion_v;\
	real3 const ion_vL = coord_lower(ion_v, pt);\
	real const ion_hTotal = (eig)->ion_hTotal;\
	real const ion_vSq = real3_dot(ion_v, ion_vL);\
	real const ion_Cs = (eig)->ion_Cs;\
	real const ion_Cs_over_nLen = ion_Cs / nLen; \
	\
	real3 const elec_v = (eig)->elec_v;\
	real3 const elec_vL = coord_lower(elec_v, pt);\
	real const elec_hTotal = (eig)->elec_hTotal;\
	real const elec_vSq = real3_dot(elec_v, elec_vL);\
	real const elec_Cs = (eig)->elec_Cs;\
	real const elec_Cs_over_nLen = elec_Cs / nLen; \
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
	real const n1x = normal_l2x(n);\
	real const n1y = normal_l2y(n);\
	real const n1z = normal_l2z(n);\
	real const n2x = normal_l3x(n);\
	real const n2y = normal_l3y(n);\
	real const n2z = normal_l3z(n);\
	real3 const nU = normal_u1(n);\
\
	real3 const ion_v_ns = normal_vecDotNs(n, ion_v);\
	real const ion_v_n = ion_v_ns.x, ion_v_n1 = ion_v_ns.y, ion_v_n2 = ion_v_ns.z;\
	real3 const elec_v_ns = normal_vecDotNs(n, elec_v);\
	real const elec_v_n = elec_v_ns.x, elec_v_n1 = elec_v_ns.y, elec_v_n2 = elec_v_ns.z;\
\
\
	if (n.side == 0) {\
<? --\
					for i,fluid in ipairs(eqn.fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] =\
			  (UX)->ptr[<?=5*i-5?>]\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-1?>];\
		(UY)->ptr[<?=5*i-4?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_v.x\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=5*i-3?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * nU.y / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_v.y\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * nU.y / nLen);\
		(UY)->ptr[<?=5*i-2?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * nU.z / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_v.z\
			+ (UX)->ptr[<?=5*i-2?>]\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * nU.z / nLen);\
		(UY)->ptr[<?=5*i-1?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.x / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.x / nLen);\
<? --\
					end --\
?>\
		/* EM */\
		(UY)->ptr[10] = ((-((UX)->ptr[10] - (UX)->ptr[16])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[11] = ((-(sqrt_eps * ((UX)->ptr[13] - (UX)->ptr[15]))) / sqrt_2);\
		(UY)->ptr[12] = ((sqrt_eps * ((UX)->ptr[12] - (UX)->ptr[14])) / sqrt_2);\
		(UY)->ptr[13] = ((-((UX)->ptr[11] - (UX)->ptr[17])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[14] = ((sqrt_mu * ((UX)->ptr[12] + (UX)->ptr[14])) / sqrt_2);\
		(UY)->ptr[15] = ((sqrt_mu * ((UX)->ptr[13] + (UX)->ptr[15])) / sqrt_2);\
		(UY)->ptr[16] = (((UX)->ptr[10] + (UX)->ptr[16]) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[17] = (((UX)->ptr[11] + (UX)->ptr[17]) / (sqrt_2 * sqrt_eps));\
	\
	} else if (n.side == 1) {\
<? --\
					for i,fluid in ipairs(eqn.fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] =\
			  (UX)->ptr[<?=5*i-5?>]\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-1?>];\
		(UY)->ptr[<?=5*i-4?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * nU.x / nLen)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_v.x\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * nU.x / nLen);\
		(UY)->ptr[<?=5*i-3?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_v.y\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=5*i-2?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * nU.z / nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_v.z\
			+ (UX)->ptr[<?=5*i-2?>]\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * nU.z / nLen);\
		(UY)->ptr[<?=5*i-1?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.y / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.y / nLen);\
<? --\
					end --\
?>\
		/* EM */\
		(UY)->ptr[10] = ((sqrt_eps * ((UX)->ptr[13] - (UX)->ptr[15])) / sqrt_2);\
		(UY)->ptr[11] = ((-((UX)->ptr[10] - (UX)->ptr[16])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[12] = ((-(sqrt_eps * ((UX)->ptr[12] - (UX)->ptr[14]))) / sqrt_2);\
		(UY)->ptr[13] = ((sqrt_mu * ((UX)->ptr[12] + (UX)->ptr[14])) / sqrt_2);\
		(UY)->ptr[14] = ((-((UX)->ptr[11] - (UX)->ptr[17])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[15] = ((sqrt_mu * ((UX)->ptr[13] + (UX)->ptr[15])) / sqrt_2);\
		(UY)->ptr[16] = (((UX)->ptr[10] + (UX)->ptr[16]) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[17] = (((UX)->ptr[11] + (UX)->ptr[17]) / (sqrt_2 * sqrt_eps));\
\
	} else if (n.side == 2) {\
<? --\
					for i,fluid in ipairs(eqn.fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] =\
			  (UX)->ptr[<?=5*i-5?>]\
			+ (UX)->ptr[<?=5*i-2?>]\
			+ (UX)->ptr[<?=5*i-1?>];\
		(UY)->ptr[<?=5*i-4?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * nU.x / nLen)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_v.x\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * nU.x / nLen);\
		(UY)->ptr[<?=5*i-3?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * nU.y / nLen)\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_v.y\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * nU.y / nLen);\
		(UY)->ptr[<?=5*i-2?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_v.z\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=5*i-1?>] =\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.z / nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.z / nLen);\
<? --\
					end --\
?>\
		/* EM */\
		(UY)->ptr[10] = ((-(sqrt_eps * ((UX)->ptr[13] - (UX)->ptr[15]))) / sqrt_2);\
		(UY)->ptr[11] = ((sqrt_eps * ((UX)->ptr[12] - (UX)->ptr[14])) / sqrt_2);\
		(UY)->ptr[12] = ((-((UX)->ptr[10] - (UX)->ptr[16])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[13] = ((sqrt_mu * ((UX)->ptr[12] + (UX)->ptr[14])) / sqrt_2);\
		(UY)->ptr[14] = ((sqrt_mu * ((UX)->ptr[13] + (UX)->ptr[15])) / sqrt_2);\
		(UY)->ptr[15] = ((-((UX)->ptr[11] - (UX)->ptr[17])) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[16] = (((UX)->ptr[10] + (UX)->ptr[16]) / (sqrt_2 * sqrt_eps));\
		(UY)->ptr[17] = (((UX)->ptr[11] + (UX)->ptr[17]) / (sqrt_2 * sqrt_eps));\
\
	}\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: units

static inline void <?=eigen_fluxTransform?>(
	<?=cons_t?> * const UY,
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> const * const eig,
	<?=cons_t?> const * const UX,
	real3 const x,
	<?=normal_t?> const n
) {
	real const nLen = normal_len(n);
	real const nLenSq = nLen * nLen;

	/* g^ij for fixed j=side */
	real3 const ion_v = (eig)->ion_v;
	real3 const ion_vL = coord_lower(ion_v, x);
	real const ion_hTotal = (eig)->ion_hTotal;
	real const ion_vSq = real3_dot(ion_v, ion_vL);
	real const ion_Cs = (eig)->ion_Cs;
	real const ion_Cs_over_nLen = ion_Cs / nLen; 
	
	real3 const elec_v = (eig)->elec_v;
	real3 const elec_vL = coord_lower(elec_v, x);
	real const elec_hTotal = (eig)->elec_hTotal;
	real const elec_vSq = real3_dot(elec_v, elec_vL);
	real const elec_Cs = (eig)->elec_Cs;
	real const elec_Cs_over_nLen = elec_Cs / nLen; 

	real const nx = normal_l1x(n);
	real const ny = normal_l1y(n);
	real const nz = normal_l1z(n);
	real const n1x = normal_l2x(n);
	real const n1y = normal_l2y(n);
	real const n1z = normal_l2z(n);
	real const n2x = normal_l3x(n);
	real const n2y = normal_l3y(n);
	real const n2z = normal_l3z(n);
	real3 const nU = normal_u1(n);

	real3 const ion_v_ns = normal_vecDotNs(n, ion_v);
	real const ion_v_n = ion_v_ns.x, ion_v_n1 = ion_v_ns.y, ion_v_n2 = ion_v_ns.z;
	real3 const elec_v_ns = normal_vecDotNs(n, elec_v);
	real const elec_v_n = elec_v_ns.x, elec_v_n1 = elec_v_ns.y, elec_v_n2 = elec_v_ns.z;


	global real const * const X = (UX)->ptr;
<?
					for i,fluid	in ipairs(eqn.fluids) do 
?>
	(UY)-><?=fluid?>_rho = X[<?=5*i-4?>] * nx 
		+ X[<?=5*i-3?>] * ny 
		+ X[<?=5*i-2?>] * nz;
	(UY)-><?=fluid?>_m.x = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.x + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * nU.x)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.x * nx - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.x * ny - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.x * nz - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nx;
	(UY)-><?=fluid?>_m.y = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.y + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * nU.y)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.y * nx - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.y * ny - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.y * nz - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * ny;
	(UY)-><?=fluid?>_m.z = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.z + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * nU.z)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.z * nx - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.z * ny - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.z * nz - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nz;
	(UY)-><?=fluid?>_ETotal = X[<?=5*i-5?>] * <?=fluid?>_v_n * ((solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq - <?=fluid?>_hTotal)
		+ X[<?=5*i-4?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.x + nx * <?=fluid?>_hTotal)
		+ X[<?=5*i-3?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.y + ny * <?=fluid?>_hTotal)
		+ X[<?=5*i-2?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.z + nz * <?=fluid?>_hTotal)
		+ X[<?=5*i-1?>] * solver->heatCapacityRatio * <?=fluid?>_v_n;
<? 
					end 
?>	
	
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;

	real3 E = real3_real_mul(UX.D, 1. / eps);
	real3 H = real3_real_mul(UX.B, 1. / mu);
	if (n.side == 0) {
		(UY)->D = _real3(solver->divPhiWavespeed / unit_m_per_s * UX.phi, H.z, -H.y);
		(UY)->B = _real3(solver->divPsiWavespeed / unit_m_per_s * UX.psi, -E.z, E.y);
	} else if (n.side == 1) {
		(UY)->D = _real3(-H.z, solver->divPhiWavespeed / unit_m_per_s * UX.phi, H.x);
		(UY)->B = _real3(E.z, solver->divPsiWavespeed / unit_m_per_s * UX.psi, -E.x);
	} else if (n.side == 2) {
		(UY)->D = _real3(H.y, -H.x, solver->divPhiWavespeed / unit_m_per_s * UX.phi);
		(UY)->B = _real3(-E.y, E.x, solver->divPsiWavespeed / unit_m_per_s * UX.psi);
	}
	(UY)->phi = solver->divPhiWavespeed / unit_m_per_s * normal_vecDotN1(n, UX.D);
	(UY)->psi = solver->divPsiWavespeed / unit_m_per_s * normal_vecDotN1(n, UX.B);
	(UY)->ePot = 0;
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: units <?=elecChargeMassRatio?> <?=primFromCons?> 

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
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
	deriv->D.x -= (U->ion_m.x * solver->ionChargeMassRatio + U->elec_m.x * elecChargeMassRatio) / unit_C_per_kg;
	deriv->D.y -= (U->ion_m.y * solver->ionChargeMassRatio + U->elec_m.y * elecChargeMassRatio) / unit_C_per_kg;
	deriv->D.z -= (U->ion_m.z * solver->ionChargeMassRatio + U->elec_m.z * elecChargeMassRatio) / unit_C_per_kg;
	
	deriv->phi += eps * (U->ion_rho * solver->ionChargeMassRatio + U->elec_rho * elecChargeMassRatio) / unit_C_per_kg * solver->divPhiWavespeed / unit_m_per_s;

<? if not require 'hydro.coord.cartesian'.is(solver.coord) then ?>
	real3 x = cellBuf[index].pos;
	/* connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system */
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	real3 const conn1_u = coord_conn_trace23(x);
	<? for _,fluid in ipairs(eqn.fluids) do ?>{
		real3 m_conn_vv = coord_conn_apply23(W.<?=fluid?>_v, U-><?=fluid?>_m, x);
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

<? for _,fluid in ipairs(eqn.fluids) do
?>	W.<?=fluid?>_rho = max((real)W.<?=fluid?>_rho, (real)solver->min_<?=fluid?>_rho);
	W.<?=fluid?>_P = max((real)W.<?=fluid?>_P, (real)solver->min_<?=fluid?>_P);
<? end
?>
	<?=consFromPrim?>(U, solver, &W, x);
}
