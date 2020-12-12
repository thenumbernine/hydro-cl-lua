//// MODULE_NAME: elecChargeMassRatio

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

//// MODULE_NAME: sqrt_2_and_1_2

#define sqrt_1_2 <?=("%.50f"):format(math.sqrt(.5))?>
#define sqrt_2 <?=("%.50f"):format(math.sqrt(2))?>

//// MODULE_NAME: primFromCons
//// MODULE_DEPENDS: <?=solver_t?> eqn.common <?=prim_t?> <?=cons_t?>

#define primFromCons(\
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
	(W)->D_g = (U)->D_g;\
	(W)->B_g = (U)->B_g;\
	(W)->psi_g = (U)->psi_g;\
	(W)->phi_g = (U)->phi_g;\
}

//// MODULE_NAME: consFromPrim
//// MODULE_DEPENDS: <?=solver_t?> eqn.common <?=prim_t?> <?=cons_t?> consFromPrim

#define consFromPrim(\
	/*<?=cons_t?> * const */U,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */x\
) {\
<? for _,fluid in ipairs(fluids) do ?>\
	(U)-><?=fluid?>_rho = (W)-><?=fluid?>_rho;\
	(U)-><?=fluid?>_m = real3_real_mul((W)-><?=fluid?>_v, (W)-><?=fluid?>_rho);\
	(U)-><?=fluid?>_ETotal = calc_<?=fluid?>_ETotal(solver, W, x);\
<? end ?>\
	(U)->D = (W)->D;\
	(U)->B = (W)->B;\
	(U)->psi = (W)->psi;\
	(U)->phi = (W)->phi;\
	(U)->D_g = (W)->D_g;\
	(U)->B_g = (W)->B_g;\
	(U)->psi_g = (W)->psi_g;\
	(U)->phi_g = (W)->phi_g;\
}

//// MODULE_NAME: apply_dU_dW
//// MODULE_DEPENDS: real3 coord_lower <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define apply_dU_dW(\
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
	(result)->B_g = (W)->B_g;\
	(result)->D_g = (W)->D_g;\
	(result)->phi_g = (W)->phi_g;\
	(result)->psi_g = (W)->psi_g;\
}

//// MODULE_NAME: apply_dW_dU
//// MODULE_DEPENDS: real3 coord_lower <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define apply_dW_dU(\
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
	(result)->B_g = (U)->B_g;\
	(result)->D_g = (U)->D_g;\
	(result)->phi_g = (U)->phi_g;\
	(result)->psi_g = (U)->psi_g;\
}

//// MODULE_NAME: eqn.common
//// MODULE_DEPENDS: units coordLenSq cartesianToCoord

#define /*real3*/ calc_EField(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
) \
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

static inline real calc_H(constant <?=solver_t?> const * const solver, real const P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
static inline real calc_h(constant <?=solver_t?> const * const solver, real const rho, real const P) { return calc_H(solver, P) / rho; }
static inline real calc_hTotal(constant <?=solver_t?> const * const solver, real const rho, real const P, real const ETotal) { return (P + ETotal) / rho; }
static inline real calc_HTotal(real const P, real const ETotal) { return P + ETotal; }

<? for _,fluid in ipairs(fluids) do ?>
static inline real calc_<?=fluid?>_eKin(<?=prim_t?> const * const W, real3 const x) { return .5 * coordLenSq((W)-><?=fluid?>_v, x); }
static inline real calc_<?=fluid?>_EKin(<?=prim_t?> const * const W, real3 const x) { return (W)-><?=fluid?>_rho * calc_<?=fluid?>_eKin(W, x); }
static inline real calc_<?=fluid?>_EInt(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) { return (W)-><?=fluid?>_P / (solver->heatCapacityRatio - 1.); }
static inline real calc_<?=fluid?>_eInt(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) { return calc_<?=fluid?>_EInt(solver, W) / (W)-><?=fluid?>_rho; }

static inline real calc_<?=fluid?>_EKin_fromCons(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	return .5 * coordLenSq((U)-><?=fluid?>_m, x) / (U)-><?=fluid?>_rho;
}

static inline real calc_<?=fluid?>_ETotal(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W, real3 const x) {
	return calc_<?=fluid?>_EKin(W, x) + calc_<?=fluid?>_EInt(solver, W);
}
static inline real calc_<?=fluid?>_Cs(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) {
	return sqrt(solver->heatCapacityRatio * W-><?=fluid?>_P / W-><?=fluid?>_rho);
}
<? end ?>

static inline real calc_EM_energy(constant <?=solver_t?> const * const solver, global <?=cons_t?> const * const U, real3 const x) {
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return .5 * (coordLenSq(U->D, x) / eps + coordLenSq(U->B, x) / mu);
}


/*
units:
eps_g = 1 / (4 pi G) 
mu_g = 4 pi G / c^2
[eps_g] = kg s^2 / m^3
[mu_g] = m / kg
[E_g] = [D_g / eps_g]
kg/m^2 * m^3/(kg s^2)
m/s^2
(rho * D_g / eps_g + m * B_g) / c
[rho * D_g / eps_g]
kg/m^3 * kg/m^2 * m^3 / (kg s^2)
kg/(m^2 s^2)
[m * B_g]
kg/(m^2 s) 1/s = kg/(m^2 s^2)
kg/m^3 * m/s^2 = kg / (m^2 s^2)
densitized force, in units of kg/(m^2 s^2)
*/
real3 calcIonGravForce(constant <?=solver_t?> const * const solver, global <?=cons_t?> const * const U, real3 const x) {
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real const eps_g = 1. / (4. * M_PI * G);
	return _real3(
		U->ion_rho * U->D_g.x / eps_g + 4. * (U->ion_m.y * U->B_g.z - U->ion_m.z * U->B_g.y),
		U->ion_rho * U->D_g.y / eps_g + 4. * (U->ion_m.z * U->B_g.x - U->ion_m.x * U->B_g.z),
		U->ion_rho * U->D_g.z / eps_g + 4. * (U->ion_m.x * U->B_g.y - U->ion_m.y * U->B_g.x));
}

real3 calcElecGravForce(constant <?=solver_t?> const * const solver, global <?=cons_t?> const * const U, real3 const x) {
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real const eps_g = 1. / (4. * M_PI * G);
	return _real3(
		U->elec_rho * U->D_g.x / eps_g + 4. * (U->elec_m.y * U->B_g.z - U->elec_m.z * U->B_g.y),
		U->elec_rho * U->D_g.y / eps_g + 4. * (U->elec_m.z * U->B_g.x - U->elec_m.x * U->B_g.z),
		U->elec_rho * U->D_g.z / eps_g + 4. * (U->elec_m.x * U->B_g.y - U->elec_m.y * U->B_g.x));
}

//// MODULE_NAME: applyInitCond

kernel void applyInitCond(
	constant <?=solver_t?> const * const solver,
	constant initCond_t const * const initCond,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	real3 const x = cellBuf[index].pos;
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

	// intel OpenCL compiler bug crashing when I initialize W with a struct assign
	<?=prim_t?> W;
<? 
if eqn.useEulerInitState then 
?>
	W.ion_rho = rho;
	W.elec_rho = rho / solver->ionElectronMassRatio;

	// "the electron pressure is taken to be elec_P = 5 ion_rho"
	// is that arbitrary?
	W.elec_P = 5. * rho;
	
	// "the ion pressure is 1/100th the electron pressure"
	// is that from the mass ratio of ion/electron?
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
	W.D_g = <?=vec3?>_zero;
	W.B_g = <?=vec3?>_zero;
	W.psi_g = 0;
	W.phi_g = 0;
	consFromPrim(UBuf + index, solver, &W, x);
}


//// MODULE_NAME: fluxFromCons
//// MODULE_DEPENDS: units normal_t

#define fluxFromCons(\
	/*<?=cons_t?> const * const */F,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=prim_t?> W;\
	primFromCons(&W, solver, U, pt);\
\
<? --\
for _,fluid in ipairs(fluids) do --\
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
?>\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
	/* taken from glm-maxwell instead of the 2014 Abgrall, Kumar */\
	/*  then replace D = epsilon E and phi' -> epsilon phi */\
	<? for _,suffix in ipairs{"", "_g"} do ?>{\
		real3 const E = real3_real_mul((U)->D<?=suffix?>, 1. / eps<?=suffix?>);\
		real3 const H = real3_real_mul((U)->B<?=suffix?>, 1. / mu<?=suffix?>);\
		if (n.side == 0) {\
			(F)->D<?=suffix?> = _real3((U)->phi<?=suffix?> * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s, H.z, -H.y);\
			(F)->B<?=suffix?> = _real3((U)->psi<?=suffix?> * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s, -E.z, E.y);\
		} else if (n.side == 1) {\
			(F)->D<?=suffix?> = _real3(-H.z, (U)->phi<?=suffix?> * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s, H.x);\
			(F)->B<?=suffix?> = _real3(E.z, (U)->psi<?=suffix?> * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s, -E.x);\
		} else if (n.side == 2) {\
			(F)->D<?=suffix?> = _real3(H.y, -H.x, (U)->phi<?=suffix?> * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s);\
			(F)->B<?=suffix?> = _real3(-E.y, E.x, (U)->psi<?=suffix?> * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s);\
		}\
		(F)->phi<?=suffix?> = normal_vecDotN1(n, (U)->D<?=suffix?>) * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s;\
		(F)->psi<?=suffix?> = normal_vecDotN1(n, (U)->B<?=suffix?>) * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s;\
	}<? end ?>\
}

//// MODULE_NAME: eigen_forInterface

#define eigen_forInterface(\
	/*<?=eigen_t?> * const */eig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=prim_t?> WL;\
	primFromCons(&WL, solver, UL, pt);\
	<?=prim_t?> WR;\
	primFromCons(&WR, solver, UR, pt);\
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
	(eig)-><?=fluid?>_rho = <?=fluid?>_sqrtRhoL * <?=fluid?>_sqrtRhoR;\
	(eig)-><?=fluid?>_v = real3_add(\
			real3_real_mul(<?=fluid?>_vL, <?=fluid?>_sqrtRhoL * <?=fluid?>_invDenom),\
			real3_real_mul(<?=fluid?>_vR, <?=fluid?>_sqrtRhoR * <?=fluid?>_invDenom));\
	(eig)-><?=fluid?>_hTotal = <?=fluid?>_invDenom * (<?=fluid?>_sqrtRhoL * <?=fluid?>_hTotalL + <?=fluid?>_sqrtRhoR * <?=fluid?>_hTotalR);\
\
	/* derived: */\
	(eig)-><?=fluid?>_vSq = coordLenSq((eig)-><?=fluid?>_v, pt);\
	real const <?=fluid?>_eKin = .5 * (eig)-><?=fluid?>_vSq;\
	real const <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * ((eig)-><?=fluid?>_hTotal - <?=fluid?>_eKin);\
	(eig)-><?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);\
\
<? end ?>\
}

//// MODULE_NAME: eigen_forCell

#define eigen_forCell(\
	/*<?=eigen_t?> * const */eig, \
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=prim_t?> W;\
	primFromCons(&W, solver, U, pt);\
<? for _,fluid in ipairs(fluids) do ?>\
	(eig)-><?=fluid?>_rho = W.<?=fluid?>_rho;\
	(eig)-><?=fluid?>_v = W.<?=fluid?>_v;\
	(eig)-><?=fluid?>_vSq = coordLenSq(W.<?=fluid?>_v, pt);\
	real const <?=fluid?>_eKin = .5 * (eig)-><?=fluid?>_vSq;\
	(eig)-><?=fluid?>_hTotal = calc_hTotal(solver, W.<?=fluid?>_rho, W.<?=fluid?>_P, (U)-><?=fluid?>_ETotal);\
	real const <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * ((eig)-><?=fluid?>_hTotal - <?=fluid?>_eKin);\
	(eig)-><?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);\
<? end ?>\
}

//// MODULE_NAME: eigen_left/rightTransform
//// MODULE_DEPENDS: sqrt_2_and_1_2

#define eigen_leftTransform(\
	/*<?=waves_t?> * const */UY,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */UX,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real* X = UX.ptr;\
\
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
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
	real const sqrt_eps = sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_mu = sqrt(mu);\
\
	real const sqrt_eps_g = sqrt(eps_g);	/*  TODO sqrt units */\
	real const sqrt_mu_g = sqrt(mu_g);\
\
	if (n.side == 0) {\
<? --\
					for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x - (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-4?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * (eig)-><?=fluid?>_vSq)\
			+ (UX)->ptr[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-3?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.x * nU.y * inv_nLenSq - (eig)-><?=fluid?>_v.y)\
			+ (UX)->ptr[<?=5*i-4?>] * -nU.y * inv_nLenSq\
			+ (UX)->ptr[<?=5*i-3?>];\
		(UY)->ptr[<?=5*i-2?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.x * nU.z * inv_nLenSq - (eig)-><?=fluid?>_v.z)\
			+ (UX)->ptr[<?=5*i-4?>] * -nU.z * inv_nLenSq\
			+ (UX)->ptr[<?=5*i-2?>];\
		(UY)->ptr[<?=5*i-1?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x + (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
					end --\
?>\
		/* EM & gravity */\
		<? for i,suffix in ipairs{"", "_g"} do ?>{\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+0?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+1?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+2?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] * sqrt_mu<?=suffix?>) + ((UX)->ptr[<?=5*#fluids+8*(i-1)+4?>] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+3?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] * sqrt_mu<?=suffix?>) - ((UX)->ptr[<?=5*#fluids+8*(i-1)+5?>] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+4?>] = ((-(((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] * sqrt_mu<?=suffix?>) - ((UX)->ptr[<?=5*#fluids+8*(i-1)+4?>] * sqrt_eps<?=suffix?>))) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+5?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] * sqrt_mu<?=suffix?>) + ((UX)->ptr[<?=5*#fluids+8*(i-1)+5?>] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+6?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+7?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>])) / sqrt_2);\
		}<? end ?>\
\
	} else if (n.side == 1) {\
<? --\
					for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y - (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-4?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.y * nU.x * inv_nLenSq - (eig)-><?=fluid?>_v.x)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-3?>] * -nU.x * inv_nLenSq;\
		(UY)->ptr[<?=5*i-3?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * (eig)-><?=fluid?>_vSq)\
			+ (UX)->ptr[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-2?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.y * nU.z * inv_nLenSq - (eig)-><?=fluid?>_v.z)\
			+ (UX)->ptr[<?=5*i-3?>] * -nU.z * inv_nLenSq\
			+ (UX)->ptr[<?=5*i-2?>];\
		(UY)->ptr[<?=5*i-1?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y + (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
					end --\
?>\
		/* EM & gravity */\
		<? for i,suffix in ipairs{"", "_g"} do ?>{\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+0?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+1?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+4?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+2?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] * sqrt_mu<?=suffix?>) - ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+3?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] * sqrt_mu<?=suffix?>) + ((UX)->ptr[<?=5*#fluids+8*(i-1)+5?>] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+4?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] * sqrt_mu<?=suffix?>) + ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+5?>] = ((-(((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] * sqrt_mu<?=suffix?>) - ((UX)->ptr[<?=5*#fluids+8*(i-1)+5?>] * sqrt_eps<?=suffix?>))) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+6?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+7?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+4?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>])) / sqrt_2);\
		}<? end ?>\
\
	} else if (n.side == 2) {\
<? --\
					for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z - (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-4?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.z * nU.x * inv_nLenSq - (eig)-><?=fluid?>_v.x)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-2?>] * -nU.x * inv_nLenSq;\
		(UY)->ptr[<?=5*i-3?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.z * nU.y * inv_nLenSq - (eig)-><?=fluid?>_v.y)\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-2?>] * -nU.y * inv_nLenSq;\
		(UY)->ptr[<?=5*i-2?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * (eig)-><?=fluid?>_vSq)\
			+ (UX)->ptr[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * -2. * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
		(UY)->ptr[<?=5*i-1?>] = (\
			  (UX)->ptr[<?=5*i-5?>] * (.5 * heatRatioMinusOne * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z + (eig)-><?=fluid?>_Cs * inv_nLen)\
			+ (UX)->ptr[<?=5*i-1?>] * heatRatioMinusOne\
		) * <?=fluid?>_invDenom;\
<? --\
					end --\
?>\
		/* EM & gravity */\
		<? for i,suffix in ipairs{"", "_g"} do ?>{\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+0?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+1?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+5?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+2?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] * sqrt_mu<?=suffix?>) + ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+3?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] * sqrt_mu<?=suffix?>) - ((UX)->ptr[<?=5*#fluids+8*(i-1)+4?>] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+4?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] * sqrt_mu<?=suffix?>) - ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+5?>] = ((((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] * sqrt_mu<?=suffix?>) + ((UX)->ptr[<?=5*#fluids+8*(i-1)+4?>] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+6?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+7?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+5?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>])) / sqrt_2);\
		}<? end ?>\
\
	}\
}

#define eigen_rightTransform(\
	/*<?=cons_t?> * const */UY,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */UX,	/* numWaves = 26 */\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
\
	/* g^ij for fixed j=side */\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const <?=fluid?>_vL = coord_lower((eig)-><?=fluid?>_v, pt);\
<? end ?>\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
	real const sqrt_eps = sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_mu = sqrt(mu);\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
	real const sqrt_eps_g = sqrt(eps_g);	/*  TODO sqrt units */\
	real const sqrt_mu_g = sqrt(mu_g);\
\
	real3 const nU = normal_u1(n);\
\
	if (n.side == 0) {\
<? --\
					for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] =\
			  (UX)->ptr[<?=5*i-5?>]\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-1?>];\
		(UY)->ptr[<?=5*i-4?>] =\
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (eig)-><?=fluid?>_v.x\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=5*i-3?>] =\
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * nU.y * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (eig)-><?=fluid?>_v.y\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * nU.y * inv_nLen);\
		(UY)->ptr[<?=5*i-2?>] =\
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * nU.z * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (eig)-><?=fluid?>_v.z\
			+ (UX)->ptr[<?=5*i-2?>]\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * nU.z * inv_nLen);\
		(UY)->ptr[<?=5*i-1?>] =\
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * (eig)-><?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.x * inv_nLen);\
<? --\
					end --\
?>\
		/* EM & gravity */\
		<? for i,suffix in ipairs{"", "_g"} do ?>{\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+0?>] = ((-((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>])) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+1?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+5?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+2?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+4?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+3?>] = ((-((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>])) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+4?>] = ((sqrt_mu<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+4?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+5?>] = ((sqrt_mu<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+5?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+6?>] = (((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>]) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+7?>] = (((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>]) / (sqrt_2 * sqrt_eps<?=suffix?>));\
		}<? end ?>\
\
	} else if (n.side == 1) {\
<? --\
					for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] =\
			  (UX)->ptr[<?=5*i-5?>]\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-1?>];\
		(UY)->ptr[<?=5*i-4?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * nU.x * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-3?>] * (eig)-><?=fluid?>_v.x\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * nU.x * inv_nLen);\
		(UY)->ptr[<?=5*i-3?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * (eig)-><?=fluid?>_v.y\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=5*i-2?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * nU.z * inv_nLen)\
			+ (UX)->ptr[<?=5*i-3?>] * (eig)-><?=fluid?>_v.z\
			+ (UX)->ptr[<?=5*i-2?>]\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * nU.z * inv_nLen);\
		(UY)->ptr[<?=5*i-1?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * (eig)-><?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=5*i-2?>] * <?=fluid?>_vL.z\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.y * inv_nLen);\
<? --\
					end --\
?>\
		/* EM & gravity */\
		<? for i,suffix in ipairs{"", "_g"} do ?>{\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+0?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+5?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+1?>] = ((-((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>])) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+2?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+4?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+3?>] = ((sqrt_mu<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+4?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+4?>] = ((-((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>])) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+5?>] = ((sqrt_mu<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+5?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+6?>] = (((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>]) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+7?>] = (((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>]) / (sqrt_2 * sqrt_eps<?=suffix?>));\
		}<? end ?>\
\
	} else if (n.side == 2) {\
<? --\
					for i,fluid in ipairs(fluids) do --\
?>\
		(UY)->ptr[<?=5*i-5?>] = \
			  (UX)->ptr[<?=5*i-5?>]\
			+ (UX)->ptr[<?=5*i-2?>]\
			+ (UX)->ptr[<?=5*i-1?>];\
		(UY)->ptr[<?=5*i-4?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * nU.x * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>]\
			+ (UX)->ptr[<?=5*i-2?>] * (eig)-><?=fluid?>_v.x\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * nU.x * inv_nLen);\
		(UY)->ptr[<?=5*i-3?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * nU.y * inv_nLen)\
			+ (UX)->ptr[<?=5*i-3?>]\
			+ (UX)->ptr[<?=5*i-2?>] * (eig)-><?=fluid?>_v.y\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * nU.y * inv_nLen);\
		(UY)->ptr[<?=5*i-2?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * nLen)\
			+ (UX)->ptr[<?=5*i-2?>] * (eig)-><?=fluid?>_v.z\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * nLen);\
		(UY)->ptr[<?=5*i-1?>] = \
			  (UX)->ptr[<?=5*i-5?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen)\
			+ (UX)->ptr[<?=5*i-4?>] * <?=fluid?>_vL.x\
			+ (UX)->ptr[<?=5*i-3?>] * <?=fluid?>_vL.y\
			+ (UX)->ptr[<?=5*i-2?>] * (eig)-><?=fluid?>_vSq / 2.\
			+ (UX)->ptr[<?=5*i-1?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_v.z * inv_nLen);\
<? --\
					end --\
?>\
		/* EM & gravity */\
		<? for i,suffix in ipairs{"", "_g"} do ?>{\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+0?>] = ((-(sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+5?>]))) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+1?>] = ((sqrt_eps<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+4?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+2?>] = ((-((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>])) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+3?>] = ((sqrt_mu<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+2?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+4?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+4?>] = ((sqrt_mu<?=suffix?> * ((UX)->ptr[<?=5*#fluids+8*(i-1)+3?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+5?>])) / sqrt_2);\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+5?>] = ((-((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] - (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>])) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+6?>] = (((UX)->ptr[<?=5*#fluids+8*(i-1)+0?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+6?>]) / (sqrt_2 * sqrt_eps<?=suffix?>));\
			(UY)->ptr[<?=5*#fluids+8*(i-1)+7?>] = (((UX)->ptr[<?=5*#fluids+8*(i-1)+1?>] + (UX)->ptr[<?=5*#fluids+8*(i-1)+7?>]) / (sqrt_2 * sqrt_eps<?=suffix?>));\
		}<? end ?>\
\
	}\
}

//// MODULE_NAME: eigen_fluxTransform

#define eigen_fluxTransform(\
	/*<?=cons_t?> const * const */UY,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */UX,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
\
	/* g^ij for fixed j=side */\
<? for _,fluid in ipairs(fluids) do ?>\
	real3 const <?=fluid?>_vL = coord_lower(<?=fluid?>_v, pt);\
	real const <?=fluid?>_v_n = normal_vecDotN1(n, <?=fluid?>_v);\
<? end ?>\
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
\
	real3 const nU = normal_u1(n);\
\
<? --\
					for i,fluid	in ipairs(fluids) do --\
?>\
	(UY)-><?=fluid?>_rho =\
		  (UX)->ptr[<?=5*i-4?>] * nx \
		+ (UX)->ptr[<?=5*i-3?>] * ny \
		+ (UX)->ptr[<?=5*i-2?>] * nz;\
	(UY)-><?=fluid?>_m.x =\
		  (UX)->ptr[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.x + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.x)\
		+ (UX)->ptr[<?=5*i-4?>] * (<?=fluid?>_v.x * nx - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)\
		+ (UX)->ptr[<?=5*i-3?>] * (<?=fluid?>_v.x * ny - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.y)\
		+ (UX)->ptr[<?=5*i-2?>] * (<?=fluid?>_v.x * nz - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.z)\
		+ (UX)->ptr[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nx;\
	(UY)-><?=fluid?>_m.y =\
		  (UX)->ptr[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.y + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.y)\
		+ (UX)->ptr[<?=5*i-4?>] * (<?=fluid?>_v.y * nx - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.x)\
		+ (UX)->ptr[<?=5*i-3?>] * (<?=fluid?>_v.y * ny - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)\
		+ (UX)->ptr[<?=5*i-2?>] * (<?=fluid?>_v.y * nz - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.z)\
		+ (UX)->ptr[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * ny;\
	(UY)-><?=fluid?>_m.z =\
		  (UX)->ptr[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.z + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.z)\
		+ (UX)->ptr[<?=5*i-4?>] * (<?=fluid?>_v.z * nx - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.x)\
		+ (UX)->ptr[<?=5*i-3?>] * (<?=fluid?>_v.z * ny - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.y)\
		+ (UX)->ptr[<?=5*i-2?>] * (<?=fluid?>_v.z * nz - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)\
		+ (UX)->ptr[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nz;\
	(UY)-><?=fluid?>_ETotal =\
		  (UX)->ptr[<?=5*i-5?>] * <?=fluid?>_v_n * ((solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=5*i-4?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.x + nx * (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=5*i-3?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.y + ny * (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=5*i-2?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.z + nz * (eig)-><?=fluid?>_hTotal)\
		+ (UX)->ptr[<?=5*i-1?>] * solver->heatCapacityRatio * <?=fluid?>_v_n;\
<? --\
					end --\
?>\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	real3 const E = real3_real_mul((UX)->D, 1. / eps);\
	real3 const H = real3_real_mul((UX)->B, 1. / mu);\
	if (n.side == 0) {\
		(UY)->D = _real3(solver->divPhiWavespeed / unit_m_per_s * (UX)->phi, H.z, -H.y);\
		(UY)->B = _real3(solver->divPsiWavespeed / unit_m_per_s * (UX)->psi, -E.z, E.y);\
	} else if (n.side == 1) {\
		(UY)->D = _real3(-H.z, solver->divPhiWavespeed / unit_m_per_s * (UX)->phi, H.x);\
		(UY)->B = _real3(E.z, solver->divPsiWavespeed / unit_m_per_s * (UX)->psi, -E.x);\
	} else if (n.side == 2) {\
		(UY)->D = _real3(H.y, -H.x, solver->divPhiWavespeed / unit_m_per_s * (UX)->phi);\
		(UY)->B = _real3(-E.y, E.x, solver->divPsiWavespeed / unit_m_per_s * (UX)->psi);\
	}\
	(UY)->phi = solver->divPhiWavespeed / unit_m_per_s * normal_vecDotN1(n, (UX)->D);\
	(UY)->psi = solver->divPsiWavespeed / unit_m_per_s * normal_vecDotN1(n, (UX)->B);\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
	real3 const E_g = real3_real_mul((UX)->D_g, 1. / eps_g);\
	real3 const H_g = real3_real_mul((UX)->B_g, 1. / mu_g);\
	if (n.side == 0) {\
		(UY)->D_g = _real3(solver->divPhiWavespeed_g / unit_m_per_s * (UX)->phi_g, H_g.z, -H_g.y);\
		(UY)->B_g = _real3(solver->divPsiWavespeed_g / unit_m_per_s * (UX)->psi_g, -E_g.z, E_g.y);\
	} else if (n.side == 1) {\
		(UY)->D_g = _real3(-H_g.z, solver->divPhiWavespeed_g / unit_m_per_s * (UX)->phi_g, H_g.x);\
		(UY)->B_g = _real3(E_g.z, solver->divPsiWavespeed_g / unit_m_per_s * (UX)->psi_g, -E_g.x);\
	} else if (n.side == 2) {\
		(UY)->D_g = _real3(H_g.y, -H_g.x, solver->divPhiWavespeed_g / unit_m_per_s * (UX)->phi_g);\
		(UY)->B_g = _real3(-E_g.y, E_g.x, solver->divPsiWavespeed_g / unit_m_per_s * (UX)->psi_g);\
	}\
	(UY)->phi_g = solver->divPhiWavespeed_g / unit_m_per_s * normal_vecDotN1(n, E_g);\
	(UY)->psi_g = solver->divPsiWavespeed_g / unit_m_per_s * normal_vecDotN1(n, H_g);\
}

//// MODULE_NAME: addSource
//// MODULE_DEPENDS: eqn.common elecChargeMassRatio

kernel void addSource(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

	/* TODO double check all units */

	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real const eps_g = 1. / (4. * M_PI * G);
	real const mu_g = 1. / (eps_g * speedOfLightSq);

	/* rho * (E + v * B) has units kg/(m^2 s^2) */
	real3 const ionGravForce = calcIonGravForce(solver, U, x);

	/*
	electric force:
	kg/(m^2 s^2) = C/kg * (kg/m^3 * C/m^2 / ((C^2 s^2)/(kg m^3))
	kg/(m^2 s^2) = C/kg * (kg^2/(C m^2 s^2)
	kg/(m^2 s^2) = kg/(m^2 s^2) 

	magnetic force
	kg/(m^2 s^2) = C/kg * (kg/(m^2 s) kg/(C s))
	kg/(m^2 s^2) = C/kg * (kg^2/(C m^2 s^2)
	kg/(m^2 s^2) = kg/(m^2 s^2) 
	
	ion_m has units kg/m^3 * m/s = kg/(m^2 s)
	ion_m,t has units kg/(m^2 s^2)
	*/
	deriv->ion_m.x += solver->ionChargeMassRatio / unit_C_per_kg * (U->ion_rho * U->D.x / eps + U->ion_m.y * U->B.z - U->ion_m.z * U->B.y) + ionGravForce.x;
	deriv->ion_m.y += solver->ionChargeMassRatio / unit_C_per_kg * (U->ion_rho * U->D.y / eps + U->ion_m.z * U->B.x - U->ion_m.x * U->B.z) + ionGravForce.y;
	deriv->ion_m.z += solver->ionChargeMassRatio / unit_C_per_kg * (U->ion_rho * U->D.z / eps + U->ion_m.x * U->B.y - U->ion_m.y * U->B.x) + ionGravForce.z;
	
	/*
	kg/(m s^3) = C/kg * C/m^2 * kg/(m^2 s) / ((C^2 s^2)/(kg m^3))
	kg/(m s^3) = C/kg * (C kg)/(m^4 s) * (kg m^3)/(C^2 s^2)
	kg/(m s^3) = C/kg * (kg^2)/(C m s^3)
	kg/(m s^3) = kg/(m s^3)
	*/
	deriv->ion_ETotal += solver->ionChargeMassRatio / unit_C_per_kg * real3_dot(U->D, U->ion_m) / eps + real3_dot(ionGravForce, U->ion_m) / U->ion_rho;
	
	real3 elecGravForce = calcElecGravForce(solver, U, x);

	deriv->elec_m.x -= elecChargeMassRatio / unit_C_per_kg * (U->elec_rho * U->D.x / eps + U->elec_m.y * U->B.z - U->elec_m.z * U->B.y) + elecGravForce.x;
	deriv->elec_m.y -= elecChargeMassRatio / unit_C_per_kg * (U->elec_rho * U->D.y / eps + U->elec_m.z * U->B.x - U->elec_m.x * U->B.z) + elecGravForce.y;
	deriv->elec_m.z -= elecChargeMassRatio / unit_C_per_kg * (U->elec_rho * U->D.z / eps + U->elec_m.x * U->B.y - U->elec_m.y * U->B.x) + elecGravForce.z;
	
	deriv->elec_ETotal -= elecChargeMassRatio / unit_C_per_kg * real3_dot(U->D, U->elec_m) / eps + real3_dot(elecGravForce, U->elec_m) / U->elec_rho;

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

	/* source of D_g is J_g is the momentum + Poynting vector */
	real3 J_g = real3_zero;
	/*  I'm symmetrizing the stress-energy */
		/* matter */
	J_g.x -= U->ion_m.x + U->elec_m.x;
	J_g.y -= U->ion_m.y + U->elec_m.y;
	J_g.z -= U->ion_m.z + U->elec_m.z;
		/* electromagnetism */
	real const v_p = 1. / (solver->sqrt_eps * solver->sqrt_mu);		/*  phase vel = sqrt(1 / (eps * mu)) */
	real const v_pSq = v_p * v_p;
	real const nSq = 1. / v_pSq;		/*  index of refraction^2 = 1 / phase vel^2 */
	real3 const S = real3_real_mul(real3_cross(U->D, U->B), 1. / nSq);
	J_g.x -= .5 * S.x * (1. + nSq) / speedOfLightSq;
	J_g.y -= .5 * S.y * (1. + nSq) / speedOfLightSq;
	J_g.z -= .5 * S.z * (1. + nSq) / speedOfLightSq;
	
	deriv->D_g.x -= J_g.x;
	deriv->D_g.y -= J_g.y;
	deriv->D_g.z -= J_g.z;

	/*  source of phi_g is T_00 is rho + .5 (D^2 + B^2) */
	real T_00_over_c2 = 0;
		/* matter */
	T_00_over_c2 += U->ion_rho + U->elec_rho;
		/* electromagnetism			 */
	T_00_over_c2 += calc_EM_energy(solver, U, x) / speedOfLightSq;
	
	deriv->phi_g += T_00_over_c2 * solver->divPhiWavespeed_g / unit_m_per_s;


<? if not require "hydro.coord.cartesian".is(solver.coord) then ?>
	/* connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system */
	<?=prim_t?> W;
	primFromCons(&W, solver, U, x);
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

//// MODULE_NAME: constrainU
//// MODULE_DEPENDS: consFromPrim

kernel void constrainU(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	global <?=cons_t?> * const U = UBuf + index;
	real3 const x = cellBuf[index].pos;
	<?=prim_t?> W;
	primFromCons(&W, solver, U, x);

<? for _,fluid in ipairs(fluids) do
?>	W.<?=fluid?>_rho = max((real)W.<?=fluid?>_rho, (real)solver->min_<?=fluid?>_rho);
	W.<?=fluid?>_P = max((real)W.<?=fluid?>_P, (real)solver->min_<?=fluid?>_P);
<? end
?>
	consFromPrim(U, solver, &W, x);
}
