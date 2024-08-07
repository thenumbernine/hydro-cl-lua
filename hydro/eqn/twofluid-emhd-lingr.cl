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

//// MODULE_NAME: <?=cons_setEB?>

#define <?=cons_setEB?>(/*cons_t & */U, /*real3 */Eval, /*real3 */Bval) {\
	(U).B = Bval;\
	(U).D = real3_real_mul(Eval, 1. / (solver->sqrt_eps * solver->sqrt_eps));\
}

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

//units are kg/(m*s^2)
static inline real calc_EM_energy(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const pt
) {
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return .5 * (coordLenSq(U->D, pt) / eps + coordLenSq(U->B, pt) / mu);
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
	(W)->D_g = (U)->D_g;\
	(W)->B_g = (U)->B_g;\
	(W)->psi_g = (U)->psi_g;\
	(W)->phi_g = (U)->phi_g;\
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
	(U)->D_g = (W)->D_g;\
	(U)->B_g = (W)->B_g;\
	(U)->psi_g = (W)->psi_g;\
	(U)->phi_g = (W)->phi_g;\
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
	(result)->B_g = (W)->B_g;\
	(result)->D_g = (W)->D_g;\
	(result)->phi_g = (W)->phi_g;\
	(result)->psi_g = (W)->psi_g;\
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
	(result)->B_g = (U)->B_g;\
	(result)->D_g = (U)->D_g;\
	(result)->phi_g = (U)->phi_g;\
	(result)->psi_g = (U)->psi_g;\
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
	W.D_g = <?=vec3?>_zero;
	W.B_g = <?=vec3?>_zero;
	W.psi_g = 0;
	W.phi_g = 0;
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
?>	real const <?=fluid?>_v_n = normal_vecDotN1(n, W.<?=fluid?>_v);\
	(resultFlux)-><?=fluid?>_rho = normal_vecDotN1(n, (U)-><?=fluid?>_m);\
	(resultFlux)-><?=fluid?>_m = real3_add(\
		real3_real_mul((U)-><?=fluid?>_m, <?=fluid?>_v_n),\
		real3_real_mul(normal_u1(n), W.<?=fluid?>_P)\
	);\
	real const <?=fluid?>_HTotal = (U)-><?=fluid?>_ETotal + W.<?=fluid?>_P;\
	(resultFlux)-><?=fluid?>_ETotal = <?=fluid?>_HTotal * <?=fluid?>_v_n;\
<? --\
end --\
?>\
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
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
\
		(resultFlux)->D.x = H.y * nz - H.z * ny + nx * (U)->phi * solver->divPhiWavespeed / unit_m_per_s;	/* F_D^i = -eps^ijk n_j H_k */\
		(resultFlux)->B.x = E.z * ny - E.y * nz + nx * (U)->psi * solver->divPsiWavespeed / unit_m_per_s;	/* F_B^i = +eps^ijk n_j B_k */\
\
		(resultFlux)->D.y = H.z * nx - H.x * nz + ny * (U)->phi * solver->divPhiWavespeed / unit_m_per_s;\
		(resultFlux)->B.y = E.x * nz - E.z * nx + ny * (U)->psi * solver->divPsiWavespeed / unit_m_per_s;\
\
		(resultFlux)->D.z = H.x * ny - H.y * nx + nz * (U)->phi * solver->divPhiWavespeed / unit_m_per_s;\
		(resultFlux)->B.z = E.y * nx - E.x * ny + nz * (U)->psi * solver->divPsiWavespeed / unit_m_per_s;\
\
		(resultFlux)->phi<?=suffix?> = normal_vecDotN1(n, (U)->D<?=suffix?>) * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s;\
		(resultFlux)->psi<?=suffix?> = normal_vecDotN1(n, (U)->B<?=suffix?>) * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s;\
	}<? end ?>\
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
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const nLen = normal_len(n);\
	real const nLenSq = nLen * nLen;\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
	real const inv_nLen = 1. / nLen;\
	real const inv_nLenSq = 1. / nLenSq;\
	real const gamma_1 = solver->heatCapacityRatio - 1.;\
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
	real const sqrt_1_eps = 1./sqrt(eps);		/* TODO sqrt units */\
	real const sqrt_1_mu = 1./sqrt(mu);\
\
	real const sqrt_1_eps_g = 1./sqrt(eps_g);	/* TODO sqrt units */\
	real const sqrt_1_mu_g = 1./sqrt(mu_g);\
\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
	{\
		real3 const v_n = normal_vecDotNs(n, (eig)-><?=fluid?>_v);\
		real3 const vL = coord_lower((eig)-><?=fluid?>_v, pt);\
		real const denom = 2. * (eig)-><?=fluid?>_Cs * (eig)-><?=fluid?>_Cs;\
		real const invDenom = 1. / denom;\
		(result)->ptr[0+<?=k?>] = (\
				(X)->ptr[0+<?=k?>] * (.5 * gamma_1 * (eig)-><?=fluid?>_vSq + (eig)-><?=fluid?>_Cs * v_n.x * inv_nLen)\
				+ (X)->ptr[1+<?=k?>] * (-gamma_1 * vL.x - (eig)-><?=fluid?>_Cs * normal_l1x_over_len(n))\
				+ (X)->ptr[2+<?=k?>] * (-gamma_1 * vL.y - (eig)-><?=fluid?>_Cs * normal_l1y_over_len(n))\
				+ (X)->ptr[3+<?=k?>] * (-gamma_1 * vL.z - (eig)-><?=fluid?>_Cs * normal_l1z_over_len(n))\
				+ (X)->ptr[4+<?=k?>] * gamma_1\
			) * invDenom;\
		(result)->ptr[1+<?=k?>] =\
			(\
				(X)->ptr[0+<?=k?>] * (denom - gamma_1 * (eig)-><?=fluid?>_vSq)\
				+ (X)->ptr[1+<?=k?>] * 2. * gamma_1 * vL.x\
				+ (X)->ptr[2+<?=k?>] * 2. * gamma_1 * vL.y\
				+ (X)->ptr[3+<?=k?>] * 2. * gamma_1 * vL.z\
				+ (X)->ptr[4+<?=k?>] * -2. * gamma_1\
			) * invDenom;\
		(result)->ptr[2+<?=k?>] =\
			(X)->ptr[0+<?=k?>] * -v_n.y\
			+ (X)->ptr[1+<?=k?>] * normal_l2x(n)\
			+ (X)->ptr[2+<?=k?>] * normal_l2y(n)\
			+ (X)->ptr[3+<?=k?>] * normal_l2z(n);\
		(result)->ptr[3+<?=k?>] =\
			(X)->ptr[0+<?=k?>] * -v_n.z\
			+ (X)->ptr[1+<?=k?>] * normal_l3x(n)\
			+ (X)->ptr[2+<?=k?>] * normal_l3y(n)\
			+ (X)->ptr[3+<?=k?>] * normal_l3z(n);\
		(result)->ptr[4+<?=k?>] =\
			(\
				(X)->ptr[0+<?=k?>] * (.5 * gamma_1 * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_Cs * v_n.x * inv_nLen)\
				+ (X)->ptr[1+<?=k?>] * (-gamma_1 * vL.x + (eig)-><?=fluid?>_Cs * normal_l1x_over_len(n))\
				+ (X)->ptr[2+<?=k?>] * (-gamma_1 * vL.y + (eig)-><?=fluid?>_Cs * normal_l1y_over_len(n))\
				+ (X)->ptr[3+<?=k?>] * (-gamma_1 * vL.z + (eig)-><?=fluid?>_Cs * normal_l1z_over_len(n))\
				+ (X)->ptr[4+<?=k?>] * gamma_1\
			) * invDenom;\
	}\
<? k=k+5 end ?>\
		/* EM & gravity */\
<? --\
for i,suffix in ipairs{"", "_g"} do --\
?>\
	{\
		real const tmp1 = 1. / sqrt_1_eps<?=suffix?>;\
		real const tmp2 = tmp1 * n_l.x;\
		real const tmp3 = 1. / 2.;\
		real const tmp4 = (X)->ptr[0+<?=k?>] * tmp2;\
		real const tmp5 = tmp3 * tmp4;\
		real const tmp7 = n_l.y * tmp1;\
		real const tmp9 = (X)->ptr[1+<?=k?>] * tmp7;\
		real const tmp11 = tmp9 * tmp3;\
		real const tmp13 = n_l.z * tmp1;\
		real const tmp15 = (X)->ptr[2+<?=k?>] * tmp13;\
		real const tmp17 = tmp15 * tmp3;\
		real const tmp20 = (X)->ptr[6+<?=k?>] * tmp1;\
		real const tmp22 = tmp20 * tmp3;\
		real const tmp26 = (X)->ptr[3+<?=k?>] * tmp2;\
		real const tmp27 = tmp26 * tmp3;\
		real const tmp31 = (X)->ptr[4+<?=k?>] * tmp7;\
		real const tmp33 = tmp31 * tmp3;\
		real const tmp37 = (X)->ptr[5+<?=k?>] * tmp13;\
		real const tmp39 = tmp37 * tmp3;\
		real const tmp42 = (X)->ptr[7+<?=k?>] * tmp1;\
		real const tmp44 = tmp42 * tmp3;\
		real const tmp45 = sqrt_1_mu<?=suffix?> * n2_l.z;\
		real const tmp47 = (X)->ptr[5+<?=k?>] * tmp45;\
		real const tmp48 = sqrt_1_mu<?=suffix?> * n2_l.y;\
		real const tmp50 = (X)->ptr[4+<?=k?>] * tmp48;\
		real const tmp51 = sqrt_1_mu<?=suffix?> * n2_l.x;\
		real const tmp53 = (X)->ptr[3+<?=k?>] * tmp51;\
		real const tmp54 = sqrt_1_eps<?=suffix?> * n3_l.z;\
		real const tmp56 = (X)->ptr[2+<?=k?>] * tmp54;\
		real const tmp57 = sqrt_1_eps<?=suffix?> * n3_l.y;\
		real const tmp59 = (X)->ptr[1+<?=k?>] * tmp57;\
		real const tmp60 = sqrt_1_eps<?=suffix?> * n3_l.x;\
		real const tmp62 = (X)->ptr[0+<?=k?>] * tmp60;\
		real const tmp63 = tmp59 * tmp3;\
		real const tmp64 = tmp62 * tmp3;\
		real const tmp65 = tmp56 * tmp3;\
		real const tmp67 = tmp53 * tmp3;\
		real const tmp69 = tmp50 * tmp3;\
		real const tmp71 = tmp47 * tmp3;\
		real const tmp73 = sqrt_1_eps<?=suffix?> * n2_l.x;\
		real const tmp75 = (X)->ptr[0+<?=k?>] * tmp73;\
		real const tmp76 = tmp75 * tmp3;\
		real const tmp77 = sqrt_1_eps<?=suffix?> * n2_l.y;\
		real const tmp79 = (X)->ptr[1+<?=k?>] * tmp77;\
		real const tmp81 = tmp79 * tmp3;\
		real const tmp82 = sqrt_1_eps<?=suffix?> * n2_l.z;\
		real const tmp84 = (X)->ptr[2+<?=k?>] * tmp82;\
		real const tmp86 = tmp84 * tmp3;\
		real const tmp87 = sqrt_1_mu<?=suffix?> * n3_l.x;\
		real const tmp89 = (X)->ptr[3+<?=k?>] * tmp87;\
		real const tmp90 = sqrt_1_mu<?=suffix?> * n3_l.y;\
		real const tmp92 = (X)->ptr[4+<?=k?>] * tmp90;\
		real const tmp93 = sqrt_1_mu<?=suffix?> * n3_l.z;\
		real const tmp95 = (X)->ptr[5+<?=k?>] * tmp93;\
		real const tmp96 = tmp92 * tmp3;\
		real const tmp97 = tmp95 * tmp3;\
		real const tmp98 = tmp89 * tmp3;\
		(result)->ptr[0+<?=k?>] = tmp22 + -tmp5 - tmp11 - tmp17;\
		(result)->ptr[1+<?=k?>] = tmp44 + -tmp27 - tmp33 - tmp39;\
		(result)->ptr[2+<?=k?>] = tmp71 + tmp69 + tmp67 + tmp65 + tmp63 + tmp64;\
		(result)->ptr[3+<?=k?>] = tmp98 + tmp96 + tmp97 + -tmp76 - tmp81 - tmp86;\
		(result)->ptr[4+<?=k?>] = tmp67 + tmp69 + tmp71 + -tmp64 - tmp63 - tmp65;\
		(result)->ptr[5+<?=k?>] = tmp97 + tmp96 + tmp98 + tmp86 + tmp81 + tmp76;\
		(result)->ptr[6+<?=k?>] = tmp22 + tmp17 + tmp11 + tmp5;\
		(result)->ptr[7+<?=k?>] = tmp44 + tmp39 + tmp33 + tmp27;\
	}\
<? k=k+8 end ?>\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: units <?=eigen_t?> <?=waves_t?> <?=coord_lower?> <?=sqrt_2_and_1_2?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,	/* numWaves = 26 */\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
\
	real3 const nU = normal_u1(n);\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
	real const sqrt_1_eps = 1./sqrt(eps);	/*  TODO sqrt units */\
	real const sqrt_1_mu = 1./sqrt(mu);\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
	real const sqrt_1_eps_g = 1./sqrt(eps_g);	/*  TODO sqrt units */\
	real const sqrt_1_mu_g = 1./sqrt(mu_g);\
\
<? --\
local k = 0 --\
for i,fluid in ipairs(fluids) do --\
?>\
	{\
		real3 const v_n = normal_vecDotNs(n, (eig)-><?=fluid?>_v);\
		(result)->ptr[0+<?=k?>] =\
			(X)->ptr[0+<?=k?>]\
			+ (X)->ptr[1+<?=k?>]\
			+ (X)->ptr[4+<?=k?>];\
		(result)->ptr[1+<?=k?>] =\
			(X)->ptr[0+<?=k?>] * ((eig)-><?=fluid?>_v.x - (eig)-><?=fluid?>_Cs * normal_u1x_over_len(n))\
			+ (X)->ptr[1+<?=k?>] * (eig)-><?=fluid?>_v.x\
			+ (X)->ptr[2+<?=k?>] * normal_u2x(n)\
			+ (X)->ptr[3+<?=k?>] * normal_u3x(n)\
			+ (X)->ptr[4+<?=k?>] * ((eig)-><?=fluid?>_v.x + (eig)-><?=fluid?>_Cs * normal_u1x_over_len(n));\
		(result)->ptr[2+<?=k?>] =\
			(X)->ptr[0+<?=k?>] * ((eig)-><?=fluid?>_v.y - (eig)-><?=fluid?>_Cs * normal_u1y_over_len(n))\
			+ (X)->ptr[1+<?=k?>] * (eig)-><?=fluid?>_v.y\
			+ (X)->ptr[2+<?=k?>] * normal_u2y(n)\
			+ (X)->ptr[3+<?=k?>] * normal_u3y(n)\
			+ (X)->ptr[4+<?=k?>] * ((eig)-><?=fluid?>_v.y + (eig)-><?=fluid?>_Cs * normal_u1y_over_len(n));\
		(result)->ptr[3+<?=k?>] =\
			(X)->ptr[0+<?=k?>] * ((eig)-><?=fluid?>_v.z - (eig)-><?=fluid?>_Cs * normal_u1z_over_len(n))\
			+ (X)->ptr[1+<?=k?>] * (eig)-><?=fluid?>_v.z\
			+ (X)->ptr[2+<?=k?>] * normal_u2z(n)\
			+ (X)->ptr[3+<?=k?>] * normal_u3z(n)\
			+ (X)->ptr[4+<?=k?>] * ((eig)-><?=fluid?>_v.z + (eig)-><?=fluid?>_Cs * normal_u1z_over_len(n));\
		(result)->ptr[4+<?=k?>] =\
			(X)->ptr[0+<?=k?>] * ((eig)-><?=fluid?>_hTotal - (eig)-><?=fluid?>_Cs * v_n.x * inv_nLen)\
			+ (X)->ptr[1+<?=k?>] * .5 * (eig)-><?=fluid?>_vSq\
			+ (X)->ptr[2+<?=k?>] * v_n.y\
			+ (X)->ptr[3+<?=k?>] * v_n.z\
			+ (X)->ptr[4+<?=k?>] * ((eig)-><?=fluid?>_hTotal + (eig)-><?=fluid?>_Cs * v_n.x * inv_nLen);\
	}\
<? k=k+5 end ?>\
		/* EM & gravity */\
<? --\
for i,suffix in ipairs{"", "_g"} do --\
?>\
	{\
		real const tmp1 = sqrt_1_eps<?=suffix?> * n_l.x;\
		real const tmp3 = 1. / sqrt_1_eps<?=suffix?>;\
		real const tmp4 = tmp3 * n3_l.x;\
		real const tmp6 = n2_l.x * tmp3;\
		real const tmp22 = sqrt_1_eps<?=suffix?> * n_l.y;\
		real const tmp25 = n3_l.y * tmp3;\
		real const tmp27 = n2_l.y * tmp3;\
		real const tmp43 = sqrt_1_eps<?=suffix?> * n_l.z;\
		real const tmp46 = n3_l.z * tmp3;\
		real const tmp48 = n2_l.z * tmp3;\
		real const tmp66 = 1. / sqrt_1_mu<?=suffix?>;\
		real const tmp67 = tmp66 * n2_l.x;\
		real const tmp69 = n3_l.x * tmp66;\
		real const tmp88 = n2_l.y * tmp66;\
		real const tmp90 = n3_l.y * tmp66;\
		real const tmp109 = n2_l.z * tmp66;\
		real const tmp111 = n3_l.z * tmp66;\
		(result)->ptr[0+<?=k?>] = (X)->ptr[6+<?=k?>] * tmp1 + (X)->ptr[5+<?=k?>] * tmp6 + (X)->ptr[2+<?=k?>] * tmp4 - (X)->ptr[3+<?=k?>] * tmp6 - (X)->ptr[4+<?=k?>] * tmp4 + -(X)->ptr[0+<?=k?>] * tmp1;\
		(result)->ptr[1+<?=k?>] = (X)->ptr[6+<?=k?>] * tmp22 + (X)->ptr[5+<?=k?>] * tmp27 + (X)->ptr[2+<?=k?>] * tmp25 - (X)->ptr[3+<?=k?>] * tmp27 - (X)->ptr[4+<?=k?>] * tmp25 + -(X)->ptr[0+<?=k?>] * tmp22;\
		(result)->ptr[2+<?=k?>] = (X)->ptr[6+<?=k?>] * tmp43 + (X)->ptr[5+<?=k?>] * tmp48 + (X)->ptr[2+<?=k?>] * tmp46 - (X)->ptr[3+<?=k?>] * tmp48 - (X)->ptr[4+<?=k?>] * tmp46 + -(X)->ptr[0+<?=k?>] * tmp43;\
		(result)->ptr[3+<?=k?>] = (X)->ptr[7+<?=k?>] * tmp1 + (X)->ptr[5+<?=k?>] * tmp69 + (X)->ptr[4+<?=k?>] * tmp67 + (X)->ptr[3+<?=k?>] * tmp69 + (X)->ptr[2+<?=k?>] * tmp67 + -(X)->ptr[1+<?=k?>] * tmp1;\
		(result)->ptr[4+<?=k?>] = (X)->ptr[7+<?=k?>] * tmp22 + (X)->ptr[5+<?=k?>] * tmp90 + (X)->ptr[4+<?=k?>] * tmp88 + (X)->ptr[3+<?=k?>] * tmp90 + (X)->ptr[2+<?=k?>] * tmp88 + -(X)->ptr[1+<?=k?>] * tmp22;\
		(result)->ptr[5+<?=k?>] = (X)->ptr[7+<?=k?>] * tmp43 + (X)->ptr[5+<?=k?>] * tmp111 + (X)->ptr[4+<?=k?>] * tmp109 + (X)->ptr[3+<?=k?>] * tmp111 + (X)->ptr[2+<?=k?>] * tmp109 + -(X)->ptr[1+<?=k?>] * tmp43;\
		(result)->ptr[6+<?=k?>] = (X)->ptr[6+<?=k?>] * sqrt_1_eps<?=suffix?> + (X)->ptr[0+<?=k?>] * sqrt_1_eps<?=suffix?>;\
		(result)->ptr[7+<?=k?>] = (X)->ptr[7+<?=k?>] * sqrt_1_eps<?=suffix?> + (X)->ptr[1+<?=k?>] * sqrt_1_eps<?=suffix?>;\
	}\
<? k=k+8 end ?>\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: units

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
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
	(result)-><?=fluid?>_rho =\
		  (X)->ptr[<?=k+1?>] * nx\
		+ (X)->ptr[<?=k+2?>] * ny\
		+ (X)->ptr[<?=k+3?>] * nz;\
	(result)-><?=fluid?>_m.x =\
		  (X)->ptr[<?=k+0?>] * (-<?=fluid?>_v_n * (eig)-><?=fluid?>_v.x + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.x)\
		+ (X)->ptr[<?=k+1?>] * ((eig)-><?=fluid?>_v.x * nx - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)\
		+ (X)->ptr[<?=k+2?>] * ((eig)-><?=fluid?>_v.x * ny - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.y)\
		+ (X)->ptr[<?=k+3?>] * ((eig)-><?=fluid?>_v.x * nz - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.z)\
		+ (X)->ptr[<?=k+4?>] * (solver->heatCapacityRatio - 1.) * nx;\
	(result)-><?=fluid?>_m.y =\
		  (X)->ptr[<?=k+0?>] * (-<?=fluid?>_v_n * (eig)-><?=fluid?>_v.y + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.y)\
		+ (X)->ptr[<?=k+1?>] * ((eig)-><?=fluid?>_v.y * nx - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.x)\
		+ (X)->ptr[<?=k+2?>] * ((eig)-><?=fluid?>_v.y * ny - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)\
		+ (X)->ptr[<?=k+3?>] * ((eig)-><?=fluid?>_v.y * nz - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.z)\
		+ (X)->ptr[<?=k+4?>] * (solver->heatCapacityRatio - 1.) * ny;\
	(result)-><?=fluid?>_m.z =\
		  (X)->ptr[<?=k+0?>] * (-<?=fluid?>_v_n * (eig)-><?=fluid?>_v.z + (solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq * nU.z)\
		+ (X)->ptr[<?=k+1?>] * ((eig)-><?=fluid?>_v.z * nx - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.x)\
		+ (X)->ptr[<?=k+2?>] * ((eig)-><?=fluid?>_v.z * ny - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.y)\
		+ (X)->ptr[<?=k+3?>] * ((eig)-><?=fluid?>_v.z * nz - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)\
		+ (X)->ptr[<?=k+4?>] * (solver->heatCapacityRatio - 1.) * nz;\
	(result)-><?=fluid?>_ETotal =\
		  (X)->ptr[<?=k+0?>] * <?=fluid?>_v_n * ((solver->heatCapacityRatio - 1.) * .5 * (eig)-><?=fluid?>_vSq - (eig)-><?=fluid?>_hTotal)\
		+ (X)->ptr[<?=k+1?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.x + nx * (eig)-><?=fluid?>_hTotal)\
		+ (X)->ptr[<?=k+2?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.y + ny * (eig)-><?=fluid?>_hTotal)\
		+ (X)->ptr[<?=k+3?>] * (-(solver->heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.z + nz * (eig)-><?=fluid?>_hTotal)\
		+ (X)->ptr[<?=k+4?>] * solver->heatCapacityRatio * <?=fluid?>_v_n;\
<? k=k+5 end ?>\
\
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
<? for _,suffix in ipairs{"", "_g"} do ?>\
	{\
		real3 const E = real3_real_mul((X)->D<?=suffix?>, 1. / eps<?=suffix?>);\
		real3 const H = real3_real_mul((X)->B<?=suffix?>, 1. / mu<?=suffix?>);\
\
		(result)->D<?=suffix?>.x = H.y * nz - H.z * ny + nx * (U)->phi * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s;	/* F_D^i = -eps^ijk n_j H_k */\
		(result)->B<?=suffix?>.x = E.z * ny - E.y * nz + nx * (U)->psi * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s;	/* F_B^i = +eps^ijk n_j B_k */\
\
		(result)->D<?=suffix?>.y = H.z * nx - H.x * nz + ny * (U)->phi * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s;\
		(result)->B<?=suffix?>.y = E.x * nz - E.z * nx + ny * (U)->psi * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s;\
\
		(result)->D<?=suffix?>.z = H.x * ny - H.y * nx + nz * (U)->phi * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s;\
		(result)->B<?=suffix?>.z = E.y * nx - E.x * ny + nz * (U)->psi * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s;\
\
		real const D_n = normal_vecDotN1(n, (U)->D<?=suffix?>);\
		real const B_n = normal_vecDotN1(n, (U)->B<?=suffix?>);\
		(result)->phi<?=suffix?> = <?=real_mul?>(D_n, solver->divPhiWavespeed<?=suffix?> / unit_m_per_s);\
		(result)->psi<?=suffix?> = <?=real_mul?>(B_n, solver->divPsiWavespeed<?=suffix?> / unit_m_per_s);\
	}\
<? end ?>\
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


<? if not require "hydro.coord.cartesian":isa(solver.coord) then ?>
	/* connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system */
//// MODULE_DEPENDS: <?=coord_conn_trace13?>
//// MODULE_DEPENDS: <?=coord_raise?>
//// MODULE_DEPENDS: <?=coord_conn_apply23?>
//// MODULE_DEPENDS: <?=coord_conn_trace23?>
//// MODULE_DEPENDS: <?=coord_conn_apply13?>
//// MODULE_DEPENDS: <?=coord_conn_apply123?>
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
