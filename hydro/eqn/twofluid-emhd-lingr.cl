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

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=solver.solver_t?> solver_t;

<? if moduleName == nil then ?>
<? elseif moduleName == "sqrt_2_and_1_2" then ?>

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>
#define sqrt_2 <?=('%.50f'):format(math.sqrt(2))?>

<? elseif moduleName == "eqn.prim-cons" then ?>

<?=eqn.prim_t?> primFromCons(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
	<? for _,fluid in ipairs(eqn.fluids) do ?>
	real <?=fluid?>_EKin = calc_<?=fluid?>_EKin_fromCons(U, x);
	real <?=fluid?>_EInt = U.<?=fluid?>_ETotal - <?=fluid?>_EKin;
	<? end ?>
	return (<?=eqn.prim_t?>){
		<? for _,fluid in ipairs(eqn.fluids) do ?>
		.<?=fluid?>_rho = U.<?=fluid?>_rho,
		.<?=fluid?>_v = real3_real_mul(U.<?=fluid?>_m, 1./U.<?=fluid?>_rho),
		.<?=fluid?>_P = (solver->heatCapacityRatio - 1.) * <?=fluid?>_EInt,
		<? end ?>
		.D = U.D,
		.B = U.B,
		.psi = U.psi,
		.phi = U.phi,
		.D_g = U.D_g,
		.B_g = U.B_g,
		.psi_g = U.psi_g,
		.phi_g = U.phi_g,
	};
}

<?=eqn.cons_t?> consFromPrim(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
<? for _,fluid in ipairs(eqn.fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_m = real3_real_mul(W.<?=fluid?>_v, W.<?=fluid?>_rho),
		.<?=fluid?>_ETotal = calc_<?=fluid?>_ETotal(solver, W, x),
<? end ?>
		.D = W.D,
		.B = W.B,
		.psi = W.psi,
		.phi = W.phi,
		.D_g = W.D_g,
		.B_g = W.B_g,
		.psi_g = W.psi_g,
		.phi_g = W.phi_g,
	};
}

<? elseif moduleName == "eqn.dU-dW" then ?>

<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
<? for _,fluid in ipairs(eqn.fluids) do ?>
	real3 WA_<?=fluid?>_vL = coord_lower(WA.<?=fluid?>_v, x);
<? end ?>
	return (<?=eqn.cons_t?>){
<? for _,fluid in ipairs(eqn.fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_m = real3_add(
			real3_real_mul(WA.<?=fluid?>_v, W.<?=fluid?>_rho), 
			real3_real_mul(W.<?=fluid?>_v, WA.<?=fluid?>_rho)),
		.<?=fluid?>_ETotal = W.<?=fluid?>_rho * .5 * real3_dot(WA.<?=fluid?>_v, WA_<?=fluid?>_vL) 
			+ WA.<?=fluid?>_rho * real3_dot(W.<?=fluid?>_v, WA_<?=fluid?>_vL)
			+ W.<?=fluid?>_P / (solver->heatCapacityRatio - 1.),
<? end ?>
		.B = W.B,
		.D = W.D,
		.phi = W.phi,
		.psi = W.psi,
		.B_g = W.B_g,
		.D_g = W.D_g,
		.phi_g = W.phi_g,
		.psi_g = W.psi_g,
	};
}

<?=eqn.prim_t?> apply_dW_dU(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
<? for _,fluid in ipairs(eqn.fluids) do ?>
	real3 WA_<?=fluid?>_vL = coord_lower(WA.<?=fluid?>_v, x);
<? end ?>
	return (<?=eqn.prim_t?>){
<? for _,fluid in ipairs(eqn.fluids) do ?>
		.<?=fluid?>_rho = U.<?=fluid?>_rho,
		.<?=fluid?>_v = real3_sub(
			real3_real_mul(U.<?=fluid?>_m, 1. / WA.<?=fluid?>_rho),
			real3_real_mul(WA.<?=fluid?>_v, U.<?=fluid?>_rho / WA.<?=fluid?>_rho)),
		.<?=fluid?>_P = (solver->heatCapacityRatio - 1.) * (
			.5 * real3_dot(WA.<?=fluid?>_v, WA_<?=fluid?>_vL) * U.<?=fluid?>_rho 
			- real3_dot(U.<?=fluid?>_m, WA_<?=fluid?>_vL)
			+ U.<?=fluid?>_ETotal),
<? end ?>
		.B = U.B,
		.D = U.D,
		.phi = U.phi,
		.psi = U.psi,
		.B_g = U.B_g,
		.D_g = U.D_g,
		.phi_g = U.phi_g,
		.psi_g = U.psi_g,
	};
}

<? elseif moduleName == "eqn.common" then ?>

real3 calc_EField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	return real3_real_mul(U.D, 1. / eps);
}
 
real3 calc_HField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) { 
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return real3_real_mul(U.B, 1. / mu);
}

real3 calc_SField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	return real3_cross(
		calc_EField(solver, U),
		calc_HField(solver, U));
}

real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
real calc_hTotal(constant <?=solver.solver_t?>* solver, real rho, real P, real ETotal) { return (P + ETotal) / rho; }
real calc_HTotal(real P, real ETotal) { return P + ETotal; }

<? for _,fluid in ipairs(eqn.fluids) do ?>
real calc_<?=fluid?>_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.<?=fluid?>_v, x); }
real calc_<?=fluid?>_EKin(<?=eqn.prim_t?> W, real3 x) { return W.<?=fluid?>_rho * calc_<?=fluid?>_eKin(W, x); }
real calc_<?=fluid?>_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.<?=fluid?>_P / (solver->heatCapacityRatio - 1.); }
real calc_<?=fluid?>_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_<?=fluid?>_EInt(solver, W) / W.<?=fluid?>_rho; }
real calc_<?=fluid?>_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.<?=fluid?>_m, x) / U.<?=fluid?>_rho; }
real calc_<?=fluid?>_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return calc_<?=fluid?>_EKin(W, x) + calc_<?=fluid?>_EInt(solver, W);
}
real calc_<?=fluid?>_Cs(constant <?=solver.solver_t?>* solver, const <?=eqn.prim_t?>* W) {
	return sqrt(solver->heatCapacityRatio * W-><?=fluid?>_P / W-><?=fluid?>_rho);
}
<? end ?>

real calc_EM_energy(constant <?=solver.solver_t?>* solver, const global <?=eqn.cons_t?>* U, real3 x) {
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
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
real3 calcIonGravForce(constant <?=solver.solver_t?>* solver, const global <?=eqn.cons_t?>* U, real3 x) {
	const real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	const real eps_g = 1. / (4. * M_PI * G);
	return _real3(
		U->ion_rho * U->D_g.x / eps_g + 4. * (U->ion_m.y * U->B_g.z - U->ion_m.z * U->B_g.y),
		U->ion_rho * U->D_g.y / eps_g + 4. * (U->ion_m.z * U->B_g.x - U->ion_m.x * U->B_g.z),
		U->ion_rho * U->D_g.z / eps_g + 4. * (U->ion_m.x * U->B_g.y - U->ion_m.y * U->B_g.x));
}

real3 calcElecGravForce(constant <?=solver.solver_t?>* solver, const global <?=eqn.cons_t?>* U, real3 x) {
	const real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	const real eps_g = 1. / (4. * M_PI * G);
	return _real3(
		U->elec_rho * U->D_g.x / eps_g + 4. * (U->elec_m.y * U->B_g.z - U->elec_m.z * U->B_g.y),
		U->elec_rho * U->D_g.y / eps_g + 4. * (U->elec_m.z * U->B_g.x - U->elec_m.x * U->B_g.z),
		U->elec_rho * U->D_g.z / eps_g + 4. * (U->elec_m.x * U->B_g.y - U->elec_m.y * U->B_g.x));
}

<? elseif moduleName == "applyInitCond" then ?>

<? 
local cons_t = eqn.cons_t
local susc_t = eqn.susc_t
local scalar = eqn.scalar
local vec3 = eqn.vec3
local zero = scalar..'_zero'
local inv = scalar..'_inv'
local fromreal = scalar..'_from_real'
local sqrt = scalar..'_sqrt'
?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = x.x < mids.x
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
	 for _,fluid in ipairs(eqn.fluids) do
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
	<?=eqn.prim_t?> W;
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
	for _,fluid in ipairs(eqn.fluids) do ?>
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
	UBuf[index] = consFromPrim(solver, W, x);
}


<? elseif moduleName == "fluxFromCons" then ?>

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	<?=eqn.cons_t?> F;

<? 
for _,fluid in ipairs(eqn.fluids) do
?>	real <?=fluid?>_vj = normal_vecDotN1(n, W.<?=fluid?>_v);
	real <?=fluid?>_HTotal = U.<?=fluid?>_ETotal + W.<?=fluid?>_P;
	
	F.<?=fluid?>_rho = normal_vecDotN1(n, U.<?=fluid?>_m);
	F.<?=fluid?>_m = real3_real_mul(U.<?=fluid?>_m, <?=fluid?>_vj);
<? 	for i,xi in ipairs(xNames) do
?>	F.<?=fluid?>_m.<?=xi?> += normal_u1<?=xi?>(n) * W.<?=fluid?>_P;
<? 	end
?>	F.<?=fluid?>_ETotal = <?=fluid?>_HTotal * <?=fluid?>_vj;
<? 
end
?>
	
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real eps_g = 1. / (4. * M_PI * G);
	real mu_g = 1. / (eps_g * speedOfLightSq);

	//taken from glm-maxwell instead of the 2014 Abgrall, Kumar
	// then replace D = epsilon E and phi' -> epsilon phi
	<? for _,suffix in ipairs{'', '_g'} do ?>{
		real3 E = real3_real_mul(U.D<?=suffix?>, 1. / eps<?=suffix?>);
		real3 H = real3_real_mul(U.B<?=suffix?>, 1. / mu<?=suffix?>);
		if (n.side == 0) {
			F.D<?=suffix?> = _real3(U.phi<?=suffix?> * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s, H.z, -H.y);
			F.B<?=suffix?> = _real3(U.psi<?=suffix?> * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s, -E.z, E.y);
		} else if (n.side == 1) {
			F.D<?=suffix?> = _real3(-H.z, U.phi<?=suffix?> * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s, H.x);
			F.B<?=suffix?> = _real3(E.z, U.psi<?=suffix?> * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s, -E.x);
		} else if (n.side == 2) {
			F.D<?=suffix?> = _real3(H.y, -H.x, U.phi<?=suffix?> * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s);
			F.B<?=suffix?> = _real3(-E.y, E.x, U.psi<?=suffix?> * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s);
		}
		F.phi<?=suffix?> = normal_vecDotN1(n, U.D<?=suffix?>) * solver->divPhiWavespeed<?=suffix?> / unit_m_per_s;
		F.psi<?=suffix?> = normal_vecDotN1(n, U.B<?=suffix?>) * solver->divPsiWavespeed<?=suffix?> / unit_m_per_s;
	}<? end ?>

	return F;
}

<? elseif moduleName == "eigen_forInterface" then ?>

typedef <?=eqn.eigen_t?> eigen_t;

eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	normal_t n
) {
	prim_t WL = primFromCons(solver, UL, x);
	prim_t WR = primFromCons(solver, UR, x);
	eigen_t eig;

<? for _,fluid in ipairs(eqn.fluids) do ?>

	real <?=fluid?>_sqrtRhoL = sqrt(WL.<?=fluid?>_rho);
	real3 <?=fluid?>_vL = WL.<?=fluid?>_v;
	real <?=fluid?>_hTotalL = calc_hTotal(solver, WL.<?=fluid?>_rho, WL.<?=fluid?>_P, UL.<?=fluid?>_ETotal);

	real <?=fluid?>_sqrtRhoR = sqrt(UR.<?=fluid?>_rho);
	real3 <?=fluid?>_vR = WR.<?=fluid?>_v;
	real <?=fluid?>_hTotalR = calc_hTotal(solver, WR.<?=fluid?>_rho, WR.<?=fluid?>_P, UR.<?=fluid?>_ETotal);

	real <?=fluid?>_invDenom = 1./(<?=fluid?>_sqrtRhoL + <?=fluid?>_sqrtRhoR);
	
	//Roe-averaged
	real <?=fluid?>_rho = <?=fluid?>_sqrtRhoL * <?=fluid?>_sqrtRhoR;
	real3 <?=fluid?>_v = real3_add(
			real3_real_mul(<?=fluid?>_vL, <?=fluid?>_sqrtRhoL * <?=fluid?>_invDenom),
			real3_real_mul(<?=fluid?>_vR, <?=fluid?>_sqrtRhoR * <?=fluid?>_invDenom));
	real <?=fluid?>_hTotal = <?=fluid?>_invDenom * (<?=fluid?>_sqrtRhoL * <?=fluid?>_hTotalL + <?=fluid?>_sqrtRhoR * <?=fluid?>_hTotalR);

	//derived:
	real <?=fluid?>_vSq = coordLenSq(<?=fluid?>_v, x);
	real <?=fluid?>_eKin = .5 * <?=fluid?>_vSq;
	real <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * (<?=fluid?>_hTotal - <?=fluid?>_eKin);
	real <?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);
	
	eig.<?=fluid?>_rho = <?=fluid?>_rho; 
	eig.<?=fluid?>_v = <?=fluid?>_v;
	eig.<?=fluid?>_hTotal = <?=fluid?>_hTotal;
	eig.<?=fluid?>_vSq = <?=fluid?>_vSq;
	eig.<?=fluid?>_Cs = <?=fluid?>_Cs;

<? end ?>
	
	return eig;
}

<? elseif moduleName == "eigen_forCell" then ?>

typedef <?=eqn.eigen_t?> eigen_t;

eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	normal_t n
) {
	prim_t W = primFromCons(solver, U, x);
<? for _,fluid in ipairs(eqn.fluids) do ?>
	real <?=fluid?>_vSq = coordLenSq(W.<?=fluid?>_v, x);
	real <?=fluid?>_eKin = .5 * <?=fluid?>_vSq;
	real <?=fluid?>_hTotal = calc_hTotal(solver, W.<?=fluid?>_rho, W.<?=fluid?>_P, U.<?=fluid?>_ETotal);
	real <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * (<?=fluid?>_hTotal - <?=fluid?>_eKin);
	real <?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);
<? end ?>	
	return (eigen_t){
<? for _,fluid in ipairs(eqn.fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_v = W.<?=fluid?>_v,
		.<?=fluid?>_hTotal = <?=fluid?>_hTotal,
		.<?=fluid?>_vSq = <?=fluid?>_vSq,
		.<?=fluid?>_Cs = <?=fluid?>_Cs,
<? end ?>	
	};
}

<? elseif moduleName == "eigen_left/rightTransform" then ?>

typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;

<?
local prefix = [[
	real nLen = normal_len(n);
	real nLenSq = nLen * nLen;

	//g^ij for fixed j=side
	real3 ion_v = eig.ion_v;
	real3 ion_vL = coord_lower(ion_v, x);
	real ion_hTotal = eig.ion_hTotal;
	real ion_vSq = real3_dot(ion_v, ion_vL);
	real ion_Cs = eig.ion_Cs;
	real ion_Cs_over_nLen = ion_Cs / nLen; 
	
	real3 elec_v = eig.elec_v;
	real3 elec_vL = coord_lower(elec_v, x);
	real elec_hTotal = eig.elec_hTotal;
	real elec_vSq = real3_dot(elec_v, elec_vL);
	real elec_Cs = eig.elec_Cs;
	real elec_Cs_over_nLen = elec_Cs / nLen; 

	real nx = normal_l1x(n);
	real ny = normal_l1y(n);
	real nz = normal_l1z(n);
	real n1x = normal_l2x(n);
	real n1y = normal_l2y(n);
	real n1z = normal_l2z(n);
	real n2x = normal_l3x(n);
	real n2y = normal_l3y(n);
	real n2z = normal_l3z(n);

	real3 nU = normal_u1(n);

	real3 ion_v_ns = normal_vecDotNs(n, ion_v);
	real ion_v_n = ion_v_ns.x, ion_v_n1 = ion_v_ns.y, ion_v_n2 = ion_v_ns.z;
	real3 elec_v_ns = normal_vecDotNs(n, elec_v);
	real elec_v_n = elec_v_ns.x, elec_v_n1 = elec_v_ns.y, elec_v_n2 = elec_v_ns.z;
 
]]
?>

waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t UX,
	real3 x,
	normal_t n
) { 
	waves_t UY;
	real* Y = UY.ptr;
	real* X = UX.ptr;

	<?=prefix?>

	real ion_denom = 2. * ion_Cs * ion_Cs;
	real ion_invDenom = 1. / ion_denom;
	
	real elec_denom = 2. * elec_Cs * elec_Cs;
	real elec_invDenom = 1. / elec_denom;

	const real heatRatioMinusOne = solver->heatCapacityRatio - 1.;

	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	real sqrt_eps = sqrt(eps);	// TODO sqrt units
	real sqrt_mu = sqrt(mu);
	real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real eps_g = 1. / (4. * M_PI * G);
	real mu_g = 1. / (eps_g * speedOfLightSq);
	real sqrt_eps_g = sqrt(eps_g);	// TODO sqrt units
	real sqrt_mu_g = sqrt(mu_g);

	if (n.side == 0) {
<?
					for i,fluid in ipairs(eqn.fluids) do
?>
		Y[<?=5*i-5?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.x / nLen)
			+ X[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x - <?=fluid?>_Cs / nLen)
			+ X[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
		Y[<?=5*i-4?>] = (X[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)
			+ X[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * -2. * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
		Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x * nU.y / nLenSq - <?=fluid?>_v.y)
			+ X[<?=5*i-4?>] * -nU.y / nLenSq
			+ X[<?=5*i-3?>];
		Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x * nU.z / nLenSq - <?=fluid?>_v.z)
			+ X[<?=5*i-4?>] * -nU.z / nLenSq
			+ X[<?=5*i-2?>];
		Y[<?=5*i-1?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.x / nLen)
			+ X[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x + <?=fluid?>_Cs / nLen)
			+ X[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
<? 
					end
?>	
		X += 2;
		Y += 2;
		//EM & gravity
		<? for i,suffix in ipairs{'', '_g'} do ?>{
			X += 8;
			Y += 8;
			Y[0] = ((-(sqrt_eps<?=suffix?> * (X[0] - X[6]))) / sqrt_2);
			Y[1] = ((-(sqrt_eps<?=suffix?> * (X[3] - X[7]))) / sqrt_2);
			Y[2] = (((X[2] * sqrt_mu<?=suffix?>) + (X[4] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>));
			Y[3] = (((X[1] * sqrt_mu<?=suffix?>) - (X[5] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));
			Y[4] = ((-((X[2] * sqrt_mu<?=suffix?>) - (X[4] * sqrt_eps<?=suffix?>))) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));
			Y[5] = (((X[1] * sqrt_mu<?=suffix?>) + (X[5] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));
			Y[6] = ((sqrt_eps<?=suffix?> * (X[0] + X[6])) / sqrt_2);
			Y[7] = ((sqrt_eps<?=suffix?> * (X[3] + X[7])) / sqrt_2);
		}<? end ?>

	} else if (n.side == 1) {
<?	
					for i,fluid in ipairs(eqn.fluids) do
?>
		Y[<?=5*i-5?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.y / nLen)
			+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y - <?=fluid?>_Cs / nLen)
			+ X[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
		Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y * nU.x / nLenSq - <?=fluid?>_v.x)
			+ X[<?=5*i-4?>]
			+ X[<?=5*i-3?>] * -nU.x / nLenSq;
		Y[<?=5*i-3?>] = (X[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)
			+ X[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * -2. * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
		Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y * nU.z / nLenSq - <?=fluid?>_v.z)
			+ X[<?=5*i-3?>] * -nU.z / nLenSq
			+ X[<?=5*i-2?>];
		Y[<?=5*i-1?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.y / nLen)
			+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y + <?=fluid?>_Cs / nLen)
			+ X[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
<?
					end 
?>
		X += 2;
		Y += 2;
		//EM & gravity
		<? for i,suffix in ipairs{'', '_g'} do ?>{
			X += 8;
			Y += 8;
			Y[0] = ((-(sqrt_eps<?=suffix?> * (X[1] - X[6]))) / sqrt_2);
			Y[1] = ((-(sqrt_eps<?=suffix?> * (X[4] - X[7]))) / sqrt_2);
			Y[2] = (((X[2] * sqrt_mu<?=suffix?>) - (X[3] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));
			Y[3] = (((X[0] * sqrt_mu<?=suffix?>) + (X[5] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>));
			Y[4] = (((X[2] * sqrt_mu<?=suffix?>) + (X[3] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));
			Y[5] = ((-((X[0] * sqrt_mu<?=suffix?>) - (X[5] * sqrt_eps<?=suffix?>))) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));
			Y[6] = ((sqrt_eps<?=suffix?> * (X[1] + X[6])) / sqrt_2);
			Y[7] = ((sqrt_eps<?=suffix?> * (X[4] + X[7])) / sqrt_2);
		}<? end ?>

	} else if (n.side == 2) {
<?
					for i,fluid in ipairs(eqn.fluids) do
?>
		Y[<?=5*i-5?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.z / nLen)
			+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z - <?=fluid?>_Cs / nLen)
			+ X[<?=5*i-1?>] * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
		Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z * nU.x / nLenSq - <?=fluid?>_v.x)
			+ X[<?=5*i-4?>]
			+ X[<?=5*i-2?>] * -nU.x / nLenSq;
		Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z * nU.y / nLenSq - <?=fluid?>_v.y)
			+ X[<?=5*i-3?>]
			+ X[<?=5*i-2?>] * -nU.y / nLenSq;
		Y[<?=5*i-2?>] = (X[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)
			+ X[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * -2. * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
		Y[<?=5*i-1?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.z / nLen)
			+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z + <?=fluid?>_Cs / nLen)
			+ X[<?=5*i-1?>] * heatRatioMinusOne
		) * <?=fluid?>_invDenom;
<? 
					end
?>
		X += 2;
		Y += 2;
		//EM & gravity
		<? for i,suffix in ipairs{'', '_g'} do ?>{
			X += 8;
			Y += 8;
			Y[0] = ((-(sqrt_eps<?=suffix?> * (X[2] - X[6]))) / sqrt_2);
			Y[1] = ((-(sqrt_eps<?=suffix?> * (X[5] - X[7]))) / sqrt_2);
			Y[2] = (((X[1] * sqrt_mu<?=suffix?>) + (X[3] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>));
			Y[3] = (((X[0] * sqrt_mu<?=suffix?>) - (X[4] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));
			Y[4] = (((X[1] * sqrt_mu<?=suffix?>) - (X[3] * sqrt_eps<?=suffix?>)) / (-(sqrt_mu<?=suffix?> * sqrt_2 * sqrt_eps<?=suffix?>)));
			Y[5] = (((X[0] * sqrt_mu<?=suffix?>) + (X[4] * sqrt_eps<?=suffix?>)) / (sqrt_mu<?=suffix?> * sqrt_eps<?=suffix?> * sqrt_2));
			Y[6] = ((sqrt_eps<?=suffix?> * (X[2] + X[6])) / sqrt_2);
			Y[7] = ((sqrt_eps<?=suffix?> * (X[5] + X[7])) / sqrt_2);
		}<? end ?>

	}

	return UY;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t UX,	//numWaves = 26
	real3 x,
	normal_t n
) {
	cons_t UY;
	real* Y = UY.ptr;
	real* X = UX.ptr;

	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	real sqrt_eps = sqrt(eps);	// TODO sqrt units
	real sqrt_mu = sqrt(mu);
	real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real eps_g = 1. / (4. * M_PI * G);
	real mu_g = 1. / (eps_g * speedOfLightSq);
	real sqrt_eps_g = sqrt(eps_g);	// TODO sqrt units
	real sqrt_mu_g = sqrt(mu_g);

	<?=prefix?>
	
	
	if (n.side == 0) {
<?
					for i,fluid in ipairs(eqn.fluids) do
?>
		Y[<?=5*i-5?>] = X[<?=5*i-5?>] + X[<?=5*i-4?>] + X[<?=5*i-1?>];
		Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * nLen)
			+ X[<?=5*i-4?>] * <?=fluid?>_v.x
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * nLen);
		Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * nU.y / nLen)
			+ X[<?=5*i-4?>] * <?=fluid?>_v.y
			+ X[<?=5*i-3?>]
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * nU.y / nLen);
		Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * nU.z / nLen)
			+ X[<?=5*i-4?>] * <?=fluid?>_v.z
			+ X[<?=5*i-2?>]
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * nU.z / nLen);
		Y[<?=5*i-1?>] = X[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.x / nLen)
			+ X[<?=5*i-4?>] * <?=fluid?>_vSq / 2.
			+ X[<?=5*i-3?>] * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.x / nLen);
<?
					end
?>
		X += 2;
		Y += 2;
		//EM & gravity
		<? for i,suffix in ipairs{'', '_g'} do ?>{
			X += 8;
			Y += 8;
			Y[0] = ((-(X[0] - X[6])) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[1] = ((-(sqrt_eps<?=suffix?> * (X[3] - X[5]))) / sqrt_2);
			Y[2] = ((sqrt_eps<?=suffix?> * (X[2] - X[4])) / sqrt_2);
			Y[3] = ((-(X[1] - X[7])) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[4] = ((sqrt_mu<?=suffix?> * (X[2] + X[4])) / sqrt_2);
			Y[5] = ((sqrt_mu<?=suffix?> * (X[3] + X[5])) / sqrt_2);
			Y[6] = ((X[0] + X[6]) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[7] = ((X[1] + X[7]) / (sqrt_2 * sqrt_eps<?=suffix?>));
		}<? end ?>

	} else if (n.side == 1) {
<?
					for i,fluid in ipairs(eqn.fluids) do
?>	
		Y[<?=5*i-5?>] = X[<?=5*i-5?>] + X[<?=5*i-3?>] + X[<?=5*i-1?>];
		Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * nU.x / nLen)
			+ X[<?=5*i-4?>]
			+ X[<?=5*i-3?>] * <?=fluid?>_v.x
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * nU.x / nLen);
		Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * nLen)
			+ X[<?=5*i-3?>] * <?=fluid?>_v.y
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * nLen);
		Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * nU.z / nLen)
			+ X[<?=5*i-3?>] * <?=fluid?>_v.z
			+ X[<?=5*i-2?>]
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * nU.z / nLen);
		Y[<?=5*i-1?>] = X[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.y / nLen)
			+ X[<?=5*i-4?>] * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * <?=fluid?>_vSq / 2.
			+ X[<?=5*i-2?>] * <?=fluid?>_vL.z
			+ X[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.y / nLen);
<?
					end
?>
		X += 2;
		Y += 2;
		//EM & gravity
		<? for i,suffix in ipairs{'', '_g'} do ?>{
			X += 8;
			Y += 8;
			Y[0] = ((sqrt_eps<?=suffix?> * (X[3] - X[5])) / sqrt_2);
			Y[1] = ((-(X[0] - X[6])) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[2] = ((-(sqrt_eps<?=suffix?> * (X[2] - X[4]))) / sqrt_2);
			Y[3] = ((sqrt_mu<?=suffix?> * (X[2] + X[4])) / sqrt_2);
			Y[4] = ((-(X[1] - X[7])) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[5] = ((sqrt_mu<?=suffix?> * (X[3] + X[5])) / sqrt_2);
			Y[6] = ((X[0] + X[6]) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[7] = ((X[1] + X[7]) / (sqrt_2 * sqrt_eps<?=suffix?>));
		}<? end ?>

	} else if (n.side == 2) {
<?
					for i,fluid in ipairs(eqn.fluids) do
?>
		Y[<?=5*i-5?>] = X[<?=5*i-5?>] + X[<?=5*i-2?>] + X[<?=5*i-1?>];
		Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * nU.x / nLen)
			+ X[<?=5*i-4?>]
			+ X[<?=5*i-2?>] * <?=fluid?>_v.x
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * nU.x / nLen);
		Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * nU.y / nLen)
			+ X[<?=5*i-3?>]
			+ X[<?=5*i-2?>] * <?=fluid?>_v.y
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * nU.y / nLen);
		Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * nLen)
			+ X[<?=5*i-2?>] * <?=fluid?>_v.z
			+ X[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * nLen);
		Y[<?=5*i-1?>] = X[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.z / nLen)
			+ X[<?=5*i-4?>] * <?=fluid?>_vL.x
			+ X[<?=5*i-3?>] * <?=fluid?>_vL.y
			+ X[<?=5*i-2?>] * <?=fluid?>_vSq / 2.
			+ X[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.z / nLen);
<?
					end
?>
		X += 2;
		Y += 2;
		//EM & gravity
		<? for i,suffix in ipairs{'', '_g'} do ?>{
			X += 8;
			Y += 8;
			Y[0] = ((-(sqrt_eps<?=suffix?> * (X[3] - X[5]))) / sqrt_2);
			Y[1] = ((sqrt_eps<?=suffix?> * (X[2] - X[4])) / sqrt_2);
			Y[2] = ((-(X[0] - X[6])) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[3] = ((sqrt_mu<?=suffix?> * (X[2] + X[4])) / sqrt_2);
			Y[4] = ((sqrt_mu<?=suffix?> * (X[3] + X[5])) / sqrt_2);
			Y[5] = ((-(X[1] - X[7])) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[6] = ((X[0] + X[6]) / (sqrt_2 * sqrt_eps<?=suffix?>));
			Y[7] = ((X[1] + X[7]) / (sqrt_2 * sqrt_eps<?=suffix?>));
		}<? end ?>

	}

	return UY;
}

<? elseif moduleName == "eigen_fluxTransform" then ?>

typedef <?=eqn.eigen_t?> eigen_t;

cons_t eigen_fluxTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t UX,
	real3 x,
	normal_t n
) {
	<?=prefix?>
	cons_t UY;
	real* X = UX.ptr;
<?
					for i,fluid	in ipairs(eqn.fluids) do 
?>
	UY.<?=fluid?>_rho = X[<?=5*i-4?>] * nx 
		+ X[<?=5*i-3?>] * ny 
		+ X[<?=5*i-2?>] * nz;
	UY.<?=fluid?>_m.x = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.x + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * nU.x)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.x * nx - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.x * ny - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.x * nz - (solver->heatCapacityRatio - 1.) * nU.x * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nx;
	UY.<?=fluid?>_m.y = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.y + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * nU.y)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.y * nx - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.y * ny - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.y * nz - (solver->heatCapacityRatio - 1.) * nU.y * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * ny;
	UY.<?=fluid?>_m.z = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.z + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * nU.z)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.z * nx - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.z * ny - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.z * nz - (solver->heatCapacityRatio - 1.) * nU.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nz;
	UY.<?=fluid?>_ETotal = X[<?=5*i-5?>] * <?=fluid?>_v_n * ((solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq - <?=fluid?>_hTotal)
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
		UY.D = _real3(solver->divPhiWavespeed / unit_m_per_s * UX.phi, H.z, -H.y);
		UY.B = _real3(solver->divPsiWavespeed / unit_m_per_s * UX.psi, -E.z, E.y);
	} else if (n.side == 1) {
		UY.D = _real3(-H.z, solver->divPhiWavespeed / unit_m_per_s * UX.phi, H.x);
		UY.B = _real3(E.z, solver->divPsiWavespeed / unit_m_per_s * UX.psi, -E.x);
	} else if (n.side == 2) {
		UY.D = _real3(H.y, -H.x, solver->divPhiWavespeed / unit_m_per_s * UX.phi);
		UY.B = _real3(-E.y, E.x, solver->divPsiWavespeed / unit_m_per_s * UX.psi);
	}
	UY.phi = solver->divPhiWavespeed / unit_m_per_s * normal_vecDotN1(n, UX.D);
	UY.psi = solver->divPsiWavespeed / unit_m_per_s * normal_vecDotN1(n, UX.B);

	real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real eps_g = 1. / (4. * M_PI * G);
	real mu_g = 1. / (eps_g * speedOfLightSq);
	
	real3 E_g = real3_real_mul(UX.D_g, 1. / eps_g);
	real3 H_g = real3_real_mul(UX.B_g, 1. / mu_g);
	if (n.side == 0) {
		UY.D_g = _real3(solver->divPhiWavespeed_g / unit_m_per_s * UX.phi_g, H_g.z, -H_g.y);
		UY.B_g = _real3(solver->divPsiWavespeed_g / unit_m_per_s * UX.psi_g, -E_g.z, E_g.y);
	} else if (n.side == 1) {
		UY.D_g = _real3(-H_g.z, solver->divPhiWavespeed_g / unit_m_per_s * UX.phi_g, H_g.x);
		UY.B_g = _real3(E_g.z, solver->divPsiWavespeed_g / unit_m_per_s * UX.psi_g, -E_g.x);
	} else if (n.side == 2) {
		UY.D_g = _real3(H_g.y, -H_g.x, solver->divPhiWavespeed_g / unit_m_per_s * UX.phi_g);
		UY.B_g = _real3(-E_g.y, E_g.x, solver->divPsiWavespeed_g / unit_m_per_s * UX.psi_g);
	}
	UY.phi_g = solver->divPhiWavespeed_g / unit_m_per_s * normal_vecDotN1(n, E_g);
	UY.psi_g = solver->divPsiWavespeed_g / unit_m_per_s * normal_vecDotN1(n, H_g);

	return UY;
}

<? elseif moduleName == "addSource" then ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cellBuf[index].pos;
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	//TODO double check all units

	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	real G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real eps_g = 1. / (4. * M_PI * G);
	real mu_g = 1. / (eps_g * speedOfLightSq);

	//rho * (E + v * B) has units kg/(m^2 s^2)
	real3 ionGravForce = calcIonGravForce(solver, U, x);

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

	//source of D is the current
	deriv->D.x -= J.x;
	deriv->D.y -= J.y;
	deriv->D.z -= J.z;
	
	//source of phi is the charge
	deriv->phi += eps * (U->ion_rho * solver->ionChargeMassRatio + U->elec_rho * elecChargeMassRatio) / unit_C_per_kg * solver->divPhiWavespeed / unit_m_per_s;

	//source of D_g is J_g is the momentum + Poynting vector
	real3 J_g = real3_zero;
	// I'm symmetrizing the stress-energy
		//matter
	J_g.x -= U->ion_m.x + U->elec_m.x;
	J_g.y -= U->ion_m.y + U->elec_m.y;
	J_g.z -= U->ion_m.z + U->elec_m.z;
		//electromagnetism
	real v_p = 1. / (solver->sqrt_eps * solver->sqrt_mu);		// phase vel = sqrt(1 / (eps * mu))
	real v_pSq = v_p * v_p;
	real nSq = 1. / v_pSq;		// index of refraction^2 = 1 / phase vel^2
	real3 S = real3_real_mul(real3_cross(U->D, U->B), 1. / nSq);
	J_g.x -= .5 * S.x * (1. + nSq) / speedOfLightSq;
	J_g.y -= .5 * S.y * (1. + nSq) / speedOfLightSq;
	J_g.z -= .5 * S.z * (1. + nSq) / speedOfLightSq;
	
	deriv->D_g.x -= J_g.x;
	deriv->D_g.y -= J_g.y;
	deriv->D_g.z -= J_g.z;

	// source of phi_g is T_00 is rho + .5 (D^2 + B^2)
	real T_00_over_c2 = 0;
		//matter
	T_00_over_c2 += U->ion_rho + U->elec_rho;
		//electromagnetism			
	T_00_over_c2 += calc_EM_energy(solver, U, x) / speedOfLightSq;
	
	deriv->phi_g += T_00_over_c2 * solver->divPhiWavespeed_g / unit_m_per_s;


<? if not require 'hydro.coord.cartesian'.is(solver.coord) then ?>
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	prim_t W = primFromCons(solver, *U, x);
	real3 conn1_u = coord_conn_trace23(x);
	<? for _,fluid in ipairs(eqn.fluids) do ?>{
		real3 m_conn_vv = coord_conn_apply23(W.<?=fluid?>_v, U-><?=fluid?>_m, x);
		deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, m_conn_vv);	//-Conn^i_jk rho v^j v^k 
		deriv-><?=fluid?>_m = real3_add(deriv-><?=fluid?>_m, real3_real_mul(coord_raise(coord_conn_trace13(x), x), W.<?=fluid?>_P));		//+Conn^j_kj g^ki P
		deriv-><?=fluid?>_m = real3_add(deriv-><?=fluid?>_m, real3_real_mul(coord_conn_apply13(W.<?=fluid?>_v, U-><?=fluid?>_m, x), (solver->heatCapacityRatio - 1.) ));	//+ (gamma-1) rho v^k v^l Gamma_kjl g^ij	
		deriv-><?=fluid?>_ETotal -= (solver->heatCapacityRatio - 1.) * coord_conn_apply123(W.<?=fluid?>_v, W.<?=fluid?>_v, U-><?=fluid?>_m, x);		//- (gamma-1) rho v^j v^k v^l Gamma_jkl
	}<? end ?>
<? end ?>
}

<? elseif moduleName == "constrainU" then ?>

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	global cons_t* U = UBuf + index;
	real3 x = cellBuf[index].pos;
	prim_t W = primFromCons(solver, *U, x);

<? for _,fluid in ipairs(eqn.fluids) do
?>	W.<?=fluid?>_rho = max((real)W.<?=fluid?>_rho, (real)solver->min_<?=fluid?>_rho);
	W.<?=fluid?>_P = max((real)W.<?=fluid?>_P, (real)solver->min_<?=fluid?>_P);
<? end
?>

	*U = consFromPrim(solver, W, x);
}

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
