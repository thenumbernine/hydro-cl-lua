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

<? elseif moduleName == "primFromCons" then 
depmod{
	"real3",
	"solver_t",
	"prim_t",
	"cons_t",
	"eqn.common",	-- calc_*
}
?>

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
		.ePot = U.ePot,
	};
}

<? elseif moduleName == "consFromPrim" then 
depmod{
	"real3",
	"solver_t",
	"prim_t",
	"cons_t",
	"eqn.common",	-- calc_*
}
?>

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
		.ePot = W.ePot,
	};
}

<? elseif moduleName == "eqn.dU-dW" then 
depmod{
	"coord_lower",
}
?>

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
		.ePot = W.ePot,
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
		.ePot = U.ePot,
	};
}

<? elseif moduleName == "eqn.common" then 
depmod{
	"units",
}
?>

real3 calc_EField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	return real3_real_mul(U.D, 1. / eps);
}
 
real3 calc_HField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) { 
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return real3_real_mul(U.B, 1. / mu);
}

real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
real calc_HTotal(real P, real ETotal) { return P + ETotal; }
real calc_hTotal(constant <?=solver.solver_t?>* solver, real rho, real P, real ETotal) { return calc_HTotal(P, ETotal) / rho; }

real calc_rho_from_U(<?=eqn.cons_t?> U) {
	real rho = 0.;
<? for _,fluid in ipairs(eqn.fluids) do 
?>	rho += U.<?=fluid?>_rho;
<? end 
?>	return rho;
}

real calc_rho_from_W(<?=eqn.prim_t?> W) {
	real rho = 0.;
<? for _,fluid in ipairs(eqn.fluids) do 
?>	rho += W.<?=fluid?>_rho;
<? end 
?>	return rho;
}

real calc_EPot(<?=eqn.cons_t?> U) { return calc_rho_from_U(U) * U.ePot; }
real calc_EPot_from_W(<?=eqn.prim_t?> W) { return calc_rho_from_W(W) * W.ePot; }

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

<? elseif moduleName == "applyInitCond" then 
depmod{
	"consFromPrim",
	"cartesianToCoord",
}
?>

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

	<?=eqn.prim_t?> W = {
<? 
if eqn.useEulerInitState then 
?>
		.ion_rho = rho,
		.elec_rho = rho / solver->ionElectronMassRatio, 

		// "the electron pressure is taken to be elec_P = 5 ion_rho"
		// is that arbitrary?
		.elec_P = 5. * rho,
		
		// "the ion pressure is 1/100th the electron pressure"
		// is that from the mass ratio of ion/electron?
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
	UBuf[index] = consFromPrim(solver, W, x);
}

<? elseif moduleName == "fluxFromCons" then 
depmod{
	"units",
	"primFromCons",
}
?>

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
?>	F.ePot = 0.;
	
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;

	//taken from glm-maxwell instead of the 2014 Abgrall, Kumar
	real3 E = real3_real_mul(U.D, 1. / eps);
	real3 H = real3_real_mul(U.B, 1. / mu);
	if (n.side == 0) {
		F.D = _real3(U.phi * solver->divPhiWavespeed / unit_m_per_s, H.z, -H.y);
		F.B = _real3(U.psi * solver->divPsiWavespeed / unit_m_per_s, -E.z, E.y);
	} else if (n.side == 1) {
		F.D = _real3(-H.z, U.phi * solver->divPhiWavespeed / unit_m_per_s, H.x);
		F.B = _real3(E.z, U.psi * solver->divPsiWavespeed / unit_m_per_s, -E.x);
	} else if (n.side == 2) {
		F.D = _real3(H.y, -H.x, U.phi * solver->divPhiWavespeed / unit_m_per_s);
		F.B = _real3(-E.y, E.x, U.psi * solver->divPsiWavespeed / unit_m_per_s);
	}
	F.phi = normal_vecDotN1(n, U.D) * solver->divPhiWavespeed / unit_m_per_s;
	F.psi = normal_vecDotN1(n, U.B) * solver->divPsiWavespeed / unit_m_per_s;

	return F;
}

<? elseif moduleName == "eigen_forInterface" then 
depmod{
	"eigen_t",
	"primFromCons",
}
?>

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

<? elseif moduleName == "eigen_forCell" then 
depmod{
	"eigen_t",
	"primFromCons",
}
?>

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

<? elseif moduleName == "eigen_left/rightTransform" then 
depmod{
	"units",
	"eigen_t",
	"waves_t",
	"coord_lower",
}
?>

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
		//EM
		X += 10;
		Y += 10;
		Y[0] = ((-(sqrt_eps * (X[0] - X[6]))) / sqrt_2);
		Y[1] = ((-(sqrt_eps * (X[3] - X[7]))) / sqrt_2);
		Y[2] = (((X[2] * sqrt_mu) + (X[4] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));
		Y[3] = (((X[1] * sqrt_mu) - (X[5] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
		Y[4] = ((-((X[2] * sqrt_mu) - (X[4] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));
		Y[5] = (((X[1] * sqrt_mu) + (X[5] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));
		Y[6] = ((sqrt_eps * (X[0] + X[6])) / sqrt_2);
		Y[7] = ((sqrt_eps * (X[3] + X[7])) / sqrt_2);

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
		//EM
		X += 10;
		Y += 10;
		Y[0] = ((-(sqrt_eps * (X[1] - X[6]))) / sqrt_2);
		Y[1] = ((-(sqrt_eps * (X[4] - X[7]))) / sqrt_2);
		Y[2] = (((X[2] * sqrt_mu) - (X[3] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
		Y[3] = (((X[0] * sqrt_mu) + (X[5] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));
		Y[4] = (((X[2] * sqrt_mu) + (X[3] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));
		Y[5] = ((-((X[0] * sqrt_mu) - (X[5] * sqrt_eps))) / (sqrt_mu * sqrt_eps * sqrt_2));
		Y[6] = ((sqrt_eps * (X[1] + X[6])) / sqrt_2);
		Y[7] = ((sqrt_eps * (X[4] + X[7])) / sqrt_2);

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
		//EM
		X += 10;
		Y += 10;
		Y[0] = ((-(sqrt_eps * (X[2] - X[6]))) / sqrt_2);
		Y[1] = ((-(sqrt_eps * (X[5] - X[7]))) / sqrt_2);
		Y[2] = (((X[1] * sqrt_mu) + (X[3] * sqrt_eps)) / (sqrt_mu * sqrt_2 * sqrt_eps));
		Y[3] = (((X[0] * sqrt_mu) - (X[4] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
		Y[4] = (((X[1] * sqrt_mu) - (X[3] * sqrt_eps)) / (-(sqrt_mu * sqrt_2 * sqrt_eps)));
		Y[5] = (((X[0] * sqrt_mu) + (X[4] * sqrt_eps)) / (sqrt_mu * sqrt_eps * sqrt_2));
		Y[6] = ((sqrt_eps * (X[2] + X[6])) / sqrt_2);
		Y[7] = ((sqrt_eps * (X[5] + X[7])) / sqrt_2);

	}

	return UY;
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t UX,	//numWaves = 16
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
		//EM
		X += 10;
		Y += 10;
		Y[0] = ((-(X[0] - X[6])) / (sqrt_2 * sqrt_eps));
		Y[1] = ((-(sqrt_eps * (X[3] - X[5]))) / sqrt_2);
		Y[2] = ((sqrt_eps * (X[2] - X[4])) / sqrt_2);
		Y[3] = ((-(X[1] - X[7])) / (sqrt_2 * sqrt_eps));
		Y[4] = ((sqrt_mu * (X[2] + X[4])) / sqrt_2);
		Y[5] = ((sqrt_mu * (X[3] + X[5])) / sqrt_2);
		Y[6] = ((X[0] + X[6]) / (sqrt_2 * sqrt_eps));
		Y[7] = ((X[1] + X[7]) / (sqrt_2 * sqrt_eps));
	
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
		//EM
		X += 10;
		Y += 10;
		Y[0] = ((sqrt_eps * (X[3] - X[5])) / sqrt_2);
		Y[1] = ((-(X[0] - X[6])) / (sqrt_2 * sqrt_eps));
		Y[2] = ((-(sqrt_eps * (X[2] - X[4]))) / sqrt_2);
		Y[3] = ((sqrt_mu * (X[2] + X[4])) / sqrt_2);
		Y[4] = ((-(X[1] - X[7])) / (sqrt_2 * sqrt_eps));
		Y[5] = ((sqrt_mu * (X[3] + X[5])) / sqrt_2);
		Y[6] = ((X[0] + X[6]) / (sqrt_2 * sqrt_eps));
		Y[7] = ((X[1] + X[7]) / (sqrt_2 * sqrt_eps));

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
		//EM
		X += 10;
		Y += 10;
		Y[0] = ((-(sqrt_eps * (X[3] - X[5]))) / sqrt_2);
		Y[1] = ((sqrt_eps * (X[2] - X[4])) / sqrt_2);
		Y[2] = ((-(X[0] - X[6])) / (sqrt_2 * sqrt_eps));
		Y[3] = ((sqrt_mu * (X[2] + X[4])) / sqrt_2);
		Y[4] = ((sqrt_mu * (X[3] + X[5])) / sqrt_2);
		Y[5] = ((-(X[1] - X[7])) / (sqrt_2 * sqrt_eps));
		Y[6] = ((X[0] + X[6]) / (sqrt_2 * sqrt_eps));
		Y[7] = ((X[1] + X[7]) / (sqrt_2 * sqrt_eps));

	}

	return UY;
}

<? elseif moduleName == "eigen_fluxTransform" then 
depmod{
	"units",
}
?>

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
	UY.ePot = 0;

	return UY;
}

<? elseif moduleName == "addSource" then 
depmod{
	"units",
	"primFromCons",
}
?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;

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

<? elseif moduleName == "constrainU" then 
depmod{
	"primFromCons",
	"consFromPrim",
}
?>

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

<? elseif moduleName == "calcDT" then 
depmod{
	"units",
	"primFromCons",
}
?>

//2014 Abgrall, Kumar eqn 2.25
// dt < sqrt(EInt_a/rho_a) sqrt(2) |lHat_r^a| / |E + v_a cross B|
//lHat_r^a = lHat_r for a=i, -lHat_r/m for a=e
kernel void calcDT(
	constant solver_t* solver,
	global real* dtBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	
	real3 x = cellBuf[index].pos;
	const global cons_t* U = UBuf + index;
	
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;

	prim_t W = primFromCons(solver, *U, x);
	real lHat_ion = normalizedIonLarmorRadius;
	real lHat_elec = lHat_ion / solver->ionElectronMassRatio;
<? for _,fluid in ipairs(eqn.fluids) do ?>
	real EInt_<?=fluid?> = calc_<?=fluid?>_EInt(solver, W);
	real LorentzForceSq_<?=fluid?> = coordLenSq(
		real3_add(
			real3_real_mul(W.D, 1. / eps),
			real3_cross(W.<?=fluid?>_v, W.B)
		), x);
	real sqrt_EInt_lHat_over_rho_<?=fluid?> = sqrt(2. * EInt_<?=fluid?> * lHat_<?=fluid?> / (W.<?=fluid?>_rho * LorentzForceSq_<?=fluid?>));
<? end ?>

	dtBuf[index] = min(
		sqrt_EInt_lHat_over_rho_ion,
		sqrt_EInt_lHat_over_rho_elec);
}

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
