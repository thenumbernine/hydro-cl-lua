<?
local clnumber = require 'cl.obj.number'
local fluids = eqn.fluids
?>

#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>
#define sqrt_2 <?=('%.50f'):format(math.sqrt(2))?>

// r_e = q_e / m_e
// r_e = m q_e / m_i
// if q_e = 1 then r_e = m / m_i
// https://en.wikipedia.org/wiki/Mass-to-charge_ratio
// q_e / m_e = -1.758820024e+11 C/kg
#define elecChargeMassRatio			(solver->ionElectronMassRatio / (solver->ionMass / unit_kg))

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

<? if solver.hasCalcDTCode then ?>
//2014 Abgrall, Kumar eqn 2.25
// dt < sqrt(EInt_a/rho_a) sqrt(2) |lHat_r^a| / |E + v_a cross B|
//lHat_r^a = lHat_r for a=i, -lHat_r/m for a=e
kernel void calcDT(
	constant solver_t* solver,
	global real* dtBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	
	real3 x = cell_x(i);
	const global cons_t* U = UBuf + index;
	
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;

	prim_t W = primFromCons(solver, *U, x);
	real lHat_ion = normalizedIonLarmorRadius;
	real lHat_elec = lHat_ion / solver->ionElectronMassRatio;
<? for _,fluid in ipairs(fluids) do ?>
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
<? end ?>

<? for side=0,solver.dim-1 do ?>
cons_t fluxFromCons_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	prim_t W = primFromCons(solver, U, x);
	cons_t F;

<? 
for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_vj = W.<?=fluid?>_v.s<?=side?>;
	real <?=fluid?>_HTotal = U.<?=fluid?>_ETotal + W.<?=fluid?>_P;
	
	F.<?=fluid?>_rho = U.<?=fluid?>_m.s<?=side?>;
	F.<?=fluid?>_m = real3_real_mul(U.<?=fluid?>_m, <?=fluid?>_vj);
<? 	for i=0,2 do
?>	F.<?=fluid?>_m.s<?=i?> += coord_gU<?=i?><?=side?>(x) * W.<?=fluid?>_P;
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
	<? if side == 0 then ?>
	F.D = _real3(U.phi * solver->divPhiWavespeed / unit_m_per_s, H.z, -H.y);
	F.B = _real3(U.psi * solver->divPsiWavespeed / unit_m_per_s, -E.z, E.y);
	<? elseif side == 1 then ?>
	F.D = _real3(-H.z, U.phi * solver->divPhiWavespeed / unit_m_per_s, H.x);
	F.B = _real3(E.z, U.psi * solver->divPsiWavespeed / unit_m_per_s, -E.x);
	<? elseif side == 2 then ?>
	F.D = _real3(H.y, -H.x, U.phi * solver->divPhiWavespeed / unit_m_per_s);
	F.B = _real3(-E.y, E.x, U.psi * solver->divPsiWavespeed / unit_m_per_s);
	<? end ?>
	F.phi = U.D.s<?=side?> * solver->divPhiWavespeed / unit_m_per_s;
	F.psi = U.B.s<?=side?> * solver->divPsiWavespeed / unit_m_per_s;

	return F;
}
<? end ?>

eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3 n
) {
	prim_t WL = primFromCons(solver, UL, x);
	prim_t WR = primFromCons(solver, UR, x);
	eigen_t eig;

<? for _,fluid in ipairs(fluids) do ?>

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

<? for side=0,solver.dim-1 do ?>
eigen_t eigen_forCell_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	prim_t W = primFromCons(solver, U, x);
<? for _,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_vSq = coordLenSq(W.<?=fluid?>_v, x);
	real <?=fluid?>_eKin = .5 * <?=fluid?>_vSq;
	real <?=fluid?>_hTotal = calc_hTotal(solver, W.<?=fluid?>_rho, W.<?=fluid?>_P, U.<?=fluid?>_ETotal);
	real <?=fluid?>_CsSq = (solver->heatCapacityRatio - 1.) * (<?=fluid?>_hTotal - <?=fluid?>_eKin);
	real <?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);
<? end ?>	
	return (eigen_t){
<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_v = W.<?=fluid?>_v,
		.<?=fluid?>_hTotal = <?=fluid?>_hTotal,
		.<?=fluid?>_vSq = <?=fluid?>_vSq,
		.<?=fluid?>_Cs = <?=fluid?>_Cs,
<? end ?>	
	};
}
<? end ?>

<? for side=0,solver.dim-1 do 
	local prefix
	if side == 0 then
		prefix = [[
	const real nx = 1, ny = 0, nz = 0;
	const real n1x = 0, n1y = 1, n1z = 0;
	const real n2x = 0, n2y = 0, n2z = 1;
	real ion_v_n = ion_v.x, ion_v_n1 = ion_v.y, ion_v_n2 = ion_v.z;
	real elec_v_n = elec_v.x, elec_v_n1 = elec_v.y, elec_v_n2 = elec_v.z;
]] 
	elseif side == 1 then
		prefix = [[
	const real nx = 0, ny = 1, nz = 0;
	const real n1x = 0, n1y = 0, n1z = 1;
	const real n2x = 1, n2y = 0, n2z = 0;
	real ion_v_n = ion_v.y, ion_v_n1 = ion_v.z, ion_v_n2 = ion_v.x;
	real elec_v_n = elec_v.y, elec_v_n1 = elec_v.z, elec_v_n2 = elec_v.x;
]] 
	elseif side == 2 then
		prefix = [[
	const real nx = 0, ny = 0, nz = 1;
	const real n1x = 1, n1y = 0, n1z = 0;
	const real n2x = 0, n2y = 1, n2z = 0;
	real ion_v_n = ion_v.z, ion_v_n1 = ion_v.x, ion_v_n2 = ion_v.y;
	real elec_v_n = elec_v.z, elec_v_n1 = elec_v.x, elec_v_n2 = elec_v.y;
]]
	end
	prefix = [[
	sym3 gU = coord_gU(x);
	real gUjj = gU.s]]..side..side..[[;
	real sqrt_gUjj = coord_sqrt_gU]]..side..side..[[(x);
	
	//g^ij for fixed j=side
	real3 ion_v = eig.ion_v;
	real3 ion_vL = coord_lower(ion_v, x);
	real ion_hTotal = eig.ion_hTotal;
	real ion_vSq = real3_dot(ion_v, ion_vL);
	real ion_Cs = eig.ion_Cs;
	real ion_Cs_over_sqrt_gUjj = ion_Cs / sqrt_gUjj; 
	
	real3 elec_v = eig.elec_v;
	real3 elec_vL = coord_lower(elec_v, x);
	real elec_hTotal = eig.elec_hTotal;
	real elec_vSq = real3_dot(elec_v, elec_vL);
	real elec_Cs = eig.elec_Cs;
	real elec_Cs_over_sqrt_gUjj = elec_Cs / sqrt_gUjj; 

]] .. prefix

	local gUdef = '\treal3 gUj = _real3(\n'
	for i=0,2 do
		gUdef = gUdef .. '\t\tcoord_gU'..side..i..'(x)'..(i<2 and ',' or '')..'\n'
	end
	gUdef = gUdef .. '\t);\n'
	prefix = gUdef .. prefix
?>

waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t UX,
	real3 x
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

<? 
				if side == 0 then 
?>
	real sqrt_gUxx = sqrt_gUjj;
<?
					for i,fluid in ipairs(fluids) do
?>
	Y[<?=5*i-5?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.x / sqrt_gUxx)
		+ X[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x - <?=fluid?>_Cs / sqrt_gUxx)
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
	Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x * gU.xy / gU.xx - <?=fluid?>_v.y)
		+ X[<?=5*i-4?>] * -gU.xy / gU.xx
		+ X[<?=5*i-3?>];
	Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x * gU.xz / gU.xx - <?=fluid?>_v.z)
		+ X[<?=5*i-4?>] * -gU.xz / gU.xx
		+ X[<?=5*i-2?>];
	Y[<?=5*i-1?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.x / sqrt_gUxx)
		+ X[<?=5*i-4?>] * (-heatRatioMinusOne * <?=fluid?>_vL.x + <?=fluid?>_Cs / sqrt_gUxx)
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

<?
				elseif side == 1 then 
?>
	real sqrt_gUyy = sqrt_gUjj;
<?	
					for i,fluid in ipairs(fluids) do
?>
	Y[<?=5*i-5?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.y / sqrt_gUyy)
		+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y - <?=fluid?>_Cs / sqrt_gUyy)
		+ X[<?=5*i-2?>] * -heatRatioMinusOne * <?=fluid?>_vL.z
		+ X[<?=5*i-1?>] * heatRatioMinusOne
	) * <?=fluid?>_invDenom;
	Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y * gU.xy / gU.yy - <?=fluid?>_v.x)
		+ X[<?=5*i-4?>]
		+ X[<?=5*i-3?>] * -gU.xy / gU.yy;
	Y[<?=5*i-3?>] = (X[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)
		+ X[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y
		+ X[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z
		+ X[<?=5*i-1?>] * -2. * heatRatioMinusOne
	) * <?=fluid?>_invDenom;
	Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y * gU.yz / gU.yy - <?=fluid?>_v.z)
		+ X[<?=5*i-3?>] * -gU.yz / gU.yy
		+ X[<?=5*i-2?>];
	Y[<?=5*i-1?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.y / sqrt_gUyy)
		+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * (-heatRatioMinusOne * <?=fluid?>_vL.y + <?=fluid?>_Cs / sqrt_gUyy)
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

<? 
				elseif side == 2 then
?>
	real sqrt_gUzz = sqrt_gUjj;
<?
					for i,fluid in ipairs(fluids) do
?>
	Y[<?=5*i-5?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq + <?=fluid?>_Cs * <?=fluid?>_v.z / sqrt_gUzz)
		+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y
		+ X[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z - <?=fluid?>_Cs / sqrt_gUzz)
		+ X[<?=5*i-1?>] * heatRatioMinusOne
	) * <?=fluid?>_invDenom;
	Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z * gU.xz / gU.zz - <?=fluid?>_v.x)
		+ X[<?=5*i-4?>]
		+ X[<?=5*i-2?>] * -gU.xz / gU.zz;
	Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z * gU.yz / gU.zz - <?=fluid?>_v.y)
		+ X[<?=5*i-3?>]
		+ X[<?=5*i-2?>] * -gU.yz / gU.zz;
	Y[<?=5*i-2?>] = (X[<?=5*i-5?>] * (<?=fluid?>_denom - heatRatioMinusOne * <?=fluid?>_vSq)
		+ X[<?=5*i-4?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.y
		+ X[<?=5*i-2?>] * 2. * heatRatioMinusOne * <?=fluid?>_vL.z
		+ X[<?=5*i-1?>] * -2. * heatRatioMinusOne
	) * <?=fluid?>_invDenom;
	Y[<?=5*i-1?>] = (X[<?=5*i-5?>] * (.5 * heatRatioMinusOne * <?=fluid?>_vSq - <?=fluid?>_Cs * <?=fluid?>_v.z / sqrt_gUzz)
		+ X[<?=5*i-4?>] * -heatRatioMinusOne * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * -heatRatioMinusOne * <?=fluid?>_vL.y
		+ X[<?=5*i-2?>] * (-heatRatioMinusOne * <?=fluid?>_vL.z + <?=fluid?>_Cs / sqrt_gUzz)
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

<?
				end 
?>
	return UY;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t UX,	//numWaves = 16
	real3 x
) {
	cons_t UY;
	real* Y = UY.ptr;
	real* X = UX.ptr;
	
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	real sqrt_eps = sqrt(eps);	// TODO sqrt units
	real sqrt_mu = sqrt(mu);

	<?=prefix?>
<? 
				if side == 0 then
?>
	real sqrt_gUxx = sqrt_gUjj;
<?
					for i,fluid in ipairs(fluids) do
?>
	Y[<?=5*i-5?>] = X[<?=5*i-5?>] + X[<?=5*i-4?>] + X[<?=5*i-1?>];
	Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * sqrt_gUxx)
		+ X[<?=5*i-4?>] * <?=fluid?>_v.x
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * sqrt_gUxx);
	Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * gU.xy / sqrt_gUxx)
		+ X[<?=5*i-4?>] * <?=fluid?>_v.y
		+ X[<?=5*i-3?>]
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * gU.xy / sqrt_gUxx);
	Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * gU.xz / sqrt_gUxx)
		+ X[<?=5*i-4?>] * <?=fluid?>_v.z
		+ X[<?=5*i-2?>]
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * gU.xz / sqrt_gUxx);
	Y[<?=5*i-1?>] = X[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.x / sqrt_gUxx)
		+ X[<?=5*i-4?>] * <?=fluid?>_vSq / 2.
		+ X[<?=5*i-3?>] * <?=fluid?>_vL.y
		+ X[<?=5*i-2?>] * <?=fluid?>_vL.z
		+ X[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.x / sqrt_gUxx);
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
<? 
				elseif side == 1 then
?>
	real sqrt_gUyy = sqrt_gUjj;
<?
					for i,fluid in ipairs(fluids) do
?>	
	Y[<?=5*i-5?>] = X[<?=5*i-5?>] + X[<?=5*i-3?>] + X[<?=5*i-1?>];
	Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * gU.xy / sqrt_gUyy)
		+ X[<?=5*i-4?>]
		+ X[<?=5*i-3?>] * <?=fluid?>_v.x
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * gU.xy / sqrt_gUyy);
	Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * sqrt_gUyy)
		+ X[<?=5*i-3?>] * <?=fluid?>_v.y
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * sqrt_gUyy);
	Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * gU.yz / sqrt_gUyy)
		+ X[<?=5*i-3?>] * <?=fluid?>_v.z
		+ X[<?=5*i-2?>]
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * gU.yz / sqrt_gUyy);
	Y[<?=5*i-1?>] = X[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.y / sqrt_gUyy)
		+ X[<?=5*i-4?>] * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * <?=fluid?>_vSq / 2.
		+ X[<?=5*i-2?>] * <?=fluid?>_vL.z
		+ X[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.y / sqrt_gUyy);
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

<?
				elseif side == 2 then
?>
	real sqrt_gUzz = sqrt_gUjj;
<?
					for i,fluid in ipairs(fluids) do
?>
	Y[<?=5*i-5?>] = X[<?=5*i-5?>] + X[<?=5*i-2?>] + X[<?=5*i-1?>];
	Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.x - <?=fluid?>_Cs * gU.xz / sqrt_gUzz)
		+ X[<?=5*i-4?>]
		+ X[<?=5*i-2?>] * <?=fluid?>_v.x
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.x + <?=fluid?>_Cs * gU.xz / sqrt_gUzz);
	Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.y - <?=fluid?>_Cs * gU.yz / sqrt_gUzz)
		+ X[<?=5*i-3?>]
		+ X[<?=5*i-2?>] * <?=fluid?>_v.y
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.y + <?=fluid?>_Cs * gU.yz / sqrt_gUzz);
	Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (<?=fluid?>_v.z - <?=fluid?>_Cs * sqrt_gUzz)
		+ X[<?=5*i-2?>] * <?=fluid?>_v.z
		+ X[<?=5*i-1?>] * (<?=fluid?>_v.z + <?=fluid?>_Cs * sqrt_gUzz);
	Y[<?=5*i-1?>] = X[<?=5*i-5?>] * (<?=fluid?>_hTotal - <?=fluid?>_Cs * <?=fluid?>_v.z / sqrt_gUzz)
		+ X[<?=5*i-4?>] * <?=fluid?>_vL.x
		+ X[<?=5*i-3?>] * <?=fluid?>_vL.y
		+ X[<?=5*i-2?>] * <?=fluid?>_vSq / 2.
		+ X[<?=5*i-1?>] * (<?=fluid?>_hTotal + <?=fluid?>_Cs * <?=fluid?>_v.z / sqrt_gUzz);
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

<?
				end
?>

	return UY;
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t UX,
	real3 x
) {
	<?=prefix?>
	cons_t UY;
	real* X = UX.ptr;
<?
					for i,fluid	in ipairs(fluids) do 
?>
	UY.<?=fluid?>_rho = X[<?=5*i-4?>] * nx 
		+ X[<?=5*i-3?>] * ny 
		+ X[<?=5*i-2?>] * nz;
	UY.<?=fluid?>_m.x = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.x + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * gUj.x)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.x * nx - (solver->heatCapacityRatio - 1.) * gUj.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.x * ny - (solver->heatCapacityRatio - 1.) * gUj.x * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.x * nz - (solver->heatCapacityRatio - 1.) * gUj.x * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * nx;
	UY.<?=fluid?>_m.y = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.y + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * gUj.y)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.y * nx - (solver->heatCapacityRatio - 1.) * gUj.y * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.y * ny - (solver->heatCapacityRatio - 1.) * gUj.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.y * nz - (solver->heatCapacityRatio - 1.) * gUj.y * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (solver->heatCapacityRatio - 1.) * ny;
	UY.<?=fluid?>_m.z = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.z + (solver->heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * gUj.z)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.z * nx - (solver->heatCapacityRatio - 1.) * gUj.z * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.z * ny - (solver->heatCapacityRatio - 1.) * gUj.z * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.z * nz - (solver->heatCapacityRatio - 1.) * gUj.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)
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
	<? if side == 0 then ?>
	UY.D = _real3(solver->divPhiWavespeed / unit_m_per_s * UX.phi, H.z, -H.y);
	UY.B = _real3(solver->divPsiWavespeed / unit_m_per_s * UX.psi, -E.z, E.y);
	<? elseif side == 1 then ?>
	UY.D = _real3(-H.z, solver->divPhiWavespeed / unit_m_per_s * UX.phi, H.x);
	UY.B = _real3(E.z, solver->divPsiWavespeed / unit_m_per_s * UX.psi, -E.x);
	<? elseif side == 2 then ?>
	UY.D = _real3(H.y, -H.x, solver->divPhiWavespeed / unit_m_per_s * UX.phi);
	UY.B = _real3(-E.y, E.x, solver->divPsiWavespeed / unit_m_per_s * UX.psi);
	<? end ?>
	UY.phi = solver->divPhiWavespeed / unit_m_per_s * UX.D.s<?=side?>;
	UY.psi = solver->divPsiWavespeed / unit_m_per_s * UX.B.s<?=side?>;
	UY.ePot = 0;

	return UY;
}
<? end ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
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
	deriv->D.x -= U->ion_m.x * solver->ionChargeMassRatio / unit_C_per_kg + U->elec_m.x * elecChargeMassRatio / unit_C_per_kg;
	deriv->D.y -= U->ion_m.y * solver->ionChargeMassRatio / unit_C_per_kg + U->elec_m.y * elecChargeMassRatio / unit_C_per_kg;
	deriv->D.z -= U->ion_m.z * solver->ionChargeMassRatio / unit_C_per_kg + U->elec_m.z * elecChargeMassRatio / unit_C_per_kg;
	
	deriv->phi += eps * (U->ion_rho * solver->ionChargeMassRatio / unit_C_per_kg + U->elec_rho * elecChargeMassRatio / unit_C_per_kg) * solver->divPhiWavespeed / unit_m_per_s;

<? if not require 'coord.cartesian'.is(solver.coord) then ?>
	real3 x = cell_x(i);
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	prim_t W = primFromCons(solver, *U, x);
	real3 conn1_u = coord_conn_trace23(x);
	<? for _,fluid in ipairs(fluids) do ?>{
		real3 m_conn_vv = coord_conn_apply23(W.<?=fluid?>_v, U-><?=fluid?>_m, x);
		deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, m_conn_vv);	//-Conn^i_jk rho v^j v^k 
		deriv-><?=fluid?>_m = real3_add(deriv-><?=fluid?>_m, real3_real_mul(coord_raise(coord_conn_trace13(x), x), W.<?=fluid?>_P));		//+Conn^j_kj g^ki P
	}<? end ?>
<? end ?>
}

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	global cons_t* U = UBuf + index;
	real3 x = cell_x(i);
	prim_t W = primFromCons(solver, *U, x);

<? for _,fluid in ipairs(fluids) do
?>	W.<?=fluid?>_rho = max((real)W.<?=fluid?>_rho, (real)solver->min_<?=fluid?>_rho);
	W.<?=fluid?>_P = max((real)W.<?=fluid?>_P, (real)solver->min_<?=fluid?>_P);
<? end
?>
	
	*U = consFromPrim(solver, W, x);
}
