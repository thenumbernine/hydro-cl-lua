<?
local clnumber = require 'cl.obj.number'
local fluids = eqn.fluids
?>

#define ionReferenceThermalVelocity (solver->ionLarmorRadius * solver->ionChargeMassRatio * solver->referenceMagneticField)
#define normalizedSpeedOfLight 		(solver->speedOfLight / ionReferenceThermalVelocity)
#define normalizedSpeedOfLightSq 	(normalizedSpeedOfLight * normalizedSpeedOfLight)
#define normalizedIonLarmorRadius 	(solver->ionLarmorRadius / solver->referenceLength)
#define normalizedIonDebyeLength	(solver->ionDebyeLength / solver->ionLarmorRadius)

// r_e = q_e / m_e
// r_e = m q_e / m_i
// if q_e = 1 then r_e = m / m_i
// https://en.wikipedia.org/wiki/Mass-to-charge_ratio
// q_e / m_e = -1.758820024e+11 C/kg
#define elecChargeMassRatio			(solver->ionElectronMassRatio / solver->ionMass)

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

<? local useGEM = true ?>

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
?>

	//taken from glm-maxwell instead of the 2014 Abgrall, Kumar
	real3 B = U.B;
	real3 E = U.E;
	<? if side == 0 then ?>
	F.E = _real3(normalizedSpeedOfLight * solver->divPhiWavespeed * U.phi, normalizedSpeedOfLight * B.z, -normalizedSpeedOfLight * B.y);
	F.B = _real3(solver->divPsiWavespeed * U.psi, -E.z, E.y);
	<? elseif side == 1 then ?>
	F.E = _real3(-normalizedSpeedOfLight * B.z, normalizedSpeedOfLight * solver->divPhiWavespeed * U.phi, normalizedSpeedOfLight * B.x);
	F.B = _real3(E.z, solver->divPsiWavespeed * U.psi, -E.x);
	<? elseif side == 2 then ?>
	F.E = _real3(normalizedSpeedOfLight * B.y, -normalizedSpeedOfLight * B.x, normalizedSpeedOfLight * solver->divPhiWavespeed * U.phi);
	F.B = _real3(-E.y, E.x, solver->divPsiWavespeed * U.psi);
	<? end ?>
	F.phi = E.s<?=side?> * solver->divPhiWavespeed;
	F.psi = B.s<?=side?> * solver->divPsiWavespeed * normalizedSpeedOfLight;

<? if useGEM then ?>
	//same deal but with E_g and B_g
	real3 B_g = U.B_g;
	real3 E_g = U.E_g;
	<? if side == 0 then ?>
	F.E_g = _real3(normalizedSpeedOfLight * solver->divPhiGWavespeed * U.phi_g, normalizedSpeedOfLight * B_g.z, -normalizedSpeedOfLight * B_g.y);
	F.B_g = _real3(solver->divPsiGWavespeed * U.psi_g, -E_g.z, E_g.y);
	<? elseif side == 1 then ?>
	F.E_g = _real3(-normalizedSpeedOfLight * B_g.z, normalizedSpeedOfLight * solver->divPhiGWavespeed * U.phi_g, normalizedSpeedOfLight * B_g.x);
	F.B_g = _real3(E_g.z, solver->divPsiGWavespeed * U.psi_g, -E_g.x);
	<? elseif side == 2 then ?>
	F.E_g = _real3(normalizedSpeedOfLight * B_g.y, -normalizedSpeedOfLight * B_g.x, normalizedSpeedOfLight * solver->divPhiGWavespeed * U_g.phi);
	F.B_g = _real3(-E_g.y, E_g.x, solver->divPsiGWavespeed * U.psi_g);
	<? end ?>
	F.phi_g = E_g.s<?=side?> * solver->divPhiGWavespeed;
	F.psi_g = B_g.s<?=side?> * solver->divPsiGWavespeed * normalizedSpeedOfLight;
<? else ?>
	F.E_g = real3_zero;
	F.B_g = real3_zero;
	F.phi_g = 0;
	F.psi_g = 0;
<? end ?>
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
	X += 2;
	Y += 2;
	//EM & gravity
	<? for i=1,2 do ?>
	X += 8;
	Y += 8;
	Y[0] = (X[0] - X[6] * normalizedSpeedOfLight) * .5;
	Y[1] = (X[3] - X[7] / normalizedSpeedOfLight) * .5;
	Y[2] = (X[1] - X[5] * normalizedSpeedOfLight) * .5;
	Y[3] = (X[4] * normalizedSpeedOfLight + X[2]) * .5;
	Y[4] = (X[5] * normalizedSpeedOfLight + X[1]) * .5;
	Y[5] = (X[2] - X[4] * normalizedSpeedOfLight) * .5;
	Y[6] = (X[7] / normalizedSpeedOfLight + X[3]) * .5;
	Y[7] = (X[6] * normalizedSpeedOfLight + X[0]) * .5;
	<? end ?>
<? if not useGEM then ?>
	for (int k = 0; k < 8; ++k) Y[k] = X[k];
<? end ?>

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
	X += 2;
	Y += 2;
	//EM & gravity
	<? for i=1,2 do ?>
	X += 8;
	Y += 8;
	Y[0] = (X[1] - X[6] * normalizedSpeedOfLight) * .5;
	Y[1] = (X[4] - X[7] / normalizedSpeedOfLight) * .5;
	Y[2] = (X[5] * normalizedSpeedOfLight + X[0]) * .5;
	Y[3] = (X[2] - X[3] * normalizedSpeedOfLight) * .5;
	Y[4] = (X[0] - X[5] * normalizedSpeedOfLight) * .5;
	Y[5] = (X[3] * normalizedSpeedOfLight + X[2]) * .5;
	Y[6] = (X[7] / normalizedSpeedOfLight + X[4]) * .5;
	Y[7] = (X[6] * normalizedSpeedOfLight + X[1]) * .5;
	<? end ?>
<? if not useGEM then ?>
	for (int k = 0; k < 8; ++k) Y[k] = X[k];
<? end ?>

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
	X += 2;
	Y += 2;
	//EM & gravity
	<? for i=1,2 do ?>
	X += 8;
	Y += 8;
	Y[0] = (X[2] - X[6] * normalizedSpeedOfLight) * .5;
	Y[1] = (X[5] - X[7] / normalizedSpeedOfLight) * .5;
	Y[2] = (X[0] - X[4] * normalizedSpeedOfLight) * .5;
	Y[3] = (X[3] * normalizedSpeedOfLight + X[1]) * .5;
	Y[4] = (X[4] * normalizedSpeedOfLight + X[0]) * .5;
	Y[5] = (X[1] - X[3] * normalizedSpeedOfLight) * .5;
	Y[6] = (X[7] / normalizedSpeedOfLight + X[5]) * .5;
	Y[7] = (X[6] * normalizedSpeedOfLight + X[2]) * .5;
	<? end ?>
<? if not useGEM then ?>
	for (int k = 0; k < 8; ++k) Y[k] = X[k];
<? end ?>

<?
				end 
?>
	return UY;
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t UX,	//numWaves = 26
	real3 x
) {
	cons_t UY;
	real* Y = UY.ptr;
	real* X = UX.ptr;

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
	X += 2;
	Y += 2;
	//EM & gravity
	<? for i=1,2 do ?>
	X += 8;
	Y += 8;
	Y[0] = X[7] + X[0];
	Y[1] = X[4] + X[2];
	Y[2] = X[5] + X[3];
	Y[3] = X[6] + X[1];
	Y[4] = (X[3] - X[5]) / normalizedSpeedOfLight;
	Y[5] = (X[4] - X[2]) / normalizedSpeedOfLight;
	Y[6] = (X[7] - X[0]) / normalizedSpeedOfLight;
	Y[7] = (X[6] - X[1]) * normalizedSpeedOfLight;
	<? end ?>
<? if not useGEM then ?>
	for (int k = 0; k < 8; ++k) Y[k] = X[k];
<? end ?>

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
	X += 2;
	Y += 2;
	//EM & gravity
	<? for i=1,2 do ?>
	X += 8;
	Y += 8;
	Y[0] = X[4] + X[2];
	Y[1] = X[7] + X[0];
	Y[2] = X[5] + X[3];
	Y[3] = (X[5] - X[3]) / normalizedSpeedOfLight;
	Y[4] = X[6] + X[1];
	Y[5] = (X[2] - X[4]) / normalizedSpeedOfLight;
	Y[6] = (X[7] - X[0]) / normalizedSpeedOfLight;
	Y[7] = (X[6] - X[1]) * normalizedSpeedOfLight;
	<? end ?>
<? if not useGEM then ?>
	for (int k = 0; k < 8; ++k) Y[k] = X[k];
<? end ?>

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
	X += 2;
	Y += 2;
	//EM & gravity
	<? for i=1,2 do ?>
	X += 8;
	Y += 8;
	Y[0] = X[4] + X[2];
	Y[1] = X[5] + X[3];
	Y[2] = X[7] + X[0];
	Y[3] = (X[3] - X[5]) / normalizedSpeedOfLight;
	Y[4] = (X[4] - X[2]) / normalizedSpeedOfLight;
	Y[5] = X[6] + X[1];
	Y[6] = (X[7] - X[0]) / normalizedSpeedOfLight;
	Y[7] = (X[6] - X[1]) * normalizedSpeedOfLight;
	<? end ?>
<? if not useGEM then ?>
	for (int k = 0; k < 8; ++k) Y[k] = X[k];
<? end ?>

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
	
	<? if side == 0 then ?>
	UY.E = _real3(normalizedSpeedOfLightSq * solver->divPhiWavespeed * UX.phi, normalizedSpeedOfLightSq * UX.B.z, -normalizedSpeedOfLightSq * UX.B.y);
	UY.B = _real3(solver->divPsiWavespeed * UX.psi, -UX.E.z, UX.E.y);
	<? elseif side == 1 then ?>
	UY.E = _real3(-normalizedSpeedOfLightSq * UX.B.z, normalizedSpeedOfLightSq * solver->divPhiWavespeed * UX.phi, normalizedSpeedOfLightSq * UX.B.x);
	UY.B = _real3(UX.E.z, solver->divPsiWavespeed * UX.psi, -UX.E.x);
	<? elseif side == 2 then ?>
	UY.E = _real3(normalizedSpeedOfLightSq * UX.B.y, -normalizedSpeedOfLightSq * UX.B.x, normalizedSpeedOfLightSq * solver->divPhiWavespeed * UX.phi);
	UY.B = _real3(-UX.E.y, UX.E.x, solver->divPsiWavespeed * UX.psi);
	<? end ?>
	UY.phi = solver->divPhiWavespeed * UX.E.s<?=side?>;
	UY.psi = normalizedSpeedOfLightSq * solver->divPsiWavespeed * UX.B.s<?=side?>;

<? if useGEM then ?>
	<? if side == 0 then ?>
	UY.E_g = _real3(normalizedSpeedOfLightSq * solver->divPhiGWavespeed * UX.phi_g, normalizedSpeedOfLightSq * UX.B_g.z, -normalizedSpeedOfLightSq * UX.B_g.y);
	UY.B_g = _real3(solver->divPsiGWavespeed * UX.psi_g, -UX.E_g.z, UX.E_g.y);
	<? elseif side == 1 then ?>
	UY.E_g = _real3(-normalizedSpeedOfLightSq * UX.B_g.z, normalizedSpeedOfLightSq * solver->divPhiGWavespeed * UX.phi_g, normalizedSpeedOfLightSq * UX.B_g.x);
	UY.B_g = _real3(UX.E_g.z, solver->divPsiGWavespeed * UX.psi_g, -UX.E_g.x);
	<? elseif side == 2 then ?>
	UY.E_g = _real3(normalizedSpeedOfLightSq * UX.B_g.y, -normalizedSpeedOfLightSq * UX.B_g.x, normalizedSpeedOfLightSq * solver->divPhiGWavespeed * UX.phi_g);
	UY.B_g = _real3(-UX.E_g.y, UX.E_g.x, solver->divPsiGWavespeed * UX.psi_g);
	<? end ?>
	UY.phi_g = solver->divPhiGWavespeed * UX.E_g.s<?=side?>;
	UY.psi_g = normalizedSpeedOfLightSq * solver->divPsiGWavespeed * UX.B_g.s<?=side?>;
<? else ?>
	UY.E_g = UX.E_g;
	UY.B_g = UX.B_g;
	UY.phi_g = UX.phi_g;
	UY.psi_g = UX.psi_g;
<? end ?>

	return UY;
}
<? end ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

<? if useGEM then ?>
	//E_g has units m^2/s^3
	//1/c E_g has units m/s^2
	//B_g has units m/s^2
	//B_g * v/c has units m/s^2
	//rho * (E + v * B)/c has units kg/(m^2 s^2)
	real3 ionGravForce = _real3(
		(U->ion_rho * U->E_g.x + 4. * (U->ion_m.y * U->B_g.z - U->ion_m.z * U->B_g.y)) / normalizedSpeedOfLight,
		(U->ion_rho * U->E_g.y + 4. * (U->ion_m.z * U->B_g.x - U->ion_m.x * U->B_g.z)) / normalizedSpeedOfLight,
		(U->ion_rho * U->E_g.z + 4. * (U->ion_m.x * U->B_g.y - U->ion_m.y * U->B_g.x)) / normalizedSpeedOfLight);

	real3 elecGravForce = _real3(
		(U->elec_rho * U->E_g.x + 4. * (U->elec_m.y * U->B_g.z - U->elec_m.z * U->B_g.y)) / normalizedSpeedOfLight,
		(U->elec_rho * U->E_g.y + 4. * (U->elec_m.z * U->B_g.x - U->elec_m.x * U->B_g.z)) / normalizedSpeedOfLight,
		(U->elec_rho * U->E_g.z + 4. * (U->elec_m.x * U->B_g.y - U->elec_m.y * U->B_g.x)) / normalizedSpeedOfLight);
	
	const real _4_pi_G = 4. * M_PI * solver->gravitationalConstant;
	const real eps_g = 1. / _4_pi_G;
	const real mu_g = _4_pi_G / normalizedSpeedOfLightSq;

	//source of E_g is T_ti is the momentum + Poynting vector
	real3 T_ti = real3_zero;
	// I'm symmetrizing the stress-energy
		//matter
	T_ti.x -= U->ion_m.x * normalizedSpeedOfLight;
	T_ti.y -= U->ion_m.y * normalizedSpeedOfLight;
	T_ti.z -= U->ion_m.z * normalizedSpeedOfLight;
		//electromagnetism
	real3 S = real3_real_mul(real3_cross(U->E, U->B), 1. / solver->mu);
	real nSq = solver->eps * solver->mu;	// index of refraction^2 = 1 / phase vel^2 = permittivity * permeability
	T_ti.x -= .5 * S.x * (1. + nSq) / normalizedSpeedOfLight;
	T_ti.y -= .5 * S.y * (1. + nSq) / normalizedSpeedOfLight;
	T_ti.z -= .5 * S.z * (1. + nSq) / normalizedSpeedOfLight;
	
	deriv->E_g.x -= T_ti.x * mu_g / normalizedSpeedOfLight;
	deriv->E_g.y -= T_ti.y * mu_g / normalizedSpeedOfLight;
	deriv->E_g.z -= T_ti.z * mu_g / normalizedSpeedOfLight;

	// source of phi_g is T_00 is rho + .5 (E^2 + B^2)
	real T_tt = 0.;	
		//matter
	T_tt += (U->ion_rho + U->elec_rho) * normalizedSpeedOfLightSq;
		//electromagnetism			
	T_tt += calc_EM_energy(solver, U, x);
	
	deriv->phi_g += T_tt * _4_pi_G / normalizedSpeedOfLight * solver->divPhiGWavespeed;

	//ion_m has units kg/m^3 * m/s = kg/(m^2 s)
	//ion_m,t has units kg/(m^2 s^2)
	deriv->ion_m.x += (1. / normalizedIonLarmorRadius) * (U->ion_rho * U->E.x + U->ion_m.y * U->B.z - U->ion_m.z * U->B.y) + ionGravForce.x;
	deriv->ion_m.y += (1. / normalizedIonLarmorRadius) * (U->ion_rho * U->E.y + U->ion_m.z * U->B.x - U->ion_m.x * U->B.z) + ionGravForce.y;
	deriv->ion_m.z += (1. / normalizedIonLarmorRadius) * (U->ion_rho * U->E.z + U->ion_m.x * U->B.y - U->ion_m.y * U->B.x) + ionGravForce.z;
	deriv->ion_ETotal += (1. / normalizedIonLarmorRadius) * real3_dot(U->E, U->ion_m) + real3_dot(U->E_g, U->ion_m);

	deriv->elec_m.x -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * (U->elec_rho * U->E.x + U->elec_m.y * U->B.z - U->elec_m.z * U->B.y) + elecGravForce.x;
	deriv->elec_m.y -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * (U->elec_rho * U->E.y + U->elec_m.z * U->B.x - U->elec_m.x * U->B.z) + elecGravForce.y;
	deriv->elec_m.z -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * (U->elec_rho * U->E.z + U->elec_m.x * U->B.y - U->elec_m.y * U->B.x) + elecGravForce.z;
	deriv->elec_ETotal -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * real3_dot(U->E, U->elec_m) + real3_dot(U->E_g, U->elec_m);

<? else ?>

	deriv->ion_m.x += (1. / normalizedIonLarmorRadius) * (U->ion_rho * U->E.x + U->ion_m.y * U->B.z - U->ion_m.z * U->B.y);
	deriv->ion_m.y += (1. / normalizedIonLarmorRadius) * (U->ion_rho * U->E.y + U->ion_m.z * U->B.x - U->ion_m.x * U->B.z);
	deriv->ion_m.z += (1. / normalizedIonLarmorRadius) * (U->ion_rho * U->E.z + U->ion_m.x * U->B.y - U->ion_m.y * U->B.x);
	deriv->ion_ETotal += (1. / normalizedIonLarmorRadius) * real3_dot(U->E, U->ion_m);
	
	deriv->elec_m.x -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * (U->elec_rho * U->E.x + U->elec_m.y * U->B.z - U->elec_m.z * U->B.y);
	deriv->elec_m.y -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * (U->elec_rho * U->E.y + U->elec_m.z * U->B.x - U->elec_m.x * U->B.z);
	deriv->elec_m.z -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * (U->elec_rho * U->E.z + U->elec_m.x * U->B.y - U->elec_m.y * U->B.x);
	deriv->elec_ETotal -= solver->ionElectronMassRatio / normalizedIonLarmorRadius * real3_dot(U->E, U->elec_m);

<? end ?>
	
	real normalizedIonDebyeLengthSq	= normalizedIonDebyeLength * normalizedIonDebyeLength;
	//source of E is the current
	deriv->E.x -= (U->ion_m.x * solver->ionChargeMassRatio + U->elec_m.x * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius);
	deriv->E.y -= (U->ion_m.y * solver->ionChargeMassRatio + U->elec_m.y * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius);
	deriv->E.z -= (U->ion_m.z * solver->ionChargeMassRatio + U->elec_m.z * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius);
	//source of phi is the charge
	deriv->phi += (U->ion_rho * solver->ionChargeMassRatio + U->elec_rho * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius) * solver->divPhiWavespeed;

<? if not require 'coord.cartesian'.is(solver.coord) then ?>
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	prim_t W = primFromCons(solver, *U, x);
	real3 conn1_u = coord_conn_trace23(x);
	<? for _,fluid in ipairs(fluids) do ?>{
		real3 m_conn_vv = coord_conn_apply23(W.<?=fluid?>_v, U-><?=fluid?>_m, x);
		deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, m_conn_vv);	//-Conn^i_jk rho v^j v^k 
		deriv-><?=fluid?>_m = real3_sub(deriv-><?=fluid?>_m, real3_real_mul(conn1_u, W.<?=fluid?>_P));		//-Conn^i_jk g^jk P
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
