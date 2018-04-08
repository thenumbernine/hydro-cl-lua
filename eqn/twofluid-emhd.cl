#define sqrt_1_2 <?=('%.50f'):format(math.sqrt(.5))?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	<?=eqn.cons_t?> F;

<? 
for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_vj = W.<?=fluid?>_v.s<?=side?>;
	real <?=fluid?>_HTotal = U.<?=fluid?>_ETotal + W.<?=fluid?>_P;
	
	F.<?=fluid?>_rho = U.<?=fluid?>_m.s<?=side?>;
	F.<?=fluid?>_m = real3_scale(U.<?=fluid?>_m, <?=fluid?>_vj);
<? 	for i=0,2 do
?>	F.<?=fluid?>_m.s<?=i?> += coord_gU<?=i?><?=side?>(x) * W.<?=fluid?>_P;
<? 	end
?>	F.<?=fluid?>_ETotal = <?=fluid?>_HTotal * <?=fluid?>_vj;
	F.<?=fluid?>_ePot = 0.;
<? 
end
?>

	real3 B = U.B;
	real3 epsE = U.epsE;
	real mu = U.mu;
	real eps = U.eps;

	<? if side == 0 then ?>
	F.epsE = _real3(0., B.z / mu, -B.y / mu);
	F.B = _real3(0., -epsE.z / eps, epsE.y / eps);
	<? elseif side == 1 then ?>
	F.epsE = _real3(-B.z / mu, 0., B.x / mu);
	F.B = _real3(epsE.z / eps, 0., -epsE.x / eps);
	<? elseif side == 2 then ?>
	F.epsE = _real3(B.y / mu, -B.x / mu, 0.);
	F.B = _real3(-epsE.y / eps, epsE.x / eps, 0.);
	<? end ?>
	F.BPot = 0.;
	F.sigma = 0.;
	F.eps = 0.;
	F.mu = 0.;

	return F;
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	
	real EM_lambda = 1. / sqrt(U->eps * U->mu);
	range_t range = (range_t){-EM_lambda, EM_lambda};

<? 
for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_Cs = calc_<?=fluid?>_Cs(&W);
	real <?=fluid?>_Cs_sqrt_gU = <?=fluid?>_Cs * coord_sqrt_gU<?=side..side?>(x);
	range.min = min(range.min, W.<?=fluid?>_v.s<?=side?> - <?=fluid?>_Cs_sqrt_gU);
	range.max = max(range.max, W.<?=fluid?>_v.s<?=side?> + <?=fluid?>_Cs_sqrt_gU);
<? 
end
?>
	return range;
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	<?=eqn.prim_t?> WL = primFromCons(*UL, x);
	<?=eqn.prim_t?> WR = primFromCons(*UR, x);
	<?=eqn.eigen_t?> eig;

	eig.eps = .5 * (UL->eps + UR->eps);
	eig.mu = .5 * (UL->mu + UR->mu);

<? for _,fluid in ipairs(fluids) do ?>

	real <?=fluid?>_sqrtRhoL = sqrt(WL.<?=fluid?>_rho);
	real3 <?=fluid?>_vL = WL.<?=fluid?>_v;
	real <?=fluid?>_hTotalL = calc_hTotal(WL.<?=fluid?>_rho, WL.<?=fluid?>_P, UL-><?=fluid?>_ETotal) - UL-><?=fluid?>_ePot;
	
	real <?=fluid?>_sqrtRhoR = sqrt(UR-><?=fluid?>_rho);
	real3 <?=fluid?>_vR = WR.<?=fluid?>_v;
	real <?=fluid?>_hTotalR = calc_hTotal(WR.<?=fluid?>_rho, WR.<?=fluid?>_P, UR-><?=fluid?>_ETotal) - UR-><?=fluid?>_ePot;

	real <?=fluid?>_invDenom = 1./(<?=fluid?>_sqrtRhoL + <?=fluid?>_sqrtRhoR);
	
	//Roe-averaged
	real <?=fluid?>_rho = <?=fluid?>_sqrtRhoL * <?=fluid?>_sqrtRhoR;
	real3 <?=fluid?>_v = real3_add(
			real3_scale(<?=fluid?>_vL, <?=fluid?>_sqrtRhoL * <?=fluid?>_invDenom),
			real3_scale(<?=fluid?>_vR, <?=fluid?>_sqrtRhoR * <?=fluid?>_invDenom));
	real <?=fluid?>_hTotal = <?=fluid?>_invDenom * (<?=fluid?>_sqrtRhoL * <?=fluid?>_hTotalL + <?=fluid?>_sqrtRhoR * <?=fluid?>_hTotalR);

	//derived:
	real <?=fluid?>_vSq = coordLenSq(<?=fluid?>_v, x);
	real <?=fluid?>_eKin = .5 * <?=fluid?>_vSq;
	real <?=fluid?>_CsSq = (heatCapacityRatio - 1.) * (<?=fluid?>_hTotal - <?=fluid?>_eKin);
	real <?=fluid?>_Cs = sqrt(<?=fluid?>_CsSq);
	
	eig.<?=fluid?>_rho = <?=fluid?>_rho; 
	eig.<?=fluid?>_v = <?=fluid?>_v;
	eig.<?=fluid?>_hTotal = <?=fluid?>_hTotal;
	eig.<?=fluid?>_vSq = <?=fluid?>_vSq;
	eig.<?=fluid?>_Cs = <?=fluid?>_Cs;

<? end ?>
	
	return eig;
}
<? end ?>

//same as in eqn/euler.cl
kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,		//[volume][dim]
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
	
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		*eig = eigen_forSide_<?=side?>(UL, UR, xInt);
	}<? end ?>
}

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,solver.dim-1 do 
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
	real3 ion_v = eig->ion_v;
	real3 ion_vL = coord_lower(ion_v, x);
	real ion_hTotal = eig->ion_hTotal;
	real ion_vSq = real3_dot(ion_v, ion_vL);
	real ion_Cs = eig->ion_Cs;
	real ion_Cs_over_sqrt_gUjj = ion_Cs / sqrt_gUjj; 
	
	real3 elec_v = eig->elec_v;
	real3 elec_vL = coord_lower(elec_v, x);
	real elec_hTotal = eig->elec_hTotal;
	real elec_vSq = real3_dot(elec_v, elec_vL);
	real elec_Cs = eig->elec_Cs;
	real elec_Cs_over_sqrt_gUjj = elec_Cs / sqrt_gUjj; 

]] .. prefix

				local gUdef = '\treal3 gUj = _real3(\n'
				for i=0,2 do
					gUdef = gUdef .. '\t\tcoord_gU'..side..i..'(x)'..(i<2 and ',' or '')..'\n'
				end
				gUdef = gUdef .. '\t);\n'
				prefix = gUdef .. prefix
?>

void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) { 
	<?=prefix?>

	real ion_denom = 2. * ion_Cs * ion_Cs;
	real ion_invDenom = 1. / ion_denom;
	
	real elec_denom = 2. * elec_Cs * elec_Cs;
	real elec_invDenom = 1. / elec_denom;

	const real heatRatioMinusOne = heatCapacityRatio - 1.;
	
	const real ise = sqrt_1_2 / sqrt(eig->eps);
	const real isu = sqrt_1_2 / sqrt(eig->mu);

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
	Y[10] = X[12] *  ise + X[14] * isu;
	Y[11] = X[11] * -ise + X[15] * isu;
	Y[12] = X[10] * -ise + X[13] * isu;
	Y[13] = X[10] *  ise + X[13] * isu;
	Y[14] = X[11] *  ise + X[15] * isu;
	Y[15] = X[12] * -ise + X[14] * isu;
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
	Y[10] = X[10] *  ise + X[15] * isu;
	Y[11] = X[12] * -ise + X[13] * isu;
	Y[12] = X[11] * -ise + X[14] * isu;
	Y[13] = X[11] *  ise + X[14] * isu;
	Y[14] = X[12] *  ise + X[13] * isu;
	Y[15] = X[10] * -ise + X[15] * isu;
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
	Y[10] = X[11] *  ise + X[13] * isu;
	Y[11] = X[10] * -ise + X[14] * isu;
	Y[12] = X[12] * -ise + X[15] * isu;
	Y[13] = X[12] *  ise + X[15] * isu;
	Y[14] = X[10] *  ise + X[14] * isu;
	Y[15] = X[11] * -ise + X[13] * isu;
<?
				end 
?>
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,			//numIntStates = 16
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,	//numWaves = 16
	real3 x
) {
	<?=prefix?>
	const real se = sqrt_1_2 * sqrt(eig->eps);
	const real su = sqrt_1_2 * sqrt(eig->mu);
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
	Y[10] = se * (-X[12] + X[13]);
	Y[11] = se * (-X[11] + X[14]);
	Y[12] = se * (X[10] + -X[15]);
	Y[13] = su * (X[12] + X[13]);
	Y[14] = su * (X[10] + X[15]);
	Y[15] = su * (X[11] + X[14]);
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
	Y[10] = se * (X[10] - X[15]);
	Y[11] = se * (-X[12] + X[13]);
	Y[12] = se * (-X[11] + X[14]);
	Y[13] = su * (X[11] + X[14]);
	Y[14] = su * (X[12] + X[13]);
	Y[15] = su * (X[10] + X[15]);

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
	Y[10] = se * (-X[11] + X[14]);
	Y[11] = se * (X[10] - X[15]);
	Y[12] = se * (-X[12] + X[13]);
	Y[13] = su * (X[10] + X[15]);
	Y[14] = su * (X[11] + X[14]);
	Y[15] = su * (X[12] + X[13]);
<?
				end
?>
}

<?
				if solver.checkFluxError then 
?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) {
	<?=prefix?>
<?
					for i,fluid	in ipairs(fluids) do 
?>
	Y[<?=5*i-5?>] = X[<?=5*i-4?>] * nx 
		+ X[<?=5*i-3?>] * ny 
		+ X[<?=5*i-2?>] * nz;
	Y[<?=5*i-4?>] = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.x + (heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * gUj.x)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.x * nx - (heatCapacityRatio - 1.) * gUj.x * <?=fluid?>_vL.x + <?=fluid?>_v_n)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.x * ny - (heatCapacityRatio - 1.) * gUj.x * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.x * nz - (heatCapacityRatio - 1.) * gUj.x * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (heatCapacityRatio - 1.) * nx;
	Y[<?=5*i-3?>] = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.y + (heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * gUj.y)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.y * nx - (heatCapacityRatio - 1.) * gUj.y * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.y * ny - (heatCapacityRatio - 1.) * gUj.y * <?=fluid?>_vL.y + <?=fluid?>_v_n)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.y * nz - (heatCapacityRatio - 1.) * gUj.y * <?=fluid?>_vL.z)
		+ X[<?=5*i-1?>] * (heatCapacityRatio - 1.) * ny;
	Y[<?=5*i-2?>] = X[<?=5*i-5?>] * (-<?=fluid?>_v_n * <?=fluid?>_v.z + (heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq * gUj.z)
		+ X[<?=5*i-4?>] * (<?=fluid?>_v.z * nx - (heatCapacityRatio - 1.) * gUj.z * <?=fluid?>_vL.x)
		+ X[<?=5*i-3?>] * (<?=fluid?>_v.z * ny - (heatCapacityRatio - 1.) * gUj.z * <?=fluid?>_vL.y)
		+ X[<?=5*i-2?>] * (<?=fluid?>_v.z * nz - (heatCapacityRatio - 1.) * gUj.z * <?=fluid?>_vL.z + <?=fluid?>_v_n)
		+ X[<?=5*i-1?>] * (heatCapacityRatio - 1.) * nz;
	Y[<?=5*i-1?>] = X[<?=5*i-5?>] * <?=fluid?>_v_n * ((heatCapacityRatio - 1.) * .5 * <?=fluid?>_vSq - <?=fluid?>_hTotal)
		+ X[<?=5*i-4?>] * (-(heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.x + nx * <?=fluid?>_hTotal)
		+ X[<?=5*i-3?>] * (-(heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.y + ny * <?=fluid?>_hTotal)
		+ X[<?=5*i-2?>] * (-(heatCapacityRatio - 1.) * <?=fluid?>_v_n * <?=fluid?>_vL.z + nz * <?=fluid?>_hTotal)
		+ X[<?=5*i-1?>] * heatCapacityRatio * <?=fluid?>_v_n;
<? 
					end 
?>
	real3 epsE = ((<?=addr2?> const <?=eqn.cons_t?>*)X)->epsE;
	real3 B = ((<?=addr2?> const <?=eqn.cons_t?>*)X)->B;
	real eps = eig->eps;
	real mu = eig->mu;
<? 
					if side==0 then 
?>
	Y[10] = 0;
	Y[11] = B.z / mu;
	Y[12] = -B.y / mu;
	Y[13] = 0;
	Y[14] = -epsE.z / eps;
	Y[15] = epsE.y / eps;
<?
					elseif side==1 then
?>
	Y[10] = -B.z / mu;
	Y[11] = 0;
	Y[12] = B.x / mu;
	Y[13] = epsE.z / eps;
	Y[14] = 0;
	Y[15] = -epsE.x / eps;
<?
					elseif side==2 then 
?>
	Y[10] = B.y / mu;
	Y[11] = -B.x / mu;
	Y[12] = 0;
	Y[13] = -epsE.y / eps;
	Y[14] = epsE.x / eps;
	Y[15] = 0;
<?
					end 
?>
}
<?
				end
			end
		end
	end
end
?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	deriv->ion_m.x += (1. / ionLarmorRadius) * (U->ion_rho * U->epsE.x / U->eps + U->ion_m.y * U->B.z - U->ion_m.z * U->B.y);
	deriv->ion_m.y += (1. / ionLarmorRadius) * (U->ion_rho * U->epsE.y / U->eps + U->ion_m.z * U->B.x - U->ion_m.x * U->B.z);
	deriv->ion_m.z += (1. / ionLarmorRadius) * (U->ion_rho * U->epsE.z / U->eps + U->ion_m.x * U->B.y - U->ion_m.y * U->B.x);
	deriv->ion_ETotal += (1. / ionLarmorRadius) * real3_dot(U->epsE, U->ion_m) / U->eps;
	
	deriv->elec_m.x -= ionElectronMassRatio / ionLarmorRadius * (U->elec_rho * U->epsE.x / U->eps + U->elec_m.y * U->B.z - U->elec_m.z * U->B.y);
	deriv->elec_m.y -= ionElectronMassRatio / ionLarmorRadius * (U->elec_rho * U->epsE.y / U->eps + U->elec_m.z * U->B.x - U->elec_m.x * U->B.z);
	deriv->elec_m.z -= ionElectronMassRatio / ionLarmorRadius * (U->elec_rho * U->epsE.z / U->eps + U->elec_m.x * U->B.y - U->elec_m.y * U->B.x);
	deriv->elec_ETotal -= ionElectronMassRatio / ionLarmorRadius * real3_dot(U->epsE, U->elec_m) / U->eps;

#define ionDebyeLengthSq (ionDebyeLength*ionDebyeLength)
	deriv->epsE.x -= (U->ion_m.x * ionChargeMassRatio
						+ U->elec_m.x * elecChargeMassRatio
					) / ionDebyeLengthSq * ionLarmorRadius;
	deriv->epsE.y -= (U->ion_m.y * ionChargeMassRatio
						+ U->elec_m.y * elecChargeMassRatio
					) / (ionDebyeLengthSq * ionLarmorRadius);
	deriv->epsE.z -= (U->ion_m.z * ionChargeMassRatio
						+ U->elec_m.z * elecChargeMassRatio
					) / (ionDebyeLengthSq * ionLarmorRadius);
}
