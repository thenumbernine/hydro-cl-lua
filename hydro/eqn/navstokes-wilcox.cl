<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real vTilde_j = W.vTilde.s<?=side?>;
	
	//this is the flux term used, but is it technically called 'HTotal' ?
	//'HTotal' = rhoBar (1/2 vTilde^2 + (1 - 2/3/gamma_1) k) + 1/gamma_1 PStar ... + PStar
	//'HTotal' = rhoBar (1/2 vTilde^2 + (1 - 2/3 C_v/R) k) + (1 + C_v/R) PStar
	//'hTotal' = 1/2 vTilde^2 + (1 - 2/3 C_v/R) k + (1 + C_v/R) PStar / rhoBar
	real HTotal = U.rhoBar_eTotalTilde + W.PStar;
	
	<?=eqn.cons_t?> F;
	
	F.rhoBar = U.rhoBar_vTilde.s<?=side?>;
	
	F.rhoBar_vTilde = real3_real_mul(U.rhoBar_vTilde, vTilde_j);

<? for i=0,2 do
?>	F.rhoBar_vTilde.s<?=i?> += coord_g_uu<?=i?><?=side?>(x) * W.PStar;
<? end
?>	
	F.rhoBar_eTotalTilde = HTotal * vTilde_j;
	
	F.rhoBar_k = 0;
	F.rhoBar_omega = 0;

	F.ePot = 0;
	
	return F;
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real Cs = calc_Cs(solver, W);
	real Cs_sqrt_gU = Cs * coord_sqrt_g_uu<?=side..side?>(x);
	return (range_t){
		.min = W.vTilde.s<?=side?> - Cs_sqrt_gU, 
		.max = W.vTilde.s<?=side?> + Cs_sqrt_gU,
	};
}
<? end ?>

<?=eqn.eigen_t?> eigen_forInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> WL = primFromCons(solver, UL, x);
	real sqrtRhoL = sqrt(WL.rhoBar);
	real3 vTildeL = WL.vTilde;
	real hTotalL = (WL.PStar + UL.rhoBar_eTotalTilde) / UL.rhoBar - UL.ePot;
	real kL = WL.k;
	real omegaL = WL.omega;

	<?=eqn.prim_t?> WR = primFromCons(solver, UR, x);
	real sqrtRhoR = sqrt(WR.rhoBar);
	real3 vTildeR = WR.vTilde;
	real hTotalR = (WR.PStar + UR.rhoBar_eTotalTilde) / UR.rhoBar - UR.ePot;
	real kR = WR.k;
	real omegaR = WR.omega;

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);
	
	//Roe-averaged
	real rhoBar = sqrtRhoL * sqrtRhoR;
	real3 vTilde = real3_add(
			real3_real_mul(vTildeL, sqrtRhoL * invDenom),
			real3_real_mul(vTildeR, sqrtRhoR * invDenom));
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);
	//I'm making this part up:
	real k = invDenom * (sqrtRhoL * kL + sqrtRhoR * kR);
	real omega = invDenom * (sqrtRhoL * omegaL + sqrtRhoR * omegaR);

	//derived:
	real vTildeSq = coordLenSq(vTilde, x);
	real eKinTilde = .5 * vTildeSq;
/*	
rhoBar eIntTilde = rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot
rhoBar TTilde C_v = rhoBar eIntTilde
rhoBar TTilde = 1/C_v ( rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot )
PBar / R = rhoBar TTilde
PBar = R/C_v ( rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot )
PBar = PStar - 2/3 rhoBar k
PStar - 2/3 rhoBar k = R/C_v ( rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot )

rhoBar eTotalTilde = rhoBar (hTotalTilde + ePot) - PStar
(R/C_v + 1) PStar / rhoBar = R/C_v (hTotalTilde - 1/2 vTilde^2 - (1 - 2/3 C_v/R) k)

Cs^2 = (R/C_v+1) PStar / rhoBar
Cs^2 = R/C_v (hTotalTilde - 1/2 vTilde^2 - (1 - 2/3 C_v/R) k)
*/	
	real CsSq = R_over_C_v * (hTotal - eKinTilde - k) + 2./3. * k;
	real Cs = sqrt(CsSq);

	return (<?=eqn.eigen_t?>){
		.rhoBar = rhoBar, 
		.vTilde = vTilde,
		.hTotal = hTotal,
		.k = k,
		.omega = omega,
		.vTildeSq = vTildeSq,
		.Cs = Cs,
	};
}

<?
local prefixes = {}
for side=0,2 do 
	local prefix	
	if side == 0 then
		prefix = [[
	const real nx = 1, ny = 0, nz = 0;
	const real n1x = 0, n1y = 1, n1z = 0;
	const real n2x = 0, n2y = 0, n2z = 1;
	real vTilde_n = vTilde.x, vTilde_n1 = vTilde.y, vTilde_n2 = vTilde.z;
]] 
	elseif side == 1 then
		prefix = [[
	const real nx = 0, ny = 1, nz = 0;
	const real n1x = 0, n1y = 0, n1z = 1;
	const real n2x = 1, n2y = 0, n2z = 0;
	real vTilde_n = vTilde.y, vTilde_n1 = vTilde.z, vTilde_n2 = vTilde.x;
]] 
	elseif side == 2 then
		prefix = [[
	const real nx = 0, ny = 0, nz = 1;
	const real n1x = 1, n1y = 0, n1z = 0;
	const real n2x = 0, n2y = 1, n2z = 0;
	real vTilde_n = vTilde.z, vTilde_n1 = vTilde.x, vTilde_n2 = vTilde.y;
]]
	end
	
	prefix = [[
	sym3 gU = coord_g_uu(x);
	real gUjj = gU.s]]..side..side..[[;
	real sqrt_gUjj = coord_sqrt_g_uu]]..side..side..[[(x);
	
	real3 vTilde = eig.vTilde;
	real3 vTildeL = coord_lower(vTilde, x);
	real hTotal = eig.hTotal;
	real vTildeSq = real3_dot(vTilde, vTildeL);
	real Cs = eig.Cs;
	real Cs_over_sqrt_gUjj = Cs / sqrt_gUjj; 
	real rhoBar = eig.rhoBar;
	real k = eig.k;
	real omega = eig.omega;
	//g^ij for fixed j=side
]] .. prefix

	local gUdef = '\treal3 gUj = _real3(\n'
	for i=0,2 do
		gUdef = gUdef .. '\t\tcoord_g_uu'..side..i..'(x)'..(i<2 and ',' or '')..'\n'
	end
	gUdef = gUdef .. '\t);\n'
	prefix = gUdef .. prefix
	prefixes[side] = prefix
end
?>

<?=eqn.waves_t?> eigen_leftTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x,
	normal_t n
) { 
	if (n.side == 0) {
		<?=prefixes[0]?>
		real denom = 2. * Cs * Cs;
		real invDenom = 1. / denom;
		real invRhoBar = 1. / eig.rhoBar;
		const real tmp = 1. - 2./3. * C_v_over_R;
	
		real sqrt_gUxx = sqrt_gUjj;
		return (<?=eqn.waves_t?>){.ptr={
			(X.ptr[0] * (.5 * R_over_C_v * vTildeSq + Cs * vTilde.x / sqrt_gUxx)
				+ X.ptr[1] * (-R_over_C_v * vTildeL.x - Cs / sqrt_gUxx)
				+ X.ptr[2] * -R_over_C_v * vTildeL.y
				+ X.ptr[3] * -R_over_C_v * vTildeL.z
				+ X.ptr[4] * R_over_C_v
				+ X.ptr[5] * (2./3. - R_over_C_v)
			) * invDenom,
			(X.ptr[0] * (denom - R_over_C_v * vTildeSq)
				+ X.ptr[1] * 2. * R_over_C_v * vTildeL.x
				+ X.ptr[2] * 2. * R_over_C_v * vTildeL.y
				+ X.ptr[3] * 2. * R_over_C_v * vTildeL.z
				+ X.ptr[4] * -2. * R_over_C_v
				+ X.ptr[5] * 2. * (R_over_C_v - 2./3.)
			) * invDenom,
			(X.ptr[0] * (vTilde.x * gU.xy / gU.xx - vTilde.y)
				+ X.ptr[1] * -gU.xy / gU.xx
				+ X.ptr[2]) * invRhoBar,
			(X.ptr[0] * (vTilde.x * gU.xz / gU.xx - vTilde.z)
				+ X.ptr[1] * -gU.xz / gU.xx
				+ X.ptr[3]) * invRhoBar,
			(X.ptr[0] * -k
				+ X.ptr[5]) * invRhoBar,
			(X.ptr[0] * -omega
				+ X.ptr[6]) * invRhoBar,
			(X.ptr[0] * (.5 * R_over_C_v * vTildeSq - Cs * vTilde.x / sqrt_gUxx)
				+ X.ptr[1] * (-R_over_C_v * vTildeL.x + Cs / sqrt_gUxx)
				+ X.ptr[2] * -R_over_C_v * vTildeL.y
				+ X.ptr[3] * -R_over_C_v * vTildeL.z
				+ X.ptr[4] * R_over_C_v
				+ X.ptr[5] * (2./3. - R_over_C_v)
			) * invDenom,
		}};
	} else if (n.side == 1) {
		<?=prefixes[1]?>
		real denom = 2. * Cs * Cs;
		real invDenom = 1. / denom;
		real invRhoBar = 1. / eig.rhoBar;
		const real tmp = 1. - 2./3. * C_v_over_R;
		
		real sqrt_gUyy = sqrt_gUjj;
		return (<?=eqn.waves_t?>){.ptr={
			(X.ptr[0] * (.5 * R_over_C_v * vTildeSq + Cs * vTilde.y / sqrt_gUyy)
				+ X.ptr[1] * -R_over_C_v * vTildeL.x
				+ X.ptr[2] * (-R_over_C_v * vTildeL.y - Cs / sqrt_gUyy)
				+ X.ptr[3] * -R_over_C_v * vTildeL.z
				+ X.ptr[4] * R_over_C_v
				+ X.ptr[5] * (2./3. - R_over_C_v)
			) * invDenom,
			(X.ptr[0] * (vTilde.y * gU.xy / gU.yy - vTilde.x)
				+ X.ptr[1]
				+ X.ptr[2] * -gU.xy / gU.yy) * invRhoBar,
			(X.ptr[0] * (denom - R_over_C_v * vTildeSq)
				+ X.ptr[1] * 2. * R_over_C_v * vTildeL.x
				+ X.ptr[2] * 2. * R_over_C_v * vTildeL.y
				+ X.ptr[3] * 2. * R_over_C_v * vTildeL.z
				+ X.ptr[4] * -2. * R_over_C_v
				+ X.ptr[5] * 2. * (R_over_C_v - 2./3.)
			) * invDenom,
			(X.ptr[0] * (vTilde.y * gU.yz / gU.yy - vTilde.z)
				+ X.ptr[2] * -gU.yz / gU.yy
				+ X.ptr[3]) * invRhoBar,
			(X.ptr[0] * -k
				+ X.ptr[5]) * invRhoBar,
			(X.ptr[0] * -omega
				+ X.ptr[6]) * invRhoBar,
			(X.ptr[0] * (.5 * R_over_C_v * vTildeSq - Cs * vTilde.y / sqrt_gUyy)
				+ X.ptr[1] * -R_over_C_v * vTildeL.x
				+ X.ptr[2] * (-R_over_C_v * vTildeL.y + Cs / sqrt_gUyy)
				+ X.ptr[3] * -R_over_C_v * vTildeL.z
				+ X.ptr[4] * R_over_C_v
				+ X.ptr[5] * (2./3. - R_over_C_v)
			) * invDenom,
		}};
	} else if (n.side == 2) {
		<?=prefixes[2]?>
		real denom = 2. * Cs * Cs;
		real invDenom = 1. / denom;
		real invRhoBar = 1. / eig.rhoBar;
		const real tmp = 1. - 2./3. * C_v_over_R;
		
		real sqrt_gUzz = sqrt_gUjj;
		return (<?=eqn.waves_t?>){.ptr={
			(X.ptr[0] * (.5 * R_over_C_v * vTildeSq + Cs * vTilde.z / sqrt_gUzz)
				+ X.ptr[1] * -R_over_C_v * vTildeL.x
				+ X.ptr[2] * -R_over_C_v * vTildeL.y
				+ X.ptr[3] * (-R_over_C_v * vTildeL.z - Cs / sqrt_gUzz)
				+ X.ptr[4] * R_over_C_v
				+ X.ptr[5] * (2./3. - R_over_C_v)
			) * invDenom,
			(X.ptr[0] * (vTilde.z * gU.xz / gU.zz - vTilde.x)
				+ X.ptr[1]
				+ X.ptr[3] * -gU.xz / gU.zz) * invRhoBar,
			(X.ptr[0] * (vTilde.z * gU.yz / gU.zz - vTilde.y)
				+ X.ptr[2]
				+ X.ptr[3] * -gU.yz / gU.zz) * invRhoBar,
			(X.ptr[0] * (denom - R_over_C_v * vTildeSq)
				+ X.ptr[1] * 2. * R_over_C_v * vTildeL.x
				+ X.ptr[2] * 2. * R_over_C_v * vTildeL.y
				+ X.ptr[3] * 2. * R_over_C_v * vTildeL.z
				+ X.ptr[4] * -2. * R_over_C_v
				+ X.ptr[5] * 2. * (R_over_C_v - 2./3.)
			) * invDenom,
			(X.ptr[0] * -k
				+ X.ptr[5]) * invRhoBar,
			(X.ptr[0] * -omega
				+ X.ptr[6]) * invRhoBar,
			(X.ptr[0] * (.5 * R_over_C_v * vTildeSq - Cs * vTilde.z / sqrt_gUzz)
				+ X.ptr[1] * -R_over_C_v * vTildeL.x
				+ X.ptr[2] * -R_over_C_v * vTildeL.y
				+ X.ptr[3] * (-R_over_C_v * vTildeL.z + Cs / sqrt_gUzz)
				+ X.ptr[4] * R_over_C_v
				+ X.ptr[5] * (2./3. - R_over_C_v)
			) * invDenom,
		}};
	}
}

<?=eqn.cons_t?> eigen_rightTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x,
	normal_t n
) {
	if (n.side == 0) {
		<?=prefixes[0]?>
		real sqrt_gUxx = sqrt_gUjj;
		return (<?=eqn.cons_t?>){.ptr={
			X.ptr[0] + X.ptr[1] + X.ptr[6],
			X.ptr[0] * (vTilde.x - Cs * sqrt_gUxx)
				+ X.ptr[1] * vTilde.x
				+ X.ptr[6] * (vTilde.x + Cs * sqrt_gUxx),
			X.ptr[0] * (vTilde.y - Cs * gU.xy / sqrt_gUxx)
				+ X.ptr[1] * vTilde.y
				+ X.ptr[2] * rhoBar
				+ X.ptr[6] * (vTilde.y + Cs * gU.xy / sqrt_gUxx),
			X.ptr[0] * (vTilde.z - Cs * gU.xz / sqrt_gUxx)
				+ X.ptr[1] * vTilde.z
				+ X.ptr[3] * rhoBar
				+ X.ptr[6] * (vTilde.z + Cs * gU.xz / sqrt_gUxx),
			X.ptr[0] * (hTotal - Cs * vTilde.x / sqrt_gUxx)
				+ X.ptr[1] * (.5 * vTildeSq + (1. - 2./3. * C_v_over_R) * k)
				+ X.ptr[2] * rhoBar * vTildeL.y
				+ X.ptr[3] * rhoBar * vTildeL.z
				+ X.ptr[4] * rhoBar * (1. - 2./3. * C_v_over_R)
				+ X.ptr[6] * (hTotal + Cs * vTilde.x / sqrt_gUxx),
			(X.ptr[0] + X.ptr[1] + X.ptr[6]) * k
				+ X.ptr[4] * rhoBar,
			(X.ptr[0] + X.ptr[1] + X.ptr[6]) * omega
				+ X.ptr[5] * rhoBar,
			0,
		}};
	} else if (n.side == 1) {
		<?=prefixes[1]?>
		real sqrt_gUyy = sqrt_gUjj;
		return (<?=eqn.cons_t?>){.ptr={
			X.ptr[0] + X.ptr[2] + X.ptr[6],
			X.ptr[0] * (vTilde.x - Cs * gU.xy / sqrt_gUyy)
				+ X.ptr[1] * rhoBar
				+ X.ptr[2] * vTilde.x
				+ X.ptr[6] * (vTilde.x + Cs * gU.xy / sqrt_gUyy),
			X.ptr[0] * (vTilde.y - Cs * sqrt_gUyy)
				+ X.ptr[2] * vTilde.y
				+ X.ptr[6] * (vTilde.y + Cs * sqrt_gUyy),
			X.ptr[0] * (vTilde.z - Cs * gU.yz / sqrt_gUyy)
				+ X.ptr[2] * vTilde.z
				+ X.ptr[3] * rhoBar
				+ X.ptr[6] * (vTilde.z + Cs * gU.yz / sqrt_gUyy),
			X.ptr[0] * (hTotal - Cs * vTilde.y / sqrt_gUyy)
				+ X.ptr[1] * rhoBar * vTildeL.x
				+ X.ptr[2] * (.5 * vTildeSq + (1. - 2./3. * C_v_over_R) * k)
				+ X.ptr[3] * rhoBar * vTildeL.z
				+ X.ptr[4] * (1. - 2./3. * C_v_over_R) * k
				+ X.ptr[6] * (hTotal + Cs * vTilde.y / sqrt_gUyy),
			(X.ptr[0] + X.ptr[2] + X.ptr[6]) * k
				+ X.ptr[4] * rhoBar,
			(X.ptr[0] + X.ptr[2] + X.ptr[6]) * omega
				+ X.ptr[5] * rhoBar,
			0,
		}};
	} else if (n.side == 2) {
		<?=prefixes[2]?>
		real sqrt_gUzz = sqrt_gUjj;
		return (<?=eqn.cons_t?>){.ptr={
			X.ptr[0] + X.ptr[3] + X.ptr[6],
			X.ptr[0] * (vTilde.x - Cs * gU.xz / sqrt_gUzz)
				+ X.ptr[1] * rhoBar
				+ X.ptr[3] * vTilde.x
				+ X.ptr[6] * (vTilde.x + Cs * gU.xz / sqrt_gUzz),
			X.ptr[0] * (vTilde.y - Cs * gU.yz / sqrt_gUzz)
				+ X.ptr[2] * rhoBar
				+ X.ptr[3] * vTilde.y
				+ X.ptr[6] * (vTilde.y + Cs * gU.yz / sqrt_gUzz),
			X.ptr[0] * (vTilde.z - Cs * sqrt_gUzz)
				+ X.ptr[3] * vTilde.z
				+ X.ptr[6] * (vTilde.z + Cs * sqrt_gUzz),
			X.ptr[0] * (hTotal - Cs * vTilde.z / sqrt_gUzz)
				+ X.ptr[1] * rhoBar * vTildeL.x
				+ X.ptr[2] * rhoBar * vTildeL.y
				+ X.ptr[3] * (.5 * vTildeSq + (1. - 2./3. * C_v_over_R) * k)
				+ X.ptr[4] * (1. - 2./3. * C_v_over_R) * k
				+ X.ptr[6] * (hTotal + Cs * vTilde.z / sqrt_gUzz),
			(X.ptr[0] + X.ptr[3] + X.ptr[6]) * k
				+ X.ptr[4] * rhoBar,
			(X.ptr[0] + X.ptr[3] + X.ptr[6]) * omega
				+ X.ptr[5] * rhoBar,
			0,
		}};
	}
}

<?=eqn.cons_t?> eigen_fluxTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x,
	normal_t n
) {
	if (n.side == 0) {
		<?=prefixes[0]?>
		real CsSq = Cs * Cs;
		real PStar = CsSq * rhoBar / (1. + R_over_C_v);
		return (<?=eqn.cons_t?>){.ptr={
			X.ptr[1],
			- X.ptr[0] * (vTilde.x * vTilde.x - .5 * gUj.x * R_over_C_v * vTildeSq)
				+ X.ptr[1] * (2 * vTilde.x - gUj.x * R_over_C_v * vTildeL.x) 
				- X.ptr[2] * gUj.x * R_over_C_v * vTildeL.y 
				- X.ptr[3] * gUj.x * R_over_C_v * vTildeL.z 
				+ X.ptr[4] * gUj.x * R_over_C_v
				- X.ptr[5] * gUj.x * (R_over_C_v - 2./.3),
			- X.ptr[0] * (vTilde.x * vTilde.y - .5 * gUj.y * R_over_C_v * vTildeSq)
				+ X.ptr[1] * (vTilde.y - gUj.y * R_over_C_v * vTildeL.x)
				+ X.ptr[2] * (vTilde.x - gUj.y * R_over_C_v * vTildeL.y) 
				- X.ptr[3] * gUj.y * R_over_C_v * vTildeL.z 
				+ X.ptr[4] * gUj.y * R_over_C_v
				- X.ptr[5] * gUj.y * (R_over_C_v - 2./3.),
			- X.ptr[0] * (vTilde.x * vTilde.z - .5 * gUj.z * R_over_C_v * vTildeSq)
				+ X.ptr[1] * (vTilde.z - gUj.z * R_over_C_v * vTildeL.x)
				- X.ptr[2] * gUj.z * R_over_C_v * vTildeL.y 
				+ X.ptr[3] * (vTilde.x - gUj.z * R_over_C_v * vTildeL.z) 
				+ X.ptr[4] * gUj.z * R_over_C_v
				- X.ptr[5] * gUj.z * (R_over_C_v - 2./3.),
			+ X.ptr[0] * vTilde.x * (.5 * R_over_C_v * vTildeSq - hTotal)
				- X.ptr[1] * (R_over_C_v * vTildeL.x * vTilde.x - hTotal)
				- X.ptr[2] * vTilde.x * R_over_C_v * vTildeL.y
				- X.ptr[3] * vTilde.x * R_over_C_v * vTildeL.z
				+ X.ptr[4] * vTilde.x * (R_over_C_v + 1.)
				- X.ptr[5] * vTilde.x * (R_over_C_v - 2./3.),
			- X.ptr[0] * k * vTilde.x
				+ X.ptr[1] * k
				+ X.ptr[5] * vTilde.x,
			- X.ptr[0] * vTilde.x * omega
				+ X.ptr[1] * omega
				+ X.ptr[6] * vTilde.x,
			0.
		}};
	} else if (n.side == 1) {
		<?=prefixes[1]?>
		real CsSq = Cs * Cs;
		real PStar = CsSq * rhoBar / (1. + R_over_C_v);
		return (<?=eqn.cons_t?>){.ptr={
			X.ptr[2],
			-X.ptr[0] * (vTilde.x * vTilde.y - .5 * gUj.x * R_over_C_v * vTildeSq)
				+ X.ptr[1] * (vTilde.y - gUj.x * R_over_C_v * vTildeL.x)
				+ X.ptr[2] * (vTilde.x - gUj.x * R_over_C_v * vTildeL.y) 
				- X.ptr[3] * gUj.x * R_over_C_v * vTildeL.z 
				+ X.ptr[4] * gUj.x * R_over_C_v
				- X.ptr[5] * gUj.x * (R_over_C_v - 2./3.),
			-X.ptr[0] * (vTilde.y * vTilde.y - .5 * gUj.y * R_over_C_v * vTildeSq)
				- X.ptr[1] * gUj.y * R_over_C_v * vTildeL.x 
				+ X.ptr[2] * (2 * vTilde.y - gUj.y * R_over_C_v * vTildeL.y) 
				- X.ptr[3] * gUj.y * R_over_C_v * vTildeL.z
				+ X.ptr[4] * gUj.y * R_over_C_v
				- X.ptr[5] * gUj.y * (R_over_C_v - 2./3.),
			-X.ptr[0] * (vTilde.y * vTilde.z - .5 * gUj.z * R_over_C_v * vTildeSq)
				- X.ptr[1] * gUj.z * R_over_C_v * vTildeL.x 
				+ X.ptr[2] * (vTilde.z - gUj.z * R_over_C_v * vTildeL.y)
				+ X.ptr[3] * (vTilde.y - gUj.z * R_over_C_v * vTildeL.z)
				+ X.ptr[4] * gUj.z * R_over_C_v
				- X.ptr[5] * gUj.z * (R_over_C_v - 2./3.),
			X.ptr[0] * vTilde.y * (.5 * R_over_C_v * vTildeSq - hTotal)
				- X.ptr[1] * vTildeL.x * vTilde.y * R_over_C_v
				- X.ptr[2] * (R_over_C_v * vTildeL.y * vTilde.y - hTotal)
				- X.ptr[3] * vTilde.y * R_over_C_v * vTildeL.z 
				+ X.ptr[4] * vTilde.y * (R_over_C_v + 1.) 
				- X.ptr[5] * vTilde.y * (R_over_C_v - 2./3.),
			-X.ptr[0] * k * vTilde.y
				+ X.ptr[2] * k
				+ X.ptr[5] * vTilde.y,
			-X.ptr[0] * omega * vTilde.y
				+ X.ptr[2] * omega
				+ X.ptr[6] * vTilde.y,
			0.,
		}};
	} else if (n.side == 2) {
		<?=prefixes[2]?>
		real CsSq = Cs * Cs;
		real PStar = CsSq * rhoBar / (1. + R_over_C_v);
		return (<?=eqn.cons_t?>){.ptr={
			X.ptr[3],
			- X.ptr[0] * (vTilde.x * vTilde.z - .5 * gUj.x * R_over_C_v * vTildeSq)
				+ X.ptr[1] * (vTilde.z - gUj.x * R_over_C_v * vTildeL.x)
				- X.ptr[2] * gUj.x * R_over_C_v * vTildeL.y 
				+ X.ptr[3] * (vTilde.x - gUj.x * R_over_C_v * vTildeL.z) 
				+ X.ptr[4] * gUj.x * R_over_C_v
				- X.ptr[5] * gUj.x * (R_over_C_v - 2./3.),
			- X.ptr[0] * (vTilde.y * vTilde.z - .5 * gUj.y * R_over_C_v * vTildeSq)
				- X.ptr[1] * gUj.y * R_over_C_v * vTildeL.x 
				+ X.ptr[2] * (vTilde.z - gUj.y * R_over_C_v * vTildeL.y)
				+ X.ptr[3] * (vTilde.y - gUj.y * R_over_C_v * vTildeL.z)
				+ X.ptr[4] * gUj.y * R_over_C_v
				- X.ptr[5] * gUj.y * (R_over_C_v - 2./3.),
			- X.ptr[0] * (vTilde.z * vTilde.z - .5 * gUj.z * R_over_C_v * vTildeSq)
				- X.ptr[1] * gUj.z * R_over_C_v * vTildeL.x 
				- X.ptr[2] * gUj.z * R_over_C_v * vTildeL.y
				+ X.ptr[3] * (2. * vTilde.z - gUj.z * R_over_C_v * vTildeL.z)
				+ X.ptr[4] * gUj.z * R_over_C_v
				- X.ptr[5] * gUj.z * (R_over_C_v - 2./3.),
			+ X.ptr[0] * vTilde.z * (.5 * R_over_C_v * vTildeSq - hTotal)
				- X.ptr[1] * vTilde.z * vTildeL.x * R_over_C_v
				- X.ptr[2] * vTilde.z * vTildeL.y * R_over_C_v
				- X.ptr[3] * (R_over_C_v * vTildeL.z * vTilde.z - hTotal)
				+ X.ptr[4] * vTilde.z * (R_over_C_v + 1.) 
				- X.ptr[5] * vTilde.z * (R_over_C_v - 2./3.),
			- X.ptr[0] * k * vTilde.z
				+ X.ptr[5] * vTilde.z
				+ X.ptr[3] * k,
			- X.ptr[0] * omega * vTilde.z
				+ X.ptr[6] * vTilde.z
				+ X.ptr[3] * omega,
			0.,
		}};
	}
}


// used by PLM


<?=eqn.eigen_t?> eigen_forCell(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real vTildeSq = coordLenSq(W.vTilde, x);
	real eKinTilde = .5 * vTildeSq;
	real hTotal = (W.PStar + U.rhoBar_eTotalTilde) / U.rhoBar - U.ePot;
	real CsSq = R_over_C_v * (hTotal - eKinTilde) + 2./3. * W.k;
	real Cs = sqrt(CsSq);
	return (<?=eqn.eigen_t?>){
		.rhoBar = W.rhoBar,
		.vTilde = W.vTilde,
		.hTotal = hTotal,
		.k = W.k,
		.omega = W.omega,
		.vTildeSq = vTildeSq,
		.Cs = Cs,
	};
}

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

<? if not require 'hydro.coord.cartesian'.is(solver.coord) then ?>
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real3 m_conn_vv = coord_conn_apply23(W.vTilde, U->rhoBar_vTilde, x);
	deriv->rhoBar_vTilde = real3_sub(deriv->rhoBar_vTilde, m_conn_vv);	//-Conn^i_jk rhoBar vTilde^j vTilde^k 
	deriv->rhoBar_vTilde = real3_add(deriv->rhoBar_vTilde, real3_real_mul(coord_raise(coord_conn_trace13(x), x), W.PStar));		//+Conn^j_kj g^ki PStar
<? end ?>
}
