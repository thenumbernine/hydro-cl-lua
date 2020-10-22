<? if moduleName == nil then ?>
<? elseif moduleName == "primFromCons" then 
depmod{
	"prim_t",
	"cons_t",
	"coordLenSq",
}
?>

<?=eqn.prim_t?> primFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 vTilde = real3_real_mul(U.rhoBar_vTilde, 1. / U.rhoBar);
	real vTildeSq = coordLenSq(vTilde, x);
	real rhoBar_eIntTilde = U.rhoBar_eTotalTilde - .5 * U.rhoBar * vTildeSq - U.rhoBar_k;
	real rhoBar_TTilde = rhoBar_eIntTilde / solver->C_v;
	real PBar = rhoBar_TTilde * solver->gasConstant;
	real PStar = PBar + 2./3. * U.rhoBar_k;
	
	return (<?=eqn.prim_t?>){
		.rhoBar = U.rhoBar,
		.vTilde = vTilde,
		.PStar = PStar,
		.k = U.rhoBar_k / U.rhoBar,
		.omega = U.rhoBar_omega / U.rhoBar,
		.ePot = U.ePot,
	};
}

<? elseif moduleName == "consFromPrim" then 
depmod{
	"prim_t",
	"cons_t",
	"coordLenSq",
}
?>

<?=eqn.cons_t?> consFromPrim(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> W,
	real3 x
) {
	real rhoBar_k = W.rhoBar * W.k;

	//eqn 6: PStar = PBar + 2/3 rhoBar k
	real PBar = W.PStar - 2./3. * rhoBar_k;

	//eqn 10: PBar = rhoBar R TTilde
	real TTilde = PBar / (W.rhoBar * solver->gasConstant);
	
	//eqn 6: eIntTilde = C_v TTilde
	real eIntTilde = solver->C_v * TTilde;
	
	//eqn 6: eTotalTilde = eIntTilde + 1/2 vTilde^2 + W.k
	//so eTotalTilde = C_v PStar / rhoBar + 1/2 vTilde^2 + (1 - 2/3 C_v / solver->gasConstant) k
	real eTotalTilde = eIntTilde + .5 * coordLenSq(W.vTilde, x) + W.k;
	
	return (<?=eqn.cons_t?>){
		.rhoBar = W.rhoBar,
		.rhoBar_vTilde = real3_real_mul(W.vTilde, W.rhoBar),
		.rhoBar_eTotalTilde = W.rhoBar * eTotalTilde,
		.rhoBar_k = rhoBar_k,
		.rhoBar_omega = W.rhoBar * W.omega,
		.ePot = W.ePot,
	};
}

<? elseif moduleName == "apply_dU_dW" then 
depmod{
	"real3",
	"solver_t",
	"prim_t",
	"cons_t",
	"coord_lower",
}
?>

<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	real3 WA_vTildeL = coord_lower(WA.vTilde, x);
	return (<?=eqn.cons_t?>){
		.rhoBar = W.rhoBar,
		.rhoBar_vTilde = real3_add(
			real3_real_mul(WA.vTilde, W.rhoBar), 
			real3_real_mul(W.vTilde, WA.rhoBar)),
		.rhoBar_eTotalTilde = W.rhoBar * (.5 * real3_dot(WA.vTilde, WA_vTildeL) 
				+ (1. - 2./3. * C_v_over_R) * WA.k)
			+ WA.rhoBar * real3_dot(W.vTilde, WA_vTildeL)
			+ W.PStar * C_v_over_R
			+ (1. - 2./3 * C_v_over_R) * WA.rhoBar * W.k,
		.rhoBar_k = WA.k * W.rhoBar + WA.rhoBar * W.k,
		.rhoBar_omega = WA.omega * W.rhoBar + WA.rhoBar * W.omega,
		.ePot = W.ePot,
	};
}

<? elseif moduleName == "apply_dW_dU" then 
depmod{
	"real3",
	"solver_t",
	"prim_t",
	"cons_t",
	"coord_lower",
}
?>

<?=eqn.prim_t?> apply_dW_dU(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 WA_vTildeL = coord_lower(WA.vTilde, x);
	return (<?=eqn.prim_t?>){
		.rhoBar = U.rhoBar,
		.vTilde = real3_sub(
			real3_real_mul(U.rhoBar_vTilde, 1. / WA.rhoBar),
			real3_real_mul(WA.vTilde, U.rhoBar / WA.rhoBar)),
		.PStar = R_over_C_v * (
				.5 * real3_dot(WA.vTilde, WA_vTildeL) * U.rhoBar 
				- real3_dot(U.rhoBar_vTilde, WA_vTildeL)
				+ U.rhoBar_eTotalTilde
			) + (2./3. * R_over_C_v - 1.) * U.rhoBar_k,
		.k = U.rhoBar_k / WA.rhoBar - WA.k / WA.rhoBar * U.rhoBar,
		.omega = U.rhoBar_omega / WA.rhoBar - WA.omega / WA.rhoBar * U.rhoBar,
		.ePot = U.ePot,
	};
}

<? elseif moduleName == "eqn.common" then ?>

#define R_over_C_v (solver->gasConstant / solver->C_v)
#define C_v_over_R (solver->C_v / solver->gasConstant)

//real calc_H(real PStar) { return PStar * ((R_over_C_v + 1.) / (R_over_C_v)); }
//real calc_h(real rhoBar, real PStar) { return calc_H(PStar) / rhoBar; }
//real calc_hTotal(real rhoBar, real PStar, real rhoBar_eTotalTilde) { return (PStar + rhoBar_eTotalTilde) / rhoBar; }
//real calc_HTotal(real PStar, real rhoBar_eTotalTilde) { return PStar + rhoBar_eTotalTilde; }
real calc_eKinTilde(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.vTilde, x); }
real calc_EKinTilde(<?=eqn.prim_t?> W, real3 x) { return W.rhoBar * calc_eKinTilde(W, x); }

//before
//real calc_EIntTilde(<?=eqn.prim_t?> W) { return W.PStar * C_v_over_R; }
//real calc_eIntTilde(<?=eqn.prim_t?> W) { return calc_EIntTilde(W) / W.rhoBar; }

//after
real calc_PBar(<?=eqn.prim_t?> W) { return W.PStar - 2./3. * W.rhoBar * W.k; }
real calc_TTilde(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_PBar(W) / (W.rhoBar * solver->gasConstant); }
real calc_eIntTilde(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return solver->C_v * calc_TTilde(solver, W); }
real calc_EIntTilde(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.rhoBar * calc_eIntTilde(solver, W); }

real calc_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.rhoBar_vTilde, x) / U.rhoBar; }
real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return calc_EKinTilde(W, x) + calc_EIntTilde(solver, W);
}

real calc_Cs(constant <?=solver.solver_t?>* solver, const <?=eqn.prim_t?> W) {
	return sqrt((R_over_C_v + 1.) * W.PStar / W.rhoBar);
}

<? elseif moduleName == "applyInitCond" then 
depmod{
	"cartesianToCoord",
	"consFromPrim",
}
?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = true
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>	&& x.<?=xi?> < mids.<?=xi?>
<?
end
?>;
	
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	
	//TODO make this B for Maxwell
	
	real3 B = real3_zero;	//set for MHD / thrown away for pure NavierStokesWilcox
	real ePot = 0;

	<?=initCode()?>

	<?=eqn.prim_t?> W = {
		.rhoBar = rho,
		.vTilde = cartesianToCoord(v, x),	//transform from cartesian to coordinate space 
		.PStar = P,
		.k = 0,
		.omega = 0,
		.ePot = ePot,
	};
	UBuf[index] = consFromPrim(solver, W, x);
}

<? elseif moduleName == "fluxFromCons" then 
depmod{
	"coord_g_uu##",
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
	real vTilde_j = W.vTilde.s[n.side];
	
	//this is the flux term used, but is it technically called 'HTotal' ?
	//'HTotal' = rhoBar (1/2 vTilde^2 + (1 - 2/3/gamma_1) k) + 1/gamma_1 PStar ... + PStar
	//'HTotal' = rhoBar (1/2 vTilde^2 + (1 - 2/3 C_v/R) k) + (1 + C_v/R) PStar
	//'hTotal' = 1/2 vTilde^2 + (1 - 2/3 C_v/R) k + (1 + C_v/R) PStar / rhoBar
	real HTotal = U.rhoBar_eTotalTilde + W.PStar;
	
	<?=eqn.cons_t?> F;
	
	F.rhoBar = U.rhoBar_vTilde.s[n.side];
	
	F.rhoBar_vTilde = real3_real_mul(U.rhoBar_vTilde, vTilde_j);

	if (false) {}
<? for side=0,2 do ?>
	else if (n.side == <?=side?>) {
<? for i=0,2 do
?>	F.rhoBar_vTilde.s<?=i?> += coord_g_uu<?=i?><?=side?>(x) * W.PStar;
<? end
?>	}
<? end ?>

	F.rhoBar_eTotalTilde = HTotal * vTilde_j;
	
	F.rhoBar_k = 0;
	F.rhoBar_omega = 0;

	F.ePot = 0;
	
	return F;
}

<? elseif moduleName == "calcCellMinMaxEigenvalues" then 
depmod{
	"coord_sqrt_g_uu##",
	"primFromCons",
	"eqn.common",
}
?>

range_t calcCellMinMaxEigenvalues(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real Cs = calc_Cs(solver, W);
	if (false) {}
<? for side=0,2 do ?>
	else if (n.side == <?=side?>) {
		real Cs_sqrt_gU = Cs * coord_sqrt_g_uu<?=side..side?>(x);
		return (range_t){
			.min = W.vTilde.s<?=side?> - Cs_sqrt_gU, 
			.max = W.vTilde.s<?=side?> + Cs_sqrt_gU,
		};
	}
<? end ?>
}

<? elseif moduleName == "eigen_forInterface" then 
depmod{
	"primFromCons",
}
?>

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

<? elseif moduleName == "eigen_left/rightTransform" then 
depmod{
	"coord_g_uu",
	"coord_g_uu##",
	"coord_sqrt_g_uu##",
	"coord_lower",
}
?>

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

<? elseif moduleName == "eigen_fluxTransform" then ?>

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

<? elseif moduleName == "eigen_forCell" then 
depmod{
	"primFromCons",
}
?>

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

<? elseif moduleName == "addSource" then 
depmod{
	"cell_x",
	"primFromCons",
}
?>

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

<? 
else
	error("unknown moduleName "..require 'ext.tolua'(moduleName))
end 
?>
