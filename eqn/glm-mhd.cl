/*
started with the mhd.lua and mhd.cl
tweaked it while looking at
2006 Xueshang et al - A 3rd Order WENO GLM-MHD Scheme for Magnetic Reconnection
2009 Mignone, Tzeferacos - A Second-Order Unsplit Godunov Scheme for Cell-Centered MHD- the CTU-GLM scheme
*/

<?
local common = require 'common'	-- xNames, symNames
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym
?>

//assumes UL and UR are already rotated so the 'x' direction is our flux direction
Roe_t calcRoeValues(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> UL, 
	<?=eqn.cons_t?> UR,
	real3 x
) {
	Roe_t W;

	// should I use Bx, or BxL/R, for calculating the PMag at the L and R states?
	<?=eqn.prim_t?> WL = primFromCons(solver, UL, x);
	real sqrtRhoL = sqrt(UL.rho);
	real PMagL = .5 * real3_lenSq(UL.B);
	real hTotalL = (UL.ETotal + WL.P + PMagL) / UL.rho;

	<?=eqn.prim_t?> WR = primFromCons(solver, UR, x);
	real sqrtRhoR = sqrt(UR.rho);
	real PMagR = .5 * real3_lenSq(UR.B);
	real hTotalR = (UR.ETotal + WR.P + PMagR) / UR.rho;
	
	real dby = WL.B.y - WR.B.y;
	real dbz = WL.B.z - WR.B.z;
	
	real invDenom = 1 / (sqrtRhoL + sqrtRhoR);
	
	W.rho  = sqrtRhoL * sqrtRhoR;
	W.v = real3_real_mul(real3_add(
		real3_real_mul(WL.v, sqrtRhoL),
		real3_real_mul(WR.v, sqrtRhoR)), invDenom);
	
	W.hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
	
	W.B.x = (sqrtRhoL * WL.B.x + sqrtRhoR * WR.B.x) * invDenom;
	// why does athena switch the weights of the By and Bz components?
	W.B.y = (sqrtRhoR * WL.B.y + sqrtRhoL * WR.B.y) * invDenom;
	W.B.z = (sqrtRhoR * WL.B.z + sqrtRhoL * WR.B.z) * invDenom;
	
	W.X = .5 * (dby * dby + dbz * dbz) * invDenom * invDenom;
	W.Y = .5 * (UL.rho + UR.rho) / W.rho;
	
	return W;
};

//assumes the vector values are x-axis aligned with the interface normal
<?=eqn.eigen_t?> eigen_forRoeAvgs(
	constant <?=solver.solver_t?>* solver,
	Roe_t roe,
	real3 x
) {
	<?=eqn.eigen_t?> eig;
	
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	const real gamma_3 = gamma - 3.;
	
	real rho = roe.rho;
	real3 v = roe.v;
	real hTotal = roe.hTotal;
	real3 B = roe.B;
	real X = roe.X;
	real Y = roe.Y;

	real _1_rho = 1. / rho;
	real vSq = real3_lenSq(v);
	real BPerpSq = B.y*B.y + B.z*B.z;
	real BStarPerpSq = (gamma_1 - gamma_2 * Y) * BPerpSq;
	real CAxSq = B.x*B.x*_1_rho;
	real CASq = CAxSq + BPerpSq * _1_rho;
	eig.hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// eig.hHydro = hTotal - CASq, CASq = EMag/rho
	// eig.hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// eig.hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(eig.hHydro - eKin) = gamma P / rho
	eig.aTildeSq = max((gamma_1 * (eig.hHydro - .5 * vSq) - gamma_2 * X), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + eig.aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - eig.aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + eig.aTildeSq * BStarPerpSq_rho);
	
	eig.CAx = sqrt(CAxSq);
	
	real CfSq = CStarSq + sqrtDiscr;
	eig.Cf = sqrt(CfSq);
	
	real CsSq = eig.aTildeSq * CAxSq / CfSq;
	eig.Cs = sqrt(CsSq);

<? if not eqn.useFixedCh then ?>
	eig.Ch = max(max(fabs(v.x) + eig.Cf, fabs(v.y) + eig.Cf), fabs(v.z) + eig.Cf);
<? end ?>

	real BPerpLen = sqrt(BPerpSq);
	eig.BStarPerpLen = sqrt(BStarPerpSq);
	
	if (BPerpLen == 0) {
		eig.betaY = 1;
		eig.betaZ = 0;
	} else {
		eig.betaY = B.y / BPerpLen;
		eig.betaZ = B.z / BPerpLen;
	}
	eig.betaStarY = eig.betaY / sqrt(gamma_1 - gamma_2*Y);
	eig.betaStarZ = eig.betaZ / sqrt(gamma_1 - gamma_2*Y);
	eig.betaStarSq = eig.betaStarY*eig.betaStarY + eig.betaStarZ*eig.betaStarZ;


	if (CfSq - CsSq == 0) {
		eig.alphaF = 1;
		eig.alphaS = 0;
	} else if (eig.aTildeSq - CsSq <= 0) {
		eig.alphaF = 0;
		eig.alphaS = 1;
	} else if (CfSq - eig.aTildeSq <= 0) {
		eig.alphaF = 1;
		eig.alphaS = 0;
	} else {
		eig.alphaF = sqrt((eig.aTildeSq - CsSq) / (CfSq - CsSq));
		eig.alphaS = sqrt((CfSq - eig.aTildeSq) / (CfSq - CsSq));
	}


	eig.sqrtRho = sqrt(rho);
	real _1_sqrtRho = 1. / eig.sqrtRho;
	eig.sbx = B.x >= 0 ? 1 : -1;
	real aTilde = sqrt(eig.aTildeSq);
	eig.Qf = eig.Cf * eig.alphaF * eig.sbx;
	eig.Qs = eig.Cs * eig.alphaS * eig.sbx;
	eig.Af = aTilde * eig.alphaF * _1_sqrtRho;
	eig.As = aTilde * eig.alphaS * _1_sqrtRho;


	//used for eigenvectors and eigenvalues
<? 	for _,var in ipairs(eqn.roeVars) do
		local name = var.name
?>	eig.<?=name?> = roe.<?=name?>;
<?	end
?>

	return eig;
}

<?=eqn.cons_t?> cons_alignFrom(<?=eqn.cons_t?> U, real3 n) {
	//U.m = real3_swap<?=side?>(U.m);
	//U.B = real3_swap<?=side?>(U.B);
	U.m = real3_rotateFrom(U.m, n);
	U.B = real3_rotateFrom(U.B, n);
	return U;
}

<?=eqn.cons_t?> cons_alignTo(<?=eqn.cons_t?> U, real3 n) {
	//U.m = real3_swap<?=side?>(U.m);
	//U.B = real3_swap<?=side?>(U.B);
	U.m = real3_rotateTo(U.m, n);
	U.B = real3_rotateTo(U.B, n);
	return U;
}

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> cons_swapFrom<?=side?>(<?=eqn.cons_t?> U) {
	return cons_alignFrom(U, normalForSide<?=side?>());
}
<?=eqn.cons_t?> cons_swapTo<?=side?>(<?=eqn.cons_t?> U) {
	return cons_alignTo(U, normalForSide<?=side?>());
}
<? end ?>

<?=eqn.eigen_t?> eigen_forInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	real3 n
) {
	//swap the sides with x here, so all the fluxes are in the 'x' direction
	<?=eqn.cons_t?> UL_ = cons_alignFrom(UL, n);
	<?=eqn.cons_t?> UR_ = cons_alignFrom(UR, n);
	Roe_t roe = calcRoeValues(solver, UL_, UR_, x);
	return eigen_forRoeAvgs(solver, roe, x);
}

<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real PMag = .5 * real3_lenSq(W.B);
	real hTotal = (U.ETotal + W.P + PMag) / W.rho;
	Roe_t roe = (Roe_t){
		.rho = W.rho,
		.v = W.v,
		.hTotal = hTotal,
		.B = W.B,
		.X = 0,
		.Y = 1,
	};
	return eigen_forRoeAvgs(solver, roe, x);
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real vj = W.v.s<?=side?>;
	real Bj = W.B.s<?=side?>;
	real BSq = real3_lenSq(W.B);
	real BDotV = real3_dot(W.B, W.v);
	real PMag = .5 * BSq;
	real PTotal = W.P + PMag;
	real HTotal = U.ETotal + PTotal;
	
	<?=eqn.cons_t?> F;
	F.rho = U.m.s<?=side?>;
	F.m = real3_sub(real3_real_mul(U.m, vj), real3_real_mul(U.B, Bj / solver->mu0));
	F.m.s<?=side?> += PTotal;
	F.B = real3_sub(real3_real_mul(U.B, vj), real3_real_mul(W.v, Bj));
	F.B.s<?=side?> += F.psi;
	F.ETotal = HTotal * vj - BDotV * Bj / solver->mu0;

<? if not eqn.useFixedCh then ?>
	//TODO don't need the whole eigen here, just the Ch
	real Ch = 0;
	<? for side=0,solver.dim-1 do ?>{
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(solver, U, x);
		Ch = max(Ch, eig.Ch);
	}<? end ?>
<? end ?>

	F.psi = Ch * Ch;
	return F;
}
<? end ?>

//called from calcDT
<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.cons_t?> U_ = cons_swapFrom<?=side?>(*U);
	<?=eqn.prim_t?> W = primFromCons(solver, U_, x);
	
#if 0
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real3 v = W.v;
	real3 B = W.B;
	
	real BSq = real3_lenSq(B);
	real invRho = 1./W.rho;
	
	real aSq = solver->heatCapacityRatio * W.P * invRho;
	real CaxSq = B.s<?=side?> * B.s<?=side?> * invRho;
	real CaSq = BSq * invRho;
	
	real CStarSq = .5 * (CaSq + aSq);
	real sqrtCfsDiscr = sqrt(max(0., CStarSq * CStarSq - aSq * CaxSq));
	
	real CfSq = CStarSq + sqrtCfsDiscr;
	real CsSq = CStarSq - sqrtCfsDiscr;

	real Cf = sqrt(CfSq);
	real Cs = sqrt(max(CsSq, 0.));
	return (range_t){.min=v.s<?=side?> - Cf, .max=v.s<?=side?> + Cf};
#else
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	
	real rho = W.rho;
	real3 v = W.v;
	real3 B = W.B;
	real hTotal = .5 * real3_lenSq(W.v) + (W.P * gamma / gamma_1 + real3_lenSq(B)) / W.rho;

	//the rest of this matches calcEigenBasis:

	real _1_rho = 1. / rho;
	real vSq = real3_lenSq(v);
	real BPerpSq = B.y*B.y + B.z*B.z;
	real BStarPerpSq = (gamma_1 - gamma_2) * BPerpSq;
	real CAxSq = B.x*B.x*_1_rho;
	real CASq = CAxSq + BPerpSq * _1_rho;
	real hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// hHydro = hTotal - CASq, CASq = EMag/rho
	// hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * BStarPerpSq_rho);
	
	real CfSq = CStarSq + sqrtDiscr;
	real Cf = sqrt(CfSq);

	real CsSq = aTildeSq * CAxSq / CfSq;
	real Cs = sqrt(CsSq);

<? if not eqn.useFixedCh then ?>
	real Ch = max(max(fabs(v.x) + Cf, fabs(v.y) + Cf), fabs(v.z) + Cf);
<? end ?>

	real lambdaFastMin = v.x - Cf;
	real lambdaFastMax = v.x + Cf;
	
	return (range_t){
		.min = min(lambdaFastMin, -Ch),
		.max = max(lambdaFastMax, Ch),
	};
#endif
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> inputU,
	real3 x
) {	
	inputU = cons_swapFrom<?=side?>(inputU);
	
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;

<? for _,var in ipairs(eqn.eigenVars) do
	local name = var.name
	local ctype = var.type
?> 	<?=ctype?> <?=name?> = eig.<?=name?>;
<? end ?>

	real vSq = real3_lenSq(v);
	
	// left eigenvectors
	real norm = .5 / eig.aTildeSq;
	real Cff = norm * alphaF * Cf;
	real Css = norm * alphaS * Cs;
	Qf = Qf * norm;
	Qs = Qs * norm;
	real AHatF = norm * Af * rho;
	real AHatS = norm * As * rho;
	real afpb = norm * Af * BStarPerpLen;
	real aspb = norm * As * BStarPerpLen;

	norm = norm * gamma_1;
	alphaF = alphaF * norm;
	alphaS = alphaS * norm;
	real QStarY = betaStarY/betaStarSq;
	real QStarZ = betaStarZ/betaStarSq;
	real vqstr = (v.y * QStarY + v.z * QStarZ);
	norm = norm * 2.;

	
	real l16 = AHatS * QStarY - alphaF * B.y;
	real l17 = AHatS * QStarZ - alphaF * B.z;
	real l21 = .5 * (v.y * betaZ - v.z * betaY);
	real l23 = .5 * betaZ;
	real l24 = .5 * betaY;
	real l26 = -.5 * sqrtRho * betaZ * sbx;
	real l27 = .5 * sqrtRho * betaY * sbx;
	real l36 = -AHatF * QStarY - alphaS * B.y;
	real l37 = -AHatF * QStarZ - alphaS * B.z;
	

	<?=eqn.waves_t?> result;
	result.ptr[0] = 
		  inputU.rho * (alphaF * (vSq - eig.hHydro) + Cff * (Cf + v.x) - Qs * vqstr - aspb) 
		+ inputU.m.x * (-alphaF * v.x - Cff)
		+ inputU.m.y * (-alphaF * v.y + Qs * QStarY)
		+ inputU.m.z * (-alphaF * v.z + Qs * QStarZ)
		+ inputU.ETotal * alphaF
		+ inputU.B.y * l16
		+ inputU.B.z * l17;
	result.ptr[1] = 
		  inputU.rho * l21
		+ inputU.m.y * l23
		+ inputU.m.z * l24
		+ inputU.B.y * l26
		+ inputU.B.z * l27;
	result.ptr[2] = 
		  inputU.rho * (alphaS * (vSq - eig.hHydro) + Css * (Cs + v.x) + Qf * vqstr + afpb)
		+ inputU.m.x * (-alphaS * v.x - Css)
		+ inputU.m.y * (-alphaS * v.y - Qf * QStarY)
		+ inputU.m.z * (-alphaS * v.z - Qf * QStarZ)
		+ inputU.ETotal * alphaS
		+ inputU.B.y * l36
		+ inputU.B.z * l37;
	result.ptr[3] = 
		  inputU.rho * (1. - norm * (.5 * vSq - gamma_2 * X / gamma_1))
		+ inputU.m.x * norm*v.x
		+ inputU.m.y * norm*v.y
		+ inputU.m.z * norm*v.z
		+ inputU.ETotal * -norm
		+ inputU.B.y * norm*B.y
		+ inputU.B.z * norm*B.z;
	result.ptr[4] = 
		  inputU.rho * (alphaS * (vSq - eig.hHydro) + Css * (Cs - v.x) - Qf * vqstr + afpb)
		+ inputU.m.x * (-alphaS * v.x + Css)
		+ inputU.m.y * (-alphaS * v.y + Qf * QStarY)
		+ inputU.m.z * (-alphaS * v.z + Qf * QStarZ)
		+ inputU.ETotal * alphaS
		+ inputU.B.y * l36
		+ inputU.B.z * l37;
	result.ptr[5] = 
		  inputU.rho * -l21
		+ inputU.m.y * -l23
		+ inputU.m.z * -l24
		+ inputU.B.y * l26
		+ inputU.B.z * l27;
	result.ptr[6] = 
		  inputU.rho * (alphaF * (vSq - eig.hHydro) + Cff * (Cf - v.x) + Qs * vqstr - aspb)
		+ inputU.m.x * (-alphaF * v.x + Cff)
		+ inputU.m.y * (-alphaF * v.y - Qs * QStarY)
		+ inputU.m.z * (-alphaF * v.z - Qs * QStarZ)
		+ inputU.ETotal * alphaF
		+ inputU.B.y * l16
		+ inputU.B.z * l17;
	result.ptr[7] = inputU.B.x * .5;
	result.ptr[8] = inputU.B.x * .5;
	if (Ch != 0) {
		result.ptr[7] += inputU.psi * -.5 / Ch;
		result.ptr[8] += inputU.psi * .5 / Ch;
	}
	return result;
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> input,
	real3 x
) {
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;

<? for _,var in ipairs(eqn.eigenVars) do
	local name = var.name
	local ctype = var.type
?> 	<?=ctype?> <?=name?> = eig.<?=name?>;
<? end ?>

	real vSq = real3_lenSq(v);
	real vDotBeta = v.y*betaStarY + v.z*betaStarZ;
	real _1_sqrtRho = 1. / eig.sqrtRho;
	real Afpbb = Af * BStarPerpLen * betaStarSq;
	real Aspbb = As * BStarPerpLen * betaStarSq;

	real lambdaFastMin = eig.v.x - eig.Cf;
	real lambdaSlowMin = eig.v.x - eig.Cs;
	real lambdaSlowMax = eig.v.x + eig.Cs;
	real lambdaFastMax = eig.v.x + eig.Cf;

	// right eigenvectors
	real qa3 = alphaF * v.y;
	real qb3 = alphaS * v.y;
	real qc3 = Qs * betaStarY;
	real qd3 = Qf * betaStarY;
	real qa4 = alphaF * v.z;
	real qb4 = alphaS * v.z;
	real qc4 = Qs * betaStarZ;
	real qd4 = Qf * betaStarZ;
	real r52 = -(v.y * betaZ - v.z * betaY);
	real r61 = As * betaStarY;
	real r62 = -betaZ * sbx * _1_sqrtRho;
	real r63 = -Af * betaStarY;
	real r71 = As * betaStarZ;
	real r72 = betaY * sbx * _1_sqrtRho;
	real r73 = -Af * betaStarZ;

	<?=eqn.cons_t?> resultU;
	resultU.rho =
		  input.ptr[0] * alphaF
		+ input.ptr[2] * alphaS
		+ input.ptr[3]
		+ input.ptr[4] * alphaS
		+ input.ptr[6] * alphaF;
	resultU.m.x =
		  input.ptr[0] * alphaF * lambdaFastMin
		+ input.ptr[2] * alphaS * lambdaSlowMin
		+ input.ptr[3] * v.x
		+ input.ptr[4] * alphaS * lambdaSlowMax
		+ input.ptr[6] * alphaF * lambdaFastMax;
	resultU.m.y =
		  input.ptr[0] * (qa3 + qc3)
		+ input.ptr[1] * -betaZ
		+ input.ptr[2] * (qb3 - qd3)
		+ input.ptr[3] * v.y
		+ input.ptr[4] * (qb3 + qd3)
		+ input.ptr[5] * betaZ
		+ input.ptr[6] * (qa3 - qc3);
	resultU.m.z =
		  input.ptr[0] * (qa4 + qc4)
		+ input.ptr[1] * betaY
		+ input.ptr[2] * (qb4 - qd4)
		+ input.ptr[3] * v.z
		+ input.ptr[4] * (qb4 + qd4)
		+ input.ptr[5] * -betaY
		+ input.ptr[6] * (qa4 - qc4);
	resultU.ETotal =
		  input.ptr[0] * (alphaF*(eig.hHydro - v.x*Cf) + Qs*vDotBeta + Aspbb)
		+ input.ptr[1] * r52
		+ input.ptr[2] * (alphaS*(eig.hHydro - v.x*Cs) - Qf*vDotBeta - Afpbb)
		+ input.ptr[3] * (.5*vSq + gamma_2*X/gamma_1)
		+ input.ptr[4] * (alphaS*(eig.hHydro + v.x*Cs) + Qf*vDotBeta - Afpbb)
		+ input.ptr[5] * -r52
		+ input.ptr[6] * (alphaF*(eig.hHydro + v.x*Cf) - Qs*vDotBeta + Aspbb);
	resultU.B.x = 
		  input.ptr[7] * .5
		+ input.ptr[8] * .5;
	resultU.B.y =
		  input.ptr[0] * r61
		+ input.ptr[1] * r62
		+ input.ptr[2] * r63
		+ input.ptr[4] * r63
		+ input.ptr[5] * r62
		+ input.ptr[6] * r61;
	resultU.B.z =
		  input.ptr[0] * r71
		+ input.ptr[1] * r72
		+ input.ptr[2] * r73
		+ input.ptr[4] * r73
		+ input.ptr[5] * r72
		+ input.ptr[6] * r71;
	resultU.psi = 
		  input.ptr[7] * -Ch
		+ input.ptr[8] * Ch;

	real3 n = real3_zero; n.s<?=side?> = 1.;
	return cons_alignTo(resultU, n);
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> inputU,
	real3 x
) {
	inputU = cons_swapFrom<?=side?>(inputU);

	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	const real gamma_3 = gamma - 3.;

<? for _,var in ipairs(eqn.eigenVars) do
	local name = var.name
	local ctype = var.type
?> 	<?=ctype?> <?=name?> = eig.<?=name?>;
<? end ?>

	real _1_rho = 1. / rho;
	real vSq = real3_lenSq(v);
	real BDotV = real3_dot(B,v);

	// dF/dU
	<?=eqn.cons_t?> resultU;
	resultU.rho = inputU.m.x;
	resultU.m.x =
		  inputU.rho * (-v.x*v.x + .5*gamma_1*vSq - gamma_2*X)
		+ inputU.m.x * -gamma_3*v.x
		+ inputU.m.y * -gamma_1*v.y
		+ inputU.m.z * -gamma_1*v.z
		+ inputU.ETotal * gamma_1
		+ inputU.B.y * -gamma_2*Y*B.y
		+ inputU.B.z * -gamma_2*Y*B.z;
	resultU.m.y =
		  inputU.rho * -v.x*v.y
		+ inputU.m.x * v.y
		+ inputU.m.y * v.x
		+ inputU.B.y * -B.x;
	resultU.m.z =
		  inputU.rho * -v.x*v.z
		+ inputU.m.x * v.z
		+ inputU.m.z * v.x
		+ inputU.B.z * -B.x;
	resultU.ETotal =
		  inputU.rho * (v.x*(.5*gamma_1*vSq - hTotal) + B.x*BDotV * _1_rho)
		+ inputU.m.x * (-gamma_1*v.x*v.x + hTotal - B.x*B.x * _1_rho)
		+ inputU.m.y * (-gamma_1*v.x*v.y - B.x*B.y * _1_rho)
		+ inputU.m.z * (-gamma_1*v.x*v.z - B.x*B.z * _1_rho)
		+ inputU.ETotal * gamma*v.x
		+ inputU.B.y * (-gamma_2*Y*B.y*v.x - B.x*v.y)
		+ inputU.B.z * (-gamma_2*Y*B.z*v.x - B.x*v.z);
	resultU.B.x = inputU.B.x;
	resultU.B.y =
		  inputU.rho * (B.x*v.y - B.y*v.x) * _1_rho
		+ inputU.m.x * B.y * _1_rho
		+ inputU.m.y * -B.x * _1_rho
		+ inputU.B.y * v.x;
	resultU.B.z =
		  inputU.rho * (B.x*v.z - B.z*v.x) * _1_rho
		+ inputU.m.x * B.z * _1_rho
		+ inputU.m.z * -B.x * _1_rho
		+ inputU.B.z * v.x;
	resultU.psi = inputU.psi;
	
	real3 n = real3_zero; n.s<?=side?> = 1.;
	return cons_alignTo(resultU, n);
}
<? end ?>

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

#if 0
	//TODO don't need the whole eigen here, just the Ch	
<? if not eqn.useFixedCh then ?>
	real Ch = 0;
	<? for side=0,solver.dim-1 do ?>{
		<?=eqn.eigen_t?> eig = eigen_forCell_<?=side?>(solver, *U, x);
		Ch = max(Ch, eig.Ch);
	}<? end ?>
<? end ?>
#endif

	real divB = 0.<? 
for i,xi in ipairs(xNames:sub(1,solver.dim)) do
?> + (U[solver->stepsize.<?=xi?>].B.<?=xi?> - U[-solver->stepsize.<?=xi?>].B.<?=xi?>) / (2. * solver->grid_dx.s<?=i?>) <?
end ?>;

	real3 grad_psi = (real3){
<? for i,xi in ipairs(xNames:sub(1,solver.dim)) do
?>		.<?=xi?> = (U[solver->stepsize.<?=xi?>].psi - U[-solver->stepsize.<?=xi?>].psi) / (2. * solver->grid_dx.s<?=i?>),
<? end 
for i=solver.dim+1,3 do
?>		.<?=xNames[i]?> = 0.,
<? end
?>	};

#if 0	//2009 Mignone eqn 3
	deriv->psi -= Ch * U->psi;
#endif
#if 1	//2002 Dedner eqn 24 -- EGLM
<? for i,xi in ipairs(xNames) do
?>	deriv->m.<?=xi?> -= divB * U->B.<?=xi?>;
<? end
?>	deriv->ETotal -= real3_dot(U->B, grad_psi);
//	deriv->psi -= Ch * Ch / (Cp * Cp) * U->psi;
#endif
#if 0	//2002 Dedner eqn 38
<? 
for i,xi in ipairs(xNames) do
?>	deriv->m.<?=xi?> -= divB * U->B.<?=xi?>;
<? 
end
for i,xi in ipairs(xNames) do
?>	deriv->B.<?=xi?> -= divB * U->m.<?=xi?> / U->rho;
<?
end
?>	deriv->ETotal -= real3_dot(U->B, grad_psi) + real3_dot(U->m, U->B) * (divB / U->rho);
	
	//update this in the second step:
	//psi[n+1] = psi[n] * exp(-dt * ch * ch / (cp * cp))
	//deriv->psi -= Ch * Ch / (Cp * Cp) * U->psi;
	deriv->psi -= real3_dot(U->m, grad_psi) / U->rho;
#endif
}

kernel void constrainU(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);

	W.rho = max(W.rho, 1e-7);
	W.P = max(W.P, 1e-7);

	*U = consFromPrim(solver, W, x);
}

//this is a temporary fix until I implement MHD's inline eigenvalue code

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);

	const global <?=eqn.cons_t?>* U = UBuf + index;

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		//use cell-centered eigenvalues
		range_t lambda = calcCellMinMaxEigenvalues_<?=side?>(solver, U, x); 
		real absLambdaMax = max(fabs(lambda.min), fabs(lambda.max));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt;
}
