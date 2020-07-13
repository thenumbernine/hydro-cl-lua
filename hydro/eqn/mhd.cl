/*
Stone et al 2008 - https://arxiv.org/pdf/0804.0402v1.pdf
based on Athena's version of eigenvectors of derivative of adiabatic MHD flux wrt primitives
ideal-mhd, divergence-free, conservative-based eigensystem
*/

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normalInfo_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real vj = normalInfo_vecDotN1(n, W.v);
	real Bj = normalInfo_vecDotN1(n, W.B);
	real BSq = coordLenSq(W.B, x);
	real BDotV = real3_dot(W.B, W.v);
	real PMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);
	real PTotal = W.P + PMag;
	real HTotal = U.ETotal + PTotal;
	
	<?=eqn.cons_t?> F;
	F.rho = normalInfo_vecDotN1(n, U.m);
	F.m = real3_sub(real3_real_mul(U.m, vj), real3_real_mul(U.B, Bj / (solver->mu0 / unit_kg_m_per_C2)));
	F.m.x += PTotal * normalInfo_u1x(n);
	F.m.y += PTotal * normalInfo_u1y(n);
	F.m.z += PTotal * normalInfo_u1z(n);
	F.B = real3_sub(real3_real_mul(U.B, vj), real3_real_mul(W.v, Bj));
	F.ETotal = HTotal * vj - BDotV * Bj / (solver->mu0 / unit_kg_m_per_C2);
	F.psi = 0;
	F.ePot = 0;
	return F;
}

//align from vector coordinates to the normal basis
<?=eqn.cons_t?> cons_rotateFrom(<?=eqn.cons_t?> U, normalInfo_t n) {
	U.m = normalInfo_vecDotNs(n, U.m);
	U.B = normalInfo_vecDotNs(n, U.B);
	return U;
}

//align from normal basis to vector coordinates
<?=eqn.cons_t?> cons_rotateTo(<?=eqn.cons_t?> U, normalInfo_t n) {
	U.m = normalInfo_vecFromNs(n, U.m);
	U.B = normalInfo_vecFromNs(n, U.B);
	return U;
}

// TODO find out where mu_0 goes in the code below

//called from calcDT
range_t calcCellMinMaxEigenvalues(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normalInfo_t n
) {
	<?=eqn.cons_t?> U_ = cons_rotateFrom(U, n);
	<?=eqn.prim_t?> W = primFromCons(solver, U_, x);
	
#if 0
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real3 v = W.v;
	real3 B = W.B;
	
	real BSq = coordLenSq(B, x);
	real invRho = 1./W.rho;
	
	real aSq = solver->heatCapacityRatio * W.P * invRho;
	real B_n = normalInfo_vecDotN1(n, B);
	real CaxSq = B_n * B_n * invRho;
	real CaSq = BSq * invRho;
	
	real CStarSq = .5 * (CaSq + aSq);
	real sqrtCfsDiscr = sqrt(max(0., CStarSq * CStarSq - aSq * CaxSq));
	
	real CfSq = CStarSq + sqrtCfsDiscr;
	real CsSq = CStarSq - sqrtCfsDiscr;

	real Cf = sqrt(CfSq);
	real Cs = sqrt(max(CsSq, 0.));
	real v_n = normalInfo_vecDotN1(n, v);
	return (range_t){.min=v_n - Cf, .max=v_n + Cf};
#else
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	
	real rho = W.rho;
	real3 v = W.v;
	real3 B = W.B;
	real hTotal = .5 * coordLenSq(W.v, x) + (W.P * gamma / gamma_1 + coordLenSq(B, x)) / W.rho;

	//the rest of this matches calcEigenBasis:

	real _1_rho = 1. / rho;
	real vSq = coordLenSq(v, x);
#warning consider g_ij	
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

	real lambdaFastMin = v.x - Cf;
	real lambdaFastMax = v.x + Cf;
	
	return (range_t){
		.min = lambdaFastMin,
		.max = lambdaFastMax,
	};
#endif
}

//assumes UL and UR are already rotated so the 'x' direction is our flux direction
<?=eqn.roe_t?> calcRoeValues(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> UL, 
	<?=eqn.cons_t?> UR,
	real3 x
) {
	<?=eqn.roe_t?> W;
	
	// should I use Bx, or BxL/R, for calculating the PMag at the L and R states?
	<?=eqn.prim_t?> WL = primFromCons(solver, UL, x);
	real sqrtRhoL = sqrt(UL.rho);
	real PMagL = .5 * coordLenSq(UL.B, x);
	real hTotalL = (UL.ETotal + WL.P + PMagL) / UL.rho;

	<?=eqn.prim_t?> WR = primFromCons(solver, UR, x);
	real sqrtRhoR = sqrt(UR.rho);
	real PMagR = .5 * coordLenSq(UR.B, x);
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
	<?=eqn.roe_t?> roe,
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
	real vSq = coordLenSq(v, x);
#warning consider g_ij	
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

	
#warning consider g_ij	
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
?>	eig.<?=var.name?> = roe.<?=var.name?>;
<?	end
?>

	return eig;
}

<?=eqn.eigen_t?> eigen_forInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	normalInfo_t n
) {
	//rotate UL and UR to be x-aligned?  that takes the normal ...

	//swap the sides with x here, so all the fluxes are in the 'x' direction
	<?=eqn.cons_t?> UL_ = cons_rotateFrom(UL, n);
	<?=eqn.cons_t?> UR_ = cons_rotateFrom(UR, n);
	<?=eqn.roe_t?> roe = calcRoeValues(solver, UL_, UR_, x);
	return eigen_forRoeAvgs(solver, roe, x);
}

<?=eqn.waves_t?> eigen_leftTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> inputU,
	real3 x,
	normalInfo_t n
) {	
	inputU = cons_rotateFrom(inputU, n);
	
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;

<? for _,var in ipairs(eqn.eigenVars) do
?> 	<?=var.type or 'real'?> <?=var.name?> = eig.<?=var.name?>;
<? end ?>

#warning you can't use coordLenSq (which uses g_ij) after rotating coordinates 
	real vSq = coordLenSq(v, x);
	
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

	return result;
}

<?=eqn.cons_t?> eigen_rightTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> input,
	real3 x,
	normalInfo_t n
) {
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;

<? for _,var in ipairs(eqn.eigenVars) do
?> 	<?=var.type or 'real'?> <?=var.name?> = eig.<?=var.name?>;
<? end ?>

	real vSq = coordLenSq(v, x);
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
	resultU.B.x = 0;
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
	resultU.psi = 0;
	return cons_rotateTo(resultU, n);
}

<?=eqn.cons_t?> eigen_fluxTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> inputU,
	real3 x,
	normalInfo_t n
) {
	inputU = cons_rotateFrom(inputU, n);

	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	const real gamma_3 = gamma - 3.;

<? for _,var in ipairs(eqn.eigenVars) do
?> 	<?=var.type or 'real'?> <?=var.name?> = eig.<?=var.name?>;
<? end ?>

	real _1_rho = 1. / rho;
#warning you can't use coordLenSq (which uses g_ij) after rotating coordinates 
	real vSq = coordLenSq(v, x);
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
	resultU.B.x = 0;
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
	resultU.psi = 0;
	return cons_rotateTo(resultU, n);
}

<?=eqn.eigen_t?> eigen_forCell(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normalInfo_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real PMag = .5 * coordLenSq(W.B, x);
	real hTotal = (U.ETotal + W.P + PMag) / W.rho;
	<?=eqn.roe_t?> roe = {
		.rho = W.rho,
		.v = W.v,
		.hTotal = hTotal,
		.B = W.B,
		.X = 0,
		.Y = 1,
	};
	return eigen_forRoeAvgs(solver, roe, x);
}

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cellBuf[index].pos;
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

<? if not require 'hydro.coord.cartesian'.is(solver.coord) then ?>
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real BSq = coordLenSq(U->B, x);
	real PMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);
	real PTotal = W.P + PMag;
	real3 m_conn_vv = coord_conn_apply23(W.v, U->m, x);
	deriv->m = real3_sub(deriv->m, m_conn_vv);	//-Conn^i_jk rho v^j v^k 
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_raise(coord_conn_trace13(x), x), PTotal));		//+Conn^j_kj g^ki PTotal
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_conn_apply23(U->B, U->B, x), 1. / (solver->mu0 / unit_kg_m_per_C2)));	//+ 1/mu0 Conn^i_jk B^j B^k
<? end ?>
}

kernel void constrainU(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);

	W.rho = max(W.rho, 1e-7);
	W.P = max(W.P, 1e-7);

	*U = consFromPrim(solver, W, x);
}
