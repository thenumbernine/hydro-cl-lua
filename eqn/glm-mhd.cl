/*
started with the mhd.lua and mhd.cl
tweaked it while looking at
2006 Xueshang et al - A 3rd Order WENO GLM-MHD Scheme for Magnetic Reconnection
2009 Mignone, Tzeferacos - A Second-Order Unsplit Godunov Scheme for Cell-Centered MHD- the CTU-GLM scheme
*/

constant real c_h = 1.;
constant real c_p = 1.;

//use Eqn.hasFluxFromCons to allow the calcFlux function to take advantage of this function
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real vj = W.v.s<?=side?>;
	real Bj = W.B.s<?=side?>;
	real BSq = real3_lenSq(W.B);
	real BDotV = real3_dot(W.B, W.v);
	real PMag = .5 * BSq;
	real PTotal = W.P + PMag;
	real HTotal = U.ETotal + PTotal;
	
	<?=eqn.cons_t?> F;
	F.rho = U.m.s<?=side?>;
	F.m = real3_sub(real3_scale(U.m, vj), real3_scale(U.B, Bj / mu0));
	F.m.s<?=side?> += PTotal;
	F.B = real3_sub(real3_scale(U.B, vj), real3_scale(W.v, Bj));
	F.B.s<?=side?> += F.BPot;
	F.ETotal = HTotal * vj - BDotV * Bj / mu0;
	F.BPot = c_h * c_h;
	return F;
}
<? end ?>

<? for side=0,2 do ?>
<?=eqn.cons_t?> cons_swapFrom<?=side?>(<?=eqn.cons_t?> U) {
	U.m = real3_swap<?=side?>(U.m);
	U.B = real3_swap<?=side?>(U.B);
	return U;
}

<?=eqn.cons_t?> cons_swapTo<?=side?>(<?=eqn.cons_t?> U) {
	U.m = real3_swap<?=side?>(U.m);
	U.B = real3_swap<?=side?>(U.B);
	return U;
}
<? end ?>

//called from calcDT
<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.cons_t?> U_ = cons_swapFrom<?=side?>(*U);
	<?=eqn.prim_t?> W = primFromCons(U_, x);
	
#if 0
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real3 v = W.v;
	real3 B = W.B;
	
	real BSq = real3_lenSq(B);
	real invRho = 1./W.rho;
	
	real aSq = heatCapacityRatio * W.P * invRho;
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
	const real gamma = heatCapacityRatio;
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

	real lambdaFastMin = v.x - Cf;
	real lambdaFastMax = v.x + Cf;
	
	return (range_t){
		.min = min(lambdaFastMin, -c_h),
		.max = max(lambdaFastMax, c_h),
	};
#endif
}
<? end ?>

//assumes UL and UR are already rotated so the 'x' direction is our flux direction
void calcRoeValues(
	Roe_t* W, 
	const <?=eqn.cons_t?>* UL, 
	const <?=eqn.cons_t?>* UR,
	real3 x
) {
	// should I use Bx, or BxL/R, for calculating the PMag at the L and R states?
	<?=eqn.prim_t?> WL = primFromCons(*UL, x);
	real sqrtRhoL = sqrt(UL->rho);
	real PMagL = .5 * real3_lenSq(UL->B);
	real hTotalL = (UL->ETotal + WL.P + PMagL) / UL->rho;

	<?=eqn.prim_t?> WR = primFromCons(*UR, x);
	real sqrtRhoR = sqrt(UR->rho);
	real PMagR = .5 * real3_lenSq(UR->B);
	real hTotalR = (UR->ETotal + WR.P + PMagR) / UR->rho;
	
	real dby = WL.B.y - WR.B.y;
	real dbz = WL.B.z - WR.B.z;
	
	real invDenom = 1 / (sqrtRhoL + sqrtRhoR);
	
	W->rho  = sqrtRhoL * sqrtRhoR;
	W->v = real3_scale(real3_add(
		real3_scale(WL.v, sqrtRhoL),
		real3_scale(WR.v, sqrtRhoR)), invDenom);
	
	W->hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
	
	W->B.x = (sqrtRhoL * WL.B.x + sqrtRhoR * WR.B.x) * invDenom;
	// why does athena switch the weights of the By and Bz components?
	W->B.y = (sqrtRhoR * WL.B.y + sqrtRhoL * WR.B.y) * invDenom;
	W->B.z = (sqrtRhoR * WL.B.z + sqrtRhoL * WR.B.z) * invDenom;
	
	W->X = .5 * (dby * dby + dbz * dbz) * invDenom * invDenom;
	W->Y = .5 * (UL->rho + UR->rho) / W->rho;
};

<? 
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do 
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr0?> real* wave,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	real3 x
) {
	wave[0] = eig->v.x - eig->Cf;
	wave[1] = eig->v.x - eig->CAx;
	wave[2] = eig->v.x - eig->Cs;
	wave[3] = eig->v.x;
	wave[4] = eig->v.x + eig->Cs;
	wave[5] = eig->v.x + eig->CAx;
	wave[6] = eig->v.x + eig->Cf;
	wave[7] = -c_h;
	wave[8] = c_h;
}
<?
		end
	end 
end

for side=0,solver.dim-1 do
	for _,addr1 in ipairs{'', 'global'} do 
?>
void eig_forSide_<?=side?>_<?=addr1?>(
	<?=addr1?> <?=eqn.eigen_t?>* eig,
	Roe_t roe,
	real3 x
) {
	const real gamma = heatCapacityRatio;
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
	eig->hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// eig->hHydro = hTotal - CASq, CASq = EMag/rho
	// eig->hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// eig->hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(eig->hHydro - eKin) = gamma P / rho
	eig->aTildeSq = max((gamma_1 * (eig->hHydro - .5 * vSq) - gamma_2 * X), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + eig->aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - eig->aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + eig->aTildeSq * BStarPerpSq_rho);
	
	eig->CAx = sqrt(CAxSq);
	
	real CfSq = CStarSq + sqrtDiscr;
	eig->Cf = sqrt(CfSq);
	
	real CsSq = eig->aTildeSq * CAxSq / CfSq;
	eig->Cs = sqrt(CsSq);

	
	real BPerpLen = sqrt(BPerpSq);
	eig->BStarPerpLen = sqrt(BStarPerpSq);
	
	if (BPerpLen == 0) {
		eig->betaY = 1;
		eig->betaZ = 0;
	} else {
		eig->betaY = B.y / BPerpLen;
		eig->betaZ = B.z / BPerpLen;
	}
	eig->betaStarY = eig->betaY / sqrt(gamma_1 - gamma_2*Y);
	eig->betaStarZ = eig->betaZ / sqrt(gamma_1 - gamma_2*Y);
	eig->betaStarSq = eig->betaStarY*eig->betaStarY + eig->betaStarZ*eig->betaStarZ;


	if (CfSq - CsSq == 0) {
		eig->alphaF = 1;
		eig->alphaS = 0;
	} else if (eig->aTildeSq - CsSq <= 0) {
		eig->alphaF = 0;
		eig->alphaS = 1;
	} else if (CfSq - eig->aTildeSq <= 0) {
		eig->alphaF = 1;
		eig->alphaS = 0;
	} else {
		eig->alphaF = sqrt((eig->aTildeSq - CsSq) / (CfSq - CsSq));
		eig->alphaS = sqrt((CfSq - eig->aTildeSq) / (CfSq - CsSq));
	}


	eig->sqrtRho = sqrt(rho);
	real _1_sqrtRho = 1. / eig->sqrtRho;
	eig->sbx = B.x >= 0 ? 1 : -1;
	real aTilde = sqrt(eig->aTildeSq);
	eig->Qf = eig->Cf * eig->alphaF * eig->sbx;
	eig->Qs = eig->Cs * eig->alphaS * eig->sbx;
	eig->Af = aTilde * eig->alphaF * _1_sqrtRho;
	eig->As = aTilde * eig->alphaS * _1_sqrtRho;


	//used for eigenvectors and eigenvalues
<? 	for _,kv in ipairs(eqn.roeVars) do
		local name = next(kv) 
?>	eig-><?=name?> = roe.<?=name?>;
<?	end
?>

}
<? 
	end 
end
?>

kernel void calcEigenBasis(
	global real* waveBuf,			//[volume][dim][numWaves]
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

		//swap the sides with x here, so all the fluxes are in the 'x' direction
		<?=eqn.cons_t?> UL_ = cons_swapFrom<?=side?>(*UL);
		<?=eqn.cons_t?> UR_ = cons_swapFrom<?=side?>(*UR);

		Roe_t roe;
		calcRoeValues(&roe, &UL_, &UR_, xInt);

		global real* wave = waveBuf + numWaves * indexInt;
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;

		eig_forSide_<?=side?>_global(eig, roe, xInt);
		
		eigen_calcWaves_<?=side?>_global_global(wave, eig, xInt);
	}<? end ?>
}

<? for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,2 do ?>
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* result,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 x
) {	
	<?=eqn.cons_t?> inputU = cons_swapTo<?=side?>(*(<?=addr2?> <?=eqn.cons_t?>*)input);
	
	const real gamma = heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;

<? for _,kv in ipairs(eqn.eigenVars) do
	local name, ctype = next(kv)
?> 	<?=ctype?> <?=name?> = eig-><?=name?>;
<? end ?>

	real vSq = real3_lenSq(v);
	
	// left eigenvectors
	real norm = .5 / eig->aTildeSq;
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
	

	result[0] = 
		  inputU.rho * (alphaF * (vSq - eig->hHydro) + Cff * (Cf + v.x) - Qs * vqstr - aspb) 
		+ inputU.m.x * (-alphaF * v.x - Cff)
		+ inputU.m.y * (-alphaF * v.y + Qs * QStarY)
		+ inputU.m.z * (-alphaF * v.z + Qs * QStarZ)
		+ inputU.ETotal * alphaF
		+ inputU.B.y * l16
		+ inputU.B.z * l17;
	result[1] = 
		  inputU.rho * l21
		+ inputU.m.y * l23
		+ inputU.m.z * l24
		+ inputU.B.y * l26
		+ inputU.B.z * l27;
	result[2] = 
		  inputU.rho * (alphaS * (vSq - eig->hHydro) + Css * (Cs + v.x) + Qf * vqstr + afpb)
		+ inputU.m.x * (-alphaS * v.x - Css)
		+ inputU.m.y * (-alphaS * v.y - Qf * QStarY)
		+ inputU.m.z * (-alphaS * v.z - Qf * QStarZ)
		+ inputU.ETotal * alphaS
		+ inputU.B.y * l36
		+ inputU.B.z * l37;
	result[3] = 
		  inputU.rho * (1. - norm * (.5 * vSq - gamma_2 * X / gamma_1))
		+ inputU.m.x * norm*v.x
		+ inputU.m.y * norm*v.y
		+ inputU.m.z * norm*v.z
		+ inputU.ETotal * -norm
		+ inputU.B.y * norm*B.y
		+ inputU.B.z * norm*B.z;
	result[4] = 
		  inputU.rho * (alphaS * (vSq - eig->hHydro) + Css * (Cs - v.x) - Qf * vqstr + afpb)
		+ inputU.m.x * (-alphaS * v.x + Css)
		+ inputU.m.y * (-alphaS * v.y + Qf * QStarY)
		+ inputU.m.z * (-alphaS * v.z + Qf * QStarZ)
		+ inputU.ETotal * alphaS
		+ inputU.B.y * l36
		+ inputU.B.z * l37;
	result[5] = 
		  inputU.rho * -l21
		+ inputU.m.y * -l23
		+ inputU.m.z * -l24
		+ inputU.B.y * l26
		+ inputU.B.z * l27;
	result[6] = 
		  inputU.rho * (alphaF * (vSq - eig->hHydro) + Cff * (Cf - v.x) + Qs * vqstr - aspb)
		+ inputU.m.x * (-alphaF * v.x + Cff)
		+ inputU.m.y * (-alphaF * v.y - Qs * QStarY)
		+ inputU.m.z * (-alphaF * v.z - Qs * QStarZ)
		+ inputU.ETotal * alphaF
		+ inputU.B.y * l16
		+ inputU.B.z * l17;
	result[7] = inputU.B.x;
	result[8] = inputU.BPot;
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* result,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 x
) {
	const real gamma = heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;

<? for _,kv in ipairs(eqn.eigenVars) do
	local name, ctype = next(kv)
?> 	<?=ctype?> <?=name?> = eig-><?=name?>;
<? end ?>

	real vSq = real3_lenSq(v);
	real vDotBeta = v.y*betaStarY + v.z*betaStarZ;
	real _1_sqrtRho = 1. / eig->sqrtRho;
	real Afpbb = Af * BStarPerpLen * betaStarSq;
	real Aspbb = As * BStarPerpLen * betaStarSq;

	real lambdaFastMin = eig->v.x - eig->Cf;
	real lambdaSlowMin = eig->v.x - eig->Cs;
	real lambdaSlowMax = eig->v.x + eig->Cs;
	real lambdaFastMax = eig->v.x + eig->Cf;

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
		  input[0] * alphaF
		+ input[2] * alphaS
		+ input[3]
		+ input[4] * alphaS
		+ input[6] * alphaF;
	resultU.m.x =
		  input[0] * alphaF * lambdaFastMin
		+ input[2] * alphaS * lambdaSlowMin
		+ input[3] * v.x
		+ input[4] * alphaS * lambdaSlowMax
		+ input[6] * alphaF * lambdaFastMax;
	resultU.m.y =
		  input[0] * (qa3 + qc3)
		+ input[1] * -betaZ
		+ input[2] * (qb3 - qd3)
		+ input[3] * v.y
		+ input[4] * (qb3 + qd3)
		+ input[5] * betaZ
		+ input[6] * (qa3 - qc3);
	resultU.m.z =
		  input[0] * (qa4 + qc4)
		+ input[1] * betaY
		+ input[2] * (qb4 - qd4)
		+ input[3] * v.z
		+ input[4] * (qb4 + qd4)
		+ input[5] * -betaY
		+ input[6] * (qa4 - qc4);
	resultU.ETotal =
		  input[0] * (alphaF*(eig->hHydro - v.x*Cf) + Qs*vDotBeta + Aspbb)
		+ input[1] * r52
		+ input[2] * (alphaS*(eig->hHydro - v.x*Cs) - Qf*vDotBeta - Afpbb)
		+ input[3] * (.5*vSq + gamma_2*X/gamma_1)
		+ input[4] * (alphaS*(eig->hHydro + v.x*Cs) + Qf*vDotBeta - Afpbb)
		+ input[5] * -r52
		+ input[6] * (alphaF*(eig->hHydro + v.x*Cf) - Qs*vDotBeta + Aspbb);
	resultU.B.x = input[7];
	resultU.B.y =
		  input[0] * r61
		+ input[1] * r62
		+ input[2] * r63
		+ input[4] * r63
		+ input[5] * r62
		+ input[6] * r61;
	resultU.B.z =
		  input[0] * r71
		+ input[1] * r72
		+ input[2] * r73
		+ input[4] * r73
		+ input[5] * r72
		+ input[6] * r71;
	resultU.BPot = input[8];
	*(<?=addr0?> <?=eqn.cons_t?>*)result = cons_swapFrom<?=side?>(resultU);
}

<? 				if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* result,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input,
	real3 x
) {
	<?=eqn.cons_t?> inputU = cons_swapTo<?=side?>(*(<?=addr2?> <?=eqn.cons_t?>*)input);

	const real gamma = heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	const real gamma_3 = gamma - 3.;

<? for _,kv in ipairs(eqn.eigenVars) do
	local name, ctype = next(kv)
?> 	<?=ctype?> <?=name?> = eig-><?=name?>;
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
	resultU.BPot = inputU.BPot;
	*(<?=addr0?> <?=eqn.cons_t?>*)result = cons_swapFrom<?=side?>(resultU);
}
<?				end
			end
		end
	end
end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

	deriv->BPot -= c_h * c_h / (c_p * c_p) * U->BPot;
}


<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	<?=eqn.eigen_t?>* eig,
	global const <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real PMag = .5 * real3_lenSq(W.B);
	real hTotal = (U->ETotal + W.P + PMag) / W.rho;
	Roe_t roe = (Roe_t){
		.rho = W.rho,
		.v = W.v,
		.hTotal = hTotal,
		.B = W.B,
		.X = 0,
		.Y = 1,
	};
	eig_forSide_<?=side?>_(eig, roe, x);
}
<? end ?>

//U = output
//WA = W components that make up the jacobian matrix
//W = input
//x = coordinate location
void apply_dU_dW(
	<?=eqn.cons_t?>* U, 
	const <?=eqn.prim_t?>* WA, 
	const <?=eqn.prim_t?>* W, 
	real3 x
) {
	*U = (<?=eqn.cons_t?>){
		.rho = W->rho,
		.m = real3_add(
			real3_scale(WA->v, W->rho),
			real3_scale(W->v, WA->rho)),
		.B = WA->B,
		.ETotal = W->rho * .5 * real3_dot(WA->v, WA->v)
			+ WA->rho * real3_dot(W->v, WA->v)
			+ real3_dot(W->B, WA->B) / mu0
			+ W->P / (heatCapacityRatio - 1.),
		.BPot = W->BPot,
	};
}

//W = output
//WA = W components that make up the jacobian matrix
//U = input
//x = coordinate location
void apply_dW_dU(
	<?=eqn.prim_t?>* W,
	const <?=eqn.prim_t?>* WA,
	const <?=eqn.cons_t?>* U,
	real3 x
) {
	*W = (<?=eqn.prim_t?>){
		.rho = U->rho,
		.v = real3_sub(
			real3_scale(U->m, 1. / WA->rho),
			real3_scale(WA->v, U->rho / WA->rho)),
		.B = U->B,
		.P = (heatCapacityRatio - 1.) *  (
			.5 * U->rho * real3_dot(WA->v, WA->v)
			- real3_dot(U->m, WA->v)
			- real3_dot(U->B, WA->B) / mu0
			+ U->ETotal),
		.BPot = U->BPot,
	};
}
