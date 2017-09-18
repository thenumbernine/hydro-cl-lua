/*
Stone et al 2008 - https://arxiv.org/pdf/0804.0402v1.pdf
based on Athena's version of eigenvectors of derivative of adiabatic MHD flux wrt primitives
*/

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
	F.ETotal = HTotal * vj - BDotV * Bj / mu0;
	return F;
}
<? end ?>

<? 

-- [[ rotating using 2D rotations
-- rotate 3D vectors into the x-is-fwd plane
local function putXFwd(dst, src, side)
	return ({
		[0] = '',
		[1] = 'dst = _real3(src.y, -src.x, src.z);',
		[2] = 'dst = _real3(src.z, src.y, -src.x);',
	})[side]:gsub('src', src):gsub('dst', dst)
end 

-- reverse the above rotation 
local function putXBack(dst, src, side)
	return ({
		[0] = '',
		[1] = 'dst = _real3(-src.y, src.x, src.z);',
		[2] = 'dst = _real3(-src.z, src.y, src.x);',
	})[side]:gsub('src', src):gsub('dst', dst)
end 
--]]
--[[ rotate using 3-permutations
-- rotate 3D vectors into the x-is-fwd plane
local function putXFwd(dst, src, side)
	return ({
		[0] = '',
		[1] = 'dst = _real3(src.y, src.z, src.x);',
		[2] = 'dst = _real3(src.z, src.x, src.y);',
	})[side]:gsub('src', src):gsub('dst', dst)
end 

-- reverse the above rotation 
local function putXBack(dst, src, side)
	return ({
		[0] = '',
		[1] = 'dst = _real3(src.z, src.x, src.y);',
		[2] = 'dst = _real3(src.y, src.z, src.x);',
	})[side]:gsub('src', src):gsub('dst', dst)
end 
--]]

local function consPutXFwd(var, side)
	return putXFwd(var..'.m', var..'.m', side)..'\n'
		..putXFwd(var..'.B', var..'.B', side)
end

local function consPutXBack(var, side)
	return putXBack(var..'.m', var..'.m', side)..'\n'
		..putXBack(var..'.B', var..'.B', side)
end

-- rotate 3D vectors of a cons_t into the x-is-fwd plane, and remove the Bx component 
local function _7to8code(addr,side)
	return [[
	]]..eqn.cons_t..[[ inputU = *(]]..addr..[[ ]]..eqn.cons_t..[[*)input_;
	]]..consPutXFwd('inputU', side)..[[
	real input[7] = { inputU.rho, inputU.m.x, inputU.m.y, inputU.m.z, inputU.ETotal, inputU.B.y, inputU.B.z }; 
]]
end

-- re-insert Bx=0 and rotate x back to its original direction
local function _8to7code(addr, side)
	return [[
	]]..eqn.cons_t..[[ resultU = { 
		.rho = result[0], 
		.m = {.x = result[1], .y = result[2], .z = result[3] }, 
		.ETotal = result[4], 
		.B = {.x = 0, .y = result[5], .z = result[6] },
	};
	]]..consPutXBack('resultU', side)..[[
	*(]]..addr..[[ ]]..eqn.cons_t..[[*)result = resultU;
]]
end
?>

//called from calcDT
<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.cons_t?> U_ = *U;
	<?=consPutXFwd('U_', side)?>
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
		.min = lambdaFastMin,
		.max = lambdaFastMax,
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

<? for _,addr in ipairs{'', 'global'} do ?>
void fill_<?=addr?>(<?=addr?> real* ptr, int step, real a, real b, real c, real d, real e, real f, real g) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
	ptr[3*step] = d;
	ptr[4*step] = e;
	ptr[5*step] = f;
	ptr[6*step] = g;
}
<? 
end 
?>

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
	wave[0] = eig->vx - eig->Cf;
	wave[1] = eig->vx - eig->CAx;
	wave[2] = eig->vx - eig->Cs;
	wave[3] = eig->vx;
	wave[4] = eig->vx + eig->Cs;
	wave[5] = eig->vx + eig->CAx;
	wave[6] = eig->vx + eig->Cf;
}
<?
		end
	end 
end

for side=0,solver.dim-1 do
	for _,addr0 in ipairs{'', 'global'} do 
		for _,addr1 in ipairs{'', 'global'} do 
?>
void calcEigenSystemForRoeValues_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr0?> real* wave,
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
	real BDotV = real3_dot(B,v);
	real BPerpSq = B.y*B.y + B.z*B.z;
	real BStarPerpSq = (gamma_1 - gamma_2 * Y) * BPerpSq;
	real CAxSq = B.x*B.x*_1_rho;
	real CASq = CAxSq + BPerpSq * _1_rho;
	real hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// hHydro = hTotal - CASq, CASq = EMag/rho
	// hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2 * X), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * BStarPerpSq_rho);
	
	real CfSq = CStarSq + sqrtDiscr;
	real Cf = sqrt(CfSq);
	
	real CsSq = aTildeSq * CAxSq / CfSq;
	real Cs = sqrt(CsSq);
	
	real BPerpLen = sqrt(BPerpSq);
	real BStarPerpLen = sqrt(BStarPerpSq);
	real betaY, betaZ;
	if (BPerpLen == 0) {
		betaY = 1;
		betaZ = 0;
	} else {
		betaY = B.y / BPerpLen;
		betaZ = B.z / BPerpLen;
	}
	real betaStarY = betaY / sqrt(gamma_1 - gamma_2*Y);
	real betaStarZ = betaZ / sqrt(gamma_1 - gamma_2*Y);
	real betaStarSq = betaStarY*betaStarY + betaStarZ*betaStarZ;
	real vDotBeta = v.y*betaStarY + v.z*betaStarZ;

	real alphaF, alphaS;
	if (CfSq - CsSq == 0) {
		alphaF = 1;
		alphaS = 0;
	} else if (aTildeSq - CsSq <= 0) {
		alphaF = 0;
		alphaS = 1;
	} else if (CfSq - aTildeSq <= 0) {
		alphaF = 1;
		alphaS = 0;
	} else {
		alphaF = sqrt((aTildeSq - CsSq) / (CfSq - CsSq));
		alphaS = sqrt((CfSq - aTildeSq) / (CfSq - CsSq));
	}

	real sqrtRho = sqrt(rho);
	real _1_sqrtRho = 1./sqrtRho;
	real sbx = B.x >= 0 ? 1 : -1;
	real aTilde = sqrt(aTildeSq);
	real Qf = Cf*alphaF*sbx;
	real Qs = Cs*alphaS*sbx;
	real Af = aTilde*alphaF*_1_sqrtRho;
	real As = aTilde*alphaS*_1_sqrtRho;
	real Afpbb = Af*BStarPerpLen*betaStarSq;
	real Aspbb = As*BStarPerpLen*betaStarSq;

	eig->roe = roe;

	eig->vx = v.x;
	eig->Cs = Cs;
	eig->CAx = sqrt(CAxSq);
	eig->Cf = Cf;

	real lambdaFastMin, lambdaSlowMin, lambdaSlowMax, lambdaFastMax;
	if (wave) {
		eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(wave, eig, x);
		lambdaFastMin = wave[0];
		lambdaSlowMin = wave[2];
		lambdaSlowMax = wave[4];
		lambdaFastMax = wave[6];
	} else {
		lambdaFastMin = eig->vx - eig->Cf;
		lambdaSlowMin = eig->vx - eig->Cs;
		lambdaSlowMax = eig->vx + eig->Cs;
		lambdaFastMax = eig->vx + eig->Cf;
	}

	// dF/dU
	<? if solver.checkFluxError then ?>
	<?=addr1?> real* A = eig->A;
	fill_<?=addr1?>(A+0,7,	0, 												1,											0,								0,								0,			0,								0								);
	fill_<?=addr1?>(A+1,7,	-v.x*v.x + .5*gamma_1*vSq - gamma_2*X,			-gamma_3*v.x,								-gamma_1*v.y,					-gamma_1*v.z,					gamma_1,	-gamma_2*Y*B.y,					-gamma_2*Y*B.z					);
	fill_<?=addr1?>(A+2,7,	-v.x*v.y,										v.y,										v.x,							0, 								0,			-B.x,							0								);
	fill_<?=addr1?>(A+3,7,	-v.x*v.z,										v.z,										0,								v.x, 							0,			0,								-B.x							);
	fill_<?=addr1?>(A+4,7,	v.x*(.5*gamma_1*vSq - hTotal) + B.x*BDotV/rho,	-gamma_1*v.x*v.x + hTotal - B.x*B.x/rho,	-gamma_1*v.x*v.y - B.x*B.y/rho,	-gamma_1*v.x*v.z - B.x*B.z/rho,	gamma*v.x,	-gamma_2*Y*B.y*v.x - B.x*v.y,	-gamma_2*Y*B.z*v.x - B.x*v.z	);
	fill_<?=addr1?>(A+5,7,	(B.x*v.y - B.y*v.x)/rho,						B.y/rho,									-B.x/rho,						0, 								0,			v.x,							0								);
	fill_<?=addr1?>(A+6,7,	(B.x*v.z - B.z*v.x)/rho,						B.z/rho,									0,								-B.x/rho, 						0,			0,								v.x								);
	<? end ?>
}
<? 
		end
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
		<?=eqn.cons_t?> UL_ = *UL;
		<?=eqn.cons_t?> UR_ = *UR;
		<?=consPutXFwd('UL_',side)?>
		<?=consPutXFwd('UR_',side)?>

		Roe_t roe;
		calcRoeValues(&roe, &UL_, &UR_, xInt);

		global real* wave = waveBuf + numWaves * indexInt;
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;

		calcEigenSystemForRoeValues_<?=side?>_global_global(wave, eig, roe, xInt);
	}<? end ?>
}

<? for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,2 do ?>
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* result,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input_,
	real3 x
) {	
	<?=_7to8code(addr2, side)?>
	
	const real gamma = heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	const real gamma_3 = gamma - 3.;
	
	real rho = eig->roe.rho;
	real3 v = eig->roe.v;
	real hTotal = eig->roe.hTotal;
	real3 B = eig->roe.B;
	real X = eig->roe.X;
	real Y = eig->roe.Y;

	real _1_rho = 1. / rho;
	real vSq = real3_lenSq(v);
	real BDotV = real3_dot(B,v);
	real BPerpSq = B.y*B.y + B.z*B.z;
	real BStarPerpSq = (gamma_1 - gamma_2 * Y) * BPerpSq;
	real CAxSq = B.x*B.x*_1_rho;
	real CASq = CAxSq + BPerpSq * _1_rho;
	real hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// hHydro = hTotal - CASq, CASq = EMag/rho
	// hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2 * X), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * BStarPerpSq_rho);
	
	real CfSq = CStarSq + sqrtDiscr;
	real Cf = sqrt(CfSq);
	
	real CsSq = aTildeSq * CAxSq / CfSq;
	real Cs = sqrt(CsSq);
	
	real BPerpLen = sqrt(BPerpSq);
	real BStarPerpLen = sqrt(BStarPerpSq);
	real betaY, betaZ;
	if (BPerpLen == 0) {
		betaY = 1;
		betaZ = 0;
	} else {
		betaY = B.y / BPerpLen;
		betaZ = B.z / BPerpLen;
	}
	real betaStarY = betaY / sqrt(gamma_1 - gamma_2*Y);
	real betaStarZ = betaZ / sqrt(gamma_1 - gamma_2*Y);
	real betaStarSq = betaStarY*betaStarY + betaStarZ*betaStarZ;
	real vDotBeta = v.y*betaStarY + v.z*betaStarZ;

	real alphaF, alphaS;
	if (CfSq - CsSq == 0) {
		alphaF = 1;
		alphaS = 0;
	} else if (aTildeSq - CsSq <= 0) {
		alphaF = 0;
		alphaS = 1;
	} else if (CfSq - aTildeSq <= 0) {
		alphaF = 1;
		alphaS = 0;
	} else {
		alphaF = sqrt((aTildeSq - CsSq) / (CfSq - CsSq));
		alphaS = sqrt((CfSq - aTildeSq) / (CfSq - CsSq));
	}

	real sqrtRho = sqrt(rho);
	real _1_sqrtRho = 1./sqrtRho;
	real sbx = B.x >= 0 ? 1 : -1;
	real aTilde = sqrt(aTildeSq);
	real Qf = Cf*alphaF*sbx;
	real Qs = Cs*alphaS*sbx;
	real Af = aTilde*alphaF*_1_sqrtRho;
	real As = aTilde*alphaS*_1_sqrtRho;
	real Afpbb = Af*BStarPerpLen*betaStarSq;
	real Aspbb = As*BStarPerpLen*betaStarSq;


	real lambdaFastMin = eig->vx - eig->Cf;
	real lambdaSlowMin = eig->vx - eig->Cs;
	real lambdaSlowMax = eig->vx + eig->Cs;
	real lambdaFastMax = eig->vx + eig->Cf;


	// left eigenvectors
	real norm = .5/aTildeSq;
	real Cff = norm*alphaF*Cf;
	real Css = norm*alphaS*Cs;
	Qf = Qf * norm;
	Qs = Qs * norm;
	real AHatF = norm*Af*rho;
	real AHatS = norm*As*rho;
	real afpb = norm*Af*BStarPerpLen;
	real aspb = norm*As*BStarPerpLen;

	norm = norm * gamma_1;
	alphaF = alphaF * norm;
	alphaS = alphaS * norm;
	real QStarY = betaStarY/betaStarSq;
	real QStarZ = betaStarZ/betaStarSq;
	real vqstr = (v.y*QStarY + v.z*QStarZ);
	norm = norm * 2.;

	
	real l16 = AHatS*QStarY - alphaF*B.y;
	real l17 = AHatS*QStarZ - alphaF*B.z;
	real l21 = .5*(v.y*betaZ - v.z*betaY);
	real l23 = .5*betaZ;
	real l24 = .5*betaY;
	real l26 = -.5*sqrtRho*betaZ*sbx;
	real l27 = .5*sqrtRho*betaY*sbx;
	real l36 = -AHatF*QStarY - alphaS*B.y;
	real l37 = -AHatF*QStarZ - alphaS*B.z;
	

	result[0] = 
		  input[0] * (alphaF*(vSq-hHydro) + Cff*(Cf+v.x) - Qs*vqstr - aspb) 
		+ input[1] * (-alphaF*v.x - Cff)
		+ input[2] * (-alphaF*v.y + Qs*QStarY)
		+ input[3] * (-alphaF*v.z + Qs*QStarZ)
		+ input[4] * alphaF
		+ input[5] * l16
		+ input[6] * l17;
	result[1] = 
		  input[0] * l21
		+ input[2] * l23
		+ input[3] * l24
		+ input[5] * l26
		+ input[6] * l27;
	result[2] = 
		  input[0] * (alphaS*(vSq-hHydro) + Css*(Cs+v.x) + Qf*vqstr + afpb)
		+ input[1] * (-alphaS*v.x - Css)
		+ input[2] * (-alphaS*v.y - Qf*QStarY)
		+ input[3] * (-alphaS*v.z - Qf*QStarZ)
		+ input[4] * alphaS
		+ input[5] * l36
		+ input[6] * l37;
	result[3] = 
		  input[0] * (1. - norm*(.5*vSq - gamma_2*X/gamma_1))
		+ input[1] * norm*v.x
		+ input[2] * norm*v.y
		+ input[3] * norm*v.z
		+ input[4] * -norm
		+ input[5] * norm*B.y
		+ input[6] * norm*B.z;
	result[4] = 
		  input[0] * (alphaS*(vSq-hHydro) + Css*(Cs-v.x) - Qf*vqstr + afpb)
		+ input[1] * (-alphaS*v.x + Css)
		+ input[2] * (-alphaS*v.y + Qf*QStarY)
		+ input[3] * (-alphaS*v.z + Qf*QStarZ)
		+ input[4] * alphaS
		+ input[5] * l36
		+ input[6] * l37;
	result[5] = 
		  input[0] * -l21
		+ input[2] * -l23
		+ input[3] * -l24
		+ input[5] * l26
		+ input[6] * l27;
	result[6] = 
		  input[0] * (alphaF*(vSq-hHydro) + Cff*(Cf-v.x) + Qs*vqstr - aspb)
		+ input[1] * (-alphaF*v.x + Cff)
		+ input[2] * (-alphaF*v.y - Qs*QStarY)
		+ input[3] * (-alphaF*v.z - Qs*QStarZ)
		+ input[4] * alphaF
		+ input[5] * l16
		+ input[6] * l17;
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
	const real gamma_3 = gamma - 3.;
	
	real rho = eig->roe.rho;
	real3 v = eig->roe.v;
	real hTotal = eig->roe.hTotal;
	real3 B = eig->roe.B;
	real X = eig->roe.X;
	real Y = eig->roe.Y;

	real _1_rho = 1. / rho;
	real vSq = real3_lenSq(v);
	real BDotV = real3_dot(B,v);
	real BPerpSq = B.y*B.y + B.z*B.z;
	real BStarPerpSq = (gamma_1 - gamma_2 * Y) * BPerpSq;
	real CAxSq = B.x*B.x*_1_rho;
	real CASq = CAxSq + BPerpSq * _1_rho;
	real hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// hHydro = hTotal - CASq, CASq = EMag/rho
	// hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2 * X), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * BStarPerpSq_rho);
	
	real CfSq = CStarSq + sqrtDiscr;
	real Cf = sqrt(CfSq);
	
	real CsSq = aTildeSq * CAxSq / CfSq;
	real Cs = sqrt(CsSq);
	
	real BPerpLen = sqrt(BPerpSq);
	real BStarPerpLen = sqrt(BStarPerpSq);
	real betaY, betaZ;
	if (BPerpLen == 0) {
		betaY = 1;
		betaZ = 0;
	} else {
		betaY = B.y / BPerpLen;
		betaZ = B.z / BPerpLen;
	}
	real betaStarY = betaY / sqrt(gamma_1 - gamma_2*Y);
	real betaStarZ = betaZ / sqrt(gamma_1 - gamma_2*Y);
	real betaStarSq = betaStarY*betaStarY + betaStarZ*betaStarZ;
	real vDotBeta = v.y*betaStarY + v.z*betaStarZ;

	real alphaF, alphaS;
	if (CfSq - CsSq == 0) {
		alphaF = 1;
		alphaS = 0;
	} else if (aTildeSq - CsSq <= 0) {
		alphaF = 0;
		alphaS = 1;
	} else if (CfSq - aTildeSq <= 0) {
		alphaF = 1;
		alphaS = 0;
	} else {
		alphaF = sqrt((aTildeSq - CsSq) / (CfSq - CsSq));
		alphaS = sqrt((CfSq - aTildeSq) / (CfSq - CsSq));
	}

	real sqrtRho = sqrt(rho);
	real _1_sqrtRho = 1./sqrtRho;
	real sbx = B.x >= 0 ? 1 : -1;
	real aTilde = sqrt(aTildeSq);
	real Qf = Cf*alphaF*sbx;
	real Qs = Cs*alphaS*sbx;
	real Af = aTilde*alphaF*_1_sqrtRho;
	real As = aTilde*alphaS*_1_sqrtRho;
	real Afpbb = Af*BStarPerpLen*betaStarSq;
	real Aspbb = As*BStarPerpLen*betaStarSq;


	real lambdaFastMin = eig->vx - eig->Cf;
	real lambdaSlowMin = eig->vx - eig->Cs;
	real lambdaSlowMax = eig->vx + eig->Cs;
	real lambdaFastMax = eig->vx + eig->Cf;


	// right eigenvectors
	real qa3 = alphaF*v.y;
	real qb3 = alphaS*v.y;
	real qc3 = Qs*betaStarY;
	real qd3 = Qf*betaStarY;
	real qa4 = alphaF*v.z;
	real qb4 = alphaS*v.z;
	real qc4 = Qs*betaStarZ;
	real qd4 = Qf*betaStarZ;
	real r52 = -(v.y*betaZ - v.z*betaY);
	real r61 = As*betaStarY;
	real r62 = -betaZ*sbx*_1_sqrtRho;
	real r63 = -Af*betaStarY;
	real r71 = As*betaStarZ;
	real r72 = betaY*sbx*_1_sqrtRho;
	real r73 = -Af*betaStarZ;
	
	result[0] =
		  input[0] * alphaF
		+ input[2] * alphaS
		+ input[3]
		+ input[4] * alphaS
		+ input[6] * alphaF;
	result[1] =
		  input[0] * alphaF*lambdaFastMin
		+ input[2] * alphaS*lambdaSlowMin
		+ input[3] * v.x
		+ input[4] * alphaS*lambdaSlowMax
		+ input[6] * alphaF*lambdaFastMax;
	result[2] =
		  input[0] * (qa3 + qc3)
		+ input[1] * -betaZ
		+ input[2] * (qb3 - qd3)
		+ input[3] * v.y
		+ input[4] * (qb3 + qd3)
		+ input[5] * betaZ
		+ input[6] * (qa3 - qc3);
	result[3] =
		  input[0] * (qa4 + qc4)
		+ input[1] * betaY
		+ input[2] * (qb4 - qd4)
		+ input[3] * v.z
		+ input[4] * (qb4 + qd4)
		+ input[5] * -betaY
		+ input[6] * (qa4 - qc4);
	result[4] =
		  input[0] * (alphaF*(hHydro - v.x*Cf) + Qs*vDotBeta + Aspbb)
		+ input[1] * r52
		+ input[2] * (alphaS*(hHydro - v.x*Cs) - Qf*vDotBeta - Afpbb)
		+ input[3] * (.5*vSq + gamma_2*X/gamma_1)
		+ input[4] * (alphaS*(hHydro + v.x*Cs) + Qf*vDotBeta - Afpbb)
		+ input[5] * -r52
		+ input[6] * (alphaF*(hHydro + v.x*Cf) - Qs*vDotBeta + Aspbb);
	result[5] =
		  input[0] * r61
		+ input[1] * r62
		+ input[2] * r63
		+ input[4] * r63
		+ input[5] * r62
		+ input[6] * r61;
	result[6] =
		  input[0] * r71
		+ input[1] * r72
		+ input[2] * r73
		+ input[4] * r73
		+ input[5] * r72
		+ input[6] * r71;


	<?=_8to7code(addr0, side)?>
}

<? 				if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* result,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* input_,
	real3 x
) {
	<?=_7to8code(addr2, side)?>
	<?=addr1?> const real* A = eig->A;
	for (int i = 0; i < 7; ++i) {
		real sum = 0;
		for (int j = 0; j < 7; ++j) {
			sum += A[i+7*j] * X[j];
		}
		Y[i] = sum;
	}
	<?=_8to7code(addr0, side)?>
}
<?				end
			end
		end
	end
end ?>

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
	calcEigenSystemForRoeValues_<?=side?>__(NULL, eig, roe, x);
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
	};
}
