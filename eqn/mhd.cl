//use Eqn.hasFluxFromCons to allow the calcFlux function to take advantage of this function
<? for side=0,solver.dim-1 do ?>
cons_t fluxFromCons_<?=side?>(cons_t U) {
	prim_t W = primFromCons(U);
	real vj = W.v.s<?=side?>;
	real bj = W.b.s<?=side?>;
	real bSq = real3_lenSq(W.b);
	real bDotV = real3_dot(W.b, W.v);
	real PMag = .5 * bSq;
	real PTotal = W.P + PMag;
	real HTotal = U.ETotal + PTotal;
	
	cons_t F;
	F.rho = U.m.s<?=side?>;
	F.m = real3_sub(real3_scale(U.m, vj), real3_scale(U.b, bj / mu0));
	F.m.s<?=side?> += PTotal;
	F.b = real3_sub(real3_scale(U.b, vj), real3_scale(W.v, bj));
	F.ETotal = HTotal * vj - bDotV * bj / mu0;
	return F;
}
<? end ?>

// https://arxiv.org/pdf/0804.0402v1.pdf
// based on Athena's version of eigenvectors of derivative of adiabatic MHD flux wrt primitives

typedef struct {
	real rho;
	real3 v;
	real hTotal;
	real3 b;
	real X, Y;
} Roe_t;

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
		..putXFwd(var..'.b', var..'.b', side)
end

local function consPutXBack(var, side)
	return putXBack(var..'.m', var..'.m', side)..'\n'
		..putXBack(var..'.b', var..'.b', side)
end

-- rotate 3D vectors of a cons_t into the x-is-fwd plane, and remove the bx component 
local function _7to8code(addr,side)
	return [[
	cons_t xU = *(]]..addr..[[ cons_t*)x_;
	]]..consPutXFwd('xU', side)..[[
	real x[7] = { xU.rho, xU.m.x, xU.m.y, xU.m.z, xU.ETotal, xU.b.y, xU.b.z }; 
]]
end

-- re-insert bx=0 and rotate x back to its original direction
local function _8to7code(addr, side)
	return [[
	cons_t yU = { 
		.rho = y[0], 
		.m = {.x = y[1], .y = y[2], .z = y[3] }, 
		.ETotal = y[4], 
		.b = {.x = 0, .y = y[5], .z = y[6] },
	};
	]]..consPutXBack('yU', side)..[[
	*(]]..addr..[[ cons_t*)y = yU;
]]
end
?>

//called from calcDT
<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global cons_t* U
) {
	cons_t U_ = *U;
	<?=consPutXFwd('U_', side)?>
	prim_t W = primFromCons(U_);
	
#if 0
	prim_t W = primFromCons(*U);
	real3 v = W.v;
	real3 b = W.b;
	
	real bSq = real3_lenSq(b);
	real invRho = 1./W.rho;
	
	real aSq = heatCapacityRatio * W.P * invRho;
	real CaxSq = b.s<?=side?> * b.s<?=side?> * invRho;
	real CaSq = bSq * invRho;
	
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
	real3 b = W.b;
	real hTotal = .5 * real3_lenSq(W.v) + (W.P * gamma / gamma_1 + real3_lenSq(b)) / W.rho;

	//the rest of this matches calcEigenBasis:

	real _1_rho = 1. / rho;
	real vSq = real3_lenSq(v);
	real bPerpSq = b.y*b.y + b.z*b.z;
	real bStarPerpSq = (gamma_1 - gamma_2) * bPerpSq;
	real CAxSq = b.x*b.x*_1_rho;
	real CASq = CAxSq + bPerpSq * _1_rho;
	real hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// hHydro = hTotal - CASq, CASq = EMag/rho
	// hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2), 1e-20);

	real bStarPerpSq_rho = bStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + bStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * bStarPerpSq_rho);
	
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
	const cons_t* UL, 
	const cons_t* UR
) {
	// should I use bx, or bxL/R, for calculating the PMag at the L and R states?
	prim_t WL = primFromCons(*UL);
	real sqrtRhoL = sqrt(UL->rho);
	real PMagL = .5 * real3_lenSq(UL->b);
	real hTotalL = (UL->ETotal + WL.P + PMagL) / UL->rho;

	prim_t WR = primFromCons(*UR);
	real sqrtRhoR = sqrt(UR->rho);
	real PMagR = .5 * real3_lenSq(UR->b);
	real hTotalR = (UR->ETotal + WR.P + PMagR) / UR->rho;
	
	real dby = WL.b.y - WR.b.y;
	real dbz = WL.b.z - WR.b.z;
	
	real invDenom = 1 / (sqrtRhoL + sqrtRhoR);
	
	W->rho  = sqrtRhoL * sqrtRhoR;
	W->v = real3_scale(real3_add(
		real3_scale(WL.v, sqrtRhoL),
		real3_scale(WR.v, sqrtRhoR)), invDenom);
	
	W->hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
	
	W->b.x = (sqrtRhoL * WL.b.x + sqrtRhoR * WR.b.x) * invDenom;
	// why does athena switch the weights of the by and bz components?
	W->b.y = (sqrtRhoR * WL.b.y + sqrtRhoL * WR.b.y) * invDenom;
	W->b.z = (sqrtRhoR * WL.b.z + sqrtRhoL * WR.b.z) * invDenom;
	
	W->X = .5 * (dby * dby + dbz * dbz) * invDenom * invDenom;
	W->Y = .5 * (UL->rho + UR->rho) / W->rho;
};

void fill(global real* ptr, int step, real a, real b, real c, real d, real e, real f, real g) {
	ptr[0*step] = a;
	ptr[1*step] = b;
	ptr[2*step] = c;
	ptr[3*step] = d;
	ptr[4*step] = e;
	ptr[5*step] = f;
	ptr[6*step] = g;
}

kernel void calcEigenBasis(
	global real* waveBuf,			//[volume][dim][numWaves]
	global eigen_t* eigenBuf,		//[volume][dim]
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(2,1);
	int indexR = index;

	const real gamma = heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	const real gamma_3 = gamma - 3.;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
		<?= solver.getULRCode ?>

		//swap the sides with x here, so all the fluxes are in the 'x' direction
		cons_t UL_ = *UL;
		cons_t UR_ = *UR;
		<?=consPutXFwd('UL_',side)?>
		<?=consPutXFwd('UR_',side)?>

		Roe_t roe;
		calcRoeValues(&roe, &UL_, &UR_);

		real rho = roe.rho;
		real3 v = roe.v;
		real hTotal = roe.hTotal;
		real3 b = roe.b;
		real X = roe.X;
		real Y = roe.Y;

		real _1_rho = 1. / rho;
		real vSq = real3_lenSq(v);
		real bDotV = real3_dot(b,v);
		real bPerpSq = b.y*b.y + b.z*b.z;
		real bStarPerpSq = (gamma_1 - gamma_2 * Y) * bPerpSq;
		real CAxSq = b.x*b.x*_1_rho;
		real CASq = CAxSq + bPerpSq * _1_rho;
		real hHydro = hTotal - CASq;
		// hTotal = (EHydro + EMag + P)/rho
		// hHydro = hTotal - CASq, CASq = EMag/rho
		// hHydro = eHydro + P/rho = eKin + eInt + P/rho
		// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
		// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
		real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2 * X), 1e-20);

		real bStarPerpSq_rho = bStarPerpSq * _1_rho;
		real CATildeSq = CAxSq + bStarPerpSq_rho;
		real CStarSq = .5 * (CATildeSq + aTildeSq);
		real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
		real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * bStarPerpSq_rho);
		
		real CfSq = CStarSq + sqrtDiscr;
		real Cf = sqrt(CfSq);
		
		real CsSq = aTildeSq * CAxSq / CfSq;
		real Cs = sqrt(CsSq);
		
		real bPerpLen = sqrt(bPerpSq);
		real bStarPerpLen = sqrt(bStarPerpSq);
		real betaY, betaZ;
		if (bPerpLen == 0) {
			betaY = 1;
			betaZ = 0;
		} else {
			betaY = b.y / bPerpLen;
			betaZ = b.z / bPerpLen;
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
		real sbx = b.x >= 0 ? 1 : -1;
		real aTilde = sqrt(aTildeSq);
		real Qf = Cf*alphaF*sbx;
		real Qs = Cs*alphaS*sbx;
		real Af = aTilde*alphaF*_1_sqrtRho;
		real As = aTilde*alphaS*_1_sqrtRho;
		real Afpbb = Af*bStarPerpLen*betaStarSq;
		real Aspbb = As*bStarPerpLen*betaStarSq;

		real CAx = sqrt(CAxSq);
		
		int intindex = side + dim * index;
		global real* wave = waveBuf + numWaves * intindex;
		
		real lambdaFastMin = v.x - Cf;
		real lambdaSlowMin = v.x - Cs;
		real lambdaSlowMax = v.x + Cs;
		real lambdaFastMax = v.x + Cf;
		fill(wave, 1, lambdaFastMin, v.x - CAx, lambdaSlowMin, v.x, lambdaSlowMax, v.x + CAx, lambdaFastMax);

		global eigen_t* eig = eigenBuf + intindex;

		// dF/dU
		<? if solver.checkFluxError then ?>
		global real* A = eig->A;
		fill(A+0,7,	0, 												1,											0,								0,								0,			0,								0								);
		fill(A+1,7,	-v.x*v.x + .5*gamma_1*vSq - gamma_2*X,			-gamma_3*v.x,								-gamma_1*v.y,					-gamma_1*v.z,					gamma_1,	-gamma_2*Y*b.y,					-gamma_2*Y*b.z					);
		fill(A+2,7,	-v.x*v.y,										v.y,										v.x,							0, 								0,			-b.x,							0								);
		fill(A+3,7,	-v.x*v.z,										v.z,										0,								v.x, 							0,			0,								-b.x							);
		fill(A+4,7,	v.x*(.5*gamma_1*vSq - hTotal) + b.x*bDotV/rho,	-gamma_1*v.x*v.x + hTotal - b.x*b.x/rho,	-gamma_1*v.x*v.y - b.x*b.y/rho,	-gamma_1*v.x*v.z - b.x*b.z/rho,	gamma*v.x,	-gamma_2*Y*b.y*v.x - b.x*v.y,	-gamma_2*Y*b.z*v.x - b.x*v.z	);
		fill(A+5,7,	(b.x*v.y - b.y*v.x)/rho,						b.y/rho,									-b.x/rho,						0, 								0,			v.x,							0								);
		fill(A+6,7,	(b.x*v.z - b.z*v.x)/rho,						b.z/rho,									0,								-b.x/rho, 						0,			0,								v.x								);
		<? end ?>

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
		//rows
		global real* evR = eig->evR;
		fill(evR+0,7, 	alphaF, 										0, 		alphaS, 										1, 							alphaS, 										0, 		alphaF											);
		fill(evR+1,7, 	alphaF*lambdaFastMin, 							0, 		alphaS*lambdaSlowMin, 							v.x, 	 	 	 	 	 	alphaS*lambdaSlowMax, 							0, 		alphaF*lambdaFastMax							);
		fill(evR+2,7, 	qa3 + qc3, 										-betaZ,	qb3 - qd3, 										v.y, 	 	 	 	 	 	qb3 + qd3, 										betaZ, 	qa3 - qc3										);
		fill(evR+3,7, 	qa4 + qc4, 										betaY,	qb4 - qd4, 										v.z, 	 	 	 	 	 	qb4 + qd4, 										-betaY, qa4 - qc4										);
		fill(evR+4,7, 	alphaF*(hHydro - v.x*Cf) + Qs*vDotBeta + Aspbb, r52,	alphaS*(hHydro - v.x*Cs) - Qf*vDotBeta - Afpbb, .5*vSq + gamma_2*X/gamma_1, alphaS*(hHydro + v.x*Cs) + Qf*vDotBeta - Afpbb, -r52, 	alphaF*(hHydro + v.x*Cf) - Qs*vDotBeta + Aspbb	);
		fill(evR+5,7, 	r61, 											r62, 	r63, 											0, 							r63, 											r62, 	r61												);
		fill(evR+6,7,	r71, 											r72, 	r73, 											0, 							r73, 											r72, 	r71												);

		// left eigenvectors
		real norm = .5/aTildeSq;
		real Cff = norm*alphaF*Cf;
		real Css = norm*alphaS*Cs;
		Qf = Qf * norm;
		Qs = Qs * norm;
		real AHatF = norm*Af*rho;
		real AHatS = norm*As*rho;
		real afpb = norm*Af*bStarPerpLen;
		real aspb = norm*As*bStarPerpLen;

		norm = norm * gamma_1;
		alphaF = alphaF * norm;
		alphaS = alphaS * norm;
		real QStarY = betaStarY/betaStarSq;
		real QStarZ = betaStarZ/betaStarSq;
		real vqstr = (v.y*QStarY + v.z*QStarZ);
		norm = norm * 2.;

#if 1	//computing beforehand
		real l16 = AHatS*QStarY - alphaF*b.y;
		real l17 = AHatS*QStarZ - alphaF*b.z;
		real l21 = .5*(v.y*betaZ - v.z*betaY);
		real l23 = .5*betaZ;
		real l24 = .5*betaY;
		real l26 = -.5*sqrtRho*betaZ*sbx;
		real l27 = .5*sqrtRho*betaY*sbx;
		real l36 = -AHatF*QStarY - alphaS*b.y;
		real l37 = -AHatF*QStarZ - alphaS*b.z;
		//rows
		global real* evL = eig->evL;
		fill(evL+0,7,	alphaF*(vSq-hHydro) + Cff*(Cf+v.x) - Qs*vqstr - aspb, -alphaF*v.x - Cff, -alphaF*v.y + Qs*QStarY, -alphaF*v.z + Qs*QStarZ, alphaF, l16, l17);
		fill(evL+1,7,	l21, 0, l23, l24, 0, l26, l27);
		fill(evL+2,7,	alphaS*(vSq-hHydro) + Css*(Cs+v.x) + Qf*vqstr + afpb, -alphaS*v.x - Css, -alphaS*v.y - Qf*QStarY, -alphaS*v.z - Qf*QStarZ, alphaS, l36, l37);
		fill(evL+3,7,	1. - norm*(.5*vSq - gamma_2*X/gamma_1) , norm*v.x, norm*v.y, norm*v.z, -norm, norm*b.y, norm*b.z);
		fill(evL+4,7,	alphaS*(vSq-hHydro) + Css*(Cs-v.x) - Qf*vqstr + afpb, -alphaS*v.x + Css, -alphaS*v.y + Qf*QStarY, -alphaS*v.z + Qf*QStarZ, alphaS, l36, l37);
		fill(evL+5,7,	-l21, 0, -l23, -l24, 0, l26, l27);
		fill(evL+6,7,	alphaF*(vSq-hHydro) + Cff*(Cf-v.x) + Qs*vqstr - aspb, -alphaF*v.x + Cff, -alphaF*v.y - Qs*QStarY, -alphaF*v.z - Qs*QStarZ, alphaF, l16, l17);
#else	//self-referencing evL
		//rows
		global real* evL = eig->evL;
		fill(evL+0,7,	alphaF*(vSq-hHydro) + Cff*(Cf+v.x) - Qs*vqstr - aspb, -alphaF*v.x - Cff, -alphaF*v.y + Qs*QStarY, -alphaF*v.z + Qs*QStarZ, alphaF, AHatS*QStarY - alphaF*b.y, AHatS*QStarZ - alphaF*b.z);
		fill(evL+1,7,	.5*(v.y*betaZ - v.z*betaY), 0, -.5*betaZ, .5*betaY, 0, -.5*sqrtRho*betaZ*sbx, .5*sqrtRho*betaY*sbx);
		fill(evL+2,7,	alphaS*(vSq-hHydro) + Css*(Cs+v.x) + Qf*vqstr + afpb, -alphaS*v.x - Css, -alphaS*v.y - Qf*QStarY, -alphaS*v.z - Qf*QStarZ, alphaS, -AHatF*QStarY - alphaS*b.y, -AHatF*QStarZ - alphaS*b.z);
		fill(evL+3,7,	1. - norm*(.5*vSq - gamma_2*X/gamma_1) , norm*v.x, norm*v.y, norm*v.z, -norm, norm*b.y, norm*b.z);
		fill(evL+4,7,	alphaS*(vSq-hHydro) + Css*(Cs-v.x) - Qf*vqstr + afpb, -alphaS*v.x + Css, -alphaS*v.y + Qf*QStarY, -alphaS*v.z + Qf*QStarZ, alphaS, evL[2+7*5], evL[2+7*6]);
		fill(evL+5,7,	-evL[1+7*0], 0, -evL[1+7*2], -evL[1+7*3], 0, evL[1+7*5], evL[1+7*6]);
		fill(evL+6,7,	alphaF*(vSq-hHydro) + Cff*(Cf-v.x) + Qs*vqstr - aspb, -alphaF*v.x + Cff, -alphaF*v.y - Qs*QStarY, -alphaF*v.z - Qs*QStarZ, alphaF, evL[0+7*5], evL[0+7*6]);
#endif
	}<? end ?>
}

<? for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,2 do ?>
void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x_
) {
	<?=_7to8code(addr2,side)?>
	<?=addr1?> const real* A = eig->evL;
	for (int i = 0; i < 7; ++i) {
		real sum = 0;
		for (int j = 0; j < 7; ++j) {
			sum += A[i+7*j] * x[j];
		}
		y[i] = sum;
	}
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	<?=addr1?> const real* A = eig->evR;
	for (int i = 0; i < 7; ++i) {
		real sum = 0;
		for (int j = 0; j < 7; ++j) {
			sum += A[i+7*j] * x[j];
		}
		y[i] = sum;
	}
	<?=_8to7code(addr0, side)?>
}

<? 				if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x_
) {
	<?=_7to8code(addr2, side)?>
	<?=addr1?> const real* A = eig->A;
	for (int i = 0; i < 7; ++i) {
		real sum = 0;
		for (int j = 0; j < 7; ++j) {
			sum += A[i+7*j] * x[j];
		}
		y[i] = sum;
	}
	<?=_8to7code(addr0, side)?>
}
<?				end
			end
		end
	end
end ?>
