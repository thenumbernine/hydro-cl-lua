//// MODULE_NAME: primFromCons
//// MODULE_DEPENDS: units real3 solver_t prim_t cons_t coordLenSq

#define primFromCons(\
	/*prim_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */x\
) {\
	(result)->rho = (U)->rho;\
	(result)->v = real3_real_mul((U)->m, 1./(U)->rho);\
	(result)->B = (U)->B;\
	real const vSq = coordLenSq((result)->v, x);\
	real const BSq = coordLenSq((result)->B, x);\
	real const EKin = .5 * (U)->rho * vSq;\
	real const EMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);\
	real const EInt = (U)->ETotal - EKin - EMag;\
	(result)->P = EInt * (solver->heatCapacityRatio - 1.);\
	(result)->P = max((result)->P, (real)1e-7);\
	(result)->rho = max((result)->rho, (real)1e-7);\
	(result)->psi = (U)->psi;\
	(result)->ePot = (U)->ePot;\
}

//// MODULE_NAME: consFromPrim
//// MODULE_DEPENDS: units real3 solver_t prim_t cons_t coordLenSq

#define consFromPrim(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*prim_t const * const */W,\
	/*real3 const */x\
) {\
	(result)->rho = (W)->rho;\
	(result)->m = real3_real_mul((W)->v, (W)->rho);\
	(result)->B = (W)->B;\
	real const vSq = coordLenSq((W)->v, x);\
	real const BSq = coordLenSq((W)->B, x);\
	real const EKin = .5 * (W)->rho * vSq;\
	real const EMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);\
	real const EInt = (W)->P / (solver->heatCapacityRatio - 1.);\
	(result)->ETotal = EInt + EKin + EMag;\
	(result)->psi = (W)->psi;\
	(result)->ePot = (W)->ePot;\
}

//// MODULE_NAME: apply_dU_dW
//// MODULE_DEPENDS: units real3 solver_t prim_t cons_t

#define apply_dU_dW(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*prim_t const * const */WA, \
	/*prim_t const * const */W, \
	/*real3 const */x\
) {\
	(result)->rho = (W)->rho;\
	(result)->m = real3_add(\
		real3_real_mul((WA)->v, (W)->rho),\
		real3_real_mul((W)->v, (WA)->rho)),\
	(result)->B = (WA)->B;\
	(result)->ETotal = (W)->rho * .5 * real3_dot((WA)->v, (WA)->v)\
		+ (WA)->rho * real3_dot((W)->v, (WA)->v)\
		+ real3_dot((W)->B, (WA)->B) / (solver->mu0 / unit_kg_m_per_C2)\
		+ (W)->P / (solver->heatCapacityRatio - 1.);\
	(result)->psi = (W)->psi;\
	(result)->ePot = (W)->ePot;\
}

//// MODULE_NAME: apply_dW_dU
//// MODULE_DEPENDS: units real3 solver_t prim_t cons_t

#define apply_dW_dU(\
	/*prim_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*prim_t const * const */WA,\
	/*cons_t const * const */U,\
	/*real3 const */x\
) {\
	(result)->rho = (U)->rho;\
	(result)->v = real3_sub(\
		real3_real_mul((U)->m, 1. / (WA)->rho),\
		real3_real_mul((WA)->v, (U)->rho / (WA)->rho));\
	(result)->B = (U)->B;\
	(result)->P = (solver->heatCapacityRatio - 1.) *  (\
		.5 * (U)->rho * real3_dot((WA)->v, (WA)->v)\
		- real3_dot((U)->m, (WA)->v)\
		- real3_dot((U)->B, (WA)->B) / (solver->mu0 / unit_kg_m_per_C2)\
		+ (U)->ETotal);\
	(result)->psi = (U)->psi;\
	(result)->ePot = (U)->ePot;\
}

//// MODULE_NAME: cons_rotateFrom
//// MODULE_DEPENDS: cons_t normal_t

//TODO don't call this.
//instead separate out the vectors and rotate those vars inline, and use the orig struct for the rest.
//but this might require picking apart the individual fields of calcRoeValues
// or at least providing separate overrides of the vector components

//align from vector coordinates to the normal basis
static inline cons_t cons_rotateFrom(cons_t U, normal_t n) {
	U.m = normal_vecDotNs(n, U.m);
	U.B = normal_vecDotNs(n, U.B);
	return U;
}

//// MODULE_NAME: prim_rotateFrom
//// MODULE_DEPENDS: prim_t normal_t

//TODO also don't call this.

//align from vector coordinates to the normal basis
static inline prim_t prim_rotateFrom(prim_t W, normal_t n) {
	W.v = normal_vecDotNs(n, W.v);
	W.B = normal_vecDotNs(n, W.B);
	return W;
}
//// MODULE_NAME: cons_rotateTo
//// MODULE_DEPENDS: cons_t normal_t

//align from normal basis to vector coordinates
cons_t cons_rotateTo(cons_t U, normal_t n) {
	U.m = normal_vecFromNs(n, U.m);
	U.B = normal_vecFromNs(n, U.B);
	return U;
}

//// MODULE_NAME: applyInitCond
//// MODULE_DEPENDS: cartesianToCoord consFromPrim

kernel void applyInitCond(
	constant solver_t const * const solver,
	constant initCond_t const * const initCond,
	global cons_t * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS(0,0);
	real3 const x = cellBuf[index].pos;
	
	global cons_t * const U = UBuf + index;
	
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool const lhs = true
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
	real3 B = real3_zero;
	real ePot = 0;
	/* ignored: */
	real3 D = real3_zero;

	<?=initCode()?>
	
	prim_t W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
		.B = cartesianToCoord(B, x),
		.psi = 0,
		.ePot = ePot,
	};
	consFromPrim(UBuf + index, solver, &W, x);
}

//// MODULE_NAME: calcRoeValues
//// MODULE_DEPENDS: primFromCons roe_t

// TODO find out where mu_0 goes in the code below

//assumes UL and UR are already rotated so the 'x' direction is our flux direction
#define calcRoeValues(\
	/*roe_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */UL, \
	/*cons_t const * const */UR,\
	/*real3 const */pt\
) {\
	/*  should I use Bx, or BxL/R, for calculating the PMag at the L and R states? */\
	prim_t WL;\
	primFromCons(&WL, solver, UL, pt);\
	real const sqrtRhoL = sqrt((UL)->rho);\
	real const PMagL = .5 * coordLenSq((UL)->B, pt);\
	real const hTotalL = ((UL)->ETotal + WL.P + PMagL) / (UL)->rho;\
\
	prim_t WR;\
	primFromCons(&WR, solver, UR, pt);\
	real const sqrtRhoR = sqrt((UR)->rho);\
	real const const PMagR = .5 * coordLenSq((UR)->B, pt); real hTotalR = ((UR)->ETotal + WR.P + PMagR) / (UR)->rho;\
\
	real const dby = WL.B.y - WR.B.y;\
	real const dbz = WL.B.z - WR.B.z;\
\
	real const invDenom = 1 / (sqrtRhoL + sqrtRhoR);\
\
	(result)->rho = sqrtRhoL * sqrtRhoR;\
	(result)->v = real3_real_mul(real3_add(\
		real3_real_mul(WL.v, sqrtRhoL),\
		real3_real_mul(WR.v, sqrtRhoR)), invDenom);\
\
	(result)->hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;\
\
	(result)->B.x = (sqrtRhoL * WL.B.x + sqrtRhoR * WR.B.x) * invDenom;\
	/*  why does athena switch the weights of the By and Bz components? */\
	(result)->B.y = (sqrtRhoR * WL.B.y + sqrtRhoL * WR.B.y) * invDenom;\
	(result)->B.z = (sqrtRhoR * WL.B.z + sqrtRhoL * WR.B.z) * invDenom;\
\
	(result)->X = .5 * (dby * dby + dbz * dbz) * invDenom * invDenom;\
	(result)->Y = .5 * ((UL)->rho + (UR)->rho) / (result)->rho;\
};

//// MODULE_NAME: eigen_forRoeAvgs
//// MODULE_DEPENDS: roe_t

//assumes the vector values are x-axis aligned with the interface normal
#define eigen_forRoeAvgs(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*roe_t const * const */roe,\
	/*real3 const */pt\
) {\
	real const gamma = solver->heatCapacityRatio;\
	real const gamma_1 = gamma - 1.;\
	real const gamma_2 = gamma - 2.;\
	real const gamma_3 = gamma - 3.;\
	\
	real const rho = (roe)->rho;\
	real3 const v = (roe)->v;\
	real const hTotal = (roe)->hTotal;\
	real3 const B = (roe)->B;\
	real const X = (roe)->X;\
	real const Y = (roe)->Y;\
\
	real const _1_rho = 1. / rho;\
	real const vSq = coordLenSq(v, pt);\
/* TODO consider g_ij */\
	real const BPerpSq = B.y*B.y + B.z*B.z;\
	real const BStarPerpSq = (gamma_1 - gamma_2 * Y) * BPerpSq;\
	real const CAxSq = B.x*B.x*_1_rho;\
	real const CASq = CAxSq + BPerpSq * _1_rho;\
	(eig)->hHydro = hTotal - CASq;\
	/*  hTotal = (EHydro + EMag + P)/rho */\
	/*  (eig)->hHydro = hTotal - CASq, CASq = EMag/rho */\
	/*  (eig)->hHydro = eHydro + P/rho = eKin + eInt + P/rho */\
	/*  (eig)->hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho */\
	/*  a^2 = (gamma-1)((eig)->hHydro - eKin) = gamma P / rho */\
	(eig)->aTildeSq = max((gamma_1 * ((eig)->hHydro - .5 * vSq) - gamma_2 * X), 1e-20);\
\
	real const BStarPerpSq_rho = BStarPerpSq * _1_rho;\
	real const CATildeSq = CAxSq + BStarPerpSq_rho;\
	real const CStarSq = .5 * (CATildeSq + (eig)->aTildeSq);\
	real const CA_a_TildeSqDiff = .5 * (CATildeSq - (eig)->aTildeSq);\
	real const sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + (eig)->aTildeSq * BStarPerpSq_rho);\
	\
	(eig)->CAx = sqrt(CAxSq);\
	\
	real const CfSq = CStarSq + sqrtDiscr;\
	(eig)->Cf = sqrt(CfSq);\
	\
	real const CsSq = (eig)->aTildeSq * CAxSq / CfSq;\
	(eig)->Cs = sqrt(CsSq);\
\
/* TODO consider g_ij */\
	real const BPerpLen = sqrt(BPerpSq);\
	(eig)->BStarPerpLen = sqrt(BStarPerpSq);\
	\
	if (BPerpLen == 0) {\
		(eig)->betaY = 1;\
		(eig)->betaZ = 0;\
	} else {\
		(eig)->betaY = B.y / BPerpLen;\
		(eig)->betaZ = B.z / BPerpLen;\
	}\
	(eig)->betaStarY = (eig)->betaY / sqrt(gamma_1 - gamma_2*Y);\
	(eig)->betaStarZ = (eig)->betaZ / sqrt(gamma_1 - gamma_2*Y);\
	(eig)->betaStarSq = (eig)->betaStarY*(eig)->betaStarY + (eig)->betaStarZ*(eig)->betaStarZ;\
\
	if (CfSq - CsSq == 0) {\
		(eig)->alphaF = 1;\
		(eig)->alphaS = 0;\
	} else if ((eig)->aTildeSq - CsSq <= 0) {\
		(eig)->alphaF = 0;\
		(eig)->alphaS = 1;\
	} else if (CfSq - (eig)->aTildeSq <= 0) {\
		(eig)->alphaF = 1;\
		(eig)->alphaS = 0;\
	} else {\
		(eig)->alphaF = sqrt(((eig)->aTildeSq - CsSq) / (CfSq - CsSq));\
		(eig)->alphaS = sqrt((CfSq - (eig)->aTildeSq) / (CfSq - CsSq));\
	}\
\
	(eig)->sqrtRho = sqrt(rho);\
	real const _1_sqrtRho = 1. / (eig)->sqrtRho;\
	(eig)->sbx = B.x >= 0 ? 1 : -1;\
	real const aTilde = sqrt((eig)->aTildeSq);\
	(eig)->Qf = (eig)->Cf * (eig)->alphaF * (eig)->sbx;\
	(eig)->Qs = (eig)->Cs * (eig)->alphaS * (eig)->sbx;\
	(eig)->Af = aTilde * (eig)->alphaF * _1_sqrtRho;\
	(eig)->As = aTilde * (eig)->alphaS * _1_sqrtRho;\
\
	/* used for eigenvectors and eigenvalues */\
<? 	for _,var in ipairs(eqn.roeVars) do --\
?>	(eig)-><?=var.name?> = (roe)-><?=var.name?>;\
<?	end --\
?>\
}

//// MODULE_NAME: eqn.common
//// MODULE_DEPENDS: units coordLenSq

static inline real calc_eKin(
	prim_t const * const W,
	real3 const x
) { 
	return .5 * coordLenSq(W->v, x);
}

static inline real calc_EKin(
	prim_t const * const W,
	real3 const x
) { 
	return W->rho * calc_eKin(W, x); 
}

static inline real calc_EInt(
	constant solver_t const * const solver,
	prim_t const * const W
) { 
	return W->P / (solver->heatCapacityRatio - 1.); 
}

static inline real calc_eInt(
	constant solver_t const * const solver,
	prim_t const * const W
) { 
	return calc_EInt(solver, W) / W->rho; 
}

//units: 
//B has units kg/(C*s)
//mu0 has units kg*m/C^2
//PMag = 1/2 B^2 / mu0 has units kg/(m*s^2)
static inline real calc_EM_energy(
	constant solver_t const * const solver,
	prim_t const * const W,
	real3 const x
) { 
	return .5 * coordLenSq(W->B, x) / (solver->mu0 / unit_kg_m_per_C2);
}

//same as calc_EM_energy
static inline real calc_PMag(
	constant solver_t const * const solver,
	prim_t const * const W,
	real3 const x
) { 
	return .5 * coordLenSq(W->B, x) / (solver->mu0 / unit_kg_m_per_C2);
}

static inline real calc_EHydro(
	constant solver_t const * const solver,
	prim_t const * const W,
	real3 const x
) { 
	return calc_EKin(W, x) + calc_EInt(solver, W); 
}

static inline real calc_eHydro(
	constant solver_t const * const solver,
	prim_t const * const W,
	real3 const x
) { 
	return calc_EHydro(solver, W, x) / W->rho; 
}

static inline real calc_ETotal(
	constant solver_t const * const solver,
	prim_t const * const W,
	real3 x
) { 
	return calc_EKin(W, x) + calc_EInt(solver, W) + calc_EM_energy(solver, W, x); 
}

static inline real calc_eTotal(
	constant solver_t const * const solver,
	prim_t const * const W,
	real3 const x
) {
	return calc_ETotal(solver, W, x) / W->rho; 
}

static inline real calc_H(
	constant solver_t const * const solver,
	real const P
) { 
	return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); 
}

static inline real calc_h(
	constant solver_t const * const solver,
	real const rho,
	real const P
) { 
	return calc_H(solver, P) / rho; 
}

static inline real calc_HTotal(
	constant solver_t const * const solver,
	prim_t const * const W,
	real const ETotal,
	real3 const x
) { 
	return W->P + calc_PMag(solver, W, x) + ETotal; 
}

static inline real calc_hTotal(
	constant solver_t const * const solver,
	prim_t const * const W,
	real const ETotal,
	real3 const x
) { 
	return calc_HTotal(solver, W, ETotal, x) / W->rho; 
}

//notice, this is speed of sound, to match the name convention of hydro/eqn/euler
//but Cs in eigen_t is the slow speed
//most the MHD papers use 'a' for the speed of sound
static inline real calc_Cs(
	constant solver_t const * const solver,
	prim_t const * const W
) { 
	return sqrt(solver->heatCapacityRatio * W->P / W->rho);
}

//CA = B/sqrt(mu0 rho)
//B has units kg/(C*s)
//mu0 has units kg*m/C^2
//rho has units kg/m^3
//CA has units m/s
static inline real3 calc_CA(
	constant solver_t const * const solver,
	cons_t const * const U
) {
	return real3_real_mul(U->B, 1./sqrt(U->rho * solver->mu0 / unit_kg_m_per_C2));
}

//// MODULE_NAME: fluxFromCons
//// MODULE_DEPENDS: units solver_t cons_t prim_t primFromCons normal_t coordLenSq

#define fluxFromCons(\
	/*cons_t * const */F,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	prim_t W;\
	primFromCons(&W, solver, U, pt);\
	real vj = normal_vecDotN1(n, W.v);\
	real Bj = normal_vecDotN1(n, W.B);\
	real BSq = coordLenSq(W.B, pt);\
	real BDotV = real3_dot(W.B, W.v);\
	real PMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);\
	real PTotal = W.P + PMag;\
	real HTotal = (U)->ETotal + PTotal;\
\
	(F)->rho = normal_vecDotN1(n, (U)->m);\
	(F)->m = real3_sub(real3_real_mul((U)->m, vj), real3_real_mul((U)->B, Bj / (solver->mu0 / unit_kg_m_per_C2)));\
	(F)->m.x += PTotal * normal_u1x(n);\
	(F)->m.y += PTotal * normal_u1y(n);\
	(F)->m.z += PTotal * normal_u1z(n);\
	(F)->B = real3_sub(real3_real_mul((U)->B, vj), real3_real_mul(W.v, Bj));\
	(F)->ETotal = HTotal * vj - BDotV * Bj / (solver->mu0 / unit_kg_m_per_C2);\
	(F)->psi = 0;\
	(F)->ePot = 0;\
}

//// MODULE_NAME: calcCellMinMaxEigenvalues
//// MODULE_DEPENDS: range_t cons_rotateFrom primFromCons

//called from calcDT
#define calcCellMinMaxEigenvalues(\
	/*range_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	cons_t U_ = cons_rotateFrom(*U, n);\
	prim_t W;\
	primFromCons(&W, solver, &U_, pt);\
\
<? if false then ?>\
	prim_t W;\
	primFromCons(&W, solver, U, pt);\
	real3 v = W.v;\
	real3 B = W.B;\
\
	real BSq = coordLenSq(B, pt);\
	real invRho = 1./W.rho;\
\
	real aSq = solver->heatCapacityRatio * W.P * invRho;\
	real B_n = normal_vecDotN1(n, B);\
	real CaxSq = B_n * B_n * invRho;\
	real CaSq = BSq * invRho;\
\
	real CStarSq = .5 * (CaSq + aSq);\
	real sqrtCfsDiscr = sqrt(max(0., CStarSq * CStarSq - aSq * CaxSq));\
\
	real CfSq = CStarSq + sqrtCfsDiscr;\
	real CsSq = CStarSq - sqrtCfsDiscr;\
\
	real Cf = sqrt(CfSq);\
	real Cs = sqrt(max(CsSq, 0.));\
	real v_n = normal_vecDotN1(n, v);\
	return (range_t){.min=v_n - Cf, .max=v_n + Cf};\
<? else ?>\
	const real gamma = solver->heatCapacityRatio;\
	const real gamma_1 = gamma - 1.;\
	const real gamma_2 = gamma - 2.;\
\
	real const rho = W.rho;\
	real3 const v = W.v;\
	real3 const B = W.B;\
	real const hTotal = .5 * coordLenSq(W.v, pt) + (W.P * gamma / gamma_1 + coordLenSq(B, pt)) / W.rho;\
\
	/* the rest of this matches calcEigenBasis: */\
\
	real const _1_rho = 1. / rho;\
	real const vSq = coordLenSq(v, pt);\
/* TODO consider g_ij */\
	real const BPerpSq = B.y*B.y + B.z*B.z;\
	real const BStarPerpSq = (gamma_1 - gamma_2) * BPerpSq;\
	real const CAxSq = B.x*B.x*_1_rho;\
	real const CASq = CAxSq + BPerpSq * _1_rho;\
	real const hHydro = hTotal - CASq;\
	/*  hTotal = (EHydro + EMag + P)/rho */\
	/*  hHydro = hTotal - CASq, CASq = EMag/rho */\
	/*  hHydro = eHydro + P/rho = eKin + eInt + P/rho */\
	/*  hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho */\
	/*  a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho */\
	real const aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2), 1e-20);\
\
	real const BStarPerpSq_rho = BStarPerpSq * _1_rho;\
	real const CATildeSq = CAxSq + BStarPerpSq_rho;\
	real const CStarSq = .5 * (CATildeSq + aTildeSq);\
	real const CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);\
	real const sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * BStarPerpSq_rho);\
\
	real const CfSq = CStarSq + sqrtDiscr;\
	real const Cf = sqrt(CfSq);\
\
	real const CsSq = aTildeSq * CAxSq / CfSq;\
	real const Cs = sqrt(CsSq);\
\
	real const lambdaFastMin = v.x - Cf;\
	real const lambdaFastMax = v.x + Cf;\
\
	(result)->min = lambdaFastMin;\
	(result)->max = lambdaFastMax;\
<? end ?>\
}

//// MODULE_NAME: eigen_forInterface
//// MODULE_DEPENDS: roe_t cons_rotateFrom calcRoeValues eigen_forRoeAvgs

#define eigen_forInterface(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */UL,\
	/*cons_t const * const */UR,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	/* rotate UL and UR to be x-aligned?  that takes the normal ... */\
	/* swap the sides with x here, so all the fluxes are in the 'x' direction */\
	cons_t UL_ = cons_rotateFrom(*UL, n);\
	cons_t UR_ = cons_rotateFrom(*UR, n);\
	roe_t roe;\
	calcRoeValues(&roe, solver, &UL_, &UR_, x);\
	eigen_forRoeAvgs(eig, solver, &roe, x);\
}

//// MODULE_NAME: eigen_left/rightTransform
//// MODULE_DEPENDS: eigen_t waves_t cons_rotateFrom cons_rotateTo

#define eigen_leftTransform(\
	/*waves_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*cons_t const * const */inputU_,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	cons_t const inputU = cons_rotateFrom(*(inputU_), n);\
	real const gamma = solver->heatCapacityRatio;\
	real const gamma_1 = gamma - 1.;\
	real const gamma_2 = gamma - 2.;\
\
<? for _,var in ipairs(eqn.eigenVars) do --\
?> 	<?=var.type or 'real'?> <?=var.name?> = (eig)-><?=var.name?>;\
<? end ?>\
\
<? print("you can't use coordLenSq (which uses g_ij) after rotating coordinates") ?>\
	real const vSq = coordLenSq(v, pt);\
\
	/*  left eigenvectors */\
	real norm = .5 / (eig)->aTildeSq;\
	real const Cff = norm * alphaF * Cf;\
	real const Css = norm * alphaS * Cs;\
	Qf = Qf * norm;\
	Qs = Qs * norm;\
	real const AHatF = norm * Af * rho;\
	real const AHatS = norm * As * rho;\
	real const afpb = norm * Af * BStarPerpLen;\
	real const aspb = norm * As * BStarPerpLen;\
\
	norm = norm * gamma_1;\
	alphaF = alphaF * norm;\
	alphaS = alphaS * norm;\
	real const QStarY = betaStarY/betaStarSq;\
	real const QStarZ = betaStarZ/betaStarSq;\
	real const vqstr = (v.y * QStarY + v.z * QStarZ);\
	norm = norm * 2.;\
\
	real const l16 = AHatS * QStarY - alphaF * B.y;\
	real const l17 = AHatS * QStarZ - alphaF * B.z;\
	real const l21 = .5 * (v.y * betaZ - v.z * betaY);\
	real const l23 = .5 * betaZ;\
	real const l24 = .5 * betaY;\
	real const l26 = -.5 * sqrtRho * betaZ * sbx;\
	real const l27 = .5 * sqrtRho * betaY * sbx;\
	real const l36 = -AHatF * QStarY - alphaS * B.y;\
	real const l37 = -AHatF * QStarZ - alphaS * B.z;\
\
	(result)->ptr[0] = \
		  inputU.rho * (alphaF * (vSq - (eig)->hHydro) + Cff * (Cf + v.x) - Qs * vqstr - aspb) \
		+ inputU.m.x * (-alphaF * v.x - Cff)\
		+ inputU.m.y * (-alphaF * v.y + Qs * QStarY)\
		+ inputU.m.z * (-alphaF * v.z + Qs * QStarZ)\
		+ inputU.ETotal * alphaF\
		+ inputU.B.y * l16\
		+ inputU.B.z * l17;\
	(result)->ptr[1] = \
		  inputU.rho * l21\
		+ inputU.m.y * l23\
		+ inputU.m.z * l24\
		+ inputU.B.y * l26\
		+ inputU.B.z * l27;\
	(result)->ptr[2] = \
		  inputU.rho * (alphaS * (vSq - (eig)->hHydro) + Css * (Cs + v.x) + Qf * vqstr + afpb)\
		+ inputU.m.x * (-alphaS * v.x - Css)\
		+ inputU.m.y * (-alphaS * v.y - Qf * QStarY)\
		+ inputU.m.z * (-alphaS * v.z - Qf * QStarZ)\
		+ inputU.ETotal * alphaS\
		+ inputU.B.y * l36\
		+ inputU.B.z * l37;\
	(result)->ptr[3] = \
		  inputU.rho * (1. - norm * (.5 * vSq - gamma_2 * X / gamma_1))\
		+ inputU.m.x * norm*v.x\
		+ inputU.m.y * norm*v.y\
		+ inputU.m.z * norm*v.z\
		+ inputU.ETotal * -norm\
		+ inputU.B.y * norm*B.y\
		+ inputU.B.z * norm*B.z;\
	(result)->ptr[4] = \
		  inputU.rho * (alphaS * (vSq - (eig)->hHydro) + Css * (Cs - v.x) - Qf * vqstr + afpb)\
		+ inputU.m.x * (-alphaS * v.x + Css)\
		+ inputU.m.y * (-alphaS * v.y + Qf * QStarY)\
		+ inputU.m.z * (-alphaS * v.z + Qf * QStarZ)\
		+ inputU.ETotal * alphaS\
		+ inputU.B.y * l36\
		+ inputU.B.z * l37;\
	(result)->ptr[5] = \
		  inputU.rho * -l21\
		+ inputU.m.y * -l23\
		+ inputU.m.z * -l24\
		+ inputU.B.y * l26\
		+ inputU.B.z * l27;\
	(result)->ptr[6] = \
		  inputU.rho * (alphaF * (vSq - (eig)->hHydro) + Cff * (Cf - v.x) + Qs * vqstr - aspb)\
		+ inputU.m.x * (-alphaF * v.x + Cff)\
		+ inputU.m.y * (-alphaF * v.y - Qs * QStarY)\
		+ inputU.m.z * (-alphaF * v.z - Qs * QStarZ)\
		+ inputU.ETotal * alphaF\
		+ inputU.B.y * l16\
		+ inputU.B.z * l17;\
}

#define eigen_rightTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*waves_t const * const */input,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real const gamma = solver->heatCapacityRatio;\
	real const gamma_1 = gamma - 1.;\
	real const gamma_2 = gamma - 2.;\
\
<? for _,var in ipairs(eqn.eigenVars) do --\
?> 	<?=var.type or "real"?> <?=var.name?> = (eig)-><?=var.name?>;\
<? end ?>\
\
	real const vSq = coordLenSq(v, pt);\
	real const vDotBeta = v.y*betaStarY + v.z*betaStarZ;\
	real const _1_sqrtRho = 1. / (eig)->sqrtRho;\
	real const Afpbb = Af * BStarPerpLen * betaStarSq;\
	real const Aspbb = As * BStarPerpLen * betaStarSq;\
\
	real const lambdaFastMin = (eig)->v.x - (eig)->Cf;\
	real const lambdaSlowMin = (eig)->v.x - (eig)->Cs;\
	real const lambdaSlowMax = (eig)->v.x + (eig)->Cs;\
	real const lambdaFastMax = (eig)->v.x + (eig)->Cf;\
\
	/*  right eigenvectors */\
	real const qa3 = alphaF * v.y;\
	real const qb3 = alphaS * v.y;\
	real const qc3 = Qs * betaStarY;\
	real const qd3 = Qf * betaStarY;\
	real const qa4 = alphaF * v.z;\
	real const qb4 = alphaS * v.z;\
	real const qc4 = Qs * betaStarZ;\
	real const qd4 = Qf * betaStarZ;\
	real const r52 = -(v.y * betaZ - v.z * betaY);\
	real const r61 = As * betaStarY;\
	real const r62 = -betaZ * sbx * _1_sqrtRho;\
	real const r63 = -Af * betaStarY;\
	real const r71 = As * betaStarZ;\
	real const r72 = betaY * sbx * _1_sqrtRho;\
	real const r73 = -Af * betaStarZ;\
\
	cons_t resultU;\
	resultU.rho =\
		  (input)->ptr[0] * alphaF\
		+ (input)->ptr[2] * alphaS\
		+ (input)->ptr[3]\
		+ (input)->ptr[4] * alphaS\
		+ (input)->ptr[6] * alphaF;\
	resultU.m.x =\
		  (input)->ptr[0] * alphaF * lambdaFastMin\
		+ (input)->ptr[2] * alphaS * lambdaSlowMin\
		+ (input)->ptr[3] * v.x\
		+ (input)->ptr[4] * alphaS * lambdaSlowMax\
		+ (input)->ptr[6] * alphaF * lambdaFastMax;\
	resultU.m.y =\
		  (input)->ptr[0] * (qa3 + qc3)\
		+ (input)->ptr[1] * -betaZ\
		+ (input)->ptr[2] * (qb3 - qd3)\
		+ (input)->ptr[3] * v.y\
		+ (input)->ptr[4] * (qb3 + qd3)\
		+ (input)->ptr[5] * betaZ\
		+ (input)->ptr[6] * (qa3 - qc3);\
	resultU.m.z =\
		  (input)->ptr[0] * (qa4 + qc4)\
		+ (input)->ptr[1] * betaY\
		+ (input)->ptr[2] * (qb4 - qd4)\
		+ (input)->ptr[3] * v.z\
		+ (input)->ptr[4] * (qb4 + qd4)\
		+ (input)->ptr[5] * -betaY\
		+ (input)->ptr[6] * (qa4 - qc4);\
	resultU.ETotal =\
		  (input)->ptr[0] * (alphaF*((eig)->hHydro - v.x*Cf) + Qs*vDotBeta + Aspbb)\
		+ (input)->ptr[1] * r52\
		+ (input)->ptr[2] * (alphaS*((eig)->hHydro - v.x*Cs) - Qf*vDotBeta - Afpbb)\
		+ (input)->ptr[3] * (.5*vSq + gamma_2*X/gamma_1)\
		+ (input)->ptr[4] * (alphaS*((eig)->hHydro + v.x*Cs) + Qf*vDotBeta - Afpbb)\
		+ (input)->ptr[5] * -r52\
		+ (input)->ptr[6] * (alphaF*((eig)->hHydro + v.x*Cf) - Qs*vDotBeta + Aspbb);\
	resultU.B.x = 0;\
	resultU.B.y =\
		  (input)->ptr[0] * r61\
		+ (input)->ptr[1] * r62\
		+ (input)->ptr[2] * r63\
		+ (input)->ptr[4] * r63\
		+ (input)->ptr[5] * r62\
		+ (input)->ptr[6] * r61;\
	resultU.B.z =\
		  (input)->ptr[0] * r71\
		+ (input)->ptr[1] * r72\
		+ (input)->ptr[2] * r73\
		+ (input)->ptr[4] * r73\
		+ (input)->ptr[5] * r72\
		+ (input)->ptr[6] * r71;\
	resultU.psi = 0;\
	*(result) = cons_rotateTo(resultU, n);\
}

//// MODULE_NAME: eigen_fluxTransform

#define eigen_fluxTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*cons_t const * const */inputU_,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	cons_t inputU = cons_rotateFrom(*(inputU_), n);\
\
	real const gamma = solver->heatCapacityRatio;\
	real const gamma_1 = gamma - 1.;\
	real const gamma_2 = gamma - 2.;\
	real const gamma_3 = gamma - 3.;\
\
<? for _,var in ipairs(eqn.eigenVars) do --\
?> 	<?=var.type or "real"?> <?=var.name?> = (eig)-><?=var.name?>;\
<? end ?>\
\
	real const _1_rho = 1. / rho;\
#warning you can't use coordLenSq (which uses g_ij) after rotating coordinates \
	real const vSq = coordLenSq(v, pt);\
	real const BDotV = real3_dot(B,v);\
\
	/*  dF/dU */\
	cons_t resultU;\
	resultU.rho = inputU.m.x;\
	resultU.m.x =\
		  inputU.rho * (-v.x*v.x + .5*gamma_1*vSq - gamma_2*X)\
		+ inputU.m.x * -gamma_3*v.x\
		+ inputU.m.y * -gamma_1*v.y\
		+ inputU.m.z * -gamma_1*v.z\
		+ inputU.ETotal * gamma_1\
		+ inputU.B.y * -gamma_2*Y*B.y\
		+ inputU.B.z * -gamma_2*Y*B.z;\
	resultU.m.y =\
		  inputU.rho * -v.x*v.y\
		+ inputU.m.x * v.y\
		+ inputU.m.y * v.x\
		+ inputU.B.y * -B.x;\
	resultU.m.z =\
		  inputU.rho * -v.x*v.z\
		+ inputU.m.x * v.z\
		+ inputU.m.z * v.x\
		+ inputU.B.z * -B.x;\
	resultU.ETotal =\
		  inputU.rho * (v.x*(.5*gamma_1*vSq - hTotal) + B.x*BDotV * _1_rho)\
		+ inputU.m.x * (-gamma_1*v.x*v.x + hTotal - B.x*B.x * _1_rho)\
		+ inputU.m.y * (-gamma_1*v.x*v.y - B.x*B.y * _1_rho)\
		+ inputU.m.z * (-gamma_1*v.x*v.z - B.x*B.z * _1_rho)\
		+ inputU.ETotal * gamma*v.x\
		+ inputU.B.y * (-gamma_2*Y*B.y*v.x - B.x*v.y)\
		+ inputU.B.z * (-gamma_2*Y*B.z*v.x - B.x*v.z);\
	resultU.B.x = 0;\
	resultU.B.y =\
		  inputU.rho * (B.x*v.y - B.y*v.x) * _1_rho\
		+ inputU.m.x * B.y * _1_rho\
		+ inputU.m.y * -B.x * _1_rho\
		+ inputU.B.y * v.x;\
	resultU.B.z =\
		  inputU.rho * (B.x*v.z - B.z*v.x) * _1_rho\
		+ inputU.m.x * B.z * _1_rho\
		+ inputU.m.z * -B.x * _1_rho\
		+ inputU.B.z * v.x;\
	resultU.psi = 0;\
	*(result) = cons_rotateTo(resultU, n);\
}

//// MODULE_NAME: eigen_forCell
//// MODULE_DEPENDS: roe_t primFromCons

#define eigen_forCell(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	prim_t W;\
	primFromCons(&W, solver, U, pt);\
	real PMag = .5 * coordLenSq(W.B, pt);\
	real hTotal = ((U)->ETotal + W.P + PMag) / W.rho;\
	roe_t roe = {\
		.rho = W.rho,\
		.v = W.v,\
		.hTotal = hTotal,\
		.B = W.B,\
		.X = 0,\
		.Y = 1,\
	};\
	eigen_forRoeAvgs(eig, solver, &roe, pt);\
}

//// MODULE_NAME: addSource
//// MODULE_DEPENDS: units primFromCons

kernel void addSource(
	constant solver_t const * const solver,
	global cons_t * const derivBuf,
	global cons_t const * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 const x = cellBuf[index].pos;
	global cons_t * const deriv = derivBuf + index;
	global cons_t const * const U = UBuf + index;

<? if not require 'hydro.coord.cartesian'.is(solver.coord) then ?>
	prim_t W;
	primFromCons(&W, solver, U, x);
	real const BSq = coordLenSq(U->B, x);
	real const PMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);
	real const PTotal = W.P + PMag;
	real3 const m_conn_vv = coord_conn_apply23(W.v, U->m, x);
	deriv->m = real3_sub(deriv->m, m_conn_vv);	/* -Conn^i_jk rho v^j v^k  */
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_raise(coord_conn_trace13(x), x), PTotal));		/* +Conn^j_kj g^ki PTotal */
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_conn_apply23(U->B, U->B, x), 1. / (solver->mu0 / unit_kg_m_per_C2)));	/* + 1/mu0 Conn^i_jk B^j B^k */
<? end ?>
}

//// MODULE_NAME: constrainU
//// MODULE_DEPENDS: consFromPrim primFromCons

kernel void constrainU(
	constant solver_t const * const solver,
	global cons_t * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS(0,0);
	real3 const x = cellBuf[index].pos;
	
	global cons_t * const U = UBuf + index;
	prim_t W;
	primFromCons(&W, solver, U, x);

	W.rho = max(W.rho, 1e-7);
	W.P = max(W.P, 1e-7);

	consFromPrim(U, solver, &W, x);
}
