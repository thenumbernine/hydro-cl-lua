//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: <?=prim_t?> <?=cons_t?> coordLenSq

#define <?=primFromCons?>(\
	/*<?=prim_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
	(result)->rhoBar = (U)->rhoBar;\
	(result)->k = (U)->rhoBar_k / (U)->rhoBar;\
	(result)->omega = (U)->rhoBar_omega / (U)->rhoBar;\
	(result)->ePot = (U)->ePot;\
	(result)->vTilde = real3_real_mul((U)->rhoBar_vTilde, 1. / (U)->rhoBar);\
	\
	real const vTildeSq = coordLenSq((result)->vTilde, pt);\
	real const rhoBar_eIntTilde = (U)->rhoBar_eTotalTilde - .5 * (U)->rhoBar * vTildeSq - (U)->rhoBar_k;\
	real const rhoBar_TTilde = rhoBar_eIntTilde / solver->C_v;\
	real const PBar = rhoBar_TTilde * solver->gasConstant;\
	(result)->PStar = PBar + 2./3. * (U)->rhoBar_k;\
}

//// MODULE_NAME: <?=consFromPrim?>
//// MODULE_DEPENDS: <?=prim_t?> <?=cons_t?> coordLenSq

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
	real const rhoBar_k = (W)->rhoBar * (W)->k;\
\
	/* eqn 6: PStar = PBar + 2/3 rhoBar k */\
	real const PBar = (W)->PStar - 2./3. * rhoBar_k;\
\
	/* eqn 10: PBar = rhoBar R TTilde */\
	real const TTilde = PBar / ((W)->rhoBar * solver->gasConstant);\
	\
	/* eqn 6: eIntTilde = C_v TTilde */\
	real const eIntTilde = solver->C_v * TTilde;\
	\
	/* eqn 6: eTotalTilde = eIntTilde + 1/2 vTilde^2 + (W)->k */\
	/* so eTotalTilde = C_v PStar / rhoBar + 1/2 vTilde^2 + (1 - 2/3 C_v / solver->gasConstant) k */\
	real const eTotalTilde = eIntTilde + .5 * coordLenSq((W)->vTilde, pt) + (W)->k;\
	\
	(result)->rhoBar = (W)->rhoBar;\
	(result)->rhoBar_vTilde = real3_real_mul((W)->vTilde, (W)->rhoBar);\
	(result)->rhoBar_eTotalTilde = (W)->rhoBar * eTotalTilde;\
	(result)->rhoBar_k = rhoBar_k;\
	(result)->rhoBar_omega = (W)->rhoBar * (W)->omega;\
	(result)->ePot = (W)->ePot;\
}

//// MODULE_NAME: <?=apply_dU_dW?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> coord_lower

#define <?=apply_dU_dW?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA, \
	/*<?=prim_t?> const * const */W, \
	/*real3 const */pt\
) {\
	real3 WA_vTildeL = coord_lower((WA)->vTilde, pt);\
	(result)->rhoBar = (W)->rhoBar;\
	(result)->rhoBar_vTilde = real3_add(\
		real3_real_mul((WA)->vTilde, (W)->rhoBar), \
		real3_real_mul((W)->vTilde, (WA)->rhoBar));\
	(result)->rhoBar_eTotalTilde = (W)->rhoBar * (.5 * real3_dot((WA)->vTilde, WA_vTildeL) \
			+ (1. - 2./3. * C_v_over_R) * (WA)->k)\
		+ (WA)->rhoBar * real3_dot((W)->vTilde, WA_vTildeL)\
		+ (W)->PStar * C_v_over_R\
		+ (1. - 2./3 * C_v_over_R) * (WA)->rhoBar * (W)->k;\
	(result)->rhoBar_k = (WA)->k * (W)->rhoBar + (WA)->rhoBar * (W)->k;\
	(result)->rhoBar_omega = (WA)->omega * (W)->rhoBar + (WA)->rhoBar * (W)->omega;\
	(result)->ePot = (W)->ePot;\
}

//// MODULE_NAME: <?=apply_dW_dU?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> coord_lower

#define <?=apply_dW_dU?>(\
	/*<?=prim_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
	real3 const WA_vTildeL = coord_lower((WA)->vTilde, pt);\
	(result)->rhoBar = (U)->rhoBar;\
	(result)->vTilde = real3_sub(\
		real3_real_mul((U)->rhoBar_vTilde, 1. / (WA)->rhoBar),\
		real3_real_mul((WA)->vTilde, (U)->rhoBar / (WA)->rhoBar));\
	(result)->PStar = R_over_C_v * (\
			.5 * real3_dot((WA)->vTilde, WA_vTildeL) * (U)->rhoBar \
			- real3_dot((U)->rhoBar_vTilde, WA_vTildeL)\
			+ (U)->rhoBar_eTotalTilde\
		) + (2./3. * R_over_C_v - 1.) * (U)->rhoBar_k;\
	(result)->k = (U)->rhoBar_k / (WA)->rhoBar - (WA)->k / (WA)->rhoBar * (U)->rhoBar;\
	(result)->omega = (U)->rhoBar_omega / (WA)->rhoBar - (WA)->omega / (WA)->rhoBar * (U)->rhoBar;\
	(result)->ePot = (U)->ePot;\
}

//// MODULE_NAME: eqn.common
//// MODULE_DEPENDS: coordLenSq

#define R_over_C_v (solver->gasConstant / solver->C_v)
#define C_v_over_R (solver->C_v / solver->gasConstant)

//real calc_H(real PStar) { return PStar * ((R_over_C_v + 1.) / (R_over_C_v)); }
//real calc_h(real rhoBar, real PStar) { return calc_H(PStar) / rhoBar; }
//real calc_hTotal(real rhoBar, real PStar, real rhoBar_eTotalTilde) { return (PStar + rhoBar_eTotalTilde) / rhoBar; }
//real calc_HTotal(real PStar, real rhoBar_eTotalTilde) { return PStar + rhoBar_eTotalTilde; }
static inline real calc_eKinTilde(<?=prim_t?> const * const W, real3 x) { return .5 * coordLenSq((W)->vTilde, x); }
static inline real calc_EKinTilde(<?=prim_t?> const * const W, real3 x) { return (W)->rhoBar * calc_eKinTilde(W, x); }

//before
//real calc_EIntTilde(<?=prim_t?> W) { return W.PStar * C_v_over_R; }
//real calc_eIntTilde(<?=prim_t?> W) { return calc_EIntTilde(W) / W.rhoBar; }

//after
static inline real calc_PBar(<?=prim_t?> const * const W) { return (W)->PStar - 2./3. * (W)->rhoBar * (W)->k; }
static inline real calc_TTilde(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) { return calc_PBar(W) / ((W)->rhoBar * solver->gasConstant); }
static inline real calc_eIntTilde(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) { return solver->C_v * calc_TTilde(solver, W); }
static inline real calc_EIntTilde(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) { return (W)->rhoBar * calc_eIntTilde(solver, W); }

static inline real calc_EKin_fromCons(<?=cons_t?> const * const U, real3 const x) { return .5 * coordLenSq((U)->rhoBar_vTilde, x) / (U)->rhoBar; }
static inline real calc_ETotal(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W, real3 const x) {
	return calc_EKinTilde(W, x) + calc_EIntTilde(solver, W);
}

static inline real calc_Cs(constant <?=solver_t?> const * const solver, <?=prim_t?> const * const W) {
	return sqrt((R_over_C_v + 1.) * (W)->PStar / (W)->rhoBar);
}

//// MODULE_NAME: <?=applyInitCond?>
//// MODULE_DEPENDS: cartesianToCoord <?=consFromPrim?>

kernel void <?=applyInitCond?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	real3 const x = cellBuf[index].pos;
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
	
	/* TODO make this B for Maxwell */
	
	real3 B = real3_zero;	/* set for MHD / thrown away for pure NavierStokesWilcox */
	real ePot = 0;

	<?=initCode()?>

	<?=prim_t?> W = {
		.rhoBar = rho,
		.vTilde = cartesianToCoord(v, x),	/* transform from cartesian to coordinate space  */
		.PStar = P,
		.k = 0,
		.omega = 0,
		.ePot = ePot,
	};
	<?=consFromPrim?>(UBuf + index, solver, &W, x);
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: coord_g_uu## <?=primFromCons?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, x);\
	real const vTilde_j = W.vTilde.s[n.side];\
	\
	/* this is the flux term used, but is it technically called 'HTotal' ? */\
	/* 'HTotal' = rhoBar (1/2 vTilde^2 + (1 - 2/3/gamma_1) k) + 1/gamma_1 PStar ... + PStar */\
	/* 'HTotal' = rhoBar (1/2 vTilde^2 + (1 - 2/3 C_v/R) k) + (1 + C_v/R) PStar */\
	/* 'hTotal' = 1/2 vTilde^2 + (1 - 2/3 C_v/R) k + (1 + C_v/R) PStar / rhoBar */\
	real const HTotal = (U)->rhoBar_eTotalTilde + W.PStar;\
	\
	(resultFlux)->rhoBar = (U)->rhoBar_vTilde.s[n.side];\
	(resultFlux)->rhoBar_vTilde = real3_real_mul((U)->rhoBar_vTilde, vTilde_j);\
	\
	if (false) {}\
<? for side=0,2 do --\
?>	else if (n.side == <?=side?>) {\
<? for i=0,2 do --\
?>	(resultFlux)->rhoBar_vTilde.s<?=i?> += coord_g_uu<?=i?><?=side?>(x) * W.PStar;\
<? end --\
?>	}\
<? end --\
?>	(resultFlux)->rhoBar_eTotalTilde = HTotal * vTilde_j;\
	(resultFlux)->rhoBar_k = 0;\
	(resultFlux)->rhoBar_omega = 0;\
	(resultFlux)->ePot = 0;\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: coord_sqrt_g_uu## <?=primFromCons?> eqn.common

#define <?=calcCellMinMaxEigenvalues?>(\
	/*range_t * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, x);\
	real Cs = calc_Cs(solver, W);\
	if (false) {}\
<? for side=0,2 do ?>\
	else if (n.side == <?=side?>) {\
		real Cs_sqrt_gU = Cs * coord_sqrt_g_uu<?=side..side?>(x);\
		(result)->min = W.vTilde.s<?=side?> - Cs_sqrt_gU,\
		(result)->max = W.vTilde.s<?=side?> + Cs_sqrt_gU,\
	}\
<? end ?>\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=primFromCons?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	<?=prim_t?> WL;\
	<?=primFromCons?>(&WL, solver, UL, x);\
	real const sqrtRhoL = sqrt(WL.rhoBar);\
	real3 const vTildeL = WL.vTilde;\
	real const hTotalL = (WL.PStar + (UL)->rhoBar_eTotalTilde) / (UL)->rhoBar - (UL)->ePot;\
	real const kL = WL.k;\
	real const omegaL = WL.omega;\
\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, x);\
	real const sqrtRhoR = sqrt(WR.rhoBar);\
	real3 const vTildeR = WR.vTilde;\
	real const hTotalR = (WR.PStar + (UR)->rhoBar_eTotalTilde) / (UR)->rhoBar - (UR)->ePot;\
	real const kR = WR.k;\
	real const omegaR = WR.omega;\
\
	real const invDenom = 1./(sqrtRhoL + sqrtRhoR);\
	\
	/* Roe-averaged */\
	(result)->rhoBar = sqrtRhoL * sqrtRhoR;\
	(result)->vTilde = real3_add(\
		real3_real_mul(vTildeL, sqrtRhoL * invDenom),\
		real3_real_mul(vTildeR, sqrtRhoR * invDenom));\
	(result)->hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);\
	/* I'm making this part up: */\
	(result)->k = invDenom * (sqrtRhoL * kL + sqrtRhoR * kR);\
	(result)->omega = invDenom * (sqrtRhoL * omegaL + sqrtRhoR * omegaR);\
\
	/* derived: */\
	(result)->vTildeSq = coordLenSq((result)->vTilde, x);\
	\
/*	\
rhoBar eIntTilde = rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot\
rhoBar TTilde C_v = rhoBar eIntTilde\
rhoBar TTilde = 1/C_v ( rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot )\
PBar / R = rhoBar TTilde\
PBar = R/C_v ( rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot )\
PBar = PStar - 2/3 rhoBar k\
PStar - 2/3 rhoBar k = R/C_v ( rhoBar eTotalTilde - 1/2 rhoBar vTilde^2 - rhoBar k - rhoBar ePot )\
\
rhoBar eTotalTilde = rhoBar (hTotalTilde + ePot) - PStar\
(R/C_v + 1) PStar / rhoBar = R/C_v (hTotalTilde - 1/2 vTilde^2 - (1 - 2/3 C_v/R) k)\
\
Cs^2 = (R/C_v+1) PStar / rhoBar\
Cs^2 = R/C_v (hTotalTilde - 1/2 vTilde^2 - (1 - 2/3 C_v/R) k)\
*/	\
	real const eKinTilde = .5 * (result)->vTildeSq;\
	real const CsSq = R_over_C_v * ((result)->hTotal - eKinTilde - (result)->k) + 2./3. * (result)->k;\
	(result)->Cs = sqrt(CsSq);\
}

<?
local prefixes = {}
for side=0,2 do 
	local prefix	
	if side == 0 then
		prefix = [[
	real const nx = 1, ny = 0, nz = 0;
	real const n1x = 0, n1y = 1, n1z = 0;
	real const n2x = 0, n2y = 0, n2z = 1;
	real const vTilde_n = vTilde.x, vTilde_n1 = vTilde.y, vTilde_n2 = vTilde.z;
]] 
	elseif side == 1 then
		prefix = [[
	real const nx = 0, ny = 1, nz = 0;
	real const n1x = 0, n1y = 0, n1z = 1;
	real const n2x = 1, n2y = 0, n2z = 0;
	real const vTilde_n = vTilde.y, vTilde_n1 = vTilde.z, vTilde_n2 = vTilde.x;
]] 
	elseif side == 2 then
		prefix = [[
	real const nx = 0, ny = 0, nz = 1;
	real const n1x = 1, n1y = 0, n1z = 0;
	real const n2x = 0, n2y = 1, n2z = 0;
	real const vTilde_n = vTilde.z, vTilde_n1 = vTilde.x, vTilde_n2 = vTilde.y;
]]
	end
	
	prefix = [[
	sym3 const gU = coord_g_uu(pt);
	real const gUjj = gU.s]]..side..side..[[;
	real const sqrt_gUjj = coord_sqrt_g_uu]]..side..side..[[(pt);
	
	real3 const vTilde = (eig)->vTilde;
	real3 const vTildeL = coord_lower(vTilde, pt);
	real const hTotal = (eig)->hTotal;
	real const vTildeSq = real3_dot(vTilde, vTildeL);
	real const Cs = (eig)->Cs;
	real const Cs_over_sqrt_gUjj = Cs / sqrt_gUjj; 
	real const rhoBar = (eig)->rhoBar;
	real const k = (eig)->k;
	real const omega = (eig)->omega;
	/* g^ij for fixed j=side */
]] .. prefix

	local gUdef = "\treal3 gUj = _real3(\n"
	for i=0,2 do
		gUdef = gUdef .. "\t\tcoord_g_uu"..side..i.."(pt)"..(i<2 and "," or "").."\n"
	end
	gUdef = gUdef .. "\t);\n"
	prefix = gUdef .. prefix
	prefixes[side] = prefix
end
?>

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: coord_g_uu coord_g_uu## coord_sqrt_g_uu## coord_lower <?=waves_t?> 

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	if (n.side == 0) {\
		<?=prefixes[0]:gsub("\n", "\\\n")?>\
		real const denom = 2. * Cs * Cs;\
		real const invDenom = 1. / denom;\
		real const invRhoBar = 1. / (eig)->rhoBar;\
		real const tmp = 1. - 2./3. * C_v_over_R;\
		real const sqrt_gUxx = sqrt_gUjj;\
		(result)->ptr[0] =\
			((X)->ptr[0] * (.5 * R_over_C_v * vTildeSq + Cs * vTilde.x / sqrt_gUxx)\
				+ (X)->ptr[1] * (-R_over_C_v * vTildeL.x - Cs / sqrt_gUxx)\
				+ (X)->ptr[2] * -R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * -R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * R_over_C_v\
				+ (X)->ptr[5] * (2./3. - R_over_C_v)\
			) * invDenom;\
		(result)->ptr[1] =\
			((X)->ptr[0] * (denom - R_over_C_v * vTildeSq)\
				+ (X)->ptr[1] * 2. * R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * 2. * R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * 2. * R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * -2. * R_over_C_v\
				+ (X)->ptr[5] * 2. * (R_over_C_v - 2./3.)\
			) * invDenom;\
		(result)->ptr[2] =\
			((X)->ptr[0] * (vTilde.x * gU.xy / gU.xx - vTilde.y)\
				+ (X)->ptr[1] * -gU.xy / gU.xx\
				+ (X)->ptr[2]) * invRhoBar;\
		(result)->ptr[3] =\
			((X)->ptr[0] * (vTilde.x * gU.xz / gU.xx - vTilde.z)\
				+ (X)->ptr[1] * -gU.xz / gU.xx\
				+ (X)->ptr[3]) * invRhoBar;\
		(result)->ptr[4] =\
			((X)->ptr[0] * -k\
				+ (X)->ptr[5]) * invRhoBar;\
		(result)->ptr[5] =\
			((X)->ptr[0] * -omega\
				+ (X)->ptr[6]) * invRhoBar;\
		(result)->ptr[6] =\
			((X)->ptr[0] * (.5 * R_over_C_v * vTildeSq - Cs * vTilde.x / sqrt_gUxx)\
				+ (X)->ptr[1] * (-R_over_C_v * vTildeL.x + Cs / sqrt_gUxx)\
				+ (X)->ptr[2] * -R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * -R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * R_over_C_v\
				+ (X)->ptr[5] * (2./3. - R_over_C_v)\
			) * invDenom;\
	} else if (n.side == 1) {\
		<?=prefixes[1]:gsub("\n", "\\\n")?>\
		real const denom = 2. * Cs * Cs;\
		real const invDenom = 1. / denom;\
		real const invRhoBar = 1. / (eig)->rhoBar;\
		real const tmp = 1. - 2./3. * C_v_over_R;\
		real const sqrt_gUyy = sqrt_gUjj;\
		(result)->ptr[0] =\
			((X)->ptr[0] * (.5 * R_over_C_v * vTildeSq + Cs * vTilde.y / sqrt_gUyy)\
				+ (X)->ptr[1] * -R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * (-R_over_C_v * vTildeL.y - Cs / sqrt_gUyy)\
				+ (X)->ptr[3] * -R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * R_over_C_v\
				+ (X)->ptr[5] * (2./3. - R_over_C_v)\
			) * invDenom;\
		(result)->ptr[1] =\
			((X)->ptr[0] * (vTilde.y * gU.xy / gU.yy - vTilde.x)\
				+ (X)->ptr[1]\
				+ (X)->ptr[2] * -gU.xy / gU.yy) * invRhoBar;\
		(result)->ptr[2] =\
			((X)->ptr[0] * (denom - R_over_C_v * vTildeSq)\
				+ (X)->ptr[1] * 2. * R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * 2. * R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * 2. * R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * -2. * R_over_C_v\
				+ (X)->ptr[5] * 2. * (R_over_C_v - 2./3.)\
			) * invDenom;\
		(result)->ptr[3] =\
			((X)->ptr[0] * (vTilde.y * gU.yz / gU.yy - vTilde.z)\
				+ (X)->ptr[2] * -gU.yz / gU.yy\
				+ (X)->ptr[3]) * invRhoBar;\
		(result)->ptr[4] =\
			((X)->ptr[0] * -k\
				+ (X)->ptr[5]) * invRhoBar;\
		(result)->ptr[5] =\
			((X)->ptr[0] * -omega\
				+ (X)->ptr[6]) * invRhoBar;\
		(result)->ptr[6] =\
			((X)->ptr[0] * (.5 * R_over_C_v * vTildeSq - Cs * vTilde.y / sqrt_gUyy)\
				+ (X)->ptr[1] * -R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * (-R_over_C_v * vTildeL.y + Cs / sqrt_gUyy)\
				+ (X)->ptr[3] * -R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * R_over_C_v\
				+ (X)->ptr[5] * (2./3. - R_over_C_v)\
			) * invDenom;\
	} else if (n.side == 2) {\
		<?=prefixes[2]:gsub("\n", "\\\n")?>\
		real const denom = 2. * Cs * Cs;\
		real const invDenom = 1. / denom;\
		real const invRhoBar = 1. / (eig)->rhoBar;\
		real const tmp = 1. - 2./3. * C_v_over_R;\
		real const sqrt_gUzz = sqrt_gUjj;\
		(result)->ptr[0] =\
			((X)->ptr[0] * (.5 * R_over_C_v * vTildeSq + Cs * vTilde.z / sqrt_gUzz)\
				+ (X)->ptr[1] * -R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * -R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * (-R_over_C_v * vTildeL.z - Cs / sqrt_gUzz)\
				+ (X)->ptr[4] * R_over_C_v\
				+ (X)->ptr[5] * (2./3. - R_over_C_v)\
			) * invDenom;\
		(result)->ptr[1] =\
			((X)->ptr[0] * (vTilde.z * gU.xz / gU.zz - vTilde.x)\
				+ (X)->ptr[1]\
				+ (X)->ptr[3] * -gU.xz / gU.zz) * invRhoBar;\
		(result)->ptr[2] =\
			((X)->ptr[0] * (vTilde.z * gU.yz / gU.zz - vTilde.y)\
				+ (X)->ptr[2]\
				+ (X)->ptr[3] * -gU.yz / gU.zz) * invRhoBar;\
		(result)->ptr[3] =\
			((X)->ptr[0] * (denom - R_over_C_v * vTildeSq)\
				+ (X)->ptr[1] * 2. * R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * 2. * R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * 2. * R_over_C_v * vTildeL.z\
				+ (X)->ptr[4] * -2. * R_over_C_v\
				+ (X)->ptr[5] * 2. * (R_over_C_v - 2./3.)\
			) * invDenom;\
		(result)->ptr[4] =\
			((X)->ptr[0] * -k\
				+ (X)->ptr[5]) * invRhoBar;\
		(result)->ptr[5] =\
			((X)->ptr[0] * -omega\
				+ (X)->ptr[6]) * invRhoBar;\
		(result)->ptr[6] =\
			((X)->ptr[0] * (.5 * R_over_C_v * vTildeSq - Cs * vTilde.z / sqrt_gUzz)\
				+ (X)->ptr[1] * -R_over_C_v * vTildeL.x\
				+ (X)->ptr[2] * -R_over_C_v * vTildeL.y\
				+ (X)->ptr[3] * (-R_over_C_v * vTildeL.z + Cs / sqrt_gUzz)\
				+ (X)->ptr[4] * R_over_C_v\
				+ (X)->ptr[5] * (2./3. - R_over_C_v)\
			) * invDenom;\
	}\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: coord_g_uu coord_g_uu## coord_sqrt_g_uu## coord_lower <?=waves_t?> 

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> const * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	if (n.side == 0) {\
		<?=prefixes[0]:gsub("\n", "\\\n")?>\
		real const sqrt_gUxx = sqrt_gUjj;\
		(result)->ptr[0] =\
			(X)->ptr[0] + (X)->ptr[1] + (X)->ptr[6];\
		(result)->ptr[1] =\
			(X)->ptr[0] * (vTilde.x - Cs * sqrt_gUxx)\
				+ (X)->ptr[1] * vTilde.x\
				+ (X)->ptr[6] * (vTilde.x + Cs * sqrt_gUxx);\
		(result)->ptr[2] =\
			(X)->ptr[0] * (vTilde.y - Cs * gU.xy / sqrt_gUxx)\
				+ (X)->ptr[1] * vTilde.y\
				+ (X)->ptr[2] * rhoBar\
				+ (X)->ptr[6] * (vTilde.y + Cs * gU.xy / sqrt_gUxx);\
		(result)->ptr[3] =\
			(X)->ptr[0] * (vTilde.z - Cs * gU.xz / sqrt_gUxx)\
				+ (X)->ptr[1] * vTilde.z\
				+ (X)->ptr[3] * rhoBar\
				+ (X)->ptr[6] * (vTilde.z + Cs * gU.xz / sqrt_gUxx);\
		(result)->ptr[4] =\
			(X)->ptr[0] * (hTotal - Cs * vTilde.x / sqrt_gUxx)\
				+ (X)->ptr[1] * (.5 * vTildeSq + (1. - 2./3. * C_v_over_R) * k)\
				+ (X)->ptr[2] * rhoBar * vTildeL.y\
				+ (X)->ptr[3] * rhoBar * vTildeL.z\
				+ (X)->ptr[4] * rhoBar * (1. - 2./3. * C_v_over_R)\
				+ (X)->ptr[6] * (hTotal + Cs * vTilde.x / sqrt_gUxx);\
		(result)->ptr[5] =\
			((X)->ptr[0] + (X)->ptr[1] + (X)->ptr[6]) * k\
				+ (X)->ptr[4] * rhoBar;\
		(result)->ptr[6] =\
			((X)->ptr[0] + (X)->ptr[1] + (X)->ptr[6]) * omega\
				+ (X)->ptr[5] * rhoBar;\
		(result)->ptr[7] =\
			0;\
	} else if (n.side == 1) {\
		<?=prefixes[1]:gsub("\n", "\\\n")?>\
		real const sqrt_gUyy = sqrt_gUjj;\
		(result)->ptr[0] =\
			(X)->ptr[0] + (X)->ptr[2] + (X)->ptr[6];\
		(result)->ptr[1] =\
			(X)->ptr[0] * (vTilde.x - Cs * gU.xy / sqrt_gUyy)\
				+ (X)->ptr[1] * rhoBar\
				+ (X)->ptr[2] * vTilde.x\
				+ (X)->ptr[6] * (vTilde.x + Cs * gU.xy / sqrt_gUyy);\
		(result)->ptr[2] =\
			(X)->ptr[0] * (vTilde.y - Cs * sqrt_gUyy)\
				+ (X)->ptr[2] * vTilde.y\
				+ (X)->ptr[6] * (vTilde.y + Cs * sqrt_gUyy);\
		(result)->ptr[3] =\
			(X)->ptr[0] * (vTilde.z - Cs * gU.yz / sqrt_gUyy)\
				+ (X)->ptr[2] * vTilde.z\
				+ (X)->ptr[3] * rhoBar\
				+ (X)->ptr[6] * (vTilde.z + Cs * gU.yz / sqrt_gUyy);\
		(result)->ptr[4] =\
			(X)->ptr[0] * (hTotal - Cs * vTilde.y / sqrt_gUyy)\
				+ (X)->ptr[1] * rhoBar * vTildeL.x\
				+ (X)->ptr[2] * (.5 * vTildeSq + (1. - 2./3. * C_v_over_R) * k)\
				+ (X)->ptr[3] * rhoBar * vTildeL.z\
				+ (X)->ptr[4] * (1. - 2./3. * C_v_over_R) * k\
				+ (X)->ptr[6] * (hTotal + Cs * vTilde.y / sqrt_gUyy);\
		(result)->ptr[5] =\
			((X)->ptr[0] + (X)->ptr[2] + (X)->ptr[6]) * k\
				+ (X)->ptr[4] * rhoBar;\
		(result)->ptr[6] =\
			((X)->ptr[0] + (X)->ptr[2] + (X)->ptr[6]) * omega\
				+ (X)->ptr[5] * rhoBar;\
		(result)->ptr[7] =\
			0;\
	} else if (n.side == 2) {\
		<?=prefixes[2]:gsub("\n", "\\\n")?>\
		real const sqrt_gUzz = sqrt_gUjj;\
		(result)->ptr[0] =\
			(X)->ptr[0] + (X)->ptr[3] + (X)->ptr[6];\
		(result)->ptr[1] =\
			(X)->ptr[0] * (vTilde.x - Cs * gU.xz / sqrt_gUzz)\
				+ (X)->ptr[1] * rhoBar\
				+ (X)->ptr[3] * vTilde.x\
				+ (X)->ptr[6] * (vTilde.x + Cs * gU.xz / sqrt_gUzz);\
		(result)->ptr[2] =\
			(X)->ptr[0] * (vTilde.y - Cs * gU.yz / sqrt_gUzz)\
				+ (X)->ptr[2] * rhoBar\
				+ (X)->ptr[3] * vTilde.y\
				+ (X)->ptr[6] * (vTilde.y + Cs * gU.yz / sqrt_gUzz);\
		(result)->ptr[3] =\
			(X)->ptr[0] * (vTilde.z - Cs * sqrt_gUzz)\
				+ (X)->ptr[3] * vTilde.z\
				+ (X)->ptr[6] * (vTilde.z + Cs * sqrt_gUzz);\
		(result)->ptr[4] =\
			(X)->ptr[0] * (hTotal - Cs * vTilde.z / sqrt_gUzz)\
				+ (X)->ptr[1] * rhoBar * vTildeL.x\
				+ (X)->ptr[2] * rhoBar * vTildeL.y\
				+ (X)->ptr[3] * (.5 * vTildeSq + (1. - 2./3. * C_v_over_R) * k)\
				+ (X)->ptr[4] * (1. - 2./3. * C_v_over_R) * k\
				+ (X)->ptr[6] * (hTotal + Cs * vTilde.z / sqrt_gUzz);\
		(result)->ptr[5] =\
			((X)->ptr[0] + (X)->ptr[3] + (X)->ptr[6]) * k\
				+ (X)->ptr[4] * rhoBar;\
		(result)->ptr[6] =\
			((X)->ptr[0] + (X)->ptr[3] + (X)->ptr[6]) * omega\
				+ (X)->ptr[5] * rhoBar;\
		(result)->ptr[7] =\
			0;\
	}\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

static inline void <?=eigen_fluxTransform?>(
	<?=cons_t?> const * const result,
	constant <?=solver_t?> const * const solver,
	<?=eigen_t?> eig,
	<?=cons_t?> X,
	real3 pt,
	normal_t n
) {
	if (n.side == 0) {
		<?=prefixes[0]:gsub("\n", "\\\n")?>
		real const CsSq = Cs * Cs;
		real const PStar = CsSq * rhoBar / (1. + R_over_C_v);
		*(result) = (<?=cons_t?>){.ptr={
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
		<?=prefixes[1]:gsub("\n", "\\\n")?>
		real const CsSq = Cs * Cs;
		real const PStar = CsSq * rhoBar / (1. + R_over_C_v);
		*(result) = (<?=cons_t?>){.ptr={
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
		<?=prefixes[2]:gsub("\n", "\\\n")?>
		real const CsSq = Cs * Cs;
		real const PStar = CsSq * rhoBar / (1. + R_over_C_v);
		*(result) = (<?=cons_t?>){.ptr={
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

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=primFromCons?>

// used by PLM

#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> const * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
	(result)->rhoBar = W.rhoBar;\
	(result)->vTilde = W.vTilde;\
	(result)->k = W.k;\
	(result)->omega = W.omega;\
	(result)->vTildeSq = coordLenSq(W.vTilde, pt);\
	(result)->hTotal = (W.PStar + (U)->rhoBar_eTotalTilde) / (U)->rhoBar - (U)->ePot;\
	real const eKinTilde = .5 * (result)->vTildeSq;\
	real const CsSq = R_over_C_v * ((result)->hTotal - eKinTilde) + 2./3. * W.k;\
	(result)->Cs = sqrt(CsSq);\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=primFromCons?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

<? if not require "hydro.coord.cartesian".is(solver.coord) then ?>
	/* connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system */
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	real3 m_conn_vv = coord_conn_apply23(W.vTilde, U->rhoBar_vTilde, x);
	deriv->rhoBar_vTilde = real3_sub(deriv->rhoBar_vTilde, m_conn_vv);	/* -Conn^i_jk rhoBar vTilde^j vTilde^k  */
	deriv->rhoBar_vTilde = real3_add(deriv->rhoBar_vTilde, real3_real_mul(coord_raise(coord_conn_trace13(x), x), W.PStar));		/* +Conn^j_kj g^ki PStar */
<? end ?>
}
