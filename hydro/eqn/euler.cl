//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=eqn_common?>
// eqn_common is for all the calc_* stuff

static inline <?=prim_t?> <?=primFromCons?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	real3 const & x
) {
	<?=prim_t?> result;
	if (U.rho < solver.rhoMin) {
		result.rho = 0.;
		result.v = {};
		result.P = 0.;
	} else {
		result.rho = U.rho;
		result.v = <?=calc_v?>(U);
		result.P = <?=calc_P?>(solver, U, x);
	}
	result.rho = U.rho;
	result.v = <?=calc_v?>(U);
	result.P = <?=calc_P?>(solver, U, x);
	result.ePot = U.ePot;
	return result;
}

//// MODULE_NAME: <?=consFromPrim?>
//// MODULE_DEPENDS: real3 <?=eqn_common?> <?=solver_t?> <?=prim_t?> <?=cons_t?>
// eqn_common is for all the calc_* stuff

static inline <?=cons_t?> <?=consFromPrim?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & W,
	real3 const pt
) {
	<?=cons_t?> result;
	result.rho = W.rho;
	result.m = W.v * W.rho;
	result.ETotal = <?=calc_ETotal?>(solver, W, pt);
	result.ePot = W.ePot;
	return result;
}

//// MODULE_NAME: <?=apply_dU_dW?>
//// MODULE_DEPENDS: <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>
// only used by PLM

static inline <?=cons_t?> <?=apply_dU_dW?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & WA,
	<?=prim_t?> const & W,
	real3 const x
) {
	<?=cons_t?> result;
	real3 const WA_vL = coord_lower(WA.v, x);
	result.rho = W.rho;
	result.m = 
		WA.v * W.rho
		+ W.v * WA.rho;
	result.ETotal = W.rho * (real).5 * dot(WA.v, WA_vL)
		+ WA.rho * dot(W.v, WA_vL)
		+ W.P / (solver.heatCapacityRatio - 1.);
	result.ePot = W.ePot;
	return result;
}

//// MODULE_NAME: <?=apply_dW_dU?>	
//// MODULE_DEPENDS: <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>
// only used by PLM

static inline <?=prim_t?> <?=apply_dW_dU?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & WA,
	<?=cons_t?> const & U,
	real3 const x
) {
	<?=prim_t?> result;
	real3 const WA_vL = coord_lower(WA.v, x);
	result.rho = U.rho;
	if (U.rho < solver.rhoMin) {
		result.v = {};
		result.P = 0.;
	} else {
		result.v = 
			U.m * (1. / WA.rho)
			- WA.v * (U.rho / WA.rho);
		result.P = (solver.heatCapacityRatio - 1.) * (
			.5 * dot(WA.v, WA_vL) * U.rho 
			- dot(U.m, WA_vL)
			+ U.ETotal);
	}
	result.ePot = U.ePot;
	return result;
}

//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: <?=coordLenSq?> <?=cons_t?> <?=prim_t?> <?=waves_t?> <?=eigen_t?> <?=eqn_guiVars_compileTime?>

static inline real <?=calc_H?>(
	constant <?=solver_t?> const & solver,
	real const P
) {
	return (P * (solver.heatCapacityRatio / (solver.heatCapacityRatio - 1.)));
}

static inline real <?=calc_h?>(
	constant <?=solver_t?> const & solver,
	real const rho,
	real const P
) {
	return <?=calc_H?>(solver, P) / rho;
}

static inline real <?=calc_HTotal?>(
	real const P,
	real const ETotal
) {
	return P + ETotal;
}

static inline real <?=calc_hTotal?>(
	real const rho,
	real const P,
	real const ETotal
) {
	return <?=calc_HTotal?>(P, ETotal) / rho;
}

static inline real <?=calc_eKin?>(
	<?=prim_t?> const & W,
	real3 const x
) {
	return (real).5 * coordLenSq(W.v, x);
}

static inline real <?=calc_EKin?>(
	<?=prim_t?> const & W,
	real3 const x
) {
	return W.rho * <?=calc_eKin?>(W, x);
}

static inline real <?=calc_EInt?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & W
) {
	return W.P / (solver.heatCapacityRatio - 1.);
}

static inline real <?=calc_eInt?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & W
) {
	return <?=calc_EInt?>(solver, W) / W.rho;
}

static inline real <?=calc_EKin_fromCons?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	real3 const x
) {
	if (U.rho < solver.rhoMin) return 0;
	return (real).5 * coordLenSq(U.m, x) / U.rho;
}

static inline real <?=calc_ETotal?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & W,
	real3 const x
) {
	return <?=calc_EKin?>(W, x) + <?=calc_EInt?>(solver, W);
}

static inline real <?=calc_P?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	real3 const x
) {
	if (U.rho < solver.rhoMin) return 0;
	return (solver.heatCapacityRatio - 1.) * 
		/*EInt=*/(U.ETotal - <?=calc_EKin_fromCons?>(solver, U, x));
}

static inline real <?=calc_Cs?>(
	constant <?=solver_t?> const & solver,
	<?=prim_t?> const & W
) {
	if (W.P <= solver.PMin) return 0.;
	if (W.rho < solver.rhoMin) return INFINITY;
	return sqrt(solver.heatCapacityRatio * W.P / W.rho);
}

static inline real <?=calc_Cs_fromCons?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	real3 const pt
) {
	real const P = <?=calc_P?>(solver, U, pt);
	if (P <= solver.PMin) {
		return 0;
	} else if (U.rho < solver.rhoMin) {
		return INFINITY;
	} else {
		return sqrt(solver.heatCapacityRatio * P / U.rho);
	}
}

static inline real <?=calc_eInt_fromCons?>(
	<?=cons_t?> const & U,
	real3 const x
) {
	return (U.ETotal - .5 * coordLenSq(U.m, x) / U.rho) / U.rho - U.ePot;
}

<? local materials = require "hydro.materials" ?>
static constexpr real C_v = <?=("%.50f"):format(materials.Air.C_v)?>;

static inline real <?=calc_T?>(
	<?=cons_t?> const & U,
	real3 const x
) {
	return <?=calc_eInt_fromCons?>(U, x) / C_v;
}

static inline real3 <?=calc_v?>(
	<?=cons_t?> const & U
) {
	return U.m * (1. / U.rho);
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=cartesianToCoord?>

/*
I've highjacked all of this.  It was a normal Euler eqn solver.
But I experimented with a curved-space solver.  
To get back to the original code,
just replace all the g_ab stuff with their constant values and simplify away.
*/
void <?=applyInitCondCell?>(
	constant <?=solver_t?> const & solver,
	constant <?=initCond_t?> const & initCond,
	global <?=cons_t?> & U,
	global <?=cell_t?> const & cell
) {
	real3 const x = cell.pos;
	real3 const mids = (solver.initCondMins + solver.initCondMaxs) * (real).5;
	bool const lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	// these are all standard for all init/euler.lua initial conditions
	real rho = 0;
	real3 v = {};
	real P = 0;
	real3 D = {};
	real3 B = {};
	real ePot = 0;
	
<?=initCode()?>
	
	<?=prim_t?> W = <?=prim_t?>()
		.set_rho(rho)
		.set_v(cartesianToCoord(v, x))
		.set_P(P)
		.set_ePot(ePot)
	;

	U = <?=consFromPrim?>(solver, W, x);
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=primFromCons?> <?=normal_t?>

static inline <?=cons_t?> <?=fluxFromCons?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	<?=cell_t?> const & cell,
	<?=normal_t?> const n
) {
	<?=cons_t?> resultF;
	<?=prim_t?> W = <?=primFromCons?>(solver, U, cell.pos);
	real const v_n = normal_vecDotN1(n, W.v);
	resultF.rho = U.rho * v_n;
	resultF.m = 
		U.m * v_n
		+ normal_u1(n) * W.P;
	real const HTotal = U.ETotal + W.P;
	resultF.ETotal = HTotal * v_n;
	resultF.ePot = 0;
	return resultF;
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: real3x3 <?=primFromCons?>
// added by request only, so I don't have to compile the real3x3 code. 
// not used at the moment

static inline <?=range_t?> <?=calcCellMinMaxEigenvalues?>(
	constant <?=solver_t?> const & solver,
	global <?=cons_t?> const & U,
	real3 const pt,
	real3x3 const nL,
	real3x3 const nU,
	real const nLen
) {
	<?=prim_t?> W = <?=primFromCons?>(solver, U, pt);
	real const v_n = dot(W.v, nL.x);
	real const Cs = <?=calc_Cs?>(solver, &W);
	real const Cs_nLen = Cs * nLen;
	<?=range_t?> result;
	result.min = v_n - Cs_nLen; 
	result.max = v_n + Cs_nLen;
	return result;
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=normal_t?> <?=coord_lower?> <?=cons_t?> <?=prim_t?> <?=eigen_t?> <?=primFromCons?> <?=eqn_common?>
// eqn_common is for all the calc_* stuff

// used by PLM
static inline <?=eigen_t?> <?=eigen_forCell?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	<?=cell_t?> const & cell,
	<?=normal_t?> const n
) {
	<?=prim_t?> W = <?=primFromCons?>(solver, U, cell.pos);
	real3 const vL = coord_lower(W.v, cell.pos);
	real const vSq = dot(W.v, vL);
	real const v_n = normal_vecDotN1(n, W.v);
	real const eKin = .5 * vSq;
	real const hTotal = <?=calc_hTotal?>(W.rho, W.P, U.ETotal);
	real const CsSq = (solver.heatCapacityRatio - 1.) * (hTotal - eKin);
	real const Cs = sqrt(CsSq);
	<?=eigen_t?> result;
	result.rho = W.rho;
	result.v = W.v;
	result.vSq = vSq;
	result.vL = vL;
	result.hTotal = hTotal;
	result.Cs = Cs;
	return result;
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=primFromCons?> <?=eigen_t?> <?=normal_t?> <?=coord_lower?>

//used by the mesh version
static inline <?=eigen_t?> <?=eigen_forInterface?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & UL,
	<?=cons_t?> const & UR,
	<?=cell_t?> const & cellL,
	<?=cell_t?> const & cellR,
	real3 const pt,
	<?=normal_t?> const n
) {
	<?=eigen_t?> result;
	real const rhoEpsilon = 1e-5;

	if (UL.rho < rhoEpsilon && UR.rho < rhoEpsilon) {
		// left and right are both vacuum:
		result.rho = 0.;
		result.v = {};
		result.vSq = 0.;
		result.vL = {};
		result.hTotal = 0;
		result.Cs = 0;
	} else {
		if (UL.rho < rhoEpsilon) {
			// left is vacuum:
			<?=prim_t?> WR = <?=primFromCons?>(solver, UR, cellR.pos);
			result.rho = UR.rho;
			result.v = WR.v;
			result.vL = coord_lower(WR.v, pt);
			result.vSq = dot(WR.v, result.vL);
			result.hTotal = <?=calc_hTotal?>(WR.rho, WR.P, UR.ETotal);
			result.Cs = <?=calc_Cs?>(solver, WR);
		} else if (UR.rho < rhoEpsilon) {
			// right is vacuum:
			<?=prim_t?> WL = <?=primFromCons?>(solver, UL, cellL.pos);
			result.rho = UL.rho;
			result.v = WL.v;
			result.vL = coord_lower(WL.v, pt);
			result.vSq = dot(WL.v, result.vL);
			result.hTotal = <?=calc_hTotal?>(WL.rho, WL.P, UL.ETotal);
			result.Cs = <?=calc_Cs?>(solver, WL);
		} else {
			<?=prim_t?> WL = <?=primFromCons?>(solver, UL, cellL.pos);
			real const sqrtRhoL = sqrt(WL.rho);
			real3 const vLeft = WL.v;
			real const hTotalL = <?=calc_hTotal?>(WL.rho, WL.P, UL.ETotal);

			<?=prim_t?> WR = <?=primFromCons?>(solver, UR, cellR.pos);
			real const sqrtRhoR = sqrt(WR.rho);
			real3 const vR = WR.v;
			real const hTotalR = <?=calc_hTotal?>(WR.rho, WR.P, UR.ETotal);

			real const invDenom = 1./(sqrtRhoL + sqrtRhoR);

			//Roe-averaged
			result.rho = sqrtRhoL * sqrtRhoR;
			real3 const v = 
					vLeft * (sqrtRhoL * invDenom)
					+ vR * (sqrtRhoR * invDenom);
			real const hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);

			//derived:
			real3 const vLower = coord_lower(v, pt);
			real const vSq = dot(v, vLower);
			real const eKin = .5 * vSq;
			real const h = hTotal - eKin;
			// TODO verify hTotal = 1/2 v^2 + Cs^2 / (gamma-1)
			if (h < rhoEpsilon) {
				result.hTotal = eKin;
				result.Cs = 0.;
			} else {
				result.hTotal = hTotal;
				real const CsSq = h < rhoEpsilon ? 0. : (solver.heatCapacityRatio - 1.) * h;
				result.Cs = sqrt(CsSq);
			}

			result.v = v;
			result.vSq = vSq;
			result.vL = vLower;
		}
	}
	return result;
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>

static inline <?=waves_t?> <?=eigen_leftTransform?>(
	constant <?=solver_t?> const & solver,
	<?=eigen_t?> const & eig,
	<?=cons_t?> const & X,
	real3 const pt,
	<?=normal_t?> n
) {
	<?=waves_t?> result;
	if (eig.rho < solver.rhoMin) {
		result.ptr[0] = X.ptr[0];
		result.ptr[1] = X.ptr[1];
		result.ptr[2] = X.ptr[2];
		result.ptr[3] = X.ptr[3];
		result.ptr[4] = X.ptr[4];
	} else {
		real3 const v_n = normal_vecDotNs(n, eig.v);
		real const nLen = normal_len(n);
		real const inv_nLen = 1. / nLen;
		real const denom = 2. * eig.Cs * eig.Cs;
		real const invDenom = 1. / denom;
		real const gamma_1 = solver.heatCapacityRatio - 1.;
		result.ptr[0] = (
				X.ptr[0] * (.5 * gamma_1 * eig.vSq + eig.Cs * v_n.x * inv_nLen)
				+ X.ptr[1] * (-gamma_1 * eig.vL.x - eig.Cs * normal_l1x_over_len(n))
				+ X.ptr[2] * (-gamma_1 * eig.vL.y - eig.Cs * normal_l1y_over_len(n))
				+ X.ptr[3] * (-gamma_1 * eig.vL.z - eig.Cs * normal_l1z_over_len(n))
				+ X.ptr[4] * gamma_1
			) * invDenom;
		result.ptr[1] =
			(
				X.ptr[0] * (denom - gamma_1 * eig.vSq)
				+ X.ptr[1] * 2. * gamma_1 * eig.vL.x
				+ X.ptr[2] * 2. * gamma_1 * eig.vL.y
				+ X.ptr[3] * 2. * gamma_1 * eig.vL.z
				+ X.ptr[4] * -2. * gamma_1
			) * invDenom;
		result.ptr[2] =
			X.ptr[0] * -v_n.y
			+ X.ptr[1] * normal_l2x(n)
			+ X.ptr[2] * normal_l2y(n)
			+ X.ptr[3] * normal_l2z(n);
		result.ptr[3] =
			X.ptr[0] * -v_n.z
			+ X.ptr[1] * normal_l3x(n)
			+ X.ptr[2] * normal_l3y(n)
			+ X.ptr[3] * normal_l3z(n);
		result.ptr[4] =
			(
				X.ptr[0] * (.5 * gamma_1 * eig.vSq - eig.Cs * v_n.x * inv_nLen)
				+ X.ptr[1] * (-gamma_1 * eig.vL.x + eig.Cs * normal_l1x_over_len(n))
				+ X.ptr[2] * (-gamma_1 * eig.vL.y + eig.Cs * normal_l1y_over_len(n))
				+ X.ptr[3] * (-gamma_1 * eig.vL.z + eig.Cs * normal_l1z_over_len(n))
				+ X.ptr[4] * gamma_1
			) * invDenom;
	}
	return result;
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>

static inline <?=cons_t?> <?=eigen_rightTransform?>(
	constant <?=solver_t?> const & solver,
	<?=eigen_t?> const & eig,
	<?=waves_t?> const & X,
	real3 const pt,
	<?=normal_t?> const n
) {
	<?=cons_t?> result;
	if (eig.rho < solver.rhoMin) {
		result.ptr[0] = X.ptr[0];
		result.ptr[1] = X.ptr[1];
		result.ptr[2] = X.ptr[2];
		result.ptr[3] = X.ptr[3];
		result.ptr[4] = X.ptr[4];
	} else {
		real3 const v_n = normal_vecDotNs(n, eig.v);
		real const nLen = normal_len(n);
		real const inv_nLen = 1. / nLen;
		result.ptr[0] =
			X.ptr[0]
			+ X.ptr[1]
			+ X.ptr[4];
		result.ptr[1] =
			X.ptr[0] * (eig.v.x - eig.Cs * normal_u1x_over_len(n))
			+ X.ptr[1] * eig.v.x
			+ X.ptr[2] * normal_u2x(n)
			+ X.ptr[3] * normal_u3x(n)
			+ X.ptr[4] * (eig.v.x + eig.Cs * normal_u1x_over_len(n));
		result.ptr[2] =
			X.ptr[0] * (eig.v.y - eig.Cs * normal_u1y_over_len(n))
			+ X.ptr[1] * eig.v.y
			+ X.ptr[2] * normal_u2y(n)
			+ X.ptr[3] * normal_u3y(n)
			+ X.ptr[4] * (eig.v.y + eig.Cs * normal_u1y_over_len(n));
		result.ptr[3] =
			X.ptr[0] * (eig.v.z - eig.Cs * normal_u1z_over_len(n))
			+ X.ptr[1] * eig.v.z
			+ X.ptr[2] * normal_u2z(n)
			+ X.ptr[3] * normal_u3z(n)
			+ X.ptr[4] * (eig.v.z + eig.Cs * normal_u1z_over_len(n));
		result.ptr[4] =
			X.ptr[0] * (eig.hTotal - eig.Cs * v_n.x * inv_nLen)
			+ X.ptr[1] * (real).5 * eig.vSq
			+ X.ptr[2] * v_n.y
			+ X.ptr[3] * v_n.z
			+ X.ptr[4] * (eig.hTotal + eig.Cs * v_n.x * inv_nLen);
	}
	result.ptr[5] = 0;
	return result;
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>
// Not used anymore.  was used by Roe, but I switched that to a <?=fluxFromCons?>.
// <?=fluxFromCons?> only matches <?=eigen_fluxTransform?> when the eig properties are derived from X_ 

static inline <?=eigen_fluxTransform?>(
	constant <?=solver_t?> const & solver,
	<?=eigen_t?> const & eig,
	<?=cons_t?> const & X,
	<?=cell_t?> const & cell,
	<?=normal_t?> const n
) {
	<?=cons_t?> result;
	if (eig.rho < solver.rhoMin) {
		result.ptr[0] = X.ptr[0];
		result.ptr[1] = X.ptr[1];
		result.ptr[2] = X.ptr[2];
		result.ptr[3] = X.ptr[3];
		result.ptr[4] = X.ptr[4];
	} else {
		real3 const v_n = normal_vecDotNs(n, eig.v);
		real const nLen = normal_len(n);
		real const gamma = solver.heatCapacityRatio;
		real const gamma_1 = gamma - 1.;
		real const gamma_2 = gamma - 2.;
	
		result.ptr[0] =
			X.ptr[1] * normal_l1x(n)
			+ X.ptr[2] * normal_l1y(n)
			+ X.ptr[3] * normal_l1z(n);
	
		result.ptr[1] =
			X.ptr[0] * (-v_n.x * eig.v.x + gamma_1 * (real).5 * eig.vSq * normal_u1x(n))
			+ X.ptr[1] * (eig.v.x * normal_l1x(n) - gamma_2 * normal_u1x(n) * eig.vL.x + v_n.x)
			+ X.ptr[2] * (eig.v.x * normal_l1y(n) - gamma_2 * normal_u1x(n) * eig.vL.y)
			+ X.ptr[3] * (eig.v.x * normal_l1z(n) - gamma_2 * normal_u1x(n) * eig.vL.z)
			+ X.ptr[4] * gamma_1 * normal_u1x(n);
	
		result.ptr[2] =
			X.ptr[0] * (-v_n.x * eig.v.y + gamma_1 * (real).5 * eig.vSq * normal_u1y(n))
			+ X.ptr[1] * (eig.v.y * normal_l1x(n) - gamma_2 * normal_u1y(n) * eig.vL.x)
			+ X.ptr[2] * (eig.v.y * normal_l1y(n) - gamma_2 * normal_u1y(n) * eig.vL.y + v_n.x)
			+ X.ptr[3] * (eig.v.y * normal_l1z(n) - gamma_2 * normal_u1y(n) * eig.vL.z)
			+ X.ptr[4] * gamma_1 * normal_u1y(n);
	
		result.ptr[3] =
			X.ptr[0] * (-v_n.x * eig.v.z + gamma_1 * (real).5 * eig.vSq * normal_u1z(n))
			+ X.ptr[1] * (eig.v.z * normal_l1x(n) - gamma_2 * normal_u1z(n) * eig.vL.x)
			+ X.ptr[2] * (eig.v.z * normal_l1y(n) - gamma_2 * normal_u1z(n) * eig.vL.y)
			+ X.ptr[3] * (eig.v.z * normal_l1z(n) - gamma_2 * normal_u1z(n) * eig.vL.z + v_n.x)
			+ X.ptr[4] * gamma_1 * normal_u1z(n);
	
		result.ptr[4] =
			X.ptr[0] * v_n.x * (.5 * gamma_1 * eig.vSq - eig.hTotal)
			+ X.ptr[1] * (normal_l1x(n) * eig.hTotal - gamma_1 * v_n.x * eig.vL.x)
			+ X.ptr[2] * (normal_l1y(n) * eig.hTotal - gamma_1 * v_n.x * eig.vL.y)
			+ X.ptr[3] * (normal_l1z(n) * eig.hTotal - gamma_1 * v_n.x * eig.vL.z)
			+ X.ptr[4] * gamma * v_n.x;
	}

	result.ptr[5] = 0;
	return result;
}

<? if false then -- TODO sort <?=addSource?> out. ?>
//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS_NOGHOST?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	constant <?=solver_t?> const & solver = *psolver;
	<?=SETBOUNDS_NOGHOST?>();

	global <?=cons_t?> & deriv = derivBuf[index];
	global <?=cons_t?> const & U = UBuf[index];
	global <?=cell_t?> const & cell = cellBuf[index];
	real3 const x = cell.pos;

<? if false 
and solver.coord.vectorComponent == "anholonomic" 
and require "hydro.coord.cylinder":isa(solver.coord) 
then ?>
<? 	if true then -- 2009 Trangenstein, p.474, 1999 Toro, p.29, eqn.1.104, 1.105 ?>
	<? for side=0,1 do ?>{
		real3 xL = x; xL.s<?=side?> -= solver.grid_dx.s<?=side?>;
		real3 xR = x; xR.s<?=side?> += solver.grid_dx.s<?=side?>;
		
		global <?=cons_t?> const & UL = U[solver.stepsize.s<?=side?>];
		global <?=cons_t?> const & UR = U[solver.stepsize.s<?=side?>];
		real const PL = <?=calc_P?>(solver, UL, xL);
		real const PR = <?=calc_P?>(solver, UR, xR);
	
		deriv.m.s<?=side?> -= (PR - PL) / (2. * solver.grid_dx.s<?=side?>);
	}<? end ?>
<?	end ?>
<?	if false then -- 1999 Toro p.28 eqn.1.102, 1.103 ?>
	<?=cons_t?> F = <?=fluxFromCons?>(solver, U, cell, normal_forSide0(x));
	deriv.rho -= F.rho / x.x;
	deriv.m.x -= F.m.x / x.x;
	deriv.m.y -= F.m.y / x.x;
	deriv.m.z -= F.m.z / x.x;
	deriv.ETotal -= F.ETotal / x.x;
<?	end ?>
<? end ?>

<? do -- if not solver.coord.vectorComponent == "anholonomic" then ?>
<? if not (require "hydro.coord.cartesian":isa(solver.coord) 
		or solver.coord.vectorComponent == "cartesian")
then ?>
//// MODULE_DEPENDS: <?=primFromCons?> <?=coord_conn_apply23?> <?=coord_conn_trace23?> <?=coord_conn_apply13?>
/*
This is working for init conds with zero velocity.
Introducing constant velocity of v=[x=1,y=1]=[r=sqrt(2),theta=pi/4] in the init cond causes some numerical errors.
However the problem isn't the terms below -- because disabling this for zero-vel init conds causes things to become unsteady.
That means that the volume gradient in calcDerivFV is causing nonzero velocities to emerge, and this is cancelling them.
Maybe for an initial constant vel as large as sqrt(2) this fails, but it works only for small perturbations?
*/
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	<?=prim_t?> W = <?=primFromCons?>(solver, U, x);
	
	//- Γ^i_jk ρ v^j v^k 
	deriv.m -= coord_conn_apply23(W.v, U.m, x);
	
	//- Γ^i_jk g^jk P
	deriv.m -= coord_conn_trace23(x) * W.P;
	
	//+ (γ-1) ρ v^k v^l Γ_kjl g^ij
	deriv.m += coord_conn_apply13(W.v, U.m, x) * (solver.heatCapacityRatio - 1.);
	
	//- (γ-1) ρ v^j v^k v^l Γ_jkl
//	deriv.ETotal -= (solver.heatCapacityRatio - 1.) * coord_conn_apply123(W.v, W.v, U.m, x);	

	//+ c_jk^k * Flux^Ij
<? 	if false and solver.coord.vectorComponent == "anholonomic" then ?>
	real3 const commTrace = coord_tr23_c(x);
	<? for i=0,solver.dim-1 do ?>{
		<?=cons_t?> flux;
		calcFluxFromCons(&F, U, x);
		for (int j = 0; j < numIntStates; ++j) {
			deriv.ptr[j] += commTrace.s<?=i?> * flux.ptr[j];
		}
	}<? end ?>
<? 	end ?>
<? end ?>
<? end -- vectorComponent == "anholonomic" ?>
}
<? end ?>

<? if false then ?>
//// MODULE_NAME: <?=constrainU?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=primFromCons?> <?=consFromPrim?>

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	constant <?=solver_t?> const & solver = *psolver;
	<?=SETBOUNDS?>(0,0);
	real3 const x = cellBuf[index].pos;

	global <?=cons_t?> & U = UBuf[index];
	<?=prim_t?> W = <?=primFromCons?>(solver, U, x);

	if (W.rho < solver.rhoMin) W.rho = solver.rhoMin;
	if (W.P < solver.PMin) W.P = solver.PMin;

	U = <?=consFromPrim?>(solver, W, x);
}
<? end ?>
