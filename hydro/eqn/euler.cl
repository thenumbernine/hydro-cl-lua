//// MODULE_NAME: <?=Equation?>
//// MODULE_DEPENDS: <?=coordLenSq?> <?=waves_t?> <?=eigen_t?> <?=eqn_guiVars_compileTime?>
// Equation is for all the calc_* stuff

<?=eqn:template(require "ext.file" "hydro/eqn/eqn.clcpp":read())?>

// TODO 'namespace Equation' then 'class <?=Equation?>'
// or TODO namespace <?=Solver?> and then class Solver, class Eqn, etc
namespace <?=Equation?> {

// TODO this is usually the <?=prim_t?> structure ...
// which I want to get rid of in favor of member-class of the equation, Cons and Prim
// but this does match the codegen struct so far except ...
union Prim {
	real ptr[6];
	struct {
		real rho;
		real3 v;
		real P;
		real ePot;
	};
//// BEGIN EXCLUDE FROM FFI_CDEF
	Prim() {}
	Prim & set_rho(real const & value_) { rho = value_; return *this; }
	Prim & set_v(real3 const & value_) { v = value_; return *this; }
	Prim & set_P(real const & value_) { P = value_; return *this; }
	Prim & set_ePot(real const & value_) { ePot = value_; return *this; }
//// END EXCLUDE FROM FFI_CDEF
};

union Cons {
	real ptr[6];
	struct {
		real rho;
		real3 m;
		real ETotal;
		real ePot;
	};
//// BEGIN EXCLUDE FROM FFI_CDEF
	Cons() {}
	Cons & set_rho(real const & value_) { rho = value_; return *this; }
	Cons & set_m(real3 const & value_) { m = value_; return *this; }
	Cons & set_ETotal(real const & value_) { ETotal = value_; return *this; }
	Cons & set_ePot(real const & value_) { ePot = value_; return *this; }
//// END EXCLUDE FROM FFI_CDEF
};
	
//// MODULE_DEPENDS: <?=normal_t?> <?=initCond_t?> <?=cell_t?>
using Solver = <?=solver_t?>;
using Normal = <?=normal_t?>;
using InitCond = <?=initCond_t?>;
using Cell = <?=cell_t?>;

// TODO eventually:
//template<typename Solver>
struct Eqn : public Hydro::Eqn<Eqn, Solver, Cons, Prim, Normal, Cell> {
	// until then just use a typedef and a namespace

	// hmm for some reason crtp parent can't seem to 'using' see into the childs member class's 'using's
	//using Cons = ::<?=Equation?>::Cons;
	//using Prim = ::<?=Equation?>::Prim;
	//using Normal = <?=normal_t?>;
	//using Cell = <?=cell_t?>;

	static inline real calc_H(
		constant Solver const & solver,
		real const P
	) {
		return (P * (solver.heatCapacityRatio / (solver.heatCapacityRatio - 1.)));
	}

	static inline real calc_h(
		constant Solver const & solver,
		real const rho,
		real const P
	) {
		return calc_H(solver, P) / rho;
	}

	static inline real calc_HTotal(
		real const P,
		real const ETotal
	) {
		return P + ETotal;
	}

	static inline real calc_hTotal(
		real const rho,
		real const P,
		real const ETotal
	) {
		return calc_HTotal(P, ETotal) / rho;
	}

	static inline real calc_eKin(
		Prim const & W,
		real3 const pt
	) {
		return (real).5 * coordLenSq(W.v, pt);
	}

	static inline real calc_EKin(
		Prim const & W,
		real3 const pt
	) {
		return W.rho * calc_eKin(W, pt);
	}

	static inline real calc_EInt(
		constant Solver const & solver,
		Prim const & W
	) {
		return W.P / (solver.heatCapacityRatio - 1.);
	}

	static inline real calc_eInt(
		constant Solver const & solver,
		Prim const & W
	) {
		return calc_EInt(solver, W) / W.rho;
	}

	static inline real calc_EKin_fromCons(
		constant Solver const & solver,
		Cons const & U,
		real3 const pt
	) {
		if (U.rho < solver.rhoMin) return 0;
		return (real).5 * coordLenSq(U.m, pt) / U.rho;
	}

	static inline real calc_ETotal(
		constant Solver const & solver,
		Prim const & W,
		real3 const pt
	) {
		return calc_EKin(W, pt) + calc_EInt(solver, W);
	}

	static inline real calc_P(
		constant Solver const & solver,
		Cons const & U,
		real3 const pt
	) {
		if (U.rho < solver.rhoMin) return 0;
		return (solver.heatCapacityRatio - 1.) * 
			/*EInt=*/(U.ETotal - calc_EKin_fromCons(solver, U, pt));
	}

	static inline real calc_Cs(
		constant Solver const & solver,
		Prim const & W
	) {
		if (W.P <= solver.PMin) return 0.;
		if (W.rho < solver.rhoMin) return INFINITY;
		return sqrt(solver.heatCapacityRatio * W.P / W.rho);
	}

	static inline real calc_Cs_fromCons(
		constant Solver const & solver,
		Cons const & U,
		real3 const pt
	) {
		real const P = calc_P(solver, U, pt);
		if (P <= solver.PMin) {
			return 0;
		} else if (U.rho < solver.rhoMin) {
			return INFINITY;
		} else {
			return sqrt(solver.heatCapacityRatio * P / U.rho);
		}
	}

	static inline real calc_eInt_fromCons(
		Cons const & U,
		real3 const pt
	) {
		return (U.ETotal - .5 * coordLenSq(U.m, pt) / U.rho) / U.rho - U.ePot;
	}

<? local materials = require "hydro.materials" ?>
	static constexpr real C_v = <?=("%.50f"):format(materials.Air.C_v)?>;

	static inline real calc_T(
		Cons const & U,
		real3 const pt
	) {
		return calc_eInt_fromCons(U, pt) / C_v;
	}

	static inline real3 calc_v(
		Cons const & U
	) {
		return U.m * (1. / U.rho);
	}

	// TODO make this a member-function
	// so it can access equation-specific vars (like heatCapacityRatio)
	// and move those from Solver to Equation
	// then there's one less argument 
	
	static inline Prim primFromCons(
		constant Solver const & solver,
		Cons const & U,
		real3 const pt
	) {
		if (U.rho < solver.rhoMin) {
			return Prim()
			.set_rho(0)
			.set_v({})
			.set_P(0)
			.set_ePot(U.ePot)
			;
		} else {
			return Prim()
			.set_rho(U.rho)
			.set_v(calc_v(U))
			.set_P(calc_P(solver, U, pt))
			.set_ePot(U.ePot)
			;
		}
	}

	static inline Cons consFromPrim(
		constant Solver const & solver,
		Prim const & W,
		real3 const pt
	) {
		return Cons()
		.set_rho(W.rho)
		.set_m(W.v * W.rho)
		.set_ETotal(calc_ETotal(solver, W, pt))
		.set_ePot(W.ePot)
		;
	}

	// only used by PLM
	static inline Cons apply_dU_dW(
		constant Solver const & solver,
		Prim const & WA,
		Prim const & W,
		real3 const pt
	) {
		Cons result;
		result.rho = W.rho;
		result.m = 
			WA.v * W.rho
			+ W.v * WA.rho;
//// MODULE_DEPENDS: <?=coord_lower?>
		real3 const WA_vL = coord_lower(WA.v, pt);
		result.ETotal = W.rho * (real).5 * dot(WA.v, WA_vL)
			+ WA.rho * dot(W.v, WA_vL)
			+ W.P / (solver.heatCapacityRatio - 1.);
		result.ePot = W.ePot;
		return result;
	}

	// only used by PLM
	static inline Prim apply_dW_dU(
		constant Solver const & solver,
		Prim const & WA,
		Cons const & U,
		real3 const pt
	) {
		Prim result;
		result.rho = U.rho;
		if (U.rho < solver.rhoMin) {
			result.v = {};
			result.P = 0.;
		} else {
			result.v = 
				U.m * (1. / WA.rho)
				- WA.v * (U.rho / WA.rho);
//// MODULE_DEPENDS: <?=coord_lower?>
			real3 const WA_vL = coord_lower(WA.v, pt);
			result.P = (solver.heatCapacityRatio - 1.) * (
				.5 * dot(WA.v, WA_vL) * U.rho 
				- dot(U.m, WA_vL)
				+ U.ETotal);
		}
		result.ePot = U.ePot;
		return result;
	}

	struct ConsWaveCodeMinMaxAllSides {
		real Cs;
		ConsWaveCodeMinMaxAllSides(
			constant Solver const & solver,
			Cons const & U,
			real3 const pt
		) : Cs(calc_Cs_fromCons(solver, U, pt))
		{}

		void operator()(
			constant Solver const & solver,
			Normal n,
			Cons const & U,
			real3 const pt,
			real * resultMin,
			real * resultMax
		) {
			real const Cs_nLen = Cs * normal_len(n);
			real const v_n = U.rho < solver.rhoMin 
				? 0. 
				: normal_vecDotN1(n, U.m) / U.rho;
			//waveCodeAssignMinMax
			if (resultMin) {
				*resultMin = v_n - Cs_nLen;
			}
			if (resultMax) {
				*resultMax = v_n + Cs_nLen;
			}
		}
	};

	/*
	I've highjacked all of this.  It was a normal Euler eqn solver.
	But I experimented with a curved-space solver.
	To get back to the original code,
	just replace all the g_ab stuff with their constant values and simplify away.
	*/
	static inline void applyInitCondCell(
		constant Solver const & solver,
		constant InitCond const & initCond,
		global Cons & U,
		global Cell const & cell
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
		Prim W = Prim()
			.set_rho(rho)
//// MODULE_DEPENDS: <?=cartesianToCoord?>
			.set_v(cartesianToCoord(v, x))
			.set_P(P)
			.set_ePot(ePot)
		;
		U = consFromPrim(solver, W, x);
	}

};

}	//namespace <?=Equation?>

using <?=cons_t?> = <?=Equation?>::Cons;
using <?=prim_t?> = <?=Equation?>::Prim;
using <?=initCond_t?> = <?=Equation?>::InitCond;


// TODO use the one in hydro/eqn/eqn.clcpp
#if 1
// TODO std::forward implementation
kernel void <?=calcDT?>(
	constant <?=solver_t?> const * const psolver,
	global real * const dtBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
<?
if require "hydro.solver.meshsolver":isa(solver) then
?>	,
	global <?=face_t?> const * const faces,
	global int const * const cellFaceIndexes
<?
end
?>
) {
	auto const & solver = *psolver;
	<?=Equation?>::Eqn::calcDT(solver, dtBuf, UBuf, cellBuf
<?
if require "hydro.solver.meshsolver":isa(solver) then
?>, faces, cellFaceIndexes
<?
end
?>	);
}


kernel void <?=applyInitCond?>(
	constant <?=solver_t?> const * const psolver,
	constant <?=initCond_t?> const * const pinitCond,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> * const cellBuf
) {
	auto const & solver = *psolver;
	constant <?=initCond_t?> const & initCond = *pinitCond;
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?> & U = UBuf[index];
	global <?=cell_t?> & cell = cellBuf[index];
	<?=Equation?>::Eqn::applyInitCondCell(solver, initCond, U, cell);
}

#endif


//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=normal_t?> <?=Equation?>

static inline <?=cons_t?> <?=fluxFromCons?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	<?=cell_t?> const & cell,
	<?=normal_t?> const n
) {
	<?=cons_t?> resultF;
	<?=prim_t?> W = <?=Equation?>::Eqn::primFromCons(solver, U, cell.pos);
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
//// MODULE_DEPENDS: real3x3 <?=Equation?>
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
	<?=prim_t?> W = <?=Equation?>::Eqn::primFromCons(solver, U, pt);
	real const v_n = dot(W.v, nL.x);
	real const Cs = <?=Equation?>::Eqn::calc_Cs(solver, &W);
	real const Cs_nLen = Cs * nLen;
	<?=range_t?> result;
	result.min = v_n - Cs_nLen; 
	result.max = v_n + Cs_nLen;
	return result;
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=normal_t?> <?=coord_lower?> <?=cons_t?> <?=prim_t?> <?=eigen_t?> <?=Equation?>
// Equation is for all the calc_* stuff

// used by PLM
static inline <?=eigen_t?> <?=eigen_forCell?>(
	constant <?=solver_t?> const & solver,
	<?=cons_t?> const & U,
	<?=cell_t?> const & cell,
	<?=normal_t?> const n
) {
	<?=prim_t?> W = <?=Equation?>::Eqn::primFromCons(solver, U, cell.pos);
	real3 const vL = coord_lower(W.v, cell.pos);
	real const vSq = dot(W.v, vL);
	real const v_n = normal_vecDotN1(n, W.v);
	real const eKin = .5 * vSq;
	real const hTotal = <?=Equation?>::Eqn::calc_hTotal(W.rho, W.P, U.ETotal);
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
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?> <?=coord_lower?> <?=Equation?>

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
			<?=prim_t?> WR = <?=Equation?>::Eqn::primFromCons(solver, UR, cellR.pos);
			result.rho = UR.rho;
			result.v = WR.v;
			result.vL = coord_lower(WR.v, pt);
			result.vSq = dot(WR.v, result.vL);
			result.hTotal = <?=Equation?>::Eqn::calc_hTotal(WR.rho, WR.P, UR.ETotal);
			result.Cs = <?=Equation?>::Eqn::calc_Cs(solver, WR);
		} else if (UR.rho < rhoEpsilon) {
			// right is vacuum:
			<?=prim_t?> WL = <?=Equation?>::Eqn::primFromCons(solver, UL, cellL.pos);
			result.rho = UL.rho;
			result.v = WL.v;
			result.vL = coord_lower(WL.v, pt);
			result.vSq = dot(WL.v, result.vL);
			result.hTotal = <?=Equation?>::Eqn::calc_hTotal(WL.rho, WL.P, UL.ETotal);
			result.Cs = <?=Equation?>::Eqn::calc_Cs(solver, WL);
		} else {
			<?=prim_t?> WL = <?=Equation?>::Eqn::primFromCons(solver, UL, cellL.pos);
			real const sqrtRhoL = sqrt(WL.rho);
			real3 const vLeft = WL.v;
			real const hTotalL = <?=Equation?>::Eqn::calc_hTotal(WL.rho, WL.P, UL.ETotal);

			<?=prim_t?> WR = <?=Equation?>::Eqn::primFromCons(solver, UR, cellR.pos);
			real const sqrtRhoR = sqrt(WR.rho);
			real3 const vR = WR.v;
			real const hTotalR = <?=Equation?>::Eqn::calc_hTotal(WR.rho, WR.P, UR.ETotal);

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
// Not used anymore.  Was used by Roe, but I switched that to a <?=fluxFromCons?>.
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
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS_NOGHOST?> <?=Equation?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	auto const & solver = *psolver;
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
		real const PL = <?=Equation?>::Eqn::calc_P?>(solver, UL, xL);
		real const PR = <?=Equation?>::Eqn::calc_P?>(solver, UR, xR);
	
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
//// MODULE_DEPENDS: <?=coord_conn_apply23?> <?=coord_conn_trace23?> <?=coord_conn_apply13?>
/*
This is working for init conds with zero velocity.
Introducing constant velocity of v=[x=1,y=1]=[r=sqrt(2),theta=pi/4] in the init cond causes some numerical errors.
However the problem isn't the terms below -- because disabling this for zero-vel init conds causes things to become unsteady.
That means that the volume gradient in calcDerivFV is causing nonzero velocities to emerge, and this is cancelling them.
Maybe for an initial constant vel as large as sqrt(2) this fails, but it works only for small perturbations?
*/
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	<?=prim_t?> W = <?=Equation?>::Eqn::primFromCons(solver, U, x);
	
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
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=Equation?>

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(0,0);
	real3 const x = cellBuf[index].pos;

	global <?=cons_t?> & U = UBuf[index];
	<?=prim_t?> W = <?=Equation?>::Eqn::primFromCons(solver, U, x);

	if (W.rho < solver.rhoMin) W.rho = solver.rhoMin;
	if (W.P < solver.PMin) W.P = solver.PMin;

	U = <?=Equation?>::Eqn::consFromPrim(solver, W, x);
}
<? end ?>
