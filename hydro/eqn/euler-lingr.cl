//// MODULE_NAME: <?=calc_H?>
#define /*real*/ <?=calc_H?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*real const */P\
)	((P) * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)))

//// MODULE_NAME: <?=calc_h?>
#define /*real*/ <?=calc_h?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*real const */rho,\
	/*real const */P\
)	(<?=calc_H?>(solver, P) / (rho))

//// MODULE_NAME: <?=calc_HTotal?>
#define /*real*/ <?=calc_HTotal?>(\
	/*real const */P,\
	/*real const */ETotal\
)	((P) + (ETotal))

//// MODULE_NAME: <?=calc_hTotal?>
#define /*real*/ <?=calc_hTotal?>(\
	/*real const */rho,\
	/*real const */P,\
	/*real const */ETotal\
)	(<?=calc_HTotal?>(P, ETotal) / (rho))

//// MODULE_NAME: <?=calc_eKin?>
//// MODULE_DEPENDS: <?=coordLenSq?>
#define /*real*/ <?=calc_eKin?>(\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
)	(.5 * coordLenSq((W)->v, pt))

//// MODULE_NAME: <?=calc_EKin?>
#define /*real*/ <?=calc_EKin?>(\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) 	((W)->rho * <?=calc_eKin?>(W, pt))

//// MODULE_NAME: <?=calc_EInt?>
#define /*real*/ <?=calc_EInt?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W\
) 	((W)->P / (solver->heatCapacityRatio - 1.))

//// MODULE_NAME: <?=calc_eInt?>
#define /*real*/ <?=calc_eInt?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W\
)	(<?=calc_EInt?>(solver, W) / (W)->rho)

//// MODULE_NAME: <?=calc_EKin_fromCons?>
//// MODULE_DEPENDS: <?=coordLenSq?>
#define /*real*/ <?=calc_EKin_fromCons?>(\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */pt\
)	(.5 * coordLenSq((U)->m, pt) / (U)->rho)

//// MODULE_NAME: <?=calc_ETotal?>
#define /*real*/ <?=calc_ETotal?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) 	(<?=calc_EKin?>(W, pt) + <?=calc_EInt?>(solver, W))

//// MODULE_NAME: <?=calc_Cs?>
#define /*real*/ <?=calc_Cs?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W\
) 	(sqrt(solver->heatCapacityRatio * (W)->P / (W)->rho))

//// MODULE_NAME: <?=calc_P?>
#define /*real*/ <?=calc_P?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */pt\
)	(solver->heatCapacityRatio - 1.) * (/*EInt=*/(U)->ETotal - /*EKin=*/<?=calc_EKin_fromCons?>(U, pt))

//// MODULE_NAME: <?=calc_Cs_fromCons?>
#define /*real*/ <?=calc_Cs_fromCons?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) 	(sqrt(solver->heatCapacityRatio * <?=calc_P?>(solver, U, pt) / (U)->rho))

//// MODULE_NAME: <?=calc_eInt_fromCons?>
//// MODULE_DEPENDS: <?=coordLenSq?>
#define /*real*/ <?=calc_eInt_fromCons?>(\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */pt\
) (((U)->ETotal - .5 * coordLenSq((U)->m, pt) / (U)->rho) / (U)->rho)

//// MODULE_NAME: <?=calc_T?>
<? local materials = require "hydro.materials" ?>
#define C_v				<?=("%.50f"):format(materials.Air.C_v)?>
#define /*real*/ <?=calc_T?>(\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */pt\
) (<?=calc_eInt_fromCons?>(U, pt) / C_v)

//// MODULE_NAME: <?=calc_v?>
#define /*real3*/ <?=calc_v?>(\
	/*<?=cons_t?> const * const*/U\
) (real3_real_mul((U)->m, 1. / (U)->rho))

//// MODULE_NAME: <?=calc_EgField?>
//// MODULE_DEPENDS: units
#define /*real3*/ <?=calc_EgField?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U\
)\
	(real3_real_mul((U)->D_g, 4. * M_PI * solver->gravitationalConstant / unit_m3_per_kg_s2))

//// MODULE_NAME: <?=calc_HgField?>
#define /*real3*/ <?=calc_HgField?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U\
)\
	(real3_real_mul((U)->B_g, solver->speedOfLight * solver->speedOfLight / (4. * M_PI * solver->gravitationalConstant) / unit_kg_per_m))

//// MODULE_NAME: <?=calc_SgField?>
#define /*real3*/ <?=calc_SgField?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U\
)\
	(real3_cross(<?=calc_EgField?>(solver, U), <?=calc_HgField?>(solver, U)))

//// MODULE_NAME: <?=calc_GEM_energy?>
//// MODULE_DEPENDS: <?=coordLenSq?>
//units are energy-per-volume: kg/(m*s^2)
static inline real <?=calc_GEM_energy?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const x
) {
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;
	real const _1_eps_g = 4. * M_PI * G;
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real const _1_mu_g = speedOfLightSq / _1_eps_g;

	return .5 * (coordLenSq(U->D_g, x) * _1_eps_g + coordLenSq(U->B_g, x) * _1_mu_g);
}

//// MODULE_NAME: <?=calcGravityForcePerVolume?>
//rho * (E + v * B) has units kg/(m^2 s^2)
//so this is a force-per-volume
//divide by rho to get m/s^2 acceleration
real3 <?=calcGravityForcePerVolume?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const pt
) {
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;
	real const _1_eps_g = 4. * M_PI * G;
	return _real3(
		U->rho * U->D_g.x * _1_eps_g + U->m.y * U->B_g.z - U->m.z * U->B_g.y,
		U->rho * U->D_g.y * _1_eps_g + U->m.z * U->B_g.x - U->m.x * U->B_g.z,
		U->rho * U->D_g.z * _1_eps_g + U->m.x * U->B_g.y - U->m.y * U->B_g.x);
}

//// MODULE_NAME: <?=primFromCons?>

#define <?=primFromCons?>(\
	/*<?=prim_t?> * const */resultW,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */pt\
) {\
	(resultW)->rho = (U)->rho;\
	(resultW)->v = <?=calc_v?>(U);\
	(resultW)->P = <?=calc_P?>(solver, U, pt);\
	(resultW)->D_g = (U)->D_g;\
	(resultW)->B_g = (U)->B_g;\
	(resultW)->psi_g = (U)->psi_g;\
	(resultW)->phi_g = (U)->phi_g;\
}

//// MODULE_NAME: <?=consFromPrim?>

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */ resultU,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
	(resultU)->rho = (W)->rho;\
	(resultU)->m = real3_real_mul((W)->v, (W)->rho);\
	(resultU)->ETotal = <?=calc_ETotal?>(solver, W, pt);\
	(resultU)->D_g = (W)->D_g;\
	(resultU)->B_g = (W)->B_g;\
	(resultU)->psi_g = (W)->psi_g;\
	(resultU)->phi_g = (W)->phi_g;\
}

//// MODULE_NAME: <?=apply_dU_dW?>
//// MODULE_DEPENDS: <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>
// only used by PLM

#define <?=apply_dU_dW?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */x\
) {\
	real3 const WA_vL = coord_lower((WA)->v, x);\
	(result)->rho = (W)->rho;\
	(result)->m = real3_add(\
		real3_real_mul((WA)->v, (W)->rho),\
		real3_real_mul((W)->v, (WA)->rho));\
	(result)->ETotal = (W)->rho * .5 * real3_dot((WA)->v, WA_vL)\
		+ (WA)->rho * real3_dot((W)->v, WA_vL)\
		+ (W)->P / (solver->heatCapacityRatio - 1.);\
	(result)->B_g = (W)->B_g;\
	(result)->D_g = (W)->D_g;\
	(result)->phi_g = (W)->phi_g;\
	(result)->psi_g = (W)->psi_g;\
}

//// MODULE_NAME: <?=apply_dW_dU?>
//// MODULE_DEPENDS: <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=coord_lower?>
// only used by PLM

#define <?=apply_dW_dU?>(\
	/*<?=prim_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x\
) {\
	real3 const WA_vL = coord_lower((WA)->v, x);\
	(result)->rho = (U)->rho;\
	(result)->v = real3_sub(\
		real3_real_mul((U)->m, 1. / (WA)->rho),\
		real3_real_mul((WA)->v, (U)->rho / (WA)->rho));\
	(result)->P = (solver->heatCapacityRatio - 1.) * (\
		.5 * real3_dot((WA)->v, WA_vL) * (U)->rho \
		- real3_dot((U)->m, WA_vL)\
		+ (U)->ETotal);\
	(result)->B_g = (U)->B_g;\
	(result)->D_g = (U)->D_g;\
	(result)->phi_g = (U)->phi_g;\
	(result)->psi_g = (U)->psi_g;\
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
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(solver->initCondMins, solver->initCondMaxs), .5);
	bool const lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	// these are all standard for all init/euler.lua initial conditions
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	
<?=initCode()?>
	
	<?=prim_t?> W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
		.D_g = real3_zero,
		.B_g = real3_zero,
		.psi_g = 0,
		.phi_g = 0,
	};

	<?=consFromPrim?>(U, solver, &W, x);
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=primFromCons?> <?=normal_t?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, (cell)->pos);\
	real const v_n = normal_vecDotN1(n, W.v);\
	(resultFlux)->rho = (U)->rho * v_n;\
	(resultFlux)->m = real3_add(\
		real3_real_mul((U)->m, v_n),\
		real3_real_mul(normal_u1(n), W.P)\
	);\
	real const HTotal = (U)->ETotal + W.P;\
	(resultFlux)->ETotal = HTotal * v_n;\
\
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;\
	real const _1_eps_g = 4. * M_PI * G;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const _1_mu_g = speedOfLightSq / _1_eps_g;\
\
	/* taken from glm-maxwell instead of the 2014 Abgrall, Kumar */\
	/*  then replace D = ε E and φ' -> ε φ */\
	real3 const E_g = real3_real_mul((U)->D_g, _1_eps_g);\
	real3 const H_g = real3_real_mul((U)->B_g, _1_mu_g);\
\
	real const divPhiWavespeed_g = solver->divPhiWavespeed_g / unit_m_per_s;\
	real const divPsiWavespeed_g = solver->divPsiWavespeed_g / unit_m_per_s;\
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
\
	(resultFlux)->D_g.x = H_g.y * nz - H_g.z * ny + nx * (U)->phi_g * divPhiWavespeed_g;	/* F_D^i = -ε^ijk n_j H_k */\
	(resultFlux)->B_g.x = E_g.z * ny - E_g.y * nz + nx * (U)->psi_g * divPsiWavespeed_g;	/* F_B^i = +ε^ijk n_j B_k */\
\
	(resultFlux)->D_g.y = H_g.z * nx - H_g.x * nz + ny * (U)->phi_g * divPhiWavespeed_g;\
	(resultFlux)->B_g.y = E_g.x * nz - E_g.z * nx + ny * (U)->psi_g * divPsiWavespeed_g;\
\
	(resultFlux)->D_g.z = H_g.x * ny - H_g.y * nx + nz * (U)->phi_g * divPhiWavespeed_g;\
	(resultFlux)->B_g.z = E_g.y * nx - E_g.x * ny + nz * (U)->psi_g * divPsiWavespeed_g;\
\
	real const D_n = normal_vecDotN1(n, (U)->D_g);\
	real const B_n = normal_vecDotN1(n, (U)->B_g);\
	(resultFlux)->phi_g = D_n * divPhiWavespeed_g;\
	(resultFlux)->psi_g = B_n * divPsiWavespeed_g;\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: real3x3 <?=primFromCons?>
// added by request only, so I don't have to compile the real3x3 code.
// not used at the moment

#define <?=calcCellMinMaxEigenvalues?>(\
	/*<?=range_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*real3x3 const */nL,\
	/*real3x3 const */nU,\
	/*real const */nLen\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
	real const v_n = real3_dot(W.v, nL.x);\
	real const Cs = <?=calc_Cs?>(solver, &W);\
	real const Cs_nLen = Cs * nLen;\
	(result)->min = v_n - Cs_nLen; \
	(result)->max = v_n + Cs_nLen;\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=coord_lower?>

// used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, (cell)->pos);\
	real3 const vL = coord_lower(W.v, (cell)->pos);\
	real const vSq = real3_dot(W.v, vL);\
	real const v_n = normal_vecDotN1(n, W.v);\
	real const eKin = .5 * vSq;\
	real const hTotal = <?=calc_hTotal?>(W.rho, W.P, (U)->ETotal);\
	real const CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);\
	real const Cs = sqrt(CsSq);\
	(result)->rho = W.rho;\
	(result)->v = W.v;\
	(result)->vSq = vSq;\
	(result)->vL = vL;\
	(result)->hTotal = hTotal;\
	(result)->Cs = Cs;\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=primFromCons?> <?=eigen_t?> <?=normal_t?> <?=coord_lower?>

//used by the mesh version
#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> WL;\
	<?=primFromCons?>(&WL, solver, UL, (cellL)->pos);\
	real const sqrtRhoL = sqrt(WL.rho);\
	real3 const vLeft = WL.v;\
	real const hTotalL = <?=calc_hTotal?>(WL.rho, WL.P, (UL)->ETotal);\
\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, (cellR)->pos);\
	real const sqrtRhoR = sqrt(WR.rho);\
	real3 const vR = WR.v;\
	real const hTotalR = <?=calc_hTotal?>(WR.rho, WR.P, (UR)->ETotal);\
\
	real const invDenom = 1./(sqrtRhoL + sqrtRhoR);\
\
	/*Roe-averaged*/\
	real const rho = sqrtRhoL * sqrtRhoR;\
	real3 const v = real3_add(\
			real3_real_mul(vLeft, sqrtRhoL * invDenom),\
			real3_real_mul(vR, sqrtRhoR * invDenom));\
	real const hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);\
\
	/*derived:*/\
	real3 const vLower = coord_lower(v, pt);\
	real const vSq = real3_dot(v, vLower);\
	real const eKin = .5 * vSq;\
	real const CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);\
	real const Cs = sqrt(CsSq);\
\
	(result)->rho = rho;\
	(result)->v = v;\
	(result)->vSq = vSq;\
	(result)->vL = vLower;\
	(result)->hTotal = hTotal;\
	(result)->Cs = Cs;\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> */n\
) {\
	real3 const v_n = normal_vecDotNs(n, (eig)->v);\
	real const inv_nLen = 1. / normal_len(n);\
	real const denom = 2. * (eig)->Cs * (eig)->Cs;\
	real const invDenom = 1. / denom;\
	real const gamma_1 = solver->heatCapacityRatio - 1.;\
\
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;\
	real const _1_eps_g = 4. * M_PI * G;\
	real const sqrt_1_eps_g = sqrt(_1_eps_g);\
	real const sqrt_1_mu_g = solver->speedOfLight / sqrt_1_eps_g;\
	real const _1_mu_g = sqrt_1_mu_g * sqrt_1_mu_g;\
\
	(Y)->ptr[0] = (\
			(X)->ptr[0] * (.5 * gamma_1 * (eig)->vSq + (eig)->Cs * v_n.x * inv_nLen)\
			+ (X)->ptr[1] * (-gamma_1 * (eig)->vL.x - (eig)->Cs * normal_l1x_over_len(n))\
			+ (X)->ptr[2] * (-gamma_1 * (eig)->vL.y - (eig)->Cs * normal_l1y_over_len(n))\
			+ (X)->ptr[3] * (-gamma_1 * (eig)->vL.z - (eig)->Cs * normal_l1z_over_len(n))\
			+ (X)->ptr[4] * gamma_1\
		) * invDenom;\
	(Y)->ptr[1] =\
		(\
			(X)->ptr[0] * (denom - gamma_1 * (eig)->vSq)\
			+ (X)->ptr[1] * 2. * gamma_1 * (eig)->vL.x\
			+ (X)->ptr[2] * 2. * gamma_1 * (eig)->vL.y\
			+ (X)->ptr[3] * 2. * gamma_1 * (eig)->vL.z\
			+ (X)->ptr[4] * -2. * gamma_1\
		) * invDenom;\
	(Y)->ptr[2] =\
		(X)->ptr[0] * -v_n.y\
		+ (X)->ptr[1] * normal_l2x(n)\
		+ (X)->ptr[2] * normal_l2y(n)\
		+ (X)->ptr[3] * normal_l2z(n);\
	(Y)->ptr[3] =\
		(X)->ptr[0] * -v_n.z\
		+ (X)->ptr[1] * normal_l3x(n)\
		+ (X)->ptr[2] * normal_l3y(n)\
		+ (X)->ptr[3] * normal_l3z(n);\
	(Y)->ptr[4] =\
		(\
			(X)->ptr[0] * (.5 * gamma_1 * (eig)->vSq - (eig)->Cs * v_n.x * inv_nLen)\
			+ (X)->ptr[1] * (-gamma_1 * (eig)->vL.x + (eig)->Cs * normal_l1x_over_len(n))\
			+ (X)->ptr[2] * (-gamma_1 * (eig)->vL.y + (eig)->Cs * normal_l1y_over_len(n))\
			+ (X)->ptr[3] * (-gamma_1 * (eig)->vL.z + (eig)->Cs * normal_l1z_over_len(n))\
			+ (X)->ptr[4] * gamma_1\
		) * invDenom;\
\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = 1. / sqrt_1_eps_g;\
	real const tmp2 = tmp1 * n_l.x;\
	real const tmp3 = 1. / 2.;\
	real const tmp4 = (X)->ptr[5+0] * tmp2;\
	real const tmp5 = tmp3 * tmp4;\
	real const tmp7 = n_l.y * tmp1;\
	real const tmp9 = (X)->ptr[5+1] * tmp7;\
	real const tmp11 = tmp9 * tmp3;\
	real const tmp13 = n_l.z * tmp1;\
	real const tmp15 = (X)->ptr[5+2] * tmp13;\
	real const tmp17 = tmp15 * tmp3;\
	real const tmp20 = (X)->ptr[5+6] * tmp1;\
	real const tmp22 = tmp20 * tmp3;\
	real const tmp26 = (X)->ptr[5+3] * tmp2;\
	real const tmp27 = tmp26 * tmp3;\
	real const tmp31 = (X)->ptr[5+4] * tmp7;\
	real const tmp33 = tmp31 * tmp3;\
	real const tmp37 = (X)->ptr[5+5] * tmp13;\
	real const tmp39 = tmp37 * tmp3;\
	real const tmp42 = (X)->ptr[5+7] * tmp1;\
	real const tmp44 = tmp42 * tmp3;\
	real const tmp45 = sqrt_1_mu_g * n2_l.z;\
	real const tmp47 = (X)->ptr[5+5] * tmp45;\
	real const tmp48 = sqrt_1_mu_g * n2_l.y;\
	real const tmp50 = (X)->ptr[5+4] * tmp48;\
	real const tmp51 = sqrt_1_mu_g * n2_l.x;\
	real const tmp53 = (X)->ptr[5+3] * tmp51;\
	real const tmp54 = sqrt_1_eps_g * n3_l.z;\
	real const tmp56 = (X)->ptr[5+2] * tmp54;\
	real const tmp57 = sqrt_1_eps_g * n3_l.y;\
	real const tmp59 = (X)->ptr[5+1] * tmp57;\
	real const tmp60 = sqrt_1_eps_g * n3_l.x;\
	real const tmp62 = (X)->ptr[5+0] * tmp60;\
	real const tmp63 = tmp59 * tmp3;\
	real const tmp64 = tmp62 * tmp3;\
	real const tmp65 = tmp56 * tmp3;\
	real const tmp67 = tmp53 * tmp3;\
	real const tmp69 = tmp50 * tmp3;\
	real const tmp71 = tmp47 * tmp3;\
	real const tmp73 = sqrt_1_eps_g * n2_l.x;\
	real const tmp75 = (X)->ptr[5+0] * tmp73;\
	real const tmp76 = tmp75 * tmp3;\
	real const tmp77 = sqrt_1_eps_g * n2_l.y;\
	real const tmp79 = (X)->ptr[5+1] * tmp77;\
	real const tmp81 = tmp79 * tmp3;\
	real const tmp82 = sqrt_1_eps_g * n2_l.z;\
	real const tmp84 = (X)->ptr[5+2] * tmp82;\
	real const tmp86 = tmp84 * tmp3;\
	real const tmp87 = sqrt_1_mu_g * n3_l.x;\
	real const tmp89 = (X)->ptr[5+3] * tmp87;\
	real const tmp90 = sqrt_1_mu_g * n3_l.y;\
	real const tmp92 = (X)->ptr[5+4] * tmp90;\
	real const tmp93 = sqrt_1_mu_g * n3_l.z;\
	real const tmp95 = (X)->ptr[5+5] * tmp93;\
	real const tmp96 = tmp92 * tmp3;\
	real const tmp97 = tmp95 * tmp3;\
	real const tmp98 = tmp89 * tmp3;\
	(Y)->ptr[5+0] = tmp22 + -tmp5 - tmp11 - tmp17;\
	(Y)->ptr[5+1] = tmp44 + -tmp27 - tmp33 - tmp39;\
	(Y)->ptr[5+2] = tmp71 + tmp69 + tmp67 + tmp65 + tmp63 + tmp64;\
	(Y)->ptr[5+3] = tmp98 + tmp96 + tmp97 + -tmp76 - tmp81 - tmp86;\
	(Y)->ptr[5+4] = tmp67 + tmp69 + tmp71 + -tmp64 - tmp63 - tmp65;\
	(Y)->ptr[5+5] = tmp97 + tmp96 + tmp98 + tmp86 + tmp81 + tmp76;\
	(Y)->ptr[5+6] = tmp22 + tmp17 + tmp11 + tmp5;\
	(Y)->ptr[5+7] = tmp44 + tmp39 + tmp33 + tmp27;\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */Y,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real3 const v_n = normal_vecDotNs(n, (eig)->v);\
	real const inv_nLen = 1. / normal_len(n);\
\
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;\
	real const _1_eps_g = 4. * M_PI * G;\
	real const sqrt_1_eps_g = sqrt(_1_eps_g);\
	real const sqrt_1_mu_g = solver->speedOfLight / sqrt_1_eps_g;\
\
	(Y)->ptr[0] =\
		(X)->ptr[0]\
		+ (X)->ptr[1]\
		+ (X)->ptr[4];\
	(Y)->ptr[1] =\
		(X)->ptr[0] * ((eig)->v.x - (eig)->Cs * normal_u1x_over_len(n))\
		+ (X)->ptr[1] * (eig)->v.x\
		+ (X)->ptr[2] * normal_u2x(n)\
		+ (X)->ptr[3] * normal_u3x(n)\
		+ (X)->ptr[4] * ((eig)->v.x + (eig)->Cs * normal_u1x_over_len(n));\
	(Y)->ptr[2] =\
		(X)->ptr[0] * ((eig)->v.y - (eig)->Cs * normal_u1y_over_len(n))\
		+ (X)->ptr[1] * (eig)->v.y\
		+ (X)->ptr[2] * normal_u2y(n)\
		+ (X)->ptr[3] * normal_u3y(n)\
		+ (X)->ptr[4] * ((eig)->v.y + (eig)->Cs * normal_u1y_over_len(n));\
	(Y)->ptr[3] =\
		(X)->ptr[0] * ((eig)->v.z - (eig)->Cs * normal_u1z_over_len(n))\
		+ (X)->ptr[1] * (eig)->v.z\
		+ (X)->ptr[2] * normal_u2z(n)\
		+ (X)->ptr[3] * normal_u3z(n)\
		+ (X)->ptr[4] * ((eig)->v.z + (eig)->Cs * normal_u1z_over_len(n));\
	(Y)->ptr[4] =\
		(X)->ptr[0] * ((eig)->hTotal - (eig)->Cs * v_n.x * inv_nLen)\
		+ (X)->ptr[1] * .5 * (eig)->vSq\
		+ (X)->ptr[2] * v_n.y\
		+ (X)->ptr[3] * v_n.z\
		+ (X)->ptr[4] * ((eig)->hTotal + (eig)->Cs * v_n.x * inv_nLen);\
\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = sqrt_1_eps_g * n_l.x;\
	real const tmp3 = 1. / sqrt_1_eps_g;\
	real const tmp4 = tmp3 * n3_l.x;\
	real const tmp6 = n2_l.x * tmp3;\
	real const tmp22 = sqrt_1_eps_g * n_l.y;\
	real const tmp25 = n3_l.y * tmp3;\
	real const tmp27 = n2_l.y * tmp3;\
	real const tmp43 = sqrt_1_eps_g * n_l.z;\
	real const tmp46 = n3_l.z * tmp3;\
	real const tmp48 = n2_l.z * tmp3;\
	real const tmp66 = 1. / sqrt_1_mu_g;\
	real const tmp67 = tmp66 * n2_l.x;\
	real const tmp69 = n3_l.x * tmp66;\
	real const tmp88 = n2_l.y * tmp66;\
	real const tmp90 = n3_l.y * tmp66;\
	real const tmp109 = n2_l.z * tmp66;\
	real const tmp111 = n3_l.z * tmp66;\
	(Y)->ptr[5+0] = (X)->ptr[5+6] * tmp1 + (X)->ptr[5+5] * tmp6 + (X)->ptr[5+2] * tmp4 - (X)->ptr[5+3] * tmp6 - (X)->ptr[5+4] * tmp4 + -(X)->ptr[5+0] * tmp1;\
	(Y)->ptr[5+1] = (X)->ptr[5+6] * tmp22 + (X)->ptr[5+5] * tmp27 + (X)->ptr[5+2] * tmp25 - (X)->ptr[5+3] * tmp27 - (X)->ptr[5+4] * tmp25 + -(X)->ptr[5+0] * tmp22;\
	(Y)->ptr[5+2] = (X)->ptr[5+6] * tmp43 + (X)->ptr[5+5] * tmp48 + (X)->ptr[5+2] * tmp46 - (X)->ptr[5+3] * tmp48 - (X)->ptr[5+4] * tmp46 + -(X)->ptr[5+0] * tmp43;\
	(Y)->ptr[5+3] = (X)->ptr[5+7] * tmp1 + (X)->ptr[5+5] * tmp69 + (X)->ptr[5+4] * tmp67 + (X)->ptr[5+3] * tmp69 + (X)->ptr[5+2] * tmp67 + -(X)->ptr[5+1] * tmp1;\
	(Y)->ptr[5+4] = (X)->ptr[5+7] * tmp22 + (X)->ptr[5+5] * tmp90 + (X)->ptr[5+4] * tmp88 + (X)->ptr[5+3] * tmp90 + (X)->ptr[5+2] * tmp88 + -(X)->ptr[5+1] * tmp22;\
	(Y)->ptr[5+5] = (X)->ptr[5+7] * tmp43 + (X)->ptr[5+5] * tmp111 + (X)->ptr[5+4] * tmp109 + (X)->ptr[5+3] * tmp111 + (X)->ptr[5+2] * tmp109 + -(X)->ptr[5+1] * tmp43;\
	(Y)->ptr[5+6] = (X)->ptr[5+6] * sqrt_1_eps_g + (X)->ptr[5+0] * sqrt_1_eps_g;\
	(Y)->ptr[5+7] = (X)->ptr[5+7] * sqrt_1_eps_g + (X)->ptr[5+1] * sqrt_1_eps_g;\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>
// Not used anymore.  was used by Roe, but I switched that to a <?=fluxFromCons?>.
// <?=fluxFromCons?> only matches <?=eigen_fluxTransform?> when the eig properties are derived from X_

#error Getting rid of this, right?
#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real3 const v_n = normal_vecDotNs(n, (eig).v);\
	real const gamma = solver->heatCapacityRatio;\
	real const gamma_1 = gamma - 1.;\
	real const gamma_2 = gamma - 2.;\
\
	(resultFlux)->ptr[0] =\
		(X)->ptr[1] * normal_l1x(n)\
		+ (X)->ptr[2] * normal_l1y(n)\
		+ (X)->ptr[3] * normal_l1z(n);\
\
	(resultFlux)->ptr[1] =\
		(X)->ptr[0] * (-v_n.x * (eig).v.x + gamma_1 * .5 * (eig).vSq * normal_u1x(n))\
		+ (X)->ptr[1] * ((eig).v.x * normal_l1x(n) - gamma_2 * normal_u1x(n) * (eig).vL.x + v_n.x)\
		+ (X)->ptr[2] * ((eig).v.x * normal_l1y(n) - gamma_2 * normal_u1x(n) * (eig).vL.y)\
		+ (X)->ptr[3] * ((eig).v.x * normal_l1z(n) - gamma_2 * normal_u1x(n) * (eig).vL.z)\
		+ (X)->ptr[4] * gamma_1 * normal_u1x(n);\
\
	(resultFlux)->ptr[2] =\
		(X)->ptr[0] * (-v_n.x * (eig).v.y + gamma_1 * .5 * (eig).vSq * normal_u1y(n))\
		+ (X)->ptr[1] * ((eig).v.y * normal_l1x(n) - gamma_2 * normal_u1y(n) * (eig).vL.x)\
		+ (X)->ptr[2] * ((eig).v.y * normal_l1y(n) - gamma_2 * normal_u1y(n) * (eig).vL.y + v_n.x)\
		+ (X)->ptr[3] * ((eig).v.y * normal_l1z(n) - gamma_2 * normal_u1y(n) * (eig).vL.z)\
		+ (X)->ptr[4] * gamma_1 * normal_u1y(n);\
\
	(resultFlux)->ptr[3] =\
		(X)->ptr[0] * (-v_n.x * (eig).v.z + gamma_1 * .5 * (eig).vSq * normal_u1z(n))\
		+ (X)->ptr[1] * ((eig).v.z * normal_l1x(n) - gamma_2 * normal_u1z(n) * (eig).vL.x)\
		+ (X)->ptr[2] * ((eig).v.z * normal_l1y(n) - gamma_2 * normal_u1z(n) * (eig).vL.y)\
		+ (X)->ptr[3] * ((eig).v.z * normal_l1z(n) - gamma_2 * normal_u1z(n) * (eig).vL.z + v_n.x)\
		+ (X)->ptr[4] * gamma_1 * normal_u1z(n);\
\
	(resultFlux)->ptr[4] =\
		(X)->ptr[0] * v_n.x * (.5 * gamma_1 * (eig).vSq - (eig).hTotal)\
		+ (X)->ptr[1] * (normal_l1x(n) * (eig).hTotal - gamma_1 * v_n.x * (eig).vL.x)\
		+ (X)->ptr[2] * (normal_l1y(n) * (eig).hTotal - gamma_1 * v_n.x * (eig).vL.y)\
		+ (X)->ptr[3] * (normal_l1z(n) * (eig).hTotal - gamma_1 * v_n.x * (eig).vL.z)\
		+ (X)->ptr[4] * gamma * v_n.x;\
\
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const _1_eps_g = 4. * M_PI * G;\
	real const _1_mu_g = speedOfLightSq / _1_eps_g;\
\
	real3 const E_g = real3_real_mul((X)->D_g, _1_eps_g);\
	real3 const H_g = real3_real_mul((X)->B_g, _1_mu_g);\
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
\
	(resultFlux)->D_g.x = H_g.y * nz - H_g.z * ny + nx * (U)->phi_g * solver->divPhiWavespeed_g / unit_m_per_s;	/* F_D^i = -ε^ijk n_j H_k */\
	(resultFlux)->B_g.x = E_g.z * ny - E_g.y * nz + nx * (U)->psi_g * solver->divPsiWavespeed_g / unit_m_per_s;	/* F_B^i = +ε^ijk n_j B_k */\
\
	(resultFlux)->D_g.y = H_g.z * nx - H_g.x * nz + ny * (U)->phi_g * solver->divPhiWavespeed_g / unit_m_per_s;\
	(resultFlux)->B_g.y = E_g.x * nz - E_g.z * nx + ny * (U)->psi_g * solver->divPsiWavespeed_g / unit_m_per_s;\
\
	(resultFlux)->D_g.z = H_g.x * ny - H_g.y * nx + nz * (U)->phi_g * solver->divPhiWavespeed_g / unit_m_per_s;\
	(resultFlux)->B_g.z = E_g.y * nx - E_g.x * ny + nz * (U)->psi_g * solver->divPsiWavespeed_g / unit_m_per_s;\
\
	real const D_n = normal_vecDotN1(n, (U)->D_g);\
	real const B_n = normal_vecDotN1(n, (U)->B_g);\
	(resultFlux)->phi_g = D_n * solver->divPhiWavespeed_g / unit_m_per_s;\
	(resultFlux)->psi_g = B_n * solver->divPsiWavespeed_g / unit_m_per_s;\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS_NOGHOST?> <?=primFromCons?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();

	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);

<? if false
and solver.coord.vectorComponent == "anholonomic"
and require "hydro.coord.cylinder":isa(solver.coord)
then ?>
<? 	if true then -- 2009 Trangenstein, p.474, 1999 Toro, p.29, eqn.1.104, 1.105 ?>
	<? for side=0,1 do ?>{
		real3 xL = x; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		real3 xR = x; xR.s<?=side?> += solver->grid_dx.s<?=side?>;
		
		global <?=cons_t?> const * const UL = U - solver->stepsize.s<?=side?>;
		global <?=cons_t?> const * const UR = U + solver->stepsize.s<?=side?>;
		real const PL = <?=calc_P?>(solver, UL, xL);
		real const PR = <?=calc_P?>(solver, UR, xR);
	
		deriv->m.s<?=side?> -= (PR - PL) / (2. * solver->grid_dx.s<?=side?>);
	}<? end ?>
<?	end ?>
<?	if false then -- 1999 Toro p.28 eqn.1.102, 1.103 ?>
	<?=cons_t?> F;
	fluxFromCons(&F, solver, U, cell, normal_forSide0(x));
	deriv->rho -= F.rho / x.x;
	deriv->m.x -= F.m.x / x.x;
	deriv->m.y -= F.m.y / x.x;
	deriv->m.z -= F.m.z / x.x;
	deriv->ETotal -= F.ETotal / x.x;
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
	
	//- Γ^i_jk ρ v^j v^k
	deriv->m = real3_sub(deriv->m, coord_conn_apply23(W.v, U->m, x));
	
	//- Γ^i_jk g^jk P
	deriv->m = real3_sub(deriv->m, real3_real_mul(coord_conn_trace23(x), W.P));
	
	//+ (γ-1) ρ v^k v^l Γ_kjl g^ij
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_conn_apply13(W.v, U->m, x), (solver->heatCapacityRatio - 1.) ));
	
	//- (γ-1) ρ v^j v^k v^l Γ_jkl
//	deriv->ETotal -= (solver->heatCapacityRatio - 1.) * coord_conn_apply123(W.v, W.v, U->m, x);

	//+ c_jk^k * Flux^Ij
<? 	if false and solver.coord.vectorComponent == "anholonomic" then ?>
	real3 const commTrace = coord_tr23_c(x);
	<? for i=0,solver.dim-1 do ?>{
		<?=cons_t?> flux;
		<?=calcFluxFromCons?>(&F, *U, x);
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] += commTrace.s<?=i?> * flux.ptr[j];
		}
	}<? end ?>
<? 	end ?>
<? end ?>
<? end -- vectorComponent == "anholonomic" ?>


//	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;
//	real const _1_eps_g = 4. * M_PI * G;
//	real const mu_g = _1_eps_g / speedOfLightSq;
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;

	real3 const gravityForce = <?=calcGravityForcePerVolume?>(solver, U, x);
	
	deriv->m = real3_add(deriv->m, gravityForce);

	// u is unitless, ~= [1, v/c]
	//T_ab = (c^2 ρ + P) u_a u_b + P g_ab
	//T_00 = c^2 ρ + P
	//T_0i = (c ρ + P/c) v_i
	//T_ij = (c^2 ρ + P) v_i v_j + P δ_ij

	// J_i = 1/c T_0i = (ρ + P/c^2) v_i
	// [T_0i] = kg/(m*s^2)
	// [J_i] = kg/(m^2*s)

	real const rho_plus_P_over_c2 = 0.
		/* matter */
		+ U->rho					/* kg/m^3 */
		+ W.P / speedOfLightSq 	/* kg/m^3 */
	;

	/* source of D_g is J_g is the momentum + Poynting vector */
	real3 J_g = real3_zero;
	/*  I'm symmetrizing the stress-energy */
		/* matter */
	J_g.x += W.v.x * rho_plus_P_over_c2;	/* [m] = kg/(m^2*s) */
	J_g.y += W.v.y * rho_plus_P_over_c2;
	J_g.z += W.v.z * rho_plus_P_over_c2;
	
	/* D_g's units: kg/m^2 */
	deriv->D_g.x += J_g.x;	/* ∫ [J_g] dt = kg/m^2 */
	deriv->D_g.y += J_g.y;
	deriv->D_g.z += J_g.z;
	
	/*  source of φ_g is T_00 is ρ + .5 (D^2 + B^2) */
	real const T_00_over_c2 = rho_plus_P_over_c2;
	deriv->phi_g -= 								/* = kg/(m^2*s) */
		T_00_over_c2 								/* kg/m^3 */
		* solver->divPhiWavespeed_g / unit_m_per_s	/* m/s */
	;	/* and the deriv is ∂/∂t, so it is integrated by seconds */
	/* after integrating the derivative wrt t it becomes in units of kg/m^2, which is φ_g's units */
}

//// MODULE_NAME: <?=constrainU?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=primFromCons?> <?=consFromPrim?>

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(0,0);
	real3 const x = cellBuf[index].pos;

	global <?=cons_t?> * const U = UBuf + index;
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);

	if (W.rho < solver->rhoMin) W.rho = solver->rhoMin;
	if (W.P < solver->PMin) W.P = solver->PMin;

	<?=consFromPrim?>(U, solver, &W, x);
}
