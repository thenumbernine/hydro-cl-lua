//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: <?=coordLenSq?> <?=cons_t?> <?=prim_t?> <?=waves_t?> <?=eigen_t?> <?=eqn_guiVars_compileTime?>

#define /*real*/ calc_H(\
	/*constant <?=solver_t?> const * const */solver,\
	/*real const */P\
)	((P) * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)))

#define /*real*/ calc_h(\
	/*constant <?=solver_t?> const * const */solver,\
	/*real const */rho,\
	/*real const */P\
)	(calc_H(solver, P) / (rho))

#define /*real*/ calc_HTotal(\
	/*real const */P,\
	/*real const */ETotal\
)	((P) + (ETotal))

#define /*real*/ calc_hTotal(\
	/*real const */rho,\
	/*real const */P,\
	/*real const */ETotal\
)	(calc_HTotal(P, ETotal) / (rho))

#define /*real*/ calc_eKin(\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */x\
)	(.5 * coordLenSq((W)->v, x))

#define /*real*/ calc_EKin(\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */x\
) 	((W)->rho * calc_eKin(W, x))

#define /*real*/ calc_EInt(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W\
) 	((W)->P / (solver->heatCapacityRatio - 1.))

#define /*real*/ calc_eInt(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W\
)	(calc_EInt(solver, W) / (W)->rho)

#define /*real*/ calc_EKin_fromCons(\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */x\
)	(.5 * coordLenSq((U)->m, x) / (U)->rho)

#define /*real*/ calc_ETotal(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */x\
) 	(calc_EKin(W, x) + calc_EInt(solver, W))

#define /*real*/ calc_Cs(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W\
) 	(sqrt(solver->heatCapacityRatio * (W)->P / (W)->rho))

#define /*real*/ calc_P(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */x\
)	(solver->heatCapacityRatio - 1.) * (/*EInt=*/(U)->ETotal - /*EKin=*/calc_EKin_fromCons(U, x))

#define /*real*/ calc_eInt_from_cons(\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */x\
) (((U)->ETotal - .5 * coordLenSq((U)->m, x) / (U)->rho) / (U)->rho)

<? local materials = require "hydro.materials" ?>
#define C_v				<?=("%.50f"):format(materials.Air.C_v)?>

#define /*real*/ calc_T(\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */x\
) (calc_eInt_from_cons(U, x) / C_v)

#define /*real3*/ calc_v(\
	/*<?=cons_t?> const * const*/U\
) (real3_real_mul((U)->m, 1. / (U)->rho))

//// MODULE_DEPENDS: units

// rho * (E + v * B) has units kg/(m^2 s^2)
// so this isn't a force, it's a kg-times-force
real3 calcGravForce(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const x
) {
	real const eps_g = 1. / (4. * M_PI * (solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2));
	return _real3(
		U->rho * U->D_g.x / eps_g + 4. * (U->m.y * U->B_g.z - U->m.z * U->B_g.y),
		U->rho * U->D_g.y / eps_g + 4. * (U->m.z * U->B_g.x - U->m.x * U->B_g.z),
		U->rho * U->D_g.z / eps_g + 4. * (U->m.x * U->B_g.y - U->m.y * U->B_g.x));
}

//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?> <?=eqn_common?>
// eqn_common is for all the calc_* stuff

#define <?=primFromCons?>(\
	/*<?=prim_t?> * const */resultW,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const*/U,\
	/*real3 const */x\
) {\
	(resultW)->rho = (U)->rho;\
	(resultW)->v = calc_v(U);\
	(resultW)->P = calc_P(solver, U, x);\
	(resultW)->D_g = (U)->D_g;\
	(resultW)->B_g = (U)->B_g;\
	(resultW)->psi_g = (U)->psi_g;\
	(resultW)->phi_g = (U)->phi_g;\
}

//// MODULE_NAME: <?=consFromPrim?>
//// MODULE_DEPENDS: real3 <?=eqn_common?> <?=solver_t?> <?=prim_t?> <?=cons_t?>
// eqn_common is for all the calc_* stuff

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */ resultU,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
	(resultU)->rho = (W)->rho;\
	(resultU)->m = real3_real_mul((W)->v, (W)->rho);\
	(resultU)->ETotal = calc_ETotal(solver, W, pt);\
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
	real const eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;\
	real const mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
	/* taken from glm-maxwell instead of the 2014 Abgrall, Kumar */\
	/*  then replace D = epsilon E and phi' -> epsilon phi */\
	real3 const E = real3_real_mul((U)->D_g, 1. / eps_g);\
	real3 const H = real3_real_mul((U)->B_g, 1. / mu_g);\
\
	real const nx = normal_l1x(n);\
	real const ny = normal_l1y(n);\
	real const nz = normal_l1z(n);\
\
	(resultFlux)->D_g.x = H_g.y * nz - H_g.z * ny + nx * (U)->phi_g * solver->divPhiWavespeed_g / unit_m_per_s;	/* F_D^i = -eps^ijk n_j H_k */\
	(resultFlux)->B_g.x = E_g.z * ny - E_g.y * nz + nx * (U)->psi_g * solver->divPsiWavespeed_g / unit_m_per_s;	/* F_B^i = +eps^ijk n_j B_k */\
\
	(resultFlux)->D_g.y = H_g.z * nx - H_g.x * nz + ny * (U)->phi_g * solver->divPhiWavespeed_g / unit_m_per_s;\
	(resultFlux)->B_g.y = E_g.x * nz - E_g.z * nx + ny * (U)->psi_g * solver->divPsiWavespeed_g / unit_m_per_s;\
\
	(resultFlux)->D_g.z = H_g.x * ny - H_g.y * nx + nz * (U)->phi_g * solver->divPhiWavespeed_g / unit_m_per_s;\
	(resultFlux)->B_g.z = E_g.y * nx - E_g.x * ny + nz * (U)->psi_g * solver->divPsiWavespeed_g / unit_m_per_s;\
\
	real const D_n = normal_vecDotN1(n, (U)->D_g);\
	real const B_n = normal_vecDotN1(n, (U)->B_g);\
	(resultFlux)->phi_g = <?=real_mul?>(D_n, solver->divPhiWavespeed_g / unit_m_per_s);\
	(resultFlux)->psi_g = <?=real_mul?>(B_n, solver->divPsiWavespeed_g / unit_m_per_s);\
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
	real const Cs = calc_Cs(solver, &W);\
	real const Cs_nLen = Cs * nLen;\
	(result)->min = v_n - Cs_nLen; \
	(result)->max = v_n + Cs_nLen;\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=normal_t?> <?=coord_lower?> <?=cons_t?> <?=prim_t?> <?=eigen_t?> <?=primFromCons?> <?=eqn_common?>
// eqn_common is for all the calc_* stuff

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
	real const hTotal = calc_hTotal(W.rho, W.P, (U)->ETotal);\
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
	real const hTotalL = calc_hTotal(WL.rho, WL.P, (UL)->ETotal);\
\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, (cellR)->pos);\
	real const sqrtRhoR = sqrt(WR.rho);\
	real3 const vR = WR.v;\
	real const hTotalR = calc_hTotal(WR.rho, WR.P, (UR)->ETotal);\
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
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
	real const denom = 2. * (eig)->Cs * (eig)->Cs;\
	real const invDenom = 1. / denom;\
	real const gamma_1 = solver->heatCapacityRatio - 1.;\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
	real const sqrt_eps_g = sqrt(eps_g);	/* TODO sqrt units */\
	real const sqrt_mu_g = sqrt(mu_g);\
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
	real const sqrt_1_eps_g = 1. / sqrt_eps_g;\
	real const sqrt_1_mu_g = 1. / sqrt_mu_g;\
\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	real const tmp1 = 1. / sqrt_1_eps_g;\
	real const tmp2 = n_l.x * (X)->ptr[5+0];\
	real const tmp3 = 1. / 2.;\
	real const tmp4 = tmp1 * tmp2;\
	real const tmp5 = tmp3 * tmp4;\
	real const tmp7 = n2_l.x * (X)->ptr[5+1];\
	real const tmp9 = tmp7 * tmp1;\
	real const tmp11 = tmp9 * tmp3;\
	real const tmp13 = n3_l.x * (X)->ptr[5+2];\
	real const tmp15 = tmp13 * tmp1;\
	real const tmp17 = tmp15 * tmp3;\
	real const tmp20 = (X)->ptr[5+6] * tmp1;\
	real const tmp22 = tmp20 * tmp3;\
	real const tmp24 = n_l.x * (X)->ptr[5+3];\
	real const tmp26 = tmp24 * tmp1;\
	real const tmp27 = tmp26 * tmp3;\
	real const tmp29 = n2_l.x * (X)->ptr[5+4];\
	real const tmp31 = tmp29 * tmp1;\
	real const tmp33 = tmp31 * tmp3;\
	real const tmp35 = n3_l.x * (X)->ptr[5+5];\
	real const tmp37 = tmp35 * tmp1;\
	real const tmp39 = tmp37 * tmp3;\
	real const tmp42 = (X)->ptr[5+7] * tmp1;\
	real const tmp44 = tmp42 * tmp3;\
	real const tmp45 = n3_l.y * (X)->ptr[5+5];\
	real const tmp47 = sqrt_1_mu_g * tmp45;\
	real const tmp48 = n2_l.y * (X)->ptr[5+4];\
	real const tmp50 = sqrt_1_mu_g * tmp48;\
	real const tmp51 = n_l.y * (X)->ptr[5+3];\
	real const tmp53 = sqrt_1_mu_g * tmp51;\
	real const tmp54 = n3_l.z * (X)->ptr[5+2];\
	real const tmp56 = sqrt_1_eps_g * tmp54;\
	real const tmp57 = n2_l.z * (X)->ptr[5+1];\
	real const tmp59 = sqrt_1_eps_g * tmp57;\
	real const tmp60 = n_l.z * (X)->ptr[5+0];\
	real const tmp62 = sqrt_1_eps_g * tmp60;\
	real const tmp63 = tmp59 * tmp3;\
	real const tmp64 = tmp62 * tmp3;\
	real const tmp65 = tmp56 * tmp3;\
	real const tmp67 = tmp53 * tmp3;\
	real const tmp69 = tmp50 * tmp3;\
	real const tmp71 = tmp47 * tmp3;\
	real const tmp73 = n_l.y * (X)->ptr[5+0];\
	real const tmp75 = sqrt_1_eps_g * tmp73;\
	real const tmp76 = tmp75 * tmp3;\
	real const tmp77 = n2_l.y * (X)->ptr[5+1];\
	real const tmp79 = sqrt_1_eps_g * tmp77;\
	real const tmp81 = tmp79 * tmp3;\
	real const tmp82 = n3_l.y * (X)->ptr[5+2];\
	real const tmp84 = sqrt_1_eps_g * tmp82;\
	real const tmp86 = tmp84 * tmp3;\
	real const tmp87 = n_l.z * (X)->ptr[5+3];\
	real const tmp89 = sqrt_1_mu_g * tmp87;\
	real const tmp90 = n2_l.z * (X)->ptr[5+4];\
	real const tmp92 = sqrt_1_mu_g * tmp90;\
	real const tmp93 = n3_l.z * (X)->ptr[5+5];\
	real const tmp95 = sqrt_1_mu_g * tmp93;\
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
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
	real const sqrt_eps_g = sqrt(eps_g);	/*  TODO sqrt units */\
	real const sqrt_mu_g = sqrt(mu_g);\
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
	real const sqrt_1_eps_g = 1. / sqrt_eps_g;\
	real const sqrt_1_mu_g = 1. / sqrt_mu_g;\
\
	real3 const n_l = normal_l1(n);\
	real3 const n2_l = normal_l2(n);\
	real3 const n3_l = normal_l3(n);\
\
	(Y)->ptr[5+0] = sqrt_1_eps_g * n_l.x * (X)->ptr[5+6] + sqrt_eps_g * n_l.y * (X)->ptr[5+5] + sqrt_eps_g * n_l.z * (X)->ptr[5+2] - sqrt_eps_g * n_l.y * (X)->ptr[5+3] - sqrt_eps_g * n_l.z * (X)->ptr[5+4] + -sqrt_1_eps_g * n_l.x * (X)->ptr[5+0];\
	(Y)->ptr[5+1] = sqrt_1_eps_g * n2_l.x * (X)->ptr[5+6] + sqrt_eps_g * n2_l.y * (X)->ptr[5+5] + sqrt_eps_g * n2_l.z * (X)->ptr[5+2] - sqrt_eps_g * n2_l.y * (X)->ptr[5+3] - sqrt_eps_g * n2_l.z * (X)->ptr[5+4] + -sqrt_1_eps_g * n2_l.x * (X)->ptr[5+0];\
	(Y)->ptr[5+2] = sqrt_1_eps_g * n3_l.x * (X)->ptr[5+6] + sqrt_eps_g * n3_l.y * (X)->ptr[5+5] + sqrt_eps_g * n3_l.z * (X)->ptr[5+2] - sqrt_eps_g * n3_l.y * (X)->ptr[5+3] - sqrt_eps_g * n3_l.z * (X)->ptr[5+4] + -sqrt_1_eps_g * n3_l.x * (X)->ptr[5+0];\
	(Y)->ptr[5+3] = sqrt_mu_g * n_l.z * (X)->ptr[5+5] + sqrt_mu_g * n_l.y * (X)->ptr[5+4] + sqrt_mu_g * n_l.z * (X)->ptr[5+3] + sqrt_mu_g * n_l.y * (X)->ptr[5+2] + sqrt_1_eps_g * n_l.x * (X)->ptr[5+7] + -sqrt_1_eps_g * n_l.x * (X)->ptr[5+1];\
	(Y)->ptr[5+4] = sqrt_mu_g * n2_l.z * (X)->ptr[5+5] + sqrt_mu_g * n2_l.y * (X)->ptr[5+4] + sqrt_mu_g * n2_l.z * (X)->ptr[5+3] + sqrt_mu_g * n2_l.y * (X)->ptr[5+2] + sqrt_1_eps_g * n2_l.x * (X)->ptr[5+7] + -sqrt_1_eps_g * n2_l.x * (X)->ptr[5+1];\
	(Y)->ptr[5+5] = sqrt_mu_g * n3_l.z * (X)->ptr[5+5] + sqrt_mu_g * n3_l.y * (X)->ptr[5+4] + sqrt_mu_g * n3_l.z * (X)->ptr[5+3] + sqrt_mu_g * n3_l.y * (X)->ptr[5+2] + sqrt_1_eps_g * n3_l.x * (X)->ptr[5+7] + -sqrt_1_eps_g * n3_l.x * (X)->ptr[5+1];\
	(Y)->ptr[5+6] = sqrt_1_eps_g * (X)->ptr[5+6] + sqrt_1_eps_g * (X)->ptr[5+0];\
	(Y)->ptr[5+7] = sqrt_1_eps_g * (X)->ptr[5+7] + sqrt_1_eps_g * (X)->ptr[5+1];\
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
	real const nLen = normal_len(n);\
	const real gamma = solver->heatCapacityRatio;\
	const real gamma_1 = gamma - 1.;\
	const real gamma_2 = gamma - 2.;\
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
	(resultFlux)->ptr[5] =\
		0;\
\
	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;\
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;\
	real const eps_g = 1. / (4. * M_PI * G);\
	real const mu_g = 1. / (eps_g * speedOfLightSq);\
\
	real3 const E_g = real3_real_mul((X)->D_g, 1. / eps_g);\
	real3 const H_g = real3_real_mul((X)->B_g, 1. / mu_g);\
	if (n.side == 0) {\
		(resultFlux)->D_g = _real3(solver->divPhiWavespeed_g / unit_m_per_s * (X)->phi_g, H_g.z, -H_g.y);\
		(resultFlux)->B_g = _real3(solver->divPsiWavespeed_g / unit_m_per_s * (X)->psi_g, -E_g.z, E_g.y);\
	} else if (n.side == 1) {\
		(resultFlux)->D_g = _real3(-H_g.z, solver->divPhiWavespeed_g / unit_m_per_s * (X)->phi_g, H_g.x);\
		(resultFlux)->B_g = _real3(E_g.z, solver->divPsiWavespeed_g / unit_m_per_s * (X)->psi_g, -E_g.x);\
	} else if (n.side == 2) {\
		(resultFlux)->D_g = _real3(H_g.y, -H_g.x, solver->divPhiWavespeed_g / unit_m_per_s * (X)->phi_g);\
		(resultFlux)->B_g = _real3(-E_g.y, E_g.x, solver->divPsiWavespeed_g / unit_m_per_s * (X)->psi_g);\
	}\
	(resultFlux)->phi_g = solver->divPhiWavespeed_g / unit_m_per_s * normal_vecDotN1(n, E_g);\
	(resultFlux)->psi_g = solver->divPsiWavespeed_g / unit_m_per_s * normal_vecDotN1(n, H_g);\
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS_NOGHOST?>

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
		real const PL = calc_P(solver, UL, xL);
		real const PR = calc_P(solver, UR, xR);
	
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
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);
	
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
		calcFluxFromCons(&F, *U, x);
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] += commTrace.s<?=i?> * flux.ptr[j];
		}
	}<? end ?>
<? 	end ?>
<? end ?>
<? end -- vectorComponent == "anholonomic" ?>


	real const G = solver->sqrt_G * solver->sqrt_G / unit_m3_per_kg_s2;
	real const speedOfLightSq = solver->speedOfLight * solver->speedOfLight / unit_m2_per_s2;
	real const eps_g = 1. / (4. * M_PI * G);
	real const mu_g = 1. / (eps_g * speedOfLightSq);

//// MODULE_DEPENDS: <?=eqn_common?>
	real3 const gravForce = calcGravForce(solver, U, x);
	
	deriv->m = real3_add(deriv->m, gravForce);

	/* source of D_g is J_g is the momentum + Poynting vector */
	real3 J_g = real3_zero;
	/*  I'm symmetrizing the stress-energy */
		/* matter */
	J_g.x -= U->m.x;
	J_g.y -= U->m.y;
	J_g.z -= U->m.z;
	
	deriv->D_g.x -= J_g.x;
	deriv->D_g.y -= J_g.y;
	deriv->D_g.z -= J_g.z;

	/*  source of phi_g is T_00 is rho + .5 (D^2 + B^2) */
	real const T_00_over_c2 = 0
		/* matter */
		+ U->rho
	;
	
	deriv->phi_g += T_00_over_c2 * solver->divPhiWavespeed_g / unit_m_per_s;
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
