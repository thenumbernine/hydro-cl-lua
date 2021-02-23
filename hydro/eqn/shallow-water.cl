//// MODULE_NAME: <?=primFromCons?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define <?=primFromCons?>(\
	/*<?=prim_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U, \
	/*real3 const */pt\
) {\
	(result)->h = (U)->h;\
	(result)->v = real3_real_mul((U)->m, 1. / (U)->h);\
	(result)->depth = (U)->depth;\
}

//// MODULE_NAME: <?=consFromPrim?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define <?=consFromPrim?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */W,\
	/*real3 const */pt\
) {\
	(result)->h = (W)->h;\
	(result)->m = real3_real_mul((W)->v, (W)->h);\
	(result)->depth = (W)->depth;\
}

//// MODULE_NAME: <?=apply_dU_dW?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define <?=apply_dU_dW?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA, \
	/*<?=prim_t?> const * const */W, \
	/*real3 const */pt\
) {\
	(result)->h = (W)->h;\
	(result)->m = real3_add(\
		real3_real_mul((WA)->v, (W)->h), \
		real3_real_mul((W)->v, (WA)->h));\
	(result)->depth = (W)->depth;\
}

//// MODULE_NAME: <?=apply_dW_dU?>
//// MODULE_DEPENDS: real3 <?=solver_t?> <?=prim_t?> <?=cons_t?>

#define <?=apply_dW_dU?>(\
	/*<?=prim_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=prim_t?> const * const */WA,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt\
) {\
	(result)->h = (U)->h;\
	(result)->v = real3_sub(\
		real3_real_mul((U)->m, 1. / (WA)->h),\
		real3_real_mul((WA)->v, (U)->h / (WA)->h));\
	(result)->depth = (U)->depth;\
}

//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?>

#define calc_C(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
)\
	sqrt(solver->gravity * (U)->h)

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: <?=cartesianToCoord?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	
	// this is in all init/euler.lua
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
	real ePot = 0;

	<?=initCode()?>

	<?=prim_t?> const W = {
		.h = rho,
		.v = cartesianToCoord(v, x),
		.depth = 0.,
	};
	<?=consFromPrim?>(U, solver, &W, x);
}

//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=primFromCons?> <?=normal_t?>

#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
	real const v_n = normal_vecDotN1(n, W.v);\
	real3 const nU = normal_u1(n);\
	(result)->h = (U)->h * v_n,\
	(result)->m = real3_add(\
		real3_real_mul((U)->m, v_n),	/*h v^i v_n*/\
		real3_real_mul(nU, .5 * solver->gravity * (U)->h * (U)->h)	/*.5 g h^2 n^i*/\
	);\
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>

// used by PLM
#define <?=calcCellMinMaxEigenvalues?>(\
	/*<?=range_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
	real const v_n = normal_vecDotN1(n, W.v);\
	real const C = calc_C(solver, U);\
	real const C_nLen = C * normal_len(n);\
	(result)->min = v_n - C_nLen;\
	(result)->max = v_n + C_nLen;\
}

//// MODULE_NAME: <?=eigen_forInterface?>

#define <?=eigen_forInterface?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */UL,\
	/*<?=cons_t?> const * const */UR,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> WL;\
	<?=primFromCons?>(&WL, solver, UL, pt);\
	real const sqrtHL = sqrt(WL.h);\
	real3 const vL = WL.v;\
	\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, pt);\
	real const sqrtHR = sqrt(WR.h);\
	real3 const vR = WR.v;\
	\
	real const invDenom = 1. / (sqrtHL + sqrtHR);\
	\
	/*Roe-averaged:*/\
	(result)->h = sqrtHL * sqrtHR;\
	(result)->v = real3_add(\
		real3_real_mul(vL, sqrtHL * invDenom),\
		real3_real_mul(vR, sqrtHR * invDenom));\
\
	/*derived:*/\
	real const CSq = solver->gravity * (result)->h;\
	(result)->C = sqrt(CSq);\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=waves_t?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X_,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) { \
	real3 const nL = normal_l1(n);\
	real const nLen = normal_len(n);\
	real const nLenSq = nLen * nLen;\
	real3 const v_ns = normal_vecDotNs(n, (eig)->v);\
	real const v_n = v_ns.x;\
	real3 const n2U = normal_u2(n);\
	real3 const n3U = normal_u3(n);\
\
	/*\
	TODO double-check that this is the correct 3D generalization...\
	 seems strange to have an upper dot an upper\
	*/\
	real const vU_dot_n2U = real3_dot((eig)->v, n2U);\
	real const vU_dot_n3U = real3_dot((eig)->v, n3U);\
\
	real const C_nLen = (eig)->C * nLen;\
	real const denom = 2. * (eig)->h * nLenSq;\
	real const invDenom = 1. / denom;\
\
	(result)->ptr[0] = (\
			  (X_)->ptr[0] * (-C_nLen - v_n) \
			+ (X_)->ptr[1] * nL.x \
			+ (X_)->ptr[2] * nL.y \
			+ (X_)->ptr[3] * nL.z\
		) * invDenom;\
	(result)->ptr[1] = (\
			  (X_)->ptr[0] * -vU_dot_n2U \
			+ (X_)->ptr[1] * n2U.x \
			+ (X_)->ptr[2] * n2U.y \
			+ (X_)->ptr[3] * n2U.z\
		) * 2. * invDenom;\
	(result)->ptr[2] = (\
			  (X_)->ptr[0] * -vU_dot_n3U \
			+ (X_)->ptr[1] * n3U.x \
			+ (X_)->ptr[2] * n3U.y \
			+ (X_)->ptr[3] * n3U.z\
		) * 2. * invDenom;\
	(result)->ptr[3] = (\
			  (X_)->ptr[0] * (C_nLen - v_n) \
			+ (X_)->ptr[1] * nL.x \
			+ (X_)->ptr[2] * nL.y \
			+ (X_)->ptr[3] * nL.z\
		) * invDenom;\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=waves_t?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X_,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real3 const nU = normal_u1(n);\
	real const nLen = normal_len(n);\
	real const invC = 1. / (eig)->C;\
\
	(result)->ptr[0] = \
		nLen * (eig)->h * invC * ((X_)->ptr[3] - (X_)->ptr[0]);\
		\
	(result)->ptr[1] = \
		((X_)->ptr[0] + (X_)->ptr[3]) * (eig)->h * nU.x \
		+ ((X_)->ptr[3] - (X_)->ptr[0]) * (nLen * (eig)->h * invC * (eig)->v.x)\
		+ (eig)->h * ((X_)->ptr[1] * normal_l2x(n) + (X_)->ptr[2] * normal_l3x(n));\
		\
	(result)->ptr[2] = \
		((X_)->ptr[0] + (X_)->ptr[3]) * (eig)->h * nU.y\
		+ ((X_)->ptr[3] - (X_)->ptr[0]) * (nLen * (eig)->h * invC * (eig)->v.y)\
		+ (eig)->h * ((X_)->ptr[1] * normal_l2y(n) + (X_)->ptr[2] * normal_l3y(n));\
\
	(result)->ptr[3] = \
		((X_)->ptr[0] + (X_)->ptr[3]) * (eig)->h * nU.z\
		+ ((X_)->ptr[3] - (X_)->ptr[0]) * (nLen * (eig)->h * invC * (eig)->v.z)\
		+ (eig)->h * ((X_)->ptr[1] * normal_l2z(n) + (X_)->ptr[2] * normal_l3z(n));\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X_,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real3 const nL = normal_l1(n);\
	real3 const nU = normal_u1(n);\
	real const CSq = (eig)->C * (eig)->C;\
	real3 const v_ns = normal_vecDotNs(n, (eig)->v);\
	real const v_n = v_ns.x;\
\
	(result)->ptr[0] = \
		(X_)->ptr[1] * nL.x \
		+ (X_)->ptr[2] * nL.y \
		+ (X_)->ptr[3] * nL.z;\
		\
	(result)->ptr[1] = \
		(X_)->ptr[0] * (CSq * nU.x - (eig)->v.x * v_n)\
		+ (X_)->ptr[1] * ((eig)->v.x * nL.x + v_n)\
		+ (X_)->ptr[2] * ((eig)->v.x * nL.y)\
		+ (X_)->ptr[3] * ((eig)->v.x * nL.z);\
		\
	(result)->ptr[2] = \
		(X_)->ptr[0] * (CSq * nU.y - (eig)->v.y * v_n)\
		+ (X_)->ptr[1] * ((eig)->v.y * nL.x)\
		+ (X_)->ptr[2] * ((eig)->v.y * nL.y + v_n)\
		+ (X_)->ptr[3] * ((eig)->v.y * nL.z);\
		\
	(result)->ptr[3] = \
		(X_)->ptr[0] * (CSq * nU.z - (eig)->v.z * v_n)\
		+ (X_)->ptr[1] * ((eig)->v.z * nL.x)\
		+ (X_)->ptr[2] * ((eig)->v.z * nL.y)\
		+ (X_)->ptr[3] * ((eig)->v.z * nL.z + v_n);\
}

//// MODULE_NAME: <?=eigen_forCell?>

// used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */pt,\
	/*real3 const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, pt);\
	real const CSq = solver->gravity * (U)->h;\
	real const C = sqrt(CSq);\
	(result)->h = W.h;\
	(result)->v = W.v;\
	(result)->C = C;\
}

//// MODULE_NAME: <?=addSource?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	real3 const pt = cellBuf[index].pos;
	
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);

//// MODULE_DEPENDS: <?=cell_volume?>
	//I'm treating conn as if it is fixed at the cell center.
	//TODO instead: treat the cell values constant at center and integrate conn across the cell.
	real const volume = cell_volume(pt);

//// MODULE_DEPENDS: <?=coord_conn_apply23?>
	//covariant derivative connection:
	//integral -vol Γ^i_jk h v^j v^k dx^3
	deriv->m = real3_sub(
		deriv->m,
		real3_real_mul(
			coord_conn_apply23(W.v, U->m, pt),
			volume
		)
	);
	
//// MODULE_DEPENDS: <?=coord_conn_trace23?>
	//covariant derivative connection:
	//integral -vol 1/2 g h^2 Γ^i_jk g^jk dx^3
	deriv->m = real3_sub(
		deriv->m,
		real3_real_mul(
			coord_conn_trace23(pt),
			.5 * volume * solver->gravity * U->h * U->h
		)
	);

// \partial_tilde{j} depth
<?=eqn:makePartial1"depth"?>	

//// MODULE_DEPENDS: <?=coord_holBasisLen_i?>
	//e_j(depth) = {e_j}^\tilde{j} \partial_\tilde{j} (depth)
	real3 const e_depth_l = _real3(
		partial_depth_l.x / coord_holBasisLen0(pt),
		partial_depth_l.y / coord_holBasisLen1(pt),	
		partial_depth_l.z / coord_holBasisLen2(pt)
	);

	//2011 Berger, eqn 1
//// MODULE_DEPENDS: <?=coord_raise?>
	//g^ij e_j(depth)
	real3 const e_depth_u = coord_raise(e_depth_l, pt);	

	//source of shallow water equations:
	//integral vol g h g^ij e_j(d) dx^3
	deriv->m = real3_sub(
		deriv->m,
		real3_real_mul(
			e_depth_u,
			volume * solver->gravity * U->h
		)
	);
	
	//2011 Berger, eqn 3
	real const drag = solver->gravity * solver->Manning * solver->Manning * coordLen(U->m, pt) * pow(U->h, -8./3.);	//|u|/h^(5/3) = |m|*h^(-8/3)
	deriv->m = real3_sub(deriv->m, U->m);
}
