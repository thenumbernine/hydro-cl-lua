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
<? if not eqn.depthInCell then ?>\
	(result)->depth = (U)->depth;\
<? end ?>\
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
<? if not eqn.depthInCell then ?>\
	(result)->depth = (W)->depth;\
<? end ?>\
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
<? if not eqn.depthInCell then ?>\
	(result)->depth = (W)->depth;\
<? end ?>\
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
<? if not eqn.depthInCell then ?>\
	(result)->depth = (U)->depth;\
<? end ?>
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
	global <?=cell_t?> * const cell
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
	//specific to shallow water problems:
	real depth = 0;	//same as 'water_H'

	<?=initCode()?>

	<?=prim_t?> W = {
		.h = rho,
		.v = cartesianToCoord(v, x),
	};

	<?=getDepth("&W", "cell")?> = depth;
	
	<?=consFromPrim?>(U, solver, &W, x);
}

<? if not eqn.useDefaultFluxFromCons then ?>
//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=primFromCons?> <?=normal_t?>

/* TODO why is there an error when I include the 1/2 in the F->m term? */
#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, (cell)->pos);\
	real const v_n = normal_vecDotN1(n, W.v);\
	real3 const nU = normal_u1(n);\
	(result)->h = (U)->h * v_n,\
	(result)->m = real3_add(\
		real3_real_mul((U)->m, v_n),	/*h v^i v_n*/\
		real3_real_mul(nU, .5 * solver->gravity * (U)->h * ((U)->h\
<? if eqn.extraTermInFlux then ?>\
			- 2. * <?=getDepth("U", "cell")?>\
<? end ?>\
		))	/*.5 g h^2 n^i*/\
	);\
}
<? end -- not eqn.useDefaultFluxFromCons ?>

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
	/*<?=eigen_t?> * const */resultEig,\
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
	real const sqrtHL = sqrt(WL.h);\
	real3 const vL = WL.v;\
\
	<?=prim_t?> WR;\
	<?=primFromCons?>(&WR, solver, UR, (cellR)->pos);\
	real const sqrtHR = sqrt(WR.h);\
	real3 const vR = WR.v;\
\
	real const invDenom = 1. / (sqrtHL + sqrtHR);\
\
	/*Roe-averaged:*/\
\
	/*(resultEig)->h = sqrtHL * sqrtHR;*/\
	(resultEig)->h = .5 * ((UL)->h + (UR)->h);\
\
	(resultEig)->v = real3_add(\
		real3_real_mul(vL, sqrtHL * invDenom),\
		real3_real_mul(vR, sqrtHR * invDenom));\
<? if not eqn.depthInCell then ?>\
	(resultEig)->depth = .5 * (<?=getDepth("UR","cellR")?> + <?=getDepth("UL","cellL")?>);\
<? end ?>\
\
	/*derived:*/\
\
	real const CSq = solver->gravity * (resultEig)->h;\
	(resultEig)->C = sqrt(CSq);\
}

//// MODULE_NAME: <?=eigen_forCell?>

// used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */resultEig,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*real3 const */n\
) {\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, (cell)->pos);\
	real const CSq = solver->gravity * (U)->h;\
	real const C = sqrt(CSq);\
	(resultEig)->h = W.h;\
	(resultEig)->v = W.v;\
	(resultEig)->C = C;\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=waves_t?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) { \
	real3 const nL = normal_l1(n);\
	real3 const n2L = normal_l2(n);\
	real3 const n3L = normal_l3(n);\
	real const nLen = normal_len(n);\
	real const v_n = normal_vecDotN1(n, (eig)->v);\
\
	real const invDenom = 1. / (2. * (eig)->h * nLen * nLen);\
\
	real const tmp1 = nLen * (X)->ptr[0];\
	real const tmp2 = (eig)->C * tmp1;\
	real const tmp3 = nL.x * (X)->ptr[1];\
	real const tmp5 = nL.y * (X)->ptr[2];\
	real const tmp7 = nL.z * (X)->ptr[3];\
	real const tmp8 = v_n * (X)->ptr[0];\
	real const tmp12 = (eig)->v.x * (X)->ptr[0];\
	real const tmp13 = (eig)->v.y * (X)->ptr[0];\
	real const tmp14 = (eig)->v.z * (X)->ptr[0];\
	(result)->ptr[0] = invDenom * -(tmp8 + tmp2 - tmp3 - tmp5 - tmp7);\
	(result)->ptr[1] = invDenom * -2. * (n2L.z * tmp14 - n2L.x * (X)->ptr[1] - n2L.y * (X)->ptr[2] - n2L.z * (X)->ptr[3] + n2L.y * tmp13 + n2L.x * tmp12);\
	(result)->ptr[2] = invDenom * -2. * (n3L.z * tmp14 - n3L.x * (X)->ptr[1] - n3L.y * (X)->ptr[2] - n3L.z * (X)->ptr[3] + n3L.y * tmp13 + n3L.x * tmp12);\
	(result)->ptr[3] = invDenom * -(tmp8 - tmp2 - tmp3 - tmp5 - tmp7);\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=waves_t?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real3 const nU = normal_u1(n);\
	real3 const n2U = normal_u2(n);\
	real3 const n3U = normal_u3(n);\
	real const nLen = normal_len(n);\
\
	real const invDenom = 1. / (eig)->C;\
\
	real const tmp5 = (eig)->C * (X)->ptr[0];\
	real const tmp6 = nLen * (X)->ptr[0];\
	real const tmp9 = (eig)->C * (X)->ptr[1];\
	real const tmp10 = (eig)->C * (X)->ptr[2];\
	real const tmp11 = (eig)->C * (X)->ptr[3];\
	real const tmp12 = nLen * (X)->ptr[3];\
	real const tmp22 = invDenom * (eig)->h;\
	(result)->ptr[0] = invDenom * -nLen * (eig)->h * ((X)->ptr[0] - (X)->ptr[3]);\
	(result)->ptr[1] = tmp22 * ((eig)->v.x * tmp12 + nU.x * tmp11 + n3U.x * tmp10 + n2U.x * tmp9 + nU.x * tmp5 - (eig)->v.x * tmp6);\
	(result)->ptr[2] = tmp22 * ((eig)->v.y * tmp12 + nU.y * tmp11 + n3U.y * tmp10 + n2U.y * tmp9 + nU.y * tmp5 - (eig)->v.y * tmp6);\
	(result)->ptr[3] = tmp22 * ((eig)->v.z * tmp12 + nU.z * tmp11 + n3U.z * tmp10 + n2U.z * tmp9 + nU.z * tmp5 - (eig)->v.z * tmp6);\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>

#define <?=eigen_fluxTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real3 const nL = normal_l1(n);\
	real3 const nU = normal_u1(n);\
	real const v_n = normal_vecDotN1(n, (eig)->v);\
\
	real const tmp5 = v_n * (X)->ptr[0];\
	real const tmp6 = (eig)->h * (X)->ptr[0];\
	real const tmp7 = solver->gravity * tmp6;\
	(result)->ptr[0] = nL.z * (X)->ptr[3] + nL.y * (X)->ptr[2] + nL.x * (X)->ptr[1];\
	(result)->ptr[1] = -((eig)->v.x * tmp5 - nU.x * tmp7 - nL.x * (eig)->v.x * (X)->ptr[1] - v_n * (X)->ptr[1] - nL.y * (eig)->v.x * (X)->ptr[2] - nL.z * (eig)->v.x * (X)->ptr[3]);\
	(result)->ptr[2] = -((eig)->v.y * tmp5 - nU.y * tmp7 - nL.x * (eig)->v.y * (X)->ptr[1] - nL.y * (eig)->v.y * (X)->ptr[2] - v_n * (X)->ptr[2] - nL.z * (eig)->v.y * (X)->ptr[3]);\
	(result)->ptr[3] = -((eig)->v.z * tmp5 - nU.z * tmp7 - nL.x * (eig)->v.z * (X)->ptr[1] - nL.y * (eig)->v.z * (X)->ptr[2] - nL.z * (eig)->v.z * (X)->ptr[3] - v_n * (X)->ptr[3]);\
}

//// MODULE_NAME: <?=addSource?>

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
	
	real3 const pt = cell->pos;
	
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, x);

//// MODULE_DEPENDS: <?=cell_volume?>
	//I'm treating conn as if it is fixed at the cell center.
	//TODO instead: treat the cell values constant at center and integrate conn across the cell.
	real const volume = cell_volume(solver, pt);

#if 0
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
#endif
#if 1
// \partial_tilde{j} depth
<?=eqn:makePartial1(
	"depth",			-- field
	"real",				-- fieldType ... needs to be manually specified if we are pulling from a field that is not in consStruct, where the inferred types are checked
	nil,				-- nameOverride
	getDepthSource()	-- srcName (cellBuf instead of UBuf)
)?>

//// MODULE_DEPENDS: <?=coord_holBasisLen_i?>
	//e_j(depth) = {e_j}^\tilde{j} \partial_\tilde{j} (depth)
	real3 const e_depth_l = _real3(
		partial_depth_l.x / coord_holBasisLen0(pt),
		partial_depth_l.y / coord_holBasisLen1(pt),	
		partial_depth_l.z / coord_holBasisLen2(pt)
	);

	//negate changes in h due to differences in water from resting depth?
	//h += v^j nabla_j H
	//deriv->h += real3_dot(U->m, e_depth_l) / U->h;
	//but then again, this should only happen proportional to the velocity, when at zero shouldn't change the height

	//2011 Berger, eqn 1
//// MODULE_DEPENDS: <?=coord_raise?>
	//g^ij e_j(depth)
	real3 const e_depth_u = coord_raise(e_depth_l, pt);	

	//source of shallow water equations:
	//integral vol g h g^ij e_j(d) dx^3
	// this value needs to match whatever is generated by the flux update of d/dt h when H varies in order for constant-(h-H) to be a steady-state.
	deriv->m = real3_add(
		deriv->m,
		real3_real_mul(
			e_depth_u,
			solver->gravity * U->h
		)
	);
#endif
#if 0
	//2011 Berger, eqn 3
	real const drag = solver->gravity * solver->Manning * solver->Manning * coordLen(U->m, pt) * pow(U->h, -8./3.);	//|u|/h^(5/3) = |m|*h^(-8/3)
	deriv->m = real3_sub(deriv->m, U->m);
#endif
}
