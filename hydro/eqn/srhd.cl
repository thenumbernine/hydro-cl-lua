//// MODULE_NAME: eqn.common
//// MODULE_DEPENDS: coordLenSq cons_only_t,prim_only_t

//pressure function for ideal gas
#define calc_P(\
	/*constant solver_t const * const */solver,\
	/*real const */rho,\
	/*real const */eInt\
)\
	((solver->heatCapacityRatio - 1.) * (rho) * (eInt))

//chi in most papers
#define calc_dP_drho(\
	/*constant solver_t const * const */solver,\
	/*real const */rho,\
	/*real const */eInt\
)\
	((solver->heatCapacityRatio - 1.) * (eInt))

//kappa in most papers
#define calc_dP_deInt(\
	/*constant solver_t const * const */solver,\
	/*real const */rho,\
	/*real const */eInt\
)\
	((solver->heatCapacityRatio - 1.) * (rho))

#define calc_eInt_from_P(\
	/*constant solver_t const * const */solver,\
	/*real const */rho,\
	/*real const */P\
)\
	((P) / ((solver->heatCapacityRatio - 1.) * (rho)))

#define calc_h(\
	/*real const */rho,\
	/*real const */P,\
	/*real const */eInt\
)\
	(1. + (eInt) + (P) / (rho))

#define consFromPrimOnly(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*<?=eqn.prim_only_t?> const * const */prim,\
	/*real3 const */x\
) {\
	real const vSq = coordLenSq((prim)->v, x);\
	real const WSq = 1. / (1. - vSq);\
	real const W = sqrt(WSq);\
	real const P = calc_P(solver, (prim)->rho, (prim)->eInt);\
	real const h = calc_h((prim)->rho, P, (prim)->eInt);\
\
	/* 2008 Font, eqn 40-42: */\
\
	/* rest-mass density = J^0 = rho u^0 */\
	(result)->D = (prim)->rho * W;	\
\
	/* momentum = T^0i = rho h u^0 u^i + P g^0i */\
	(result)->S = real3_real_mul((prim)->v, (prim)->rho * h * WSq);\
\
	/* energy = T^00 = rho h u^0 u^0 + P g^00 */\
	(result)->tau = (prim)->rho * h * WSq - (result)->D - P;\
\
	(result)->rho = (prim)->rho;\
	(result)->v = (prim)->v;\
	(result)->eInt = (prim)->eInt;\
	(result)->ePot = 0;\
}

//build the cons_only_t from the cons_t's prim_only_t fields
//used for checking the error between cons_only_t and its prim-reconstructed-from-cons_only_t
#define consOnlyFromPrim(\
	/*<?=eqn.cons_only_t?> * const */result,\
	/*constant solver_t const * const */solver, \
	/*cons_t const * const */U,\
	/*real3 const */x\
) {\
	real const vSq = coordLenSq((U)->v, x);\
	real const WSq = 1. / (1. - vSq);\
	real const W = sqrt(WSq);\
	real const P = calc_P(solver, (U)->rho, (U)->eInt);\
	real const h = calc_h((U)->rho, P, (U)->eInt);\
\
	/* 2008 Font, eqn 40-42: */\
\
	/* rest-mass density = J^0 = rho u^0 */\
	(result)->D = (U)->rho * W;	\
\
	/* momentum = T^0i = rho h u^0 u^i + P g^0i */\
	(result)->S = real3_real_mul((U)->v, (U)->rho * h * WSq);\
\
	/* energy = T^00 = rho h u^0 u^0 + P g^00 */\
	(result)->tau = (U)->rho * h * WSq - (result)->D - P;\
}

#define primOnlyFromCons(\
	/*<?=eqn.prim_only_t?> * const */result,\
	/*constant solver_t const * const */solver, \
	/*cons_t const * const */U,\
	/*real3 const */x\
) {\
	(result)->rho = (U)->rho;\
	(result)->v = (U)->v;\
	(result)->eInt = (U)->eInt;\
}

//PLM uses prim_only_t and cons_t, esp using the 'numIntStates' reals that they start with
//...and PLM uses consFromPrim and primFromCons

//// MODULE_NAME: applyInitCond
//// MODULE_DEPENDS: eqn.common

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
	bool const lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	/* ignored: */
	real3 D = real3_zero;
	real3 B = real3_zero;

	<?=initCode()?>
	
	<?=eqn.prim_only_t?> prim = {
		.rho = rho,
		.v = v,
		.eInt = calc_eInt_from_P(solver, rho, P),
	};
	consFromPrimOnly(UBuf + index, solver, &prim, x);
}


//// MODULE_NAME: fluxFromCons
//// MODULE_DEPENDS: solver_t cons_t normal_t eqn.common

#define fluxFromCons(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
	real const v_n = normal_vecDotN1(n, (U)->v);\
	real const P = calc_P(solver, (U)->rho, (U)->eInt);\
	(result)->D = (U)->D * v_n;\
	(result)->S = real3_add(\
		real3_real_mul((U)->S, v_n),\
		_real3(\
			normal_u1x(n) * P,\
			normal_u1y(n) * P,\
			normal_u1z(n) * P\
		)\
	);\
	(result)->tau = (U)->tau * v_n + P * v_n;\
	(result)->rho = 0;\
	(result)->v = real3_zero;\
	(result)->eInt = 0;\
	(result)->ePot = 0;\
}

//// MODULE_NAME: calcDT
//// MODULE_DEPENDS: SETBOUNDS coordLenSq eqn.common normal_t

//everything matches the default except the params passed through to calcCellMinMaxEigenvalues
kernel void calcDT(
	constant solver_t const * const solver,
	global real * const dtBuf,
	global cons_t const * const UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 const x = cellBuf[index].pos;

	global cons_t const * const U = UBuf + index;
	real const rho = (U)->rho;
	real const eInt = (U)->eInt;
	real const vSq = coordLenSq((U)->v, x);
	real const P = calc_P(solver, rho, eInt);
	real const h = calc_h(rho, P, eInt);
	real const csSq = solver->heatCapacityRatio * P / (rho * h);
	real const cs = sqrt(csSq);
	
	real dt = INFINITY;
	/* for (int side = 0; side < dim; ++side) { */
	<? for side=0,solver.dim-1 do ?>{
		normal_t const n = normal_forSide<?=side?>(x);
		/* for the particular direction */
		real const vi = normal_vecDotN1(n, (U)->v);
		real const viSq = vi * vi;
		
		/*  Marti 1998 eqn 19 */
		/*  also Marti & Muller 2008 eqn 68 */
		/*  also Font 2008 eqn 106 */
		real const discr = sqrt((1. - vSq) * (1. - vSq * csSq - viSq * (1. - csSq)));
		real const lambdaMin = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq);
		real const lambdaMax = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq);
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
	
	dtBuf[index] = dt; 
}

//// MODULE_NAME: eigen_forInterface
//// MODULE_DEPENDS: coord_lower

#define eigen_forInterface(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */UL,\
	/*cons_t const * const */UR,\
	/*real3 const */xInt,\
	/*normal_t const */n\
) {\
<? if true then -- arithmetic averaging ?>\
	prim_only_t avg = {\
		.rho = .5 * ((UL)->rho + (UR)->rho),\
		.v = real3_real_mul(real3_add((UL)->v, (UR)->v), .5),\
		.eInt = .5 * ((UL)->eInt + (UR)->eInt),\
	};\
<? -- else -- Roe-averaging, Font 2008 eqn 38 ?>\
	/*  weight by k = sqrt(sqrt(g) rho h) */\
<? end ?>\
\
	/*  rotate avg.v into n */\
	avg.v = normal_vecDotNs(n, avg.v);\
\
	real const rho = avg.rho;\
	real3 const v = avg.v;\
	real const eInt = avg.eInt;\
\
	/* TODO NOTE if you're swapping vector components, you have to swap metric components too  */\
	real3 const vL = coord_lower(v, xInt);\
	real const vSq = real3_dot(v, vL);\
	real const oneOverW2 = 1. - vSq;\
	real const oneOverW = sqrt(oneOverW2);\
	real const W = 1. / oneOverW;\
	real const W2 = 1. / oneOverW2;\
	real const P = (solver->heatCapacityRatio - 1.) * rho * eInt;\
	real const h = 1. + eInt + P / rho;\
\
	real const hW = h * W;\
\
	/* just after 2008 Font eqn 107: */\
	/* h cs^2 = chi + P / rho^2 kappa = dp/drho + p / rho^2 dp/deInt */\
	/*  = (gamma-1) eInt + P/rho^2 (gamma-1) rho  for an ideal gas */\
	/*  = (gamma-1) (eInt + P/rho) */\
	/*  = 1/rho ( (gamma-1) rho eInt + (gamma-1) P ) */\
	/*  = 1/rho ( P + (gamma-1) P) */\
	/*  = gamma P / rho */\
	real const vxSq = v.x * v.x;\
	real const csSq = solver->heatCapacityRatio * P / (rho * h);\
	real const cs = sqrt(csSq);\
\
	real const discr = sqrt((1. - vSq) * ((1. - vSq * csSq) - vxSq * (1. - csSq)));\
	real const lambdaMin = (v.x * (1. - csSq) - cs * discr) / (1. - vSq * csSq);\
	real const lambdaMax = (v.x * (1. - csSq) + cs * discr) / (1. - vSq * csSq);\
\
	/* used by evL and evR */\
	real const ATildeMinus = (1. - vxSq) / (1. - v.x * lambdaMin);	/* 2008 Font eqn 113 */\
	real const ATildePlus  = (1. - vxSq) / (1. - v.x * lambdaMax);	/* 2008 Font eqn 113 */\
\
	/* used by evL */\
	real const VMinus = (v.x - lambdaMin) / (1. - v.x * lambdaMin);	/* 2008 Font eqn 113 */\
	real const VPlus = (v.x - lambdaMax) / (1. - v.x * lambdaMax);	/* 2008 Font eqn 113 */\
\
	/* used by evL and evR */\
	real const CMinus = vL.x - VMinus;	/* 2008 Font eqn 112 */\
	real const CPlus = vL.x - VPlus;	/* 2008 Font eqn 112 */\
\
	real const kappa = calc_dP_deInt(solver, rho, eInt);	/* 2008 Font note just after eqn 107 */\
	real const kappaTilde = kappa / rho;	/* 2008 Font eqn 112.   */\
	/* used by evL and evR */\
	real const Kappa = kappaTilde / (kappaTilde - csSq);	/* 2008 Font eqn 112.   */\
	/* Kappa = h;	/ * approx for ideal gas */\
\
<? for _,var in ipairs(eqn.eigenVars) do --\
?>	(eig)-><?=var.name?> = <?=var.name?>;\
<? end --\
?>\
}

//// MODULE_NAME: eigen_forCell

//used by PLM
#define eigen_forCell(\
	/*eigen_t * const */eig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */x,\
	/*normal_t const */n\
)\
	(eigen_forInterface(eig, solver, U, U, x, n))

//// MODULE_NAME: eigen_left/rightTransform 
//// MODULE_DEPENDS: waves_t

<? -- create code to initialize local vars of all the eig vars
local eigVarCode = require "ext.table".map(eqn.eigenVars, function(var)
	return "\t"..var.type.." "..var.name.." = (eig)->"..var.name..";\n"
end):concat()
?>

#define eigen_leftTransform(\
	/*waves_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*cons_t const * const */X_,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	/* rotate incoming v's in X */\
	cons_t X = *X_;\
	X.S = normal_vecDotNs(n, X.S);\
\
	<?=eigVarCode:gsub("\n", "\\\n")?>\
\
	real3 const vL = coord_lower(v, pt);\
	real const vxSq = v.x * v.x;\
	real const hSq = h * h;\
	real const hW = h * W;\
	real const W2 = W * W;\
\
	real const xi = 1. - vxSq;	/* 2008 Font eqn 121 */\
	real const Delta = hSq * hW * (Kappa - 1.) * (CPlus - CMinus) * xi;	/* 2008 Font eqn 121 */\
\
	{\
		/* min row	2008 Font eqn 118 */\
		real const scale = hSq / Delta;\
		real const l5minus = (1 - Kappa) * (-v.x + VPlus * (W2 * xi - 1.)) - Kappa * W2 * VPlus * xi;\
		(result)->ptr[0] = (\
			X.ptr[0] * (hW * VPlus * xi + l5minus)\
			+ X.ptr[1] * (1 - Kappa * ATildePlus + (2. * Kappa - 1.) * VPlus * (W2 * v.x * xi - v.x))\
			+ X.ptr[2] * ((2. * Kappa - 1.) * VPlus * W2 * v.y * xi)\
			+ X.ptr[3] * ((2. * Kappa - 1.) * VPlus * W2 * v.z * xi)\
			+ X.ptr[4] * l5minus\
		) * scale;\
	}\
	{\
		/* mid normal row	2008 Font eqn 115 */\
		real const scale = W / (Kappa - 1.);\
		(result)->ptr[1] = (\
			X.ptr[0] * (h - W) \
			+ X.ptr[1] * (W * v.x) \
			+ X.ptr[2] * (W * v.y) \
			+ X.ptr[3] * (W * v.z) \
			+ X.ptr[4] * (-W)\
		) * scale;\
	}\
	{\
		/* mid tangent A row	2008 Font eqn 116 */\
		real const scale = 1. / (h * xi);\
		(result)->ptr[2] = (\
			X.ptr[0] * -vL.y\
			+ X.ptr[1] * v.x * vL.y\
			+ X.ptr[2] * ((1. - v.x * vL.x))\
			+ X.ptr[4] * -vL.y\
		) * scale;\
		/* mid tangent B row	2008 Font eqn 117 */\
		(result)->ptr[3] = (\
			X.ptr[0] * -vL.z\
			+ X.ptr[1] * v.x * vL.z\
			+ X.ptr[3] * (1. - vL.x * v.x)\
			+ X.ptr[4] * -vL.z\
		) * scale;\
	}\
	{\
		/* max row	2008 Font eqn 118 */\
		real const scale = -hSq / Delta;\
		real const l5plus = (1 - Kappa) * (-v.x + VMinus * (W2 * xi - 1.)) - Kappa * W2 * VMinus * xi;\
		(result)->ptr[4] = (\
			X.ptr[0] * (h * W * VMinus * xi + l5plus)\
			+ X.ptr[1] * (1 - Kappa * ATildeMinus + (2. * Kappa - 1.) * VMinus * (W2 * v.x * xi - v.x))\
			+ X.ptr[2] * ((2. * Kappa - 1.) * VMinus * W2 * v.y * xi)\
			+ X.ptr[3] * ((2. * Kappa - 1.) * VMinus * W2 * v.z * xi)\
			+ X.ptr[4] * l5plus\
		) * scale;\
	}\
}

#define eigen_rightTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*waves_t const * const */X,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	<?=eigVarCode:gsub("\n", "\\\n")?>\
\
	real3 const vL = coord_lower(v, pt);\
	real const hW = h * W;\
	real const W2 = W * W;\
\
	cons_t Y;\
	/* 2008 Font eqns 108-111 */\
	Y.ptr[0] = (X)->ptr[0]\
		+ (X)->ptr[1] * (Kappa / hW)\
		+ (X)->ptr[2] * (W * vL.y)\
		+ (X)->ptr[3] * (W * vL.z)\
		+ (X)->ptr[4];\
	Y.ptr[1] = (X)->ptr[0] * (hW * CMinus)\
		+ (X)->ptr[1] * (vL.x)\
		+ (X)->ptr[2] * (h * 2. * W2 * vL.y * vL.x)\
		+ (X)->ptr[3] * (h * 2. * W2 * vL.x * vL.z)\
		+ (X)->ptr[4] * (hW * CPlus);\
	Y.ptr[2] = (X)->ptr[0] * (hW * vL.y)\
		+ (X)->ptr[1] * (vL.y)\
		+ (X)->ptr[2] * (h * (1. + 2. * W2 * vL.y * vL.y))\
		+ (X)->ptr[3] * (h * (2. * W2 * vL.y * vL.z))\
		+ (X)->ptr[4] * (hW * vL.y);\
	Y.ptr[3] = (X)->ptr[0] * (hW * vL.z)\
		+ (X)->ptr[1] * (vL.z)\
		+ (X)->ptr[2] * (h * (2. * W2 * vL.y * vL.z))\
		+ (X)->ptr[3] * (h * (1. + 2. * W2 * vL.z * vL.z))\
		+ (X)->ptr[4] * (hW * vL.z);\
	Y.ptr[4] = (X)->ptr[0] * (hW * ATildeMinus - 1.)\
		+ (X)->ptr[1] * (1. - Kappa / hW)\
		+ (X)->ptr[2] * (W * vL.y * (2. * hW - 1.))\
		+ (X)->ptr[3] * (W * vL.z * (2. * hW - 1.))\
		+ (X)->ptr[4] * (hW * ATildePlus - 1.);\
\
	/* rotate outgoing y's x's into side */\
	Y.S = normal_vecFromNs(n, Y.S);\
\
	for (int i = numWaves; i < numStates; ++i) {\
		Y.ptr[i] = 0;\
	}\
\
	*(result) = Y;\
}

//// MODULE_NAME: eigen_fluxTransform

#define eigen_fluxTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const * const */eig,\
	/*cons_t const * const */X_,\
	/*real3 const */x,\
	/*normal_t const */n\
) {\
<? if false then ?>\
	/* rotate incoming v's in x */\
	X.S = normal_vecDotNs(n, X.S);\
\
	/* TODO do the matrix multiply here */\
\
	/* rotate outgoing y's x's into side */\
	X.S = normal_vecFromNs(n, X.S);\
<? else ?>\
	/* default */\
	eigen_leftTransform(&waves, solver, eig, X_, x, n);\
	<?=eqn:eigenWaveCodePrefix("n", "eig", "x")?>\
<? for j=0,eqn.numWaves-1 do --\
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode("n", "eig", "x", j)?>;\
<? end --\
?>	eigen_rightTransform(result, solver, eig, waves, x, n);\
<? end ?>\
}

//// MODULE_NAME: constrainU
//// MODULE_DEPENDS: coordLen eqn.guiVars.compileTime

kernel void constrainU(
	constant solver_t const * const solver,
	global cons_t * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 const x = cellBuf[index].pos;

	global cons_t * const U = UBuf + index;

	U->D = max(U->D, (real)solver->DMin);
	U->tau = max(U->tau, (real)solver->tauMin);

	U->D = min(U->D, (real)solver->DMax);
	U->tau = min(U->tau, (real)solver->tauMax);
	
	real D = U->D;
	real3 S = U->S;
	real tau = U->tau;

	real3 v = U->v;

	real SLen = coordLen(S, x);
	real PMin = max(SLen - tau - D + SLen * solver->solvePrimVelEpsilon, solver->solvePrimPMinEpsilon);
	real PMax = (solver->heatCapacityRatio - 1.) * tau;
	PMax = max(PMax, PMin);
	real P = .5 * (PMin + PMax);

	for (int iter = 0; iter < solvePrimMaxIter; ++iter) {
		real vLen = SLen / (tau + D + P);
		real vSq = vLen * vLen;
		real W = 1. / sqrt(1. - vSq);
		real eInt = (tau + D * (1. - W) + P * (1. - W*W)) / (D * W);
		real rho = D / W;
		real f = (solver->heatCapacityRatio - 1.) * rho * eInt - P;
		real csSq = (solver->heatCapacityRatio - 1.) * (tau + D * (1. - W) + P) / (tau + D + P);
		real df_dP = vSq * csSq - 1.;
		real newP = P - f / df_dP;
		newP = max(newP, PMin);
		real PError = fabs(1. - newP / P);
		P = newP;
		if (PError < solver->solvePrimStopEpsilon) {
			v = real3_real_mul(S, 1. / (tau + D + P));
			vSq = coordLenSq(v, x);
			W = 1. / sqrt(1. - vSq);
			rho = D / W;
			rho = max(rho, (real)solver->rhoMin);
			rho = min(rho, (real)solver->rhoMax);
			eInt = P / (rho * (solver->heatCapacityRatio - 1.));
			eInt = min(eInt, (real)solver->eIntMax);
			U->rho = rho;
			U->v = v;
			U->eInt = eInt;
/* printf("cell %d finished with prims = %f %f %f\n", index, rho, v.x, eInt); */
			return;
		}
	}
}

//// MODULE_NAME: addSource

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

<? if not require "hydro.coord.cartesian".is(solver.coord) then ?>
	/* connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system */
	/* TODO calculate this according to SRHD flux.  I'm winging it right now. */
	real const P = calc_P(solver, U->rho, U->eInt);
	real3 const Ftau = real3_sub(U->S, real3_real_mul(U->v, U->D));

	/* - Conn^i_jk S^j v^k  */
	deriv->S = real3_sub(deriv->S, coord_conn_apply23(U->S, U->v, x));	
	
	/* - Conn^i_jk g^jk P */
	deriv->S = real3_sub(deriv->S, real3_real_mul(coord_conn_trace23(x), P));
	
	/* + (gamma-1) rho v^k v^l Gamma_kjl g^ij */
	deriv->S = real3_add(deriv->S, real3_real_mul(coord_conn_apply13(U->v, U->S, x), (solver->heatCapacityRatio - 1.) ));	
	
	/* - (gamma-1) rho v^j v^k v^l Gamma_jkl */
/* 	deriv->ETotal -= (solver->heatCapacityRatio - 1.) * coord_conn_apply123(U->v, U->v, U->S, x);	 */
<? end ?>
}
