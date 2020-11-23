//// MODULE_NAME: primFromCons
//// MODULE_DEPENDS: real3 solver_t prim_t cons_t eqn.common
// eqn.common is for all the calc_* stuff

#define primFromCons(\
	/*prim_t * const */resultW,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const*/U,\
	/*real3 const */x\
) {\
	(resultW)->rho = (U)->rho;\
	(resultW)->v = real3_real_mul((U)->m, 1./(U)->rho);\
	(resultW)->P = calc_P(solver, U, x);\
	(resultW)->ePot = (U)->ePot;\
}

//// MODULE_NAME: consFromPrim
//// MODULE_DEPENDS: real3 solver_t prim_t cons_t eqn.common
// eqn.common is for all the calc_* stuff

#define consFromPrim(\
	/*cons_t * const */ result,\
	/*constant solver_t const * const */solver,\
	/*prim_t const * const */W,\
	/*real3 const */pt\
) {\
	(result)->rho = (W)->rho;\
	(result)->m = real3_real_mul((W)->v, (W)->rho);\
	(result)->ETotal = calc_ETotal(solver, W, pt);\
	(result)->ePot = (W)->ePot;\
}

//// MODULE_NAME: apply_dU_dW	
//// MODULE_DEPENDS: solver_t prim_t cons_t coord_lower
// only used by PLM

<?=eqn.cons_t?> apply_dU_dW(
	constant solver_t const * const solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	real3 WA_vL = coord_lower(WA.v, x);
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_add(
			real3_real_mul(WA.v, W.rho), 
			real3_real_mul(W.v, WA.rho)),
		.ETotal = W.rho * .5 * real3_dot(WA.v, WA_vL) 
			+ WA.rho * real3_dot(W.v, WA_vL)
			+ W.P / (solver->heatCapacityRatio - 1.),
		.ePot = W.ePot,
	};
}

//// MODULE_NAME: apply_dW_dU	
//// MODULE_DEPENDS: solver_t prim_t cons_t coord_lower
// only used by PLM

<?=eqn.prim_t?> apply_dW_dU(
	constant solver_t const * const solver,
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 WA_vL = coord_lower(WA.v, x);
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_sub(
			real3_real_mul(U.m, 1. / WA.rho),
			real3_real_mul(WA.v, U.rho / WA.rho)),
		.P = (solver->heatCapacityRatio - 1.) * (
			.5 * real3_dot(WA.v, WA_vL) * U.rho 
			- real3_dot(U.m, WA_vL)
			+ U.ETotal),
		.ePot = U.ePot,
	};
}

//// MODULE_NAME: eqn.common
//// MODULE_DEPENDS: cons_t prim_t waves_t eigen_t coordLenSq

#define /*real*/ calc_H(/*constant solver_t const * const */solver, /*real const */P) 										((P) * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)))
#define /*real*/ calc_h(/*constant solver_t const * const */solver, /*real const */rho, /*real const */P) 					(calc_H(solver, P) / (rho))
#define /*real*/ calc_HTotal(/*real const */P, /*real const */ETotal) 														((P) + (ETotal))
#define /*real*/ calc_hTotal(/*real const */rho, /*real const */P, /*real const */ETotal) 									(calc_HTotal(P, ETotal) / (rho))
#define /*real*/ calc_eKin(/*prim_t const * const */W, /*real3 const */x) 													(.5 * coordLenSq((W)->v, x))
#define /*real*/ calc_EKin(/*prim_t const * const */W, /*real3 const */x) 													((W)->rho * calc_eKin(W, x))
#define /*real*/ calc_EInt(/*constant solver_t const * const */solver, /*prim_t const * const */W) 							((W)->P / (solver->heatCapacityRatio - 1.))
#define /*real*/ calc_eInt(/*constant solver_t const * const */solver, /*prim_t const * const */W) 							(calc_EInt(solver, W) / (W)->rho)
#define /*real*/ calc_EKin_fromCons(/*cons_t const * const*/U, /*real3 const */x) 											(.5 * coordLenSq((U)->m, x) / (U)->rho)
#define /*real*/ calc_ETotal(/*constant solver_t const * const */solver, /*prim_t const * const */W, /*real3 const */x) 	(calc_EKin(W, x) + calc_EInt(solver, W))
#define /*real*/ calc_Cs(/*constant solver_t const * const */solver, /*prim_t const * const */W) 							(sqrt(solver->heatCapacityRatio * (W)->P / (W)->rho))
#define /*real*/ calc_P(/*constant solver_t const * const */solver, /*cons_t const * const*/U, /*real3 */x) 				(solver->heatCapacityRatio - 1.) * (/*EInt=*/(U)->ETotal - /*EKin=*/calc_EKin_fromCons(U, x))\


//// MODULE_NAME: applyInitCond
//// MODULE_DEPENDS: consFromPrim cartesianToCoord

/*
I've highjacked all of this.  It was a normal Euler eqn solver.
But I experimented with a curved-space solver.  
To get back to the original code,
just replace all the g_ab stuff with their constant values and simplify away.
*/

kernel void applyInitCond(
	constant solver_t const * const solver,
	constant initCond_t const * const initCond,
	global cons_t* UBuf,
	const global cell_t* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real3 mids = real3_real_mul(real3_add(solver->initCondMins, solver->initCondMaxs), .5);
	bool lhs = true<?
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
	
	<?=eqn.prim_t?> W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
		.ePot = ePot,
	};

	consFromPrim(U, solver, &W, x);
}

//// MODULE_NAME: fluxFromCons
//// MODULE_DEPENDS: solver_t primFromCons normal_t

#define fluxFromCons(\
	/*<?=eqn.cons_t?> * const */ resultF,\
	/*constant solver_t const * const*/ solver,\
	/*cons_t const * const */ U,\
	/*real3*/ x,\
	/*normal_t*/ n\
) {\
	<?=eqn.prim_t?> W;\
	primFromCons(&W, solver, U, x);\
	real const v_n = normal_vecDotN1(n, W.v);\
	(resultF)->rho = (U)->rho * v_n;\
	(resultF)->m = real3_add(\
		real3_real_mul((U)->m, v_n),\
		_real3(\
			normal_u1x(n) * W.P,\
			normal_u1y(n) * W.P,\
			normal_u1z(n) * W.P\
		)\
	);\
	real const HTotal = (U)->ETotal + W.P;\
	(resultF)->ETotal = HTotal * v_n;\
	(resultF)->ePot = 0;\
}

//// MODULE_NAME: calcCellMinMaxEigenvalues
//// MODULE_DEPENDS: real3x3 primFromCons
// added by request only, so I don't have to compile the real3x3 code. 
// not used at the moment

#define calcCellMinMaxEigenvalues(\
	/*range_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*global cons_t const * const */U,\
	/*real3 const */pt,\
	/*real3x3 const */nL,\
	/*real3x3 const */nU,\
	/*real const */nLen\
) {\
	prim_t W;\
	primFromCons(&W, solver, U, pt);\
	real const v_n = real3_dot(W.v, nL.x);\
	real const Cs = calc_Cs(solver, &W);\
	real const Cs_nLen = Cs * nLen;\
	(result)->min = v_n - Cs_nLen; \
	(result)->max = v_n + Cs_nLen;\
}

//// MODULE_NAME: eigen_forCell
//// MODULE_DEPENDS: normal_t coord_lower cons_t prim_t eigen_t primFromCons eqn.common
// eqn.common is for all the calc_* stuff

// used by PLM
#define eigen_forCell(\
	/*eigen_t * const */resultEig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */U,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	prim_t W;\
	primFromCons(&W, solver, U, pt);\
	real3 const vL = coord_lower(W.v, pt);\
	real const vSq = real3_dot(W.v, vL);\
	real const v_n = normal_vecDotN1(n, W.v);\
	real const eKin = .5 * vSq;\
	real const hTotal = calc_hTotal(W.rho, W.P, (U)->ETotal);\
	real const CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);\
	real const Cs = sqrt(CsSq);\
	(resultEig)->rho = W.rho;\
	(resultEig)->v = W.v;\
	(resultEig)->vSq = vSq;\
	(resultEig)->vL = vL;\
	(resultEig)->hTotal = hTotal;\
	(resultEig)->Cs = Cs;\
}

//// MODULE_NAME: eigen_forInterface
//// MODULE_DEPENDS: primFromCons eigen_t normal_t coord_lower

//used by the mesh version
#define eigen_forInterface(\
	/*eigen_t * const */resultEig,\
	/*constant solver_t const * const */solver,\
	/*cons_t const * const */UL,\
	/*cons_t const * const */UR,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	prim_t WL;\
	primFromCons(&WL, solver, UL, pt);\
	real sqrtRhoL = sqrt(WL.rho);\
	real3 vLeft = WL.v;\
	real hTotalL = calc_hTotal(WL.rho, WL.P, (UL)->ETotal);\
\
	prim_t WR;\
	primFromCons(&WR, solver, UR, pt);\
	real sqrtRhoR = sqrt(WR.rho);\
	real3 vR = WR.v;\
	real hTotalR = calc_hTotal(WR.rho, WR.P, (UR)->ETotal);\
\
	real invDenom = 1./(sqrtRhoL + sqrtRhoR);\
\
	/*Roe-averaged*/\
	real rho = sqrtRhoL * sqrtRhoR;\
	real3 v = real3_add(\
			real3_real_mul(vLeft, sqrtRhoL * invDenom),\
			real3_real_mul(vR, sqrtRhoR * invDenom));\
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);\
\
	/*derived:*/\
	real3 vLower = coord_lower(v, pt);\
	real vSq = real3_dot(v, vLower);\
	real eKin = .5 * vSq;\
	real CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);\
	real Cs = sqrt(CsSq);\
\
	(resultEig)->rho = rho;\
	(resultEig)->v = v;\
	(resultEig)->vSq = vSq;\
	(resultEig)->vL = vLower;\
	(resultEig)->hTotal = hTotal;\
	(resultEig)->Cs = Cs;\
}

//// MODULE_NAME: eigen_left/rightTransform
//// MODULE_DEPENDS: eigen_t normal_t

#define eigen_leftTransform(\
	/*waves_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const */eig,\
	/*cons_t const */X_,\
	/*real3 const */pt,\
	/*normal_t */n\
) {\
	real3 const v_n = normal_vecDotNs(n, (eig).v);\
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
	real const denom = 2. * (eig).Cs * (eig).Cs;\
	real const invDenom = 1. / denom;\
\
	real const gamma_1 = solver->heatCapacityRatio - 1.;\
	(result)->ptr[0] = (\
			(X_).ptr[0] * (.5 * gamma_1 * (eig).vSq + (eig).Cs * v_n.x * inv_nLen)\
			+ (X_).ptr[1] * (-gamma_1 * (eig).vL.x - (eig).Cs * normal_l1x_over_len(n))\
			+ (X_).ptr[2] * (-gamma_1 * (eig).vL.y - (eig).Cs * normal_l1y_over_len(n))\
			+ (X_).ptr[3] * (-gamma_1 * (eig).vL.z - (eig).Cs * normal_l1z_over_len(n))\
			+ (X_).ptr[4] * gamma_1\
		) * invDenom;\
	(result)->ptr[1] =\
		(\
			(X_).ptr[0] * (denom - gamma_1 * (eig).vSq)\
			+ (X_).ptr[1] * 2. * gamma_1 * (eig).vL.x\
			+ (X_).ptr[2] * 2. * gamma_1 * (eig).vL.y\
			+ (X_).ptr[3] * 2. * gamma_1 * (eig).vL.z\
			+ (X_).ptr[4] * -2. * gamma_1\
		) * invDenom;\
	(result)->ptr[2] =\
		(X_).ptr[0] * -v_n.y\
		+ (X_).ptr[1] * normal_l2x(n)\
		+ (X_).ptr[2] * normal_l2y(n)\
		+ (X_).ptr[3] * normal_l2z(n);\
	(result)->ptr[3] =\
		(X_).ptr[0] * -v_n.z\
		+ (X_).ptr[1] * normal_l3x(n)\
		+ (X_).ptr[2] * normal_l3y(n)\
		+ (X_).ptr[3] * normal_l3z(n);\
	(result)->ptr[4] =\
		(\
			(X_).ptr[0] * (.5 * gamma_1 * (eig).vSq - (eig).Cs * v_n.x * inv_nLen)\
			+ (X_).ptr[1] * (-gamma_1 * (eig).vL.x + (eig).Cs * normal_l1x_over_len(n))\
			+ (X_).ptr[2] * (-gamma_1 * (eig).vL.y + (eig).Cs * normal_l1y_over_len(n))\
			+ (X_).ptr[3] * (-gamma_1 * (eig).vL.z + (eig).Cs * normal_l1z_over_len(n))\
			+ (X_).ptr[4] * gamma_1\
		) * invDenom;\
}

#define eigen_rightTransform(\
	/*cons_t * const */result,\
	/*constant solver_t const * const */solver,\
	/*eigen_t const */eig,\
	/*waves_t const */X_,\
	/*real3 const */pt,\
	/*normal_t const */n\
) {\
	real3 const v_n = normal_vecDotNs(n, (eig).v);\
	real const nLen = normal_len(n);\
	real const inv_nLen = 1. / nLen;\
	(result)->ptr[0] =\
		(X_).ptr[0]\
		+ (X_).ptr[1]\
		+ (X_).ptr[4];\
	(result)->ptr[1] =\
		(X_).ptr[0] * ((eig).v.x - (eig).Cs * normal_u1x_over_len(n))\
		+ (X_).ptr[1] * (eig).v.x\
		+ (X_).ptr[2] * normal_u2x(n)\
		+ (X_).ptr[3] * normal_u3x(n)\
		+ (X_).ptr[4] * ((eig).v.x + (eig).Cs * normal_u1x_over_len(n));\
	(result)->ptr[2] =\
		(X_).ptr[0] * ((eig).v.y - (eig).Cs * normal_u1y_over_len(n))\
		+ (X_).ptr[1] * (eig).v.y\
		+ (X_).ptr[2] * normal_u2y(n)\
		+ (X_).ptr[3] * normal_u3y(n)\
		+ (X_).ptr[4] * ((eig).v.y + (eig).Cs * normal_u1y_over_len(n));\
	(result)->ptr[3] =\
		(X_).ptr[0] * ((eig).v.z - (eig).Cs * normal_u1z_over_len(n))\
		+ (X_).ptr[1] * (eig).v.z\
		+ (X_).ptr[2] * normal_u2z(n)\
		+ (X_).ptr[3] * normal_u3z(n)\
		+ (X_).ptr[4] * ((eig).v.z + (eig).Cs * normal_u1z_over_len(n));\
	(result)->ptr[4] =\
		(X_).ptr[0] * ((eig).hTotal - (eig).Cs * v_n.x * inv_nLen)\
		+ (X_).ptr[1] * .5 * (eig).vSq\
		+ (X_).ptr[2] * v_n.y\
		+ (X_).ptr[3] * v_n.z\
		+ (X_).ptr[4] * ((eig).hTotal + (eig).Cs * v_n.x * inv_nLen);\
	(result)->ptr[0] =\
		0;\
}

//// MODULE_NAME: eigen_fluxTransform
//// MODULE_DEPENDS: eigen_t normal_t
// not used anymore.  was used by Roe, but I switched that to a fluxFromCons.
// is this needed anymore?
// fluxFromCons only matches eigen_fluxTransform when the eig properties are derived from X_ 

cons_t eigen_fluxTransform(
	constant solver_t const * const solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	normal_t n
) {
	real* X = X_.ptr;
	real3 const v_n = normal_vecDotNs(n, (eig).v);\
	real const nLen = normal_len(n);\
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	return (cons_t){.ptr={
		X[1] * normal_l1x(n)
		+ X[2] * normal_l1y(n)
		+ X[3] * normal_l1z(n),
		
		X[0] * (-v_n.x * (eig).v.x + gamma_1 * .5 * (eig).vSq * normal_u1x(n))
		+ X[1] * ((eig).v.x * normal_l1x(n) - gamma_2 * normal_u1x(n) * (eig).vL.x + v_n.x)
		+ X[2] * ((eig).v.x * normal_l1y(n) - gamma_2 * normal_u1x(n) * (eig).vL.y)
		+ X[3] * ((eig).v.x * normal_l1z(n) - gamma_2 * normal_u1x(n) * (eig).vL.z)
		+ X[4] * gamma_1 * normal_u1x(n),
		
		X[0] * (-v_n.x * (eig).v.y + gamma_1 * .5 * (eig).vSq * normal_u1y(n))
		+ X[1] * ((eig).v.y * normal_l1x(n) - gamma_2 * normal_u1y(n) * (eig).vL.x)
		+ X[2] * ((eig).v.y * normal_l1y(n) - gamma_2 * normal_u1y(n) * (eig).vL.y + v_n.x)
		+ X[3] * ((eig).v.y * normal_l1z(n) - gamma_2 * normal_u1y(n) * (eig).vL.z)
		+ X[4] * gamma_1 * normal_u1y(n),
		
		X[0] * (-v_n.x * (eig).v.z + gamma_1 * .5 * (eig).vSq * normal_u1z(n))
		+ X[1] * ((eig).v.z * normal_l1x(n) - gamma_2 * normal_u1z(n) * (eig).vL.x)
		+ X[2] * ((eig).v.z * normal_l1y(n) - gamma_2 * normal_u1z(n) * (eig).vL.y)
		+ X[3] * ((eig).v.z * normal_l1z(n) - gamma_2 * normal_u1z(n) * (eig).vL.z + v_n.x)
		+ X[4] * gamma_1 * normal_u1z(n),
		
		X[0] * v_n.x * (.5 * gamma_1 * (eig).vSq - (eig).hTotal)
		+ X[1] * (normal_l1x(n) * (eig).hTotal - gamma_1 * v_n.x * (eig).vL.x)
		+ X[2] * (normal_l1y(n) * (eig).hTotal - gamma_1 * v_n.x * (eig).vL.y)
		+ X[3] * (normal_l1z(n) * (eig).hTotal - gamma_1 * v_n.x * (eig).vL.z)
		+ X[4] * gamma * v_n.x,
		
		0,
	}};
}

<? if false then -- TODO sort addSource out. ?>
//// MODULE_NAME: addSource
//// MODULE_DEPENDS: solver_t cons_t cell_t SETBOUNDS_NOGHOST

kernel void addSource(
	constant solver_t const * const solver,
	global cons_t * const derivBuf,
	global cons_t const * const UBuf,
	global cell_t const * const cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 const x = cellBuf[index].pos;

	global cons_t const * const deriv = derivBuf + index;
	global cons_t const * const U = UBuf + index;

<? if false 
and solver.coord.vectorComponent == 'anholonomic' 
and require 'hydro.coord.cylinder'.is(solver.coord) 
then ?>
<? 	if true then -- 2009 Trangenstein, p.474, 1999 Toro, p.29, eqn.1.104, 1.105 ?>
	<? for side=0,1 do ?>{
		real3 xL = x; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		real3 xR = x; xR.s<?=side?> += solver->grid_dx.s<?=side?>;
		
		global cons_t const * const UL = U - solver->stepsize.s<?=side?>;
		global cons_t const * const UR = U + solver->stepsize.s<?=side?>;
		real const PL = calc_P(solver, UL, xL);
		real const PR = calc_P(solver, UR, xR);
	
		deriv->m.s<?=side?> -= (PR - PL) / (2. * solver->grid_dx.s<?=side?>);
	}<? end ?>
<?	end ?>
<?	if false then -- 1999 Toro p.28 eqn.1.102, 1.103 ?>
	cons_t F;
	fluxFromConsForSide0(&F, solver, U, x);
	deriv->rho -= F.rho / x.x;
	deriv->m.x -= F.m.x / x.x;
	deriv->m.y -= F.m.y / x.x;
	deriv->m.z -= F.m.z / x.x;
	deriv->ETotal -= F.ETotal / x.x;
<?	end ?>
<? end ?>

<? do -- if not solver.coord.vectorComponent == 'anholonomic' then ?>
<? if not (require 'hydro.coord.cartesian'.is(solver.coord) 
		or solver.coord.vectorComponent == 'cartesian')
then ?>
//// MODULE_DEPENDS: primFromCons coord_conn_apply23 coord_conn_trace23 coord_conn_apply13
/*
This is working for init conds with zero velocity.
Introducing constant velocity of v=[x=1,y=1]=[r=sqrt(2),theta=pi/4] in the init cond causes some numerical errors.
However the problem isn't the terms below -- because disabling this for zero-vel init conds causes things to become unsteady.
That means that the volume gradient in calcDerivFV is causing nonzero velocities to emerge, and this is cancelling them.
Maybe for an initial constant vel as large as sqrt(2) this fails, but it works only for small perturbations?
*/
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	prim_t W;
	primFromCons(&W, solver, U, x);
	
	//- Conn^i_jk rho v^j v^k 
	deriv->m = real3_sub(deriv->m, coord_conn_apply23(W.v, U->m, x));	
	
	//- Conn^i_jk g^jk P
	deriv->m = real3_sub(deriv->m, real3_real_mul(coord_conn_trace23(x), W.P));		
	
	//+ (gamma-1) rho v^k v^l Gamma_kjl g^ij
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_conn_apply13(W.v, U->m, x), (solver->heatCapacityRatio - 1.) ));	
	
	//- (gamma-1) rho v^j v^k v^l Gamma_jkl
//	deriv->ETotal -= (solver->heatCapacityRatio - 1.) * coord_conn_apply123(W.v, W.v, U->m, x);	

	//+ c_jk^k * Flux^Ij
<? 	if false and solver.coord.vectorComponent == 'anholonomic' then ?>
	real3 const commTrace = coord_tr23_c(x);
	<? for i=0,solver.dim-1 do ?>{
		cons_t flux;
		calcFluxFromCons(&F, *U, x);
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] += commTrace.s<?=i?> * flux.ptr[j];
		}
	}<? end ?>
<? 	end ?>
<? end ?>
<? end -- vectorComponent == 'anholonomic' ?>
}
<? end ?>

//// MODULE_NAME: constrainU
//// MODULE_DEPENDS: solver_t cons_t cell_t primFromCons consFromPrim

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

	if (W.rho < solver->rhoMin) W.rho = solver->rhoMin;
	if (W.P < solver->PMin) W.P = solver->PMin;

	consFromPrim(U, solver, &W, x);
}
