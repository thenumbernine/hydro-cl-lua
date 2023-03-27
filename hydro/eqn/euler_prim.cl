//// MODULE_NAME: <?=eqn_common?>
//// MODULE_DEPENDS: <?=coordLenSq?> <?=cons_t?> <?=waves_t?> <?=eigen_t?> <?=eqn_guiVars_compileTime?>

#define /*real*/ <?=calc_H?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*real const */P\
)	((P) * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)))

#define /*real*/ <?=calc_h?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*real const */rho,\
	/*real const */P\
)	(<?=calc_H?>(solver, P) / (rho))

#define /*real*/ <?=calc_HTotal?>(\
	/*real const */P,\
	/*real const */ETotal\
)	((P) + (ETotal))

#define /*real*/ <?=calc_hTotal?>(\
	/*real const */rho,\
	/*real const */P,\
	/*real const */ETotal\
)	(<?=calc_HTotal?>(P, ETotal) / (rho))

#define /*real*/ <?=calc_eKin?>(\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x\
)	(.5 * coordLenSq((U)->v, x))

#define /*real*/ <?=calc_EKin?>(\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x\
) 	((U)->rho * <?=calc_eKin?>(U, x))

#define /*real*/ <?=calc_EInt?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
) 	((U)->P / (solver->heatCapacityRatio - 1.))

#define /*real*/ <?=calc_eInt?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
)	(<?=calc_EInt?>(solver, U) / (U)->rho)

#define /*real*/ <?=calc_ETotal?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*real3 const */x\
) 	(<?=calc_EKin?>(U, x) + <?=calc_EInt?>(solver, U))

#define /*real*/ <?=calc_Cs?>(\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U\
) (sqrt(solver->heatCapacityRatio * (U)->P / (U)->rho))

<? local materials = require "hydro.materials" ?>
#define C_v				<?=("%.50f"):format(materials.Air.C_v)?>

#define /*real*/ <?=calc_T?>(\
	/*<?=solver_t?> const * const*/solver,\
	/*<?=cons_t?> const * const*/U\
) (<?=calc_eInt?>(solver, U) / C_v)

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
	
	*U = (<?=cons_t?>){
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
	};
}

//// MODULE_NAME: <?=calcCellMinMaxEigenvalues?>
//// MODULE_DEPENDS: real3x3
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
	real const v_n = (U)->v.dot(nL.x);\
	real const Cs = <?=calc_Cs?>(solver, U);\
	real const Cs_nLen = Cs * nLen;\
	(result)->min = v_n - Cs_nLen; \
	(result)->max = v_n + Cs_nLen;\
}

//// MODULE_NAME: <?=eigen_forCell?>
//// MODULE_DEPENDS: <?=normal_t?> <?=coord_lower?> <?=cons_t?> <?=eigen_t?> <?=eqn_common?>
// eqn_common is for all the calc_* stuff

// used by PLM
#define <?=eigen_forCell?>(\
	/*<?=eigen_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
	real3 const vL = coord_lower((U)->v, (cell)->pos);\
	real const vSq = (U)->v.dot(vL);\
	real const v_n = normal_vecDotN1(n, (U)->v);\
	real const eKin = .5 * vSq;\
	real const ETotal = <?=calc_ETotal?>(solver, U, cell->pos);\
	real const hTotal = <?=calc_hTotal?>((U)->rho, (U)->P, ETotal);\
	real const CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);\
	real const Cs = sqrt(CsSq);\
	(result)->rho = (U)->rho;\
	(result)->v = (U)->v;\
	(result)->vSq = vSq;\
	(result)->vL = vL;\
	(result)->hTotal = hTotal;\
	(result)->Cs = Cs;\
}

//// MODULE_NAME: <?=eigen_forInterface?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?> <?=coord_lower?>

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
	real const sqrtRhoL = sqrt((UL)->rho);\
	real3 const vLeft = (UL)->v;\
	real const ETotalL = <?=calc_ETotal?>(solver, UL, cellL->pos);\
	real const hTotalL = <?=calc_hTotal?>((UL)->rho, (UL)->P, ETotalL);\
\
	real const sqrtRhoR = sqrt((UR)->rho);\
	real3 const vR = (UR)->v;\
	real const ETotalR = <?=calc_ETotal?>(solver, UR, cellR->pos);\
	real const hTotalR = <?=calc_hTotal?>((UR)->rho, (UR)->P, ETotalR);\
\
	real const invDenom = 1./(sqrtRhoL + sqrtRhoR);\
\
	/*Roe-averaged*/\
	(result)->rho = sqrtRhoL * sqrtRhoR;\
	real3 const v = real3_add(\
			real3_real_mul(vLeft, sqrtRhoL * invDenom),\
			real3_real_mul(vR, sqrtRhoR * invDenom));\
	real const hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);\
\
	/*derived:*/\
	real3 const vLower = coord_lower(v, pt);\
	real const vSq = v.dot(vLower);\
	real const eKin = .5 * vSq;\
	real const h = hTotal - eKin;\
	/* TODO verify hTotal = 1/2 v^2 + Cs^2 / (gamma-1) */\
	(result)->hTotal = hTotal;\
	real const CsSq = (solver->heatCapacityRatio - 1.) * h;\
	(result)->Cs = sqrt(CsSq);\
\
	(result)->v = v;\
	(result)->vSq = vSq;\
	(result)->vL = vLower;\
}

//// MODULE_NAME: <?=eigen_leftTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>

#define <?=eigen_leftTransform?>(\
	/*<?=waves_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=cons_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> */n\
) {\
	real const nLen = normal_len(n);\
	real const nLenSq = normal_lenSq(n);\
	real const CsSq = (eig)->Cs * (eig)->Cs;\
	(result)->ptr[0] = (\
			  (X)->ptr[1] * -(eig)->rho * normal_l1x_over_len(n)\
			+ (X)->ptr[2] * -(eig)->rho * normal_l1y_over_len(n)\
			+ (X)->ptr[3] * -(eig)->rho * normal_l1z_over_len(n)\
			+ (X)->ptr[4] / (eig)->Cs\
		) * .5 / ((eig)->Cs * nLen);\
	(result)->ptr[1] =\
			  (X)->ptr[0]\
			+ (X)->ptr[4] * -1. / CsSq;\
	(result)->ptr[2] =(\
			  (X)->ptr[1] * -(nLenSq - normal_l1y(n) * normal_u1y(n)) / normal_l1x(n)\
			+ (X)->ptr[2] * normal_u1y(n) / normal_l1x(n)\
			+ (X)->ptr[3] * -normal_l1z(n)\
		) * (eig)->rho / nLenSq;\
	(result)->ptr[3] =(\
			  (X)->ptr[1] * -normal_u1z(n)\
			+ (X)->ptr[2] * -normal_l1y(n) * normal_u1z(n) / normal_l1x(n)\
			+ (X)->ptr[3] * (nLenSq - normal_l1z(n) * normal_u1z(n)) / normal_l1x(n)\
		) * (eig)->rho / nLenSq;\
	(result)->ptr[4] =(\
			  (X)->ptr[1] * (eig)->rho * normal_l1x_over_len(n)\
			+ (X)->ptr[2] * (eig)->rho * normal_l1y_over_len(n)\
			+ (X)->ptr[3] * (eig)->rho * normal_l1z_over_len(n)\
			+ (X)->ptr[4] / (eig)->Cs\
		) * .5 / ((eig)->Cs * nLen);\
}

//// MODULE_NAME: <?=eigen_rightTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>

#define <?=eigen_rightTransform?>(\
	/*<?=cons_t?> * const */result,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=eigen_t?> const * const */eig,\
	/*<?=waves_t?> const * const */X,\
	/*real3 const */pt,\
	/*<?=normal_t?> const */n\
) {\
	real const nLen = normal_len(n);\
	real const invRho = 1. / (eig)->rho;\
	(result)->ptr[0] =\
		  (X)->ptr[0] * nLen\
		+ (X)->ptr[1]\
		+ (X)->ptr[4] * nLen;\
	(result)->ptr[1] =(\
		  (X)->ptr[0] * -(eig)->Cs * normal_u1x(n)\
		+ (X)->ptr[2] * -normal_l1y(n)\
		+ (X)->ptr[3] * -normal_l1z(n)\
		+ (X)->ptr[4] * (eig)->Cs * normal_u1x(n)\
	) * invRho;\
	(result)->ptr[2] =(\
		  (X)->ptr[0] * -(eig)->Cs * normal_u1y(n)\
		+ (X)->ptr[2] * normal_l1x(n)\
		+ (X)->ptr[4] * (eig)->Cs * normal_u1y(n)\
	) * invRho;\
	(result)->ptr[3] =(\
		  (X)->ptr[0] * -(eig)->Cs * normal_u1z(n)\
		+ (X)->ptr[3] * normal_l1x(n)\
		+ (X)->ptr[4] * (eig)->Cs * normal_u1z(n)\
	) * invRho;\
	(result)->ptr[4] =(\
		  (X)->ptr[0]\
		+ (X)->ptr[4]\
	) * nLen * (eig)->Cs * (eig)->Cs;\
}

//// MODULE_NAME: <?=eigen_fluxTransform?>
//// MODULE_DEPENDS: <?=eigen_t?> <?=normal_t?>
// Not used anymore.  was used by Roe, but I switched that to a <?=fluxFromCons?>.
// <?=fluxFromCons?> only matches <?=eigen_fluxTransform?> when the eig properties are derived from X_ 
// TODO keep including it?  is all of the Roe solver based on the idea that F = dF/dU * U = R * Lambda * L * U?  I think so.
//  in that case ... this can't work anyways.


<? if false then -- TODO sort <?=addSource?> out. ?>
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
		real const PL = UL->P;
		real const PR = UR->P;
	
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

<? if false then -- if not solver.coord.vectorComponent == "anholonomic" then ?>
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
	
	//- Γ^i_jk ρ v^j v^k 
	deriv->m = real3_sub(deriv->m, coord_conn_apply23(U->v, U->m, x));	
	
	//- Γ^i_jk g^jk P
	deriv->m = real3_sub(deriv->m, real3_real_mul(coord_conn_trace23(x), U->P));		
	
	//+ (γ-1) ρ v^k v^l Γ_kjl g^ij
	deriv->m = real3_add(deriv->m, real3_real_mul(coord_conn_apply13(U->v, U->m, x), (solver->heatCapacityRatio - 1.) ));	
	
	//- (γ-1) ρ v^j v^k v^l Γ_jkl
//	deriv->ETotal -= (solver->heatCapacityRatio - 1.) * coord_conn_apply123(U->v, U->v, U->m, x);	

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
}
<? end ?>
