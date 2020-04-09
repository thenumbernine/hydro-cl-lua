/*
I've highjacked all of this.  It was a normal Euler eqn solver.
But I experimented with a curved-space solver.  
To get back to the original code,
just replace all the g_ab stuff with their constant values and simplify away.
*/
<? local solver = eqn.solver ?>

//as long as no subsequent solvers' codes are tacked onto the end of this,
//I can safely do this:
typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

cons_t fluxFromConsForNormal(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	real3x3 nL,	//covariant, row-vector basis
	real3x3 nU,
	real nLen
) {
	prim_t W = primFromCons(solver, U, x);
	real v_n = real3_dot(W.v, nL.x);
	real HTotal = U.ETotal + W.P;
	
	return (cons_t){
		.rho = U.rho * v_n,
		.m = real3_add(
			real3_real_mul(U.m, v_n),
			real3_real_mul(nU.x, W.P)
		),
		.ETotal = HTotal * v_n,
		.ePot = 0,
	};
}


// used by PLM
range_t calcCellMinMaxEigenvaluesForNormal(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) {
	prim_t W = primFromCons(solver, *U, x);
	real v_n = real3_dot(W.v, nL.x);
	real Cs = calc_Cs(solver, &W);
	real Cs_nLen = Cs * nLen;
	return (range_t){
		.min = v_n - Cs_nLen, 
		.max = v_n + Cs_nLen,
	};
}


//used by the mesh version
eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) {
	prim_t WL = primFromCons(solver, UL, x);
	real sqrtRhoL = sqrt(WL.rho);
	real3 vLeft = WL.v;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL.ETotal);
	
	prim_t WR = primFromCons(solver, UR, x);
	real sqrtRhoR = sqrt(WR.rho);
	real3 vR = WR.v;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR.ETotal);

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);
	
	//Roe-averaged
	real rho = sqrtRhoL * sqrtRhoR;
	real3 v = real3_add(
			real3_real_mul(vLeft, sqrtRhoL * invDenom),
			real3_real_mul(vR, sqrtRhoR * invDenom));
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);

	//derived:
	real3 vLower = coord_lower(v, x);
	real vSq = real3_dot(v, vLower);
	real eKin = .5 * vSq;
	real CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);

	//n is row one-forms, v is column vector
	real3 v_n = real3x3_real3_mul(nL, v);

	return (eigen_t){
		.rho = rho, 
		.v = v,
		.hTotal = hTotal,
		.vSq = vSq,
		.Cs = Cs,
		.vL = vLower,
		.v_n = v_n,
	};
}

// used by PLM
eigen_t eigen_forCellForNormal(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) {
	prim_t W = primFromCons(solver, U, x);
	real3 vL = coord_lower(W.v, x);
	real vSq = real3_dot(W.v, vL);
	real v_n = real3_dot(W.v, nL.x);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
	real CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	return (eigen_t){
		.rho = W.rho,
		.v = W.v,
		.hTotal = hTotal,
		.vSq = vSq,
		.Cs = Cs,
		.vL = vL,
		.v_n = v_n,
	};
}



<?
local prefix = [[
	real3 v = eig.v;
	real hTotal = eig.hTotal;
	real Cs = eig.Cs;
	real vSq = eig.vSq;
	real3 v_n = eig.v_n;
	real3 vL = eig.vL;
]]
?>

waves_t eigen_leftTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) { 
	real* X = X_.ptr;
	<?=prefix?>
	
	real denom = 2. * Cs * Cs;
	real invDenom = 1. / denom;

	const real gamma_1 = solver->heatCapacityRatio - 1.;
	return (waves_t){.ptr={
		(
			X[0] * (.5 * gamma_1 * vSq + Cs * v_n.x / nLen)
			+ X[1] * (-gamma_1 * vL.x - Cs * nL.x.x / nLen)
			+ X[2] * (-gamma_1 * vL.y - Cs * nL.x.y / nLen)
			+ X[3] * (-gamma_1 * vL.z - Cs * nL.x.y / nLen)
			+ X[4] * gamma_1
		) * invDenom,
		
		(
			X[0] * (denom - gamma_1 * vSq)
			+ X[1] * 2. * gamma_1 * vL.x
			+ X[2] * 2. * gamma_1 * vL.y
			+ X[3] * 2. * gamma_1 * vL.z
			+ X[4] * -2. * gamma_1
		) * invDenom,
		
		X[0] * -v_n.y
		+ X[1] * nL.y.x
		+ X[2] * nL.y.y
		+ X[3] * nL.y.z,
		
		X[0] * -v_n.z
		+ X[1] * nL.z.x
		+ X[2] * nL.z.y
		+ X[3] * nL.z.z,
		
		(
			X[0] * (.5 * gamma_1 * vSq - Cs * v_n.x / nLen)
			+ X[1] * (-gamma_1 * vL.x + Cs * nL.x.x / nLen)
			+ X[2] * (-gamma_1 * vL.y + Cs * nL.x.y / nLen)
			+ X[3] * (-gamma_1 * vL.z + Cs * nL.x.z / nLen)
			+ X[4] * gamma_1
		) * invDenom,
	}};
}

cons_t eigen_rightTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) {
	real* X = X_.ptr;
	<?=prefix?>
	return (cons_t){.ptr={
		X[0] 
		+ X[1] 
		+ X[4],
		
		X[0] * (v.x - Cs * nU.x.x / nLen)
		+ X[1] * v.x
		+ X[2] * nU.y.y
		+ X[3] * nU.z.z
		+ X[4] * (v.x + Cs * nU.x.x / nLen),
		
		X[0] * (v.y - Cs * nU.x.y / nLen)
		+ X[1] * v.y
		+ X[2] * nU.y.y
		+ X[3] * nU.z.y
		+ X[4] * (v.y + Cs * nU.x.y / nLen),
		
		X[0] * (v.z - Cs * nU.x.z / nLen)
		+ X[1] * v.z
		+ X[2] * nU.y.z
		+ X[3] * nU.z.z
		+ X[4] * (v.z + Cs * nU.x.z / nLen),
		
		X[0] * (hTotal - Cs * v_n.x / nLen)
		+ X[1] * .5 * vSq
		+ X[2] * v_n.y
		+ X[3] * v_n.z
		+ X[4] * (hTotal + Cs * v_n.x / nLen),
		
		0,
	}};
}

cons_t eigen_fluxTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) {
	real* X = X_.ptr;
	<?=prefix?>
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	return (cons_t){.ptr={
		X[1] * nL.x.x
		+ X[2] * nL.x.y
		+ X[3] * nL.x.z,
		
		X[0] * (-v_n.x * v.x + gamma_1 * .5 * vSq * nU.x.x)
		+ X[1] * (v.x * nL.x.x - gamma_2 * nU.x.x * vL.x + v_n.x)
		+ X[2] * (v.x * nL.x.y - gamma_2 * nU.x.x * vL.y)
		+ X[3] * (v.x * nL.x.z - gamma_2 * nU.x.x * vL.z)
		+ X[4] * gamma_1 * nU.x.x,
		
		X[0] * (-v_n.x * v.y + gamma_1 * .5 * vSq * nU.x.y)
		+ X[1] * (v.y * nL.x.x - gamma_2 * nU.x.y * vL.x)
		+ X[2] * (v.y * nL.x.y - gamma_2 * nU.x.y * vL.y + v_n.x)
		+ X[3] * (v.y * nL.x.z - gamma_2 * nU.x.y * vL.z)
		+ X[4] * gamma_1 * nU.x.y,
		
		X[0] * (-v_n.x * v.z + gamma_1 * .5 * vSq * nU.x.z)
		+ X[1] * (v.z * nL.x.x - gamma_2 * nU.x.z * vL.x)
		+ X[2] * (v.z * nL.x.y - gamma_2 * nU.x.z * vL.y)
		+ X[3] * (v.z * nL.x.z - gamma_2 * nU.x.z * vL.z + v_n.x)
		+ X[4] * gamma_1 * nU.x.z,
		
		X[0] * v_n.x * (.5 * gamma_1 * vSq - hTotal)
		+ X[1] * (nL.x.x * hTotal - gamma_1 * v_n.x * vL.x)
		+ X[2] * (nL.x.y * hTotal - gamma_1 * v_n.x * vL.y)
		+ X[3] * (nL.x.z * hTotal - gamma_1 * v_n.x * vL.z)
		+ X[4] * gamma * v_n.x,
		
		0,
	}};
}


kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

<? if false 
and solver.coord.vectorComponent == 'anholonomic' 
and require 'coord.cylinder'.is(solver.coord) 
then ?>
<? 	if true then -- 2009 Trangenstein, p.474, 1999 Toro, p.29, eqn.1.104, 1.105 ?>
	<? for side=0,1 do ?>{
		real3 xL = x; xL.s<?=side?> -= solver->grid_dx.s<?=side?>;
		real3 xR = x; xR.s<?=side?> += solver->grid_dx.s<?=side?>;
		
		real PL = calc_P(solver, U[-solver->stepsize.s<?=side?>], xL);
		real PR = calc_P(solver, U[solver->stepsize.s<?=side?>], xR);
	
		deriv->m.s<?=side?> -= (PR - PL) / (2. * solver->grid_dx.s<?=side?>);
	}<? end ?>
<?	end ?>
<?	if false then -- 1999 Toro p.28 eqn.1.102, 1.103 ?>
	cons_t F = fluxFromConsForSide0(solver, *U, x);
	deriv->rho -= F.rho / x.x;
	deriv->m.x -= F.m.x / x.x;
	deriv->m.y -= F.m.y / x.x;
	deriv->m.z -= F.m.z / x.x;
	deriv->ETotal -= F.ETotal / x.x;
<?	end ?>
<? end ?>

<? do -- if not solver.coord.vectorComponent == 'anholonomic' then ?>
<? if not require 'coord.cartesian'.is(solver.coord) then ?>
/*
This is working for init conds with zero velocity.
Introducing constant velocity of v=[x=1,y=1]=[r=sqrt(2),theta=pi/4] in the init cond causes some numerical errors.
However the problem isn't the terms below -- because disabling this for zero-vel init conds causes things to become unsteady.
That means that the volume gradient in calcDerivFV is causing nonzero velocities to emerge, and this is cancelling them.
Maybe for an initial constant vel as large as sqrt(2) this fails, but it works only for small perturbations?
*/
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	prim_t W = primFromCons(solver, *U, x);
	
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
	real3 commTrace = coord_tr23_c(x);
	<? for i=0,solver.dim-1 do ?>{
		cons_t flux = calcFluxFromCons(*U, x);
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] += commTrace.s<?=i?> * flux.ptr[j];
		}
	}<? end ?>
<? 	end ?>
<? end ?>
<? end -- vectorComponent == 'anholonomic' ?>
}

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf
) {
	SETBOUNDS(numGhost,numGhost);	
	global cons_t* U = UBuf + index;
	real3 x = cell_x(i);
	prim_t W = primFromCons(solver, *U, x);

	if (W.rho < solver->rhoMin) W.rho = solver->rhoMin;
	if (W.P < solver->PMin) W.P = solver->PMin;

	*U = consFromPrim(solver, W, x);
}


<? for side=0,solver.dim-1 do ?>
#define fluxFromConsForSide<?=side?>(solver, U, x) 				fluxFromConsForNormal(solver, U, x, normalBasisForSide<?=side?>, normalBasisForSide<?=side?>, 1)
#define calcCellMinMaxEigenvalues_<?=side?>(solver, U, x) 		calcCellMinMaxEigenvaluesForNormal(solver, U, x, normalBasisForSide<?=side?>, normalBasisForSide<?=side?>, 1)
#define eigen_leftTransformForSide<?=side?>(solver, eig, X, x) 	eigen_leftTransformForNormal(solver, eig, X, x, normalBasisForSide<?=side?>, normalBasisForSide<?=side?>, 1)
#define eigen_rightTransformForSide<?=side?>(solver, eig, X, x) eigen_rightTransformForNormal(solver, eig, X, x, normalBasisForSide<?=side?>, normalBasisForSide<?=side?>, 1)
#define eigen_fluxTransformForSide<?=side?>(solver, eig, X, x) 	eigen_fluxTransformForNormal(solver, eig, X, x, normalBasisForSide<?=side?>, normalBasisForSide<?=side?>, 1)
#define eigen_forCell_<?=side?>(solver, U, x) 					eigen_forCellForNormal(solver, U, x, normalBasisForSide<?=side?>, normalBasisForSide<?=side?>, 1)
<? end ?>
