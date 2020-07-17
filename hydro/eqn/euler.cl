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

cons_t fluxFromCons(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	normalInfo_t n
) {
	prim_t W = primFromCons(solver, U, x);
	real v_n = normalInfo_vecDotN1(n, W.v);
	real HTotal = U.ETotal + W.P;
	
	return (cons_t){
		.rho = U.rho * v_n,
		.m = real3_add(
			real3_real_mul(U.m, v_n),
			_real3(
				normalInfo_u1x(n) * W.P,
				normalInfo_u1y(n) * W.P,
				normalInfo_u1z(n) * W.P
			)
		),
		.ETotal = HTotal * v_n,
		.ePot = 0,
	};
}


// used by PLM
range_t calcCellMinMaxEigenvalues(
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
	normalInfo_t n
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


	return (eigen_t){
		.rho = rho, 
		.v = v,
		.vSq = vSq,
		.vL = vLower,
		.hTotal = hTotal,
		.Cs = Cs,
	};
}

// used by PLM
eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	normalInfo_t n
) {
	prim_t W = primFromCons(solver, U, x);
	real3 vL = coord_lower(W.v, x);
	real vSq = real3_dot(W.v, vL);
	real v_n = normalInfo_vecDotN1(n, W.v);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
	real CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	return (eigen_t){
		.rho = W.rho,
		.v = W.v,
		.vSq = vSq,
		.vL = vL,
		.hTotal = hTotal,
		.Cs = Cs,
	};
}



<?
local prefix = [[
	real3 v = eig.v;
	real3 vL = eig.vL;
	real vSq = eig.vSq;
	real3 v_n = normalInfo_vecDotNs(n, v);
	real hTotal = eig.hTotal;
	real Cs = eig.Cs;
	real nLen = normalInfo_len(n);
]]
?>

waves_t eigen_leftTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	normalInfo_t n
) { 
	real* X = X_.ptr;
	<?=prefix?>
	
	real denom = 2. * Cs * Cs;
	real invDenom = 1. / denom;

	const real gamma_1 = solver->heatCapacityRatio - 1.;
	return (waves_t){.ptr={
		(
			X[0] * (.5 * gamma_1 * vSq + Cs * v_n.x / nLen)
			+ X[1] * (-gamma_1 * vL.x - Cs * normalInfo_l1x_over_len(n))
			+ X[2] * (-gamma_1 * vL.y - Cs * normalInfo_l1y_over_len(n))
			+ X[3] * (-gamma_1 * vL.z - Cs * normalInfo_l1z_over_len(n))
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
		+ X[1] * normalInfo_l2x(n)
		+ X[2] * normalInfo_l2y(n)
		+ X[3] * normalInfo_l2z(n),
		
		X[0] * -v_n.z
		+ X[1] * normalInfo_l3x(n)
		+ X[2] * normalInfo_l3y(n)
		+ X[3] * normalInfo_l3z(n),
		
		(
			X[0] * (.5 * gamma_1 * vSq - Cs * v_n.x / nLen)
			+ X[1] * (-gamma_1 * vL.x + Cs * normalInfo_l1x_over_len(n))
			+ X[2] * (-gamma_1 * vL.y + Cs * normalInfo_l1y_over_len(n))
			+ X[3] * (-gamma_1 * vL.z + Cs * normalInfo_l1z_over_len(n))
			+ X[4] * gamma_1
		) * invDenom,
	}};
}

cons_t eigen_rightTransform(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x,
	normalInfo_t n
) {
	real* X = X_.ptr;
	<?=prefix?>
	return (cons_t){.ptr={
		X[0] 
		+ X[1] 
		+ X[4],
		
		X[0] * (v.x - Cs * normalInfo_u1x_over_len(n))
		+ X[1] * v.x
		+ X[2] * normalInfo_u2x(n)
		+ X[3] * normalInfo_u3x(n)
		+ X[4] * (v.x + Cs * normalInfo_u1x_over_len(n)),
		
		X[0] * (v.y - Cs * normalInfo_u1y_over_len(n))
		+ X[1] * v.y
		+ X[2] * normalInfo_u2y(n)
		+ X[3] * normalInfo_u3y(n)
		+ X[4] * (v.y + Cs * normalInfo_u1y_over_len(n)),
		
		X[0] * (v.z - Cs * normalInfo_u1z_over_len(n))
		+ X[1] * v.z
		+ X[2] * normalInfo_u2z(n)
		+ X[3] * normalInfo_u3z(n)
		+ X[4] * (v.z + Cs * normalInfo_u1z_over_len(n)),
		
		X[0] * (hTotal - Cs * v_n.x / nLen)
		+ X[1] * .5 * vSq
		+ X[2] * v_n.y
		+ X[3] * v_n.z
		+ X[4] * (hTotal + Cs * v_n.x / nLen),
		
		0,
	}};
}

cons_t eigen_fluxTransform(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	normalInfo_t n
) {
	real* X = X_.ptr;
	<?=prefix?>
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	return (cons_t){.ptr={
		X[1] * normalInfo_l1x(n)
		+ X[2] * normalInfo_l1y(n)
		+ X[3] * normalInfo_l1z(n),
		
		X[0] * (-v_n.x * v.x + gamma_1 * .5 * vSq * normalInfo_u1x(n))
		+ X[1] * (v.x * normalInfo_l1x(n) - gamma_2 * normalInfo_u1x(n) * vL.x + v_n.x)
		+ X[2] * (v.x * normalInfo_l1y(n) - gamma_2 * normalInfo_u1x(n) * vL.y)
		+ X[3] * (v.x * normalInfo_l1z(n) - gamma_2 * normalInfo_u1x(n) * vL.z)
		+ X[4] * gamma_1 * normalInfo_u1x(n),
		
		X[0] * (-v_n.x * v.y + gamma_1 * .5 * vSq * normalInfo_u1y(n))
		+ X[1] * (v.y * normalInfo_l1x(n) - gamma_2 * normalInfo_u1y(n) * vL.x)
		+ X[2] * (v.y * normalInfo_l1y(n) - gamma_2 * normalInfo_u1y(n) * vL.y + v_n.x)
		+ X[3] * (v.y * normalInfo_l1z(n) - gamma_2 * normalInfo_u1y(n) * vL.z)
		+ X[4] * gamma_1 * normalInfo_u1y(n),
		
		X[0] * (-v_n.x * v.z + gamma_1 * .5 * vSq * normalInfo_u1z(n))
		+ X[1] * (v.z * normalInfo_l1x(n) - gamma_2 * normalInfo_u1z(n) * vL.x)
		+ X[2] * (v.z * normalInfo_l1y(n) - gamma_2 * normalInfo_u1z(n) * vL.y)
		+ X[3] * (v.z * normalInfo_l1z(n) - gamma_2 * normalInfo_u1z(n) * vL.z + v_n.x)
		+ X[4] * gamma_1 * normalInfo_u1z(n),
		
		X[0] * v_n.x * (.5 * gamma_1 * vSq - hTotal)
		+ X[1] * (normalInfo_l1x(n) * hTotal - gamma_1 * v_n.x * vL.x)
		+ X[2] * (normalInfo_l1y(n) * hTotal - gamma_1 * v_n.x * vL.y)
		+ X[3] * (normalInfo_l1z(n) * hTotal - gamma_1 * v_n.x * vL.z)
		+ X[4] * gamma * v_n.x,
		
		0,
	}};
}


kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cellBuf[index].pos;

	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

<? if false 
and solver.coord.vectorComponent == 'anholonomic' 
and require 'hydro.coord.cylinder'.is(solver.coord) 
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
<? if not require 'hydro.coord.cartesian'.is(solver.coord) then ?>
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
	global cons_t* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;

	global cons_t* U = UBuf + index;
	prim_t W = primFromCons(solver, *U, x);

	if (W.rho < solver->rhoMin) W.rho = solver->rhoMin;
	if (W.P < solver->PMin) W.P = solver->PMin;

	*U = consFromPrim(solver, W, x);
}
