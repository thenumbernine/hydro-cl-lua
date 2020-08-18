typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

// used by PLM
range_t calcCellMinMaxEigenvalues(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x,
	normalInfo_t n
) {
	prim_t W = primFromCons(solver, *U, x);
	real v_n = normalInfo_vecDotN1(n, W.v);
	real C = calc_C(solver, *U);
	real C_nLen = C * normalInfo_len(n);
	return (range_t){
		.min = v_n - C_nLen, 
		.max = v_n + C_nLen,
	};
}

eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	normalInfo_t n
) {
	prim_t WL = primFromCons(solver, UL, x);
	real sqrtHL = sqrt(WL.h);
	real3 vL = WL.v;
	
	prim_t WR = primFromCons(solver, UR, x);
	real sqrtHR = sqrt(WR.h);
	real3 vR = WR.v;

	real invDenom = 1. / (sqrtHL + sqrtHR);
	
	//Roe-averaged
	real h = sqrtHL * sqrtHR;
	real3 v = real3_add(
			real3_real_mul(vL, sqrtHL * invDenom),
			real3_real_mul(vR, sqrtHR * invDenom));

	//derived:
	real CSq = solver->gravity * h;
	real C = sqrt(CSq);

	return (eigen_t){
		.h = h, 
		.v = v,
		.C = C,
	};
}

<?
local prefix = [[
	real3 nL = normalInfo_l1(n);
	real3 nU = normalInfo_u1(n);
	real nLen = normalInfo_len(n);
	real nLenSq = nLen * nLen;

	real3 v = eig.v;
	real h = eig.h;
	real C = eig.C;

	real3 v_ns = normalInfo_vecDotNs(n, v);
	real v_n = v_ns.x;
	real v_n2 = v_ns.y;
	real v_n3 = v_ns.z;
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

	real3 n2U = normalInfo_u2(n);
	real3 n3U = normalInfo_u3(n);

	//TODO double-check that this is the correct 3D generalization...
	// seems strange to have an upper dot an upper
	real vU_dot_n2U = real3_dot(v, n2U);
	real vU_dot_n3U = real3_dot(v, n3U);

	real C_nLen = C * nLen;
	real denom = 2. * eig.h * nLenSq;
	real invDenom = 1. / denom;

	return (waves_t){.ptr={
		(X[0] * (-C_nLen - v_n) + X[1] * nL.x + X[2] * nL.y + X[3] * nL.z) * invDenom,
		(X[0] * -vU_dot_n2U + X[1] * n2U.x + X[2] * n2U.y + X[3] * n2U.z) * 2. * invDenom,
		(X[0] * -vU_dot_n3U + X[1] * n3U.x + X[2] * n3U.y + X[3] * n3U.z) * 2. * invDenom,
		(X[0] * (C_nLen - v_n) + X[1] * nL.x + X[2] * nL.y + X[3] * nL.z) * invDenom,
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

	real invC = 1. / C;

	return (cons_t){.ptr={
		nLen * h * invC * (X[3] - X[0]),
		
		(X[0] + X[3]) * h * nU.x 
		+ (X[3] - X[0]) * (nLen * h * invC * v.x)
		+ h * (X[1] * normalInfo_l2x(n) + X[2] * normalInfo_l3x(n)),
		
		(X[0] + X[3]) * h * nU.y
		+ (X[3] - X[0]) * (nLen * h * invC * v.y)
		+ h * (X[1] * normalInfo_l2y(n) + X[2] * normalInfo_l3y(n)),

		(X[0] + X[3]) * h * nU.z
		+ (X[3] - X[0]) * (nLen * h * invC * v.z)
		+ h * (X[1] * normalInfo_l2z(n) + X[2] * normalInfo_l3z(n)),
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
	return (cons_t){.ptr={
		X[1] * nL.x 
		+ X[2] * nL.y 
		+ X[3] * nL.z,
		
		X[0] * (C * C * nU.x - v.x * v_n)
		+ X[1] * (v.x * nL.x + v_n)
		+ X[2] * (v.x * nL.y)
		+ X[3] * (v.x * nL.z),
		
		X[0] * (C * C * nU.y - v.y * v_n)
		+ X[1] * (v.y * nL.x)
		+ X[2] * (v.y * nL.y + v_n)
		+ X[3] * (v.y * nL.z),
		
		X[0] * (C * C * nU.z - v.z * v_n)
		+ X[1] * (v.z * nL.x)
		+ X[2] * (v.z * nL.y)
		+ X[3] * (v.z * nL.z + v_n),
	}};
}


// used by PLM
eigen_t eigen_forCell(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	real3 n
) {
	prim_t W = primFromCons(solver, U, x);
	real CSq = solver->gravity * U.h;
	real C = sqrt(CSq);
	return (eigen_t){
		.h = W.h,
		.v = W.v,
		.C = C,
	};
}
