<? local solver = eqn.solver ?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

cons_t fluxFromConsForNormal(
	constant solver_t* solver,
	cons_t U,
	real3 x,
	real3 n	//covariant components
) {
	prim_t W = primFromCons(solver, U, x);
	real3 nU = coord_raise(n, x);
	real v_n = real3_dot(W.v, n);
	return (cons_t){
		.h = U.h * v_n,
		.m = real3_add(
			real3_real_mul(U.m, v_n),	//h v^i v_n
			real3_real_mul(nU, .5 * solver->gravity * U.h * U.h)	//.5 g h^2 n^i
		),
	};
}

// used by PLM
<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	constant solver_t* solver,
	const global cons_t* U,
	real3 x
) {
	prim_t W = primFromCons(solver, *U, x);
	real C = calc_C(solver, *U);
	//for n_i = partial_i, n^i = g^ij n_j, |n| = sqrt(n^i n_i) ... when n_i is aligned to the coordinate basis then |n| = sqrt(g^ii)
	real nLen = coord_sqrt_g_uu<?=side..side?>(x);
	real C_nLen = C * nLen; 
	return (range_t){
		.min = W.v.s<?=side?> - C_nLen, 
		.max = W.v.s<?=side?> + C_nLen,
	};
}
<? end ?>

//used by the mesh version
//I'm winging this from Euler equations
eigen_t eigen_forInterface(
	constant solver_t* solver,
	cons_t UL,
	cons_t UR,
	real3 x,
	real3 n
) {
	prim_t WL = primFromCons(solver, UL, x);
	real sqrtHL = sqrt(WL.h);
	real3 vL = WL.v;
	
	prim_t WR = primFromCons(solver, UR, x);
	real sqrtHR = sqrt(WR.h);
	real3 vR = WR.v;

	real invDenom = 1./(sqrtHL + sqrtHR);
	
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
	real3 n2, n3;
	getPerpendicularBasis(n, &n2, &n3);

	real3 nU = coord_raise(n, x);
	real nLenSq = real3_dot(n, nU);
	real nLen = sqrt(nLenSq);

	real3 v = eig.v;
	real h = eig.h;
	real C = eig.C;

	real v_n = real3_dot(v, n);
	real v_n2 = real3_dot(v, n2);
	real v_n3 = real3_dot(v, n3);
]]
?>

waves_t eigen_leftTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	real3 n	//covariant
) { 
	real* X = X_.ptr;
	<?=prefix?>

	real3 n2U = coord_raise(n2, x);
	real3 n3U = coord_raise(n3, x);

	//TODO double-check that this is the correct 3D generalization...
	real vU_dot_n2U = real3_dot(v, n2U);
	real vU_dot_n3U = real3_dot(v, n3U);

	real C_nLen = C * nLen;
	real denom = 2. * eig.h * nLenSq;
	real invDenom = 1. / denom;

	return (waves_t){.ptr={
		(X[0] * (-C_nLen - v_n) + X[1] * n.x + X[2] * n.y + X[3] * n.z) * invDenom,
		(X[0] * -vU_dot_n2U + X[1] * n2U.x + X[2] * n2U.y + X[3] * n2U.z) * 2. * invDenom,
		(X[0] * -vU_dot_n3U + X[1] * n3U.x + X[2] * n3U.y + X[3] * n3U.z) * 2. * invDenom,
		(X[0] * (C_nLen - v_n) + X[1] * n.x + X[2] * n.y + X[3] * n.z) * invDenom,
	}};
}

cons_t eigen_rightTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X_,
	real3 x,
	real3 n	//covariant
) {
	real* X = X_.ptr;
	<?=prefix?>

	real invC = 1. / C;

	return (cons_t){.ptr={
		nLen * h * invC * (X[3] - X[0]),
		
		(X[0] + X[3]) * h * nU.x 
		+ (X[3] - X[0]) * (nLen * h * invC * v.x)
		+ h * (X[1] * n2.x + X[2] * n3.x),
		
		(X[0] + X[3]) * h * nU.y
		+ (X[3] - X[0]) * (nLen * h * invC * v.y)
		+ h * (X[1] * n2.y + X[2] * n3.y),

		(X[0] + X[3]) * h * nU.z
		+ (X[3] - X[0]) * (nLen * h * invC * v.z)
		+ h * (X[1] * n2.z + X[2] * n3.z),
	}};
}

cons_t eigen_fluxTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X_,
	real3 x,
	real3 n
) {
	real* X = X_.ptr;
	<?=prefix?>
	return (cons_t){.ptr={
		X[1] * n.x 
		+ X[2] * n.y 
		+ X[3] * n.z,
		
		X[0] * (C * C * nU.x - v.x * v_n)
		+ X[1] * (v.x * n.x + v_n)
		+ X[2] * (v.x * n.y)
		+ X[3] * (v.x * n.z),
		
		X[0] * (C * C * nU.y - v.y * v_n)
		+ X[1] * (v.y * n.x)
		+ X[2] * (v.y * n.y + v_n)
		+ X[3] * (v.y * n.z),
		
		X[0] * (C * C * nU.z - v.z * v_n)
		+ X[1] * (v.z * n.x)
		+ X[2] * (v.z * n.y)
		+ X[3] * (v.z * n.z + v_n),
	}};
}


// used by PLM


eigen_t eigen_forCellForNormal(
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


<? for side=0,solver.dim-1 do ?>
#define fluxFromCons_<?=side?>(solver, U, x) 				fluxFromConsForNormal(solver, U, x, normalForSide<?=side?>())
#define calcCellMinMaxEigenvalues_<?=side?>(solver, U, x) 	calcCellMinMaxEigenvaluesForNormal(solver, U, x, normalForSide<?=side?>())
#define eigen_leftTransform_<?=side?>(solver, eig, X, x) 	eigen_leftTransformForNormal(solver, eig, X, x, normalForSide<?=side?>())
#define eigen_rightTransform_<?=side?>(solver, eig, X, x) 	eigen_rightTransformForNormal(solver, eig, X, x, normalForSide<?=side?>())
#define eigen_fluxTransform_<?=side?>(solver, eig, X, x) 	eigen_fluxTransformForNormal(solver, eig, X, x, normalForSide<?=side?>())
#define eigen_forCell_<?=side?>(solver, U, x) 				eigen_forCellForNormal(solver, U, x, normalForSide<?=side?>())
<? end ?>
