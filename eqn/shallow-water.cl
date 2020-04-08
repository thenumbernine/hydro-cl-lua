<? local solver = eqn.solver ?>

typedef <?=eqn.prim_t?> prim_t;
typedef <?=eqn.cons_t?> cons_t;
typedef <?=eqn.eigen_t?> eigen_t;
typedef <?=eqn.waves_t?> waves_t;
typedef <?=solver.solver_t?> solver_t;

<? for side=0,solver.dim-1 do ?>
cons_t fluxFromCons_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
) {
	prim_t W = primFromCons(solver, U, x);
	real vj = W.v.s<?=side?>;	//v^j
	
	cons_t F;
	
	F.h = U.m.s<?=side?>;
	F.m = real3_real_mul(U.m, vj);	//h v^i v^j

<? for i=0,2 do
?>	F.m.s<?=i?> += coord_g_uu<?=i?><?=side?>(x) * .5 * solver->gravity * U.h * U.h;
<? end
?>	
	
	return F;
}
<? end ?>

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
	real3 n = normalForSide<?=side?>();
	real3 n2, n3;
	getPerpendicularBasis(n, &n2, &n3);
	
	const real nx = n.x, ny = n.y, nz = n.z;
	const real n2x = n2.x, n2y = n2.y, n2z = n2.z;
	const real n3x = n3.x, n3y = n3.y, n3z = n3.z;
	real v_n = real3_dot(v, n);
	real v_n2 = real3_dot(v, n2);
	real v_n3 = real3_dot(v, n3);

	real3 nU = coord_raise(n, x);
	real nLenSq = real3_dot(n, nU);
	real nLen = sqrt(nLenSq);

	real h = eig.h;
	real3 v = eig.v;
	real C = eig.C;
]]
?>

waves_t eigen_leftTransformForNormal(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x,
	real3 n	//covariant
) { 
	<?=prefix?>

	real3 n2 = _real3(n2x, n2y, n2z);
	real3 n2U = coord_raise(n2, x);
	
	real3 n3 = _real3(n3x, n3y, n3z);
	real3 n3U = coord_raise(n3, x);

	//TODO double-check that this is the correct 3D generalization...
	real vU_dot_n2U = real3_dot(v, n2U);
	real vU_dot_n3U = real3_dot(v, n3U);

	real nlen = sqrt_gUjj;	//sqrt(n^i n_i)
	real nsq = gU.s<?=side..side?>;
	real C_n = C * nlen;
	real denom = 2. * eig.h * nsq;
	real invDenom = 1. / denom;

	return (waves_t){.ptr={
		(X.ptr[0] * (-C_n - v_n) + X.ptr[<?=side+1?>]) * invDenom,
		(X.ptr[0] * -vU_dot_n2U + X.ptr[1] * n2U.x + X.ptr[2] * n2U.y + X.ptr[3] * n2U.z) * 2. * invDenom,
		(X.ptr[0] * -vU_dot_n3U + X.ptr[1] * n3U.x + X.ptr[2] * n3U.y + X.ptr[3] * n3U.z) * 2. * invDenom,
		(X.ptr[0] * (C_n - v_n) + X.ptr[<?=side+1?>]) * invDenom,
	}};
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x
) {
	<?=prefix?>

	real3 n = _real3(nx, ny, nz);
	real3 nU = gUj;//coord_raise(n, x);
	real nlen = sqrt_gUjj;	//sqrt(n^i n_i)
	
	real invC = 1. / C;

	return (cons_t){.ptr={
		nlen * h * invC * (X.ptr[3] - X.ptr[0]),
		
		(X.ptr[0] + X.ptr[3]) * h * nU.x 
		+ (X.ptr[3] - X.ptr[0]) * (nlen * h * invC * v.x)
		+ h * (X.ptr[1] * n2x + X.ptr[2] * n3x),
		
		(X.ptr[0] + X.ptr[3]) * h * nU.y
		+ (X.ptr[3] - X.ptr[0]) * (nlen * h * invC * v.y)
		+ h * (X.ptr[1] * n2y + X.ptr[2] * n3y),

		(X.ptr[0] + X.ptr[3]) * h * nU.z
		+ (X.ptr[3] - X.ptr[0]) * (nlen * h * invC * v.z)
		+ h * (X.ptr[1] * n2z + X.ptr[2] * n3z),
	}};
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	<?=prefix?>
	real3 nU = gUj;
	return (cons_t){.ptr={
		X.ptr[1] * nx 
		+ X.ptr[2] * ny 
		+ X.ptr[3] * nz,
		
		X.ptr[0] * (C * C * nU.x - v.x * v_n)
		+ X.ptr[1] * (v.x * nx + v_n)
		+ X.ptr[2] * (v.x * ny)
		+ X.ptr[3] * (v.x * nz),
		
		X.ptr[0] * (C * C * nU.y - v.y * v_n)
		+ X.ptr[1] * (v.y * nx)
		+ X.ptr[2] * (v.y * ny + v_n)
		+ X.ptr[3] * (v.y * nz),
		
		X.ptr[0] * (C * C * nU.z - v.z * v_n)
		+ X.ptr[1] * (v.z * nx)
		+ X.ptr[2] * (v.z * ny)
		+ X.ptr[3] * (v.z * nz + v_n),
	}};
}
<? end ?>


// used by PLM


<? for side=0,solver.dim-1 do ?>
eigen_t eigen_forCell_<?=side?>(
	constant solver_t* solver,
	cons_t U,
	real3 x
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
<? end ?>
