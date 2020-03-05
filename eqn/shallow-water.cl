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
	real C_sqrt_gU = C * coord_sqrt_g_uu<?=side..side?>(x);
	return (range_t){
		.min = W.v.s<?=side?> - C_sqrt_gU, 
		.max = W.v.s<?=side?> + C_sqrt_gU,
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
for side=0,solver.dim-1 do 
	local prefix
	
	if side == 0 then
					
		prefix = [[
	const real nx = 1, ny = 0, nz = 0;
	const real n1x = 0, n1y = 1, n1z = 0;
	const real n2x = 0, n2y = 0, n2z = 1;
	real v_n = v.x, v_n1 = v.y, v_n2 = v.z;
]] 
	elseif side == 1 then
		
		prefix = [[
	const real nx = 0, ny = 1, nz = 0;
	const real n1x = 0, n1y = 0, n1z = 1;
	const real n2x = 1, n2y = 0, n2z = 0;
	real v_n = v.y, v_n1 = v.z, v_n2 = v.x;
]] 
	elseif side == 2 then
		
		prefix = [[
	const real nx = 0, ny = 0, nz = 1;
	const real n1x = 1, n1y = 0, n1z = 0;
	const real n2x = 0, n2y = 1, n2z = 0;
	real v_n = v.z, v_n1 = v.x, v_n2 = v.y;
]]
	end
	
	prefix = [[
	sym3 gU = coord_g_uu(x);
	real sqrt_gUjj = coord_sqrt_g_uu]]..side..side..[[(x);
	
	real h = eig.h;
	real3 v = eig.v;
	real C = eig.C;
	//g^ij for fixed j=side
]] .. prefix

	local gUdef = '\treal3 gUj = _real3(\n'
	for i=0,2 do
		gUdef = gUdef .. '\t\tcoord_g_uu'..side..i..'(x)'..(i<2 and ',' or '')..'\n'
	end
	gUdef = gUdef .. '\t);\n'
	prefix = gUdef .. prefix
?>

waves_t eigen_leftTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) { 
	<?=prefix?>

	real3 n1 = _real3(n1x, n1y, n1z);
	real3 n1U = coord_raise(n1, x);
	
	real3 n2 = _real3(n2x, n2y, n2z);
	real3 n2U = coord_raise(n2, x);

	//TODO double-check that this is the correct 3D generalization...
	real vU_dot_n1U = real3_dot(v, n1U);
	real vU_dot_n2U = real3_dot(v, n2U);

	real nlen = sqrt_gUjj;	//sqrt(n^i n_i)
	real nsq = gU.s<?=side..side?>;
	real C_n = C * nlen;
	real denom = 2. * eig.h * nsq;
	real invDenom = 1. / denom;

	return (waves_t){.ptr={
		(X.ptr[0] * (-C_n - v_n) + X.ptr[<?=side+1?>]) * invDenom,
		(X.ptr[0] * -vU_dot_n1U + X.ptr[1] * n1U.x + X.ptr[2] * n1U.y + X.ptr[3] * n1U.z) * 2. * invDenom,
		(X.ptr[0] * -vU_dot_n2U + X.ptr[1] * n2U.x + X.ptr[2] * n2U.y + X.ptr[3] * n2U.z) * 2. * invDenom,
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
		+ h * (X.ptr[1] * n1x + X.ptr[2] * n2x),
		
		(X.ptr[0] + X.ptr[3]) * h * nU.y
		+ (X.ptr[3] - X.ptr[0]) * (nlen * h * invC * v.y)
		+ h * (X.ptr[1] * n1y + X.ptr[2] * n2y),

		(X.ptr[0] + X.ptr[3]) * h * nU.z
		+ (X.ptr[3] - X.ptr[0]) * (nlen * h * invC * v.z)
		+ h * (X.ptr[1] * n1z + X.ptr[2] * n2z),
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
