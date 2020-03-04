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
	real vj = W.v.s<?=side?>;
	
	cons_t F;
	
	F.h = U.m.s<?=side?>;
	F.m = real3_real_mul(U.m, vj);

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
	real C = calc_C(solver, W);
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
	real gUjj = gU.s]]..side..side..[[;
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
	
	real denom = 2. * eig.h;
	real invDenom = 1. / denom;

	//TODO the orthogonal normals should be raised by the metric 
<? if side == 0 then ?>
	real sqrt_gUxx = sqrt_gUjj;
	return (waves_t){.ptr={
		(X.ptr[0] * (-C - v.x)
			+ X.ptr[1]
		) * invDenom,
		(X.ptr[0] * (-2. * v.x)
			+ X.ptr[2] / sqrt_gUxx
		) * invDenom,
		(X.ptr[0] * (-2. * v.x)
			+ X.ptr[3] / sqrt_gUxx
		) * invDenom,
		(X.ptr[0] * (C - v.x)
			+ X.ptr[1]
		) * invDenom,
	}};
<? elseif side == 1 then ?>
	real sqrt_gUyy = sqrt_gUjj;
	return (waves_t){.ptr={
		(X.ptr[0] * (-C - v.y)
			+ X.ptr[2]
		) * invDenom,
		(X.ptr[0] * (-2. * v.y)
			+ X.ptr[3] / sqrt_gUyy
		) * invDenom,
		(X.ptr[0] * (-2. * v.y)
			+ X.ptr[1] / sqrt_gUyy
		) * invDenom,
		(X.ptr[0] * (C - v.y)
			+ X.ptr[2]
		) * invDenom,
	}};
<? elseif side == 2 then ?>
	real sqrt_gUzz = sqrt_gUjj;
	return (waves_t){.ptr={
		(X.ptr[0] * (-C - v.z)
			+ X.ptr[3]
		) * invDenom,
		(X.ptr[0] * (-2. * v.z)
			+ X.ptr[1] / sqrt_gUzz
		) * invDenom,
		(X.ptr[0] * (-2. * v.z)
			+ X.ptr[2] / sqrt_gUzz
		) * invDenom,
		(X.ptr[0] * (C - v.z)
			+ X.ptr[3]
		) * invDenom,
	}};
<? end ?>
}

cons_t eigen_rightTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	waves_t X,
	real3 x
) {
	<?=prefix?>
	real invC = 1. / C;
<? if side == 0 then ?>	
	real sqrt_gUxx = sqrt_gUjj;
	return (cons_t){.ptr={
		h * invC * (X.ptr[3] - X.ptr[0]),
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.xx / sqrt_gUxx - invC * v.x),
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.xy / sqrt_gUxx - invC * v.y)
			+ X.ptr[1] * h,	
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.xz / sqrt_gUxx - invC * v.z)
			+ X.ptr[2] * h,
	}};
<? elseif side == 1 then ?>	
	real sqrt_gUyy = sqrt_gUjj;
	return (cons_t){.ptr={
		h * invC * (X.ptr[3] - X.ptr[0]),
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.xy / sqrt_gUyy - invC * v.x)
			+ X.ptr[2] * h,	
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.yy / sqrt_gUyy - invC * v.y),
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.yz / sqrt_gUyy - invC * v.z)
			+ X.ptr[1] * h,
	}};
<? elseif side == 2 then ?>	
	real sqrt_gUzz = sqrt_gUjj;
	return (cons_t){.ptr={
		h * invC * (X.ptr[3] - X.ptr[0]),
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.xz / sqrt_gUzz - invC * v.x)
			+ X.ptr[1] * h,	
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.yz / sqrt_gUzz - invC * v.y)
			+ X.ptr[2] * h,
		
		(X.ptr[0] + X.ptr[3]) * h * (gU.zz / sqrt_gUzz - invC * v.z),
	}};
<? end ?>
}

cons_t eigen_fluxTransform_<?=side?>(
	constant solver_t* solver,
	eigen_t eig,
	cons_t X,
	real3 x
) {
	<?=prefix?>
	return (cons_t){.ptr={
		X.ptr[1] * nx 
			+ X.ptr[2] * ny 
			+ X.ptr[3] * nz,
		X.ptr[0] * (C * C * gUj.x - v.x * v_n)
			+ X.ptr[1] * (v.x * nx + v_n)
			+ X.ptr[2] * (v.x * ny)
			+ X.ptr[3] * (v.x * nz),
		X.ptr[0] * (C * C * gUj.y - v.y * v_n)
			+ X.ptr[1] * (v.y * nx)
			+ X.ptr[2] * (v.y * ny + v_n)
			+ X.ptr[3] * (v.y * nz),
		X.ptr[0] * (C * C * gUj.z - v.z * v_n)
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
