/*
I've highjacked all of this.  It was a normal Euler eqn solver.
But I experimented with a curved-space solver.  
To get back to the original code,
just replace all the g_ab stuff with their constant values and simplify away.
*/
<? local solver = eqn.solver ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real vj = W.v.s<?=side?>;
	real HTotal = U.ETotal + W.P;
	
	<?=eqn.cons_t?> F;
	
	F.rho = U.m.s<?=side?>;
	
	F.m = real3_real_mul(U.m, vj);

<? for i=0,2 do
?>	F.m.s<?=i?> += coord_gU<?=i?><?=side?>(x) * W.P;
<? end
?>	F.ETotal = HTotal * vj;
	
	F.ePot = 0;
	
	return F;
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real Cs = calc_Cs(&W);
	real Cs_sqrt_gU = Cs * coord_sqrt_gU<?=side..side?>(x);
	return (range_t){
		.min = W.v.s<?=side?> - Cs_sqrt_gU, 
		.max = W.v.s<?=side?> + Cs_sqrt_gU,
	};
}
<? end ?>

//used by the mesh version
<?=eqn.eigen_t?> eigen_forInterface(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	real3 n
) {
	<?=eqn.prim_t?> WL = primFromCons(UL, x);
	real sqrtRhoL = sqrt(WL.rho);
	real3 vL = WL.v;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL.ETotal) - UL.ePot;
	
	<?=eqn.prim_t?> WR = primFromCons(UR, x);
	real sqrtRhoR = sqrt(WR.rho);
	real3 vR = WR.v;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR.ETotal) - UR.ePot;

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);
	
	//Roe-averaged
	real rho = sqrtRhoL * sqrtRhoR;
	real3 v = real3_add(
			real3_real_mul(vL, sqrtRhoL * invDenom),
			real3_real_mul(vR, sqrtRhoR * invDenom));
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);

	//derived:
	real vSq = coordLenSq(v, x);
	real eKin = .5 * vSq;
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);

	return (<?=eqn.eigen_t?>){
		.rho = rho, 
		.v = v,
		.hTotal = hTotal,
		.vSq = vSq,
		.Cs = Cs,
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
	sym3 gU = coord_gU(x);
	real gUjj = gU.s]]..side..side..[[;
	real sqrt_gUjj = coord_sqrt_gU]]..side..side..[[(x);
	
	real3 v = eig.v;
	real3 vL = coord_lower(v, x);
	real hTotal = eig.hTotal;
	real vSq = real3_dot(v, vL);
	real Cs = eig.Cs;
	real Cs_over_sqrt_gUjj = Cs / sqrt_gUjj; 
	//g^ij for fixed j=side
]] .. prefix

	local gUdef = '\treal3 gUj = _real3(\n'
	for i=0,2 do
		gUdef = gUdef .. '\t\tcoord_gU'..side..i..'(x)'..(i<2 and ',' or '')..'\n'
	end
	gUdef = gUdef .. '\t);\n'
	prefix = gUdef .. prefix
?>

<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) { 
	<?=prefix?>
	
	real denom = 2. * Cs * Cs;
	real invDenom = 1. / denom;

#if 0	//works but isn't correct for curvilinear coordinates (which make use of g_ij)
	Y[0] = (X.ptr[0] * ((heatCapacityRatio - 1.) * .5 * vSq + Cs_over_sqrt_gUjj * v_n)
		+ X.ptr[1] * -(nx * Cs_over_sqrt_gUjj + (heatCapacityRatio - 1.) * vL.x) 
		+ X.ptr[2] * -(ny * Cs_over_sqrt_gUjj + (heatCapacityRatio - 1.) * vL.y)
		+ X.ptr[3] * -(nz * Cs_over_sqrt_gUjj + (heatCapacityRatio - 1.) * vL.z)
		+ X.ptr[4] * (heatCapacityRatio - 1.)
	) * invDenom;
	Y[1] = (X.ptr[0] * (2.*Cs*Cs - (heatCapacityRatio - 1.) * vSq)
		+ X.ptr[1] * (heatCapacityRatio - 1.) * v.x * 2
		+ X.ptr[2] * (heatCapacityRatio - 1.) * v.y * 2
		+ X.ptr[3] * (heatCapacityRatio - 1.) * v.z * 2
		+ X.ptr[4] * -(heatCapacityRatio - 1.) * 2
	) * invDenom;
	Y[2] = X.ptr[0] * -v_n1 + X.ptr[1] * n1x + X.ptr[2] * n1y + X.ptr[3] * n1z;
	Y[3] = X.ptr[0] * -v_n2 + X.ptr[1] * n2x + X.ptr[2] * n2y + X.ptr[3] * n2z;
	Y[4] = (X.ptr[0] * ((heatCapacityRatio - 1.) * .5 * vSq - Cs_over_sqrt_gUjj * v_n) 
		+ X.ptr[1] * (nx * Cs_over_sqrt_gUjj - (heatCapacityRatio - 1.) * vL.x) 
		+ X.ptr[2] * (ny * Cs_over_sqrt_gUjj - (heatCapacityRatio - 1.) * vL.y) 
		+ X.ptr[3] * (nz * Cs_over_sqrt_gUjj - (heatCapacityRatio - 1.) * vL.z) 
		+ X.ptr[4] * (heatCapacityRatio - 1.)
	) * invDenom;
#else
	const real heatRatioMinusOne = heatCapacityRatio - 1.;
<? if side == 0 then ?>
	real sqrt_gUxx = sqrt_gUjj;
	return (<?=eqn.waves_t?>){.ptr={
		(X.ptr[0] * (.5 * heatRatioMinusOne * vSq + Cs * v.x / sqrt_gUxx)
			+ X.ptr[1] * (-heatRatioMinusOne * vL.x - Cs / sqrt_gUxx)
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
		(X.ptr[0] * (denom - heatRatioMinusOne * vSq)
			+ X.ptr[1] * 2. * heatRatioMinusOne * vL.x
			+ X.ptr[2] * 2. * heatRatioMinusOne * vL.y
			+ X.ptr[3] * 2. * heatRatioMinusOne * vL.z
			+ X.ptr[4] * -2. * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (v.x * gU.xy / gU.xx - v.y)
			+ X.ptr[1] * -gU.xy / gU.xx
			+ X.ptr[2],
		X.ptr[0] * (v.x * gU.xz / gU.xx - v.z)
			+ X.ptr[1] * -gU.xz / gU.xx
			+ X.ptr[3],
		(X.ptr[0] * (.5 * heatRatioMinusOne * vSq - Cs * v.x / sqrt_gUxx)
			+ X.ptr[1] * (-heatRatioMinusOne * vL.x + Cs / sqrt_gUxx)
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
	}};
<? elseif side == 1 then ?>
	real sqrt_gUyy = sqrt_gUjj;
	return (<?=eqn.waves_t?>){.ptr={
		(X.ptr[0] * (.5 * heatRatioMinusOne * vSq + Cs * v.y / sqrt_gUyy)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * (-heatRatioMinusOne * vL.y - Cs / sqrt_gUyy)
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (v.y * gU.xy / gU.yy - v.x)
			+ X.ptr[1]
			+ X.ptr[2] * -gU.xy / gU.yy,
		(X.ptr[0] * (denom - heatRatioMinusOne * vSq)
			+ X.ptr[1] * 2. * heatRatioMinusOne * vL.x
			+ X.ptr[2] * 2. * heatRatioMinusOne * vL.y
			+ X.ptr[3] * 2. * heatRatioMinusOne * vL.z
			+ X.ptr[4] * -2. * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (v.y * gU.yz / gU.yy - v.z)
			+ X.ptr[2] * -gU.yz / gU.yy
			+ X.ptr[3],
		(X.ptr[0] * (.5 * heatRatioMinusOne * vSq - Cs * v.y / sqrt_gUyy)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * (-heatRatioMinusOne * vL.y + Cs / sqrt_gUyy)
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
	}};
<? elseif side == 2 then ?>
	real sqrt_gUzz = sqrt_gUjj;
	return (<?=eqn.waves_t?>){.ptr={
		(X.ptr[0] * (.5 * heatRatioMinusOne * vSq + Cs * v.z / sqrt_gUzz)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * (-heatRatioMinusOne * vL.z - Cs / sqrt_gUzz)
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (v.z * gU.xz / gU.zz - v.x)
			+ X.ptr[1]
			+ X.ptr[3] * -gU.xz / gU.zz,
		X.ptr[0] * (v.z * gU.yz / gU.zz - v.y)
			+ X.ptr[2]
			+ X.ptr[3] * -gU.yz / gU.zz,
		(X.ptr[0] * (denom - heatRatioMinusOne * vSq)
			+ X.ptr[1] * 2. * heatRatioMinusOne * vL.x
			+ X.ptr[2] * 2. * heatRatioMinusOne * vL.y
			+ X.ptr[3] * 2. * heatRatioMinusOne * vL.z
			+ X.ptr[4] * -2. * heatRatioMinusOne
		) * invDenom,
		(X.ptr[0] * (.5 * heatRatioMinusOne * vSq - Cs * v.z / sqrt_gUzz)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * (-heatRatioMinusOne * vL.z + Cs / sqrt_gUzz)
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
	}};
<? end ?>
#endif
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=prefix?>
#if 0	//works but doesn't account for curved space
	Y[0] = X.ptr[0] + X.ptr[1] + X.ptr[4];
	Y[1] = X.ptr[0] * (v.x - gUj.x * Cs_over_sqrt_gUjj) 
		+ X.ptr[1] * v.x 
		+ X.ptr[2] * n1x 
		+ X.ptr[3] * n2x 
		+ X.ptr[4] * (v.x + gUj.x * Cs_over_sqrt_gUjj);
	Y[2] = X.ptr[0] * (v.y - gUj.y * Cs_over_sqrt_gUjj) 
		+ X.ptr[1] * v.y 
		+ X.ptr[2] * n1y 
		+ X.ptr[3] * n2y 
		+ X.ptr[4] * (v.y + gUj.y * Cs_over_sqrt_gUjj);
	Y[3] = X.ptr[0] * (v.z - gUj.z * Cs_over_sqrt_gUjj) 
		+ X.ptr[1] * v.z 
		+ X.ptr[2] * n1z 
		+ X.ptr[3] * n2z 
		+ X.ptr[4] * (v.z + gUj.z * Cs_over_sqrt_gUjj);
	Y[4] = X.ptr[0] * (hTotal - v_n * Cs_over_sqrt_gUjj) 
		+ X.ptr[1] * .5 * vSq
		+ X.ptr[2] * v_n1 
		+ X.ptr[3] * v_n2 
		+ X.ptr[4] * (hTotal + v_n * Cs_over_sqrt_gUjj);
#else	//works for contravariant formulation of Euler fluid equations
<? if side == 0 then ?>	
	real sqrt_gUxx = sqrt_gUjj;
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[0] + X.ptr[1] + X.ptr[4],
		X.ptr[0] * (v.x - Cs * sqrt_gUxx)
			+ X.ptr[1] * v.x
			+ X.ptr[4] * (v.x + Cs * sqrt_gUxx),
		X.ptr[0] * (v.y - Cs * gU.xy / sqrt_gUxx)
			+ X.ptr[1] * v.y
			+ X.ptr[2]
			+ X.ptr[4] * (v.y + Cs * gU.xy / sqrt_gUxx),
		X.ptr[0] * (v.z - Cs * gU.xz / sqrt_gUxx)
			+ X.ptr[1] * v.z
			+ X.ptr[3]
			+ X.ptr[4] * (v.z + Cs * gU.xz / sqrt_gUxx),
		X.ptr[0] * (hTotal - Cs * v.x / sqrt_gUxx)
			+ X.ptr[1] * vSq / 2.
			+ X.ptr[2] * vL.y
			+ X.ptr[3] * vL.z
			+ X.ptr[4] * (hTotal + Cs * v.x / sqrt_gUxx),
		0,
	}};
<? elseif side == 1 then ?>	
	real sqrt_gUyy = sqrt_gUjj;
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[0] + X.ptr[2] + X.ptr[4],
		X.ptr[0] * (v.x - Cs * gU.xy / sqrt_gUyy)
			+ X.ptr[1]
			+ X.ptr[2] * v.x
			+ X.ptr[4] * (v.x + Cs * gU.xy / sqrt_gUyy),
		X.ptr[0] * (v.y - Cs * sqrt_gUyy)
			+ X.ptr[2] * v.y
			+ X.ptr[4] * (v.y + Cs * sqrt_gUyy),
		X.ptr[0] * (v.z - Cs * gU.yz / sqrt_gUyy)
			+ X.ptr[2] * v.z
			+ X.ptr[3]
			+ X.ptr[4] * (v.z + Cs * gU.yz / sqrt_gUyy),
		X.ptr[0] * (hTotal - Cs * v.y / sqrt_gUyy)
			+ X.ptr[1] * vL.x
			+ X.ptr[2] * vSq / 2.
			+ X.ptr[3] * vL.z
			+ X.ptr[4] * (hTotal + Cs * v.y / sqrt_gUyy),
		0,
	}};
<? elseif side == 2 then ?>	
	real sqrt_gUzz = sqrt_gUjj;
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[0] + X.ptr[3] + X.ptr[4],
		X.ptr[0] * (v.x - Cs * gU.xz / sqrt_gUzz)
			+ X.ptr[1]
			+ X.ptr[3] * v.x
			+ X.ptr[4] * (v.x + Cs * gU.xz / sqrt_gUzz),
		X.ptr[0] * (v.y - Cs * gU.yz / sqrt_gUzz)
			+ X.ptr[2]
			+ X.ptr[3] * v.y
			+ X.ptr[4] * (v.y + Cs * gU.yz / sqrt_gUzz),
		X.ptr[0] * (v.z - Cs * sqrt_gUzz)
			+ X.ptr[3] * v.z
			+ X.ptr[4] * (v.z + Cs * sqrt_gUzz),
		X.ptr[0] * (hTotal - Cs * v.z / sqrt_gUzz)
			+ X.ptr[1] * vL.x
			+ X.ptr[2] * vL.y
			+ X.ptr[3] * vSq / 2.
			+ X.ptr[4] * (hTotal + Cs * v.z / sqrt_gUzz),
		0,
	}};
<? end ?>
#endif
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=prefix?>
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[1] * nx 
			+ X.ptr[2] * ny 
			+ X.ptr[3] * nz,
		X.ptr[0] * (-v_n * v.x + (heatCapacityRatio - 1.) * .5 * vSq * gUj.x)
			+ X.ptr[1] * (v.x * nx - (heatCapacityRatio - 1.) * gUj.x * vL.x + v_n)
			+ X.ptr[2] * (v.x * ny - (heatCapacityRatio - 1.) * gUj.x * vL.y)
			+ X.ptr[3] * (v.x * nz - (heatCapacityRatio - 1.) * gUj.x * vL.z)
			+ X.ptr[4] * (heatCapacityRatio - 1.) * nx,
		X.ptr[0] * (-v_n * v.y + (heatCapacityRatio - 1.) * .5 * vSq * gUj.y)
			+ X.ptr[1] * (v.y * nx - (heatCapacityRatio - 1.) * gUj.y * vL.x)
			+ X.ptr[2] * (v.y * ny - (heatCapacityRatio - 1.) * gUj.y * vL.y + v_n)
			+ X.ptr[3] * (v.y * nz - (heatCapacityRatio - 1.) * gUj.y * vL.z)
			+ X.ptr[4] * (heatCapacityRatio - 1.) * ny,
		X.ptr[0] * (-v_n * v.z + (heatCapacityRatio - 1.) * .5 * vSq * gUj.z)
			+ X.ptr[1] * (v.z * nx - (heatCapacityRatio - 1.) * gUj.z * vL.x)
			+ X.ptr[2] * (v.z * ny - (heatCapacityRatio - 1.) * gUj.z * vL.y)
			+ X.ptr[3] * (v.z * nz - (heatCapacityRatio - 1.) * gUj.z * vL.z + v_n)
			+ X.ptr[4] * (heatCapacityRatio - 1.) * nz,
		X.ptr[0] * v_n * ((heatCapacityRatio - 1.) * .5 * vSq - hTotal)
			+ X.ptr[1] * (-(heatCapacityRatio - 1.) * v_n * vL.x + nx * hTotal)
			+ X.ptr[2] * (-(heatCapacityRatio - 1.) * v_n * vL.y + ny * hTotal)
			+ X.ptr[3] * (-(heatCapacityRatio - 1.) * v_n * vL.z + nz * hTotal)
			+ X.ptr[4] * heatCapacityRatio * v_n,
		0,
	}};
}
<? end ?>


// used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real vSq = coordLenSq(W.v, x);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	return (<?=eqn.eigen_t?>){
		.rho = W.rho,
		.v = W.v,
		.hTotal = hTotal,
		.vSq = vSq,
		.Cs = Cs,
	};
}
<? end ?>

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

<? if not require 'coord.cartesian'.is(solver.coord) then ?>
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real3 m_conn_vv = coord_conn_apply23(W.v, U->m, x);
	deriv->m = real3_sub(deriv->m, m_conn_vv);	//-Conn^i_jk rho v^j v^k 
	deriv->m = real3_sub(deriv->m, real3_real_mul(coord_conn_trace23(x), W.P));		//-Conn^i_jk g^jk P
<? end ?>
}
