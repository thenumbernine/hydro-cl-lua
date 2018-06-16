<? local solver = eqn.solver ?>

<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real vj = W.vTilde.s<?=side?>;
	real HTotal = U.rhoBar_eTotalTilde + W.PStar;
	
	<?=eqn.cons_t?> F;
	
	F.rhoBar = U.rhoBar_vTilde.s<?=side?>;
	
	F.rhoBar_vTilde = real3_scale(U.rhoBar_vTilde, vj);

<? for i=0,2 do
?>	F.rhoBar_vTilde.s<?=i?> += coord_gU<?=i?><?=side?>(x) * W.PStar;
<? end
?>	F.rhoBar_eTotalTilde = HTotal * vj;

	F.rhoBar_k = 0;
	F.rhoBar_omega = 0;

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
	real Cs = calc_Cs(W);
	real Cs_sqrt_gU = Cs * coord_sqrt_gU<?=side..side?>(x);
	return (range_t){
		.min = W.vTilde.s<?=side?> - Cs_sqrt_gU, 
		.max = W.vTilde.s<?=side?> + Cs_sqrt_gU,
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
	real sqrtRhoL = sqrt(WL.rhoBar);
	real3 vL = WL.vTilde;
	//real hTotalL = calc_hTotal(WL.rhoBar, WL.PStar, UL.rhoBar_eTotalTilde) - UL.ePot;
	real hTotalL = (WL.PStar + UL.rhoBar_eTotalTilde) / UL.rhoBar;

	<?=eqn.prim_t?> WR = primFromCons(UR, x);
	real sqrtRhoR = sqrt(WR.rhoBar);
	real3 vR = WR.vTilde;
	//real hTotalR = calc_hTotal(WR.rhoBar, WR.PStar, UR.rhoBar_eTotalTilde) - UR.ePot;
	real hTotalR = (WR.PStar + UR.rhoBar_eTotalTilde) / UR.rhoBar;

	real invDenom = 1./(sqrtRhoL + sqrtRhoR);
	
	//Roe-averaged
	real rhoBar = sqrtRhoL * sqrtRhoR;
	real3 vTilde = real3_add(
			real3_scale(vL, sqrtRhoL * invDenom),
			real3_scale(vR, sqrtRhoR * invDenom));
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);

	//derived:
	real vTildeSq = coordLenSq(vTilde, x);
	real eKin = .5 * vTildeSq;
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);

	return (<?=eqn.eigen_t?>){
		.rhoBar = rhoBar, 
		.vTilde = vTilde,
		.hTotal = hTotal,
		.k = 0,
		.omega = 0,
		.vTildeSq = vTildeSq,
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
	real v_n = vTilde.x, v_n1 = vTilde.y, v_n2 = vTilde.z;
]] 
	elseif side == 1 then
		
		prefix = [[
	const real nx = 0, ny = 1, nz = 0;
	const real n1x = 0, n1y = 0, n1z = 1;
	const real n2x = 1, n2y = 0, n2z = 0;
	real v_n = vTilde.y, v_n1 = vTilde.z, v_n2 = vTilde.x;
]] 
	elseif side == 2 then
		
		prefix = [[
	const real nx = 0, ny = 0, nz = 1;
	const real n1x = 1, n1y = 0, n1z = 0;
	const real n2x = 0, n2y = 1, n2z = 0;
	real v_n = vTilde.z, v_n1 = vTilde.x, v_n2 = vTilde.y;
]]
	end
	
	prefix = [[
	sym3 gU = coord_gU(x);
	real gUjj = gU.s]]..side..side..[[;
	real sqrt_gUjj = coord_sqrt_gU]]..side..side..[[(x);
	
	real3 vTilde = eig.vTilde;
	real3 vL = coord_lower(vTilde, x);
	real hTotal = eig.hTotal;
	real vTildeSq = real3_dot(vTilde, vL);
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
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) { 
	<?=prefix?>
	
	real denom = 2. * Cs * Cs;
	real invDenom = 1. / denom;

	const real heatRatioMinusOne = heatCapacityRatio - 1.;
<? if side == 0 then ?>
	real sqrt_gUxx = sqrt_gUjj;
	return (<?=eqn.waves_t?>){.ptr={
		(X.ptr[0] * (.5 * heatRatioMinusOne * vTildeSq + Cs * vTilde.x / sqrt_gUxx)
			+ X.ptr[1] * (-heatRatioMinusOne * vL.x - Cs / sqrt_gUxx)
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
		(X.ptr[0] * (denom - heatRatioMinusOne * vTildeSq)
			+ X.ptr[1] * 2. * heatRatioMinusOne * vL.x
			+ X.ptr[2] * 2. * heatRatioMinusOne * vL.y
			+ X.ptr[3] * 2. * heatRatioMinusOne * vL.z
			+ X.ptr[4] * -2. * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (vTilde.x * gU.xy / gU.xx - vTilde.y)
			+ X.ptr[1] * -gU.xy / gU.xx
			+ X.ptr[2],
		X.ptr[0] * (vTilde.x * gU.xz / gU.xx - vTilde.z)
			+ X.ptr[1] * -gU.xz / gU.xx
			+ X.ptr[3],
		0,
		0,
		(X.ptr[0] * (.5 * heatRatioMinusOne * vTildeSq - Cs * vTilde.x / sqrt_gUxx)
			+ X.ptr[1] * (-heatRatioMinusOne * vL.x + Cs / sqrt_gUxx)
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
	}};
<? elseif side == 1 then ?>
	real sqrt_gUyy = sqrt_gUjj;
	return (<?=eqn.waves_t?>){.ptr={
		(X.ptr[0] * (.5 * heatRatioMinusOne * vTildeSq + Cs * vTilde.y / sqrt_gUyy)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * (-heatRatioMinusOne * vL.y - Cs / sqrt_gUyy)
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (vTilde.y * gU.xy / gU.yy - vTilde.x)
			+ X.ptr[1]
			+ X.ptr[2] * -gU.xy / gU.yy,
		(X.ptr[0] * (denom - heatRatioMinusOne * vTildeSq)
			+ X.ptr[1] * 2. * heatRatioMinusOne * vL.x
			+ X.ptr[2] * 2. * heatRatioMinusOne * vL.y
			+ X.ptr[3] * 2. * heatRatioMinusOne * vL.z
			+ X.ptr[4] * -2. * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (vTilde.y * gU.yz / gU.yy - vTilde.z)
			+ X.ptr[2] * -gU.yz / gU.yy
			+ X.ptr[3],
		0,
		0,
		(X.ptr[0] * (.5 * heatRatioMinusOne * vTildeSq - Cs * vTilde.y / sqrt_gUyy)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * (-heatRatioMinusOne * vL.y + Cs / sqrt_gUyy)
			+ X.ptr[3] * -heatRatioMinusOne * vL.z
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
	}};
<? elseif side == 2 then ?>
	real sqrt_gUzz = sqrt_gUjj;
	return (<?=eqn.waves_t?>){.ptr={
		(X.ptr[0] * (.5 * heatRatioMinusOne * vTildeSq + Cs * vTilde.z / sqrt_gUzz)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * (-heatRatioMinusOne * vL.z - Cs / sqrt_gUzz)
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
		X.ptr[0] * (vTilde.z * gU.xz / gU.zz - vTilde.x)
			+ X.ptr[1]
			+ X.ptr[3] * -gU.xz / gU.zz,
		X.ptr[0] * (vTilde.z * gU.yz / gU.zz - vTilde.y)
			+ X.ptr[2]
			+ X.ptr[3] * -gU.yz / gU.zz,
		(X.ptr[0] * (denom - heatRatioMinusOne * vTildeSq)
			+ X.ptr[1] * 2. * heatRatioMinusOne * vL.x
			+ X.ptr[2] * 2. * heatRatioMinusOne * vL.y
			+ X.ptr[3] * 2. * heatRatioMinusOne * vL.z
			+ X.ptr[4] * -2. * heatRatioMinusOne
		) * invDenom,
		0,
		0,
		(X.ptr[0] * (.5 * heatRatioMinusOne * vTildeSq - Cs * vTilde.z / sqrt_gUzz)
			+ X.ptr[1] * -heatRatioMinusOne * vL.x
			+ X.ptr[2] * -heatRatioMinusOne * vL.y
			+ X.ptr[3] * (-heatRatioMinusOne * vL.z + Cs / sqrt_gUzz)
			+ X.ptr[4] * heatRatioMinusOne
		) * invDenom,
	}};
<? end ?>
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> X,
	real3 x
) {
	<?=prefix?>

<? if side == 0 then ?>	
	real sqrt_gUxx = sqrt_gUjj;
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[0] + X.ptr[1] + X.ptr[6],
		X.ptr[0] * (vTilde.x - Cs * sqrt_gUxx)
			+ X.ptr[1] * vTilde.x
			+ X.ptr[6] * (vTilde.x + Cs * sqrt_gUxx),
		X.ptr[0] * (vTilde.y - Cs * gU.xy / sqrt_gUxx)
			+ X.ptr[1] * vTilde.y
			+ X.ptr[2]
			+ X.ptr[6] * (vTilde.y + Cs * gU.xy / sqrt_gUxx),
		X.ptr[0] * (vTilde.z - Cs * gU.xz / sqrt_gUxx)
			+ X.ptr[1] * vTilde.z
			+ X.ptr[3]
			+ X.ptr[6] * (vTilde.z + Cs * gU.xz / sqrt_gUxx),
		X.ptr[0] * (hTotal - Cs * vTilde.x / sqrt_gUxx)
			+ X.ptr[1] * vTildeSq / 2.
			+ X.ptr[2] * vL.y
			+ X.ptr[3] * vL.z
			+ X.ptr[6] * (hTotal + Cs * vTilde.x / sqrt_gUxx),
		0,
	}};
<? elseif side == 1 then ?>	
	real sqrt_gUyy = sqrt_gUjj;
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[0] + X.ptr[2] + X.ptr[6],
		X.ptr[0] * (vTilde.x - Cs * gU.xy / sqrt_gUyy)
			+ X.ptr[1]
			+ X.ptr[2] * vTilde.x
			+ X.ptr[6] * (vTilde.x + Cs * gU.xy / sqrt_gUyy),
		X.ptr[0] * (vTilde.y - Cs * sqrt_gUyy)
			+ X.ptr[2] * vTilde.y
			+ X.ptr[6] * (vTilde.y + Cs * sqrt_gUyy),
		X.ptr[0] * (vTilde.z - Cs * gU.yz / sqrt_gUyy)
			+ X.ptr[2] * vTilde.z
			+ X.ptr[3]
			+ X.ptr[6] * (vTilde.z + Cs * gU.yz / sqrt_gUyy),
		X.ptr[0] * (hTotal - Cs * vTilde.y / sqrt_gUyy)
			+ X.ptr[1] * vL.x
			+ X.ptr[2] * vTildeSq / 2.
			+ X.ptr[3] * vL.z
			+ X.ptr[6] * (hTotal + Cs * vTilde.y / sqrt_gUyy),
		0,
	}};
<? elseif side == 2 then ?>	
	real sqrt_gUzz = sqrt_gUjj;
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[0] + X.ptr[3] + X.ptr[6],
		X.ptr[0] * (vTilde.x - Cs * gU.xz / sqrt_gUzz)
			+ X.ptr[1]
			+ X.ptr[3] * vTilde.x
			+ X.ptr[6] * (vTilde.x + Cs * gU.xz / sqrt_gUzz),
		X.ptr[0] * (vTilde.y - Cs * gU.yz / sqrt_gUzz)
			+ X.ptr[2]
			+ X.ptr[3] * vTilde.y
			+ X.ptr[6] * (vTilde.y + Cs * gU.yz / sqrt_gUzz),
		X.ptr[0] * (vTilde.z - Cs * sqrt_gUzz)
			+ X.ptr[3] * vTilde.z
			+ X.ptr[6] * (vTilde.z + Cs * sqrt_gUzz),
		X.ptr[0] * (hTotal - Cs * vTilde.z / sqrt_gUzz)
			+ X.ptr[1] * vL.x
			+ X.ptr[2] * vL.y
			+ X.ptr[3] * vTildeSq / 2.
			+ X.ptr[6] * (hTotal + Cs * vTilde.z / sqrt_gUzz),
		0,
	}};
<? end ?>
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) {
	<?=prefix?>
	return (<?=eqn.cons_t?>){.ptr={
		X.ptr[1] * nx 
			+ X.ptr[2] * ny 
			+ X.ptr[3] * nz,
		X.ptr[0] * (-v_n * vTilde.x + (heatCapacityRatio - 1.) * .5 * vTildeSq * gUj.x)
			+ X.ptr[1] * (vTilde.x * nx - (heatCapacityRatio - 1.) * gUj.x * vL.x + v_n)
			+ X.ptr[2] * (vTilde.x * ny - (heatCapacityRatio - 1.) * gUj.x * vL.y)
			+ X.ptr[3] * (vTilde.x * nz - (heatCapacityRatio - 1.) * gUj.x * vL.z)
			+ X.ptr[4] * (heatCapacityRatio - 1.) * nx,
		X.ptr[0] * (-v_n * vTilde.y + (heatCapacityRatio - 1.) * .5 * vTildeSq * gUj.y)
			+ X.ptr[1] * (vTilde.y * nx - (heatCapacityRatio - 1.) * gUj.y * vL.x)
			+ X.ptr[2] * (vTilde.y * ny - (heatCapacityRatio - 1.) * gUj.y * vL.y + v_n)
			+ X.ptr[3] * (vTilde.y * nz - (heatCapacityRatio - 1.) * gUj.y * vL.z)
			+ X.ptr[4] * (heatCapacityRatio - 1.) * ny,
		X.ptr[0] * (-v_n * vTilde.z + (heatCapacityRatio - 1.) * .5 * vTildeSq * gUj.z)
			+ X.ptr[1] * (vTilde.z * nx - (heatCapacityRatio - 1.) * gUj.z * vL.x)
			+ X.ptr[2] * (vTilde.z * ny - (heatCapacityRatio - 1.) * gUj.z * vL.y)
			+ X.ptr[3] * (vTilde.z * nz - (heatCapacityRatio - 1.) * gUj.z * vL.z + v_n)
			+ X.ptr[4] * (heatCapacityRatio - 1.) * nz,
		X.ptr[0] * v_n * ((heatCapacityRatio - 1.) * .5 * vTildeSq - hTotal)
			+ X.ptr[1] * (-(heatCapacityRatio - 1.) * v_n * vL.x + nx * hTotal)
			+ X.ptr[2] * (-(heatCapacityRatio - 1.) * v_n * vL.y + ny * hTotal)
			+ X.ptr[3] * (-(heatCapacityRatio - 1.) * v_n * vL.z + nz * hTotal)
			+ X.ptr[4] * heatCapacityRatio * v_n,
		0,
		0,
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
	real vTildeSq = coordLenSq(W.vTilde, x);
	real eKin = .5 * vTildeSq;
	//real hTotal = calc_hTotal(W.rhoBar, W.PStar, U.rhoBar_eTotalTilde);
	real hTotal = (W.PStar + U.rhoBar_eTotalTilde) / U.rhoBar;
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	return (<?=eqn.eigen_t?>){
		.rhoBar = W.rhoBar,
		.vTilde = W.vTilde,
		.hTotal = hTotal,
		.k = 0,
		.omega = 0,
		.vTildeSq = vTildeSq,
		.Cs = Cs,
	};
}
<? end ?>

kernel void addSource(
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
	real3 m_conn_vv = coord_conn_apply23(W.vTilde, U->rhoBar_vTilde, x);
	deriv->rhoBar_vTilde = real3_sub(deriv->rhoBar_vTilde, m_conn_vv);	//-Conn^i_jk rhoBar vTilde^j vTilde^k 
	deriv->rhoBar_vTilde = real3_sub(deriv->rhoBar_vTilde, real3_scale(coord_conn_trace23(x), W.PStar));		//-Conn^i_jk g^jk PStar
<? end ?>
}
