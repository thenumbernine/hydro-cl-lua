//needs Equation.hasFluxFromCons=true
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(<?=eqn.cons_t?> U, real3 x) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real vj = W.v.s<?=side?>;
	real HTotal = U.ETotal + W.P;
	
	<?=eqn.cons_t?> F;
	F.rho = U.m.s<?=side?>;
	F.m = real3_scale(U.m, vj);
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

//used for interface eigen basis
void eigen_forSide(
	global <?=eqn.eigen_t?>* eig,
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	<?=eqn.prim_t?> WL = primFromCons(*UL, x);
	real sqrtRhoL = sqrt(WL.rho);
	real3 vL = WL.v;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL->ETotal);
	
	<?=eqn.prim_t?> WR = primFromCons(*UR, x);
	real sqrtRhoR = sqrt(UR->rho);
	real3 vR = WR.v;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR->ETotal);
	
	real invDenom = 1./(sqrtRhoL + sqrtRhoR);
	
	//Roe-averaged
	real rho = sqrtRhoL * sqrtRhoR;
	real3 v = real3_add(
			real3_scale(vL, sqrtRhoL * invDenom),
			real3_scale(vR, sqrtRhoR * invDenom));
	real hTotal = invDenom * (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR);
	
	//derived:
	real vSq = coordLenSq(v, x);
	real eKin = .5 * vSq;
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	
	eig->rho = rho; 
	eig->v = v;
	eig->hTotal = hTotal;
	eig->vSq = vSq;
	eig->Cs = Cs;
}

//this is also used by PLM for cell-centered waves
//that's why it is split out
//its use there could replace calcCellMinMaxEigenvalues
//which is called by the default calcDT
<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr0?> real* wave,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	real3 x
) {
	real Cs_sqrt_gU = eig->Cs * coord_sqrt_gU<?=side..side?>(x);
	real v_n = eig->v.s[<?=side?>];
	wave[0] = v_n - Cs_sqrt_gU;
	wave[1] = v_n;
	wave[2] = v_n;
	wave[3] = v_n;
	wave[4] = v_n + Cs_sqrt_gU;
}
<?		end
	end
end
?>

kernel void calcEigenBasis(
	global real* waveBuf,			//[volume][dim][numWaves]
	global <?=eqn.eigen_t?>* eigenBuf,		//[volume][dim]
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(2,1);
	real3 x = cell_x(i);
	
	int indexR = index;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize[side];
	
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;

		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		eigen_forSide(eig, UL, UR, xInt);
		
		global real* wave = waveBuf + numWaves * indexInt;
		eigen_calcWaves_<?=side?>_global_global(wave, eig, xInt);
	}<? end ?>
}

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
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
	
	real3 v = eig->v;
	real3 vL = coord_lower(v, x);
	real hTotal = eig->hTotal;
	real vSq = real3_dot(v, vL);
	real Cs = eig->Cs;
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

void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) { 
	<?=prefix?>

	real denom = 2. * Cs * Cs;
	real invDenom = 1. / denom;

#if 0	//works but isn't correct for curved space
	Y[0] = (X[0] * ((heatCapacityRatio - 1.) * .5 * vSq + Cs_over_sqrt_gUjj * v_n)
		+ X[1] * -(nx * Cs_over_sqrt_gUjj + (heatCapacityRatio - 1.) * vL.x) 
		+ X[2] * -(ny * Cs_over_sqrt_gUjj + (heatCapacityRatio - 1.) * vL.y)
		+ X[3] * -(nz * Cs_over_sqrt_gUjj + (heatCapacityRatio - 1.) * vL.z)
		+ X[4] * (heatCapacityRatio - 1.)
	) * invDenom;
	Y[1] = (X[0] * (2.*Cs*Cs - (heatCapacityRatio - 1.) * vSq)
		+ X[1] * (heatCapacityRatio - 1.) * v.x * 2
		+ X[2] * (heatCapacityRatio - 1.) * v.y * 2
		+ X[3] * (heatCapacityRatio - 1.) * v.z * 2
		+ X[4] * -(heatCapacityRatio - 1.) * 2
	) * invDenom;
	Y[2] = X[0] * -v_n1 + X[1] * n1x + X[2] * n1y + X[3] * n1z;
	Y[3] = X[0] * -v_n2 + X[1] * n2x + X[2] * n2y + X[3] * n2z;
	Y[4] = (X[0] * ((heatCapacityRatio - 1.) * .5 * vSq - Cs_over_sqrt_gUjj * v_n) 
		+ X[1] * (nx * Cs_over_sqrt_gUjj - (heatCapacityRatio - 1.) * vL.x) 
		+ X[2] * (ny * Cs_over_sqrt_gUjj - (heatCapacityRatio - 1.) * vL.y) 
		+ X[3] * (nz * Cs_over_sqrt_gUjj - (heatCapacityRatio - 1.) * vL.z) 
		+ X[4] * (heatCapacityRatio - 1.)
	) * invDenom;
#else	//matches math

//these two blocks of code are identical
// except the first one stalls the kernel for seconds
// and the second one is just a mathematical simplification of the first - and works fine
#if 0
<? if side == 0 then ?>
	real CsSq = Cs * Cs;
	real sqrt_gUxx = sqrt_gUjj;
	const real heatRatioMinusOne = heatCapacityRatio - 1.;
	Y[0] = X[0] * (heatRatioMinusOne * vSq / (4. * CsSq) + v.x / (2. * Cs * sqrt_gUxx))
		+ X[1] * (-(heatRatioMinusOne * vL.x) / (2. * CsSq) - 1. / (2. * Cs * sqrt_gUxx))
		+ X[2] * -heatRatioMinusOne * vL.y / (2. * CsSq)
		+ X[3] * -heatRatioMinusOne * vL.z / (2. * CsSq)
		+ X[4] * heatRatioMinusOne / (2. * CsSq);
	Y[1] = X[0] * (1. - heatRatioMinusOne * vSq / (2. * CsSq))
		+ X[1] * heatRatioMinusOne * vL.x / CsSq
		+ X[2] * heatRatioMinusOne * vL.y / CsSq
		+ X[3] * heatRatioMinusOne * vL.z / CsSq
		+ X[4] * -heatRatioMinusOne / CsSq;
	Y[2] = X[0] * (v.x * gU.xy / gU.xx - v.y)
		+ X[1] * -gU.xy / gU.xx
		+ X[2];
	Y[3] = X[0] * (v.x * gU.xz / gU.xx - v.z)
		+ X[1] * -gU.xz / gU.xx
		+ X[3];
	Y[4] = X[0] * (heatRatioMinusOne * vSq / (4. * CsSq) - v.x / (2. * Cs * sqrt_gUxx))
		+ X[1] * (-heatRatioMinusOne * vL.x / (2. * CsSq) + 1. / (2. * Cs * sqrt_gUxx))
		+ X[2] * -heatRatioMinusOne * vL.y / (2. * CsSq)
		+ X[3] * -heatRatioMinusOne * vL.z / (2. * CsSq)
		+ X[4] * heatRatioMinusOne / (2. * CsSq);
<? elseif side == 1 then ?>
<? elseif side == 2 then ?>
<? end ?>
#else	//good
	const real heatRatioMinusOne = heatCapacityRatio - 1.;
<? if side == 0 then ?>
	real sqrt_gUxx = sqrt_gUjj;
	Y[0] = (X[0] * (.5 * heatRatioMinusOne * vSq + Cs * v.x / sqrt_gUxx)
		+ X[1] * (-heatRatioMinusOne * vL.x - Cs / sqrt_gUxx)
		+ X[2] * -heatRatioMinusOne * vL.y
		+ X[3] * -heatRatioMinusOne * vL.z
		+ X[4] * heatRatioMinusOne
	) * invDenom;
	Y[1] = (X[0] * (denom - heatRatioMinusOne * vSq)
		+ X[1] * 2. * heatRatioMinusOne * vL.x
		+ X[2] * 2. * heatRatioMinusOne * vL.y
		+ X[3] * 2. * heatRatioMinusOne * vL.z
		+ X[4] * -2. * heatRatioMinusOne
	) * invDenom;
	Y[2] = X[0] * (v.x * gU.xy / gU.xx - v.y)
		+ X[1] * -gU.xy / gU.xx
		+ X[2];
	Y[3] = X[0] * (v.x * gU.xz / gU.xx - v.z)
		+ X[1] * -gU.xz / gU.xx
		+ X[3];
	Y[4] = (X[0] * (.5 * heatRatioMinusOne * vSq - Cs * v.x / sqrt_gUxx)
		+ X[1] * (-heatRatioMinusOne * vL.x + Cs / sqrt_gUxx)
		+ X[2] * -heatRatioMinusOne * vL.y
		+ X[3] * -heatRatioMinusOne * vL.z
		+ X[4] * heatRatioMinusOne
	) * invDenom;
<? elseif side == 1 then ?>
	real sqrt_gUyy = sqrt_gUjj;
	Y[0] = (X[0] * (.5 * heatRatioMinusOne * vSq + Cs * v.y / sqrt_gUyy)
		+ X[1] * -heatRatioMinusOne * vL.x
		+ X[2] * (-heatRatioMinusOne * vL.y - Cs / sqrt_gUyy)
		+ X[3] * -heatRatioMinusOne * vL.z
		+ X[4] * heatRatioMinusOne
	) * invDenom;
	Y[1] = X[0] * (v.y * gU.xy / gU.yy - v.x)
		+ X[1]
		+ X[2] * -gU.xy / gU.yy;
	Y[2] = (X[0] * (denom - heatRatioMinusOne * vSq)
		+ X[1] * 2. * heatRatioMinusOne * vL.x
		+ X[2] * 2. * heatRatioMinusOne * vL.y
		+ X[3] * 2. * heatRatioMinusOne * vL.z
		+ X[4] * -2. * heatRatioMinusOne
	) * invDenom;
	Y[3] = X[0] * (v.y * gU.yz / gU.yy - v.z)
		+ X[2] * -gU.yz / gU.yy
		+ X[3];
	Y[4] = (X[0] * (.5 * heatRatioMinusOne * vSq - Cs * v.y / sqrt_gUyy)
		+ X[1] * -heatRatioMinusOne * vL.x
		+ X[2] * (-heatRatioMinusOne * vL.y + Cs / sqrt_gUyy)
		+ X[3] * -heatRatioMinusOne * vL.z
		+ X[4] * heatRatioMinusOne
	) * invDenom;
<? elseif side == 2 then ?>
	real sqrt_gUzz = sqrt_gUjj;
	Y[0] = (X[0] * (.5 * heatRatioMinusOne * vSq + Cs * v.z / sqrt_gUzz)
		+ X[1] * -heatRatioMinusOne * vL.x
		+ X[2] * -heatRatioMinusOne * vL.y
		+ X[3] * (-heatRatioMinusOne * vL.z - Cs / sqrt_gUzz)
		+ X[4] * heatRatioMinusOne
	) * invDenom;
	Y[1] = X[0] * (v.z * gU.xz / gU.zz - v.x)
		+ X[1]
		+ X[3] * -gU.xz / gU.zz;
	Y[2] = X[0] * (v.z * gU.yz / gU.zz - v.y)
		+ X[2]
		+ X[3] * -gU.yz / gU.zz;
	Y[3] = (X[0] * (denom - heatRatioMinusOne * vSq)
		+ X[1] * 2. * heatRatioMinusOne * vL.x
		+ X[2] * 2. * heatRatioMinusOne * vL.y
		+ X[3] * 2. * heatRatioMinusOne * vL.z
		+ X[4] * -2. * heatRatioMinusOne
	) * invDenom;
	Y[4] = (X[0] * (.5 * heatRatioMinusOne * vSq - Cs * v.z / sqrt_gUzz)
		+ X[1] * -heatRatioMinusOne * vL.x
		+ X[2] * -heatRatioMinusOne * vL.y
		+ X[3] * (-heatRatioMinusOne * vL.z + Cs / sqrt_gUzz)
		+ X[4] * heatRatioMinusOne
	) * invDenom;
<? end ?>
#endif

#endif
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) {
	<?=prefix?>
#if 0	//works but doesn't account for curved space
	Y[0] = X[0] + X[1] + X[4];
	Y[1] = X[0] * (v.x - gUj.x * Cs_over_sqrt_gUjj) 
		+ X[1] * v.x 
		+ X[2] * n1x 
		+ X[3] * n2x 
		+ X[4] * (v.x + gUj.x * Cs_over_sqrt_gUjj);
	Y[2] = X[0] * (v.y - gUj.y * Cs_over_sqrt_gUjj) 
		+ X[1] * v.y 
		+ X[2] * n1y 
		+ X[3] * n2y 
		+ X[4] * (v.y + gUj.y * Cs_over_sqrt_gUjj);
	Y[3] = X[0] * (v.z - gUj.z * Cs_over_sqrt_gUjj) 
		+ X[1] * v.z 
		+ X[2] * n1z 
		+ X[3] * n2z 
		+ X[4] * (v.z + gUj.z * Cs_over_sqrt_gUjj);
	Y[4] = X[0] * (hTotal - v_n * Cs_over_sqrt_gUjj) 
		+ X[1] * .5 * vSq
		+ X[2] * v_n1 
		+ X[3] * v_n2 
		+ X[4] * (hTotal + v_n * Cs_over_sqrt_gUjj);
#else	//math
<? if side == 0 then ?>	
	real sqrt_gUxx = sqrt_gUjj;
	Y[0] = X[0] + X[1] + X[4];
	Y[1] = X[0] * (v.x - Cs * sqrt_gUxx)
		+ X[1] * v.x
		+ X[4] * (v.x + Cs * sqrt_gUxx);
	Y[2] = X[0] * (v.y - Cs * gU.xy / sqrt_gUxx)
		+ X[1] * v.y
		+ X[2]
		+ X[4] * (v.y + Cs * gU.xy / sqrt_gUxx);
	Y[3] = X[0] * (v.z - Cs * gU.xz / sqrt_gUxx)
		+ X[1] * v.z
		+ X[3]
		+ X[4] * (v.z + Cs * gU.xz / sqrt_gUxx);
	Y[4] = X[0] * (hTotal - Cs * v.x / sqrt_gUxx)
		+ X[1] * vSq / 2.
		+ X[2] * vL.y
		+ X[3] * vL.z
		+ X[4] * (hTotal + Cs * v.x / sqrt_gUxx);
<? elseif side == 1 then ?>	
	real sqrt_gUyy = sqrt_gUjj;
	Y[0] = X[0] + X[2] + X[4];
	Y[1] = X[0] * (v.x - Cs * gU.xy / sqrt_gUyy)
		+ X[1]
		+ X[2] * v.x
		+ X[4] * (v.x + Cs * gU.xy / sqrt_gUyy);
	Y[2] = X[0] * (v.y - Cs * sqrt_gUyy)
		+ X[2] * v.y
		+ X[4] * (v.y + Cs * sqrt_gUyy);
	Y[3] = X[0] * (v.z - Cs * gU.yz / sqrt_gUyy)
		+ X[2] * v.z
		+ X[3]
		+ X[4] * (v.z + Cs * gU.yz / sqrt_gUyy);
	Y[4] = X[0] * (hTotal - Cs * v.y / sqrt_gUyy)
		+ X[1] * vL.x
		+ X[2] * vSq / 2.
		+ X[3] * vL.z
		+ X[4] * (hTotal + Cs * v.y / sqrt_gUyy);
<? elseif side == 2 then ?>	
	real sqrt_gUzz = sqrt_gUjj;
	Y[0] = X[0] + X[3] + X[4];
	Y[1] = X[0] * (v.x - Cs * gU.xz / sqrt_gUzz)
		+ X[1]
		+ X[3] * v.x
		+ X[4] * (v.x + Cs * gU.xz / sqrt_gUzz);
	Y[2] = X[0] * (v.y - Cs * gU.yz / sqrt_gUzz)
		+ X[2]
		+ X[3] * v.y
		+ X[4] * (v.y + Cs * gU.yz / sqrt_gUzz);
	Y[3] = X[0] * (v.z - Cs * sqrt_gUzz)
		+ X[3] * v.z
		+ X[4] * (v.z + Cs * sqrt_gUzz);
	Y[4] = X[0] * (hTotal - Cs * v.z / sqrt_gUzz)
		+ X[1] * vL.x
		+ X[2] * vL.y
		+ X[3] * vSq / 2.
		+ X[4] * (hTotal + Cs * v.z / sqrt_gUzz);
<? end ?>
#endif
}

<?	if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* Y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* X,
	real3 x
) {
	<?=prefix?>
	
	Y[0] = X[1] * nx 
		+ X[2] * ny 
		+ X[3] * nz;
	Y[1] = X[0] * (-v_n * v.x + (heatCapacityRatio - 1.) * .5 * vSq * gUj.x)
		+ X[1] * (v.x * nx - (heatCapacityRatio - 1.) * gUj.x * vL.x + v_n)
		+ X[2] * (v.x * ny - (heatCapacityRatio - 1.) * gUj.x * vL.y)
		+ X[3] * (v.x * nz - (heatCapacityRatio - 1.) * gUj.x * vL.z)
		+ X[4] * (heatCapacityRatio - 1.) * nx;
	Y[2] = X[0] * (-v_n * v.y + (heatCapacityRatio - 1.) * .5 * vSq * gUj.y)
		+ X[1] * (v.y * nx - (heatCapacityRatio - 1.) * gUj.y * vL.x)
		+ X[2] * (v.y * ny - (heatCapacityRatio - 1.) * gUj.y * vL.y + v_n)
		+ X[3] * (v.y * nz - (heatCapacityRatio - 1.) * gUj.y * vL.z)
		+ X[4] * (heatCapacityRatio - 1.) * ny;
	Y[3] = X[0] * (-v_n * v.z + (heatCapacityRatio - 1.) * .5 * vSq * gUj.z)
		+ X[1] * (v.z * nx - (heatCapacityRatio - 1.) * gUj.z * vL.x)
		+ X[2] * (v.z * ny - (heatCapacityRatio - 1.) * gUj.z * vL.y)
		+ X[3] * (v.z * nz - (heatCapacityRatio - 1.) * gUj.z * vL.z + v_n)
		+ X[4] * (heatCapacityRatio - 1.) * nz;
	Y[4] = X[0] * v_n * ((heatCapacityRatio - 1.) * .5 * vSq - hTotal)
		+ X[1] * (-(heatCapacityRatio - 1.) * v_n * vL.x + nx * hTotal)
		+ X[2] * (-(heatCapacityRatio - 1.) * v_n * vL.y + ny * hTotal)
		+ X[3] * (-(heatCapacityRatio - 1.) * v_n * vL.z + nz * hTotal)
		+ X[4] * heatCapacityRatio * v_n;
}
<?
				end
			end
		end
	end
end
?>


// used by PLM


<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxForCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real mi = U.m.s<?=side?>;
	return (<?=eqn.cons_t?>){
		.rho = mi,
		.m = (real3){
			.x = mi * W.v.x + coord_gU0<?=side?>(x) * W.P,
			.y = mi * W.v.y + coord_gU1<?=side?>(x) * W.P,
			.z = mi * W.v.z + coord_gU2<?=side?>(x) * W.P,
		},
		.ETotal = (U.ETotal + W.P) * W.v.s<?=side?>,
		.ePot = 0.,
	};
}
<? end ?>

<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	<?=eqn.eigen_t?>* eig,
	global const <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real vSq = coordLenSq(W.v, x);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U->ETotal);
	real CsSq = (heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	eig->rho = W.rho;
	eig->v = W.v;
	eig->hTotal = hTotal;
	eig->vSq = vSq;
	eig->Cs = Cs;
}
<? end ?>

//TODO where does potential energy belong in this Jacobian?
// does it belong here at all?

//U = output
//WA = W components that make up the jacobian matrix
//W = input
//x = coordinate location
void apply_dU_dW(
	<?=eqn.cons_t?>* U, 
	const <?=eqn.prim_t?>* WA, 
	const <?=eqn.prim_t?>* W, 
	real3 x
) {
	real3 WA_vL = coord_lower(WA->v, x);
	*U = (<?=eqn.cons_t?>){
		.rho = W->rho,
		.m = real3_add(
			real3_scale(WA->v, W->rho), 
			real3_scale(W->v, WA->rho)),
		.ETotal = W->rho * .5 * real3_dot(WA->v, WA_vL) 
			+ WA->rho * real3_dot(W->v, WA_vL)
			+ W->P / (heatCapacityRatio - 1.)
			+ WA->rho * W->ePot,
		.ePot = W->ePot,
	};
}

//W = output
//WA = W components that make up the jacobian matrix
//U = input
//x = coordinate location
void apply_dW_dU(
	<?=eqn.prim_t?>* W,
	const <?=eqn.prim_t?>* WA,
	const <?=eqn.cons_t?>* U,
	real3 x
) {
	real3 WA_vL = coord_lower(WA->v, x);
	*W = (<?=eqn.prim_t?>){
		.rho = U->rho,
		.v = real3_sub(
			real3_scale(U->m, 1. / WA->rho),
			real3_scale(WA->v, U->rho / WA->rho)),
		.P = (heatCapacityRatio - 1.) * (
			.5 * real3_dot(WA->v, WA_vL) * U->rho 
			- real3_dot(U->m, WA_vL)
			+ U->ETotal 
			- WA->rho * U->ePot),
		.ePot = U->ePot,
	};
}

//set Equation.useSourceTerm=true
//this has the connection terms that aren't absorbed in the change-of-volum
//for Euclidian this is -Conn^i_jk rho v^j v^k
//hmm, doesn't seem to help
kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	deriv->m = real3_sub(deriv->m, real3_scale(coord_conn(U->m, x), 1. / U->rho));
}
