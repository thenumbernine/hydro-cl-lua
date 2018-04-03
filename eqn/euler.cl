/*
I've highjacked all of this.  It was a normal Euler eqn solver.
But I experimented with a curved-space solver.  
To get back to the original code,
just replace all the g_ab stuff with their constant values and simplify away.
*/

//needs Equation.hasFluxFromCons=true
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
<?=eqn.eigen_t?> eigen_forSide(
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR,
	real3 x
) {
	/*
	hmm, where should ePot be removed from eTotal?
	before Roe averaging, or after?
	I'm doing it after in my previous framework and it works fine
	but here, removing after causes numerical errors
	however in the previous framework, I never offset ePot to be always-positive like I'm doing here ... 
	...tadaa! that's the issue.
	so if we do offset potential to be positive, and don't remove potential energy at all from hTotal, the boundary stays static.
	if we do offset potential to be positive and remove it after Roe averaging then the system explodes 
		(probably because positive potential means subtracting out a positive number, means nonphysical cases ...
		... even though as we offset the potential to be positive, we are also adding it into the total energy ...
		... but somewhere through numerical error it might be getting too close to zero non-potential total energy?)
	
	if we do offset potential to be positive and remove ePot before Roe averaging then nothing special happens

	but if we DON'T offset potential (so keep it negative) and then we remove it after Roe averaging
		then we get this nice Rayleigh-Taylor ... or Jeans ... instability from self-gravitation.
	and if we DON'T offset potential (keep negative) and remove ePot before Roe averaging, same thing

	so in summary if ePot is positive then the system either blows up (remove after) or stays boring (remove before)
	but if ePot is negative then it behaves and displays the correct(?) instability patterns 
	
	but which behavior is correct?
	is the instability at the surface of my self-gravitating demos Jeans, or is it errors?
	
	--- 

	let me take that back
	with forward Euler integrator on a grid of 256x256 ...
	... even with offsetting the potential energy to be positive ...
	... it is working fine (when subtracting potential energy pre-averaging)
	... even with radius same as yesterday, gridsize, ... idk why it was exploding before ...
	
	but still,
	offsetting the ePot to be positive does reduce the gravitational turbulence a lot
	... so the initial range of ePot influences the turbulence
	*/

	<?=eqn.prim_t?> WL = primFromCons(*UL, x);
	real sqrtRhoL = sqrt(WL.rho);
	real3 vL = WL.v;
	real hTotalL = calc_hTotal(WL.rho, WL.P, UL->ETotal) - UL->ePot;
	
	<?=eqn.prim_t?> WR = primFromCons(*UR, x);
	real sqrtRhoR = sqrt(WR.rho);
	real3 vR = WR.v;
	real hTotalR = calc_hTotal(WR.rho, WR.P, UR->ETotal) - UR->ePot;

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

	return (<?=eqn.eigen_t?>){
		.rho = rho, 
		.v = v,
		.hTotal = hTotal,
		.vSq = vSq,
		.Cs = Cs,
	};
}


//this is also used by PLM for cell-centered waves
//that's why it is split out
//its use there could replace calcCellMinMaxEigenvalues
//which is called by the default calcDT
<?
for side=0,solver.dim-1 do
	for _,addrs in ipairs{
		{'', ''},
		{'global', 'global'},
	} do
		local addr0, addr1 = table.unpack(addrs)
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

<?	end
end
?>

//this routine is pretty standard.
//why not move it to Roe or somewhere similar?
kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,		//[volume][dim]
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	
	int indexR = index;

	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		
		int indexL = index - stepsize.s<?=side?>;
		
		<?= solver.getULRCode ?>
		
		int indexInt = side + dim * index;	
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + indexInt;
		*eig = eigen_forSide(UL, UR, xInt);
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
#else
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
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	global const <?=eqn.cons_t?>* U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real vSq = coordLenSq(W.v, x);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U->ETotal);
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

/*
set Equation.useSourceTerm=true
This has the connection terms that aren't absorbed in the change-of-volume
for Euclidian this is -Conn^i_jk rho v^j v^k

rho_,t + (rho v^i)_;i = 0
(rho v^j)_,t + (rho v^i v^j + P g^ij)_;i = 0
ETotal_,t + (HTotal u^i)_;i = 0

rho_,t + m^i_;i = 0
m^j_,t + (m^i v^j + P g^ij)_;i = 0
ETotal_,t + (HTotal u^i)_;i = 0

rho_,t + m^i_,i + Conn^i_ki m^k = 0
m^j_,t + (m^i v^j + P g^ij)_,i + Conn^i_ki (m^k v^j + P g^kj) + Conn^j_ki (m^i v^k + P g^ik) = 0
ETotal_,t + (HTotal u^i)_,i + Conn^i_ki (HTotal u^k) = 0

rho_,t + m^i_,i + 1/e e_,k m^k = 0
m^j_,t + (m^i v^j + P g^ij)_,i + 1/e e_,k (m^k v^j + P g^kj) + Conn^j_ki (m^i v^k + P g^ik) = 0
ETotal_,t + (HTotal u^i)_,i + 1/e e_,k (HTotal u^k) = 0
... for e = sqrt(g), for g = det(g_ij)

rho_,t + 1/e (e m^i)_,i = 0
m^j_,t + 1/e (e (m^i v^j + P g^ij))_,i = -Conn^j_ki (m^i v^k + P g^ik)
ETotal_,t + 1/e (e (HTotal u^i))_,i = 0

... source terms:
rho_,t += 0
m^j_,t += -Conn^j_ki (rho v^k v^i + P g^ki)
ETotal_,t += 0

for cylindrical holonomic:
Conn^phi_phi_r = Conn^phi_r_phi = 1/r
Conn^r_phi_phi = -r
...and the source term becomes...
m^r_,t += -Conn^r_phi_phi (rho v^phi v^phi + P g^phi^phi)
   ... += r (rho v^phi v^phi + P / r^2)
   ... += r rho v^phi v^phi + P / r
m^phi_,t += -(Conn^phi_phi_r + Conn^phi_r_phi) (rho v^r v^phi + P g^r^phi)
     ... += -2 Conn^phi_phi_r rho v^r v^phi
     ... += -2/r rho v^r v^phi


for cylindrical anholonomic normalized:
Conn^phi_r_phi = 1/r
Conn^r_phi_phi = -1/r
...and the source term becomes...
m^r_,t += -Conn^r_phi_phi (rho v^phi v^phi + P g^phi^phi)
   ... += 1/r (rho v^phi v^phi + P)
m^phi_,t += -1/r (rho v^r v^phi)


*/
kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;

#if defined(geometry_cylinder)

//these two work with holonomic geometry
//they maintain constant initial conditions with zero velocity
//it doesn't work so well with nonzero velocity
#if 1	// holonomic coriolis force alone		
	//m^j_,t += -Conn^j_ki (rho v^k v^i + P g^ki)
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	deriv->m = real3_sub(deriv->m, real3_scale(coord_conn(W.v, x), U->rho));	//-Conn^i_jk v^j v^k rho
	deriv->m = real3_sub(deriv->m, real3_scale(coord_connTrace(x), W.P));		//-Conn^i_jk g^jk P = -Conn^i P
#elif 0	// holonomic coriolis force alone ... explicitly written out for cylindrical geometry ...
	real r = x.x, theta = x.y, z = x.z;
	<?=eqn.prim_t?> W = primFromCons(*U, x); 
	deriv->m.x += r * W.rho * W.v.y * W.v.y + W.P / r;
	deriv->m.y -= 2. / r * W.rho * W.v.x * W.v.y;

#elif 0	// holonomic: all covariant derivative terms (including the others that should be absorbed into the finite-volume computations)
	real connTrace = coord_connTrace(x);
	deriv->rho -= U->rho * connTrace;
	deriv->m = real3_sub(
		real3_sub(deriv->m, real3_scale(coord_conn(U->m, x), 1. / U->rho)),
		real3_scale(U->m, connTrace));
	deriv->ETotal -= U->ETotal * connTrace;
#elif 0	//anholonomic
	real r = x.x, theta = x.y, z = x.z;
	<?=eqn.prim_t?> W = primFromCons(*U, x); 
	real HTotal = calc_HTotal(W.P, U->ETotal);
	//anholonomic: Conn^theta_r_theta = -Conn^r_theta_theta = 1/r
	deriv->rho -= 1/r * U->m.x;
	deriv->m.x -= 1/r * (U->m.x * W.v.x);
	deriv->m.y -= 1/r * (U->m.x * W.v.y);
	deriv->m.z -= 1/r * (U->m.x * W.v.z);
	deriv->ETotal -= 1/r * W.v.x * HTotal;
#elif 0	//anholonomic coriolis force alone ... explicitly written out
		// doesn't work with anholonomic geometry ... maybe because the volume terms in the flux shouldn't be there / are messed up?
	real r = x.x, theta = x.y, z = x.z;
	<?=eqn.prim_t?> W = primFromCons(*U, x); 
	deriv->m.x += 1./r * (W.rho * W.v.y * W.v.y + W.P);
	deriv->m.y -= 1./r * W.rho * W.v.x * W.v.y;
#endif
#else	//all other geometry -- general case for holonomic coordinates
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	deriv->m = real3_sub(deriv->m, real3_scale(coord_conn(W.v, x), U->rho));	//-Conn^i_jk v^j v^k rho
	deriv->m = real3_sub(deriv->m, real3_scale(coord_connTrace(x), W.P));		//-Conn^i_jk g^jk P = -Conn^i P
#endif

/*
Navier-Stokes FANS source term: 2005 Uygun, Kirkkopru
tau_ij,j
for tau_ij = mu_T (v_i,j + v_j,i - 2/3 delta_ij v_k,k) - 2/3 rho K delta_ij
so tau_ij,j = (mu_T (v_i,j + v_j,i - 2/3 delta_ij v_k,k) - 2/3 rho K delta_ij),j
	= mu_T (v_i,j + v_j,i - 2/3 delta_ij v_k,k),j - (2/3 rho K delta_ij),j
	= mu_T (v_i,jj + 1/3 v_j,ji) - 2/3 K rho,i
	= mu_T (rho m)_i,jj + 1/3 mu_T (rho m)_j,ji - 2/3 K rho,i
	= mu_T ( (rho,j m_i + rho m_i,j),j + 1/3 (rho,j m_j + rho m_j,j),i ) - 2/3 K rho,i
	= mu_T ( 
		rho,jj m_i 
		+ 1/3 rho,ij m_j 
		+ rho,j (2 m_i,j + 1/3 m_j,i)
		+ rho (m_i,jj + 1/3 m_j,ij)
	)
	+ rho,i (mu_T m_j,j - 2 K) / 3
*/
<? if solver.useNavierStokesViscosityTerm then ?>
<?=makePartial('rho', 'real');?>
<?=makePartial2('rho', 'real');?>
<?=makePartial('m', 'real3');?>
<?=makePartial2('m', 'real3');?>

<? for i,xi in ipairs(xNames) do
?>	deriv->m.<?=xi?> += 
	- 2./3. * K
<? end
?>;
<? end ?>
}
