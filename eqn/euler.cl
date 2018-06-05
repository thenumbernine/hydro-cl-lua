/*
I've highjacked all of this.  It was a normal Euler eqn solver.
But I experimented with a curved-space solver.  
To get back to the original code,
just replace all the g_ab stuff with their constant values and simplify away.
*/
<?
local solver = eqn.solver

local clnumber = require 'cl.obj.number'
local makePartials = require 'eqn.makepartial'

local derivOrder, makePartial, makePartial2
if eqn.guiVars.useNavierStokesViscosityTerm.value then 
	derivOrder = 2 * solver.numGhost
	makePartial = function(...) return makePartials.makePartial(derivOrder, solver, ...) end
	makePartial2 = function(...) return makePartials.makePartial2(derivOrder, solver, ...) end
end
?>

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
<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forSide_<?=side?>(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
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
<? end ?>

//this routine is pretty standard.
//why not move it to Roe or somewhere similar?
kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,		//[numCells][dim]
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize.s<?=side?>;
		
		<?=solver:getULRCode()?>
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		int indexInt = side + dim * index;	
		eigenBuf[indexInt] = eigen_forSide_<?=side?>(*UL, *UR, xInt);
	}<? end ?>
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
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> X,
	real3 x
) { 
	<?=prefix?>
	
	real denom = 2. * Cs * Cs;
	real invDenom = 1. / denom;

#if 0	//works but isn't correct for curved space
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
#else	//works for covariant formulation of Euler fluid equations
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

<? if not require 'coord.cartesian'.is(solver.coord) then ?>
	//connection coefficient source terms of covariant derivative w/contravariant velocity vectors in a holonomic coordinate system
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real3 m_conn_vv = coord_conn_apply23(W.v, U->m, x);
	deriv->m = real3_sub(deriv->m, m_conn_vv);	//-Conn^i_jk rho v^j v^k 
	deriv->m = real3_sub(deriv->m, real3_scale(coord_conn_trace23(x), W.P));		//-Conn^i_jk g^jk P
	//deriv->m = real3_sub(deriv->m, real3_scale(coord_conn_apply12(W.v, U->m, x), heatCapacityRatio - 1.));	//-(gamma-1) rho v^j v^k Conn_jk^i
	//deriv->ETotal -= (heatCapacityRatio - 1.) * real3_dot(coord_lower(U->m, x), m_conn_vv);
<? end ?>

/*
Navier-Stokes FANS source term: 2005 Uygun, Kirkkopru

NOTICE this is only for flat space!

i'th term in the j'th direction
F[m]_ij += tauTilde_ij
F[ETotal]_j += ThetaTilde_j

for tau_ij = mu_T (v_i,j + v_j,i) - 2/3 delta_ij (mu_T v_k,k + K rho) (eqn 3)
Theta_j = tau_ij v_i + (kBar_L + k_T) * TTilde,j
... but eqn 3 isn't used numerically ... eqn 11 is ... and that doesn't look so straightforward as the eqn above ... ( why does the tauTilde_xx term have a vTilde_y,y term?)



so tau_ij,j = (mu_T (v_i,j + v_j,i - 2/3 delta_ij v_k,k) - 2/3 rho K delta_ij),j
	= mu_T (v_i,j + v_j,i - 2/3 delta_ij v_k,k),j - (2/3 rho K delta_ij),j
	= mu_T (v_i,jj + 1/3 v_j,ji) - 2/3 K rho,i
	= mu_T (m_i / rho),jj + 1/3 mu_T (m_j / rho),ji - 2/3 K rho,i
	= mu_T (m_i,j / rho - m_i / rho^2 rho,j),j + 1/3 mu_T (m_j,j / rho - m_j / rho^2 rho,j),i - 2/3 K rho,i
	= mu_T (
		m_i,jj / rho 
		- m_i,j / rho^2 rho,j
		- m_i,j / rho^2 rho,j
		+ 2 m_i / rho^3 (rho,j)^2
		- m_i / rho^2 rho,jj
	) + 1/3 mu_T (
		m_j,ji / rho 
		- m_j,j / rho^2 rho,i
		- m_j,i / rho^2 rho,j
		+ 2 m_j / rho^3 rho,i rho,j
		- m_j / rho^2 rho,ij
	) - 2/3 K rho,i
	
	= mu_T (
		(m_i,jj + 1/3 m_j,ji) / rho 
		- ((2 m_i,j + 1/3 m_j,i) rho,j + 1/3 m_j,j rho,i) / rho^2
		+ (
			m_i (2 (rho,j)^2 - rho rho,jj)
			+ m_j (2 rho,i rho,j - rho rho,ij)
		) / rho^3
	) - 2/3 K rho,i

	= 
		- 2/3 K rho,i
		+ mu_T (
			+ m_i,jj / rho
			+ 1/3 m_j,ji / rho
			- 2 m_i,j rho,j / rho^2 
			- 1/3 m_j,i rho,j / rho^2 
			- 1/3 m_j,j rho,i / rho^2
			+ m_i 2 (rho,j)^2 / rho^3 
			- m_i rho,jj / rho^2
			+ 2 m_j rho,i rho,j / rho^3 
			- m_j rho,ij / rho^2
		)
*/
<? if eqn.guiVars.useNavierStokesViscosityTerm.value then 
?>
	const real K = <?=clnumber(eqn.guiVars.viscosity_K.value)?>;
	const real mu_T = <?=clnumber(eqn.guiVars.viscosity_mu_T.value)?>;

<?=makePartial('rho', 'real')?>
<?=makePartial2('rho', 'real')?>
<?=makePartial('m', 'real3')?>		//m_i,j = partial_m_l[j].i
<?=makePartial2('m', 'real3')?>		//m_i,jk = partial2_m_ll[jk].i

	real invRho = 1. / U->rho;
	real invRho2 = invRho * invRho;
	real invRho3 = invRho * invRho2;

<? for i,xi in ipairs(xNames) do
?>	deriv->m.<?=xi?> += 
		-2./3. * K * partial_rho_l[<?=i-1?>]
	+ mu_T * (0.
<?	for j,xj in ipairs(xNames) do
		local ij = from3x3to6(i,j)
		local jj = from3x3to6(j,j)
?>		
		+ invRho * (
			+ partial2_m_ll[<?=jj-1?>].<?=xi?>
			+ 1./3. * partial2_m_ll[<?=ij-1?>].<?=xj?>
			+ invRho * (
				- 2. * partial_m_l[<?=j-1?>].<?=xi?> * partial_rho_l[<?=j-1?>]
				- 1./3. * partial_m_l[<?=i-1?>].<?=xj?> * partial_rho_l[<?=j-1?>]
				- 1./3. * partial_m_l[<?=j-1?>].<?=xj?> * partial_rho_l[<?=i-1?>]
				- U->m.<?=xi?> * partial2_rho_ll[<?=jj-1?>]
				- U->m.<?=xj?> * partial2_rho_ll[<?=ij-1?>]
				+ invRho * 2. * partial_rho_l[<?=j-1?>] * (
					+ U->m.<?=xi?> * partial_rho_l[<?=j-1?>] 
					+ U->m.<?=xj?> * partial_rho_l[<?=i-1?>]
				)
			)
		)
<? 	end
?>	);
<? end
?>	
<? end ?>
}
