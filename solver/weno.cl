/*
WENO solver:
sources:
1998 Shu "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws"
"A hybrid approach for the regularized long wave-Burgers equation"
2016 Rathan, Raju "An improved Non-linear Weights for Seventh-Order WENO Scheme"
https://github.com/jzrake/Mara for weno5 examples
https://github.com/wme7/WENO7-Z/blob/master/WENO7ZresAdv1d.m for weno7 examples
https://github.com/python-hydro/hydro_examples/blob/master/compressible/weno_coefficients.py likewise
*/

inline real sqr(real x) { return x * x; }

<? 
local clnumber = require 'cl.obj.number'

local stencilSize = solver.stencilSize

local coeffs = {
	[3] = {
		c = {
			{11/6, -7/6,  2/6},
			{ 2/6,  5/6, -1/6},
			{-1/6,  5/6,  2/6},
		},
		d = {1/10, 6/10, 3/10},
		betaCoeffs = {
			{
				{4/3},
				{-19/3,25/3},
				{11/3,-31/3,10/3},
			}, {
				{4/3},
				{-13/3,13/3},
				{5/3,-13/3,4/3},
			}, {
				{10/3},
				{-31/3,25/3},
				{11/3,-19/3,4/3},
			},
		},
	},
	[4] = {
		c = {	
			{ 25/12, -23/12,  13/12, -3/12, },
			{  3/12,  13/12,  -5/12,  1/12, },
			{ -1/12,   7/12,   7/12, -1/12, },
			{  1/12,  -5/12,  13/12,  3/12, },
		},
		d = {1/35, 12/35, 18/35, 4/35},
		betaCoeffs = {

--[[
			{
				{ 2107/240},
				{-9402/240,  11003/240},
				{ 7042/240, -17246/240,   7043/240},
				{-1854/240,   4642/240,  -3882/240,    547/240},
			},
			{
				{  547/240},
				{-2522/240,   3443/240},
				{ 1922/240,  -5966/240,   2843/240},
				{ -494/240,   1602/240,  -1642/240,    267/240},
			},
			{
				{  267/240},
				{-1642/240,   2843/240},
				{ 1602/240,  -5966/240,   3443/240},
				{ -494/240,   1922/240,  -2522/240,    547/240},
			},
			{
				{  547/240},
				{-3882/240,   7043/240},
				{ 4642/240, -17246/240,  11003/240},
				{-1854/240,   7042/240,  -9402/240,   2107/240},
			},
--]]		
-- [[ 2016 Rathan, Raju after eqn 12 ... 
-- https://github.com/wme7/WENO7-Z/blob/master/WENO7ZresAdv1d.m
			{
				{  6649.},
				{-47214.,	  85641.},
				{ 56694.,	-210282.,	134241.},
				{-22778.,	  86214.,	-114894.,	25729.},
			},
			{
				{  3169.},
				{-19374.,	 33441.},
				{ 19014.,	-70602.,	 41001.},
				{ -5978.,	 23094.,	-30414.,	6649.},
			},
			{
				{  6649.},
				{-30414.,	 41001.},
				{ 23094.,	-70602.,	 33441.},
				{ -5978.,	 19014.,	-19374.,	3169.},
			},
			{
				{  25729.},
				{-114894.,	134241.},
				{  86214.,	-210282.,	 85641.},
				{ -22778.,	  56694.,	-47214.,	6649.},
			},
--]]	
		},
	},
}

local coeff = coeffs[stencilSize]
local c = coeff.c
local d = coeff.d
local betaCoeffs = coeff.betaCoeffs

for _,l_or_r in ipairs{'l', 'r'} do
	local ci0 = l_or_r == 'l' and stencilSize or 1
	local cid = l_or_r == 'l' and -1 or 1
	local cj0 = l_or_r == 'l' and 1 or stencilSize
	local cjd = l_or_r == 'l' and 1 or -1
	local d0 = l_or_r == 'l' and stencilSize or 1
	local dd = l_or_r == 'l' and -1 or 1
?>
<?=eqn.waves_t?> weno_<?=l_or_r?>(const <?=eqn.waves_t?>* v) {
	<?=eqn.waves_t?> result;
	for (int k = 0; k < numWaves; ++k) {
<? 	for j=0,stencilSize-1 do 
?>		real beta<?=j?> = 0.<?
		for m=0,stencilSize-1 do
			for n=0,m do
			?> + <?=clnumber(betaCoeffs[j+1][m+1][n+1])?> * v[<?=j+m?>].ptr[k] * v[<?=j+n?>].ptr[k]<?
			end
		end?>;
<?	end
	
	if solver.wenoMethod == '1996 Jiang Shu' then 		-- WENO-JS
		--local epsilon = clnumber(1e-14)
		local epsilon = clnumber(1e-6)
?>
<? 		for i=0,stencilSize-1 do
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> / sqr(<?=epsilon?> + beta<?=i?>);
<? 		end 

	elseif solver.wenoMethod == '2008 Borges' then 	-- WENO-Z
		local epsilon = clnumber(1e-6)
		if false then 	--if stencilSize == 4 then -- for 2016 Rathan, it suggests these:
?>		real tau = fabs(beta0 + 3. * beta1 - 3. * beta2 - beta3);
<?		else	-- for weno7, 2018 Zeytinoglu suggests this:
?>		real tau = fabs(beta0 - beta<?=stencilSize-1?>);
<?		end
 		for i=0,stencilSize-1 do 
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> * (1. + (tau / (beta<?=i?> + <?=epsilon?>)));
<? 		end

	elseif solver.wenoMethod == '2010 Shen Zha' then -- WENO-BS?
		local epsilon = clnumber(1e-10)
		local shen_zha_A = clnumber(50)	-- 0-100
?>		real minB = min(min(beta0, beta1), beta2);
		real maxB = max(max(beta0, beta1), beta2);
		real R0 = minB / (maxB + <?=epsilon?>);
<? 		for i=0,stencilSize-1 do 
?>		beta<?=i?> += R0 * <?=shen_zha_A?> * minB;
<? 		end 
 		for i=0,stencilSize-1 do 
?>		real alpha<?=i?> = <?=clnumber(d[d0 + i * dd])?> / sqr(<?=epsilon?> + beta<?=i?>);
<?	 	end 
 
	else
		error("unknown wenoMethod "..tostring(solver.wenoMethod))
	end

	for i=0,stencilSize-1 do 
?>		real vs<?=i?> = 0.<?
		for j=0,stencilSize-1 do
		?> + <?=clnumber(c[ci0+cid*i][cj0+cjd*j])?> * v[<?=i+j?>].ptr[k]<?
		end ?>;
<? 
	end 
?>		real alphasum= 0.<? for i=0,stencilSize-1 do ?> + alpha<?=i?><? end ?>;
		result.ptr[k] = (0.<?
	for i=0,stencilSize-1 do
		?> + alpha<?=i?> * vs<?=i?><?
	end ?>) / alphasum;
	}
	return result;
}
<? end ?>

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost-1);
	
	real3 xR = cell_x(i);
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?
for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - solver->stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;
	
		int indexR = index;
		const global <?=eqn.cons_t?>* UR = UBuf + indexR;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, normalForSide<?=side?>());

		real maxAbsLambda = 0.;
		for (int j = 0; j < <?=2*stencilSize?>; ++j) {
			const global <?=eqn.cons_t?>* Uj = U + (j - <?=stencilSize?>) * solver->stepsize.s<?=side?>;
			
			<?=eqn:consWaveCodePrefix(side, '*Uj', 'xInt'):gsub('\n', '\n\t\t')?>
			
			real lambdaMin = <?=eqn:consMinWaveCode(side, '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMin));
			
			real lambdaMax = <?=eqn:consMaxWaveCode(side, '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMax));
		}

	<? for j=0,2*stencilSize-1 do ?>
		const global <?=eqn.cons_t?>* U<?=j?> = U + <?=j - stencilSize?> * solver->stepsize.s<?=side?>;
		<?=eqn.cons_t?> F<?=j?> = fluxFromCons_<?=side?>(solver, *U<?=j?>, xInt);
<?
		if j < 2*stencilSize-1 then
?>		<?=eqn.cons_t?> Fp<?=j?>;
<?			for k=0,eqn.numStates-1 do 
?>		Fp<?=j?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] + maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
<?			end
		end
		if j > 0 then
?>		<?=eqn.cons_t?> Fm<?=j?>;
<?			for k=0,eqn.numStates-1 do 
?>		Fm<?=j?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] - maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
<? 			end
		end
	end ?>

		<?=eqn.waves_t?> fp[<?=2*stencilSize-1?>], fm[<?=2*stencilSize-1?>];
<? 	for j=0,2*stencilSize-2 do
?>		fp[<?=j?>] = eigen_leftTransform_<?=side?>(solver, eig, Fp<?=j?>, xInt);
		fm[<?=j?>] = eigen_leftTransform_<?=side?>(solver, eig, Fm<?=j+1?>, xInt);
<? 	end ?>
	
		<?=eqn.waves_t?> wf;
		<?=eqn.waves_t?> wfp = weno_r(fp);
		<?=eqn.waves_t?> wfm = weno_l(fm);
		for (int j = 0; j < numWaves; ++j) {
			wf.ptr[j] = wfp.ptr[j] + wfm.ptr[j];
		}
		
		int fluxIndexInt = side + dim * index;
		global <?=eqn.cons_t?>* flux = fluxBuf + fluxIndexInt;
		*flux = eigen_rightTransform_<?=side?>(solver, eig, wf, xInt);
	}<? 
end ?>
}
