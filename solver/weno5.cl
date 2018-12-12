// WENO5 solver:
// courtesy of Mara
// and some matlab code

inline real sqr(real x) { return x * x; }

<? 
local clnumber = require 'cl.obj.number'

local c = {
	{11/6, -7/6,  1/3},
	{ 1/3,  5/6, -1/6},
	{-1/6,  5/6,  1/3},
	{ 1/3, -7/6, 11/6},
}
local d = { 3/10, 3/5, 1/10, }

for _,l_or_r in ipairs{'l', 'r'} do
	local cofs = l_or_r == 'l' and 1 or 2
	local d0 = l_or_r == 'l' and 3 or 1
	local dd = l_or_r == 'l' and -1 or 1
?>
<?=eqn.waves_t?> weno5<?=l_or_r?>(const <?=eqn.waves_t?>* v) {
	<?=eqn.waves_t?> result;
	for (int k = 0; k < numWaves; ++k) {
		real beta0 = (13./12.)*sqr( v[2].ptr[k] - 2*v[3].ptr[k] + v[4].ptr[k]) + (1./4.)*sqr(3*v[2].ptr[k] - 4*v[3].ptr[k] +   v[4].ptr[k]);
		real beta1 = (13./12.)*sqr( v[1].ptr[k] - 2*v[2].ptr[k] + v[3].ptr[k]) + (1./4.)*sqr(  v[1].ptr[k]                 -   v[3].ptr[k]);
		real beta2 = (13./12.)*sqr( v[0].ptr[k] - 2*v[1].ptr[k] + v[2].ptr[k]) + (1./4.)*sqr(  v[0].ptr[k] - 4*v[1].ptr[k] + 3*v[2].ptr[k]);

<?
if solver.weno5method == '1996 Jiang Shu' then 
	local epsilon = clnumber(1e-14)
?>
		real w0 = <?=clnumber(d[d0 + 0 * dd])?> / sqr(<?=epsilon?> + beta0);
		real w1 = <?=clnumber(d[d0 + 1 * dd])?> / sqr(<?=epsilon?> + beta1);
		real w2 = <?=clnumber(d[d0 + 2 * dd])?> / sqr(<?=epsilon?> + beta2);
<? elseif solver.weno5method == '2008 Borges' then 
	local epsilon = clnumber(1e-14)
?>
		real tau5 = fabs(beta0 - beta2);
		real w0 = <?=clnumber(d[d0 + 0 * dd])?> * (1 + (tau5 / (beta0 + <?=epsilon?>)));
		real w1 = <?=clnumber(d[d0 + 1 * dd])?> * (1 + (tau5 / (beta1 + <?=epsilon?>)));
		real w2 = <?=clnumber(d[d0 + 2 * dd])?> * (1 + (tau5 / (beta2 + <?=epsilon?>)));
<? elseif solver.weno5method == '2010 Shen Zha' then 
	local epsilon = clnumber(1e-10)
	local shen_zha_A = clnumber(50)	-- 0-100
?>
		real minB = min(min(beta0, beta1), beta2);
		real maxB = max(max(beta0, beta1), beta2);
		real R0 = minB / (maxB + <?=epsilon?>);
		beta0 += R0 * <?=shen_zha_A?> * minB;
		beta1 += R0 * <?=shen_zha_A?> * minB;
		beta2 += R0 * <?=shen_zha_A?> * minB;
		real w0 = <?=clnumber(d[d0 + 0 * dd])?> / sqr(<?=epsilon?> + beta0);
		real w1 = <?=clnumber(d[d0 + 1 * dd])?> / sqr(<?=epsilon?> + beta1);
		real w2 = <?=clnumber(d[d0 + 2 * dd])?> / sqr(<?=epsilon?> + beta2);
<? else
	error("unknown weno5method "..tostring(solver.weno5method))
end ?>

		real vs0 = <?=clnumber(c[cofs+0][1])?> * v[2].ptr[k] + <?=clnumber(c[cofs+0][2])?> * v[3].ptr[k] + <?=clnumber(c[cofs+0][3])?> * v[4].ptr[k];
		real vs1 = <?=clnumber(c[cofs+1][1])?> * v[1].ptr[k] + <?=clnumber(c[cofs+1][2])?> * v[2].ptr[k] + <?=clnumber(c[cofs+1][3])?> * v[3].ptr[k];
		real vs2 = <?=clnumber(c[cofs+2][1])?> * v[0].ptr[k] + <?=clnumber(c[cofs+2][2])?> * v[1].ptr[k] + <?=clnumber(c[cofs+2][3])?> * v[2].ptr[k];
	
		real wtot = w0 + w1 + w2;
		result.ptr[k] = (w0 * vs0 + w1 * vs1 + w2 * vs2) / wtot;
	}
	return result;
}
<? end ?>

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost-1,numGhost-1);
	
	real3 xR = cell_x(i);
	int indexR = index;
	const global <?=eqn.cons_t?>* UR = UBuf + indexR;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		
		int indexL = index - solver->stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;

		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;

		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, normalForSide<?=side?>());

		real maxAbsLambda = 0.;
		for (int j = 0; j < 6; ++j) {
			const global <?=eqn.cons_t?>* Uj = U + (j-2) * solver->stepsize.s<?=side?>;
			
			<?=eqn:consWaveCodePrefix(side, '*Uj', 'xInt'):gsub('\n', '\n\t\t')?>
			
			real lambdaMin = <?=eqn:consMinWaveCode(side, '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMin));
			
			real lambdaMax = <?=eqn:consMaxWaveCode(side, '*Uj', 'xInt')?>;
			maxAbsLambda = max(maxAbsLambda, fabs(lambdaMax));
		}

<? if false then ?>
		<?=eqn.cons_t?> Fp[6], Fm[6];
		for (int j = 0; j < 6; ++j) {
			int indexj = index + (j-2) * solver->stepsize.s<?=side?>;
			const global <?=eqn.cons_t?>* Uj = UBuf + index;
			<?=eqn.cons_t?> F = fluxFromCons_<?=side?>(solver, *U, x);
			for (int k = 0; k < numStates; ++k) {
				Fp[j].ptr[k] = (F.ptr[k] + maxAbsLambda * Uj->ptr[k]) * .5;
				Fm[j].ptr[k] = (F.ptr[k] - maxAbsLambda * Uj->ptr[k]) * .5;
			}
		}

		<?=eqn.waves_t?> fp[5], fm[5];
		for (int j = 0; j < 5; ++j) {
			fp[j] = eigen_leftTransform_<?=side?>(solver, eig, Fp[j], xInt);
			fm[j] = eigen_leftTransform_<?=side?>(solver, eig, Fm[j+1], xInt);
		}

<? else ?>
	<? for j=0,5 do ?>
		const global <?=eqn.cons_t?>* U<?=j?> = U + <?=j-2?> * solver->stepsize.s<?=side?>;
		<?=eqn.cons_t?> F<?=j?> = fluxFromCons_<?=side?>(solver, *U<?=j?>, xInt);
		<?=eqn.cons_t?> Fp<?=j?>, Fm<?=j?>;
		<? for k=0,eqn.numStates-1 do ?>
		Fp<?=j?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] + maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
		Fm<?=j?>.ptr[<?=k?>] = (F<?=j?>.ptr[<?=k?>] - maxAbsLambda * U<?=j?>->ptr[<?=k?>]) * .5;
		<? end ?>
	<? end ?>

		<?=eqn.waves_t?> fp[5], fm[5];
	<? for j=0,4 do ?>
		fp[<?=j?>] = eigen_leftTransform_<?=side?>(solver, eig, Fp<?=j?>, xInt);
		fm[<?=j?>] = eigen_leftTransform_<?=side?>(solver, eig, Fm<?=j+1?>, xInt);
	<? end ?>

<? end ?>
	
		<?=eqn.waves_t?> wfp = weno5r(fp);
		<?=eqn.waves_t?> wfm = weno5l(fm);
		<?=eqn.waves_t?> wf;
		for (int j = 0; j < numWaves; ++j) {
			wf.ptr[j] = wfp.ptr[j] + wfm.ptr[j];
		}
		
		int fluxIndexInt = side + dim * (index
			+ solver->stepsize.s<?=side?>
		);
		global <?=eqn.cons_t?>* flux = fluxBuf + fluxIndexInt;
		*flux = eigen_rightTransform_<?=side?>(solver, eig, wf, xInt);
	}<? end ?>
}
