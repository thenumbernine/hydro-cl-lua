// WENO5 solver:
// courtesy of Mara

inline real sqr(real x) { return x * x; }

<? 
local clnumber = require 'cl.obj.number'
for name,info in pairs{
	l = {
		c = {
			{11./6., -7./6.,  1./3. },
			{ 1./3.,  5./6., -1./6. },
			{-1./6.,  5./6.,  1./3. } 
		},
		d = { 0.1, 0.6, 0.3 },
	},
	r = {
		c = {
			{ 1./3.,  5./6., -1./6. },
			{-1./6.,  5./6.,  1./3. },
			{ 1./3., -7./6., 11./6. }
		},
		d = { 0.3, 0.6, 0.1 },
	},
} do 
	local c, d = info.c, info.d
?>
<?=eqn.waves_t?> weno5<?=name?>(const <?=eqn.waves_t?>* v) {
	const real eps = 1e-6;
	<?=eqn.waves_t?> result;
	for (int k = 0; k < numWaves; ++k) {
		real B0 = 
			 (13./12.)*sqr(  v[2].ptr[k] - 2*v[3].ptr[k] +   v[4].ptr[k]) +
			 ( 1./ 4.)*sqr(3*v[2].ptr[k] - 4*v[3].ptr[k] +   v[4].ptr[k]);
		real B1 =
			 (13./12.)*sqr(  v[1].ptr[k] - 2*v[2].ptr[k] +   v[3].ptr[k]) +
			 ( 1./ 4.)*sqr(  v[1].ptr[k]                 -   v[3].ptr[k]);
		real B2 = 
			 (13./12.)*sqr(  v[0].ptr[k] - 2*v[1].ptr[k] +   v[2].ptr[k]) +
			 ( 1./ 4.)*sqr(  v[0].ptr[k] - 4*v[1].ptr[k] + 3*v[2].ptr[k]);

		real vs0 = <?=clnumber(c[1][1])?> * v[2].ptr[k] + <?=clnumber(c[1][2])?> * v[3].ptr[k] + <?=clnumber(c[1][3])?> * v[4].ptr[k];
		real vs1 = <?=clnumber(c[2][1])?> * v[1].ptr[k] + <?=clnumber(c[2][2])?> * v[2].ptr[k] + <?=clnumber(c[2][3])?> * v[3].ptr[k];
		real vs2 = <?=clnumber(c[3][1])?> * v[0].ptr[k] + <?=clnumber(c[3][2])?> * v[1].ptr[k] + <?=clnumber(c[3][3])?> * v[2].ptr[k];

		real w0 = <?=clnumber(d[1])?> / sqr(eps + B0);
		real w1 = <?=clnumber(d[2])?> / sqr(eps + B1);
		real w2 = <?=clnumber(d[3])?> / sqr(eps + B2);
		
		real wtot = w0 + w1 + w2;
		result.ptr[k] = (w0*vs0 + w1*vs1 + w2*vs2)/wtot;
	}
	return result;
}
<? end ?>

kernel void calcFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* fluxBuf,
	<?= solver.getULRArg ?>,
	realparam dt
) {
	SETBOUNDS(numGhost,numGhost);
	real3 xR = cell_x(i);
	int indexR = index;
	const global <?=eqn.cons_t?>* UR = UBuf + indexR;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;	
		real dt_dx = dt / solver->grid_dx.s<?=side?>;//dx<?=side?>_at(i);
		int indexL = index - solver->stepsize.s<?=side?>;
		const global <?=eqn.cons_t?>* UL = UBuf + indexL;
		
		real3 xInt = xR;
		xInt.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		
		<?=eqn.eigen_t?> eig = eigen_forInterface(solver, *UL, *UR, xInt, normalForSide<?=side?>());

		real maxAbsLambda = 0.;
		for (int j = 0; j < 6; ++j) {
			const global <?=eqn.cons_t?>* U = UBuf + indexR + (j-3) * solver->stepsize.s<?=side?>;
			<?=eqn:consWaveCodePrefix(side, '*U', 'xInt')?>
			for (int k = 0; k < numWaves; ++k) {
				real lambdaMin = <?=eqn:consMinWaveCode(side, '*U', 'xInt')?>;
				real lambdaMax = <?=eqn:consMaxWaveCode(side, '*U', 'xInt')?>;
				lambdaMin = fabs(lambdaMin);
				lambdaMax = fabs(lambdaMax);
				real lambda = max(lambdaMin, lambdaMax);
				maxAbsLambda = max(maxAbsLambda, lambda);
			}
		}

		<?=eqn.cons_t?> Fp[6], Fm[6];
		for (int j = 0; j < 6; ++j) {
			const global <?=eqn.cons_t?>* U = UBuf + indexR + (j-3) * solver->stepsize.s<?=side?>;
			<?=eqn.cons_t?> F = fluxFromCons_<?=side?>(solver, *U, xInt);
			for (int k = 0; k < numStates; ++k) {
				Fp[j].ptr[k] = (F.ptr[k] + maxAbsLambda * U->ptr[k]) * .5;
				Fm[j].ptr[k] = (F.ptr[k] - maxAbsLambda * U->ptr[k]) * .5;
			}
		}

		<?=eqn.waves_t?> fp[6], fm[6];
		for (int j = 0; j < 6; ++j) {
			fp[j] = eigen_leftTransform_<?=side?>(solver, eig, Fp[j], xInt);
			fm[j] = eigen_leftTransform_<?=side?>(solver, eig, Fm[j], xInt);
		}

		<?=eqn.waves_t?> wfp = weno5r(fp);
		<?=eqn.waves_t?> wfm = weno5l(fm+1);
		<?=eqn.waves_t?> wf;
		for (int j = 0; j < numWaves; ++j) {
			wf.ptr[j] = .5 * (wfp.ptr[j] + wfm.ptr[j]);
		}

		int indexInt = side + dim * index;
		global <?=eqn.cons_t?>* flux = fluxBuf + indexInt;
		*flux = eigen_rightTransform_<?=side?>(solver, eig, wf, xInt);
	}<? end ?>
}
