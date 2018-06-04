// used by all the finite volume solvers
<?
local eqn = solver.eqn
?>

kernel void calcDerivFromFlux(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	
	real3 x = cell_x(i);
<? if eqn.weightFluxByGridVolume then ?>	
	real sqrt_det_g = sqrt_det_g_grid(x);
<? else ?>
	const real sqrt_det_g = 1.;
<? end ?>

	<? for side=0,solver.dim-1 do ?>{
		int indexIntL = <?=side?> + dim * index;
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * stepsize.s<?=side?>; 
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + indexIntR;

		//This is the covariant finite volume code that that represents the gradient of the metric determinant 
		//All other covariant terms should be accounted for in the equation source update
		//U^i_;t + F^ij_;j  = 0
		//U^i_,t + F^ij_,j + Gamma^j_kj F^ik + Gamma^i1_kj F^i1^k + ... + Gamma^in_kj F^in^k = 0
		//					(metric det gradient) 

<? if eqn.weightFluxByGridVolume then ?>	
		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real sqrt_det_g_intL = sqrt_det_g_grid(xIntL);
		real areaL = sqrt_det_g_intL / grid_dx<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		real sqrt_det_g_intR = sqrt_det_g_grid(xIntR);
		real areaR = sqrt_det_g_intR / grid_dx<?=side?>;
<? else ?>
		const real areaL = 1.;
		const real areaR = 1.;
<? end ?>

		<?=eqn.cons_t?> flux;
		for (int j = 0; j < numIntStates; ++j) {
			flux.ptr[j] = (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / sqrt_det_g;
		}

		{
<?=eqn.postComputeFluxCode or ''?>
		}

		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= flux.ptr[j];
		}
	
	}<? end ?>
}
