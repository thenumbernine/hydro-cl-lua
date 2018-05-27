// used by all the finite volume solvers

<?
-- this should add the connection terms Conn^k_jk u^I_,j (for the I'th conserved quantity u^I)
local weightFluxByGridVolume = true
?>

kernel void calcDerivFromFlux(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	
	real3 x = cell_x(i);
<? if weightFluxByGridVolume then ?>	
	real volume = volume_at(x);
<? else ?>
	const real volume = 1.;
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

<? if weightFluxByGridVolume then ?>	
		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * grid_dx<?=side?>;
		real volumeIntL = volume_at(xIntL);
		real areaL = volumeIntL / grid_dx<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * grid_dx<?=side?>;
		real volumeIntR = volume_at(xIntR);
		real areaR = volumeIntR / grid_dx<?=side?>;
<? else ?>
		const real areaL = 1.;
		const real areaR = 1.;
<? end ?>

		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
	}<? end ?>
}
