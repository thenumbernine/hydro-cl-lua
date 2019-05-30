// used by all the finite volume solvers
<?
local eqn = solver.eqn
?>

kernel void calcDerivFromFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	
	real3 x = cell_x(i);
<? if eqn.weightFluxByGridVolume then ?>
	real volume = coord_volume(x);
<? else ?>
	const real volume = 1.<? for i=0,solver.dim-1 do ?> * solver->grid_dx.s<?=i?><? end ?>;
<? end ?>

	<? for side=0,solver.dim-1 do ?>{
		int indexIntL = <?=side?> + dim * index;
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * solver->stepsize.s<?=side?>; 
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + indexIntR;

		//This is the covariant finite volume code that that represents the gradient of the metric determinant 
		//All other covariant terms should be accounted for in the equation source update
		//U^i_;t + F^ij_;j  = 0
		//U^i_,t + F^ij_,j + Gamma^j_kj F^ik + Gamma^i1_kj F^i1^k + ... + Gamma^in_kj F^in^k = 0
		//					(metric det gradient) 

<? if eqn.weightFluxByGridVolume then ?>	
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	
		real areaL = coord_volume(xIntL);
		real areaR = coord_volume(xIntR);

<? else ?>
		real areaL = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
		real areaR = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
<? end ?>

<? if not eqn.postComputeFluxCode then -- would the compiler know to optimize this? ?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / (volume * solver->grid_dx.s<?=side?>);
		}
<? else ?>
		<?=eqn.cons_t?> flux;
		for (int j = 0; j < numIntStates; ++j) {
			flux.ptr[j] = (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / (volume * solver->grid_dx.s<?=side?>);
		}

		{
<?=eqn.postComputeFluxCode or ''?>
		}

		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= flux.ptr[j];
		}
<? end ?>

	}<? end ?>
}
