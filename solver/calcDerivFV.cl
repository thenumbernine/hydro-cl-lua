// used by all the finite volume solvers
<?
local eqn = solver.eqn
?>

kernel void calcDerivFromFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf
) {
	typedef <?=eqn.cons_t?> cons_t;
	
	SETBOUNDS_NOGHOST();
	global cons_t* deriv = derivBuf + index;
	real3 x = cell_x(i);

<? if solver.coord.vectorComponent == 'anholonomic' then ?>
	//real volume = cell_volume(x);
<? elseif solver.coord.vectorComponent == 'cartesian' then ?>

	//TODO FIXME, this isn't correct
	real volume = cell_volume(x);
	//TODO rethink this flag
<? 	if eqn.weightFluxByGridVolume then ?>
//	real volume = coord_sqrt_det_g(x);
<? 	else ?>
//	const real volume = 1.<? 
	for i=0,solver.dim-1 do 
		?> * solver->grid_dx.s<?=i?><? 
	end ?>;
<? 	end ?>


<? else -- vectorComponent ~= 'anholonomic' ?>

<? 	if eqn.weightFluxByGridVolume then ?>
	real volume = coord_sqrt_det_g(x);
<? 	else ?>
	const real volume = 1.<? 
	for i=0,solver.dim-1 do 
		?> * solver->grid_dx.s<?=i?><? 
	end ?>;
<? 	end ?>

<? end -- vectorComponent == 'anholonomic' ?>

	<? for side=0,solver.dim-1 do ?>{
		int indexIntL = <?=side?> + dim * index;
		const global cons_t* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * solver->stepsize.s<?=side?>; 
		const global cons_t* fluxR = fluxBuf + indexIntR;
		
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

		//This is the covariant finite volume code that that represents the gradient of the metric determinant 
		//All other covariant terms should be accounted for in the equation source update
		//U^i_;t + F^ij_;j  = 0
		//U^i_,t + F^ij_,j + Gamma^j_kj F^ik + Gamma^i1_kj F^i1^k + ... + Gamma^in_kj F^in^k = 0
		//					(metric det gradient) 
<? if solver.coord.vectorComponent == 'anholonomic' then ?>
		//real areaL = cell_area<?=side?>(xIntL);
		//real areaR = cell_area<?=side?>(xIntR);
<? elseif solver.coord.vectorComponent == 'cartesian' then ?>
		real areaL = cell_area<?=side?>(xIntL);
		real areaR = cell_area<?=side?>(xIntR);
<? else -- vectorComponent ~= 'anholonomic' ?>
<? if eqn.weightFluxByGridVolume then ?>
		real areaL;
		if (coord_sqrt_det_g(xIntL) == 0) {
			areaL = 0;
		} else {
			areaL = coord_sqrt_det_g(xIntL) / solver->grid_dx.s<?=side?>;
		}
		real areaR;
		if (coord_sqrt_det_g(xIntR) == 0) {
			areaR = 0;
		} else {
			areaR = coord_sqrt_det_g(xIntR) / solver->grid_dx.s<?=side?>;
		}
<? else ?>
		real areaL, areaR;
		if (volume == 0) {
			areaL = areaR = 0;
		} else {
			areaL = volume / solver->grid_dx.s<?=side?>;
			areaR = volume / solver->grid_dx.s<?=side?>;
		}
<? end ?>
<? end	-- vectorComponent == 'anholonomic' ?>

<? if solver.coord.vectorComponent == 'anholonomic' then ?>

		//divide by integral of volume form
		real volInt = cell_volume(x);
		//scale flux by volume form times area
		real areaR = cell_area0(xIntR);
		real areaL = cell_area0(xIntL);
		
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				areaR * fluxR->ptr[j] 
				- areaL * fluxL->ptr[j]
			) / volInt;
		}
<? elseif not eqn.postComputeFluxCode then -- would the compiler know to optimize this?
?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
<? else ?>
		cons_t flux;
		for (int j = 0; j < numIntStates; ++j) {
			flux.ptr[j] = (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
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
