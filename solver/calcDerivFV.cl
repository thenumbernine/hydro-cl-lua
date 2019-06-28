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

<? if solver.coord.anholonomic then ?>
<? else -- anholonomic ?>
<? if eqn.weightFluxByGridVolume then ?>
	real volume = coord_sqrt_det_g(x);
<? else ?>
	const real volume = 1.<? 
	for i=0,solver.dim-1 do 
		?> * solver->grid_dx.s<?=i?><? 
	end ?>;
<? end ?>
<? end -- anholonomic ?>

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
<? if solver.coord.anholonomic then ?>
<? else -- anholonomic ?>
<? if eqn.weightFluxByGridVolume then ?>	
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

		real areaL = coord_sqrt_det_g(xIntL) / solver->grid_dx.s<?=side?>;
		real areaR = coord_sqrt_det_g(xIntR) / solver->grid_dx.s<?=side?>;
<? else ?>
		real areaL = volume / solver->grid_dx.s<?=side?>;
		real areaR = volume / solver->grid_dx.s<?=side?>;
<? end ?>
<? end	-- anholonomic ?>

<? if solver.coord.anholonomic then ?>
<? 	if require 'coord.cylinder'.is(solver.coord) then ?>
<?		if true then -- 2018 Shadab et al ?>
<? 			if side == 0 then ?>
		real rR = x.x + .5 * solver->grid_dx.x;
		real rL = x.x - .5 * solver->grid_dx.x;
		real volume = .5 * (rR * rR - rL * rL);
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (fluxR->ptr[j] * rR - fluxL->ptr[j] * rL) / volume;
		}
<? 			elseif side == 1 then ?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (fluxR->ptr[j] - fluxL->ptr[j]) / (x.x * solver->grid_dx.y);
		}
<? 			elseif side == 2 then ?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (fluxR->ptr[j] - fluxL->ptr[j]) / (solver->grid_dx.z);
		}
<? 			end ?>
<?		else	-- my attempt ?>
		real rR = x.x + .5 * solver->grid_dx.x;
		real rL = x.x - .5 * solver->grid_dx.x;
		real deltaVolume = .5 * (rR * rR - rL * rL) * solver->grid_dx.y * solver->grid_dx.z;
<? 			if side == 0 then ?>
		real deltaAreaR = rR * solver->grid_dx.y * solver->grid_dx.z;
		real deltaAreaL = rL * solver->grid_dx.y * solver->grid_dx.z;
<? 			elseif side == 1 then ?>
		real deltaArea = solver->grid_dx.x * solver->grid_dx.z;
		real deltaAreaR = deltaArea;
		real deltaAreaL = deltaArea;
<? 			elseif side == 2 then ?>
		real deltaArea = .5 * (rR * rR - rL * rL) * solver->grid_dx.y;
		real deltaAreaR = deltaArea;
		real deltaAreaL = deltaArea;
<? 			end ?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (deltaAreaR * fluxR->ptr[j] - deltaAreaL * fluxL->ptr[j]) / deltaVolume;
		}
<?		end ?>
<? 	else ?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
<? 	end ?>
<? elseif not eqn.postComputeFluxCode then -- would the compiler know to optimize this?
?>
		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
<? else ?>
		<?=eqn.cons_t?> flux;
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
