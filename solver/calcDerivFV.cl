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
	//real volume = cell_volume(x);
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
		
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

		//This is the covariant finite volume code that that represents the gradient of the metric determinant 
		//All other covariant terms should be accounted for in the equation source update
		//U^i_;t + F^ij_;j  = 0
		//U^i_,t + F^ij_,j + Gamma^j_kj F^ik + Gamma^i1_kj F^i1^k + ... + Gamma^in_kj F^in^k = 0
		//					(metric det gradient) 
<? if solver.coord.anholonomic then ?>
		//real areaL = cell_area<?=side?>(xIntL);
		//real areaR = cell_area<?=side?>(xIntR);
<? else -- anholonomic ?>
<? if eqn.weightFluxByGridVolume then ?>
		real areaL = coord_sqrt_det_g(xIntL) / solver->grid_dx.s<?=side?>;
		real areaR = coord_sqrt_det_g(xIntR) / solver->grid_dx.s<?=side?>;
<? else ?>
		real areaL = volume / solver->grid_dx.s<?=side?>;
		real areaR = volume / solver->grid_dx.s<?=side?>;
<? end ?>
<? end	-- anholonomic ?>

<? if solver.coord.anholonomic then ?>

<?	if require 'coord.cylinder'.is(solver.coord) then ?>
		real rR = x.x + .5 * solver->grid_dx.x;
		real rL = x.x - .5 * solver->grid_dx.x;
		// integral of volume form across cell
		real volInt = .5 * (rR * rR - rL * rL) * solver->grid_dx.y * solver->grid_dx.z;
		//real volInt = x.x * solver->grid_dx.x * solver->grid_dx.y * solver->grid_dx.z;

		// integral of volume form across surface
<? 		if side == 0 then ?>
		// right volume form 'r' times integral across other coords phi, z (which don't integrate the volume form)
		real areaR = rR * solver->grid_dx.y * solver->grid_dx.z;
		// left volume form 'r' times integral across other coords phi, z (which don't integrate the volume form)
		real areaL = rL * solver->grid_dx.y * solver->grid_dx.z;
<? 		elseif side == 1 then ?>
		// integral of volume form 'r' across coords r, z
		real area = .5 * (rR * rR - rL * rL) * solver->grid_dx.z;
		real areaR = area;
		real areaL = area;
<? 		elseif side == 2 then ?>
		// integral of volume form 'r' across coords r, phi
		real area = .5 * (rR * rR - rL * rL) * solver->grid_dx.y;
		real areaR = area;
		real areaL = area;
<?		end ?>

<?	else -- automatic ?>
		
		//divide by integral of volume form
		real volInt = cell_volume(x);
		//scale flux by volume form times area
		real areaR = cell_area0(xIntR);
		real areaL = cell_area0(xIntL);

<?	end ?>
		
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
