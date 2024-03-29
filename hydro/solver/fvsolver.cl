//// MODULE_NAME: <?=calcDerivFromFlux?>
//// MODULE_DEPENDS: <?=solver_macros?>

// used by all the finite volume solvers

kernel void <?=calcDerivFromFlux?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const fluxBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;

/*<?--[[
volume vs area ...
in Cartesian grids the volume is dx * dy * dz  = int dx dy dz
	volume = int |J| |dy/dx| dx0 dx1 dx2
	area_j = int |S| |dy/dx| dx-not-j

in curvilinear holonomic (coordinate)
	- basis is given by the int (coordinate volume form) dx dy dz
	- 
in curvilinear anholonomic basis 
	- basis is given by the int (coordinate volume form) dx dy dz
	- some flux rescaling for non-scalars is needed
in curvilinear coords, cartesian basis: 
	- basis is given by the int (coordinate volume form) dx dy dz
	- 

--]]?>*/
<? if solver.coord.vectorComponent == "holonomic"
or require "hydro.coord.cartesian":isa(solver.coord)
then ?>
	real const volume = 1.<?
	for i=0,solver.dim-1 do
		?> * solver->grid_dx.s<?=i?><?
	end
?>;
<? else ?>
	real const volume = cell->volume;
<? end ?>

	<? for side=0,solver.dim-1 do ?>{
		int const indexIntL = <?=side?> + dim * index;
		global <?=cons_t?> const * const fluxL = fluxBuf + indexIntL;
		
		int const indexIntR = indexIntL + dim * solver->stepsize.s<?=side?>; 
		global <?=cons_t?> const * const fluxR = fluxBuf + indexIntR;
		
		real3 xIntL = x; xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real3 xIntR = x; xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;

		//This is the covariant finite volume code that that represents the gradient of the metric determinant 
		//All other covariant terms should be accounted for in the equation source update
		//U^i_;t + F^ij_;j  = 0
		//U^i_,t + F^ij_,j + Gamma^j_kj F^ik + Gamma^i1_kj F^i1^k + ... + Gamma^in_kj F^in^k = 0
		//					(metric det gradient) 
<? if solver.coord.vectorComponent == "holonomic"
or require "hydro.coord.cartesian":isa(solver.coord)
then ?>
		real areaL = 1.<?
	for i=0,solver.dim-1 do
		if i ~= side then
			?> * solver->grid_dx.s<?=i?><?
		end
	end
?>;
		real areaR = areaL;
<? else ?>
//// MODULE_DEPENDS: <?=cell_area_i?>
		real areaL = cell_area<?=side?>(xIntL);
		real areaR = cell_area<?=side?>(xIntR);
<? end ?>

		/*
		for
		gridSize={512,512}
		mins={-1,1}
		maxs={1,1}
		grid_dx = {2/512,2/512}
		volume=1/65536
		invVolume=65536
		areaL = 256
		areaR = 256
		somewhere half is failing to do this math
		
		ok it turns out 65519 is the largest number stored in half precision
		any larger and it just maps to inf
		thats why 512x512 doesn't work, 512x511 doesn't work, but 512x510 does work
		the limit is the volume
		
		how to circumvent this?  for higher grid resolutions, maintain volume by increasing mins/maxs of the initcond domain.
		*/
		if (volume > 1e-7) {
			real const invVolume = 1. / volume;
			if (areaL <= 1e-7) areaL = 0.;
			if (areaR <= 1e-7) areaR = 0.;
			areaL *= invVolume;
			areaR *= invVolume;
<? -- TODO get rid of this, it's only used by the maxwell and glm-maxwell eqns
if eqn.postComputeFluxCode then ?>
			<?=cons_t?> flux = {.ptr={0}};
			for (int j = 0; j < numIntStates; ++j) {
				flux.ptr[j] = 
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL;
			}
<?=eqn:postComputeFluxCode()?>
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] -= flux.ptr[j];
			}

<? else -- not postComputeFluxCode ?>
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] -= 
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL;
			}
<? end -- postComputeFluxCode ?>	
		}
	}<? end ?>
}
