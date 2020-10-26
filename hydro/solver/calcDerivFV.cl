//// MODULE_NAME: calcDerivFromFlux
//// MODULE_DEPENDS: solver_t cons_t cell_t solver.macros SETBOUNDS_NOGHOST

<? if not require 'hydro.solver.meshsolver'.is(solver) then ?>

// used by all the finite volume solvers

kernel void calcDerivFromFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	typedef <?=eqn.cons_t?> cons_t;
	
	SETBOUNDS_NOGHOST();
	global cons_t* deriv = derivBuf + index;
	real3 x = cellBuf[index].pos;

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
<? if solver.coord.vectorComponent == 'holonomic'
or require 'hydro.coord.cartesian'.is(solver.coord)
then ?>
	real volume = 1.<?
	for i=0,solver.dim-1 do
		?> * solver->grid_dx.s<?=i?><?
	end
?>;
<? else ?>
//// MODULE_DEPENDS: cell_volume 
	real volume = cell_volume(x);
<? end ?>

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
<? if solver.coord.vectorComponent == 'holonomic'
or require 'hydro.coord.cartesian'.is(solver.coord)
then ?>
		real areaL, areaR;
		areaL = areaR = 1.<?
	for i=0,solver.dim-1 do
		if i ~= side then
			?> * solver->grid_dx.s<?=i?><?
		end
	end
?>;
<? else ?>
//// MODULE_DEPENDS: cell_area#
		real areaL = cell_area<?=side?>(xIntL);
		real areaR = cell_area<?=side?>(xIntR);
<? end ?>

		if (volume > 1e-7) {
			real invVolume = 1. / volume;
			if (areaL <= 1e-7) areaL = 0.;
			if (areaR <= 1e-7) areaR = 0.;
			areaL *= invVolume;
			areaR *= invVolume;
<? -- TODO get rid of this, it's only used by the maxwell and glm-maxwell eqns
if eqn.postComputeFluxCode then ?>
			cons_t flux = {.ptr={0}};
			for (int j = 0; j < numIntStates; ++j) {
				flux.ptr[j] = 
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL;
			}
<?=eqn.postComputeFluxCode?>
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


<? else -- meshsolver ?>


kernel void calcDerivFromFlux(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* fluxBuf,
//mesh-specific parameters:	
	const global <?=solver.coord.cell_t?>* cells,			//[numCells]
	const global <?=solver.coord.face_t?>* faces,			//[numFaces]
	const global int* cellFaceIndexes	//[numCellFaceIndexes]
) {
	int cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	const global <?=solver.coord.cell_t?>* cell = cells + cellIndex;
	
	global <?=eqn.cons_t?>* deriv = derivBuf + cellIndex;
	
	for (int i = 0; i < cell->faceCount; ++i) {
		int ei = cellFaceIndexes[i + cell->faceOffset];
		const global <?=solver.coord.face_t?>* e = faces + ei;
		
		const global <?=eqn.cons_t?>* flux = fluxBuf + ei;
		real areaOverVolume = e->area / cell->volume;
		
		if (cellIndex == e->cells.s0) {
//std::cout << " ... - " << *e << std::endl;
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] -= flux->ptr[j] * areaOverVolume;
			}
		} else {
//std::cout << " ... + " << *e << std::endl;
			for (int j = 0; j < numIntStates; ++j) {
				deriv->ptr[j] += flux->ptr[j] * areaOverVolume;
			}
		}
	}
}


<? end -- mesh vs grid solver ?>
