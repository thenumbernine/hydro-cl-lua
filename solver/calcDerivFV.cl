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
<? if solver.coord.vectorComponent == 'cartesian' 
	or solver.coord.vectorComponent == 'anholonomic'
then ?>
	real volume = cell_volume(x);
<? else ?>
	real volume = 1.<?
	for i=1,solver.dim do
		?> * solver->grid_dx.s<?=i-1?><?
	end
?>;
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
<? if solver.coord.vectorComponent == 'cartesian' 
	or solver.coord.vectorComponent == 'anholonomic'
then ?>
		real areaL = cell_area<?=side?>(xIntL);
		real areaR = cell_area<?=side?>(xIntR);
<? else ?>
		real areaL, areaR;
		areaL = areaR = 1.<?
	for i=1,solver.dim do
		if i ~= side then
			?> * solver->grid_dx.s<?=i-1?><?	
		end
	end
?>;
<? end ?>

		for (int j = 0; j < numIntStates; ++j) {
			deriv->ptr[j] -= (
				fluxR->ptr[j] * areaR
				- fluxL->ptr[j] * areaL
			) / volume;
		}
	}<? end ?>
}
