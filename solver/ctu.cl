/*
using the interface flux
update the LR states by integrating dt/2 times the interface flux
*/
kernel void updateCTU(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.consLR_t?>* ULRBuf,
	const global <?=eqn.cons_t?>* fluxBuf,
	real dt
) {
	SETBOUNDS(0,1);
	real3 x = cell_x(i);
	real sqrt_det_g = sqrt_det_g_grid(x);

	<? 
for side=0,solver.dim-1 do 
	?>{
		const int side = <?=side?>;
		
		int indexIntL = side + dim * index;
		const global <?=eqn.cons_t?>* fluxL = fluxBuf + indexIntL;
		
		int indexIntR = indexIntL + dim * solver->stepsize.s<?=side?>; 
		const global <?=eqn.cons_t?>* fluxR = fluxBuf + indexIntR;
			
		real3 xIntL = x;
		xIntL.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
		real sqrt_det_g_intL = sqrt_det_g_grid(xIntL);
		real areaL = sqrt_det_g_intL / solver->grid_dx.s<?=side?>;
	
		real3 xIntR = x;
		xIntR.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
		real sqrt_det_g_intR = sqrt_det_g_grid(xIntR);
		real areaR = sqrt_det_g_intR / solver->grid_dx.s<?=side?>;
	
		<?
	for side2=0,solver.dim-1 do
		if side2 ~= side then
		?>{
			int indexForSide2 = side + dim * index;
			global <?=eqn.consLR_t?>* ULR = ULRBuf + indexForSide2;
			
			for (int j = 0; j < numIntStates; ++j) {
				real dF_dx = (
					fluxR->ptr[j] * areaR
					- fluxL->ptr[j] * areaL
				) / sqrt_det_g;

				ULR->L.ptr[j] -= .5 * dt * dF_dx;
				ULR->R.ptr[j] -= .5 * dt * dF_dx;
			}
		}<? 
		end
	end ?>
	}<? 
end ?>
}
