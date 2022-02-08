//// MODULE_NAME: <?=updateCTU?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cell_t?> <?=cons_t?> <?=SETBOUNDS?> <?=cell_volume?> <?=solver_macros?> realparam
<? if solver.usePLM then ?>
//// MODULE_DEPENDS: <?=consLR_t?> <?=normal_t?> <?=fluxFromCons?>
<? end ?>

/*
2017 Zingale book
also 1990 Collela "Multidimensional Upwind Methods for Hyperbolic Conservation Laws"
also Mara
using the interface flux
update the LR states by integrating dt/2 times the interface flux
*/
kernel void <?=updateCTU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cell_t?> const * const cellBuf,
	global <?=solver.getULRArg?>,
	global <?=cons_t?> const * const fluxBuf,
	realparam const dt
) {
	<?=SETBOUNDS?>(0,1);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
<? if eqn.weightFluxByGridVolume then ?>
	real const volume = cell_volume(solver, x);
<? else ?>
	real const volume = 1.<? for i=0,solver.dim-1 do ?> * solver->grid_dx.s<?=i?><? end ?>;
<? end ?>

	//these are the inter-cell fluxes, computed from the previously calculated calcLR states
<? for side=0,solver.dim-1 do
?>	int const indexIntL<?=side?> = <?=side?> + dim * index;
	global <?=cons_t?> const * const fluxL<?=side?> = fluxBuf + indexIntL<?=side?>;
	
	int const indexIntR<?=side?> = indexIntL<?=side?> + dim * solver->stepsize.s<?=side?>;
	global <?=cons_t?> const * const fluxR<?=side?> = fluxBuf + indexIntR<?=side?>;
<? end
?>

	// calculate and save the cell side areas for fv update
<? 
for side=0,solver.dim-1 do
?>
<?	if eqn.weightFluxByGridVolume then 
?>	real3 xIntL<?=side?> = x;
	xIntL<?=side?>.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real const volume_intL<?=side?> = cell_volume(solver, xIntL<?=side?>);
	real const areaL<?=side?> = volume_intL<?=side?> / solver->grid_dx.s<?=side?>;

	real3 xIntR<?=side?> = x;
	xIntR<?=side?>.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	real const volume_intR<?=side?> = cell_volume(solver, xIntR<?=side?>);
	real const areaR<?=side?> = volume_intR<?=side?> / solver->grid_dx.s<?=side?>;
<? 
	else 
?>	real const areaL<?=side?> = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
	real const areaR<?=side?> = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
<? 
	end
end
?>

	<?
for side=0,solver.dim-1 do
	?>{
<?
	if solver.getULRArg == consLR_t.."* ULRBuf" then
?>		
		//these are the cell left and right edge states (from calcLR)'s fluxes
		int const indexForSide = <?=side?> + dim * index;
		global <?=consLR_t?> * const ULR = ULRBuf + indexForSide;
		
		//calc flux from ULR state here before modifying ULR
		<?=cons_t?> fluxCellL, fluxCellR;
		<?=normal_t?> n = normal_forSide<?=side?>(x);
		<?=fluxFromCons?>(&fluxCellL, solver, &(ULR->L), cell, n);
		<?=fluxFromCons?>(&fluxCellR, solver, &(ULR->R), cell, n);
<? 
	elseif solver.getULRArg == cons_t.."* UBuf" then
?>		global <?=cons_t?> * const U = UBuf + index;
<? 
	end
?>

		for (int j = 0; j < numIntStates; ++j) {
				
			<?
	for side2=0,solver.dim-1 do
		-- for non-PLM, i.e. constant slope, the flux dif across the cell is zero, so just skip this part
		if not (
			side == side2
			and solver.getULRArg == cons_t.."* UBuf" 
		) then
			?>{
<?			local fluxL, fluxR	-- code for source of fluxL/fluxR.  for side==side2 this is F(ULR), for otherwise it is intercell flux
			if side2 == side then
				-- by here we must be using ULRBuf ...
				fluxL = "fluxCellL"
				fluxR = "fluxCellR"
			else
				fluxL = "(*fluxL"..side2..")"
				fluxR = "(*fluxR"..side2..")"
			end

			if solver.getULRArg == consLR_t.."* ULRBuf" then
				local dF_dx = "dF_dx"..side2
?>			
			
				real const <?=dF_dx?> = (
					<?=fluxR?>.ptr[j] * areaR<?=side2?>
					- <?=fluxL?>.ptr[j] * areaL<?=side2?>
				) / volume;

				ULR->L.ptr[j] -= .5 * dt * <?=dF_dx?>;
				ULR->R.ptr[j] -= .5 * dt * <?=dF_dx?>;
			
			
<? 			elseif solver.getULRArg == cons_t.."* UBuf" then
?>			
				
				real const <?=dF_dx?> = (
					<?=fluxR?>.ptr[j] * areaR<?=side2?>
					- <?=fluxL?>.ptr[j] * areaL<?=side2?>
				) / volume;
				U->ptr[j] -= .5 * dt * <?=dF_dx?>;
			

<? 			else 
				error("can't handle getULRArg "..tostring(solver.getULRArg)) 
			end
?>			}<?
		end
	end ?>
		}
	}<?
end ?>
}

<?
--[[
TODO alternatively, if you want to do the CTU in primitive space instead of conserved: 
(a) record/recalculate the cell prim slope dW/dx_k 
(b) HERE only adjust the center cell cons value U
(c) convert adjusted cons U to prim W, then re-project slope to edges WLR = W +- 1/2 dW/dx_k
(d) then recalculate ULR at edges based on projected WLR
--]]
?>
