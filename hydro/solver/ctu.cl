<?
local table = require "ext.table"
?>
//// MODULE_NAME: <?=updateCTU?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cell_t?> <?=cons_t?> <?=SETBOUNDS?> <?=cell_volume?> <?=solver_macros?> realparam
<? if solver.usePLM then ?>
//// MODULE_DEPENDS: <?=consLR_t?> <?=normal_t?> <?=fluxFromCons?>
<? end ?>

/*
2017 Zingale "Introduction to Computational Astrophysics"
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
	real const invVolume = 1. / volume;

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
	if eqn.weightFluxByGridVolume then
?>
	real3 xIntL<?=side?> = x;
	xIntL<?=side?>.s<?=side?> -= .5 * solver->grid_dx.s<?=side?>;
	real const volume_intL<?=side?> = cell_volume(solver, xIntL<?=side?>);
	real const areaL<?=side?> = volume_intL<?=side?> / solver->grid_dx.s<?=side?>;

	real3 xIntR<?=side?> = x;
	xIntR<?=side?>.s<?=side?> += .5 * solver->grid_dx.s<?=side?>;
	real const volume_intR<?=side?> = cell_volume(solver, xIntR<?=side?>);
	real const areaR<?=side?> = volume_intR<?=side?> / solver->grid_dx.s<?=side?>;
<?
	else
?>
	real const areaL<?=side?> = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
	real const areaR<?=side?> = 1.<? for i=0,solver.dim-1 do if i ~= side then ?> * solver->grid_dx.s<?=i?><? end end ?>;
<?
	end
end
?>

	<?
for side=0,solver.dim-1 do
	?>{		// correct along direction <?=side?>
<?
	local updateFields = table()
	if solver.getULRBufType == eqn.symbols.consLR_t then
?>
		//these are the cell left and right edge states (from calcLR)'s fluxes
		int const indexForSide = <?=side?> + dim * index;
		global <?=consLR_t?> * const ULR = ULRBuf + indexForSide;
		
		//calc flux from ULR state here before modifying ULR
		<?=cons_t?> fluxCellL, fluxCellR;
		<?=normal_t?> n = normal_forSide<?=side?>(x);
		<?=fluxFromCons?>(&fluxCellL, solver, &ULR->L, cell, n);
		<?=fluxFromCons?>(&fluxCellR, solver, &ULR->R, cell, n);
<?
		updateFields:insert"ULR->L"
		updateFields:insert"ULR->R"
	elseif solver.getULRBufType == eqn.symbols.cons_t then
		-- imporatnt detail about running non-PLM with CTU...
		-- (you shouldn't but who am I to stop you)
		-- the PLM version will not update the UBuf ... so after correcting the flux, the subsequent UBuf integration step will be applied to the original UBuf
		-- however the non-PLM version *will* update the UBuf, so the new UBuf won't only contribute to the adjusted flux, but it will also contribute to the next UBuf integration, and will probably add to the error
?>		global <?=cons_t?> * const U = UBuf + index;
<?
		updateFields:insert"(*U)"
	else
		error("can't handle getULRBufType "..tostring(solver.getULRBufType))
	end
?>

		for (int j = 0; j < numIntStates; ++j) {
<?	for side2=0,solver.dim-1 do
		-- for non-PLM, i.e. constant slope, since the state is constant across the sell, the flux dif across the cell is zero, so just skip this side
		if not (
			side == side2
			and solver.getULRBufType == eqn.symbols.cons_t
		) then
			local fluxL, fluxR	-- code for source of fluxL/fluxR.  for side==side2 this is F(ULR), for otherwise it is intercell flux
			if side2 == side then
				-- by here we must be using ULRBuf ...
				fluxL = "fluxCellL"
				fluxR = "fluxCellR"
			else
				fluxL = "(*fluxL"..side2..")"
				fluxR = "(*fluxR"..side2..")"
			end

			local dF_dx = "dF_dx"..side2
?>
			real const <?=dF_dx?> = (
				<?=fluxR?>.ptr[j] * areaR<?=side2?>
				- <?=fluxL?>.ptr[j] * areaL<?=side2?>
			) * invVolume;
<?
			for _,updateField in ipairs(updateFields) do
?>			<?=updateField?>.ptr[j] -= .5 * dt * <?=dF_dx?>;
<?
			end
		end
	end
?>		}
	}<?
end ?>
}

<?
--[[
TODO alternatively, if you want to do the CTU in primitive space instead of conserved:
(a) when adjusting U, read from cell-centered U, write to ULR ... the L and R states will match, but they will differ per-dir
(b) re-run calcLR, but this time accept ULRBuf as the input instead of UBuf 
but to do this I'd have to pass in the UBuf here even when writing to ULRBuf
--]]
?>
