/*
This routine is pretty standard.

Also, to generalize this to meshes, how about 
storing the # of interfaces in the solver:
	numFaces = numCells * dim
and then cycle across that,
and - for grid-based solvers - deduce which interface to use based on the global index

Or maybe we can just make this routine for grids, and another for meshes, 
and both can call eigen_forSide...
...of course eigen_forSide is axis-aligned, so the mesh one will need to call a linear combination of all sides based on the normal ...
*/

<?
local solver = eqn.solver
?>

<? if require 'solver.gridsolver'.is(solver) then ?>

kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,		//[numCells][dim]
	<?=solver.getULRArg?>
) {
	SETBOUNDS(numGhost,numGhost-1);
	real3 x = cell_x(i);
	int indexR = index;
	
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize.s<?=side?>;
		
		<?=solver:getULRCode()?>
		
		real3 xInt = x;
		xInt.s<?=side?> -= .5 * grid_dx<?=side?>;
		
		int indexInt = side + dim * index;	
		eigenBuf[indexInt] = eigen_forSide_<?=side?>(*UL, *UR, xInt);
	}<? end ?>
}

<? else -- mesh solver ?>

kernel void calcEigenBasis(
	global <?=eqn.eigen_t?>* eigenBuf,	//[numInterfaces]
	global <?=eqn.cons_t?>* UBuf,		//[numCells]
	const global cell_t* cells,			//[numCells]
	const global iface_t* ifaces		//[numInterfaces]
) {
	int ifaceIndex = get_global_id(0);
	if (ifaceIndex >= get_global_size(0)) return;

	const global iface_t* iface = ifaces + ifaceIndex;

	const global <?=eqn.cons_t?>* UL = UBuf + iface->cellIndex[0];
	const global <?=eqn.cons_t?>* UR = UBuf + iface->cellIndex[1];

	eigenBuf[ifaceIndex] = eigen_forInterface(*UL, *UR, iface->x);
}

<? end -- mesh vs grid solver ?>
