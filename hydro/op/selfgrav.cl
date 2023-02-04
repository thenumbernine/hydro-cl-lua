//// MODULE_NAME: <?=calcGravityAccel?>

/*
This is going to give back phi_,i ... in grid coordinates ... which might not be vector coordinates
and if they're not?  convert
*/
template<int dim_>
static inline real3 <?=calcGravityAccel?>(
	constant <?=solver_t?> const & solver,
	global <?=cons_t?> const * private const U,
	real3 const pt
) {
	real3 accel_g;

	for (int side = 0; side < dim_; ++side) {
		// m/s^2
		// TODO grid coordinate influence?
		accel_g[side] = (
			U[solver.stepsize[side]].<?=op.potentialField?>
			- U[-solver.stepsize[side]].<?=op.potentialField?>
		) / (2. * solver.grid_dx[side]);
	}

<? if coord.vectorComponent == "cartesian" then ?>
//// MODULE_DEPENDS: <?=cartesianFromCoord?>
	accel_g = coord_cartesianFromCoord(accel_g, pt);
<? elseif coord.vectorComponent == "anholonomic" then
	for i=0,solver.dim-1 do
?>
//// MODULE_DEPENDS: <?=coord_dxs?>
	accel_g[<?=i?>] /= coord_dxs<<?=i?>>(pt);
<?
	end
elseif coord.vectorComponent == "holonomic" then
	-- nothing
else 
	error "here" 
end
?>
	return accel_g;
}

//// MODULE_NAME: <?=calcGravityDeriv?>
//// MODULE_DEPENDS: units realparam <?=calcGravityAccel?>

kernel void <?=calcGravityDeriv?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const derivBuffer,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(solver.numGhost, solver.numGhost);
	real3 const pt = cellBuf[index].pos;

	global <?=cons_t?> * const deriv = derivBuffer + index;
	global <?=cons_t?> const * const U = UBuf + index;

	real3 accel_g = <?=calcGravityAccel?><dim>(solver, U, pt);

	// kg/(m^2 s) = kg/m^3 * m/s^2
	deriv->m -= accel_g * U->rho;

	// kg/(m s^2) = (kg m^2 / s) * m/s^2
	deriv->ETotal -= dot(U->m, accel_g);
}

//// MODULE_NAME: <?=copyPotentialToReduce?>

//TODO just use the display var kernels
kernel void <?=copyPotentialToReduce?>(
	constant <?=solver_t?> const * const psolver,
	global real* reduceBuf,
	global const <?=cons_t?>* UBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}

//// MODULE_NAME: <?=offsetPotential?>

//keep potential energy negative
kernel void <?=offsetPotential?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?>* UBuf,
	realparam ePotMax
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> -= ePotMax;
}

