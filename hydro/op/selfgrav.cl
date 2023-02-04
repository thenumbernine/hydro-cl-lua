//// MODULE_NAME: <?=SelfGrav?>
//// MODULE_DEPENDS: <?=Solver?>

namespace <?=Solver?> {

struct SelfGrav {
	/*
	This is going to give back phi_,i ... in grid coordinates ... which might not be vector coordinates
	and if they're not?  convert
	*/
	template<int dim_>
	static inline real3 calcGravityAccel(
		constant Solver const & solver,
		global Cons const * private const U,
		real3 const pt
	) {
		real3 accel_g;

		for (int i = 0; i < dim_; ++i) {
			// m/s^2
			// TODO grid coordinate influence?
			accel_g[i] = (
				U[solver.stepsize[i]].<?=op.potentialField?>
				- U[-solver.stepsize[i]].<?=op.potentialField?>
			) / (2. * solver.grid_dx[i]);
		}

<? if coord.vectorComponent == "cartesian" then ?>
//// MODULE_DEPENDS: <?=cartesianFromCoord?>
		accel_g = coord_cartesianFromCoord(accel_g, pt);
<? elseif coord.vectorComponent == "anholonomic" then ?>
//// MODULE_DEPENDS: <?=coord_dxs?>
		accel_g /= coord_dx_vec(solver, pt);
<?
elseif coord.vectorComponent == "holonomic" then
		-- nothing
else 
	error "here" 
end
?>
		return accel_g;
	}
};

}


//// MODULE_NAME: <?=calcGravityDeriv?>
//// MODULE_DEPENDS: <?=SelfGrav?>
//// MODULE_DEPENDS: units realparam

kernel void <?=calcGravityDeriv?>(
	constant <?=solver_t?> const * const psolver,
	global <?=cons_t?> * const derivBuffer,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	auto const & solver = *psolver;
	<?=SETBOUNDS?>(solver.numGhost, solver.numGhost);
	real3 const pt = cellBuf[index].pos;
	auto * const deriv = derivBuffer + index;
	auto const * const U = UBuf + index;
	real3 accel_g = <?=Solver?>::SelfGrav::calcGravityAccel<dim>(solver, U, pt);
	deriv->m -= accel_g * U->rho;			// kg/(m^2 s) = kg/m^3 * m/s^2
	deriv->ETotal -= dot(U->m, accel_g);	// kg/(m s^2) = (kg m^2 / s) * m/s^2
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

