//// MODULE_NAME: <?=calcGravityAccel?>

void <?=calcGravityAccel?>(
	/* output: */
	real * const W,
	real3 * const u,
	real * const dW_dt,
	real3 * const du_dt,
	/* input: */
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U
) {
<? 
for i,xi in ipairs(xNames) do
	if i <= solver.dim then
?>	(du_dt)-><?=xi?> = (
		U[solver->stepsize.<?=xi?>].<?=op.potentialField?> 
		- U[-solver->stepsize.<?=xi?>].<?=op.potentialField?>
	) / (2. * solver->grid_dx.<?=xi?>);
<?	else
?>	(du_dt)-><?=xi?> = 0;
<?	end
end
?>

	real const Phi = U-><?=op.potentialField?>;

	//u = W v
	*(W) = U->D / U->rho;
	*(u) = real3_real_mul(U->v, *(W));

	real const uSq = real3_dot(*(u), *(u)) / (1. - 2. * Phi);
	real3 const dv_dt = real3_add(
		real3_real_mul(*(du_dt), 1. / *(W)),
		real3_real_mul(*(u), real3_dot(*(u), *(du_dt)) / (*(W) * (1. - 2. * Phi) * (1. + uSq)))
	);

	//W_,t = W^3 (v^i_,t v_i + v^i v^j gamma_ij,t / 2)
	//W_,t = W^3 (v^i_,t v^i) / (1 - 2 Phi)
	*(dW_dt) = *(W) * *(W) * *(W) * real3_dot(U->v, dv_dt) / (1. - 2. * Phi);
}

//// MODULE_NAME: <?=calcGravityDeriv?>
//// MODULE_DEPENDS: units

kernel void <?=calcGravityDeriv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuffer,
	global <?=cons_t?> const * const UBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const deriv = derivBuffer + index;
	global <?=cons_t?> const * const U = UBuf + index;

	real W, dW_dt;
	real3 u, du_dt;
	<?=calcGravityAccel?>(
		&W,
		&u,
		&dW_dt,
		&du_dt,
		solver,
		U
	);

	real const h = 1. + solver->heatCapacityRatio * U->eInt;

	//why am I integrating negative again?
	//why does "Hydrodynamics II" say to integrate negative for the Euler equations?

	//D = W rho
	//D,t = W,t rho
	deriv->D -= dW_dt * U->rho;

	//S = rho h W^2 v = rho h W u
	//assuming rho and h are constant ... 
	//S,t = rho h (W,t u + W u,t)
	deriv->S = real3_sub(deriv->S,
		real3_add(
			real3_real_mul(u, U->rho * h * dW_dt),
			real3_real_mul(du_dt, U->rho * h * W)
		)
	);
	
	//tau = rho h W^2 - p - rho W
	//tau,t = rho h (2 W W,t) - rho W,t
	//tau,t = rho W,t (2 h W - 1)
	deriv->tau -= U->rho * dW_dt * (2. * h * W - 1.);
}

