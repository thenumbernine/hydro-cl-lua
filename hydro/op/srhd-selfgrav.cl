//// MODULE_NAME: notes

/*
This is all using a Cartesian background iirc
TODO make it work for any grid metric
TODO TODO this looks like it's done for u^i when my srhd variables say they use u_i ... so pick your variance


u_i = W v_i

u_i,t = W_,t v_i + W v_i,t
W = 1/sqrt(1 - v^2)
in non-Cartesian grid metrics, v^2 = v^i v^j γ_ij
ADM formalism: γ_ij = g_ij
W_,t = -1/2 1/(1 - v^2)^(3/2) * (-2 v^i v_i,t)
W_,t = W^3 (v^i v_i,t)
W_,t v_i + W v_i,t = W^3 v_i (v^j v_j,t) + W v_i,t
u_i,t = W (W^2 v_i (v^j v_j,t) + v_i,t)
u_i,t = W (u_i u^j + δ_i^j) v_j,t

(u outer u + I) (I - α u outer u)
= u outer u - α u^2 u outer u + I - α u outer u
= I + u outer u (1 - α (1 + u^2))
is the inverse for 1 - α (1 + u^2) = 0
 i.e. α = 1 / (1 + u^2)
so the inverse of (I + u outer u) is (I - u outer u / (1 + u^2))

v,t = 1/W (I - u outer u / (1 + u^2)) u,t
v,t = (u,t - u (u dot u,t) / (1 + u^2)) / W


in general relativity
using a scalar metric g_ab = η_ab - 2 Φ δ_ab so that the g_at = 0
NOTICE this assumes a Cartesian background.
TODO rewrite Newtonian limit of GR for arbitrary spatial background.

W_,t = -1/2 1/(1 - v^m v^n γ_mn)^(3/2) * -(2 v^i_,t v_i + v^i v^j γ_ij,t)
W_,t = W^3 (v^i_,t v_i + v^i v^j γ_ij,t / 2)
W_,t v^j + W v^j_,t = v^j W^3 (v^i_,t v_i + v^i v^k γ_ik,t / 2) + W v^j_,t
u^j_,t = W (W^2 v^j v_i + δ^j_i) v^i_,t + v^j W^3 v^i v^k γ_ik,t / 2
	v_i = γ_ij v^j by v_a = g_ab v^b = (g_at v^t + g_aj v^j) 
u^j_,t - u^i u^j u^k γ_ik,t / 2 = (u_i u^j + δ_i^j) W v^i_,t 

δ_i^k = (δ_i^j + u_i u^j) (δ_j^k - u_j u^k / (1 + u^m u_m))
= δ_i^k + u_i u^k - u_i u^k / (1 + u^m u_m) - (u_i u^k (u^j u_j) / (1 + u^m u_m)
= δ_i^k + u_i u^k - u_i u^k (1 + u_j u^j) / (1 + u_m u^m)
= δ_i^k

(u^j_,t - u^i u^j u^l γ_il,t / 2) (δ_j^k - u_j u^k / (1 + u^m u_m))
	= (δ_j^k - u_j u^k / (1 + u^m u_m)) (u_i u^j + δ_i^j) W v^i_,t 
v^k_,t = (δ^k_j - u^k u_j / (1 + u^m u_m)) (u^j_,t - u^i u^j u^l γ_il,t / 2) / W
v^i_,t = (δ^i_j - u^i u_j / (1 + u^m u_m)) (u^j_,t - u^j u^k u^l γ_jk,t / 2) / W
	ignore γ_ij,t for now <-> grid unchanging in time
v^i_,t = (δ^i_j - u^i u_j / (1 + u^m u_m)) u^j_,t / W
v^i_,t = u^i_,t - u^i u_j / (1 + u^m u_m) u^j_,t / W
	using the scalar metric
v^i_,t = u^i_,t / W - u^i u_j u^j_,t / (W (1 - 2 Φ) (1 + u^m u_m))
*/

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

	//W_,t = W^3 (v^i_,t v_i + v^i v^j γ_ij,t / 2)
	//W_,t = W^3 (v^i_,t v^i) / (1 - 2 Φ)
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

	//D = W ρ
	//D_,t = W_,t ρ
	deriv->D -= dW_dt * U->rho;

	//S_i = ρ h W^2 v_i = ρ h W u_i
	//assuming ρ and h are constant ... 
	//S_i,t = ρ h (W_,t u_i + W u_i,t)
	deriv->S = real3_sub(deriv->S,
		real3_add(
			real3_real_mul(u, U->rho * h * dW_dt),
			real3_real_mul(du_dt, U->rho * h * W)
		)
	);

	//τ = ρ h W^2 - p - ρ W
	//τ,t = ρ h (2 W W_,t) - ρ W_,t
	//τ,t = ρ W_,t (2 h W - 1)
	deriv->tau -= U->rho * dW_dt * (2. * h * W - 1.);
}

