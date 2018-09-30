<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real f = calc_f(U.alpha);
	return (<?=eqn.cons_t?>){
		.alpha = 0,
		.gamma_xx = 0,
		.a_x = U.alpha * U.K_xx * f / U.gamma_xx,
		.d_xxx = U.alpha * U.K_xx,
		.K_xx = U.alpha * U.a_x,
	};
}
<? end ?>

//used by PLM
<? for side=0,solver.dim-1 do ?>
<?=eqn.eigen_t?> eigen_forCell_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real f = calc_f(U.alpha);
	return (<?=eqn.eigen_t?>){
		.alpha = U.alpha,
		.sqrt_f_over_gamma_xx = sqrt(f / U.gamma_xx),
	};
}
<? end ?>

//used for interface eigen basis
<?=eqn.eigen_t?> eigen_forInterface(
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	real3 n
) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	real gamma_xx = .5 * (UL.gamma_xx + UR.gamma_xx);
	real f = calc_f(alpha);
	return (<?=eqn.eigen_t?>){
		.alpha = alpha, 
		.sqrt_f_over_gamma_xx = sqrt(f / gamma_xx),
	};
}

<? for side=0,solver.dim-1 do ?>
<?=eqn.waves_t?> eigen_leftTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> x,
	real3 pt
) {
	real gamma_xx_over_f = 1. / (eig.sqrt_f_over_gamma_xx * eig.sqrt_f_over_gamma_xx);
	return (<?=eqn.waves_t?>){.ptr={
		.5 * (x.ptr[2] / eig.sqrt_f_over_gamma_xx - x.ptr[4]),
		x.ptr[3] - x.ptr[2] * gamma_xx_over_f,
		.5 * (x.ptr[2] / eig.sqrt_f_over_gamma_xx + x.ptr[4]),
	}};
}

<?=eqn.cons_t?> eigen_rightTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> x,
	real3 pt
) {
	return (<?=eqn.cons_t?>){.ptr={
		0,
		0,
		(x.ptr[0] + x.ptr[2]) * eig.sqrt_f_over_gamma_xx,
		(x.ptr[0] + x.ptr[2]) / eig.sqrt_f_over_gamma_xx + x.ptr[1],
		x.ptr[2] - x.ptr[0],
	}};
}

<?=eqn.cons_t?> eigen_fluxTransform_<?=side?>(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> x,
	real3 pt 
) {
	real f_over_gamma_xx = eig.sqrt_f_over_gamma_xx * eig.sqrt_f_over_gamma_xx;
	return (<?=eqn.cons_t?>){.ptr={
		0,
		0,
		x.ptr[4] * eig.alpha * f_over_gamma_xx,
		x.ptr[4] * eig.alpha,
		x.ptr[2] * eig.alpha,
	}};
}
<? end ?>

kernel void addSource(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	real alpha = U->alpha;
	real gamma_xx = U->gamma_xx;
	real a_x = U->a_x;
	real d_xxx = U->d_xxx;
	real K_xx = U->K_xx;
	real K = K_xx / gamma_xx;	
	real f = calc_f(alpha);
	real dalpha_f = calc_dalpha_f(alpha);
	
	deriv->alpha -= alpha * alpha * f * K;
	deriv->gamma_xx -= 2. * alpha * K_xx;
	deriv->K_xx += alpha / gamma_xx * (a_x * d_xxx - K_xx * K_xx);
// terms that mysteriously disappear when you compare the linearized flux matrix terms moved to source, vs the source that Alcubierre uses in his 1997 paper
// adding these neglected terms back in make things blow up
#if 0 
	deriv->a_x += ((2. * d_xxx / gamma_xx - a_x) * f - alpha * dalpha_f * a_x) * alpha * K;
	deriv->d_xxx -= alpha * a_x * K_xx;
	deriv->K_xx -= alpha * a_x * a_x; 
#endif

	// and now for the first-order constraints
	
<? if eqn.guiVars.a_x_convCoeff.value ~= 0 then ?>
	// a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
	real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
	deriv->a_x += a_x_convCoeff * (dx_alpha / alpha - a_x);
<? end -- eqn.guiVars.a_x_convCoeff.value  ?>
	
<? if eqn.guiVars.d_xxx_convCoeff.value ~= 0 then ?>
	// d_xxx = .5 gamma_xx,x <=> d_xxx += eta (.5 gamma_xx,x - d_xxx)
	real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
	deriv->d_xxx += d_xxx_convCoeff * (.5 * dx_gamma_xx - d_xxx);
<? end -- eqn.guiVars.d_xxx_convCoeff.value  ?>

	//Kreiss-Oligar diffusion, for stability's sake?

}
