<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global <?=eqn.cons_t?>* U
) {
	real f = calc_f(U->alpha);
	real lambda = U->alpha * sqrt(f / U->gamma_xx);
	return (range_t){.min=-lambda, .max=lambda};
}
<? end ?>

//used by PLM
<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	<?=eqn.eigen_t?>* eig,
	const global <?=eqn.cons_t?>* U
) {
	real f = calc_f(U->alpha);
	eig->alpha = U->alpha;
	eig->sqrt_f_over_gamma_xx = sqrt(f / U->gamma_xx);
}
<? end ?>

//used by PLM
<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr0?> real* wave,
	<?=addr1?> const <?=eqn.eigen_t?>* eig
) {
	real lambda = eig->alpha * eig->sqrt_f_over_gamma_xx;
	wave[0] = -lambda;
	wave[1] = 0;
	wave[2] = lambda;	
}
<?		end
	end
end
?>

//used for interface eigen basis
void eigen_forSide(
	global <?=eqn.eigen_t?>* eig,
	global const <?=eqn.cons_t?>* UL,
	global const <?=eqn.cons_t?>* UR
) {
	eig->alpha = .5 * (UL->alpha + UR->alpha);
	real gamma_xx = .5 * (UL->gamma_xx + UR->gamma_xx);
	real f = calc_f(eig->alpha);
	eig->sqrt_f_over_gamma_xx = sqrt(f / gamma_xx);
}

kernel void calcEigenBasis(
	global real* waveBuf,
	global <?=eqn.eigen_t?>* eigenBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(2,1);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
	
		<?= solver.getULRCode ?>
		
		int intindex = side + dim * index;	
		
		global <?=eqn.eigen_t?>* eig = eigenBuf + intindex;
		eigen_forSide(eig, UL, UR);

		global real* wave = waveBuf + numWaves * intindex;
		eigen_calcWaves_<?=side?>_global_global(wave, eig);
	}<? end ?>
}

<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for _,addr2 in ipairs{'', 'global'} do
			for side=0,2 do
?>

void eigen_leftTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x
) {
	real gamma_xx_over_f = 1. / (eig->sqrt_f_over_gamma_xx * eig->sqrt_f_over_gamma_xx);
	y[0] = .5 * (x[2] / eig->sqrt_f_over_gamma_xx - x[4]);
	y[1] = x[3] - x[2] * gamma_xx_over_f;
	y[2] = .5 * (x[2] / eig->sqrt_f_over_gamma_xx + x[4]);
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x
) {
	y[0] = 0;
	y[1] = 0;
	y[2] = (x[0] + x[2]) * eig->sqrt_f_over_gamma_xx;
	y[3] = (x[0] + x[2]) / eig->sqrt_f_over_gamma_xx + x[1]; 
	y[4] = x[2] - x[0];
}

<?			if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const <?=eqn.eigen_t?>* eig,
	<?=addr2?> const real* x
) {
	real f_over_gamma_xx = eig->sqrt_f_over_gamma_xx * eig->sqrt_f_over_gamma_xx;
	
	y[0] = 0;
	y[1] = 0;
	y[2] = x[4] * eig->alpha * f_over_gamma_xx;
	y[3] = x[4] * eig->alpha; 
	y[4] = x[2] * eig->alpha;
}
<?				end
			end
		end
	end
end ?>

kernel void addSource(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
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
}
