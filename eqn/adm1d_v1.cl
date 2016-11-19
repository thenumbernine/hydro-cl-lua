<? for side=0,solver.dim-1 do ?>
range_t calcCellMinMaxEigenvalues_<?=side?>(
	const global cons_t* U
) {
	real f = calc_f(U->alpha);
	real lambda = U->alpha * sqrt(f / U->gamma_xx);
	return (range_t){.min=-lambda, .max=lambda};
}
<? end ?>

//used by PLM
<? for side=0,solver.dim-1 do ?>
void eigen_forCell_<?=side?>(
	eigen_t* eig,
	const global cons_t* U
) {
	eig->alpha = U->alpha;
	eig->gamma_xx = U->gamma_xx;
	eig->f = calc_f(U->alpha);
}
<? end ?>

//used by PLM
<?
for _,addr0 in ipairs{'', 'global'} do
	for _,addr1 in ipairs{'', 'global'} do
		for side=0,solver.dim-1 do
?>
void eigen_calcWaves_<?=side?>_<?=addr0?>_<?=addr1?>(
	<?=addr1?> real* wave,
	const <?=addr0?> eigen_t* eig
) {
	real lambda = eig->alpha * sqrt(eig->f / eig->gamma_xx);
	wave[0] = -lambda;
	wave[1] = 0;
	wave[2] = lambda;
}
<?		end
	end
end
?>

void eigen_forSide(
	global eigen_t* eig,
	const global cons_t* UL,
	const global cons_t* UR
) {
	eig->alpha = .5 * (UL->alpha + UR->alpha);
	eig->gamma_xx = .5 * (UL->gamma_xx + UR->gamma_xx);
	eig->f = calc_f(eig->alpha);
}

kernel void calcEigenBasis(
	global real* waveBuf,
	global eigen_t* eigenBuf,
	<?= solver.getULRArg ?>
) {
	SETBOUNDS(2,1);
	int indexR = index;
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - stepsize[side];
	
		<?= solver.getULRCode ?>
		
		int intindex = side + dim * index;
		global eigen_t* eig = eigenBuf + intindex;
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
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	real sqrt_f = sqrt(eig->f);
	y[0] = (x[2] / eig->f - x[4] / sqrt_f) / 2.;
	y[1] = -2. * x[2] / eig->f + x[3];
	y[2] = (x[2] / eig->f + x[4] / sqrt_f) / 2.;
}

void eigen_rightTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	y[0] = 0;
	y[1] = 0;
	y[2] = (x[0] + x[2]) * eig->f;
	y[3] = 2. * (x[0] + x[2]) + x[1];
	y[4] = sqrt(eig->f) * (x[2] - x[0]);
}

<?				if solver.checkFluxError then ?>
void eigen_fluxTransform_<?=side?>_<?=addr0?>_<?=addr1?>_<?=addr2?>(
	<?=addr0?> real* y,
	<?=addr1?> const eigen_t* eig,
	<?=addr2?> const real* x
) {
	real alpha_sqrt_gamma_xx = eig->alpha / sqrt(eig->gamma_xx);
	y[0] = 0;
	y[1] = 0;
	y[2] = x[4] * eig->f * alpha_sqrt_gamma_xx;
	y[3] = x[4] * 2. * alpha_sqrt_gamma_xx;
	y[4] = x[2] * alpha_sqrt_gamma_xx;
}
<? 				end
			end
		end
	end
end
?>

kernel void addSource(
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(2,2);
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;
	
	real alpha = U->alpha;
	real gamma_xx = U->gamma_xx;
	real a_x = U->a_x;
	real D_g = U->D_g;
	real KTilde_xx = U->KTilde_xx;
	
	real sqrt_gamma_xx = sqrt(gamma_xx);
	real K_xx = KTilde_xx / sqrt_gamma_xx;
	
	real f = calc_f(alpha);
	real dalpha_f = calc_dalpha_f(alpha);
	
	deriv->alpha -= alpha * alpha * f * K_xx / gamma_xx;
	deriv->gamma_xx -= 2. * alpha * K_xx;
	deriv->a_x -= alpha * K_xx * (f * (.5 * D_g - a_x) - a_x * alpha * dalpha_f);
	deriv->D_g -= 2. * alpha * K_xx * (.5 * D_g - a_x);
	deriv->KTilde_xx -= alpha * a_x / sqrt_gamma_xx * (.5 * D_g - a_x);
}
