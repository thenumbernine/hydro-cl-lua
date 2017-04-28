kernel void calcDeriv(
	global <?=eqn.cons_t?>* derivBuf,
	const global <?=eqn.cons_t?>* UBuf,
) {
	SETBOUNDS(2,2);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
		
	const global <?=eqn.cons_t?>* Up[dim];
	const global <?=eqn.cons_t?>* Um[dim];
	for (int i = 0; i < dim; ++i) {
		Up[i] = U + index + stepsize[i];
		Um[i] = U + index - stepsize[i];
	}

	real exp_4phi = exp(4. * U->phi);
	real exp_neg4phi = 1. / exp_4phi;

	//gamma_ij = exp(4 phi) gammaTilde_ij
	sym3 gamma = sym3_scale(U->gammaTilde, exp_4phi);

	//gammaTilde^ij = inv gammaTilde_ij
	sym3 gammaTildeInv = sym3_inv(U->gammaTilde);

	//gamma^ij = inv gamma_ij	
	//sym3 gammaU = sym3_inv(gamma);
	//gamma^ij = exp(-4 phi) gammaTilde^ij
	sym3 gammaInv = sym3_scale(gammaTildeInv, exp_neg4phi);

	//tr K = K_ij gamma^ij
	real tr_K = sym3_dot(U->K, gammaU);

	//Q = f(alpha) tr K
	real Q = calc_f(alpha) * tr_K;
	
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -alpha * alpha * Q + real3_dot(dx_alpha, beta);

	//manuall update elsewhere?
	//deriv->beta += _real3(0,0,0);

	real3 dx_phi;			//phi,i = dx_phi[i]
	real3 dx_beta[3];		//beta^i_,j = dx_beta[j][i]
	real tr_dx_beta = 0.;	//beta^i_,i
	real3 dx_gammaTilde[3];
<? for i=0,dim-1 do
?>	dx_phi[<?=i?>] = (Up[<?=i?>]->phi - Um[<?=i?>]->phi) / (2 * dx<?=i?>_at(i);
	dx_beta[<?=i?>] = real3_scale( real3_sub(Up[<?=i?>]->beta, Um[<?=i?>]->beta), 1. / (2. * dx<?=i?>_at(i) ));
	dx_gammaTilde[<?=i?>] = sym3_scale( sym3_sub(Up[<?=i?>]->gammaTilde, Um[<?=i?>]->gammaTilde, 1. / (2. * dx<?=i?>_at(i) ));
	tr_dx_beta += dx_beta[<?=i?>][<?=i?>];
<? end
for i=dim,2 do
?>	dx_phi[<?=i?>] = 0.;
	dx_beta[<?=i?>] = _real3(0,0,0);
	dx_gammaTilde[<?=i?>] = (sym3){0,0,0,0,0,0};
<? end
?>
	//B&S 11.50
	//phi,t = -1/6 alpha tr K + beta^i phi,i + 1/6 beta^i_,i
	deriv->phi += -U->alpha * tr_K / 6. + real3_dot(U->beta, dx_phi) + tr_dx_beta / 6.;

	//B&S 11.51
	//gammaTilde_ij,t = -2 alpha ATilde_ij + beta^k gammaTilde_ij,k + gammaTilde_ik beta^k_,j + gammaTilde_kj beta^k_,i - 2/3 gammaTilde_ij beta^k_,k
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->gammaTilde.<?=xij?> += -2 * U->alpha * U->ATilde.<?=xij?>	//-2 alpha ATilde_ij 
<? 	for k,xk in ipairs(xNames) do
		local ik = from3x3to6(i,k)
		local jk = from3x3to6(j,k)
?>		+ dx_gammaTilde[<?=k?>].<?=xij?> * U->beta.<?=xk?>		//+ beta^k gammaTilde_ij,k 
		+ U->gammaTilde.s<?=ik?> * dx_beta[<?=j+1?>].<?=xk?> 	//+ gammaTilde_ik beta^k_,j
		+ U->gammaTilde.s<?=jk?> * dx_beta[<?=j+1?>].<?=xi?>	//+ gammaTilde_kj beta^k_,i 
<? 	end
?>		- 2./3. * U->gammaTilde.s[<?=ij?>] * tr_dx_beta;		//- 2/3 gammaTilde_ij beta^k_,k
<?
end
?>
}
