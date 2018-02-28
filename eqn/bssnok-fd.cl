/*
Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer" 2010
Alcubierre "Introduction to Numerical Relativity" 2008
*/

<? local calcConstraints = true ?>

/*
TF(K_ij) = K_ij - 1/3 gamma_ij gamma^kl K_kl

tr(A_ij)
= tr(K_ij - 1/3 gamma_ij K)
= gamma^ij K_ij - 1/3 gamma^ij gamma_ij K
= K - 1/3 3 K
= 0

tr(ATilde_ij) = exp(-4 phi) tr(A_ij) 
= exp(-4 phi) * 0
= 0

TFBar(K_ij) = K_ij - 1/3 gammaBar_ij gammaBar^kl K_kl 
	= K_ij - 1/3 gamma_ij gamma^kl K_kl
	= TF(K_ij)
*/
sym3 tracefree(sym3 A_ll, sym3 gamma_ll, sym3 gamma_uu) {
	real tr_A = sym3_dot(A_ll, gamma_uu);
	return sym3_sub(A_ll, sym3_scale(gamma_ll, tr_A / 3.));
}

kernel void constrainU(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	/*
	det(gammaBar_ij) 
	= det(gamma^-1/3 gamma_ij)
	= gamma^-1 gamma
	= 1
	
	det(a * 1/det(a)^(1/n) )
	det(a) * 1/det(a)^(1/n)^n
	det(a) * 1/det(a)
	1
	*/
<? if eqn.guiVars.constrain_det_gammaBar_ll.value then ?>
	real det_gammaBar = sym3_det(U->gammaBar_ll);
	real _1_cbrt_det_gammaBar = 1./cbrt(det_gammaBar);
<? 	for ij,xij in ipairs(symNames) do
?>	U->gammaBar_ll.<?=xij?> *= _1_cbrt_det_gammaBar;
<? 	end
end ?>

	//here we refresh gammaBar_uu
	U->gammaBar_uu = sym3_inv(U->gammaBar_ll, 1.);
	
	//in Buchman's paper it says he doesn't do this
	//likewise in my own experiences, this can tend A to grow out of control 
<? if eqn.guiVars.constrain_tr_ATilde_ll.value then ?>
	U->ATilde_ll = tracefree(U->ATilde_ll, U->gammaBar_ll, U->gammaBar_uu);
<? end ?>
}

kernel void calcDeriv(
	global <?=eqn.cons_t?>* derivBuf,
	<?=calcConstraints and '' or 'const '?>global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* deriv = derivBuf + index;

	<?=calcConstraints and '' or 'const '?>global <?=eqn.cons_t?>* U = UBuf + index;

<?=makePartial('alpha', 'real')?>		//partial_alpha_l[i] := alpha_,i
<? if eqn.useChi then ?>
<?=makePartial('chi', 'real')?>			//partial_chi_l[i] := chi_,i 
<?=makePartial('gammaBar_uu', 'sym3')?>	//partial_gammaBar_uul[k].ij = gamma^ij_,k
<? else ?>
<?=makePartial('phi', 'real')?>			//partial_phi_l[i] := phi_,i 
<? end ?>
<?=makePartial('K', 'real')	?>			//partial_K_l[i] := K,i
<?=makePartial('beta_u', 'real3')?>		//partial_beta_ul[j].i := beta^i_,j
<?=makePartial('connBar_u', 'real3')?>	//partial_connBar_ul[j].i := connBar^i_,j
<?=makePartial('gammaBar_ll', 'sym3')?>	//partial_gammaBar[k].ij := gammaBar_ij,k
<?=makePartial('ATilde_ll', 'sym3')?>	//partial_ATilde_lll[k].ij = ATilde_ij,k
<? if eqn.useHypGammaDriver then ?>
<?=makePartial('B_u', 'real3')?>
<? end ?>
	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

<?=makePartial2('alpha', 'real')?>			//partial2_alpha_ll.ij := alpha_,ij
<? if eqn.useChi then ?>
<?=makePartial2('chi', 'real')?>			//partial2_chi_ll.ij := chi_,ij
<? else ?>
<?=makePartial2('phi', 'real')?>			//partial2_phi_ll.ij := phi_,ij
<? end ?>
<?=makePartial2('gammaBar_ll', 'sym3')?>	//partial2_gammaBar_llll[kl].ij = gammaBar_ij,kl
<?=makePartial2('beta_u', 'real3')?>		//partial2_beta_ull[jk].i = beta^i_,jk

	real exp_neg4phi = calc_exp_neg4phi(U);
	//TODO minimize using these 
	real exp_4phi = 1. / exp_neg4phi;

<? if eqn.useChi then ?>
	real partial_phi_l[3] = {
<? for i,xi in ipairs(xNames) do
?>		-partial_chi_l[<?=i-1?>] / (4. * U->chi),
<? end
?>	};

	real partial2_phi_ll[6] = {
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.25 * (-partial2_chi_ll[<?=ij-1?>] + partial_chi_l[<?=i-1?>] * partial_chi_l[<?=j-1?>] / U->chi) / U->chi,
<? end
?>	};
<? end ?>
	
	//gamma_ij = exp(4 phi) gammaBar_ij
	//only used in K_ll calculation for H constraint
	// exp_4phi could be singular
	sym3 gamma_ll = sym3_scale(U->gammaBar_ll, exp_4phi);

	//gammaBar_ij = exp(-4 phi) gamma_ij
	//gammaBar^ij = exp(4 phi) gamma^ij
	//gamma^ij = exp(-4 phi) gammaBar^ij
	sym3 gamma_uu = sym3_scale(U->gammaBar_uu, exp_neg4phi);

	//connBar_lll[i].jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll.<?=xi?>.<?=xjk?> = .5 * (
		partial_gammaBar_lll[<?=k-1?>].<?=sym(i,j)?> 
		+ partial_gammaBar_lll[<?=j-1?>].<?=sym(i,k)?> 
		- partial_gammaBar_lll[<?=i-1?>].<?=xjk?>);
<?	end
end
?>	
	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	_3sym3 connBar_ull = sym3_3sym3_mul(U->gammaBar_uu, connBar_lll);	

	//from 2006 Campanelli "connBar^i is replaced by -gammaBar^ij_j wherever it is not differentiated"
	// TODO ... but do that, I should track gammaBar^ij ...
<? 
if eqn.useChi then
?>	real3 connBar_u = _real3(0,0,0);
	//connBar^i = -gammaBar^ij_,j
<?	for i,xi in ipairs(xNames) do
?>	connBar_u.<?=xi?> =<?
		for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
		end
?>;
<?
	end
else	-- not eqn.useChi
?>	real3 connBar_u = U->connBar_u;
<? 
end
?>
	//DBar^i phi = gammaBar^ij phi_,j
	real3 DBar_phi_u = sym3_real3_mul(U->gammaBar_uu, *(real3*)partial_phi_l);

	//conn_ull[i].jk := conn^i_jk
	//Alcubierre 2.8.14:
	//conn^i_jk = connBar^i_jk + 2 (delta^i_j phi_,k + delta^i_k phi_,j - gamma_jk gamma^il phi_,l)
	//B&S 3.7:
	//conn^i_jk = connBar^i_jk + 2 (delta^i_j phi_,k + delta^i_k phi_,j - gammaBar_jk gammaBar^il phi_,l)
	//conn^i_jk = connBar^i_jk + 2 (delta^i_j phi_,k + delta^i_k phi_,j - gammaBar_jk DBar^i phi)
	_3sym3 conn_ull;	
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	conn_ull.<?=xi?>.<?=xjk?> = connBar_ull.<?=xi?>.<?=xjk?> - 2 * U->gammaBar_ll.<?=xjk?> * DBar_phi_u.<?=xi?><?
		if i==j then
?> + 2 * partial_phi_l[<?=k-1?>]<?
		end
		if i==k then
?> + 2 * partial_phi_l[<?=j-1?>]<?
		end
?>;
<?	end
end
?>

<? 
if eqn.useChi then
?>	
	real partial_chi_bardot_partial_alpha = real3_weightedDot( *(real3*)partial_chi_l, *(real3*)partial_alpha_l, U->gammaBar_uu);
	
	//chi D_i D_j alpha = alpha,ij - connBar^k_ij alpha,k + 1/(2 chi) (alpha,i chi,j + alpha,j chi,i - gammaBar_ij gammaBar^kl alpha,k chi,l)
	sym3 chi_D2_alpha_ll = (sym3){
<? 	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>		.<?=xij?> = U->chi * (partial2_alpha_ll[<?=ij-1?>] <?
		for k,xk in ipairs(xNames) do
			?> - connBar_ull.<?=xk?>.<?=xij?> * partial_alpha_l[<?=k-1?>]<?
		end
			?>) + .5 * partial_alpha_l[<?=i-1?>] * partial_chi_l[<?=j-1?>] <?
			?> + .5 * partial_alpha_l[<?=i-1?>] * partial_chi_l[<?=j-1?>] <?
			?> - .5 * U->gammaBar_ll.<?=xij?> * partial_chi_bardot_partial_alpha,
<? 	end
?>	};

<? 	if false then -- why bother when I need to use chi D_i D_j alpha?  why not just take the trace?
?>
	//chi gammaTilde^ij alpha,ij - chi connBar^i alpha,i - 1/2 gammaTilde^ij chi,i alpha,j
	real tr_D2_alpha = U->chi * (
			sym3_dot(U->gammaBar_uu, *(sym3*)partial2_alpha_ll)
			- real3_dot(U->connBar_u, *(real3*)partial_alpha_l)
		)
		- .5 * real3_weightedDot(
			*(real3*)partial_chi_l,
			*(real3*)partial_alpha_l,
			U->gammaBar_uu);
<? 	else 
?>	real tr_D2_alpha = sym3_dot(U->gammaBar_uu, chi_D2_alpha_ll);

<?	end
else	-- not eqn.useChi
?>	
	
	//D2_alpha_ll.ij = D_i D_j alpha = partial_i partial_j alpha - conn^k_ij partial_k alpha
	sym3 D2_alpha_ll = _sym3(0,0,0,0,0,0);
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	D2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]<?
	for k,xk in ipairs(xNames) do 
?> - conn_ull.<?=xk?>.<?=xij?> * partial_alpha_l[<?=k-1?>]<?
	end ?>;
<? 
end
?>

	//gamma^ij D_i D_j alpha
	real tr_D2_alpha = sym3_dot(gamma_uu, D2_alpha_ll);
<? end	-- eqn.useChi
?>

	//Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	real Q = calc_f(U->alpha) * U->K;
	
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q + real3_dot(*(real3*)partial_alpha_l, U->beta_u);

<? if eqn.useChi then ?>
	//2006 Campanelli eqn 2
	//chi,t = 2/3 chi (alpha K - beta^i_,i) + beta^i chi_,i
	deriv->chi += 2./3. * U->chi * (U->alpha * U->K - tr_partial_beta) + real3_dot(U->beta_u, *(real3*)partial_chi_l); 
<? else ?>
	//B&S 11.50
	//Alcubierre 2.8.10
	//phi,t = -1/6 alpha K + beta^i phi,i + 1/6 beta^i_,i
	deriv->phi += -U->alpha * U->K / 6. + real3_dot(U->beta_u, *(real3*)partial_phi_l) + tr_partial_beta / 6.;
<? end ?>

	//B&S 11.51
	//Alcubierre 2.8.9
	//gammaBar_ij,t = -2 alpha ATilde_ij + beta^k gammaBar_ij,k + gammaBar_ik beta^k_,j + gammaBar_kj beta^k_,i - 2/3 gammaBar_ij beta^k_,k
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->gammaBar_ll.<?=xij?> += -2 * U->alpha * U->ATilde_ll.<?=xij?>	//-2 alpha ATilde_ij 
<? 	for k,xk in ipairs(xNames) do
?>		+ partial_gammaBar_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>		//+ beta^k gammaBar_ij,k 
		+ U->gammaBar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>	//+ gammaBar_jk beta^k_,i 
		+ U->gammaBar_ll.<?=sym(k,i)?> * partial_beta_ul[<?=j-1?>].<?=xk?> 	//+ gammaBar_ik beta^k_,j
<? 	end
?>		- 2./3. * U->gammaBar_ll.<?=xij?> * tr_partial_beta;				//- 2/3 gammaBar_ij beta^k_,k
<? end
?>
	mat3 ATilde_ul = sym3_sym3_mul(U->gammaBar_uu, U->ATilde_ll);		//ATilde^i_j = gammaBar^kl ATilde_kj
	sym3 ATilde_uu = mat3_sym3_to_sym3_mul(ATilde_ul, U->gammaBar_uu);	//ATilde^ij = gammaBar^ik ATilde_kl gammaBar^lj
	real tr_ATilde_sq = sym3_dot(U->ATilde_ll, ATilde_uu);			//tr_ATilde_sq := tr(ATilde^2) = ATilde_ij ATilde^ji
	
	real S = sym3_dot(U->S_ll, gamma_uu);
	
	//B&S 11.52
	//Alcubierre 2.8.12
	//K_,t = -gamma^ij D_i D_j alpha + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	deriv->K += -tr_D2_alpha
		+ U->alpha * (tr_ATilde_sq + U->K * U->K / 3.) 
		+ 4. * M_PI * U->alpha * (U->rho + S) 
		+ real3_dot(U->beta_u, *(real3*)partial_K_l);

	//tr_partial2_gammaBar_ll.ij = gammaBar^kl gammaBar_ij,kl
	sym3 tr_partial2_gammaBar_ll;	
<? for ij,xij in ipairs(symNames) do
?>	tr_partial2_gammaBar_ll.<?=xij?> = 0. <?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?> + U->gammaBar_uu.<?=sym(k,l)?> * partial2_gammaBar_llll[<?=from3x3to6(k,l)-1?>].<?=xij?><?
		end
	end
?>;
<? end
?>
	//B&S 11.54
	//Alcubierre eqn 2.8.17
	//RBar_ij = -1/2 gammaBar^lm gammaBar_ij,lm 
	//		+ 1/2 gammaBar_ki connBar^k_,j
	//		+ 1/2 gammaBar_kj connBar^k_,i 
	//		+ 1/2 connBar^k (connBar_ijk + connBar_jik)
	// 		+ gammaBar^lm (
	//			connBar^k_li connBar_jkm
	//			+ connBar^k_lj connBar_ikm
	//			+ connBar^k_im connBar_klj)
	sym3 RBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	RBar_ll.<?=xij?> = -.5 * tr_partial2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		+ .5 * U->gammaBar_ll.<?=sym(k,i)?> * partial_connBar_ul[<?=j-1?>].<?=xk?>
		+ .5 * U->gammaBar_ll.<?=sym(k,j)?> * partial_connBar_ul[<?=i-1?>].<?=xk?>
		+ .5 * connBar_u.<?=xk?> * (connBar_lll.<?=xi?>.<?=sym(j,k)?> + connBar_lll.<?=xj?>.<?=sym(i,k)?>)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>		+ U->gammaBar_uu.<?=sym(k,m)?> * (
			+ connBar_ull.<?=xk?>.<?=sym(l,i)?> * connBar_lll.<?=xj?>.<?=sym(k,m)?>
			+ connBar_ull.<?=xk?>.<?=sym(l,j)?> * connBar_lll.<?=xi?>.<?=sym(k,m)?>
			+ connBar_ull.<?=xk?>.<?=sym(i,m)?> * connBar_lll.<?=xk?>.<?=sym(l,j)?>
		)
<?			end
		end
	end
?>	;
<? end
?>

	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll[<?=ij-1?>] <?
	for k,xk in ipairs(xNames) do
?> - connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l[<?=k-1?>]<?
	end
?>;
<? end
?>
	real tr_DBar2_phi = sym3_dot(U->gammaBar_uu, DBar2_phi_ll);
	real DBar_phi_norm = real3_dot(*(real3*)partial_phi_l, DBar_phi_u);

	//Baumgarte & Shapiro p.57 eqn 3.10
	//R_ll(i,j) := R_ij = RBar_ij - 2 (DBar_i DBar_j ln(psi) + gammaBar_ij gammaBar^lm DBar_l DBar_m ln(psi)) + 4((DBar_i ln(psi)) (DBar_j ln(psi)) - gammaBar_ij gammaBar^lm (DBar_l ln(psi)) (DBar_m ln(psi)))
	//Then Baumgarte & Shapiro on p.390 say RPhi_ij is the same as p.57 substituting phi for ln(psi)
	// ... but I thought phi was ln(psi)?  Then why would you need to separate R_ij = RBar_ij + RPhi_ij ?  I thought the substitution showed that R_ij was RPhi_ij?
	//phi = ln(psi), so DBar_i ln psi = partial_phi_i
	//Alcubierre 2.8.18
	//RPhi_ll.xij := -2 DTilde_i DTilde_j phi - 2 gammaTilde_ij gammaTilde^kl DTilde_k DTilde_l phi + 4 DTilde_i phi DTilde_j phi - 4 gammaTilde_ij DTilde^k phi DTilde_k phi
	//	= -2 (DTilde_i DTilde_j phi)
	//		- 2 gammaTilde_ij gammaTilde^kl (DTilde_k DTilde_l phi)
	//		+ 4 phi_,i phi_,j 
	//		- 4 gammaTilde_ij gammaTilde^kl phi_,k phi_,l
	//it looks like Alcubierre agrees with Baumgarte & Shapiro, except without the extra RBar_ij ...
	sym3 RPhi_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	RPhi_ll.<?=xij?> = 2. * (
		- DBar2_phi_ll.<?=xij?> 
		- U->gammaBar_uu.<?=xij?> * tr_DBar2_phi 
		+ 2. * (partial_phi_l[<?=i-1?>] * partial_phi_l[<?=j-1?>] 
			- U->gammaBar_ll.<?=xij?> * DBar_phi_norm));
<? end 
?>
	
	sym3 R_ll = sym3_add(RPhi_ll, RBar_ll);

<? if eqn.useChi then
?>
	//TODO instead of R_ll, add in chi_tf_RPhi and chi_tf_RBar_ll
	//then factor out one of the chis 
	
	//(chi alpha (R_ij - 8 pi S_ij))^TF
	sym3 tf_chi_alpha_R_minus_S = sym3_scale(
		sym3_add(R_ll, sym3_scale(U->S_ll, -8. * M_PI)), 
		U->chi * U->alpha);
	tf_chi_alpha_R_minus_S = tracefree(tf_chi_alpha_R_minus_S, U->gammaBar_ll, U->gammaBar_uu);

	//(chi D_i D_j alpha)^TF
	//= chi D_i D_j alpha + gammaBar_ij ( 1/3 chi (connBar^k alpha,k - gammaBar^kl alpha,kl) + 1/6 gammaBar^kl chi,k alpha,l)
	sym3 chi_tf_D2_alpha = sym3_add( 
		chi_D2_alpha_ll,
		sym3_scale(U->gammaBar_ll, 
			1./3. * U->chi * (
				real3_dot(U->connBar_u, *(real3*)partial_alpha_l)
				- sym3_dot(U->gammaBar_uu, *(sym3*)partial2_alpha_ll)
			) + 1./6. * real3_weightedDot(
				*(real3*)partial_alpha_l,
				*(real3*)partial_chi_l,
				U->gammaBar_uu
			)
		)
	);

	sym3 chi_tracelessPart_ll = sym3_add(tf_chi_alpha_R_minus_S, chi_tf_D2_alpha);

<? else
?>
	//traceless portion of -D^2 alpha + alpha (R_ij - 8 pi S_ij)
	//...times chi = exp(-4 phi)
#if 1	//all at once
	sym3 tracelessPart_ll = sym3_sub(
		sym3_scale(
			sym3_add(R_ll, sym3_scale(U->S_ll, -8. * M_PI)), 
			U->alpha),
		D2_alpha_ll);
	tracelessPart_ll = tracefree(tracelessPart_ll, U->gammaBar_ll, U->gammaBar_uu);
#else	//each term separately
	sym3 tracelessPart_ll = sym3_sub(
		sym3_scale(
			sym3_add(
				tracefree(R_ll, gamma_ll, gamma_uu),
				sym3_scale(
					tracefree(U->S_ll, gamma_ll, gamma_uu), 
					-8. * M_PI)), 
			U->alpha),
		tracefree(D2_alpha_ll, gamma_ll, gamma_uu)
	);
#endif
	sym3 chi_tracelessPart_ll = sym3_scale(tracelessPart_ll, exp_neg4phi);
<? end
?>

	//B&S 11.53
	//Alcubierre 2.8.11
	//ATilde_ij,t = 
	//	exp(-4phi) (-(D_ij alpha) + alpha (R_ij - 8 pi S_ij) )^TF
	//	+ alpha (K ATilde_ij - 2 ATilde_il ATilde^l_j)
	//	+ beta^k ATilde_ij,k 
	//	+ ATilde_ik beta^k_,j 
	//	+ ATidle_kj beta^k_,i 
	//	- 2/3 ATilde_ij beta^k_,k
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->ATilde_ll.<?=xij?> += chi_tracelessPart_ll.<?=xij?>
		+ U->alpha * U->K * U->ATilde_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do
?>		- 2. * U->alpha * U->ATilde_ll.<?=sym(i,k)?> * ATilde_ul.<?=xk?>.<?=xj?>
		+ partial_ATilde_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>
		+ U->ATilde_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
		+ U->ATilde_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
<?	end
?>		- 2./3. * U->ATilde_ll.<?=xij?> * tr_partial_beta;
<? end
?>
	
	//connBar^i is the connection function / connection coefficient iteration with Hamiltonian constraint baked in (Baumgarte & Shapiro p.389, Alcubierre p.86).
	//B&S 11.55
	//Alcubierre 2.8.25
	//partial_t connBar^i = 
	//	-2 ATilde^ij alpha_,j
	//	+ 2 alpha (
	//		connBar^i_jk ATilde^kj 
	//		- 2/3 gammaBar^ij K_,j
	//		- 8 pi gammaBar^ij S_j 
	//		+ 6 ATilde^ij phi_,j
	//	)
	//	+ beta^j connBar^i_,j
	//	- connBar^j beta^i_,j
	//	+ 2/3 connBar^i beta^j_,j
	//	+ 1/3 gammaBar^ki beta^j_,jk
	//	+ gammaBar^kj beta^i_,jk
	real3 dt_connBar_u;
<? for i,xi in ipairs(xNames) do
?>	dt_connBar_u.<?=xi?> =
		2./3. * connBar_u.<?=xi?> * tr_partial_beta
		- 16. * M_PI * exp_4phi * U->alpha * U->S_u.<?=xi?> 
<?	for j,xj in ipairs(xNames) do
		local xij = sym(i,j)
		local jj = from3x3to6(j,j)
?>		- 2. * ATilde_uu.<?=xij?> * partial_alpha_l[<?=j-1?>]
		+ 2. * U->alpha * (
			-2./3. * U->gammaBar_uu.<?=xij?> * partial_K_l[<?=j-1?>] 
			+ 6. * ATilde_uu.<?=xij?> * partial_phi_l[<?=j-1?>])
		+ U->beta_u.<?=xi?> * partial_connBar_ul[<?=j-1?>].<?=xi?>
		- connBar_u.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?		for k,xk in ipairs(xNames) do		
			local xik = sym(i,k)
			local jk = from3x3to6(j,k)
			local xjk = symNames[jk]
?>		+ 2. * U->alpha * connBar_ull.<?=xi?>.<?=xjk?> * ATilde_uu.<?=xjk?>
		+ 1./3. * U->gammaBar_uu.<?=xik?> * partial2_beta_ull[<?=jk-1?>].<?=xj?>
		+ U->gammaBar_uu.<?=xjk?> * partial2_beta_ull[<?=jk-1?>].<?=xi?>
<?		end
	end
?>	;
<? end
?>
	deriv->connBar_u = real3_add(deriv->connBar_u, dt_connBar_u);

<? 
if eqn.guiVars.useGammaDriver.value
and eqn.useHypGammaDriver
then
	error("you can only enable either useGammaDriver or useHypGammaDriver.  I should make this a combo gui option")
end
?>
<? if eqn.guiVars.useGammaDriver.value then ?>
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k (connBar^i_,t + eta connBar^i)
	const real k = 0.;//3./4.;
	const real eta = 0.;	// 1 / (2 M), for total mass M
	deriv->beta_u = real3_add(deriv->beta_u,
		real3_add(
			real3_scale(dt_connBar_u, k),
			real3_scale(connBar_u, eta)));
<? end ?>
<? if eqn.useHypGammaDriver then ?>
	//hyperbolic Gamma driver 
	//B&S 4.83 
	//beta^i_,t = 3/4 B^i (Rezolla says 3/4 alpha B^i, 2006 Campanelli says B^i but absorbs 3/4 into B^i)
	//B^i_,t = connBar^i_,t - eta B^i ... maybe + beta^j B^i_,j (Rezolla says - beta^j B^i_,j)
	const real eta = 0.;
	deriv->beta_u = real3_add(deriv->beta_u,
		real3_scale(U->B_u, 3./4. * U->alpha));

	deriv->B_u = real3_add(deriv->B_u, dt_connBar_u);
	deriv->B_u = real3_sub(deriv->B_u, real3_scale(U->B_u, eta));
<? 	for i,xi in ipairs(xNames) do
		for j,xj in ipairs(xNames) do
?>	deriv->B_u.<?=xi?> -= U->beta_u.<?=xj?> * partial_B_ul[<?=i-1?>].<?=xj?>;
<? 		end
	end
end ?>

<? if calcConstraints then ?>
	
#if 0
	real RBar = sym3_dot(U->gammaBar_uu, RBar_ll);
	real exp_phi = exp(U->phi);
	
	//B&S 11.48
	//
	//exp(phi)_,ij = partial2_exp_phi_ll.ij
	//= (phi_,i exp(phi))_;j
	//= exp(phi) (phi_,ij + phi_,i phi_,j)
	//
	//DBar_j DBar_i exp(phi)
	//= DBar_j exp(phi)_,i 
	//= DBar_j (phi_,i exp(phi))
	//= exp(phi) (phi_,ij - connBar^k_ij phi_,k + phi_,i phi_,j)
	//= exp(phi) (DBar2_ij phi + phi_,i phi_,j)
	//
	//gammaBar^ij DBar_i DBar_j exp(phi)
	// = exp(phi) (gammaBar^ij DBar_i DBar_j phi + gammaBar^ij phi_,i phi_,j)
	// = exp(phi) (tr_DBar2_phi + DBar_phi_norm)

	//H = gammaBar^ij DBar_i DBar_j e^phi 
	//		- 1/8 e^phi RBar 
	//		+ 1/8 e^(5 phi) ATilde_ij ATilde^ij 
	//		- 1/12 e^(5 phi) K^2 
	//		+ 2 pi e^(5 phi) rho
	//= e^(phi) (
	//		tr_DBar2_phi 
	//		+ DBar_phi_norm
	//		- 1/8 RBar 
	//		+ e^(4phi) (
	//			+ 1/8 ATilde_ij ATilde^ij 
	//			- 1/12 K^2 
	//			+ 2 pi rho
	//	)	)
	U->H = exp_phi * (
		tr_DBar2_phi 
		+ DBar_phi_norm
		- 1./8. * RBar
		+ exp_4phi * (
			+ 1./8. * tr_ATilde_sq
			- 1./12. * U->K * U->K
			+ 2. * M_PI * U->rho
		)
	);
#endif
#if 1
	real R = sym3_dot(gamma_uu, R_ll);	//R = gamma^ij R_ij
	
	//K_ij = A_ij + 1/3 gamma_ij K 
	//A_ij = exp(4 phi) ATilde_ij
	//K_ij = exp(4 phi) ATilde_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_scale(U->ATilde_ll, exp_4phi),
		sym3_scale(gamma_ll, U->K/3.));
	mat3 K_ul = sym3_sym3_mul(gamma_uu, K_ll);
	sym3 K_uu = mat3_sym3_to_sym3_mul(K_ul, gamma_uu);
	
	//Alcubierre 3.1.1
	U->H = R + U->K * U->K - sym3_dot(K_uu, K_ll) - 16. * M_PI * U->rho;
#endif

#if 0
	//DBar_j (e^(6 phi) ATilde^ij)
	//= DBar_j (e^(6 phi)) ATilde^ij + e^(6 phi) DBar_j ATilde^ij
	//= 6 phi_,j e^(6 phi) ATilde^ij + e^(6 phi) (ATilde^ij_,j + connBar^i_kj ATilde^kj + connBar^j_kj ATilde^ik)
	//= exp(6 phi) (6 ATilde^ij phi_,j + (gammaBar^ik ATilde_kl gammaBar^lj)_,j + connBar^i_jk ATilde^jk) ... plus (ln det gammaBar_ll)_,k which is zero, right?
	//= exp(6 phi) (
	//		+ 6 ATilde^ij phi_,j
	//		+ gammaBar^ik_,j ATilde_k^j
	//		+ ATilde^i_l gammaBar^lj_,j
	//		+ gammaBar^ik ATilde_kl,j gammaBar^lj
	//		+ connBar^i_jk ATilde^jk)
   	//= exp(6 phi) (
	//		+ 6 ATilde^ij phi_,j
	//		- ATilde^kj gammaBar^il gammaBar_lk,j
	//		- ATilde^ik gammaBar^mj gammaBar_km,j
	//		+ gammaBar^ik gammaBar^lj ATilde_kl,j
	//		+ connBar^i_jk ATilde^jk)

	//B&S 11.49
	//0 = M^i = DBar_j (e^(6 phi) ATilde^ij) - 2/3 e^(6 phi) DBar^i K - 8 pi e^(6 phi) S^i
#endif

<? end	--calcConstraints ?>
}

kernel void addSource(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS_NOGHOST();

	global <?=eqn.cons_t?>* U = UBuf + index;

	//Kreiss-Oligar dissipation
	//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
	//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
	for (int i = 0; i < numIntStates; ++i) {
<?=makePartial2('ptr[i]', 'real', 'partial2_Ui_ll')?>	
		real lap = 0<?
for j,xj in ipairs(xNames) do
	local jj = from3x3to6(j,j)
?> + partial2_Ui_ll[<?=jj-1?>]<?
end
?>;
		U->ptr[i] -= gui_diffuseSigma/16. * lap;
	}
}
