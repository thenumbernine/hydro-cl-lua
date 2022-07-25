//// MODULE_NAME: <?=eqn_common?>

//gammaBar_ij = gammaHat_ij + epsilon_ij
sym3 <?=calc_gammaBar_ll?>(global <?=cons_t?> const * const U, real3 const x) {
	sym3 gammaHat_ll = coord_gHol_ll(x);
	return sym3_add(gammaHat_ll, U->epsilon_ll);
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2017 Ruchlin
real calc_det_gammaBar_ll(real3 const x) {
	return coord_det_gHol(x);
}

void <?=setFlatSpace?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const U,
	real3 const x
) {
	U->alpha = 1.;
	U->beta_u = real3_zero;
	U->epsilon_ll = sym3_zero;
	U->chi = 1;
	U->KHat = 0;
	U->Theta = 0;
	U->ABar_ll = sym3_ident;
	U->Delta_u = real3_zero;
<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_u = real3_zero;
<? end
?>	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	real det_gammaBar_ll = calc_det_gammaBar_ll(x);
	U->gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar_ll);

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_u = real3_zero;
}

#define <?=calc_exp_neg4phi?>(U) ((U)->chi)

//det(gamma_ij) = exp(12 phi) det(gammaBar_ij)
//				= det(gammaHat_ij) / (exp(-4 phi)^3) 
real calc_det_gamma_ll(global <?=cons_t?> const * const U, real3 const x) {
	real const exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	return calc_det_gammaBar_ll(x) / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
}

sym3 <?=calc_gamma_uu?>(global <?=cons_t?> const * const U, real3 const x) {
	real const exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	sym3 const gamma_uu = sym3_real_mul(U->gammaBar_uu, exp_neg4phi);
	return gamma_uu;
}

sym3 <?=calc_gamma_ll?>(global <?=cons_t?> const * const U, real3 const x) {
	sym3 const gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	real const exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);
	sym3 const gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	return gamma_ll;
}

/*
TF(K_ij) = K_ij - 1/3 gamma_ij gamma^kl K_kl

tr(A_ij)
= tr(K_ij - 1/3 gamma_ij K)
= gamma^ij K_ij - 1/3 gamma^ij gamma_ij K
= K - 1/3 3 K
= 0

tr(ABar_ij) = exp(-4 phi) tr(A_ij) 
= exp(-4 phi) * 0
= 0

using gammaBar metrics shouldn't make a difference (since the weights of the metric and inverse should cancel)
TFBar(K_ij) = K_ij - 1/3 gammaBar_ij gammaBar^kl K_kl 
	= K_ij - 1/3 gamma_ij gamma^kl K_kl
	= TF(K_ij)


tr(TF(K_ij)) = tr(K_ij) - 1/3 tr(gamma_ij K)
	= K - 1/3 * K * gamma^ij gamma_ij
	= K - 1/3 * K * 3 = K - K = 0

trace-free is linear:
tr(TF(A_ij + B_ij)) = gamma^ij ( A_ij + B_ij - 1/3 gamma_ij gamma^kl (A_kl + B_kl) )
	= A + B - 1/3 * 3 * (A + B)
	= A + B - A - B
	= 0
*/
sym3 tracefree(sym3 const A_ll, sym3 const gamma_ll, sym3 const gamma_uu) {
	real const tr_A = sym3_dot(A_ll, gamma_uu);
	return sym3_sub(A_ll, sym3_real_mul(gamma_ll, tr_A / 3.));
}

//// MODULE_NAME: <?=applyInitCondCell?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	
	sym3 gammaHat_ll = coord_gHol_ll(x);
	sym3 gamma_ll = gammaHat_ll;
	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

	<?=initCond()?>

	U->alpha = alpha;
	U->beta_u = beta_u;

	real det_gamma_ll = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma_ll);

	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar_ll = calc_det_gammaBar_ll(x); 
	
	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = cbrt(det_gammaBar_ll / det_gamma_ll);
	U->chi = exp_neg4phi;
	
	sym3 gammaBar_ll = sym3_real_mul(gamma_ll, exp_neg4phi);
	U->epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar_ll);

<? if false then ?>
<? for _,x in ipairs(xNames) do
?>	U->a.<?=x?> = calc_a_<?=x?>(x.x, x.y, x.z);
<? end ?>	
<? end ?>

	U->Theta = 0.;	//TODO ... Theta = -Z^mu n_mu = alpha * Z^t ... which is?

	real K = sym3_dot(K_ll, gamma_uu);
	U->KHat = K - 2. * U->Theta;
	
	sym3 A_ll = sym3_sub(K_ll, sym3_real_mul(gamma_ll, 1./3. * K));
	U->ABar_ll = sym3_real_mul(A_ll, exp_neg4phi);
	
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	
	U->H = 0.;
	U->M_u = real3_zero;
}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

//after popularing gammaBar_ll, use its finite-difference derivative to initialize connBar_u
kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf;
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const U = UBuf + index;
	
<?=eqn:makePartial1'gammaBar_uu'?>

	//connBar^i = connBar^i_jk gammaBar^jk
	// TODO is this still true?
	//= -gammaBar^ij_,j + 2 gammaBar^ij Z_j
	//= gammaBar_jk,l gammaBar^ij gammaBar^kl + 2 gammaBar^ij Z_j
	real3 connBar_u;
<? for i,xi in ipairs(xNames) do
?>	connBar_u.<?=xi?> =<?
	for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
	end ?>;
<? end ?>
	
	//Delta^i = gammaBar^jk Delta^i_jk = gammaBar^jk (connBar^i_jk - connHat^i_jk)
	//= gammaBar^jk Delta^i_jk = connBar^i - connHat^i
	_3sym3 connHat_ull = coord_connHol_ull(x);
	real3 connHat_u = _3sym3_sym3_dot23(connHat_ull, U->gammaBar_uu);
	U->Delta_u = real3_sub(connBar_u, connHat_u);
}


//// MODULE_NAME: <?=constrainU?>

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const U = UBuf + index;

<? 
if eqn.guiVars.constrain_det_gammaBar_ll.value 
or eqn.guiVars.constrain_tr_ABar_ll.value 
then 
?>
	sym3 gammaHat_ll = coord_gHol_ll(x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, U->epsilon_ll);

	/*
	if the background is flat ...
	then we need to force det(gammaBar_ij) = 1
	and we can do so with gammaBar_ij := gammaBar_ij / det(gammaBar_ij)^(1/3)
	so we find
	det(gammaBar_ij / det(gammaBar_ij)^(1/3))
	= det(gammaBar_ij)) / det(gammaBar_ij)
	= 1
	
	if the background is gammaHat_ij
	then we need to force det(gammaBar_ij) = det(gammaHat_ij)
	do so with gammaBar_ij := gammaBar_ij * ( det(gammaHat_ij) / det(gammaBar_ij) )^(1/3)
	*/
<? 	if eqn.guiVars.constrain_det_gammaBar_ll.value then ?>
	real det_gammaHat_ll = coord_det_gHol(x);
	real det_gammaBar_ll = sym3_det(gammaBar_ll);
	real rescaleMetric = cbrt(det_gammaHat_ll/det_gammaBar_ll);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_ll.<?=xij?> *= rescaleMetric;
<? 		end ?>
	U->epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
<? 	end ?>

	//here we refresh gammaBar_uu
	U->gammaBar_uu = sym3_inv(gammaBar_ll, 1.);
	
	//in Buchman's paper it says he doesn't do this
	//likewise in my own experiences, this can tend A to grow out of control 
<? 	if eqn.guiVars.constrain_tr_ABar_ll.value then ?>
	U->ABar_ll = tracefree(U->ABar_ll, gammaBar_ll, U->gammaBar_uu);
<? 	end
end
?>

	//TODO 'if constrain chi...'
	real const chiMin = 1e-3;
	U->chi = max(U->chi, chiMin);

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);
}

//// MODULE_NAME: <?=calcDeriv?>

kernel void <?=calcDeriv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const deriv = derivBuf + index;

	global <?=cons_t?> const * const U = UBuf + index;

<?=makePartial('alpha', 'real')?>		//partial_alpha_l[i] := alpha_,i
<?=makePartial('chi', 'real')?>			//partial_chi_l[i] := chi_,i 
<?=makePartial('Theta', 'real')?>		//partial_Theta_l[i] := Theta_,i 
<?=makePartial('gammaBar_uu', 'sym3')?>	//partial_gammaBar_uul[k].ij = gamma^ij_,k
<?=makePartial('KHat', 'real')	?>		//partial_KHat_l[i] := KHat,i
<?=makePartial('beta_u', 'real3')?>		//partial_beta_ul[j].i := beta^i_,j
<?=makePartial('Delta_u', 'real3')?>	//partial_Delta_ul[j].i := Delta^i_,j
<?=makePartial('epsilon_ll', 'sym3')?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k for static grids
<?=makePartial('ABar_ll', 'sym3')?>	//partial_ABar_lll[k].ij = ABar_ij,k
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=makePartial('B_u', 'real3')?>
<? end ?>
	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

<?=eqn:makePartial2'alpha'?>			//partial2_alpha_ll.ij := alpha_,ij
<?=eqn:makePartial2'chi'?>			//partial2_chi_ll.ij := chi_,ij
<?=eqn:makePartial2'epsilon_ll'?>		//partial2_epsilon_llll[kl].ij = epsilon_ij,kl = gammaBar_ij,kl for static grids
<?=eqn:makePartial2'beta_u'?>		//partial2_beta_ull[jk].i = beta^i_,jk

	real exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	//TODO minimize using these 
	real exp_4phi = 1. / exp_neg4phi;

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

	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);

	//gamma_ij = exp(4 phi) gammaBar_ij
	//only used in K_ll calculation for H constraint
	// exp_4phi could be singular
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);

	//gammaBar_ij = exp(-4 phi) gamma_ij
	//gammaBar^ij = exp(4 phi) gamma^ij
	//gamma^ij = exp(-4 phi) gammaBar^ij
	sym3 gamma_uu = sym3_real_mul(U->gammaBar_uu, exp_neg4phi);

	//connBar_lll[i].jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll.<?=xi?>.<?=xjk?> = .5 * (
		partial_epsilon_lll[<?=k-1?>].<?=sym(i,j)?> 
		+ partial_epsilon_lll[<?=j-1?>].<?=sym(i,k)?> 
		- partial_epsilon_lll[<?=i-1?>].<?=xjk?>);
<?	end
end
?>	
	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	_3sym3 connBar_ull = sym3_3sym3_mul(U->gammaBar_uu, connBar_lll);	

	//from 2006 Campanelli "connBar^i is replaced by -gammaBar^ij_j wherever it is not differentiated"
	// TODO ... but do that, I should track gammaBar^ij ...
 
	real3 connBar_u = real3_zero;
	//connBar^i = -gammaBar^ij_,j
<?	for i,xi in ipairs(xNames) do
?>	connBar_u.<?=xi?> =<?
		for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
		end
?>;
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
?>	conn_ull.<?=xi?>.<?=xjk?> = connBar_ull.<?=xi?>.<?=xjk?> - 2 * gammaBar_ll.<?=xjk?> * DBar_phi_u.<?=xi?><?
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

#if 0	//where is this equation from?
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
			?> - .5 * gammaBar_ll.<?=xij?> * partial_chi_bardot_partial_alpha,
<? 	end
?>	};
#else
	//D_i D_j alpha 
	// = D_i alpha_,j 
	// = alpha_,ij - conn^k_ij alpha_,k
	sym3 D2_alpha_ll = (sym3){
<?	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>		.<?=xij?> = partial2_alpha_ll[<?=ij-1?>] 
<?		for k,xk in ipairs(xNames) do
?>			- conn_ull.<?=xk?>.<?=xij?> * partial_alpha_l[<?=k-1?>]
<?		end ?>,
<?	end
?>	};
#endif


#if 0 	//why bother when I need to use chi D_i D_j alpha?  why not just take the trace?
	//chi gammaBar^ij alpha,ij - chi connBar^i alpha,i - 1/2 gammaBar^ij chi,i alpha,j
	real tr_D2_alpha = U->chi * (
			sym3_dot(U->gammaBar_uu, *(sym3*)partial2_alpha_ll)
			- real3_dot(connBar_u, *(real3*)partial_alpha_l)
		)
		- .5 * real3_weightedDot(
			*(real3*)partial_chi_l,
			*(real3*)partial_alpha_l,
			U->gammaBar_uu);
#elif 0
	real tr_D2_alpha = sym3_dot(U->gammaBar_uu, chi_D2_alpha_ll);
#else
	real tr_D2_alpha = sym3_dot(gamma_uu, D2_alpha_ll);
#endif

	//2011 Cao eqn 11
	//KHat = K - 2 Theta
	//K = KHat + 2 Theta
	real K = U->KHat + 2. * U->Theta;

	//Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	//2011 Cao eqn 22 ... where f == mu_L
	//2011 Cau before eqn 24 ... mu_L = 2 / alpha
	real Q = calc_f(U->alpha) * U->KHat;

	//TODO shouldn't there be a Theta in here?  Check 2011 Cao plz
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q + real3_dot(*(real3*)partial_alpha_l, U->beta_u);

	//2011 Cao eqn 12
	//chi,t = 2/3 chi (alpha (KHat + 2 Theta) - beta^i_,i) + chi_,j beta^j
	deriv->chi += 2./3. * U->chi * (U->alpha * K - tr_partial_beta) + real3_dot(*(real3*)partial_chi_l, U->beta_u);

	//B&S 11.51
	//Alcubierre 2.8.9
	//gammaBar_ij,t = -2 alpha ABar_ij + beta^k gammaBar_ij,k + gammaBar_ik beta^k_,j + gammaBar_kj beta^k_,i - 2/3 gammaBar_ij beta^k_,k
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->epsilon_ll.<?=xij?> += -2 * U->alpha * U->ABar_ll.<?=xij?>	//-2 alpha ABar_ij 
<? 	for k,xk in ipairs(xNames) do
?>		+ partial_epsilon_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>		//+ beta^k gammaBar_ij,k 
		+ gammaBar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>	//+ gammaBar_jk beta^k_,i 
		+ gammaBar_ll.<?=sym(k,i)?> * partial_beta_ul[<?=j-1?>].<?=xk?> 	//+ gammaBar_ik beta^k_,j
<? 	end
?>		- 2./3. * gammaBar_ll.<?=xij?> * tr_partial_beta;				//- 2/3 gammaBar_ij beta^k_,k
<? end
?>
	real3x3 ABar_ul = sym3_sym3_mul(U->gammaBar_uu, U->ABar_ll);		//ABar^i_j = gammaBar^kl ABar_kj
	sym3 ABar_uu = real3x3_sym3_to_sym3_mul(ABar_ul, U->gammaBar_uu);	//ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	real tr_ABar_sq = sym3_dot(U->ABar_ll, ABar_uu);			//tr_ABar_sq := tr(ABar^2) = ABar_ij ABar^ji
	
	real S = sym3_dot(U->S_ll, gamma_uu);

	real const kappa1 = 1.;
	real const kappa2 = 0.;
	//2011 Cao eqn 14
	//KHat_,t = -gamma^ij D_i D_j alpha 
	//			+ alpha (ABar_ij ABar^ij + 1/3 (KBar + 2 Theta)^2
	//				+ kappa1 (1 - kappa2) Theta)
	//			+ 4 pi alpha (S + rho_ADM) + beta^i KHat_,i
	deriv->KHat += -tr_D2_alpha
		+ U->alpha * (tr_ABar_sq + K * K / 3.) 
		+ kappa1 * (1. - kappa2) * U->Theta
		+ 4. * M_PI * U->alpha * (U->rho + S) 
		+ real3_dot(U->beta_u, *(real3*)partial_KHat_l);

	//tr_partial2_epsilon_ll.ij = gammaBar^kl gammaBar_ij,kl
	sym3 tr_partial2_epsilon_ll;	
<? for ij,xij in ipairs(symNames) do
?>	tr_partial2_epsilon_ll.<?=xij?> = 0. <?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?> + U->gammaBar_uu.<?=sym(k,l)?> * partial2_epsilon_llll[<?=from3x3to6(k,l)-1?>].<?=xij?><?
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
?>	RBar_ll.<?=xij?> = -.5 * tr_partial2_epsilon_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		+ .5 * gammaBar_ll.<?=sym(k,i)?> * partial_Delta_ul[<?=j-1?>].<?=xk?>
		+ .5 * gammaBar_ll.<?=sym(k,j)?> * partial_Delta_ul[<?=i-1?>].<?=xk?>
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
	//RPhi_ll.xij := -2 DBar_i DBar_j phi - 2 gammaBar_ij gammaBar^kl DBar_k DBar_l phi + 4 DBar_i phi DBar_j phi - 4 gammaBar_ij DBar^k phi DBar_k phi
	//	= -2 (DBar_i DBar_j phi)
	//		- 2 gammaBar_ij gammaBar^kl (DBar_k DBar_l phi)
	//		+ 4 phi_,i phi_,j 
	//		- 4 gammaBar_ij gammaBar^kl phi_,k phi_,l
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
			- gammaBar_ll.<?=xij?> * DBar_phi_norm));
<? end 
?>
	
	sym3 R_ll = sym3_add(RPhi_ll, RBar_ll);

	//TODO instead of R_ll, add in chi_tf_RPhi and chi_tf_RBar_ll
	//then factor out one of the chis 
	
	//(alpha (R_ij - 8 pi S_ij))^TF
	sym3 alpha_R_minus_S_ll = sym3_real_mul(
		sym3_add(R_ll, sym3_real_mul(U->S_ll, -8. * M_PI)), 
		U->alpha);

#if 0	//where is this from?
	//(chi D_i D_j alpha)^TF
	//= chi D_i D_j alpha + gammaBar_ij ( 1/3 chi (connBar^k alpha,k - gammaBar^kl alpha,kl) + 1/6 gammaBar^kl chi,k alpha,l)
	sym3 chi_D2_alpha = sym3_add( 
		chi_D2_alpha_ll,
		sym3_real_mul(gammaBar_ll, 
			1./3. * U->chi * (
				real3_dot(connBar_u, *(real3*)partial_alpha_l)
				- sym3_dot(U->gammaBar_uu, *(sym3*)partial2_alpha_ll)
			) + 1./6. * real3_weightedDot(
				*(real3*)partial_alpha_l,
				*(real3*)partial_chi_l,
				U->gammaBar_uu
			)
		)
	);
#endif

	//chi (alpha (R_ij - 8 pi S_ij) - D_i D_j alpha)^TF
	sym3 dABar_tracelessPart_ll = sym3_sub(alpha_R_minus_S_ll, D2_alpha_ll);
	dABar_tracelessPart_ll = tracefree(dABar_tracelessPart_ll, gammaBar_ll, U->gammaBar_uu);

#if 1	//there's something in here skewed towards the x=y direction...
	//B&S 11.53
	//Alcubierre 2.8.11
	//ABar_ij,t = 
	//	exp(-4phi) (-(D_i D_j alpha) + alpha (R_ij - 8 pi S_ij) )^TF
	//	+ alpha K ABar_ij 
	//	- 2 alpha ABar_ik ABar^k_j
	//	+ beta^k ABar_ij,k 
	//	+ ABar_ik beta^k_,j 
	//	+ ABar_kj beta^k_,i 
	//	- 2/3 ABar_ij beta^k_,k
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->ABar_ll.<?=xij?> += U->chi * dABar_tracelessPart_ll.<?=xij?>
		+ U->alpha * K * U->ABar_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do
?>		- 2. * U->alpha * U->ABar_ll.<?=sym(i,k)?> * ABar_ul.<?=xk?>.<?=xj?>
		+ partial_ABar_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>
		+ U->ABar_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
		+ U->ABar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
<?	end
?>		- 2./3. * U->ABar_ll.<?=xij?> * tr_partial_beta;
<? end
?>
#endif

	//K,i = KHat__,i + 2 Theta_,i
	real3 partial_K_l = real3_add(
		*(real3*)partial_KHat_l,
		real3_real_mul(
			*(real3*)partial_Theta_l,
			2.));

	//connBar^i is the connection function / connection coefficient iteration with Hamiltonian constraint baked in (Baumgarte & Shapiro p.389, Alcubierre p.86).
	//B&S 11.55
	//Alcubierre 2.8.25
	//partial_t connBar^i = 
	//	-2 ABar^ij alpha_,j
	//	+ 2 alpha (
	//		connBar^i_jk ABar^kj 
	//		- 2/3 gammaBar^ij K_,j
	//		- 8 pi gammaBar^ij S_j 
	//		+ 6 ABar^ij phi_,j
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
?>		- 2. * ABar_uu.<?=xij?> * partial_alpha_l[<?=j-1?>]
		+ 2. * U->alpha * (
			-2./3. * U->gammaBar_uu.<?=xij?> * partial_K_l.<?=xj?> 
			+ 6. * ABar_uu.<?=xij?> * partial_phi_l[<?=j-1?>])
		+ U->beta_u.<?=xi?> * partial_Delta_ul[<?=j-1?>].<?=xi?>
		- connBar_u.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?		for k,xk in ipairs(xNames) do		
			local xik = sym(i,k)
			local jk = from3x3to6(j,k)
			local xjk = symNames[jk]
?>		+ 2. * U->alpha * connBar_ull.<?=xi?>.<?=xjk?> * ABar_uu.<?=xjk?>
		+ 1./3. * U->gammaBar_uu.<?=xik?> * partial2_beta_ull[<?=jk-1?>].<?=xj?>
		+ U->gammaBar_uu.<?=xjk?> * partial2_beta_ull[<?=jk-1?>].<?=xi?>
<?		end
	end
?>	;
<? end
?>
	deriv->Delta_u = real3_add(deriv->Delta_u, dt_connBar_u);

<? if eqn.useShift == 'GammaDriver' then ?>
	//Gamma-driver
	//2011 Cao eqn 23
	//beta^i_,t = mu_S alpha^2 connBar^i - eta beta^i + beta^j beta^i_,j
	real mu_S = 1. / (U->alpha * U->alpha);
	real eta = .01;
<? for i,xi in ipairs(xNames) do
?>	deriv->beta_u.<?=xi?> += connBar_u.<?=xi?> * mu_S * U->alpha * U->alpha
		- eta * U->beta_u.<?=xi?>
<?	for j,xj in ipairs(xNames) do
?>		+ partial_beta_ul[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
<?	end ?>;
<? end ?>
	
<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>
	
	//hyperbolic Gamma driver 
	//B&S 4.83 
	//beta^i_,t = 3/4 B^i (Rezolla says 3/4 alpha B^i, 2006 Campanelli says B^i but absorbs 3/4 into B^i)
	//B^i_,t = connBar^i_,t - eta B^i ... maybe + beta^j B^i_,j (Rezolla says - beta^j B^i_,j)
	real const eta = 0.;
	deriv->beta_u = real3_add(deriv->beta_u,
		real3_real_mul(U->B_u, 3./4. * U->alpha));

	deriv->B_u = real3_add(deriv->B_u, dt_connBar_u);
	deriv->B_u = real3_sub(deriv->B_u, real3_real_mul(U->B_u, eta));
<? 	for i,xi in ipairs(xNames) do
		for j,xj in ipairs(xNames) do
?>	deriv->B_u.<?=xi?> -= U->beta_u.<?=xj?> * partial_B_ul[<?=i-1?>].<?=xj?>;
<? 		end
	end
end ?>

//TODO move U->H calculation into another function, or into constrainU
// for some reason, keeping it here, my intel GPU doesn't write U->H
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
	//		+ 1/8 e^(5 phi) ABar_ij ABar^ij 
	//		- 1/12 e^(5 phi) K^2 
	//		+ 2 pi e^(5 phi) rho
	//= e^(phi) (
	//		tr_DBar2_phi 
	//		+ DBar_phi_norm
	//		- 1/8 RBar 
	//		+ e^(4phi) (
	//			+ 1/8 ABar_ij ABar^ij 
	//			- 1/12 K^2 
	//			+ 2 pi rho
	//	)	)
	U->H = exp_phi * (
		tr_DBar2_phi 
		+ DBar_phi_norm
		- 1./8. * RBar
		+ exp_4phi * (
			+ 1./8. * tr_ABar_sq
			- 1./12. * K * K
			+ 2. * M_PI * U->rho
		)
	);
#endif
#if 0
	real R = sym3_dot(gamma_uu, R_ll);	//R = gamma^ij R_ij

	//2011 Cao eqn 16
	deriv->Theta += .5 * U->alpha * (
		R - tr_ABar_sq + 2. / 3. * K - 16. * M_PI * U->rho 
		- 2. * kappa1 * (2. + kappa2) * U->Theta
	) + real3_dot(*(real3*)partial_Theta_l, U->beta_u);

	//K_ij = A_ij + 1/3 gamma_ij K 
	//A_ij = exp(4 phi) ABar_ij
	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_real_mul(U->ABar_ll, exp_4phi),
		sym3_real_mul(gamma_ll, K/3.));
	real3x3 K_ul = sym3_sym3_mul(gamma_uu, K_ll);
	sym3 K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);
	
	//Alcubierre 3.1.1
	U->H = R + K * K - sym3_dot(K_uu, K_ll) - 16. * M_PI * U->rho;
#endif

#if 0
	//DBar_j (e^(6 phi) ABar^ij)
	//= DBar_j (e^(6 phi)) ABar^ij + e^(6 phi) DBar_j ABar^ij
	//= 6 phi_,j e^(6 phi) ABar^ij + e^(6 phi) (ABar^ij_,j + connBar^i_kj ABar^kj + connBar^j_kj ABar^ik)
	//= exp(6 phi) (6 ABar^ij phi_,j + (gammaBar^ik ABar_kl gammaBar^lj)_,j + connBar^i_jk ABar^jk) ... plus (ln det gammaBar_ll)_,k which is zero, right?
	//= exp(6 phi) (
	//		+ 6 ABar^ij phi_,j
	//		+ gammaBar^ik_,j ABar_k^j
	//		+ ABar^i_l gammaBar^lj_,j
	//		+ gammaBar^ik ABar_kl,j gammaBar^lj
	//		+ connBar^i_jk ABar^jk)
   	//= exp(6 phi) (
	//		+ 6 ABar^ij phi_,j
	//		- ABar^kj gammaBar^il gammaBar_lk,j
	//		- ABar^ik gammaBar^mj gammaBar_km,j
	//		+ gammaBar^ik gammaBar^lj ABar_kl,j
	//		+ connBar^i_jk ABar^jk)

	//B&S 11.49
	//0 = M^i = DBar_j (e^(6 phi) ABar^ij) - 2/3 e^(6 phi) DBar^i K - 8 pi e^(6 phi) S^i
#endif

}

//// MODULE_NAME: <?=addSource?>

kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if false then ?>
	<?=SETBOUNDS_NOGHOST?>();

	global <?=cons_t?> const * const U = UBuf + index;

	//Kreiss-Oligar dissipation
	//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
	//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
	for (int i = 0; i < numIntStates; ++i) {
<?=eqn:makePartial2('ptr[i]', 'real', 'partial2_Ui_ll')?>	
		real lap = 0<?
for j,xj in ipairs(xNames) do
	local jj = from3x3to6(j,j)
?> + partial2_Ui_ll[<?=jj-1?>]<?
end
?>;
		deriv->ptr[i] -= solver->diffuseSigma/16. * lap;
	}
<? end ?>
}

//// MODULE_NAME: <?=calcDTCell?>

void <?=calcDTCell?>(
	global real * const dtBuf,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	sym3 const gamma_uu = <?=calc_gamma_uu?>(U, x);

	<? for side=0,solver.dim-1 do ?>{
		real absLambdaMax = U->alpha * sqrt(gamma_uu.<?=sym(side+1,side+1)?>);
		real dx = solver->grid_dx.s<?=side?>;
		*(dt) = (real)min(*(dt), dx / absLambdaMax);
	}<? end ?>
}
