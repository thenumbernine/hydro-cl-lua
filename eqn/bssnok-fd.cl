<? 
local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym

local makePartials = require 'eqn.makepartial'
local derivOrder = 2 * solver.numGhost
local makePartial = function(...) return makePartials.makePartial(derivOrder, solver, ...) end
local makePartial2 = function(...) return makePartials.makePartial2(derivOrder, solver, ...) end
?>

typedef <?=eqn.cons_t?> cons_t;
typedef <?=solver.solver_t?> solver_t;

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

TFBar(K_ij) = K_ij - 1/3 gammaBar_ij gammaBar^kl K_kl 
	= K_ij - 1/3 gamma_ij gamma^kl K_kl
	= TF(K_ij)
*/
sym3 tracefree(sym3 A_ll, sym3 gamma_ll, sym3 gamma_uu) {
	real tr_A = sym3_dot(A_ll, gamma_uu);
	return sym3_sub(A_ll, sym3_real_mul(gamma_ll, tr_A / 3.));
}

//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
kernel void calcDeriv(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

<?=makePartial('alpha', 'real')?>		//partial_alpha_l[i] := alpha_,i
<? if eqn.useChi then ?>
<?=makePartial('chi', 'real')?>			//partial_chi_l[i] := chi_,i 
<?=makePartial('gammaBar_uu', 'sym3')?>	//partial_gammaBar_uul[k].ij = gamma^ij_,k
<? else ?>
<?=makePartial('phi', 'real')?>			//partial_phi_l[i] := phi_,i 
<? end ?>
<?=makePartial('K', 'real')	?>			//partial_K_l[i] := K,i
<?=makePartial('beta_u', 'real3')?>		//partial_beta_ul[j].i := beta^i_,j
<?=makePartial('LambdaBar_u', 'real3')?>	//partial_LambdaBar_ul[j].i := connBar^i_,j
<?=makePartial('epsilon_ll', 'sym3')?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=makePartial('ABar_ll', 'sym3')?>	//partial_ABar_lll[k].ij = ABar_ij,k
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
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
<?=makePartial2('epsilon_ll', 'sym3')?>		//partial2_epsilon_llll[kl].ij = epsilon_ij,kl = gammaBar_ij,kl
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
	
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	
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
<? 
if eqn.useChi then
?>	real3 LambdaBar_u = real3_zero;
	//connBar^i = -gammaBar^ij_,j
<?	for i,xi in ipairs(xNames) do
?>	LambdaBar_u.<?=xi?> =<?
		for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
		end
?>;
<?
	end
else	-- not eqn.useChi
?>	real3 LambdaBar_u = U->LambdaBar_u;
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
			?> - .5 * gammaBar_ll.<?=xij?> * partial_chi_bardot_partial_alpha,
<? 	end
?>	};

<? 	if false then -- why bother when I need to use chi D_i D_j alpha?  why not just take the trace?
?>
	//chi gammaBar^ij alpha,ij - chi connBar^i alpha,i - 1/2 gammaBar^ij chi,i alpha,j
	real tr_D2_alpha = U->chi * (
			sym3_dot(U->gammaBar_uu, *(sym3*)partial2_alpha_ll)
			- real3_dot(U->LambdaBar_u, *(real3*)partial_alpha_l)
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
	sym3 D2_alpha_ll = sym3_zero;
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

	// should be zero, but 2017 Ruchlin et al inserts it into the d/dt epsilon_ij to constrain it to zero
	real tr_ABar_ll = sym3_dot(U->ABar_ll, U->gammaBar_uu);

	//B&S 11.51
	//Alcubierre 2.8.9
	//gammaBar_ij,t = -2/3 gammaBar_ij beta^k_,k - 2 alpha ABar_ij + beta^k gammaBar_ij,k + gammaBar_ik beta^k_,j + gammaBar_kj beta^k_,i
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
?>	deriv->epsilon_ll.<?=xij?> += -2 * U->alpha * U->ABar_ll.<?=xij?>	//-2 alpha ABar_ij 
<? 	for k,xk in ipairs(xNames) do
?>		+ partial_epsilon_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>		//+ beta^k gammaBar_ij,k 
		+ gammaBar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>	//+ gammaBar_jk beta^k_,i 
		+ gammaBar_ll.<?=sym(k,i)?> * partial_beta_ul[<?=j-1?>].<?=xk?> 	//+ gammaBar_ik beta^k_,j
<? 	end
?>
		// TODO 2017 Ruchlin et al has DBar_k beta^k instead of partial_k beta^k
		- 2./3. * gammaBar_ll.<?=xij?> * tr_partial_beta;				//- 2/3 gammaBar_ij beta^k_,k
		+ 2./3. * gammaBar_ll.<?=xij?> * tr_ABar_ll;
<? end
?>
	real3x3 ABar_ul = sym3_sym3_mul(U->gammaBar_uu, U->ABar_ll);		//ABar^i_j = gammaBar^kl ABar_kj
	sym3 ABar_uu = real3x3_sym3_to_sym3_mul(ABar_ul, U->gammaBar_uu);	//ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	real tr_ABar_sq = sym3_dot(U->ABar_ll, ABar_uu);			//tr_ABar_sq := tr(ABar^2) = ABar_ij ABar^ji
	
	real S = sym3_dot(U->S_ll, gamma_uu);
	
	//B&S 11.52
	//Alcubierre 2.8.12
	//K_,t = -gamma^ij D_i D_j alpha + alpha (ABar_ij ABar^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	deriv->K += -tr_D2_alpha
		+ U->alpha * (tr_ABar_sq + U->K * U->K / 3.) 
		+ 4. * M_PI * U->alpha * (U->rho + S) 
		+ real3_dot(U->beta_u, *(real3*)partial_K_l);

	//tr_partial2_gammaBar_ll.ij = gammaBar^kl gammaBar_ij,kl
	sym3 tr_partial2_gammaBar_ll;	
<? for ij,xij in ipairs(symNames) do
?>	tr_partial2_gammaBar_ll.<?=xij?> = 0. <?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?> + U->gammaBar_uu.<?=sym(k,l)?> * partial2_epsilon_llll[<?=from3x3to6(k,l)-1?>].<?=xij?><?
		end
	end
?>;
<? end
?>
	//2017 Ruchlin eqn 12
	// where connHatBar^i_jk = connBar^i_jk raised by gammaBar_ij 
	//RBar_ij = 1/2 (
	//	- gammaBar^kl gammaBar_ij,kl 
	//	+ 2 gammaBar^kl gammaBar_im,k connHat^m_jl
	//	+ 2 gammaBar^kl gammaBar_jm,k connHat^m_il 
	//	+ gammaBar^kl connHat^m_ik,l gammaBar_jm 
	//	+ gammaBar^kl connHat^m_jk,l gammaBar_im
	//	+ gammaBar_ij,k connHatBar^k
	//	- 2 connHatBar_mki connHatBar^mk_j
	//	- connHatBar_kli connHatBar_j^kl 
	//	- connHatBar_ikl connHatBar^kl_j
	//	+ 2 Delta^kl_i Delta_jkl
	//	+ 2 Delta^kl_j Delta_ikl
	//	+ 2 Delta^kl_i Delta_klj
	//	+ gammaBar_ki LambdaBar^k_,j 
	//	+ gammaBar_kj LambdaBar^k_,i 
	//	+ connHatBar_ijk LambdaBar^k
	//	+ connHatBar_jik LambdaBar^k
	//	- connHatBar_ijk connHatBar^k
	//	- connHatBar_jik connHatBar^k
	//	+ Delta_ijk Delta^k
	//	+ Delta_jik Delta^k
	//)
	sym3 RBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	RBar_ll.<?=xij?> = 0.
		-.5 * tr_partial2_gammaBar_ll.<?=xij?>																	//diverging for spherical coordinates
<?	for k,xk in ipairs(xNames) do
?>
		+ .5 * gammaBar_ll.<?=sym(k,i)?> * partial_LambdaBar_ul[<?=j-1?>].<?=xk?>									//diverging for spherical coordinates
		+ .5 * gammaBar_ll.<?=sym(k,j)?> * partial_LambdaBar_ul[<?=i-1?>].<?=xk?>									//diverging for spherical coordinates
		+ .5 * LambdaBar_u.<?=xk?> * (connBar_lll.<?=xi?>.<?=sym(j,k)?> + connBar_lll.<?=xj?>.<?=sym(i,k)?>)		//diverging for spherical coordinates
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>		+ U->gammaBar_uu.<?=sym(k,m)?> * (0.																	//diverging for spherical coordinates
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
?>	RPhi_ll.<?=xij?> = -2. * (0.
		+ DBar2_phi_ll.<?=xij?> 									//diverging for spherical coordinates
		+ U->gammaBar_uu.<?=xij?> * tr_DBar2_phi 					//diverging for spherical coordinates
		- 2. * partial_phi_l[<?=i-1?>] * partial_phi_l[<?=j-1?>] 	//diverging for spherical coordinates
		+ 2. * gammaBar_ll.<?=xij?> * DBar_phi_norm					//diverging for spherical coordinates
	);
<? end 
?>

	sym3 R_ll = sym3_add(RPhi_ll, RBar_ll);							//diverging for spherical coordinates




<? if eqn.useChi then
?>
	//TODO instead of R_ll, add in chi_tf_RPhi and chi_tf_RBar_ll
	//then factor out one of the chis 
	
	//(chi alpha (R_ij - 8 pi S_ij))^TF
	sym3 tf_chi_alpha_R_minus_S = sym3_real_mul(
		sym3_add(
			R_ll,								//diverging for spherical coordinates
			sym3_real_mul(U->S_ll, -8. * M_PI)	//safe
		), 
		U->chi * U->alpha);
	tf_chi_alpha_R_minus_S = tracefree(tf_chi_alpha_R_minus_S, gamma_ll, gamma_uu);

	//(chi D_i D_j alpha)^TF
	//= chi D_i D_j alpha + gammaBar_ij ( 1/3 chi (connBar^k alpha,k - gammaBar^kl alpha,kl) + 1/6 gammaBar^kl chi,k alpha,l)
	sym3 chi_tf_D2_alpha = sym3_add( 
		chi_D2_alpha_ll,
		sym3_real_mul(gammaBar_ll, 
			1./3. * U->chi * (
				real3_dot(U->LambdaBar_u, *(real3*)partial_alpha_l)
				- sym3_dot(U->gammaBar_uu, *(sym3*)partial2_alpha_ll)
			) + 1./6. * real3_weightedDot(
				*(real3*)partial_alpha_l,
				*(real3*)partial_chi_l,
				U->gammaBar_uu
			)
		)
	);
	sym3 chi_tracelessPart_ll = sym3_sub(
		tf_chi_alpha_R_minus_S,				//diverging for spherical coordinates
		chi_tf_D2_alpha						//safe
	);
<? else
?>
	//traceless portion of -D^2 alpha + alpha (R_ij - 8 pi S_ij)
	//...times chi = exp(-4 phi)
#if 1	//all at once
	sym3 tracelessPart_ll = sym3_sub(
		sym3_real_mul(
			sym3_add(R_ll, sym3_real_mul(U->S_ll, -8. * M_PI)), 
			U->alpha),
		D2_alpha_ll);
	tracelessPart_ll = tracefree(tracelessPart_ll, gamma_ll, gamma_uu);
#else	//each term separately
	sym3 tracelessPart_ll = sym3_sub(
		sym3_real_mul(
			sym3_add(
				tracefree(R_ll, gamma_ll, gamma_uu),
				sym3_real_mul(
					tracefree(U->S_ll, gamma_ll, gamma_uu), 
					-8. * M_PI)), 
			U->alpha),
		tracefree(D2_alpha_ll, gamma_ll, gamma_uu)
	);
#endif
	sym3 chi_tracelessPart_ll = sym3_real_mul(tracelessPart_ll, exp_neg4phi);
<? end
?>

	//TODO this is working in cartesian 1D-3D but for spherical 1D this is getting bad values on xx, xy, xz, yy, zz (only not yz)
	//originally in  B&S 11.53, Alcubierre 2.8.11
	//ABar_ij,t = 
	//	exp(-4phi) (-(D_i D_j alpha) + alpha R_ij - 8 alpha pi S_ij )^TF
	
	//	- 2/3 ABar_ij beta^k_,k
	//	- 2 alpha ABar_il ABar^l_j
	//	+ alpha ABar_ij K
	
	//	+ ABar_kj beta^k_,i 
	//	+ ABar_ik beta^k_,j 
	//	+ beta^k ABar_ij,k 
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->ABar_ll.<?=xij?> += 0.
		+ chi_tracelessPart_ll.<?=xij?>											//diverging for spherical coordinates
		+ U->alpha * U->K * U->ABar_ll.<?=xij?> 								//safe
<?	for k,xk in ipairs(xNames) do
?>		- 2. * U->alpha * U->ABar_ll.<?=sym(i,k)?> * ABar_ul.<?=xk?>.<?=xj?>	//safe
		+ partial_ABar_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>				//safe
		+ U->ABar_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>			//safe
		+ U->ABar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>			//safe
<?	end
?>
		- 2./3. * U->ABar_ll.<?=xij?> * tr_partial_beta							//safe
	;
<? end
?>

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
		2./3. * LambdaBar_u.<?=xi?> * tr_partial_beta
		- 16. * M_PI * exp_4phi * U->alpha * U->S_u.<?=xi?> 
<?	for j,xj in ipairs(xNames) do
		local xij = sym(i,j)
		local jj = from3x3to6(j,j)
?>		- 2. * ABar_uu.<?=xij?> * partial_alpha_l[<?=j-1?>]
		+ 2. * U->alpha * (
			-2./3. * U->gammaBar_uu.<?=xij?> * partial_K_l[<?=j-1?>] 
			+ 6. * ABar_uu.<?=xij?> * partial_phi_l[<?=j-1?>])
		+ U->beta_u.<?=xi?> * partial_LambdaBar_ul[<?=j-1?>].<?=xi?>
		- LambdaBar_u.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
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
	deriv->LambdaBar_u = real3_add(deriv->LambdaBar_u, dt_connBar_u);

<? if eqn.useShift == 'GammaDriver' then ?>
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k (connBar^i_,t + eta connBar^i)
	const real k = 0.;//3./4.;
	const real eta = 0.;	// 1 / (2 M), for total mass M
	deriv->beta_u = real3_add(deriv->beta_u,
		real3_add(
			real3_real_mul(dt_connBar_u, k),
			real3_real_mul(LambdaBar_u, eta)));
<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>
	//hyperbolic Gamma driver 
	//B&S 4.83 
	//beta^i_,t = 3/4 B^i (Rezolla says 3/4 alpha B^i, 2006 Campanelli says B^i but absorbs 3/4 into B^i)
	//B^i_,t = connBar^i_,t - eta B^i ... maybe + beta^j B^i_,j (Rezolla says - beta^j B^i_,j)
	const real eta = 0.;
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
}


kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global cons_t* U = UBuf + index;

	sym3 gammaHat_ll = coord_g(x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, U->epsilon_ll);
<? 
if eqn.guiVars.constrain_det_gammaBar_ll.value 
or eqn.guiVars.constrain_tr_ABar_ll.value 
then 
?>
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
<?	if eqn.guiVars.constrain_det_gammaBar_ll.value then ?>
	real det_gammaHat_ll = sqrt_det_g_grid(x);
	real det_gammaBar_ll = sym3_det(gammaBar_ll);
	real rescaleMetric = cbrt(det_gammaHat_ll/det_gammaBar_ll);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_ll.<?=xij?> *= rescaleMetric;
<? 		end ?>
	U->epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
<?	end ?>

	//here we refresh gammaBar_uu
	U->gammaBar_uu = sym3_inv(gammaBar_ll, 1.);
	
	//in Buchman's paper it says he doesn't do this
	//likewise in my own experiences, this can tend A to grow out of control 
<? if eqn.guiVars.constrain_tr_ABar_ll.value then ?>
	U->ABar_ll = tracefree(U->ABar_ll, gammaBar_ll, U->gammaBar_uu);
<? end 
end
?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

<?=makePartial('epsilon_ll', 'sym3')?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=makePartial('ABar_ll', 'sym3')?>		//partial_ABar_lll[k].ij = ABar_ij,k
<?=makePartial('K', 'real')	?>			//partial_K_l[i] := K,i
<? if eqn.useChi then ?>
<?=makePartial('chi', 'real')?>			//partial_chi_l[i] := chi_,i 
<?=makePartial('gammaBar_uu', 'sym3')?>	//partial_gammaBar_uul[k].ij = gamma^ij_,k
<? else ?>
<?=makePartial('phi', 'real')?>			//partial_phi_l[i] := phi_,i 
<? end ?>
<?=makePartial('LambdaBar_u', 'real3')?>	//partial_LambdaBar_ul[j].i := connBar^i_,j

<? if eqn.useChi then ?>
<?=makePartial2('chi', 'real')?>			//partial2_chi_ll.ij := chi_,ij
<? else ?>
<?=makePartial2('phi', 'real')?>			//partial2_phi_ll.ij := phi_,ij
<? end ?>
<?=makePartial2('epsilon_ll', 'sym3')?>		//partial2_epsilon_llll[kl].ij = epsilon_ij,kl = gammaBar_ij,kl

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
<? 
if eqn.useChi then
?>	real3 LambdaBar_u = real3_zero;
	//connBar^i = -gammaBar^ij_,j
<?	for i,xi in ipairs(xNames) do
?>	LambdaBar_u.<?=xi?> =<?
		for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
		end
?>;
<?
	end
else	-- not eqn.useChi
?>	real3 LambdaBar_u = U->LambdaBar_u;
<? 
end
?>

	real exp_neg4phi = calc_exp_neg4phi(U);
	//TODO minimize using these 
	real exp_4phi = 1. / exp_neg4phi;
	
	sym3 gamma_uu = sym3_real_mul(U->gammaBar_uu, exp_neg4phi);
	
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);

	real3x3 ABar_ul = sym3_sym3_mul(U->gammaBar_uu, U->ABar_ll);		//ABar^i_j = gammaBar^kl ABar_kj
	sym3 ABar_uu = real3x3_sym3_to_sym3_mul(ABar_ul, U->gammaBar_uu);	//ABar^ij = gammaBar^ik ABar_kl gammaBar^lj

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
	real3 DBar_phi_u = sym3_real3_mul(U->gammaBar_uu, *(real3*)partial_phi_l);
	
	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	_3sym3 connBar_ull = sym3_3sym3_mul(U->gammaBar_uu, connBar_lll);

	//tr_partial2_gammaBar_ll.ij = gammaBar^kl gammaBar_ij,kl
	sym3 tr_partial2_gammaBar_ll;	
<? for ij,xij in ipairs(symNames) do
?>	tr_partial2_gammaBar_ll.<?=xij?> = 0. <?
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
?>	RBar_ll.<?=xij?> = -.5 * tr_partial2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		+ .5 * gammaBar_ll.<?=sym(k,i)?> * partial_LambdaBar_ul[<?=j-1?>].<?=xk?>
		+ .5 * gammaBar_ll.<?=sym(k,j)?> * partial_LambdaBar_ul[<?=i-1?>].<?=xk?>
		+ .5 * LambdaBar_u.<?=xk?> * (connBar_lll.<?=xi?>.<?=sym(j,k)?> + connBar_lll.<?=xj?>.<?=sym(i,k)?>)
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
			- 1./12. * U->K * U->K
			+ 2. * M_PI * U->rho
		)
	);
#endif
#if 1
	real R = sym3_dot(gamma_uu, R_ll);	//R = gamma^ij R_ij
	
	//K_ij = A_ij + 1/3 gamma_ij K 
	//A_ij = exp(4 phi) ABar_ij
	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_real_mul(U->ABar_ll, exp_4phi),
		sym3_real_mul(gamma_ll, U->K/3.));
	real3x3 K_ul = sym3_sym3_mul(gamma_uu, K_ll);
	sym3 K_uu = real3x3_sym3_to_sym3_mul(K_ul, gamma_uu);
	
	//Alcubierre 3.1.1
	U->H = R + U->K * U->K - sym3_dot(K_uu, K_ll) - 16. * M_PI * U->rho;
#endif
	
#if 1
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
	//M^i = exp(6 phi) (
	//			+ 6 ABar^ij phi_,j
	//			+ connBar^i_jk ABar^jk
	//			- ABar^kj gammaBar^li gammaBar_kl,j
	//			- ABar^ki gammaBar^lj gammaBar_kl,j
	//			+ gammaBar^ik gammaBar^lj ABar_kl,j
	//			- 2/3 K_,j gammaBar^ij 
	//			- 8 pi S^i
	//		)
	real exp_6phi = sqrt(exp_4phi * exp_4phi * exp_4phi);
<? for i,xi in ipairs(xNames) do
?>	U->M_u.<?=xi?> = exp_6phi * (
		- 8. * M_PI * U->S_u.<?=xi?>
<?	for j,xj in ipairs(xNames) do
?>		+ 6. * ABar_uu.<?=sym(i,j)?> * partial_phi_l[<?=j-1?>]
		- 2./3. * exp_6phi * U->gammaBar_uu.<?=sym(i,j)?> * partial_K_l[<?=j-1?>]
<?		for k,xk in ipairs(xNames) do
?>		- connBar_ull.<?=xi?>.<?=sym(j,k)?> * ABar_uu.<?=sym(j,k)?>
<?			for l,xl in ipairs(xNames) do
?>		- ABar_uu.<?=sym(k,j)?> * U->gammaBar_uu.<?=sym(l,i)?> * partial_epsilon_lll[<?=j-1?>].<?=sym(k,l)?>
		- ABar_uu.<?=sym(k,i)?> * U->gammaBar_uu.<?=sym(l,j)?> * partial_epsilon_lll[<?=j-1?>].<?=sym(k,l)?>
		+ U->gammaBar_uu.<?=sym(k,i)?> * U->gammaBar_uu.<?=sym(l,j)?> * partial_ABar_lll[<?=j-1?>].<?=sym(k,l)?>
<?			end
		end
	end
?>	);
<? end ?>
#endif
}
<? end	-- calc_H_and_M ?>

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
<? if false then ?>
	SETBOUNDS_NOGHOST();
	
	global cons_t* deriv = derivBuf + index;

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
		deriv->ptr[i] -= solver->diffuseSigma/16. * lap;
	}
<? end ?>
}
