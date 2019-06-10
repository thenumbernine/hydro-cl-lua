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

<?=makePartial('alpha', 'real')?>			//partial_alpha_l[i] := alpha_,i
<?=makePartial('W', 'real')?>				//partial_W_l[i] := W_,i 
<?=makePartial('K', 'real')	?>				//partial_K_l[i] := K,i
<?=makePartial('beta_u', 'real3')?>			//partial_beta_ul[j].i := beta^i_,j
<?=makePartial('LambdaBar_u', 'real3')?>	//partial_LambdaBar_ul[j].i := connBar^i_,j
<?=makePartial('epsilon_ll', 'sym3')?>		//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=makePartial('ABar_ll', 'sym3')?>			//partial_ABar_lll[k].ij = ABar_ij,k
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=makePartial('B_u', 'real3')?>
<? end ?>
	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

<?=makePartial2('alpha', 'real')?>			//partial2_alpha_ll.ij := alpha_,ij
<?=makePartial2('W', 'real')?>				//partial2_W_ll.ij := W_,ij
<?=makePartial2('epsilon_ll', 'sym3')?>		//partial2_epsilon_llll[kl].ij = epsilon_ij,kl = gammaBar_ij,kl
<?=makePartial2('beta_u', 'real3')?>		//partial2_beta_ull[jk].i = beta^i_,jk

	real exp_neg4phi = calc_exp_neg4phi(U);
	//TODO minimize using these 
	real exp_4phi = 1. / exp_neg4phi;

	real partial_phi_l[3] = {
<? for i,xi in ipairs(xNames) do
?>		-partial_W_l[<?=i-1?>] / (2. * U->W),
<? end
?>	};

	real partial2_phi_ll[6] = {
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.5 * (
			-partial2_W_ll[<?=ij-1?>] 
			+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
		) / U->W,
<? end
?>	};
	
	sym3 gammaHat_uu = coord_gU(x);
	_3sym3 connHat_lll = coord_conn_lll(x);
	_3sym3 connHat_ull = coord_conn_ull(x);
	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar_ll = calc_det_gammaBar_ll(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar_ll);
	_3sym3 partial_gammaBar_lll = _3sym3_add(*(_3sym3*)partial_epsilon_lll, partial_gammaHat_lll);

	//gamma_ij = exp(4 phi) gammaBar_ij
	//only used in K_ll calculation for H constraint
	// exp_4phi could be singular
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);

	//gammaBar_ij = exp(-4 phi) gamma_ij
	//gammaBar^ij = exp(4 phi) gamma^ij
	//gamma^ij = exp(-4 phi) gammaBar^ij
	sym3 gamma_uu = sym3_real_mul(gammaBar_uu, exp_neg4phi);

	//connBar_lll[i].jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll.<?=xi?>.<?=xjk?> = .5 * (
			partial_epsilon_lll[<?=k-1?>].<?=sym(i,j)?>
			+ partial_epsilon_lll[<?=j-1?>].<?=sym(i,k)?> 
			- partial_epsilon_lll[<?=i-1?>].<?=xjk?>
		) + connHat_lll.<?=xi?>.<?=xjk?>;
<?	end
end
?>	
	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	_3sym3 connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);

	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);

	real3 Delta_u = (real3){
<? for i,xi in ipairs(xNames) do 
?>		.<?=xi?> = sym3_dot(gammaBar_uu, Delta_ull.<?=xi?>),
<? end
?>	};


	/*
	partial_connHat_ulll[l].i.jk = connHat^i_jk,l
	= (gammaHat^im connHat_mjk),l
	= gammaHat^im_,l connHat_mjk + gammaHat^im connHat_mjk,l
	= 
		gammaHat^im (
			1/2 ( 
				gammaHat_mj,kl
				+ gammaHat_mk,jl
				- gammaHat_kj,ml
			)
			- gammaHat_mn,l connHat^n_jk 
		)
	*/
	_3sym3 partial_connHat_ulll[3];
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
		for l,xl in ipairs(xNames) do
?>	partial_connHat_ulll[<?=l-1?>].<?=xi?>.<?=xjk?> = 0.
<?			for m,xm in ipairs(xNames) do
?>		+ gammaHat_uu.<?=sym(i,m)?> * (
			.5 * (
				partial2_gammaHat_llll.<?=sym(m,j)?>.<?=sym(k,l)?>
				+ partial2_gammaHat_llll.<?=sym(m,k)?>.<?=sym(j,l)?>
				- partial2_gammaHat_llll.<?=sym(k,j)?>.<?=sym(m,l)?>
			)
<?				for n,xn in ipairs(xNames) do
?>			- partial_gammaHat_lll.<?=xl?>.<?=sym(m,n)?> * connHat_ull.<?=xn?>.<?=sym(j,k)?>
<?				end
?>		
		)
<?			end
?>	;
<?
		end
	end
end
?>

	_3sym3 partial_connBar_ulll[3];
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
		for l,xl in ipairs(xNames) do
?>	partial_connBar_ulll[<?=l-1?>].<?=xi?>.<?=xjk?> = 0.
<?			for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(i,m)?> * (
			.5 * (
				partial2_epsilon_llll[<?=from3x3to6(m,j)-1?>].<?=sym(k,l)?>
				+ partial2_epsilon_llll[<?=from3x3to6(m,k)-1?>].<?=sym(j,l)?>
				- partial2_epsilon_llll[<?=from3x3to6(k,j)-1?>].<?=sym(m,l)?>
				
				+ partial2_gammaHat_llll.<?=sym(m,j)?>.<?=sym(k,l)?>
				+ partial2_gammaHat_llll.<?=sym(m,k)?>.<?=sym(j,l)?>
				- partial2_gammaHat_llll.<?=sym(k,j)?>.<?=sym(m,l)?>
			)
<?				for n,xn in ipairs(xNames) do
?>			- partial_gammaBar_lll.<?=xl?>.<?=sym(m,n)?> * connBar_ull.<?=xn?>.<?=sym(j,k)?>
<?				end
?>		
		)
<?			end
?>	;
<?
		end
	end
end
?>


	//DBar^i phi = gammaBar^ij phi_,j
	real3 DBar_phi_u = sym3_real3_mul(gammaBar_uu, *(real3*)partial_phi_l);


	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>		.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]<?
	for k,xk in ipairs(xNames) do
		?> - connBar_ull.<?=xk?>.<?=xij?> * partial_alpha_l[<?=k-1?>]<?
	end
		?>,
<? end
?>	};

	real tr_DBar2_alpha = sym3_dot(gammaBar_uu, DBar2_alpha_ll);


	real3 tr13_connBar_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
?> + connBar_ull.<?=xj?>.<?=sym(i,j)?><?
	end
?>,
<? end
?>	};


	//Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	real Q = calc_f(U->alpha) * U->K;
	
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q + real3_dot(*(real3*)partial_alpha_l, U->beta_u);


	//connBar^j_kj beta^k
	real tr13_connBar_beta = real3_dot(tr13_connBar_l, U->beta_u);

	//2018 Ruchlin et al eqn 11c
	//W,t = 1/3 W (alpha K - beta^k connBar^j_kj - beta^k_,k) + beta^k W_,k
	deriv->W += (1. / 3.) * U->W * (U->alpha * U->K - tr13_connBar_beta - tr_partial_beta) + real3_dot(U->beta_u, *(real3*)partial_W_l);



	// should be zero, but 2018 Ruchlin et al inserts it into the d/dt epsilon_ij to constrain it to zero
	real tr_ABar_ll = sym3_dot(U->ABar_ll, gammaBar_uu);

	/*
	2018 Ruchlin et al, eqn 11a
	epsilon_ij,t = 2/3 gammaBar_ij (alpha ABar^k_k - DBar_k beta^k) + DHat_i beta_j + DHat_j beta_i - 2 alpha ABar_ij + epsilon_ij,k beta^k + epsilon_ik beta^k_,j + epsilon_kj beta^k_,i
	= 	
		- 2 alpha ABar_ij 
		+ 2/3 gammaBar_ij (
			alpha ABar^k_k 
			- beta^k_,k 
			- connBar^k_lk beta^l
		) 
		+ gammaBar_kj beta^k_,i
		+ gammaBar_ki beta^k_,j 
		+ epsilon_ij,k beta^k 
		+ gammaHat_ij,k beta^k 
	*/
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->epsilon_ll.<?=xij?> += 
		- 2 * U->alpha * U->ABar_ll.<?=xij?>
		+ 2. / 3. * gammaBar_ll.<?=xij?> * (
			U->alpha * tr_ABar_ll
			- tr_partial_beta
			- tr13_connBar_beta
		)
<? 	for k,xk in ipairs(xNames) do
?>		
		+ gammaBar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
		+ gammaBar_ll.<?=sym(k,i)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
		+ (
			// connHat_jki + connHat_ikj = gammaHat_ij,k
			connHat_lll.<?=xj?>.<?=sym(i,k)?>
			+ connHat_lll.<?=xi?>.<?=sym(j,k)?>
			// gammaHat_ij,k + epsilon_ij,k = gammaBar_ij,k
			+ partial_epsilon_lll[<?=k-1?>].<?=xij?>
		) * U->beta_u.<?=xk?>
<? 	end
?>	;
<? end
?>
	
	
	real3x3 ABar_ul = sym3_sym3_mul(gammaBar_uu, U->ABar_ll);		//ABar^i_j = gammaBar^kl ABar_kj
	sym3 ABar_uu = real3x3_sym3_to_sym3_mul(ABar_ul, gammaBar_uu);	//ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	real tr_ABar_sq = sym3_dot(U->ABar_ll, ABar_uu);				//tr_ABar_sq := tr(ABar^2) = ABar_ij ABar^ji
	
	real S = sym3_dot(U->S_ll, gamma_uu);


	/*
	B&S 11.52
	Alcubierre 2.8.12
	K_,t = -gamma^ij D_i D_j alpha + alpha (ABar_ij ABar^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	Ruchlin et al 2018
	K_,t = 
		1/3 alpha K^2 
		+ alpha ABar_ij ABar^ij 
		- exp(-4 phi) (
			DBar_i DBar^i alpha 
			+ 2 gammaBar^ij alpha_,i phi_,j
		) 
		+ K_,i beta^i
		+ 4 pi alpha (rho + S)
	*/
	deriv->K += 
		+ U->alpha * U->K * U->K / 3.
		+ U->alpha * tr_ABar_sq
		- exp_neg4phi * (
			tr_DBar2_alpha
			+ real3_weightedDot(
				*(real3*)partial_phi_l,
				*(real3*)partial_alpha_l,
				gammaBar_uu
			)
		)
		+ real3_dot(U->beta_u, *(real3*)partial_K_l)
		+ 4. * M_PI * U->alpha * (U->rho + S) 
	;


	//tr_partial2_gammaBar_ll.ij = gammaBar^kl gammaBar_ij,kl
	sym3 tr_partial2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
?>	tr_partial2_gammaBar_ll.<?=xij?> = 0. <?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?> + gammaBar_uu.<?=sym(k,l)?> * partial2_epsilon_llll[<?=from3x3to6(k,l)-1?>].<?=xij?><?
		end
	end
?>;
<? end
?>

	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);

	//DHat_LambdaBar_ul[i].xj = DHat_i LambdaBar^j = LambdaBar^j_,i + connHat^j_ki LambdaBar^k
	real3x3 DHat_LambdaBar_ul = (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = partial_LambdaBar_ul[<?=i-1?>].<?=xj?><?
		for k,xk in ipairs(xNames) do
			?> + connHat_ull.<?=xj?>.<?=sym(i,k)?> * U->LambdaBar_u.<?=xk?><?
		end	?>,
<?	end
?>		},
<? end
?>	};


	//gammaHat^ik gammaBar_kj 
	real3x3 gammaHat_to_gammaBar_ul = (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
?> 				+ gammaHat_uu.<?=sym(i,k)?> * gammaBar_ll.<?=sym(k,j)?>
<?		end ?>,
<?	end
?>		},
<? end
?>	};

	/*
	tr_DHat2_gammaBar_ll.ij =
	gammaBar^kl DHat_k DHat_l gammaBar_ij
	= gammaBar^kl DHat_k (gammaBar_ij,l - connHat^m_il gammaBar_mj - connHat^m_jl gammaBar_im)
	= gammaBar^kl DHat_k (
		gammaBar_ij,l 
		- connHat^m_il gammaBar_mj 
		- connHat^m_jl gammaBar_im
	)
	= gammaBar^kl (
		(
			gammaBar_ij,l 
			- connHat^m_il gammaBar_mj 
			- connHat^m_jl gammaBar_im
		),k 
		- connHat^n_ik (
			gammaBar_nj,l 
			- connHat^m_nl gammaBar_mj 
			- connHat^m_jl gammaBar_nm
		)
		- connHat^n_jk (
			gammaBar_in,l
			- connHat^m_il gammaBar_mn 
			- connHat^m_nl gammaBar_im
		)
		- connHat^n_lk (
			gammaBar_ij,n 
			- connHat^m_in gammaBar_mj 
			- connHat^m_jn gammaBar_im
		)
	)
	= gammaBar^kl (
		gammaBar_ij,lk 
		
		- connHat^m_il,k gammaBar_mj 
		- connHat^m_il gammaBar_mj,k 
		
		- connHat^m_jl,k gammaBar_im
		- connHat^m_jl gammaBar_im,k
		
		- connHat^n_ik (
			gammaBar_nj,l 
			- connHat^m_nl gammaBar_mj 
			- connHat^m_jl gammaBar_nm
		)
		- connHat^n_jk (
			gammaBar_in,l
			- connHat^m_il gammaBar_mn 
			- connHat^m_nl gammaBar_im
		)
		- connHat^n_lk (
			gammaBar_ij,n 
			- connHat^m_in gammaBar_mj 
			- connHat^m_jn gammaBar_im
		)
	)
	*/

	sym3 tr_DHat2_gammaBar_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.<?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do 
			local kl = from3x3to6(k,l)
?> 			+ gammaBar_uu.<?=sym(k,l)?> * (
				partial2_epsilon_llll[<?=kl-1?>].<?=xij?>
				+ partial2_gammaHat_llll.<?=sym(k,l)?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>				
				- partial_connHat_ulll[<?=k-1?>].<?=xm?>.<?=sym(i,l)?> * gammaBar_ll.<?=sym(m,j)?>
				- connHat_ull.<?=xm?>.<?=sym(i,l)?> * partial_gammaBar_lll.<?=xk?>.<?=sym(j,m)?>
				
				- partial_connHat_ulll[<?=k-1?>].<?=xm?>.<?=sym(j,l)?> * gammaBar_ll.<?=sym(m,i)?>
				- connHat_ull.<?=xm?>.<?=sym(j,l)?> * partial_gammaBar_lll.<?=xk?>.<?=sym(i,m)?>

<?				for n,xn in ipairs(xNames) do
?>
				- connHat_ull.<?=xn?>.<?=sym(i,k)?> * (
					partial_gammaBar_lll.<?=xl?>.<?=sym(n,j)?>
					- connHat_ull.<?=xm?>.<?=sym(n,l)?> * gammaBar_ll.<?=sym(m,j)?>
					- connHat_ull.<?=xm?>.<?=sym(j,l)?> * gammaBar_ll.<?=sym(n,m)?>
				)
				- connHat_ull.<?=xn?>.<?=sym(j,k)?> * (
					partial_gammaBar_lll.<?=xl?>.<?=sym(i,n)?>
					- connHat_ull.<?=xm?>.<?=sym(i,l)?> * gammaBar_ll.<?=sym(m,n)?>
					- connHat_ull.<?=xm?>.<?=sym(n,l)?> * gammaBar_ll.<?=sym(i,m)?>
				)
				- connHat_ull.<?=xn?>.<?=sym(l,k)?> * (
					partial_gammaBar_lll.<?=xn?>.<?=sym(i,j)?>
					- connHat_ull.<?=xm?>.<?=sym(i,n)?> * gammaBar_ll.<?=sym(m,j)?>
					- connHat_ull.<?=xm?>.<?=sym(j,n)?> * gammaBar_ll.<?=sym(i,m)?>
				)
<?				end
			end
?>			)
<?		end
	end	?>,
<? end
?>	};

	/*
	2018 Ruchlin eqn 12
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ijk
		+ 1/2 Delta^k Delta_jik
		+ gammaBar^kl (
			Delta^m_ki Delta_jml
			+ Delta^m_kj Delta_iml
			+ Delta^m_ik Delta_mjl
		)
	*/
	sym3 RBar_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.
			-.5 * tr_DHat2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * gammaBar_ll.<?=sym(k,i)?> * DHat_LambdaBar_ul.<?=xj?>.<?=xk?>
			+ .5 * gammaBar_ll.<?=sym(k,j)?> * DHat_LambdaBar_ul.<?=xi?>.<?=xk?>
			+ .5 * Delta_u.<?=xk?> * (Delta_lll.<?=xi?>.<?=sym(j,k)?> + Delta_lll.<?=xj?>.<?=sym(i,k)?>)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>				+ gammaBar_uu.<?=sym(k,l)?> * (
				Delta_ull.<?=xm?>.<?=sym(k,i)?> * Delta_lll.<?=xj?>.<?=sym(m,l)?>
				+ Delta_ull.<?=xm?>.<?=sym(k,j)?> * Delta_lll.<?=xi?>.<?=sym(m,l)?>
				+ Delta_ull.<?=xm?>.<?=sym(i,k)?> * Delta_lll.<?=xm?>.<?=sym(j,l)?>
			)<?
			end
		end
	end
?>,
<? end
?>	};


	sym3 DBar2_phi_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = partial2_phi_ll[<?=ij-1?>]<?
	for k,xk in ipairs(xNames) do
		?> - connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l[<?=k-1?>]<?
	end
?>,
<? end
?>	};

	real tr_DBar_beta = tr_partial_beta + real3_dot(tr13_connBar_l, U->beta_u);

	/*
	2018 Ruchlin et al eqn 11b
	traceless portion of ...
	exp(-4 phi) (
		-2 alpha DBar_i DBar_j phi 
		+ 4 alpha DBar_i phi DBar_j phi 
		+ 2 DBar_i phi DBar_j alpha 
		+ 2 DBar_i alpha DBar_j phi 
		- DBar_i DBar_j alpha 
		+ alpha RBar_ij 
		- 8 pi alpha S_ij
	)
	*/
	sym3 exp_neg4phi_tracelessPart_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = exp_neg4phi * (
			2. * partial_phi_l[<?=i-1?>] * partial_alpha_l[<?=j-1?>]
			+ 2. * partial_phi_l[<?=j-1?>] * partial_alpha_l[<?=i-1?>]
			- DBar2_alpha_ll.<?=xij?>
			+ U->alpha * (
				- 2. * DBar2_phi_ll.<?=xij?>
				+ 4. * partial_phi_l[<?=i-1?>] * partial_phi_l[<?=j-1?>]
				+ RBar_ll.<?=xij?>
				- 8. * M_PI * U->S_ll.<?=xij?>
			)
		),
<? end
?>	};
	
	exp_neg4phi_tracelessPart_ll = tracefree(exp_neg4phi_tracelessPart_ll, gammaBar_ll, gammaBar_uu);

	/*
	2018 Ruchlin et al, eqn. 11b
	ABar_ij,t = 
		+ exp(-4 phi) (trace-free part above)_ij
		+ alpha K ABar_ij
		- 2/3 ABar_ij DBar_k beta^k
		- 2 alpha ABar_ik ABar^k_j
		+ ABar_ij,k beta^k
		+ ABar_kj beta^k_,i
		+ ABar_ki beta^k_,j
	*/
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi = xNames[i]
	local xj = xNames[j]
?>	deriv->ABar_ll.<?=xij?> += 
		exp_neg4phi_tracelessPart_ll.<?=xij?>
		+ U->alpha * U->K * U->ABar_ll.<?=xij?>
		- 2. / 3. * U->ABar_ll.<?=xij?> * tr_DBar_beta
<?	for k,xk in ipairs(xNames) do
?>		- 2. * U->alpha * U->ABar_ll.<?=sym(i,k)?> * ABar_ul.<?=xk?>.<?=xj?>
		+ partial_ABar_lll[<?=k-1?>].<?=xij?> * U->beta_u.<?=xk?>
		+ U->ABar_ll.<?=sym(i,k)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
		+ U->ABar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
<?	end
?>		- 2./3. * U->ABar_ll.<?=xij?> * tr_partial_beta;
<? end
?>


	/*
	DHat2_beta_u.i.jk = DHat_k DHat_j beta^i
	= DHat_k (beta^i_,j + connHat^i_lj beta^l)
	= beta^i_,jk 
		+ connHat^i_lj,k beta^l 
		+ connHat^i_lj beta^l_,k
	*/
	_3sym3 DHat2_beta_u;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
?>	DHat2_beta_u.<?=xi?>.<?=xjk?> = 
		partial2_beta_ull[<?=jk-1?>].<?=xi?>
<?		for l,xl in ipairs(xNames) do
?>		+ partial_connHat_ulll[<?=k-1?>].<?=xi?>.<?=sym(l,j)?> * U->beta_u.<?=xl?>
		+ connHat_ull.<?=xi?>.<?=sym(l,j)?> * partial_beta_ul[<?=k-1?>].<?=xl?>
<?		end
?>	;
<?	end 
end
?>

	/*
	tr_gammaBar_DHat2_beta_u.i = gammaBar^jk DHat_j DHat_k beta^i
	*/
	real3 tr_gammaBar_DHat2_beta_u = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_dot(gammaBar_uu, DHat2_beta_u.<?=xi?>),
<? end
?>	};

	/*
	DBar_tr_DBar_beta_u.i = DBar^i DBar_j beta^j
	= gammaBar^ik (beta^j_,j + connBar^j_lj beta^l)_,k
	= gammaBar^ij (
		beta^k_,kj 
		+ connBar^k_lk,j beta^l 
		+ connBar^k_lk beta^l_,j
	)
	*/
	real3 DBar_tr_DBar_beta_u = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
			local jk = from3x3to6(j,k)
?>			+ gammaBar_uu.<?=sym(i,j)?> * (
				partial2_beta_ull[<?=jk-1?>].<?=xk?>
<?			for l,xl in ipairs(xNames) do
?>				+ partial_connBar_ulll[<?=j-1?>].<?=xk?>.<?=sym(l,k)?> * U->beta_u.<?=xl?>
				+ connBar_ull.<?=xk?>.<?=sym(l,k)?> * partial_beta_ul[<?=j-1?>].<?=xl?>
<?			end
?>			)
<?		end
	end
?>		,
<? end
?>	}; 

	/*
	LambdaBar^i_,t = 
		gammaBar^jk DHat_j DHat_k beta^i
		+ 2/3 Delta^i DBar_j beta^j
		+ 1/3 DBar^i DBar_j beta^j
		- 16 pi exp(4 phi) alpha S^i
		- 2 ABar^ij (alpha_,j - 6 phi_,j)
		+ 2 ABar^jk Delta^i_jk
		
		- 4/3 alpha gammaBar^ij K_,j
			// Lie derivative terms
		+ LambdaBar^i_,k beta^k
		- LambdaBar^k beta^i_,k
	*/
	real3 dt_LambdaBar_u;
<? for i,xi in ipairs(xNames) do
?>	dt_LambdaBar_u.<?=xi?> =
		tr_gammaBar_DHat2_beta_u.<?=xi?>
		+ 2. / 3. * Delta_u.<?=xi?> * tr_DBar_beta
		+ 1. / 3. * DBar_tr_DBar_beta_u.<?=xi?>
		- 16. * M_PI * exp_4phi * U->alpha * U->S_u.<?=xi?> 
<?	for j,xj in ipairs(xNames) do
		local xij = sym(i,j)
		local jj = from3x3to6(j,j)
?>		- 2. * ABar_uu.<?=xij?> * (
			partial_alpha_l[<?=j-1?>]
			- 6. * partial_phi_l[<?=j-1?>]
		)
		- 4. / 3. * U->alpha * gammaBar_uu.<?=xij?> * partial_K_l[<?=j-1?>] 
		+ U->beta_u.<?=xi?> * partial_LambdaBar_ul[<?=j-1?>].<?=xi?>
		- U->LambdaBar_u.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?		for k,xk in ipairs(xNames) do		
			local xik = sym(i,k)
			local jk = from3x3to6(j,k)
			local xjk = symNames[jk]
?>		+ 2. * Delta_ull.<?=xi?>.<?=xjk?> * ABar_uu.<?=xjk?>
<?		end
	end
?>	;
<? end
?>
	deriv->LambdaBar_u = real3_add(deriv->LambdaBar_u, dt_LambdaBar_u);


<? if eqn.useShift == 'GammaDriver' then ?>
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k (connBar^i_,t + eta connBar^i)
	const real k = 3./4.;
	const real eta = 1.;	// 1 / (2 M), for total mass M
	deriv->beta_u = real3_add(deriv->beta_u,
		real3_add(
			real3_real_mul(dt_LambdaBar_u, k),
			real3_real_mul(U->LambdaBar_u, eta)));
<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>
	//hyperbolic Gamma driver 
	//2018 Ruchlin et al, eqn 14a, 14b
	//beta^i_,t = B^i
	//B^i_,t = 3/4 (LambdaBar^i_,t - LambdaBar^i_,j beta^j + LambdaBar^j beta^i_,j) - eta B^i + B^i_,j beta^j - B^j beta^i_,j
	deriv->beta_u = real3_add(deriv->beta_u, U->B_u);

	const real eta = 0.;
<?	for i,xi in ipairs(xNames) do
?>	deriv->B_u.<?=xi?> += .75 * (
			dt_LambdaBar_u.<?=xi?>
<? 		for j,xj in ipairs(xNames) do
?>			- partial_LambdaBar_ul[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
			+ partial_beta_ul[<?=j-1?>].<?=xi?> * U->LambdaBar_u.<?=xj?>
<?		end
?>		)
		- eta * U->B_u.<?=xi?>
<?		for j,xj in ipairs(xNames) do
?>		+ partial_B_ul[<?=j-1?>].<?=xi?> * U->beta_u.<?=xj?>
		- partial_beta_ul[<?=j-1?>].<?=xi?> * U->B_u.<?=xj?>
<?		end
?>	;
<? 	end 
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
	real det_gammaHat_ll = coord_det_g(x);
	real det_gammaBar_ll = sym3_det(gammaBar_ll);
	real rescaleMetric = cbrt(det_gammaHat_ll/det_gammaBar_ll);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_ll.<?=xij?> *= rescaleMetric;
<? 		end ?>
	U->epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
<?	end ?>

	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaHat_ll);
	
	//in Buchman's paper it says he doesn't do this
	//likewise in my own experiences, this can tend A to grow out of control 
<? if eqn.guiVars.constrain_tr_ABar_ll.value then ?>
	U->ABar_ll = tracefree(U->ABar_ll, gammaBar_ll, gammaBar_uu);
<? end 
else
?>
	real det_gammaHat_ll = coord_det_g(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaHat_ll);
<?
end
?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

<?=makePartial('epsilon_ll', 'sym3')?>		//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=makePartial('ABar_ll', 'sym3')?>			//partial_ABar_lll[k].ij = ABar_ij,k
<?=makePartial('K', 'real')	?>				//partial_K_l[i] := K,i
<?=makePartial('W', 'real')?>				//partial_W_l[i] := phi_,i 
<?=makePartial('LambdaBar_u', 'real3')?>	//partial_LambdaBar_ul[j].i := connBar^i_,j

<?=makePartial2('W', 'real')?>				//partial2_W_ll.ij := phi_,ij
<?=makePartial2('epsilon_ll', 'sym3')?>		//partial2_epsilon_llll[kl].ij = epsilon_ij,kl = gammaBar_ij,kl

	real partial_phi_l[3] = {
<? for i,xi in ipairs(xNames) do
?>		-partial_W_l[<?=i-1?>] / (2. * U->W),
<? end
?>	};

	real partial2_phi_ll[6] = {
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.5 * (
			-partial2_W_ll[<?=ij-1?>] 
			+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
		) / U->W,
<? end
?>	};
	
	real exp_neg4phi = calc_exp_neg4phi(U);
	//TODO minimize using these 
	real exp_4phi = 1. / exp_neg4phi;
	
	sym3 gamma_uu = sym3_real_mul(gammaBar_uu, exp_neg4phi);	
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	
	_3sym3 connHat_lll = coord_conn_lll(x);
	
	//connBar_lll[i].jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll.<?=xi?>.<?=xjk?> = .5 * (
			partial_epsilon_lll[<?=k-1?>].<?=sym(i,j)?> 
			+ partial_epsilon_lll[<?=j-1?>].<?=sym(i,k)?> 
			- partial_epsilon_lll[<?=i-1?>].<?=xjk?>
		) + connHat_lll.<?=xi?>.<?=xjk?>;
<?	end
end
?>	

	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	_3sym3 connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);


	real3x3 ABar_ul = sym3_sym3_mul(gammaBar_uu, U->ABar_ll);		//ABar^i_j = gammaBar^kl ABar_kj
	sym3 ABar_uu = real3x3_sym3_to_sym3_mul(ABar_ul, gammaBar_uu);	//ABar^ij = gammaBar^ik ABar_kl gammaBar^lj


	_3sym3 connHat_ull = coord_conn_ull(x);


	//tr_partial2_gammaBar_ll.ij = gammaBar^kl gammaBar_ij,kl
	sym3 tr_partial2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
?>	tr_partial2_gammaBar_ll.<?=xij?> = 0. <?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?> + gammaBar_uu.<?=sym(k,l)?> * partial2_epsilon_llll[<?=from3x3to6(k,l)-1?>].<?=xij?><?
		end
	end
?>;
<? end
?>

	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);

	real3 Delta_u = (real3){
<? for i,xi in ipairs(xNames) do 
?>		.<?=xi?> = sym3_dot(gammaBar_uu, Delta_ull.<?=xi?>),
<? end
?>	};

	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);

	//DHat_LambdaBar_ul[i].xj = DHat_i LambdaBar^j = LambdaBar^j_,i + connHat^j_ki LambdaBar^k
	real3x3 DHat_LambdaBar_ul = (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = partial_LambdaBar_ul[<?=i-1?>].<?=xj?><?
		for k,xk in ipairs(xNames) do
			?> + connHat_ull.<?=xj?>.<?=sym(i,k)?> * U->LambdaBar_u.<?=xk?><?
		end	?>,
<?	end
?>		},
<? end
?>	};


	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	sym3 gammaHat_uu = coord_gU(x);

	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	_3sym3 partial_gammaBar_lll = _3sym3_add(*(_3sym3*)partial_epsilon_lll, partial_gammaHat_lll);

	//gammaHat^ik gammaBar_kj 
	real3x3 gammaHat_to_gammaBar_ul = (real3x3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (real3){
<?	for j,xj in ipairs(xNames) do
?>			.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
?> 				+ gammaHat_uu.<?=sym(i,k)?> * gammaBar_ll.<?=sym(k,j)?>
<?		end ?>,
<?	end
?>		},
<? end
?>	};

	/*
	tr_DHat2_gammaBar_ll.ij =
	gammaBar^kl DHat_k DHat_l gammaBar_ij
	= gammaBar^kl DHat_k (gammaBar_ij,l - connHat^m_il gammaBar_mj - connHat^m_jl gammaBar_im)
	= gammaBar^kl DHat_k (
		gammaBar_ij,l 
		- connHat^m_il gammaBar_mj 
		- connHat^m_jl gammaBar_im
	)
	= gammaBar^kl (
		(
			gammaBar_ij,l 
			- connHat^m_il gammaBar_mj 
			- connHat^m_jl gammaBar_im
		),k 
		- connHat^n_ik (
			gammaBar_nj,l 
			- connHat^m_nl gammaBar_mj 
			- connHat^m_jl gammaBar_nm
		)
		- connHat^n_jk (
			gammaBar_in,lk
			- connHat^m_il gammaBar_mn 
			- connHat^m_nl gammaBar_im
		)
		- connHat^n_lk (
			gammaBar_ij,n 
			- connHat^m_in gammaBar_mj 
			- connHat^m_jn gammaBar_im
		)
	)
	= gammaBar^kl (
		epsilon_ij,lk
		gammaHat_ij,lk
		
		- 1/2 (gammaBar_iq gammaHat^qm) (gammaHat_mj,lk + gammaHat_ml,jk - gammaHat_jl,mk)
		- 1/2 (gammaBar_jq gammaHat^qm) (gammaHat_mi,lk + gammaHat_ml,ik - gammaHat_il,mk) 
		
		+ (gammaBar_iq gammaHat^qm) gammaHat_mn,k connHat^n_jl 
		+ (gammaBar_jq gammaHat^qm) gammaHat_mn,k connHat^n_il 
		
		- 2 connHat^m_ik gammaBar_mj,l 
		- 2 connHat^m_jk gammaBar_im,l
		- connHat^m_lk gammaBar_ij,m 
		
		+ connHat^n_ik connHat^m_nl gammaBar_mj 
		+ connHat^n_ik connHat^m_jl gammaBar_nm
		
		+ connHat^n_jk connHat^m_il gammaBar_mn 
		+ connHat^n_jk connHat^m_nl gammaBar_im
		
		+ connHat^n_lk connHat^m_in gammaBar_mj 
		+ connHat^n_lk connHat^m_jn gammaBar_im
	)
	*/

	sym3 tr_DHat2_gammaBar_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.<?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do 
			local kl = from3x3to6(k,l)
?> 
			- gammaBar_uu.<?=sym(k,l)?> * (
				partial2_epsilon_llll[<?=kl-1?>].<?=xij?>
				+ partial2_gammaHat_llll.<?=sym(k,l)?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>				
				- .5 * gammaHat_to_gammaBar_ul.<?=xm?>.<?=xi?> * (
					partial2_gammaHat_llll.<?=sym(m,j)?>.<?=sym(k,l)?>
					+ partial2_gammaHat_llll.<?=sym(m,l)?>.<?=sym(j,k)?>
					- partial2_gammaHat_llll.<?=sym(j,l)?>.<?=sym(m,k)?>
				)
				
				- .5 * gammaHat_to_gammaBar_ul.<?=xm?>.<?=xj?> * (
					partial2_gammaHat_llll.<?=sym(m,i)?>.<?=sym(k,l)?>
					+ partial2_gammaHat_llll.<?=sym(m,l)?>.<?=sym(i,k)?>
					- partial2_gammaHat_llll.<?=sym(i,l)?>.<?=sym(m,k)?>
				)
				
				+ gammaHat_to_gammaBar_ul.<?=xm?>.<?=xi?> * (
<?				for n,xn in ipairs(xNames) do
?>					+ partial_gammaHat_lll.<?=xk?>.<?=sym(m,n)?> * connHat_ull.<?=xn?>.<?=sym(j,l)?>
<?				end
?>				)
				
				+ gammaHat_to_gammaBar_ul.<?=xm?>.<?=xj?> * (
<?				for n,xn in ipairs(xNames) do
?>					+ partial_gammaHat_lll.<?=xk?>.<?=sym(m,n)?> * connHat_ull.<?=xn?>.<?=sym(i,l)?>
<?				end
?>				)
				
				- 2. * connHat_ull.<?=xm?>.<?=sym(i,l)?> * partial_gammaBar_lll.<?=xk?>.<?=sym(m,j)?>
				- 2. * connHat_ull.<?=xm?>.<?=sym(j,l)?> * partial_gammaBar_lll.<?=xk?>.<?=sym(m,i)?>
				- connHat_ull.<?=xm?>.<?=sym(k,l)?> * partial_gammaBar_lll.<?=xm?>.<?=sym(i,j)?>

<?				for n,xn in ipairs(xNames) do
?>				
				+ connHat_ull.<?=xn?>.<?=sym(i,k)?> * connHat_ull.<?=xm?>.<?=sym(n,l)?> * gammaBar_ll.<?=sym(m,j)?>
				+ connHat_ull.<?=xn?>.<?=sym(i,k)?> * connHat_ull.<?=xm?>.<?=sym(j,l)?> * gammaBar_ll.<?=sym(m,n)?>
				
				+ connHat_ull.<?=xn?>.<?=sym(j,k)?> * connHat_ull.<?=xm?>.<?=sym(i,l)?> * gammaBar_ll.<?=sym(m,n)?>
				+ connHat_ull.<?=xn?>.<?=sym(j,k)?> * connHat_ull.<?=xm?>.<?=sym(n,l)?> * gammaBar_ll.<?=sym(i,n)?>
				
				+ connHat_ull.<?=xn?>.<?=sym(l,k)?> * connHat_ull.<?=xm?>.<?=sym(i,n)?> * gammaBar_ll.<?=sym(m,j)?>
				+ connHat_ull.<?=xn?>.<?=sym(l,k)?> * connHat_ull.<?=xm?>.<?=sym(j,l)?> * gammaBar_ll.<?=sym(m,i)?>

<?				end
?>

<?			end
?>
			
			)
<?		end
	end	?>,
<? end
?>	};

	/*
	2018 Ruchlin eqn 12
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ijk
		+ 1/2 Delta^k Delta_jik
		+ gammaBar^kl (
			Delta^m_ki Delta_jml
			+ Delta^m_kj Delta_iml
			+ Delta^m_ik Delta_mjl
		)
	*/
	sym3 RBar_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 0.
			-.5 * tr_DHat2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * gammaBar_ll.<?=sym(k,i)?> * DHat_LambdaBar_ul.<?=xj?>.<?=xk?>
			+ .5 * gammaBar_ll.<?=sym(k,j)?> * DHat_LambdaBar_ul.<?=xi?>.<?=xk?>
			+ .5 * Delta_u.<?=xk?> * (Delta_lll.<?=xi?>.<?=sym(j,k)?> + Delta_lll.<?=xj?>.<?=sym(i,k)?>)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>				+ gammaBar_uu.<?=sym(k,l)?> * (
				Delta_ull.<?=xm?>.<?=sym(k,i)?> * Delta_lll.<?=xj?>.<?=sym(m,l)?>
				+ Delta_ull.<?=xm?>.<?=sym(k,j)?> * Delta_lll.<?=xi?>.<?=sym(m,l)?>
				+ Delta_ull.<?=xm?>.<?=sym(i,k)?> * Delta_lll.<?=xm?>.<?=sym(j,l)?>
			)<?
			end
		end
	end
?>,
<? end
?>	};


	sym3 DBar2_phi_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = partial2_phi_ll[<?=ij-1?>]<?
	for k,xk in ipairs(xNames) do
		?> - connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l[<?=k-1?>]<?
	end
?>,
<? end
?>	};

	real tr_DBar2_phi = sym3_dot(gammaBar_uu, DBar2_phi_ll);

	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	sym3 RPhi_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = 
			- 2. * DBar2_phi_ll.<?=xij?>
			+ 4. * partial_phi_l[<?=i-1?>] * partial_phi_l[<?=j-1?>]
			+ gammaBar_ll.<?=xij?> * (
				- 2. * tr_DBar2_phi
				- 4. * real3_weightedLenSq(
					*(real3*)partial_phi_l,
					gammaBar_uu
				)
			),
<? end
?>	};

	sym3 R_ll = sym3_add(RPhi_ll, RBar_ll);


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
	/*
	DBar_j (e^(6 phi) ABar^ij)
	= DBar_j (e^(6 phi)) ABar^ij + e^(6 phi) DBar_j ABar^ij
	= 6 phi_,j e^(6 phi) ABar^ij + e^(6 phi) (ABar^ij_,j + connBar^i_kj ABar^kj + connBar^j_kj ABar^ik)
	= exp(6 phi) (6 ABar^ij phi_,j + (gammaBar^ik ABar_kl gammaBar^lj)_,j + connBar^i_jk ABar^jk) ... plus (ln det gammaBar_ll)_,k which is zero, right?
	= exp(6 phi) (
			+ 6 ABar^ij phi_,j
			+ gammaBar^ik_,j ABar_k^j
			+ ABar^i_l gammaBar^lj_,j
			+ gammaBar^ik ABar_kl,j gammaBar^lj
			+ connBar^i_jk ABar^jk)
   	= exp(6 phi) (
			+ 6 ABar^ij phi_,j
			- ABar^kj gammaBar^il gammaBar_lk,j
			- ABar^ik gammaBar^mj gammaBar_km,j
			+ gammaBar^ik gammaBar^lj ABar_kl,j
			+ connBar^i_jk ABar^jk)

	B&S 11.49
	0 = M^i = DBar_j (e^(6 phi) ABar^ij) - 2/3 e^(6 phi) DBar^i K - 8 pi e^(6 phi) S^i
	M^i = exp(6 phi) (
				+ 6 ABar^ij phi_,j
				+ connBar^i_jk ABar^jk
				- ABar^kj gammaBar^li gammaBar_kl,j
				- ABar^ki gammaBar^lj gammaBar_kl,j
				+ gammaBar^ik gammaBar^lj ABar_kl,j
				- 2/3 K_,j gammaBar^ij 
				- 8 pi S^i
			)
	*/
	real exp_6phi = 1. / (U->W * U->W * U->W);
<? for i,xi in ipairs(xNames) do
?>	U->M_u.<?=xi?> = exp_6phi * (
		- 8. * M_PI * U->S_u.<?=xi?>
<?	for j,xj in ipairs(xNames) do
?>		+ 6. * ABar_uu.<?=sym(i,j)?> * partial_phi_l[<?=j-1?>]
		- 2./3. * exp_6phi * gammaBar_uu.<?=sym(i,j)?> * partial_K_l[<?=j-1?>]
<?		for k,xk in ipairs(xNames) do
?>		- connBar_ull.<?=xi?>.<?=sym(j,k)?> * ABar_uu.<?=sym(j,k)?>
<?			for l,xl in ipairs(xNames) do
?>		- ABar_uu.<?=sym(k,j)?> * gammaBar_uu.<?=sym(l,i)?> * partial_epsilon_lll[<?=j-1?>].<?=sym(k,l)?>
		- ABar_uu.<?=sym(k,i)?> * gammaBar_uu.<?=sym(l,j)?> * partial_epsilon_lll[<?=j-1?>].<?=sym(k,l)?>
		+ gammaBar_uu.<?=sym(k,i)?> * gammaBar_uu.<?=sym(l,j)?> * partial_ABar_lll[<?=j-1?>].<?=sym(k,l)?>
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
