<? 
local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym

local makePartial = function(...) return eqn:makePartial(...) end
local makePartial2 = function(...) return eqn:makePartial2(...) end


-- integrates whatsoever.
local useCalcDeriv = true
local useCalcDeriv_alpha = true
local useCalcDeriv_W = true
local useCalcDeriv_epsilon_LL = true
local useCalcDeriv_K = true
local useCalcDeriv_ABar_LL = true
local useCalcDeriv_LambdaBar_U = true
local useCalcDeriv_beta_U = true

-- constrains det gammaBar_ij = det gammaHat_ij, ABar^i_i = 0, and calculates H and M^i ... if the associated flags are set
local useConstrainU = true

-- does Kreiss-Oligar dissipation
local useAddSource = false
?>

/*
I through I would save on loc by replacing expanded summation with vector/matrix operators, but it just made things worse.
The operators are pass-by-value (because I have read for C/C++ at least that this is faster for float3's,
 and for GPU's it will probably be faster for anything that fits in a float4/float8 or whatever the architecture is designed for).
However I've also read that Intel OpenCL used to have bugs with pass-by-value.
Then there's the possibility that, by moving the calc_RBar_ll code into its own function, maybe it uses too many arguments?
Then there's the possibility that I'm using too many locals and the compiler has a bug with stack allocations beyond a certain size.
*/

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

notice that gamma'_ij -> f gamma_ij; gamma'^ij -> 1/f gamma'^ij will produce the same result
so feel free to use gammaBar, gammaHat, etc
*/
sym3 tracefree(sym3 A_ll, sym3 g_ll, sym3 g_uu) {
	real tr_A = sym3_dot(A_ll, g_uu);
	return sym3_sub(A_ll, sym3_real_mul(g_ll, tr_A / 3.));
}

//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
kernel void calcDeriv(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
<? if useCalcDeriv then ?>
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global cons_t* deriv = derivBuf + index;
	const global cons_t* U = UBuf + index;

//NOTICE I'm scaling these back into coordinate form before calculating the partial derivative, while the SENR code applies the chain rule and adds the partial of the non-coordatein form times the scale transform plus the partial of the scale transform times the non-coordinate form.
//will that affect my accuracy?
<?=eqn:makePartial'alpha'?>			//partial_alpha_l[i] := alpha_,i
<?=eqn:makePartial'beta_U'?>		//partial_beta_ul[j].i := beta^i_,j
<?=eqn:makePartial'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial'W'?>				//partial_W_l[i] := W_,i 
<?=eqn:makePartial'K'	?>			//partial_K_l[i] := K,i
<?=eqn:makePartial'ABar_LL'?>		//partial_ABar_lll[k].ij = ABar_ij,k
<?=eqn:makePartial'LambdaBar_U'?>	//partial_LambdaBar_ul[j].i := connBar^i_,j
<?=eqn:makePartial2'alpha'?>		//partial2_alpha_ll.ij := alpha_,ij
<?=eqn:makePartial2'beta_U'?>		//partial2_beta_ull[jk].i = beta^i_,jk
<?=eqn:makePartial2'epsilon_LL'?>	//partial2_epsilon_llll[kl].ij = epsilon_ij,kl = gammaBar_ij,kl
<?=eqn:makePartial2'W'?>			//partial2_W_ll.ij := W_,ij
	//partial_B[i] := B^i_,t
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=eqn:makePartial'B_U'?>
<? end ?>


	//gammaHat_ll.ij := gammaHat_ij
	sym3 gammaHat_ll = coord_g_ll(x);
	
	//gammaHat_uu.ij := gammaHat^ij
	sym3 gammaHat_uu = coord_g_uu(x);
	
	//connHat_ull.i.jk := connHat^i_jk
	_3sym3 connHat_ull = coord_conn_ull(x);

	//partial_connHat_ulll[l].i.jk = connHat^i_jk,l
	_3sym3 partial_connHat_ulll[3];
	coord_partial_conn_ulll(partial_connHat_ulll, x);

	//epsilon_ll.ij := epsilon_ij = e_i^I e_j^I epsilon_IJ
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);

	//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);

	//partial_gammaHat_lll.k.ij := gammaHat_ij,k
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for ij,xij in ipairs(symNames) do
	for k,xk in ipairs(xNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

	real det_gammaHat = calc_det_gammaHat(x);
	real3 partial_det_gammaHat_l = coord_partial_det_g(x);
	sym3 partial2_det_gammaHat_ll = coord_partial2_det_g(x);

	/*
	Etienne's SENR Mathematica notebook has '*  detg'...
	 ... but the paper says det gammaBar_ij = det gammaHat_ij
	 ... so what is this mysterious 'detg' factor?
	 TODO find where detg comes from 
	
	From 2017 Ruchlin: 
	
	"We choose ˆγ ≡ det(ˆγij ) = ¯γ in the initial data for
	all applications in this paper. Since both determinants
	remain independent of time, they remain equal to each
	other throughout the evolution."
	
	... so why is there a distinction between det gammaBar_ij and det gammaHat_ij? 
	TODO detg ...
	*/
	real detg = 1.;
	real3 partial_detg_l = real3_zero;
	sym3 partial2_detg_ll = sym3_zero;
	real det_gammaBar = det_gammaHat * detg;
	real3 partial_det_gammaBar_l = real3_add(
		real3_real_mul(partial_det_gammaHat_l, detg),
		real3_real_mul(partial_detg_l, det_gammaHat));
	sym3 partial2_det_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	partial2_det_gammaBar_ll.<?=xij?> = 
		partial2_det_gammaHat_ll.<?=xij?> * detg
		+ partial_det_gammaHat_l.<?=xi?> * partial_detg_l.<?=xj?>
		+ partial_det_gammaHat_l.<?=xj?> * partial_detg_l.<?=xi?>
		+ det_gammaHat * partial2_detg_ll.<?=xij?>;
<? end
?>

	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

	real3 beta_u = real3_rescaleToCoord_U(U->beta_U, x);

	//////////////////////////////// alpha_,t //////////////////////////////// 
<? if useCalcDeriv_alpha then ?>

	//Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	real Q = calc_f(U->alpha) * U->K;

	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 Q + alpha,i beta^i
	deriv->alpha += -U->alpha * U->alpha * Q
<? for i,xi in ipairs(xNames) do
?>		+ partial_alpha_l[<?=i-1?>] * beta_u.<?=xi?>
<? end
?>	;

<? end	-- useCalcDeriv_alpha ?>
	//////////////////////////////// W_,t //////////////////////////////// 
<? if useCalcDeriv_W then ?>

	sym3 gammaBar_LL = sym3_rescaleFromCoord_ll(gammaBar_ll, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);

	_3sym3 connHat_lll = coord_conn_lll(x);

<?=eqn:getCode_connBar_ull()?>

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	//I'm reversing the index on this
	real3x3 partial_beta_UL;
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	partial_beta_UL.<?=xi?>.<?=xj?> = partial_beta_ul[<?=j-1?>].<?=xi?>
			/ coord_dx<?=j-1?>(x) * coord_dx<?=i-1?>(x);
<?	end
end ?>

	real3 tr12_connBar_l = _3sym3_tr12(connBar_ull);

	//connBar^j_kj beta^k
	real tr12_connBar_dot_beta = real3_dot(tr12_connBar_l, beta_u);

	//2017 Ruchlin et al eqn 11c
	//W,t = 1/3 W (alpha K - beta^k connBar^j_kj - beta^k_,k) + beta^k W_,k
	deriv->W += (1. / 3.) * U->W * (
			U->alpha * U->K 
			- tr12_connBar_dot_beta 
			- tr_partial_beta
		)
<? for i,xi in ipairs(xNames) do
?>		+ beta_u.<?=xi?> * partial_W_l[<?=i-1?>]
<? end
?>	;

<? end	-- useCalcDeriv_W ?>
	//////////////////////////////// epsilon_ij,t //////////////////////////////// 
	
	/*
	tr_DBar_beta := DBar_j beta^j = beta^j_,j + connBar^j_kj beta^k
	Etienne's SENR Mathematica notebook uses this instead:
	tr_DBar_beta := beta^j_,j + beta^j gammaBar_,j / (2 gammaBar)
	*/
	//real tr_DBar_beta = tr_partial_beta + tr12_connBar_dot_beta;
	real tr_DBar_beta = tr_partial_beta + real3_dot(beta_u, partial_det_gammaBar_l);

<? if useCalcDeriv_epsilon_LL then ?>

	/*
	2017 Ruchlin et al, eqn 11a
	epsilon_ij,t = 2/3 gammaBar_ij (alpha ABar^k_k - DBar_k beta^k) + DHat_i beta_j + DHat_j beta_i - 2 alpha ABar_ij + epsilon_ij,k beta^k + epsilon_ik beta^k_,j + epsilon_kj beta^k_,i
	...using DBar_(i beta_j) = DHat_(i beta_j) + epsilon_k(j beta^k_,i) + 1/2 epsilon_ij,k beta^k
	= 	
		+ 2/3 gammaBar_ij (
			alpha ABar^k_k 
			- beta^k_,k 
			- connBar^k_lk beta^l
		) 
		- 2 alpha ABar_ij 
		+ .5 DBar_i beta_j
		+ .5 DBar_j beta_i
	

	Etienne's SENR Mathematica notebook:
	= 
		// Lie derivative terms
		beta^k gammaBar_ij,k
		+ gammaBar_ki beta^k_,j
		+ gammaBar_kj beta^k_,i

		- 2/3 gammaBar_ij DBar_k beta^k
		(notice that "2/3 alpha ABar^k_k" is excluded)
		- 2 alpha ABar_ij

	Looks like the paper's notation DBar_i beta_j implies DBar_j (gammaBar_ik beta^k) = gammaBar_ki DBar_j beta^k
	...which is not DBar_j ( gamma_ik beta^k ), which was my interpretation of the definition of beta_i := gamma_ij beta^k
	*/
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	deriv->epsilon_LL.<?=xij?> += 
		- 2. / 3. * gammaBar_LL.<?=xij?> * tr_DBar_beta
		- 2 * U->alpha * U->ABar_LL.<?=xij?>
		+ (0.
<? 	for k,xk in ipairs(xNames) do
?>			+ partial_gammaBar_lll.<?=xk?>.<?=xij?> * beta_u.<?=xk?>
			+ gammaBar_ll.<?=sym(k,i)?> * partial_beta_ul[<?=j-1?>].<?=xk?>
			+ gammaBar_ll.<?=sym(k,j)?> * partial_beta_ul[<?=i-1?>].<?=xk?>
<? 	end
?>		) / coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x);
<? end
?>

<? end	-- useCalcDeriv_epsilon_LL ?>
	//////////////////////////////// K_,t //////////////////////////////// 

	//ABar^i_j = gammaBar^ik ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);

<? if useCalcDeriv_K or useCalcDeriv_ABar_LL then ?>
	
	//exp(-4 phi)
	real exp_neg4phi = calc_exp_neg4phi(U);

	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}


	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]
<?	for k,xk in ipairs(xNames) do
?>		- partial_alpha_l[<?=k-1?>] * connBar_ull.<?=xk?>.<?=xij?>
<?	end
?>	;
<? end
?>

<? end -- useCalcDeriv_K or useCalcDeriv_ABar_LL ?>
<? if useCalcDeriv_K then ?>

	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	//ABar^IJ = ABar^I_K gammaBar^KJ

	//tr_DBar2_alpha := gammaBar^ij DBar_i DBar_j alpha
	real tr_DBar2_alpha = sym3_dot(gammaBar_uu, DBar2_alpha_ll);
	
	//tr_ABarSq := ABar_ij ABar^ij = ABar_ij ABar_kl gammaBar^ik gammaBar^jl
	real tr_ABarSq = sym3_dot(U->ABar_LL, ABar_UU);
	
	/*
	gammaBar_ij = exp(-4 phi) gamma_ij
	gammaBar^ij = exp(4 phi) gamma^ij
	gamma^ij = exp(-4 phi) gammaBar^ij
	S := S_ij gamma^ij = exp(-4 phi) S_ij gammaBar^ij 
	*/
	real S = exp_neg4phi * sym3_dot(U->S_ll, gammaBar_uu);
	
	/*
	B&S 11.52
	Alcubierre 2.8.12
	K_,t = -gamma^ij D_i D_j alpha + alpha (ABar_ij ABar^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	2017 Ruchlin et al
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
	deriv->K += 0.
		+ U->alpha * U->K * U->K / 3.
		+ U->alpha * tr_ABarSq
		- exp_neg4phi * (0.
			+ tr_DBar2_alpha
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>			+ 2. * partial_phi_l.<?=xi?> * partial_alpha_l[<?=j-1?>] * gammaBar_uu.<?=sym(i,j)?>
<?	end
end
?>		)
<? for i,xi in ipairs(xNames) do
?>		+ beta_u.<?=xi?> * partial_K_l[<?=i-1?>]
<? end
?>		+ 4. * M_PI * U->alpha * (U->rho + S)
	;

<? end	-- useCalcDeriv_K ?>
	//////////////////////////////// ABar_ij_,t //////////////////////////////// 
<? if useCalcDeriv_ABar_LL or useCalcDeriv_LambdaBar_U then ?>

	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	
	//Delta_ull[i].jk := Delta^i_jk = connBar^i_jk + connHat^i_jk
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);

	
<? end -- useCalcDeriv_ABar_LL or useCalcDeriv_LambdaBar_U ?>
<? if useCalcDeriv_ABar_LL then ?>
	
	//Delta_u.i := Delta^i = LambdaBar^i - C^i
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);
	real3 Delta_u = real3_rescaleToCoord_U(Delta_U, x);

	real3 LambdaBar_u = real3_rescaleToCoord_U(U->LambdaBar_U, x);

	sym3 RBar_ll;
	{
		_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);
		<?=eqn:getCode_RBar_ll()?> 
	}
	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?	end
?>	;
<? end
?>
	
	/*
	2017 Ruchlin et al eqn 11b
	traceless portion of ...
		-2 alpha DBar_i DBar_j phi 
		+ 4 alpha DBar_i phi DBar_j phi 
		+ 2 DBar_i phi DBar_j alpha 
		+ 2 DBar_i alpha DBar_j phi 
		- DBar_i DBar_j alpha 
		+ alpha RBar_ij 
		- 8 pi alpha S_ij
	
	SENR does trace-free on a few terms individually:
	*) -DBar_i DBar_j alpha
	*) alpha RBar_ij
	*) everything else

	*/
	sym3 TF_DBar2_alpha_ll = tracefree(DBar2_alpha_ll, gammaBar_ll, gammaBar_uu);

	sym3 TF_RBar_ll = tracefree(RBar_ll, gammaBar_ll, gammaBar_uu);
	
	sym3 tracelessPart_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	tracelessPart_ll.<?=xij?> = 0.
			+ 2. * partial_phi_l.<?=xi?> * partial_alpha_l[<?=j-1?>]
			+ 2. * partial_phi_l.<?=xj?> * partial_alpha_l[<?=i-1?>]
			+ U->alpha * (0.
				- 2. * DBar2_phi_ll.<?=xij?>
				+ 4. * partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?>
				- 8. * M_PI * U->S_ll.<?=xij?>
			)
		;
<? end
?>
	tracelessPart_ll = tracefree(tracelessPart_ll, gammaBar_ll, gammaBar_uu);

	/*
	2017 Ruchlin et al, eqn. 11b
	ABar_ij,t = 
		- 2/3 ABar_ij DBar_k beta^k
		- 2 alpha ABar_ik ABar^k_j
		+ alpha ABar_ij K
		+ exp(-4 phi) (trace-free part above)_ij
		+ ABar_ij,k beta^k
		+ beta^k_,i ABar_jk
		+ beta^k_,j ABar_ik
	*/
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	deriv->ABar_LL.<?=xij?> += 0.
		- 2. / 3. * U->ABar_LL.<?=xij?> * tr_DBar_beta
<?	for k,xk in ipairs(xNames) do
?>		- 2. * U->alpha * U->ABar_LL.<?=sym(i,k)?> * ABar_UL.<?=xk?>.<?=xj?>
<?	end
?>		+ U->alpha * U->K * U->ABar_LL.<?=xij?>
		+ exp_neg4phi * (0.
			+ tracelessPart_ll.<?=xij?>
			- TF_DBar2_alpha_ll.<?=xij?>
			+ U->alpha * TF_RBar_ll.<?=xij?>//diverging near r=0
		) / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x))
<?	for k,xk in ipairs(xNames) do
?>		+ partial_ABar_lll[<?=k-1?>].<?=xij?> * beta_u.<?=xk?>
			/ (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x))
		+ U->ABar_LL.<?=sym(j,k)?> * partial_beta_UL.<?=xk?>.<?=xi?>
		+ U->ABar_LL.<?=sym(i,k)?> * partial_beta_UL.<?=xk?>.<?=xj?>
<? 	end
?>	;
<? end
?>

<? end	-- useCalcDeriv_ABar_LL ?>
	//////////////////////////////// LambdaBar^i_,t //////////////////////////////// 
<? if useCalcDeriv_LambdaBar_U then ?>
	
	/*
	DHat2_beta_u.i.j.k = DHat_k DHat_j beta^i
	= DHat_k (beta^i_,j + connHat^i_lj beta^l)
	= (beta^i_,j + connHat^i_lj beta^l)_,k
		+ connHat^i_mk (beta^m_,j + connHat^m_lj beta^l)
		- connHat^m_jk (beta^i_,m + connHat^i_lm beta^l)
	= beta^i_,jk 
		+ connHat^i_lj,k beta^l
		+ connHat^i_lj beta^l_,k
		+ connHat^i_lk beta^l_,j
		- connHat^l_jk beta^i_,l
		+ connHat^i_mk connHat^m_lj beta^l
		- connHat^m_jk connHat^i_lm beta^l
	(This is breaking convention.  usually derivatives are stored outer-most.)
	*/
	real3x3x3 DHat2_beta_u;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
			local jk = from3x3to6(j,k)
?>	DHat2_beta_u.<?=xi?>.<?=xj?>.<?=xk?> = 
		partial2_beta_ull[<?=jk-1?>].<?=xi?>
<?		for l,xl in ipairs(xNames) do
?>		+ partial_connHat_ulll[<?=k-1?>].<?=xi?>.<?=sym(l,j)?> * beta_u.<?=xl?>
		+ connHat_ull.<?=xi?>.<?=sym(l,j)?> * partial_beta_ul[<?=k-1?>].<?=xl?>
		+ connHat_ull.<?=xi?>.<?=sym(l,k)?> * partial_beta_ul[<?=j-1?>].<?=xl?>
		- connHat_ull.<?=xl?>.<?=sym(j,k)?> * partial_beta_ul[<?=l-1?>].<?=xi?>
<?			for m,xm in ipairs(xNames) do
?>		+ connHat_ull.<?=xi?>.<?=sym(m,k)?> * connHat_ull.<?=xm?>.<?=sym(l,j)?> * beta_u.<?=xl?>
		- connHat_ull.<?=xi?>.<?=sym(m,l)?> * connHat_ull.<?=xm?>.<?=sym(k,j)?> * beta_u.<?=xl?>
<?			end
		end
?>	;
<?		end
	end 
end
?>

	/*
	tr12_partial2_beta_l.i := beta^j_,ji
	beta^i_,jk is a sym3 of 3's ... so I don't have that struct yet ... 
	 What name could I use? sym3x3?  how about real3s3x3 where we have 's' for symmetric and 'x' for cartesian product.
	 Then sym3 turns into real3s3 and _3sym3 turns into real3x3s3.
	*/
	real3 tr12_partial2_beta_l = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		local ij = from3x3to6(i,j)
		local xij = symNames[ij]
?> + partial2_beta_ull[<?=ij-1?>].<?=xj?><?
	 end ?>,
<? end
?>	};

	/*
	DBar_tr_DBar_beta_l.i = DBar_i DBar_j beta^j
	= DBar_i (beta^j_,j + connBar^j_lj beta^l)
	= (beta^j_,j + connBar^j_lj beta^l)_,i
		+ connBar^j_mi (beta^m_,j + connBar^m_lj beta^l)
		- connBar^m_ji (beta^j_,m + connBar^j_lm beta^l)
	= beta^j_,ji 
		+ connBar^j_lj,i beta^l
		+ connBar^j_lj beta^l_,i
	
	Etienne's SENR uses this alternative formulation: 
	= beta^j_,ji
		+ 1/2 (
			gammaBar_,ij beta^j
			- gammaBar_,i gammaBar_,j beta^j / gammaBar
			+ gammaBar_,j beta^j_,i
		) / gammaBar
	*/
	real3 DBar_tr_DBar_beta_l;
<? for i,xi in ipairs(xNames) do
?>	DBar_tr_DBar_beta_l.<?=xi?> = 
		tr12_partial2_beta_l.<?=xi?>
<? 	for j,xj in ipairs(xNames) do
		local ij = from3x3to6(i,j)
		local xij = symNames[ij]
?>		+ .5 / det_gammaBar * (
			partial2_det_gammaBar_ll.<?=xij?> * beta_u.<?=xj?>
			- partial_det_gammaBar_l.<?=xi?> * partial_det_gammaBar_l.<?=xj?> * beta_u.<?=xj?> / det_gammaBar
			+ partial_det_gammaBar_l.<?=xj?> * partial_beta_ul[<?=i-1?>].<?=xj?>
		)
<?	end
?>	;
<? end
?>
	//DBar_tr_DBar_beta_u.i = DBar^i DBar_k beta^k = gammaBar^ij DBar_j DBar_k beta^k
	real3 DBar_tr_DBar_beta_u = sym3_real3_mul(gammaBar_uu, DBar_tr_DBar_beta_l);

	//tr_gammaBar_DHat2_beta_u.i = gammaBar^jk DHat_j DHat_k beta^i
	real3 tr_gammaBar_DHat2_beta_u = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_real3x3_dot(gammaBar_uu, DHat2_beta_u.<?=xi?>),
<? end
?>	};
	
	sym3 ABar_uu = sym3_rescaleToCoord_UU(ABar_UU, x);
	
	/*
	LambdaBar^i_,t = 
		gammaBar^jk DHat_j DHat_k beta^i
		+ 2/3 Delta^i DBar_j beta^j
		+ 1/3 DBar^i DBar_j beta^j
		- 2 ABar^ij (alpha_,j - 6 phi_,j)
		+ 2 alpha ABar^jk Delta^i_jk		-- NOTICE the paper doesn't have the extra alpha, but the Etienne SENR Mathematica notebook does have it.  Maybe I'm reading an earlier version of the paper?
		- 4/3 alpha gammaBar^ij K_,j
		- 16 pi exp(4 phi) alpha S^i
		+ LambdaBar^i_,k beta^k
		- LambdaBar^k beta^i_,k
	*/
	real3 dt_LambdaBar_U;
<? for i,xi in ipairs(xNames) do
?>	dt_LambdaBar_U.<?=xi?> = (0.
		+ tr_gammaBar_DHat2_beta_u.<?=xi?>
		+ 2. / 3. * Delta_u.<?=xi?> * tr_DBar_beta
		+ 1. / 3. * DBar_tr_DBar_beta_u.<?=xi?>
<?	for j,xj in ipairs(xNames) do
		local xij = sym(i,j)
		local jj = from3x3to6(j,j)
?>		- 2. * ABar_uu.<?=xij?> * (0.
			+ partial_alpha_l[<?=j-1?>]
			- 6. * partial_phi_l.<?=xj?>
		)
<?		for k,xk in ipairs(xNames) do		
			local xik = sym(i,k)
			local jk = from3x3to6(j,k)
			local xjk = symNames[jk]
?>		+ 2. * U->alpha * Delta_ull.<?=xi?>.<?=xjk?> * ABar_uu.<?=xjk?>
<?		end
?>		- 4. / 3. * U->alpha * gammaBar_uu.<?=xij?> * partial_K_l[<?=j-1?>] 
		+ beta_u.<?=xj?> * partial_LambdaBar_ul[<?=j-1?>].<?=xi?>
		- LambdaBar_u.<?=xj?> * partial_beta_ul[<?=j-1?>].<?=xi?>
<?	end
?>		- 16. * M_PI * U->alpha * U->S_u.<?=xi?> / exp_neg4phi
	) * coord_dx<?=i-1?>(x);
<? end
?>

	deriv->LambdaBar_U = real3_add(deriv->LambdaBar_U, dt_LambdaBar_U);

<? end	-- useCalcDeriv_LambdaBar_U ?>
	//////////////////////////////// beta^i_,t and B^i_,t //////////////////////////////// 
<? if useCalcDeriv_beta_U then ?>

<? if eqn.useShift == 'GammaDriver' then ?>
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k (connBar^i_,t + eta connBar^i)
	const real k = 3. / 4.;
	const real eta = 1.;	//1.;	// 1 / (2 M), for total mass M
<? for i,xi in ipairs(xNames) do
?>	deriv->beta_U.<?=xi?> += k * dt_LambdaBar_U.<?=xi?> + eta * U->LambdaBar_U.<?=xi?>;
<? end
?>
<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>
	/*
	hyperbolic Gamma driver 
	2017 Ruchlin et al, eqn 14a, 14b
	beta^i_,t = B^i + beta^i_,j beta^j - beta^i_,j beta^j = B^i
	B^i_,t = 3/4 (
			LambdaBar^i_,t 
			- LambdaBar^i_,j beta^j 
			+ LambdaBar^j beta^i_,j
		) - eta B^i 
		+ B^i_,j beta^j 
		- B^j beta^i_,j
	*/
<? for i,xi in ipairs(xNames) do
?>	deriv->beta_U.<?=xi?> += U->B_U.<?=xi?>;
<? end
?>
	const real eta = 1.;
<? for i,xi in ipairs(xNames) do
?>	deriv->B_U.<?=xi?> += .75 * (
			dt_LambdaBar_U.<?=xi?>
<? 	for j,xj in ipairs(xNames) do
?>			- partial_LambdaBar_ul[<?=j-1?>].<?=xi?> * beta_u.<?=xj?> * coord_dx<?=i-1?>(x)
			+ partial_beta_UL.<?=xi?>.<?=xj?> * U->LambdaBar_U.<?=xj?>
<?	end
?>		)
		- eta * U->B_U.<?=xi?>
<?	for j,xj in ipairs(xNames) do
?>		+ partial_B_ul[<?=j-1?>].<?=xi?> * beta_u.<?=xj?> * coord_dx<?=i-1?>(x)
		- partial_beta_UL.<?=xi?>.<?=xj?> * U->B_U.<?=xj?>
<?	end
?>	;
<? end ?>

<? end	-- eqn.useShift ?>
<? end	-- useCalcDeriv_beta_U ?>
<? end 	-- useCalcDeriv ?>
}

kernel void constrainU(
	constant solver_t* solver,
	global cons_t* UBuf
) {
<? if useConstrainU then ?>	
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global cons_t* U = UBuf + index;

	sym3 gammaHat_ll = coord_g_ll(x);
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);
<? 
if eqn.guiVars.constrain_det_gammaBar.value 
or eqn.guiVars.constrain_tr_ABar.value 
then 
?>
	/*
	we need to force det(gammaBar_ij) = det(gammaHat_ij)
	and we can do so with gammaBar_ij := gammaBar_ij * (det(gammaHat_ij) / det(gammaBar_ij))^(1/3)
	so we find
	det(gammaBar_ij * (det(gammaHat_ij) / det(gammaBar_ij))^(1/3) )
	= det(gammaBar_ij) * (det(gammaHat_ij) / det(gammaBar_ij))
	= det(gammaHat_ij)
	*/
<?	if eqn.guiVars.constrain_det_gammaBar.value then ?>
	real det_gammaHat = calc_det_gammaHat(x);
	real det_gammaBar = sym3_det(gammaBar_ll);
	real rescaleMetric = cbrt(det_gammaHat/det_gammaBar);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_ll.<?=xij?> *= rescaleMetric;
<? 		end ?>
	epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);
<?	end	-- constrain_det_gammaBar ?>

	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	
	sym3 gammaBar_LL = sym3_rescaleFromCoord_ll(gammaBar_ll, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
	
	//in Buchman's paper it says he doesn't do this
	//and in the new arbitrary-coord formalism, there is a tr ABar_ij term
<? if eqn.guiVars.constrain_tr_ABar.value then ?>
	U->ABar_LL = tracefree(U->ABar_LL, gammaBar_LL, gammaBar_UU);
<? end	-- constrain_tr_ABar ?>

<? else -- constrain_det_gammaBar or constrain_tr_ABar ?>
	real det_gammaHat = calc_det_gammaHat(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaHat);
	
	sym3 gammaBar_LL = sym3_rescaleFromCoord_ll(gammaBar_ll, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
<? end -- constrain_det_gammaBar or constrain_tr_ABar ?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

//TODO these need to be pre-scaled back to coordinates before computing the weighted finite difference
<?=eqn:makePartial'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial'ABar_LL'?>		//partial_ABar_lll[k].ij = ABar_ij,k
<?=eqn:makePartial'K'	?>			//partial_K_l[i] := K,i
<?=eqn:makePartial'W'?>				//partial_W_l[i] := phi_,i 
<?=eqn:makePartial'LambdaBar_U'?>	//partial_LambdaBar_ul[j].i := connBar^i_,j
<?=eqn:makePartial2'W'?>			//partial2_W_ll.ij := phi_,ij
<?=eqn:makePartial2'epsilon_LL'?>

	real3 partial_phi_l;
<? for i,xi in ipairs(xNames) do
?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

	sym3 partial2_phi_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_phi_ll.<?=xij?> = .5 * (
			-partial2_W_ll[<?=ij-1?>] 
			+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
		) / U->W;
<? end ?>
	
	real exp_neg4phi = calc_exp_neg4phi(U);
	
	_3sym3 connHat_lll = coord_conn_lll(x);
	_3sym3 connHat_ull = coord_conn_ull(x);

	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>
<?=eqn:getCode_connBar_ull()?>

	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar_uu.ij := ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	
	
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);
	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);
	
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);
	real3 Delta_u = real3_rescaleToCoord_U(Delta_U, x);
	
	real3 LambdaBar_u = real3_rescaleToCoord_U(U->LambdaBar_U, x);

	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	sym3 gammaHat_uu = coord_g_uu(x);
	
	_3sym3 partial_connHat_ulll[3];
	coord_partial_conn_ulll(partial_connHat_ulll, x);

	sym3 RBar_ll;
	{
		<?=eqn:getCode_RBar_ll()?> 
	}	
	//RBar := RBar_ij gammaBar^ij
	real RBar = sym3_dot(gammaBar_uu, RBar_ll);

	//the old way actually cached DBar2_phi_ll for calculating RPhi_ll, but that wasn't even being used, so I deleted it
	//tr_DBar2_phi := gammaBar^ij DBar_i DBar_j phi = gammaBar^ij phi_,ij - connBar^k phi_,k
	real tr_DBar2_phi = 0.
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		local ij = from3x3to6(i,j)
		local xij = symNames[ij]
?>		+ gammaBar_uu.<?=sym(i,j)?> * (
			partial2_phi_ll.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?		end
?>		)
<?	end
end
?>	;

	//2017 Ruchlin et al, eqn 46
	//H = 2/3 K^2 - ABar^ij ABar_ij + exp(-4 phi) (RBar - 8 DBar^i phi DBar_i phi - 8 gammaBar^ij DBar_i DBar_j phi)
	U->H = 2. / 3. * U->K * U->K
		- sym3_dot(U->ABar_LL, ABar_UU)
		+ exp_neg4phi * (
			RBar
			- 8. * real3_weightedLenSq(partial_phi_l, gammaBar_uu) 
			- 8. * tr_DBar2_phi
		)
		- 16. * M_PI * U->rho;

#if 1
	sym3 ABar_uu = sym3_rescaleToCoord_UU(ABar_UU, x);
	
	/*
	DBar_j (e^(6 phi) ABar^ij)
	= DBar_j (e^(6 phi)) ABar^ij + e^(6 phi) DBar_j ABar^ij
	= 6 phi_,j e^(6 phi) ABar^ij + e^(6 phi) (ABar^ij_,j + connBar^i_kj ABar^kj + connBar^j_kj ABar^ik)
	= exp(6 phi) (6 ABar^ij phi_,j + (gammaBar^ik ABar_kl gammaBar^lj)_,j + connBar^i_jk ABar^jk) ... plus (ln det gammaBar)_,k which is zero, right?
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
				- ABar^jk gammaBar^li gammaBar_kl,j
				- ABar^ik gammaBar^lj gammaBar_kl,j
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
?>		+ 6. * ABar_uu.<?=sym(i,j)?> * partial_phi_l.<?=xj?>
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
#else
	//2017 Ruchlin et al, eqn 47
	//M^i = exp(-4 phi) (DHat_j ABar^ij + 2 ABar^k(i Delta^j)_jk + 6 ABar^ij phi_,j - 2/3 gammaBar^ij K_,j)
#endif
<? end	-- calc_H_and_M ?>
<? end	-- useConstrainU ?>
}

kernel void addSource(
	constant solver_t* solver,
	global cons_t* derivBuf,
	const global cons_t* UBuf
) {
<? if useAddSource then ?>
	SETBOUNDS_NOGHOST();
	
	global cons_t* deriv = derivBuf + index;

	//Kreiss-Oligar dissipation
	//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
	//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
	for (int i = 0; i < numIntStates; ++i) {
<?=eqn:makePartial2('ptr[i]', 'real', 'partial2_Ui_ll', false)?>
		real lap = 0<?
for j,xj in ipairs(xNames) do
	local jj = from3x3to6(j,j)
?> + partial2_Ui_ll[<?=jj-1?>]<?
end
?>;
		deriv->ptr[i] -= solver->diffuseSigma/16. * lap;
	}
<? end -- addSource ?>
}

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = UBuf + index;

<? 
-- seems all the hyperbolic formalisms listed in Alcubierre's book use alpha sqrt(gamma^ii) for the speed-of-light wavespeed
-- however the 2017 Ruchlin paper says to use gamma_ij
local useGammaInvForDT = false 	
?>

<? if useGammaInvForDT then ?>	
	sym3 gamma_uu = calc_gamma_uu(U, x);
<? else ?>	
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
<? end ?>

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
<? if useGammaInvForDT then ?>
		//this is asserting alpha and W >0, which they should be
		real absLambdaMax = U->alpha * sqrt(gamma_uu.<?=sym(side+1,side+1)?>);
<? else ?>		
		//2017 Ruchlin, eqn. 53
		real absLambdaMax = U->alpha * sqrt(gammaBar_ll.<?=sym(side+1,side+1)?>);
<? end ?>

<? if false then -- hmm, do we base our CFL on delta in coordinate, or delta in Cartesian? ?>
		real dx = cell_dx<?=side?>(x); 
<? else ?>
		real dx = solver->grid_dx.s<?=side?>;
<? end ?>	
		dt = (real)min(dt, dx / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt;
}
