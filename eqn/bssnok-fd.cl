<? 
-- [[ do this every time you use the env
local env = eqn:getEnv()
local oldEnv = getfenv()
getmetatable(env).__index = oldEnv
setfenv(1, env)
--]]

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

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

//NOTICE I'm scaling these back into coordinate form before calculating the partial derivative, while the SENR code applies the chain rule and adds the partial of the non-coordatein form times the scale transform plus the partial of the scale transform times the non-coordinate form.
//will that affect my accuracy?
<?=eqn:makePartial'alpha'?>			//partial_alpha_l[i] := alpha_,i
<?=eqn:makePartial'beta_U'?>			//partial_beta_Ul[j].I := beta^I_,j
<?=eqn:makePartial'epsilon_LL'?>		//partial_epsilon_LLl[k].IJ := epsilon_IJ,k
<?=eqn:makePartial'W'?>				//partial_W_l[i] := W_,i 
<?=eqn:makePartial'K'?>				//partial_K_l[i] := K,i
<?=eqn:makePartial'ABar_LL'?>		//partial_ABar_LLl[k].IJ = ABar_IJ,k
<?=eqn:makePartial'LambdaBar_U'?>	//partial_LambdaBar_Ul[j].I := LambdaBar^I_,j
<?=eqn:makePartial2'alpha'?>			//partial2_alpha_ll.ij := alpha_,ij
<?=eqn:makePartial2'beta_U'?>		//partial2_beta_Ull[jk].I = beta^I_,jk
<?=eqn:makePartial2'W'?>				//partial2_W_ll.ij := W_,ij
	//partial_B[i] := B^i_,t
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=eqn:makePartial'B_U'?>			//partial_B_Ul[j].I := B^I_,j
<? end ?>

<?=assign_sym3'gammaHat_ll'?>
<?=assign'det_gammaHat'?>
<?=assign_sym3'gammaHat_uu'?>	
<?=assign_3sym3'connHat_ull'?>
<?=assign_real3'partial_det_gammaHat_l'?>
<?=assign_3sym3x3'partial_connHat_ulll'?>
<?=assign_sym3'partial2_det_gammaHat_ll'?>

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
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
	real3 partial_det_gammaBar_l = partial_det_gammaHat_l;
	sym3 partial2_det_gammaBar_ll = partial2_det_gammaHat_ll;
<?=assign_sym3'gammaBar_uu'?>


	//////////////////////////////// alpha_,t //////////////////////////////// 
<? if useCalcDeriv_alpha then ?>

<?=assign'dt_alpha'?>
	deriv->alpha += dt_alpha;

<? end	-- useCalcDeriv_alpha ?>
	//////////////////////////////// W_,t //////////////////////////////// 
<? if useCalcDeriv_W then ?>

<?=assign_real3'tr_connBar_l'?>
<?=assign'tr_DBar_beta'?>
<?=assign'dt_W'?>
	deriv->W += dt_W;

<? end	-- useCalcDeriv_W ?>
	//////////////////////////////// K_,t //////////////////////////////// 
	
	//exp(-4 phi)
	real exp_neg4phi = calc_exp_neg4phi(U);

<? if useCalcDeriv_K then ?>
	/*
	gammaBar_ij = exp(-4 phi) gamma_ij
	gammaBar^ij = exp(4 phi) gamma^ij
	gamma^ij = exp(-4 phi) gammaBar^ij
	S := S_ij gamma^ij = exp(-4 phi) S_ij gammaBar^ij 
	*/
	real S = exp_neg4phi * sym3_dot(U->S_ll, gammaBar_uu);

<?=assign_real3x3'ABar_ul'?>
<?=assign_real3x3'ABarSq_ul'?>
<?=assign'tr_ABarSq'?>
<?=assign_3sym3'connBar_ull'?>
<?=assign_sym3'DBar2_alpha_ll'?> 
<?=assign'trBar_DBar2_alpha'?>
<?=assign_real3'partial_phi_l'?>
<?=assign'dt_K'?>
	deriv->K += dt_K;

<? end	-- useCalcDeriv_K ?>
	//////////////////////////////// epsilon_ij,t //////////////////////////////// 
<? if useCalcDeriv_epsilon_LL then ?>

<?=assign_sym3'Lbeta_gammaBar_LL'?>
<?=assign_sym3'dt_epsilon_LL'?>
<? for ij,xij in ipairs(symNames) do
?>	deriv->epsilon_LL.<?=xij?> += dt_epsilon_LL.<?=xij?>;
<? end
?>

<? end	-- useCalcDeriv_epsilon_LL ?>
	//////////////////////////////// ABar_ij_,t //////////////////////////////// 
<? if useCalcDeriv_ABar_LL or useCalcDeriv_LambdaBar_U then ?>

<?=assign_sym3'partial2_phi_ll'?>
<?=assign_sym3sym3'partial2_gammaHat_llll'?> 
<?=assign_3sym3'Delta_ull'?>	
	
<? end -- useCalcDeriv_ABar_LL or useCalcDeriv_LambdaBar_U ?>
<? if useCalcDeriv_ABar_LL then ?>

<?=assign_real3'Delta_u'?>
<?=assign_real3'LambdaBar_u'?>
<?=assign_real3x3'partial_LambdaBar_ul'?>

<?=assign_sym3'gammaBar_ll'?>
	sym3 RBar_ll;
	{
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>
<?=assign_3sym3'Delta_lll'?>
		<?=eqn:getCode_RBar_ll()?> 
	}

<?=assign_sym3'DBar2_phi_ll'?>

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

<?=assign_3sym3('partial_ABar_lll', partial_ABar_lll:permute'_kij')?>

<?=assign_real3x3'partial_beta_ul'?>

	//I'm reversing the index on this
	real3x3 partial_beta_UL;
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	partial_beta_UL.<?=xi?>.<?=xj?> = partial_beta_ul.<?=xi?>.<?=xj?>
			/ coord_dx<?=j-1?>(x) * coord_dx<?=i-1?>(x);
<?	end
end ?>

	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);

	//ABar^i_j = gammaBar^ik ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = real3x3_trace(partial_beta_ul);

	real3 tr12_connBar_l = _3sym3_tr12(connBar_ull);

	real3 beta_u = real3_rescaleToCoord_U(U->beta_U, x);

	//connBar^j_kj beta^k
	real tr12_connBar_dot_beta = real3_dot(tr12_connBar_l, beta_u);

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
			+ U->alpha * TF_RBar_ll.<?=xij?>
		) / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x))
<?	for k,xk in ipairs(xNames) do
?>		+ partial_ABar_lll.<?=xk?>.<?=xij?> * beta_u.<?=xk?>
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

<?=assign_3sym3'partial2_beta_ull'?>

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
		partial2_beta_ull.<?=xi?>.<?=sym(j,k)?>
<?		for l,xl in ipairs(xNames) do
?>		+ partial_connHat_ulll[<?=k-1?>].<?=xi?>.<?=sym(l,j)?> * beta_u.<?=xl?>
		+ connHat_ull.<?=xi?>.<?=sym(l,j)?> * partial_beta_ul.<?=xl?>.<?=xk?>
		+ connHat_ull.<?=xi?>.<?=sym(l,k)?> * partial_beta_ul.<?=xl?>.<?=xj?>
		- connHat_ull.<?=xl?>.<?=sym(j,k)?> * partial_beta_ul.<?=xi?>.<?=xl?>
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
?> + partial2_beta_ull.<?=xj?>.<?=sym(j,i)?><?
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
			+ partial_det_gammaBar_l.<?=xj?> * partial_beta_ul.<?=xj?>.<?=xi?>
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
	
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	//ABar^IJ = ABar^I_K gammaBar^KJ
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
		+ beta_u.<?=xj?> * partial_LambdaBar_ul.<?=xi?>.<?=xj?>
		- LambdaBar_u.<?=xj?> * partial_beta_ul.<?=xi?>.<?=xj?>
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

<?=assign_real3x3'partial_B_ul'?>
	
	const real eta = 1.;
<? for i,xi in ipairs(xNames) do
?>	deriv->B_U.<?=xi?> += .75 * (
			dt_LambdaBar_U.<?=xi?>
<? 	for j,xj in ipairs(xNames) do
?>			- partial_LambdaBar_ul.<?=xi?>.<?=xj?> * beta_u.<?=xj?> * coord_dx<?=i-1?>(x)
			+ partial_beta_UL.<?=xi?>.<?=xj?> * U->LambdaBar_U.<?=xj?>
<?	end
?>		)
		- eta * U->B_U.<?=xi?>
<?	for j,xj in ipairs(xNames) do
?>		+ partial_B_ul.<?=xi?>.<?=xj?> * beta_u.<?=xj?> * coord_dx<?=i-1?>(x)
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

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>
<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>

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
	real rescaleMetric = cbrt(1. / det_gammaBar_over_det_gammaHat);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_ll.<?=xij?> *= rescaleMetric;
<? 		end ?>
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);
<?	end	-- constrain_det_gammaBar ?>

	//these are now based on the adjusted epsilon_LL:
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	sym3 gammaBar_LL = sym3_rescaleFromCoord_ll(gammaBar_ll, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
	
	//in Buchman's paper it says he doesn't do this
	//and in the new arbitrary-coord formalism, there is a tr ABar_ij term
<? if eqn.guiVars.constrain_tr_ABar.value then ?>
	U->ABar_LL = tracefree(U->ABar_LL, gammaBar_LL, gammaBar_UU);
<? end	-- constrain_tr_ABar ?>

<? else -- constrain_det_gammaBar or constrain_tr_ABar ?>
<?=assign_sym3'gammaBar_uu'?>
	
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
<? end -- constrain_det_gammaBar or constrain_tr_ABar ?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

//TODO these need to be pre-scaled back to coordinates before computing the weighted finite difference
<?=eqn:makePartial'ABar_LL'?>		//partial_ABar_LLl[k].IJ = ABar_IJ,k
<?=eqn:makePartial'K'?>				//partial_K_l[i] := K,i
<?=eqn:makePartial'W'?>				//partial_W_l[i] := phi_,i 
<?=eqn:makePartial'LambdaBar_U'?>	//partial_LambdaBar_Ul[j].I := LambdaBar^I_,j
<?=eqn:makePartial2'W'?>				//partial2_W_ll[ij] := phi_,ij

<?=assign_real3'partial_phi_l'?>
<?=assign_sym3'partial2_phi_ll'?>
	
	real exp_neg4phi = calc_exp_neg4phi(U);
	
<?=assign_3sym3'connHat_ull'?>
	
<?=eqn:makePartial'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=assign_3sym3'connBar_ull'?>

	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar_uu.ij := ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	
	
<?=assign_3sym3'Delta_ull'?>	
<?=assign_3sym3'Delta_lll'?>
<?=assign_real3'Delta_u'?>
<?=assign_real3'LambdaBar_u'?>	
<?=assign_sym3sym3'partial2_gammaHat_llll'?>
<?=assign_sym3'gammaHat_uu'?>	
<?=assign_real3'partial_det_gammaHat_l'?>
<?=assign_3sym3x3'partial_connHat_ulll'?>

<?=assign_real3x3'partial_LambdaBar_ul'?>

	sym3 RBar_ll;
	{
<?=eqn:makePartial2'epsilon_LL'?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>
		
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
	
<?=assign_3sym3('partial_ABar_lll', partial_ABar_lll:permute'_kij')?>
<?=assign_3sym3('partial_epsilon_lll', partial_epsilon_lll:permute'_kij')?>
	
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
?>		- ABar_uu.<?=sym(k,j)?> * gammaBar_uu.<?=sym(l,i)?> * partial_epsilon_lll.<?=xj?>.<?=sym(k,l)?>
		- ABar_uu.<?=sym(k,i)?> * gammaBar_uu.<?=sym(l,j)?> * partial_epsilon_lll.<?=xj?>.<?=sym(k,l)?>
		+ gammaBar_uu.<?=sym(k,i)?> * gammaBar_uu.<?=sym(l,j)?> * partial_ABar_lll.<?=xj?>.<?=sym(k,l)?>
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
<?
-- [[ do this every time you stop using the env
setfenv(1, oldEnv)
--]]
?>
