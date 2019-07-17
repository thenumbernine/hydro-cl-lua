<? 
-- constrains det gammaBar_ij = det gammaHat_ij, ABar^i_i = 0, and calculates H and M^i ... if the associated flags are set
local useConstrainU = true

-- does Kreiss-Oligar dissipation
local useAddSource = true
?>

/*
I through I would save on loc by replacing expanded summation with vector/matrix operators, but it just made things worse.
The operators are pass-by-value (because I have read for C/C++ at least that this is faster for float3's,
 and for GPU's it will probably be faster for anything that fits in a float4/float8 or whatever the architecture is designed for).
However I've also read that Intel OpenCL used to have bugs with pass-by-value.
Then there's the possibility that, by moving the getCode_RBar_LL code into its own function, maybe it uses too many arguments?
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

/*
returns the pointer to the lowerleft cell to compute upwind differencing from 
*/
const global cons_t* getUpwind(
	constant solver_t* solver,
	const global cons_t* U
) {
	const global real3* beta_U = &U->beta_U;	//don't lose track of our original beta_U
<? for i=1,solver.dim do
	local xi = xNames[i]
?>	if (beta_U-><?=xi?> < 0) {
		U -= solver->stepsize.<?=xi?>;
	}
<? end
?>	return U;
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
	const global cons_t* Uup = getUpwind(solver, U);

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

<?=eqn:makePartial1'alpha'?>			//partial_alpha_l.i := alpha_,i
<?=eqn:makePartial1'beta_U'?>			//partial_beta_Ul[j].I := beta^I_,j
<?=eqn:makePartial1'epsilon_LL'?>		//partial_epsilon_LLl[k].IJ := epsilon_IJ,k
<?=eqn:makePartial1'W'?>				//partial_W_l.i := W_,i 
<?=eqn:makePartial1'K'?>				//partial_K_l.i := K,i
<?=eqn:makePartial1'ABar_LL'?>		//partial_ABar_LLl[k].IJ = ABar_IJ,k
<?=eqn:makePartial1'LambdaBar_U'?>	//partial_LambdaBar_Ul[j].I := LambdaBar^I_,j
<?=eqn:makePartial2'alpha'?>			//partial2_alpha_ll.ij := alpha_,ij
<?=eqn:makePartial2'beta_U'?>		//partial2_beta_Ull[jk].I = beta^I_,jk
<?=eqn:makePartial2'W'?>				//partial2_W_ll.ij := W_,ij
	//partial_B[i] := B^i_,t
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=eqn:makePartial1'B_U'?>			//partial_B_Ul[j].I := B^I_,j
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
	
	because in 2013 Baumgarte et al, IIB last paragraph, they say they relax this constraint.
	
	TODO detg ...
	*/
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_LL'?>
<?=assign_sym3'gammaBar_UU'?>
<?=assign_sym3'gammaBar_uu'?>


	//////////////////////////////// alpha_,t //////////////////////////////// 

<?=eqn:makePartialUpwind'alpha'?>;
<?=assign'dt_alpha'?>
	deriv->alpha += dt_alpha;
	
	//////////////////////////////// W_,t //////////////////////////////// 

<?=eqn:makePartialUpwind'W'?>;

<?=assign_real3'partial_det_gammaBar_over_det_gammaHat_l'?>
<?=assign_real3'tr_connBar_l'?>

<?=assign'tr_DBar_beta'?>
<?=assign'dt_W'?>
	deriv->W += dt_W;

	//////////////////////////////// K_,t //////////////////////////////// 
	
	//exp(-4 phi)
	real exp_neg4phi = calc_exp_neg4phi(U);

	/*
	gammaBar_ij = exp(-4 phi) gamma_ij
	gammaBar^ij = exp(4 phi) gamma^ij
	gamma^ij = exp(-4 phi) gammaBar^ij
	S := S_ij gamma^ij = exp(-4 phi) gammaBar^ij S_ij 
	*/
	real S = exp_neg4phi * sym3_dot(U->S_ll, gammaBar_uu);

<? if true then ?>
<?=assign_real3x3'ABar_UL'?>
<?=assign_sym3'ABarSq_LL'?>
<?=assign_sym3'ABar_UU'?>
<?=assign'tr_ABarSq'?>
<? end ?>
<? if false then ?>
<?=assign_real3x3'ABar_ul'?>
<?=assign_sym3'ABarSq_ll'?>
<?=assign'tr_ABarSq'?>
<? end ?>
<? if false then ?>
<?=assign_real3x3'ABar_ul'?>
<?=assign_real3x3'ABarSq_ul'?>
<?=assign'tr_ABarSq'?>
<? end ?>

<?=eqn:makePartialUpwind'K'?>;
<?=assign_3sym3'connBar_lll'?>
<?=assign_3sym3'connBar_ULL'?>
<?=assign_sym3'DBar2_alpha_LL'?> 
<?=assign'trBar_DBar2_alpha'?>
<?=assign_real3'partial_phi_l'?>
<?=assign'dt_K'?>
	deriv->K += dt_K;

	//////////////////////////////// epsilon_ij,t //////////////////////////////// 

<?=eqn:makePartialUpwind'epsilon_LL'?>;
<?=assign_sym3'Lbeta_gammaBar_LL'?>
<?=assign_sym3'dt_epsilon_LL'?>
<? for ij,xij in ipairs(symNames) do
?>	deriv->epsilon_LL.<?=xij?> += dt_epsilon_LL.<?=xij?>;
<? end
?>

	//////////////////////////////// ABar_ij,t //////////////////////////////// 

<?=assign_sym3'partial2_phi_ll'?>
<?=assign_sym3sym3'partial2_gammaHat_llll'?> 
<?=assign_3sym3'Delta_ULL'?>
<?=assign_real3'Delta_U'?>
<?=assign_real3'LambdaBar_u'?>
<?=assign_real3x3'partial_LambdaBar_ul'?>

<?=assign_sym3'gammaBar_ll'?>
<?=eqn:makePartial2'epsilon_LL'?>
	sym3 RBar_LL;
	{
<?=eqn:makePartial1'epsilon_LL'?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>
<?=assign_3sym3'Delta_LLL'?>
		<?=eqn:getCode_RBar_LL()?> 
	}

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
	sym3 TF_DBar2_alpha_LL = tracefree(DBar2_alpha_LL, gammaBar_LL, gammaBar_UU);

	sym3 TF_RBar_LL = tracefree(RBar_LL, gammaBar_LL, gammaBar_UU);
	
<?=assign_sym3'DBar2_phi_LL'?>
<?=assign_real3'partial_phi_L'?>
<?=assign_real3'partial_alpha_L'?>
<?=assign_real3'partial_K_L'?>
<?=assign_sym3'tracelessPart_LL'?>
	tracelessPart_LL = tracefree(tracelessPart_LL, gammaBar_LL, gammaBar_UU);

<?=eqn:makePartialUpwind'ABar_LL'?>;
<?=assign_3sym3('partial_ABar_lll_upwind', partial_ABar_lll_upwind:permute'_kij')?>
<?=assign_real3x3'partial_beta_ul'?>
<?=assign_sym3'Lbeta_ABar_LL'?>

<?=assign_sym3'dt_ABar_LL'?>
<? for ij,xij in ipairs(symNames) do
?>	deriv->ABar_LL.<?=xij?> += dt_ABar_LL.<?=xij?>;
<? end
?>

	//////////////////////////////// LambdaBar^i_,t //////////////////////////////// 

	real3 beta_u = real3_rescaleToCoord_U(U->beta_U, x);
<?=assign_real3x3x3'DHat2_beta_ull'?>

<?=assign_3sym3'partial2_beta_ull'?>
<?=assign_real3'tr12_partial2_beta_l'?>

<?=assign_sym3'partial2_det_gammaBar_over_det_gammaHat_ll'?>
<?=assign_sym3'partial_tr_connBar_ll'?>
<?=assign_real3'DBar_tr_DBar_beta_l'?>
<?=assign_real3'DBar_tr_DBar_beta_u'?>
<?=assign_real3'tr_gammaBar_DHat2_beta_u'?>

<?=eqn:makePartialUpwind'LambdaBar_U'?>;
<?=assign_real3'Lbeta_LambaBar_U'?>
<?=assign_real3'dt_LambdaBar_U'?>

<? for i,xi in ipairs(xNames) do
?>	deriv->LambdaBar_U.<?=xi?> += dt_LambdaBar_U.<?=xi?>;
<? end
?>

	//////////////////////////////// beta^i_,t and B^i_,t //////////////////////////////// 

<? if eqn.useShift == 'GammaDriver' then ?>

	const real k = 3. / 4.;
<?=assign_real3'dt_beta_U_GammaDriver'?>
<? for i,xi in ipairs(xNames) do
?>	deriv->beta_U.<?=xi?> += dt_beta_U_GammaDriver.<?=xi?>;
<? end
?>

<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>

<?=eqn:makePartialUpwind'beta_U'?>;
<?=assign_real3'dt_beta_U_HyperbolicGammaDriver'?>
<? for i,xi in ipairs(xNames) do
?>	deriv->beta_U.<?=xi?> += dt_beta_U_HyperbolicGammaDriver.<?=xi?>;
<? end
?>

<?=eqn:makePartialUpwind'B_U'?>;
<?=assign_real3'dt_B_U_HyperbolicGammaDriver'?>
<? for i,xi in ipairs(xNames) do
?>	deriv->B_U.<?=xi?> += dt_B_U_HyperbolicGammaDriver.<?=xi?>;
<? end
?>

<? end	-- eqn.useShift ?>
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
<?=assign_sym3'gammaBar_LL'?>

<? 
if eqn.guiVars.constrain_det_gammaBar.value 
or eqn.guiVars.constrain_tr_ABar.value 
then 
?>

<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
	
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

<?=assign_sym3'gammaHat_ll'?>
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);
<?	end	-- constrain_det_gammaBar ?>

	//these are now based on the adjusted epsilon_LL:
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	gammaBar_LL = sym3_rescaleFromCoord_ll(gammaBar_ll, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
	
	//in Buchman's paper it says he doesn't do this
	//and in the new arbitrary-coord formalism, there is a tr ABar_ij term
<? if eqn.guiVars.constrain_tr_ABar.value then ?>
	U->ABar_LL = tracefree(U->ABar_LL, gammaBar_LL, gammaBar_UU);
<? end	-- constrain_tr_ABar ?>

<? else -- constrain_det_gammaBar or constrain_tr_ABar ?>

<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>
<?=assign_sym3'gammaBar_uu'?>
	
<? end -- constrain_det_gammaBar or constrain_tr_ABar ?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

//TODO these need to be pre-scaled back to coordinates before computing the weighted finite difference
<?=eqn:makePartial1'ABar_LL'?>			//partial_ABar_LLl[k].IJ = ABar_IJ,k
<?=eqn:makePartial1'alpha'?>			//partial_alpha_l.i := alpha_,i
<?=eqn:makePartial1'K'?>				//partial_K_l.i := K_,i
<?=eqn:makePartial1'W'?>				//partial_W_l.i := phi_,i 
<?=eqn:makePartial1'LambdaBar_U'?>		//partial_LambdaBar_Ul[j].I := LambdaBar^I_,j
<?=eqn:makePartial2'W'?>				//partial2_W_ll.ij := phi_,ij

<?=assign_real3'partial_phi_l'?>
<?=assign_sym3'partial2_phi_ll'?>
	
	real exp_neg4phi = calc_exp_neg4phi(U);
	
<?=assign_3sym3'connHat_ull'?>
	
<?=eqn:makePartial1'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=assign_3sym3'connBar_lll'?>
<?=assign_3sym3'connBar_ULL'?>
<?=assign_real3x3'ABar_UL'?>
<?=assign_sym3'ABar_UU'?>
	
<?=assign_3sym3'Delta_ULL'?>
<?=assign_3sym3'Delta_LLL'?>
<?=assign_real3'Delta_U'?>
<?=assign_real3'LambdaBar_u'?>	
<?=assign_sym3sym3'partial2_gammaHat_llll'?>
<?=assign_sym3'gammaHat_uu'?>
<?=assign_real3'partial_det_gammaHat_l'?>
<?=assign_3sym3x3'partial_connHat_ulll'?>

<?=assign_real3x3'partial_LambdaBar_ul'?>

	sym3 RBar_LL;
	{
<?=eqn:makePartial2'epsilon_LL'?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>
		
		<?=eqn:getCode_RBar_LL()?> 
	}	
	//RBar := RBar_ij gammaBar^ij
	real RBar = sym3_dot(gammaBar_UU, RBar_LL);

<?=assign_sym3'DBar2_phi_LL'?>
<?=assign'tr_DBar2_phi'?>

<?=assign_real3'partial_phi_L'?>
<?=assign_real3'partial_alpha_L'?>
	//2017 Ruchlin et al, eqn 46
	//H = 2/3 K^2 - ABar^ij ABar_ij + exp(-4 phi) (RBar - 8 DBar^i phi DBar_i phi - 8 gammaBar^ij DBar_i DBar_j phi)
	U->H = 2. / 3. * U->K * U->K
		- sym3_dot(U->ABar_LL, ABar_UU)
		+ exp_neg4phi * (
			RBar
			- 8. * real3_weightedLenSq(partial_phi_L, gammaBar_UU) 
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
		- 2./3. * exp_6phi * gammaBar_uu.<?=sym(i,j)?> * partial_K_l.<?=xj?>
<?		for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xi?>.<?=sym(j,k)?> * ABar_UU.<?=sym(j,k)?> * coord_dx<?=i-1?>(x)
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
	
	const global <?=eqn.cons_t?>* U = UBuf + index;
	global cons_t* deriv = derivBuf + index;

	if (solver->diffuseCoeff != 0.) { 
		//Kreiss-Oligar dissipation
		//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
		//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
		for (int i = 0; i < numIntStates; ++i) {
<?=require'eqn.makepartial'.makePartialRank1(4, 4, solver, 'ptr[i]', 'real', 'partial4_Ui_ll')?>
		real lap = 0<?
for j,xj in ipairs(xNames) do
?> + partial4_Ui_ll.<?=xj?><?
end
?>;
			deriv->ptr[i] += solver->diffuseCoeff * lap;
		}
	}
<? end -- useAddSource ?>
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
--local cflMethod = '2008 Alcubierre'
--local cflMethod = '2013 Baumgarte et al, eqn 32'
local cflMethod = '2017 Ruchlin et al, eqn 53'
?>

<? if cflMethod == '2008 Alcubierre' then
?>	sym3 gamma_uu = calc_gamma_uu(U, x);
<? elseif cflMethod == '2017 Ruchlin et al, eqn 53' then
?>	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
<? end 
?>
	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
<? 
if cflMethod == '2013 Baumgarte et al, eqn 32' then
	if side == 0 then 
?>		dt = (real)min(dt, solver->grid_dx.x);
<?	elseif side == 1 then 
?>		dt = (real)min(dt, .5 * solver->grid_dx.x * solver->grid_dx.y);
<? 	elseif side == 2 then 
?>		dt = (real)min(dt, .5 * solver->grid_dx.x * sin(.5 * solver->grid_dx.y) * solver->grid_dx.z);
<? 	end 
else
	if cflMethod == '2017 Alcubierre' then 
?>		//this is asserting alpha and W >0, which they should be
		real absLambdaMax = U->alpha * sqrt(gamma_uu.<?=sym(side+1,side+1)?>);
<? 	elseif cflMethod == '2017 Ruchlin et al, eqn 53' then 
?>		real absLambdaMax = sqrt(gammaBar_ll.<?=sym(side+1,side+1)?>);
<? 	end 

	if false then -- hmm, do we base our CFL on delta in coordinate, or delta in Cartesian?
?>		real dx = cell_dx<?=side?>(x); 
<? 	else
?>		real dx = solver->grid_dx.s<?=side?>;
<? 	end 
?>		dt = (real)min(dt, dx / absLambdaMax);
<?
end 
?>	}<? end ?>
	dtBuf[index] = dt;
}
