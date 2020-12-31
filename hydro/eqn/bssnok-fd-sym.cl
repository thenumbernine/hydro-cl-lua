//// MODULE_NAME: <?=calc_gammaHat_ll?>

sym3 <?=calc_gammaHat_ll?>(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>
	return gammaHat_ll;
}

//// MODULE_NAME: <?=calc_gammaHat_uu?>

sym3 <?=calc_gammaHat_uu?>(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_uu'?>
	return gammaHat_uu;
}

//// MODULE_NAME: <?=calc_det_gammaHat?>

real <?=calc_det_gammaHat?>(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaHat'?>
	return det_gammaHat;
}

//// MODULE_NAME: <?=calc_gammaHat_LL?>
// also in parent
#define <?=calc_gammaHat_LL?>(x) (sym3_ident)

//// MODULE_NAME: <?=calc_gammaHat_UU?>
// also in parent
#define <?=calc_gammaHat_UU?>(x) (sym3_ident)

//// MODULE_NAME: <?=calc_gammaBar_LL?>
//// MODULE_DEPENDS: <?=cons_t?>

sym3 <?=calc_gammaBar_LL?>(global const <?=cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaBar_LL'?>
	return gammaBar_LL;
}

//// MODULE_NAME: <?=calc_det_gammaBarLL?>
// also in parent ... at least the #if 1 part is
/*
det(epsilon_IJ + gammaHat_IJ) 
= det(epsilon_IJ + delta_IJ) 
= det(e^i_I e^j_J (gammaHat_ij + epsilon_ij))
= det(e^i_I) det(e^j_J) det(gammaHat_ij + epsilon_ij)
= det(gammaHat^ij) det(gammaBar_ij)
= det(gammaBar_ij) / det(gammaHat_ij)
= 1
TODO detg ... unless we want to change the constraint
*/
#if 0	//use the value
real <?=calc_det_gammaBarLL?>(global const <?=cons_t?>* U, ral3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
	return det_gammaBar_over_det_gammaHat;
}
#else	//use the constraint
#define <?=calc_det_gammaBarLL?>(x) 1.
#endif

//// MODULE_NAME: <?=calc_gammaBar_UU?>
//// MODULE_DEPENDS: <?=cons_t?>

sym3 <?=calc_gammaBar_UU?>(global const <?=cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<? -- assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>
	return gammaBar_UU;
}

//// MODULE_NAME: <?=calc_gammaBar_ll?>
//// MODULE_DEPENDS: <?=cons_t?>

//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 <?=calc_gammaBar_ll?>(global const <?=cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaBar_ll'?>
	return gammaBar_ll;
}

//// MODULE_NAME: <?=calc_det_gammaBar?>
//// MODULE_DEPENDS: <?=calc_det_gammaHat?>

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2017 Ruchlin
real <?=calc_det_gammaBar?>(real3 x) {
	//TODO detg ...
	real det_gammaHat = <?=calc_det_gammaHat?>(x);
	//real detg = 1.;
	//real det_gammaBar = det_gammaHat * detg;
	//return det_gammaBar;
	return det_gammaHat;
}

//// MODULE_NAME: <?=calc_exp_neg4phi?>
#define <?=calc_exp_neg4phi?>(U) ((U)->W * (U)->W)

//// MODULE_NAME: <?=calc_gammaBar_uu?>
//// MODULE_DEPENDS: <?=cons_t?> <?=calc_gammaBar_ll?> <?=calc_det_gammaBar?>

// also in parent
sym3 <?=calc_gammaBar_uu?>(global const <?=cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	real det_gammaBar = <?=calc_det_gammaBar?>(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}

//// MODULE_NAME: <?=calc_gamma_ll?>
//// MODULE_DEPENDS: <?=cons_t?> <?=calc_gammaBar_ll?> <?=calc_exp_neg4phi?>

// also in parent
sym3 <?=calc_gamma_ll?>(
	global <?=cons_t?> const * const U,
	real3 const x
) {
	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	return gamma_ll;
}
	
//// MODULE_NAME: <?=calc_gamma_uu?>
//// MODULE_DEPENDS: <?=cons_t?> <?=calc_gammaBar_ll?> <?=calc_exp_neg4phi?> <?=calc_det_gammaBar?>
// also in parent

sym3 <?=calc_gamma_uu?>(global const <?=cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	real det_gamma = <?=calc_det_gammaBar?>(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma); 
	return gamma_uu;
}

//// MODULE_NAME: <?=eqn_common?>

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
Then there's the possibility that, by moving the calc_RBar_LL code into its own function, maybe it uses too many arguments?
Then there's the possibility that I'm using too many locals and the compiler has a bug with stack allocations beyond a certain size.
*/

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
Returns the step coefficients, [+-1, +-1, +-1]
based on vector 'v'
to compute upwind differencing from.
Same thing as (int4)sgn(v)
*/
const int4 getUpwind(real3 v) {
	return (int4)(
		v.x >= 0 ? 1 : -1,
		v.y >= 0 ? 1 : -1,
		v.z >= 0 ? 1 : -1,
		0
	);
}

//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
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
	int4 const updir = getUpwind(U->beta_U);

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

<?=eqn:makePartial1'alpha'?>			//partial_alpha_l.i := alpha_,i
<?=eqn:makePartial1'beta_U'?>			//partial_beta_Ul.j.I := beta^I_,j
<?=eqn:makePartial1'epsilon_LL'?>		//partial_epsilon_LLl[k].IJ := epsilon_IJ,k
<?=eqn:makePartial1'W'?>				//partial_W_l.i := W_,i 
<?=eqn:makePartial1'K'?>				//partial_K_l.i := K,i
<?=eqn:makePartial1'ABar_LL'?>		//partial_ABar_LLl[k].IJ = ABar_IJ,k
<?=eqn:makePartial1'LambdaBar_U'?>	//partial_LambdaBar_Ul.j.I := LambdaBar^I_,j
<?=eqn:makePartial2'alpha'?>			//partial2_alpha_ll.ij := alpha_,ij
<?=eqn:makePartial2'beta_U'?>		//partial2_beta_Ull[jk].I = beta^I_,jk
<?=eqn:makePartial2'W'?>				//partial2_W_ll.ij := W_,ij
	//partial_B[i] := B^i_,t
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=eqn:makePartial1'B_U'?>			//partial_B_Ul.j.I := B^I_,j
<? end ?>

<?=assign_sym3'gammaHat_ll'?>

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
<? -- assign'det_gammaBar_over_det_gammaHat'?>


	//////////////////////////////// alpha_,t //////////////////////////////// 

<?=eqn:makePartialUpwind'alpha'?>
<?=assign'dt_alpha'?>
	deriv->alpha += dt_alpha;
	
	//////////////////////////////// W_,t //////////////////////////////// 

<?=eqn:makePartialUpwind'W'?>
<?=assign_real3'partial_det_gammaBar_over_det_gammaHat_l'?>
<?=assign'dt_W'?>
	deriv->W += dt_W;

	//////////////////////////////// K_,t //////////////////////////////// 
	
<?=eqn:makePartialUpwind'K'?>
<?=assign'dt_K'?>
	deriv->K += dt_K;

	//////////////////////////////// epsilon_ij,t //////////////////////////////// 

<?=eqn:makePartialUpwind'epsilon_LL'?>
<?=assign_sym3'dt_epsilon_LL'?>
	deriv->epsilon_LL = sym3_add(deriv->epsilon_LL, dt_epsilon_LL);

	//////////////////////////////// ABar_ij,t //////////////////////////////// 

<?=eqn:makePartial2'epsilon_LL'?>
<?=eqn:makePartialUpwind'ABar_LL'?>

<?=assign_sym3'Delta_LLL'?>

<?=assign_sym3'dt_ABar_LL'?>
	deriv->ABar_LL = sym3_add(deriv->ABar_LL, dt_ABar_LL);

	//////////////////////////////// LambdaBar^i_,t //////////////////////////////// 

<?=eqn:makePartialUpwind'LambdaBar_U'?>
<?=assign_real3'dt_LambdaBar_U'?>
	deriv->LambdaBar_U = real3_add(deriv->LambdaBar_U, dt_LambdaBar_U);

	//////////////////////////////// beta^i_,t and B^i_,t //////////////////////////////// 

<? if eqn.useShift == 'GammaDriver' then ?>

	const real k = 3. / 4.;
<?=assign_real3'dt_beta_U_GammaDriver'?>
	deriv->beta_U = real3_add(deriv->beta_U, dt_beta_U_GammaDriver);

<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>

<?=eqn:makePartialUpwind'beta_U'?>
<?=assign_real3'dt_beta_U_HyperbolicGammaDriver'?>
	deriv->beta_U = real3_add(deriv->beta_U, dt_beta_U_HyperbolicGammaDriver);

<?=eqn:makePartialUpwind'B_U'?>
<?=assign_real3'dt_B_U_HyperbolicGammaDriver'?>
	deriv->B_U = real3_add(deriv->B_U, dt_B_U_HyperbolicGammaDriver);

<? end	-- eqn.useShift ?>
}

kernel void constrainU(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useConstrainU then ?>	
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const U = UBuf + index;

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

<? 
if eqn.guiVars.constrain_det_gammaBar.value 
or eqn.guiVars.constrain_tr_ABar.value 
then 
?>

<?=assign'det_gammaBar_over_det_gammaHat'?>
	
	/*
	we need to force det(gammaBar_ij) = det(gammaHat_ij)
	and we can do so with gammaBar_ij := gammaBar_ij * (det(gammaHat_ij) / det(gammaBar_ij))^(1/3)
	so we find
	det(gammaBar_ij * (det(gammaHat_ij) / det(gammaBar_ij))^(1/3) )
	= det(gammaBar_ij) * (det(gammaHat_ij) / det(gammaBar_ij))
	= det(gammaHat_ij)
	*/
<?	if eqn.guiVars.constrain_det_gammaBar.value then ?>
<?=assign_sym3'gammaBar_ll'?>
	real rescaleMetric = cbrt(1. / det_gammaBar_over_det_gammaHat);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_ll.<?=xij?> *= rescaleMetric;
<? 		end ?>

<?=assign_sym3'gammaHat_ll'?>
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);
<?	end	-- constrain_det_gammaBar ?>

	//in Buchman's paper it says he doesn't do this
	//and in the new arbitrary-coord formalism, there is a tr ABar_ij term
<? if eqn.guiVars.constrain_tr_ABar.value then ?>
<?=assign_sym3'gammaBar_LL'?>
<?=assign_sym3'gammaBar_UU'?>
	U->ABar_LL = tracefree(U->ABar_LL, gammaBar_LL, gammaBar_UU);
<? end	-- constrain_tr_ABar ?>

<? else -- constrain_det_gammaBar or constrain_tr_ABar ?>

<? -- assign'det_gammaBar_over_det_gammaHat'?>
	
<? end -- constrain_det_gammaBar or constrain_tr_ABar ?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

//TODO these need to be pre-scaled back to coordinates before computing the weighted finite difference
<?=eqn:makePartial1'ABar_LL'?>			//partial_ABar_LLl[k].IJ = ABar_IJ,k
<?=eqn:makePartial1'alpha'?>			//partial_alpha_l.i := alpha_,i
<?=eqn:makePartial1'K'?>				//partial_K_l.i := K_,i
<?=eqn:makePartial1'W'?>				//partial_W_l.i := phi_,i 
<?=eqn:makePartial1'LambdaBar_U'?>		//partial_LambdaBar_Ul.j.I := LambdaBar^I_,j
<?=eqn:makePartial2'W'?>				//partial2_W_ll.ij := phi_,ij
<?=eqn:makePartial1'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial2'epsilon_LL'?>

<?=assign_3sym3'connBar_LLL'?>

<?=assign'H_def'?>
	U->H = H_def;

<?=assign_real3'M_U_def'?>
	U->M_U = M_U_def;

<? end	-- calc_H_and_M ?>
<? end	-- useConstrainU ?>
}

kernel void addSource(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useAddSource then ?>
	<?=SETBOUNDS_NOGHOST?>();
	
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const deriv = derivBuf + index;

	if (solver->diffuseCoeff != 0.) { 
		//Kreiss-Oligar dissipation
		//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
		//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
		for (int i = 0; i < numIntStates; ++i) {
<?=require 'hydro.eqn.makepartial'.makePartialRank1(4, 4, solver, 'ptr[i]', 'real', 'partial4_Ui_ll')?>
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
