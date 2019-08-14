<? 
-- constrains det gammaBar_ij = det gammaHat_ij, ABar^i_i = 0, and calculates H and M^i ... if the associated flags are set
local useConstrainU = true

-- does Kreiss-Oligar dissipation
local useAddSource = true
?>

<? if getCommonCode then ?>

<? local dim = solver.dim ?>

#if 1	
//rescaling, used for bssn finite-difference, but I am tempted to try it with other coordinate systems with singularities
//TODO for the initial conditions do this symbolically instead of numerically

//apply this to lower indexes to convert from coordinate metric to better metric
//apply this to upper indexes to convert from better metric to coordinate metric
real3 real3_rescaleFromCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x / coord_dx0(x),
		.y = v.y / coord_dx1(x),
		.z = v.z / coord_dx2(x),
	};
}
#define real3_rescaleToCoord_U real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_L(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y * coord_dx1(x),
		.z = v.z * coord_dx2(x),
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_L

sym3 sym3_rescaleFromCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleToCoord_UU sym3_rescaleFromCoord_ll

sym3 sym3_rescaleToCoord_LL(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleFromCoord_uu sym3_rescaleToCoord_LL

_3sym3 _3sym3_rescaleFromCoord_lll(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleToCoord_UUU _3sym3_rescaleFromCoord_lll

_3sym3 _3sym3_rescaleToCoord_LLL(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleFromCoord_uuu _3sym3_rescaleToCoord_LLL

_3sym3 _3sym3_rescaleFromCoord_ull(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * coord_dx<?=i-1?>(x) / (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}

_3sym3 _3sym3_rescaleToCoord_ULL(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)) / coord_dx<?=i-1?>(x),
<?	end
?>		},
<? end
?>	};
}

sym3sym3 sym3sym3_rescaleFromCoord_lll(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleToCoord_UUUU sym3sym3_rescaleFromCoord_llll

sym3sym3 sym3sym3_rescaleToCoord_LLLL(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleFromCoord_uuuu sym3sym3_rescaleToCoord_LLLL

#else	//debugging -- turning it off

#define real3_rescaleFromCoord_l(a,x) a
#define real3_rescaleToCoord_U(a,x) a
#define real3_rescaleToCoord_L(a,x) a
#define real3_rescaleFromCoord_u(a,x) a
#define sym3_rescaleFromCoord_ll(a,x) a
#define sym3_rescaleToCoord_UU(a,x) a
#define sym3_rescaleToCoord_LL(a,x) a
#define sym3_rescaleFromCoord_uu(a,x) a
#define _3sym3_rescaleFromCoord_lll(a,x) a
#define _3sym3_rescaleToCoord_UUU(a,x) a
#define _3sym3_rescaleToCoord_LLL(a,x) a
#define _3sym3_rescaleFromCoord_uuu(a,x) a
#define sym3sym3_rescaleFromCoord_lll(a,x) a
#define sym3sym3_rescaleToCoord_UUUU(a,x) a
#define sym3sym3_rescaleToCoord_LLLL(a,x) a
#define sym3sym3_rescaleFromCoord_uuuu (a,x) a

#endif


//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero


	// gammaHat_ij and co


sym3 calc_gammaHat_ll(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>
	return gammaHat_ll;
}

real calc_det_gammaHat(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaHat'?>
	return det_gammaHat;
}

sym3 calc_gammaHat_uu(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_uu'?>
	return gammaHat_uu;
}


	// gammaBar_IJ and co


static inline sym3 calc_gammaHat_LL(real3 x) {
	return sym3_ident;
}

static inline sym3 calc_gammaHat_UU(real3 x) {
	return sym3_ident;
}

sym3 calc_gammaBar_LL(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaBar_LL'?>
	return gammaBar_LL;
}

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
real calc_det_gammaBarLL(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
	return det_gammaBar_over_det_gammaHat;
}
#else	//use the constraint
#define calc_det_gammaBarLL(x) 1.
#endif

sym3 calc_gammaBar_UU(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>
	return gammaBar_UU;
}
	
	
	// gammaBar_ij and co


//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaBar_ll'?>
	return gammaBar_ll;
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2017 Ruchlin
real calc_det_gammaBar(real3 x) {
	//TODO detg ...
	real det_gammaHat = calc_det_gammaHat(x);
	//real detg = 1.;
	//real det_gammaBar = det_gammaHat * detg;
	//return det_gammaBar;
	return det_gammaHat;
}

sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

sym3 calc_gamma_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	return gamma_ll;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gamma_ll = calc_gamma_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	real det_gamma = calc_det_gammaBar(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma); 
	return gamma_uu;
}

void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->beta_U = real3_zero;
	U->epsilon_LL = sym3_zero;
	U->W = 1;
	U->K = 0;
	U->ABar_LL = sym3_zero;

	//LambdaBar^i = Delta^i + C^i = Delta^i_jk gammaBar^jk = (connBar^i_jk - connHat^i_jk) gammaBar^jk + C^i
	//but when space is flat we have connBar^i_jk = connHat^i_jk and therefore Delta^i_jk = 0, Delta^i = 0, and LambdaBar^i = 0
	U->LambdaBar_U = mystery_C_U;

<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_U = real3_zero;
<? end
?>

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_U = real3_zero;
}



<? else	-- getCommonCode ?>

/*
I through I would save on loc by replacing expanded summation with vector/matrix operators, but it just made things worse.
The operators are pass-by-value (because I have read for C/C++ at least that this is faster for float3's,
 and for GPU's it will probably be faster for anything that fits in a float4/float8 or whatever the architecture is designed for).
However I've also read that Intel OpenCL used to have bugs with pass-by-value.
Then there's the possibility that, by moving the calc_RBar_LL code into its own function, maybe it uses too many arguments?
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
<?=assign'det_gammaBar_over_det_gammaHat'?>


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
	constant solver_t* solver,
	global cons_t* UBuf
) {
<? if useConstrainU then ?>	
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global cons_t* U = UBuf + index;

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

<?=assign'det_gammaBar_over_det_gammaHat'?>
	
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

<?=assign'H_def'?>
	U->H = H_def;

<?=assign_real3'M_U_def'?>
	U->M_U = M_U_def;

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
local cflMethod = '2008 Alcubierre'
--local cflMethod = '2013 Baumgarte et al, eqn 32'
--local cflMethod = '2017 Ruchlin et al, eqn 53'
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
	if cflMethod == '2008 Alcubierre' then 
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

<? end -- getCommonCode ?>
