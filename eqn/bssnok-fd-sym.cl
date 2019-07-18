<? 
-- constrains det gammaBar_ij = det gammaHat_ij, ABar^i_i = 0, and calculates H and M^i ... if the associated flags are set
local useConstrainU = true

-- does Kreiss-Oligar dissipation
local useAddSource = true
?>

<? if getCommonCode then ?>

static sym3 calc_RBar_LL(
	const global <?=eqn.cons_t?>* U,
	real3 x,
	const sym3* gammaBar_ll,
	const sym3* gammaBar_uu,
	const sym3* gammaBar_UU,
	const _3sym3* connHat_ull,
	const _3sym3* partial_gammaBar_lll,
	const sym3* trBar_partial2_gammaBar_ll,
	const real3x3* partial_LambdaBar_ul,
	const real3* Delta_U,
	const _3sym3* Delta_ULL,
	const _3sym3* Delta_LLL,
	const _3sym3 partial_connHat_ulll[3],
	const real3* LambdaBar_u
) {
	//DHat_gammaBar_lll.k.ij = DHat_k gammaBar_ij 
	// = gammaBar_ij,k - connHat^l_ki gammaBar_lj - connHat^l_kj gammaBar_il
	_3sym3 DHat_gammaBar_lll;
<?
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>	DHat_gammaBar_lll.<?=xk?>.<?=xij?> = partial_gammaBar_lll-><?=xk?>.<?=xij?>
<?		for l,xl in ipairs(xNames) do
?>		- connHat_ull-><?=xl?>.<?=sym(k,i)?> * gammaBar_ll-><?=sym(l,j)?>
		- connHat_ull-><?=xl?>.<?=sym(k,j)?> * gammaBar_ll-><?=sym(l,i)?>
<?		end
?>	;
<?	end
end
?>

	/*
	partial_DHat_gammaBar_without_partial2_gammaBar_llll[l].k.ij := partial_l DHat_k gammaBar_ij
		= (gammaBar_ij,k - connHat^m_ki gammaBar_mj - connHat^m_kj gammaBar_mi)_,l
		= gammaBar_ij,kl 
			- connHat^m_ki,l gammaBar_mj 
			- connHat^m_kj,l gammaBar_mi
			- connHat^m_ki gammaBar_mj,l
			- connHat^m_kj gammaBar_mi,l
	*/
	_3sym3 partial_DHat_gammaBar_without_partial2_gammaBar_llll[3];
<? 
for k,xk in ipairs(xNames) do
	for l,xl in ipairs(xNames) do
		for ij,xij in ipairs(symNames) do
			local i,j,xi,xj = from6to3x3(ij)
?>	partial_DHat_gammaBar_without_partial2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
//		+ partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>	// moved into the above gen'd code 
<?			for m,xm in ipairs(xNames) do
?>		- partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,i)?> * gammaBar_ll-><?=sym(m,j)?>		//diverging: RBar_ij = diag(0, -.5, -.75)
		- partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,j)?> * gammaBar_ll-><?=sym(m,i)?>		//diverging: RBar_ij = diag(0, -.5, -.75)
		- connHat_ull-><?=xm?>.<?=sym(k,i)?> * partial_gammaBar_lll-><?=xl?>.<?=sym(m,j)?>		//diverging: RBar_ij = diag(0, 1, 1)
		- connHat_ull-><?=xm?>.<?=sym(k,j)?> * partial_gammaBar_lll-><?=xl?>.<?=sym(m,i)?>		//diverging: RBar_ij = diag(0, 1, 1)
<?			end
?>	;
<?		end
	end
end
?>

	/*
	DHat2_gammaBar_llll_minus_one_term[l].k.ij = DHat_l DHat_k gammaBar_ij
		= partial_l DHat_k gammaBar_ij
			- connHat^m_lk DHat_m gammaBar_ij
			- connHat^m_li DHat_k gammaBar_mj
			- connHat^m_lj DHat_k gammaBar_im
	*/
	_3sym3 DHat2_gammaBar_llll_minus_one_term[3];
<?
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>	DHat2_gammaBar_llll_minus_one_term[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
		+ partial_DHat_gammaBar_without_partial2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>	//diverging, such that in spherical vacuum the trace of these is diag(0, 1, .5)
<?			for m,xm in ipairs(xNames) do
?>		- connHat_ull-><?=xm?>.<?=sym(l,k)?> * DHat_gammaBar_lll.<?=xm?>.<?=sym(i,j)?>
		- connHat_ull-><?=xm?>.<?=sym(l,i)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(m,j)?>
		- connHat_ull-><?=xm?>.<?=sym(l,j)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(i,m)?>
<?			end
?>	;
<?		end
	end
end
?>
	
	//trBar_DHat2_gammaBar_ll.ij := gammaBar^kl DHat_k DHat_l gammaBar_ij
	sym3 trBar_DHat2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	trBar_DHat2_gammaBar_ll.<?=xij?> = 0.
		+ trBar_partial2_gammaBar_ll-><?=xij?>
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu-><?=sym(k,l)?> * DHat2_gammaBar_llll_minus_one_term[<?=l-1?>].<?=xk?>.<?=xij?>
<?		end
	end
?>	;
<? end
?>

	//derivative is the last index, unlike the partial_*'s
	//DHat_LambdaBar_ul.i.j := DHat_j LambdaBar^i = LambdaBar^i_,j + connHat^i_jk LambdaBar^k
	real3x3 DHat_LambdaBar_ul;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	DHat_LambdaBar_ul.<?=xi?>.<?=xj?> = partial_LambdaBar_ul-><?=xi?>.<?=xj?>
<?		for k,xk in ipairs(xNames) do
?>		+ connHat_ull-><?=xi?>.<?=sym(j,k)?> * LambdaBar_u-><?=xk?>
<?		end
?>	;
<?	end
end
?>

	/*
	2017 Ruchlin eqn 12
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
	
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ikj
		+ 1/2 Delta^k Delta_jki
		+ Delta^m_ki Delta_jm^k
		+ Delta^m_kj Delta_im^k
		+ Delta^m_ik Delta_mj^k
	*/
	sym3 RBar_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	RBar_LL.<?=xij?> = (0.			
			- .5 * trBar_DHat2_gammaBar_ll.<?=xij?>	
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * gammaBar_ll-><?=sym(i,k)?> * DHat_LambdaBar_ul.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_ll-><?=sym(j,k)?> * DHat_LambdaBar_ul.<?=xk?>.<?=xi?>
<?	end
?>	) / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x))
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * Delta_U-><?=xk?> * (
				Delta_LLL-><?=xi?>.<?=sym(k,j)?> 
				+ Delta_LLL-><?=xj?>.<?=sym(k,i)?>
			)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>			+ gammaBar_UU-><?=sym(k,l)?> * (0.
				+ Delta_ULL-><?=xm?>.<?=sym(k,i)?> * Delta_LLL-><?=xj?>.<?=sym(m,l)?>
				+ Delta_ULL-><?=xm?>.<?=sym(k,j)?> * Delta_LLL-><?=xi?>.<?=sym(m,l)?>
				+ Delta_ULL-><?=xm?>.<?=sym(i,k)?> * Delta_LLL-><?=xm?>.<?=sym(j,l)?>
			)
<?			end
		end
	end
?>	;
<? end ?>
	return RBar_LL;
}

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
<?=assign'det_gammaHat'?>
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
<?=assign'det_gammaHat'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
	return det_gammaBar_over_det_gammaHat;
}
#else	//use the constraint
#define calc_det_gammaBarLL(x) 1.
#endif

sym3 calc_gammaBar_UU(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaHat'?>
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
	U->M_u = real3_zero;
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
		RBar_LL = calc_RBar_LL(
			U,
			x,
			&gammaBar_ll,
			&gammaBar_uu,
			&gammaBar_UU,
			&connHat_ull,
			&partial_gammaBar_lll,
			&trBar_partial2_gammaBar_ll,
			&partial_LambdaBar_ul,
			&Delta_U,
			&Delta_ULL,
			&Delta_LLL,
			partial_connHat_ulll,
			&LambdaBar_u);
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

<?=assign_3sym3'partial2_beta_ull'?>
<?=assign_real3'tr12_partial2_beta_l'?>

<?=assign_sym3'partial2_det_gammaBar_over_det_gammaHat_ll'?>
<?=assign_sym3'partial_tr_connBar_ll'?>
<?=assign_real3'DBar_tr_DBar_beta_l'?>
<?=assign_real3'DBar_tr_DBar_beta_u'?>
<?=assign_real3'trBar_DHat2_beta_u'?>

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

<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>
	
	sym3 RBar_LL;
	{
<?=eqn:makePartial2'epsilon_LL'?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
		
		RBar_LL = calc_RBar_LL(
			U,
			x,
			&gammaBar_ll,
			&gammaBar_uu,
			&gammaBar_UU,
			&connHat_ull,
			&partial_gammaBar_lll,
			&trBar_partial2_gammaBar_ll,
			&partial_LambdaBar_ul,
			&Delta_U,
			&Delta_ULL,
			&Delta_LLL,
			partial_connHat_ulll,
			&LambdaBar_u);
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

<?=assign_sym3'ABar_ll'?>

	/*
	trBar_DHat_ABar_u.i := gammabar^ik DHat_j ABar_kl gammaBar^lj
	= (
		ABar_kl,j
		- ABar_lm connHat^m_kj
		- ABar_km connHat^m_lj
	) gammaBar^ik gammaBar^lj
	*/	
	real3 trBar_DHat_ABar_u;
<? for i,xi in ipairs(xNames) do
?>	trBar_DHat_ABar_u.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(i,k)?> * gammaBar_uu.<?=sym(j,l)?> * (0.
			+ partial_ABar_lll.<?=xj?>.<?=sym(k,l)?> 
<?				for m,xm in ipairs(xNames) do
?>			- ABar_ll.<?=sym(l,m)?> * connHat_ull.<?=xm?>.<?=sym(k,j)?>
			- ABar_ll.<?=sym(k,m)?> * connHat_ull.<?=xm?>.<?=sym(l,j)?>
<?				end
?>		)
<?			end
		end
	end
?>	;
<? end
?>

	//also used in calc_RBar_LL
	//DHat_gammaBar_lll.k.ij = DHat_k gammaBar_ij 
	// = gammaBar_ij,k - connHat^l_ki gammaBar_lj - connHat^l_kj gammaBar_il
	_3sym3 DHat_gammaBar_lll;
<?
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>	DHat_gammaBar_lll.<?=xk?>.<?=xij?> = partial_gammaBar_lll.<?=xk?>.<?=xij?>
<?		for l,xl in ipairs(xNames) do
?>		- connHat_ull.<?=xl?>.<?=sym(k,i)?> * gammaBar_ll.<?=sym(l,j)?>
		- connHat_ull.<?=xl?>.<?=sym(k,j)?> * gammaBar_ll.<?=sym(l,i)?>
<?		end
?>	;
<?	end
end
?>

<?=assign_real3'partial_det_gammaBar_over_det_gammaHat_l'?>
<?=assign_real3'tr_connBar_l'?>
	_3sym3 Delta_ull = _3sym3_rescaleToCoord_ULL(Delta_ULL, x);
	
	/*
	2013 Baumgarte et al, eqn 14
	M^i = exp(-4 phi) (
		1/sqrt(det gammaBar) DHat_j (sqrt(det gammaBar) ABar^ij)
		+ 6 ABar^ij phi_,j
		- 2/3 gammaBar^ij K_,j
		ABar^jk Delta^i_jk
	) - 8 pi S^i
	= exp(-4 phi) (
		connBar^k_jk ABar^ij
		+ DHat_j ABar_kl gammaBar^ik gammaBar^lj
		- (DHat_j gammaBar_mn) (
			gammaBar^im gammaBar^nk ABar_kl gammaBar^lj
			+ gammaBar^lm gammaBar^nj ABar_kl gammaBar^ki
		)
		+ 6 ABar^ij phi_,j
		- 2/3 gammaBar^ij K_,j
		ABar^jk Delta^i_jk
	) - 8 pi S^i
	= exp(-4 phi) (
		+ gammaBar^jk DHat_j ABar_kl gammaBar^li
		
		+ ABar^ij (connBar^k_jk + 6 phi_,j)
		- 2/3 gammaBar^ij K_,j
		
		+ Delta^i_jk ABar^jk
		
		- DHat_j gammaBar_kl (
			gammaBar^ik ABar^lj
			+ gammaBar^kj ABar^li
		)
	) - 8 pi S^i
	*/
<? for i,xi in ipairs(xNames) do
?>	U->M_u.<?=xi?> = 
		- 8. * M_PI * U->S_u.<?=xi?>
		+ exp_neg4phi * (
			+ trBar_DHat_ABar_u.<?=xi?>
<?	for j,xj in ipairs(xNames) do
?>			+ ABar_uu.<?=sym(i,j)?> * (tr_connBar_l.<?=xj?> + partial_phi_l.<?=xj?>)
			- 2./3. * gammaBar_uu.<?=sym(i,j)?> * partial_K_l.<?=xj?>
<?		for k,xk in ipairs(xNames) do
?>			+ Delta_ull.<?=xi?>.<?=sym(j,k)?> * ABar_uu.<?=sym(j,k)?> 
<?			for l,xl in ipairs(xNames) do
?>			+ DHat_gammaBar_lll.<?=xj?>.<?=sym(k,l)?> * (
				gammaBar_uu.<?=sym(i,k)?> * ABar_uu.<?=sym(l,j)?>
				+ gammaBar_uu.<?=sym(j,k)?> * ABar_uu.<?=sym(l,i)?>
			)
<?			end
		end
	end
?>		);
<? end ?>
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

<? end -- getCommonCode ?>
