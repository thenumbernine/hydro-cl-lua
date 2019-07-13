<? 

-- integrates whatsoever.
local useCalcDeriv = true
local useCalcDeriv_alpha = true
local useCalcDeriv_W = true
local useCalcDeriv_K = true
local useCalcDeriv_epsilon_LL = true
local useCalcDeriv_ABar_LL = true
local useCalcDeriv_LambdaBar_U = true
local useCalcDeriv_beta_U = true

-- constrains det gammaBar_ij = det gammaHat_ij, ABar^i_i = 0, and calculates H and M^i ... if the associated flags are set
local useConstrainU = true

-- does Kreiss-Oligar dissipation
local useAddSource = true
?>

#define real3_add3(a,b,c)	real3_add(real3_add(a,b),c)
#define real3_add4(a,b,c,d)	real3_add(real3_add(a,b),real3_add(c,d))

#define sym3_add3(a,b,c)	sym3_add(sym3_add(a,b),c)
#define sym3_add4(a,b,c,d)	sym3_add(sym3_add(a,b),sym3_add(c,d))	

<?
-- this block is shared with other things like initCond
-- it will be pasted above the !getCommonCode block below, despite that being in an exclusive condition to this
-- I'm moving it here from the .lua inline multiline string so I can get .cl syntax highlighting
if getCommonCode then	
?>

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


	// alpha^2 f(alpha)

<?
local symmath = require 'symmath'
local Tensor = symmath.Tensor
do
	local alpha = symmath.var'alpha'
	local fGuiVar = solver.eqn.guiVars.f_eqn
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	local f = assert(loadstring([=[
local alpha, symmath = ...
local log = symmath.log
return alpha^2 * ]=]..fLuaCode))(alpha, symmath)
	f = symmath.clone(f)()
?>real calc_f_times_alphaSq(real alpha) {
	return <?=symmath.export.C(f)?>;
}
<?
end
?>

	// gammaHat_ij and co

/*
e_i^I = delta_i^I f_i is a diagonal matrix with f_i indexed function.  for spherical, f_i = diag(1,r,r sin(theta))
e^i_I = delta^i_I f^i is the inverse, so f^i = 1/(f_i)
I'm trying to keep the differentiations to an absolute minimum in the bssnok-fd-num files
coord_dx#(x) is the same as f_# 
*/
<?
local coord = solver.coord

local partial_len_ll = Tensor'_ij'
partial_len_ll['_ij'] = coord.lenExprs'_i,j'()				-- derivative is last
local partial2_len_lll = Tensor'_ijk'
partial2_len_lll['_ijk'] = partial_len_ll'_ij,k'():factorDivision()	-- derivative is last

-- this is used often enough
local partial_len_over_len_lll = Tensor('_ijk', function(i,j,k)
	return (partial_len_ll[i][j] / coord.lenExprs[k])()
end)
?>

#define calc_gammaHat_ll	coord_g_ll
#define calc_det_gammaHat 	coord_det_g
#define calc_gammaHat_uu 	coord_g_uu

<?
for i,xi in ipairs(xNames) do
?>#define calc_len_<?=xi?>	coord_dx<?=i-1?>
<? 
end

for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
-- TODO? make this based on (pt) param (just like coord_dx), and get eqn:compile to cooperate
?>#define calc_partial_len_<?=xi..xj?>(pt)	(<?=eqn:compile(partial_len_ll[i][j]):gsub('x%.', 'pt.')?>)
<? 	end
end

for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>#define calc_partial_len_over_len_<?=xi..xj..xk?>(pt)	(<?=eqn:compile(partial_len_over_len_lll[i][j][k]):gsub('x%.', 'pt.')?>)
<? 		end
	end
end

for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>#define calc_partial2_len_<?=xi..xj..xk?>(pt)	(<?=eqn:compile(partial2_len_lll[i][j][k]):gsub('x%.', 'pt.')?>)
<? 		end
	end
end
?>

/*
derivative index is last
e_i^I (T^M e^i_M)_,j e^j_J
= (T^I_,j - T^I f_I,j / f_I) f_j delta^j_J
*/
real3x3 real3x3_partial_rescaleFromCoord_Ul(real3 T_U, const real3 partial_T_Ul[3], real3 x) {
	real3x3 partial_T_UL;
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	partial_T_UL.<?=xi?>.<?=xj?> = (0.
		+ partial_T_Ul[<?=j-1?>].<?=xi?>
		- T_U.<?=xi?> * calc_partial_len_<?=xi..xj?>(x)  / calc_len_<?=xi?>(x) 
	) / calc_len_<?=xj?>(x);
<?	end
end ?>
	return partial_T_UL;
}

//calc_partial_connHat_Ulll_ijkl := e_i^I connHat^i_jk,l
<? 
local partial_connHat_ulll = Tensor'^i_jkl'
partial_connHat_ulll['^i_jkl'] = coord.Gamma_ull'^i_jk,l'()

local partial_connHat_Ulll = Tensor('^I_jkl', function(i,j,k,l)
	return (partial_connHat_ulll[i][j][k][l] * coord.lenExprs[i])()
end)
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		for l,xl in ipairs(xNames) do
?>#define calc_partial_connHat_Ulll_<?=xi..xjk..xl?>(pt) (<?=eqn:compile(partial_connHat_Ulll[i][j][k][l]):gsub('x%.', 'pt.')?>)
<?		end
	end
end

?>

<?
local det_gammaHat = coord.det_g 
local partial_det_gammaHat_l = Tensor('_i', function(i)
	return det_gammaHat:diff(coord.coords[i])()
end)
local partial_det_gammaHat_over_det_gammaHat_L = Tensor('_i', function(i)
	return (partial_det_gammaHat_l[i] / (det_gammaHat * coord.lenExprs[i]))()
end)
?>
real3 calc_partial_det_gammaHat_over_det_gammaHat_L(real3 x) {
	real3 partial_det_gammaHat_over_det_gammaHat_L;
<? 
for i,xi in ipairs(xNames) do
?>	partial_det_gammaHat_over_det_gammaHat_L.<?=xi?> = <?=eqn:compile(partial_det_gammaHat_over_det_gammaHat_L[i])?>;
<? end
?>	return partial_det_gammaHat_over_det_gammaHat_L;
}

<?
local partial2_det_gammaHat_ll = partial_det_gammaHat_l'_i,j'():factorDivision()
local partial2_det_gammaHat_over_det_gammaHat_LL = Tensor('_ij', function(i,j)
	return (partial2_det_gammaHat_ll[i][j] / (det_gammaHat * coord.lenExprs[i] * coord.lenExprs[j]))()
end)
?>
sym3 calc_partial2_det_gammaHat_over_det_gammaHat_LL(real3 x) {
	sym3 partial2_det_gammaHat_over_det_gammaHat_LL;
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_det_gammaHat_over_det_gammaHat_LL.<?=xij?> = <?=eqn:compile(partial2_det_gammaHat_over_det_gammaHat_LL[i][j])?>;
<?
end
?>	return partial2_det_gammaHat_over_det_gammaHat_LL;
}


	// gammaBar_IJ and co


static inline sym3 calc_gammaHat_LL(real3 x) {
	return sym3_ident;
}

static inline sym3 calc_gammaHat_UU(real3 x) {
	return sym3_ident;
}

sym3 calc_gammaBar_LL(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_LL = calc_gammaHat_LL(x);
	sym3 gammaBar_LL = sym3_add(gammaHat_LL, U->epsilon_LL);
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
#define calc_det_gammaBarLL(x) 1.

sym3 calc_gammaBar_UU(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	return gammaBar_UU;
}


	// gammaBar_ij and co


//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_ll = calc_gammaHat_ll(x);
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);
	return gammaBar_ll;
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//...except sometimes, according to 2012 Baumgarte et al, last paragraph of II B
real calc_det_gammaBar(real3 x) {
	//TODO detg ...
	real det_gammaHat = calc_det_gammaHat(x);
	real detg = 1.;
	real det_gammaBar = det_gammaHat * detg;
	return det_gammaBar;
}

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}

sym3 calc_gamma_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	return gamma_ll;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
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

static void calc_partial_ABar_LLL(
	_3sym3* partial_ABar_LLL,
	const global <?=eqn.cons_t?>* U,
	real3 x,
	const sym3 partial_ABar_LLl[3]
) {
	//derivative first
	//partial_ABar_lll.k.ij := ABar_ij,k
	//partial_ABar_LLL.K.IJ := e_i^I e_j^J e_k^K ABar_ij,k
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_ABar_LLL-><?=xk?>.<?=xij?> = (0.
		+ partial_ABar_LLl[<?=k-1?>].<?=xij?>
		+ (
			calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			+ calc_partial_len_over_len_<?=xj..xk..xj?>(x)
		) * U->ABar_LL.<?=xij?>
	) / calc_len_<?=xk?>(x);
<?	end
end
?>
}

static void calc_partial_gammaBar_LLL(
	_3sym3 *partial_gammaBar_LLL,
	const global <?=eqn.cons_t?>* U,
	real3 x,
	const sym3 partial_epsilon_LLl[3]
) {
	//derivative first
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	//partial_gammaBar_LLL.k.IJ := e_i^I e_j^J e_k^K gammaBar_ij,k
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_gammaBar_LLL-><?=xk?>.<?=xij?> = (0.
		+ partial_epsilon_LLl[<?=k-1?>].<?=xij?>
		+ (
			calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			+ calc_partial_len_over_len_<?=xj..xk..xj?>(x)
		) * (U->epsilon_LL.<?=xij?><?=i==j and ' + 1.' or ''?>)
	) / calc_len_<?=xk?>(x);
<?	end
end
?>
}

static void calc_connHat_LLL_and_ULL(
	_3sym3* connHat_LLL,
	_3sym3* connHat_ULL,
	const global <?=eqn.cons_t?>* U,
	real3 x
) {
	// TODO inline the partial_gammaHat_LLL's
	//derivative first
	//partial_gammaHat_lll.k.ij := gammaHat_ij,k
	//partial_gammaHat_LLL.K.IJ := e_i^I e_j^J e_k^K gammaHat_ij,k
	// = 2 (f_I,k) / (f_I f_k) delta_IJ delta^k_K
	_3sym3 partial_gammaHat_LLL;
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
		if i==j then
?>	partial_gammaHat_LLL.<?=xk?>.<?=xij?> = 2. * calc_partial_len_over_len_<?=xi..xk..xi?>(x) / calc_len_<?=xk?>(x);
<?		else
?>	partial_gammaHat_LLL.<?=xk?>.<?=xij?> = 0.;
<?		end
	end
end
?>

	/*
	connHat_lll.i.jk := connHat_ijk = 1/2 (gammaHat_ij,k + gammaHat_ik,j - gammaHat_jk,i)
	connHat_LLL.I.JK := e^i_I e^j_J e^k_K connHat_ijk
	= 
		(f_I,k) / (f_I f_k) delta_IJ delta^k_K
		+ (f_I,j) / (f_I f_j) delta_IK delta^j_J
		- (f_J,i) / (f_I f_j) delta_JK delta^i_I
	*/
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>	connHat_LLL-><?=xi?>.<?=xjk?> = .5 * (0.

<?		if i==j then 
?>			+ partial_gammaHat_LLL.<?=xk?>.<?=sym(i,j)?>
<?		end
		if i==k then
?>			+ partial_gammaHat_LLL.<?=xj?>.<?=sym(i,k)?>
<?		end
		if j==k then
?>			- partial_gammaHat_LLL.<?=xi?>.<?=sym(j,k)?>
<?		end
?>		);
<?	end
end
?>	
	//connHat_ull[i].jk := connHat^i_jk = gammaHat^il connHat_ljk
	// raised by gammaHat_IJ, which equals delta^IJ
	*connHat_ULL = *connHat_LLL;
}

static void calc_connBar_ULL(
	_3sym3* connBar_ULL,
	const _3sym3* partial_gammaBar_LLL,
	const sym3* gammaBar_UU
) {
	//connBar_lll.i.jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	//connBar_LLL.I.JK := e^i_I e^j_J e^k_K connBar_ijk
	_3sym3 connBar_LLL;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>	connBar_LLL.<?=xi?>.<?=xjk?> = .5 * (0.
		+ partial_gammaBar_LLL-><?=xk?>.<?=sym(i,j)?>
		+ partial_gammaBar_LLL-><?=xj?>.<?=sym(i,k)?>
		- partial_gammaBar_LLL-><?=xi?>.<?=sym(j,k)?>
	);
<?	end
end
?>	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	*connBar_ULL = sym3_3sym3_mul(*gammaBar_UU, connBar_LLL);
}

static sym3 calc_trBar_partial2_gammaBar_LL(
	const global <?=eqn.cons_t?>* U,
	real3 x,
	const sym3* gammaBar_UU,
	const sym3 partial_epsilon_LLl[3],
	const sym3 partial2_epsilon_LLll[6]
) {
	//gives accuracy errors if I try to inline this
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(*gammaBar_UU, x); 
	
	/*
	e^i_I e^j_J gammaBar_ij,kl gammaBar^kl
	= (
		(delta_IJ + epsilon_IJ) (
			f_I,kl f^I
			+ 2 f_I,k f_J,l f^I f^J
			+ f_J,kl f^J
		)
		+ 2 epsilon_IJ,l (f_I,k f^I + f_J,k f^J)
		+ epsilon_IJ,kl
	) f^K f^L gammaBar^KL
	*/
	sym3 trBar_partial2_gammaBar_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	trBar_partial2_gammaBar_LL.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
			local kl,xkl = from3x3to6(k,l)
?>		+ (
			+ (U->epsilon_LL.<?=xij?><?=i==j and ' + 1.' or ''?>) * (0.
				+ calc_partial2_len_<?=xi..xk..xl?>(x) / calc_len_<?=xi?>(x)
				+ calc_partial2_len_<?=xj..xk..xl?>(x) / calc_len_<?=xj?>(x)
				+ 2. * calc_partial_len_over_len_<?=xi..xk..xi?>(x) * calc_partial_len_over_len_<?=xj..xl..xj?>(x) 
			)
			+ 2. * partial_epsilon_LLl[<?=l-1?>].<?=xij?> * (0.
				+ calc_partial_len_over_len_<?=xi..xk..xi?>(x)
				+ calc_partial_len_over_len_<?=xj..xk..xj?>(x)
			)
			+ partial2_epsilon_LLll[<?=kl-1?>].<?=xij?>
		) * gammaBar_uu.<?=xkl?>
<?		end
	end
?>	;
<? end
?>
	return trBar_partial2_gammaBar_LL;
}

static sym3 calc_RBar_LL(
	const global <?=eqn.cons_t?>* U,
	real3 x,
	const sym3* gammaBar_LL,
	const sym3* gammaBar_UU,
	const _3sym3* connHat_ULL,
	const _3sym3* partial_gammaBar_LLL,
	const sym3* trBar_partial2_gammaBar_LL,
	const real3x3* partial_LambdaBar_UL,
	const real3* Delta_U,
	const _3sym3* Delta_ULL,
	const _3sym3* Delta_LLL
) {
	
	//DHat_gammaBar_lll.k.ij = DHat_k gammaBar_ij 
	// = gammaBar_ij,k - connHat^l_ki gammaBar_lj - connHat^l_kj gammaBar_il
	_3sym3 DHat_gammaBar_LLL;
<?
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>	DHat_gammaBar_LLL.<?=xk?>.<?=xij?> = 0.
		+ partial_gammaBar_LLL-><?=xk?>.<?=xij?>
<?		for l,xl in ipairs(xNames) do
?>		- connHat_ULL-><?=xl?>.<?=sym(k,i)?> * gammaBar_LL-><?=sym(l,j)?>
		- connHat_ULL-><?=xl?>.<?=sym(k,j)?> * gammaBar_LL-><?=sym(l,i)?>
<?		end
?>	;
<?	end
end
?>

	_3sym3 connHat_ull = _3sym3_rescaleToCoord_ULL(*connHat_ULL, x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	_3sym3 partial_gammaBar_lll = _3sym3_rescaleToCoord_LLL(*partial_gammaBar_LLL, x);

	//connHat_times_partial_gammaBar_llll[l].k.ij := connHat^m_ij gammaBar_mk,l
	_3sym3 connHat_times_partial_gammaBar_llll[3];
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>	connHat_times_partial_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
<?			for m,xm in ipairs(xNames) do
?>		+ connHat_ull.<?=xm?>.<?=xij?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,k)?>
<?			end
?>	;
<?		end
	end
end
?>	
	
	//gammaBar_times_partial_connHat_llll[l].i.jk = gammaBar_im connHat^m_jk,l
	_3sym3 gammaBar_times_partial_connHat_llll[3];
<? for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
		for l,xl in ipairs(xNames) do
?>	gammaBar_times_partial_connHat_llll[<?=l-1?>].<?=xi?>.<?=xjk?> = 0.
<?			for m,xm in ipairs(xNames) do
?>		+ gammaBar_LL-><?=sym(i,m)?> * calc_len_<?=xi?>(x) * calc_partial_connHat_Ulll_<?=xm..xjk..xl?>(x)
<?			end
?>	;
<?		end
	end
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
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>	partial_DHat_gammaBar_without_partial2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = (0.
		- gammaBar_times_partial_connHat_llll[<?=l-1?>].<?=xj?>.<?=sym(k,i)?>
		- gammaBar_times_partial_connHat_llll[<?=l-1?>].<?=xi?>.<?=sym(k,j)?>
		- connHat_times_partial_gammaBar_llll[<?=l-1?>].<?=xj?>.<?=sym(i,k)?>
		- connHat_times_partial_gammaBar_llll[<?=l-1?>].<?=xi?>.<?=sym(j,k)?>
	) * (gammaBar_UU-><?=sym(k,l)?> / (calc_len_<?=xk?>(x) * calc_len_<?=xl?>(x)));
<?		end
	end
end
?>

	//connHat^i_jk gammaBar^jk
	real3 trBar_connHat_U = _3sym3_sym3_dot23(*connHat_ULL, *gammaBar_UU);

	_3sym3 DHat_gammaBar_lll = _3sym3_rescaleToCoord_LLL(DHat_gammaBar_LLL, x);

	//DHat_m gammaBar_ij * connHat^m_kl
	sym3 DHat_gammaBar_dot1_connHat_ll = sym3_rescaleToCoord_LL(real3_3sym3_dot1(
		trBar_connHat_U,
		DHat_gammaBar_LLL), x);

	//DHat_gammaBar_dot12_connHat_llll[i].j.kl
	//DHat_i gammaBar_jm * connHat^m_kl
	_3sym3 DHat_gammaBar_dot12_connHat_llll[3];
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for kl,xkl in ipairs(symNames) do
			local k,l,xk,xl = from6to3x3(kl)
?>	DHat_gammaBar_dot12_connHat_llll[<?=i-1?>].<?=xj?>.<?=xkl?> = 0.
<?			for m,xm in ipairs(xNames) do
?>		+ DHat_gammaBar_lll.<?=xi?>.<?=sym(j,m)?> * connHat_ull.<?=xm?>.<?=xkl?>
<?			end
?>	;
<?		end
	end
end
?>
	
	/*
	trBar_DHat2_gammaBar_without_partial2_gammaBar_ll.ij = gammaBar^kl DHat_l DHat_k gammaBar_ij
		= gammaBar^kl (partial_l DHat_k gammaBar_ij
			- connHat^m_lk DHat_m gammaBar_ij
			- connHat^m_li DHat_k gammaBar_mj
			- connHat^m_lj DHat_k gammaBar_im)
	*/
	sym3 trBar_DHat2_gammaBar_without_partial2_gammaBar_LL;
<?
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	trBar_DHat2_gammaBar_without_partial2_gammaBar_LL.<?=xij?> = (0.
		//Less accurate if I rescale this in advance.
		//Less accurate if I rescale this and subtract it later, outside this variable
		- DHat_gammaBar_dot1_connHat_ll.<?=sym(i,j)?>
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		//Less accurate to store this in advance 
		+ partial_DHat_gammaBar_without_partial2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>
		//Less accurate if I combine these two separate gammaBar^ij * ... statements 
		
		//Less accurate to convert all these to non-coord.
		//I think that is because connHat^I_JK has a few more 1/r's than connHat^i_jk
		//Also slightly less accurate to store this in advance.
		+ (0.
			- DHat_gammaBar_dot12_connHat_llll[<?=k-1?>].<?=xi?>.<?=sym(j,l)?>
			- DHat_gammaBar_dot12_connHat_llll[<?=k-1?>].<?=xj?>.<?=sym(i,l)?>
		) * (gammaBar_UU-><?=sym(k,l)?> / (calc_len_<?=xk?>(x) * calc_len_<?=xl?>(x)))
<?		end
	end
?>	) / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x));
<?
end
?>

	//trBar_DHat2_gammaBar_ll.ij := gammaBar^kl DHat_k DHat_l gammaBar_ij
	sym3 trBar_DHat2_gammaBar_LL = sym3_add(
		*trBar_partial2_gammaBar_LL,
		trBar_DHat2_gammaBar_without_partial2_gammaBar_LL);

	//derivative is the last index, unlike the partial_*'s
	//DHat_LambdaBar_ul.i.j := DHat_j LambdaBar^i = LambdaBar^i_,j + connHat^i_jk LambdaBar^k
	real3x3 DHat_LambdaBar_UL = real3x3_add(
		*partial_LambdaBar_UL,
		real3_3sym3_dot2(
			U->LambdaBar_U,
			*connHat_ULL
	));

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
?>	RBar_LL.<?=xij?> = 0.
			- .5 * trBar_DHat2_gammaBar_LL.<?=xij?>	
<?	for k,xk in ipairs(xNames) do
?>			//Less accurate when I defer this
			+ .5 * gammaBar_LL-><?=sym(i,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_LL-><?=sym(j,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xi?>
			
			//Less accurate when deferring these calculations
			+ .5 * Delta_U-><?=xk?> * (
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
<? end
?>	return RBar_LL;
}

<?
else	-- getCommonCode -- this block is for the scheme 
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

notice that gamma'_ij -> f gamma_ij; gamma'^ij -> 1/f gamma'^ij will produce the same result
so feel free to use gammaBar, gammaHat, etc
*/
sym3 tracefree(sym3 A_ll, sym3 g_ll, sym3 g_uu) {
	real tr_A = sym3_dot(A_ll, g_uu);
	return sym3_sub(A_ll, sym3_real_mul(g_ll, tr_A / 3.));
}

static void calcDeriv_alpha(
	constant solver_t* solver,
	global cons_t* deriv,
	const global cons_t* U,
	const real3* partial_alpha_L
) {
	//Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 f(alpha) K + alpha,i beta^i
	deriv->alpha += -calc_f_times_alphaSq(U->alpha) * U->K
		+ real3_dot(*partial_alpha_L, U->beta_U);
}

static void calcDeriv_W(
	constant solver_t* solver,
	global cons_t* deriv,
	const global cons_t* U,
	real3 x,
	real tr_DBar_beta,
	const real partial_W_l[3]
) {
	real3 partial_W_L;
<? for i,xi in ipairs(xNames) do
?>	partial_W_L.<?=xi?> = partial_W_l[<?=i-1?>] / calc_len_<?=xi?>(x);
<? end
?>
	
	//2017 Ruchlin et al eqn 11c
	//W,t = 1/3 W (alpha K - beta^k connBar^j_kj - beta^k_,k) + beta^k W_,k
	deriv->W += (1. / 3.) * U->W * (U->alpha * U->K - tr_DBar_beta) 
		+ real3_dot(partial_W_L, U->beta_U);
}

static void calcDeriv_K(
	constant solver_t* solver,
	global cons_t* deriv,
	const global cons_t* U,
	const sym3* gammaBar_UU,
	const sym3* ABar_UU,
	const sym3* DBar2_alpha_LL,
	const real3* partial_alpha_L,
	const real3* partial_phi_L,
	const real3* partial_K_L,
	real exp_neg4phi,
	real S
) {
	//tr_ABarSq := ABar_ij ABar^ij = ABar_ij ABar_kl gammaBar^ik gammaBar^jl
	real tr_ABarSq = sym3_dot(U->ABar_LL, *ABar_UU);
	
	//tr_DBar2_alpha := gammaBar^ij DBar_i DBar_j alpha
	real tr_DBar2_alpha = sym3_dot(*gammaBar_UU, *DBar2_alpha_LL);
	
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
	deriv->K += 
		U->alpha * (U->K * U->K / 3. + tr_ABarSq)
		- exp_neg4phi * (
			tr_DBar2_alpha
			+ 2. * real3_weightedDot(
				*partial_phi_L,
				*partial_alpha_L,
				*gammaBar_UU)
		)
		+ real3_dot(*partial_K_L, U->beta_U)
		+ 4. * M_PI * U->alpha * (U->rho + S);
}

static sym3 sym3_Lbeta_LL(
	const sym3 T_LL,
	const _3sym3* partial_T_LLL,
	const global real3* beta_U,
	const real3x3* partial_beta_UL
) {
	sym3 beta_times_partial_T_LL = real3_3sym3_dot1(*beta_U, *partial_T_LLL);
	real3x3 T_times_partial_beta_LL = sym3_real3x3_mul(T_LL, *partial_beta_UL);
	sym3 sym_T_times_partial_beta_LL = sym3_from_real3x3(T_times_partial_beta_LL);
	sym3 Lbeta_T_LL = sym3_add(beta_times_partial_T_LL, sym_T_times_partial_beta_LL);
	return Lbeta_T_LL;
}

static void calcDeriv_epsilon_LL(
	constant solver_t* solver,
	global cons_t* deriv,
	const global cons_t* U,
	const sym3* gammaBar_LL,
	const _3sym3* partial_gammaBar_LLL,
	const real3x3* partial_beta_UL,
	const real3x3* ABar_UL,
	real tr_DBar_beta
) {
	real tr_ABar = real3x3_trace(*ABar_UL);

	sym3 Lbeta_gammaBar_LL = sym3_Lbeta_LL(*gammaBar_LL, partial_gammaBar_LLL, &U->beta_U, partial_beta_UL);

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
		+ 2/3 alpha ABar^k_k
		- 2 alpha ABar_ij

	Looks like the paper's notation DBar_i beta_j implies DBar_j (gammaBar_ik beta^k) = gammaBar_ki DBar_j beta^k
	...which is not DBar_j ( gamma_ik beta^k ), which was my interpretation of the definition of beta_i := gamma_ij beta^k
	*/
	deriv->epsilon_LL = sym3_add4(
		deriv->epsilon_LL,
		sym3_real_mul(*gammaBar_LL, 2. / 3. * (U->alpha * tr_ABar - tr_DBar_beta)),
		sym3_real_mul(U->ABar_LL, -2. * U->alpha), 
		Lbeta_gammaBar_LL);
}

#define sym3_add7(a,b,c,d,e,f,g) \
	sym3_add( \
		sym3_add( \
			sym3_add(a,b), \
			sym3_add(c,d) \
		), \
		sym3_add( \
			sym3_add(e,f), \
			g \
		) \
	)

static void calcDeriv_ABar_LL(
	constant solver_t* solver,
	global cons_t* deriv,
	const global cons_t* U,
	real3 x,
	const sym3* gammaBar_LL,
	const sym3* gammaBar_UU,
	const sym3* DBar2_alpha_LL,
	const real3* partial_alpha_L, 
	const real3* partial_phi_L, 
	const _3sym3* connBar_ULL,
	const real3x3* ABar_UL,
	const real3x3* partial_beta_UL,
	const real partial_W_l[3],
	real exp_neg4phi,
	real tr_DBar_beta,
	//for RBar_IJ:
	const sym3 partial_epsilon_LLl[3],
	const _3sym3* Delta_ULL,
	const real3* Delta_U,
	const _3sym3* connHat_ULL,
	const _3sym3* partial_gammaBar_LLL,
	const real3x3* partial_LambdaBar_UL
) {
<?=eqn:makePartial2'epsilon_LL'?>

	sym3 trBar_partial2_gammaBar_LL = calc_trBar_partial2_gammaBar_LL(
		U, 
		x, 
		gammaBar_UU, 
		partial_epsilon_LLl, 
		partial2_epsilon_LLll);

	_3sym3 Delta_LLL = sym3_3sym3_mul(*gammaBar_LL, *Delta_ULL);

	sym3 RBar_LL = calc_RBar_LL(
		U,
		x,
		gammaBar_LL,
		gammaBar_UU,
		connHat_ULL,
		partial_gammaBar_LLL,
		&trBar_partial2_gammaBar_LL,
		partial_LambdaBar_UL,
		Delta_U,
		Delta_ULL,
		&Delta_LLL);

<?=eqn:makePartial'ABar_LL'?>		//partial_ABar_lll[k].ij = ABar_ij,k
	_3sym3 partial_ABar_LLL;
	calc_partial_ABar_LLL(&partial_ABar_LLL, U, x, partial_ABar_LLl);

<?=eqn:makePartial2'W'?>			//partial2_W_ll.ij := W_,ij
	//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
	sym3 partial2_phi_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	partial2_phi_LL.<?=xij?> = .5 * (
		-partial2_W_ll[<?=ij-1?>] / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
		+ partial_W_l[<?=i-1?>] / calc_len_<?=xi?>(x)
			* partial_W_l[<?=j-1?>] / calc_len_<?=xj?>(x)
			/ U->W
	) / U->W;
<? end ?>
	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_LL = sym3_sub(partial2_phi_LL, real3_3sym3_dot1(*partial_phi_L, *connBar_ULL));

	sym3 TF_DBar2_alpha_LL = tracefree(*DBar2_alpha_LL, *gammaBar_LL, *gammaBar_UU);

	sym3 TF_RBar_LL = tracefree(RBar_LL, *gammaBar_LL, *gammaBar_UU);

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
	sym3 tracelessPart_LL = tracefree(
		sym3_add4(
			sym3_real_mul(sym3_from_real3x3(real3_real3_outer(*partial_phi_L, *partial_alpha_L)), 2.),
			sym3_real_mul(DBar2_phi_LL, -2. * U->alpha),
			sym3_real_mul(real3_outer(*partial_phi_L), 4. * U->alpha),
			sym3_real_mul(sym3_rescaleFromCoord_ll(U->S_ll, x), -8. * M_PI * U->alpha)
		), 
		*gammaBar_LL, 
		*gammaBar_UU);

	sym3 ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, *ABar_UL);

	sym3 Lbeta_ABar_LL = sym3_Lbeta_LL(
		U->ABar_LL,
		&partial_ABar_LLL,
		&U->beta_U,
		partial_beta_UL);

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
	deriv->ABar_LL = sym3_add7(
		deriv->ABar_LL,
		sym3_real_mul(U->ABar_LL, U->alpha * U->K - 2. / 3. * tr_DBar_beta),
		sym3_real_mul(ABarSq_LL, -2. * U->alpha),
		sym3_real_mul(tracelessPart_LL, exp_neg4phi),
		sym3_real_mul(TF_DBar2_alpha_LL, -exp_neg4phi), 
		sym3_real_mul(TF_RBar_LL, U->alpha * exp_neg4phi),
		Lbeta_ABar_LL
	);
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

	//////////////////////////////// alpha_,t //////////////////////////////// 

	real3 partial_alpha_L;
	{
<?=eqn:makePartial'alpha'?>			//partial_alpha_l[i] := alpha_,i
<? for i,xi in ipairs(xNames) do
?>		partial_alpha_L.<?=xi?> = partial_alpha_l[<?=i-1?>] / calc_len_<?=xi?>(x);
<? end
?>	}

<? if useCalcDeriv_alpha then ?>
	calcDeriv_alpha(
		solver,
		deriv,
		U,
		&partial_alpha_L);
<? end	-- useCalcDeriv_alpha ?>

	//////////////////////////////// W_,t //////////////////////////////// 

<?=eqn:makePartial'beta_U'?>		//partial_beta_Ul[j].i := beta^I_,j
<?=eqn:makePartial'W'?>				//partial_W_l[i] := W_,i 

	/*
	notice: 
	det(gammaBar_ij)/det(gammaHat_ij) = det(gammaBar_IJ)
	*/
	real3 partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);

	/*
	Etienne's SENR Mathematica notebook has '*  detg'...
	 ... but the paper says det gammaBar_ij = det gammaHat_ij
	 ... so what is this mysterious 'detg' factor?
	 TODO find where detg comes from 
	
	From 2017 Ruchlin: 
	
	"We choose Ë†Î³ â‰¡ det(Ë†Î³ij ) = Â¯Î³ in the initial data for
	all applications in this paper. Since both determinants
	remain independent of time, they remain equal to each
	other throughout the evolution."
	
	... so why is there a distinction between det gammaBar_ij and det gammaHat_ij? 
	
	because in 2013 Baumgarte et al, IIB last paragraph, they say they relax this constraint.
	
	TODO detg ...
	*/
	real detg = 1.;
	real3 partial_detg_L = real3_zero;
	
	real3 partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);
	
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	real det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 

	//partial_beta_UL.I.J := e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = real3x3_trace(partial_beta_UL);

	/*
	tr_DBar_beta := DBar_j beta^j = beta^j_,j + connBar^j_kj beta^k
	Etienne's SENR Mathematica notebook uses this instead:
	tr_DBar_beta := beta^j_,j + beta^j gammaBar_,j / (2 gammaBar)
	*/
	real tr_DBar_beta = tr_partial_beta + real3_dot(U->beta_U, partial_det_gammaBar_over_det_gammaBar_L) * .5;

<? if useCalcDeriv_W then ?>
	calcDeriv_W(
		solver,
		deriv,
		U,
		x,
		tr_DBar_beta,
		partial_W_l);
<? end	-- useCalcDeriv_W ?>

	//////////////////////////////// K_,t //////////////////////////////// 

<?=eqn:makePartial'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial2'alpha'?>		//partial2_alpha_ll.ij := alpha_,ij

	real3 partial_K_L;
	{
<?=eqn:makePartial'K'?>			//partial_K_l[i] := K,i
<? for i,xi in ipairs(xNames) do
?>		partial_K_L.<?=xi?> = partial_K_l[<?=i-1?>] / calc_len_<?=xi?>(x);
<? end
?>	}

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

	//exp(-4 phi)
	real exp_neg4phi = calc_exp_neg4phi(U);

	real3 partial_phi_l;
	{
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>
	}
	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);

	_3sym3 partial_gammaBar_LLL;
	calc_partial_gammaBar_LLL(&partial_gammaBar_LLL, U, x, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	sym3 DBar2_alpha_LL;
	{
		sym3 partial2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		partial2_alpha_LL.<?=xij?> = partial2_alpha_ll[<?=ij-1?>] / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x));
<? end
?>
		//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
		//used below for TF_DBar2_alpha_LL
		DBar2_alpha_LL = sym3_sub(partial2_alpha_LL, real3_3sym3_dot1(partial_alpha_L, connBar_ULL));
	}

	//ABar^i_j := gammaBar^ik ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar^ij := ABar^i_k gammaBar^kj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	//ABar^IJ = ABar^I_K gammaBar^KJ

<? if useCalcDeriv_K then ?>
	
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);
	
	/*
	gammaBar_ij = exp(-4 phi) gamma_ij
	gammaBar^ij = exp(4 phi) gamma^ij
	gamma^ij = exp(-4 phi) gammaBar^ij
	S := S_ij gamma^ij = exp(-4 phi) S_ij gammaBar^ij 
	*/
	real S = exp_neg4phi * sym3_dot(U->S_ll, gammaBar_uu);

	calcDeriv_K(
		solver,
		deriv,
		U,
		&gammaBar_UU,
		&ABar_UU,
		&DBar2_alpha_LL,
		&partial_alpha_L,
		&partial_phi_L,
		&partial_K_L,
		exp_neg4phi,
		S
	);

<? end	-- useCalcDeriv_K ?>

	//////////////////////////////// epsilon_ij,t //////////////////////////////// 

<? if useCalcDeriv_epsilon_LL then ?>
	calcDeriv_epsilon_LL(
		solver,
		deriv,
		U,
		&gammaBar_LL,
		&partial_gammaBar_LLL,
		&partial_beta_UL,
		&ABar_UL,
		tr_DBar_beta
	);
<? end	-- useCalcDeriv_epsilon_LL ?>
	
	//////////////////////////////// ABar_ij_,t //////////////////////////////// 

	real3x3 partial_LambdaBar_UL;
	{
		//partial_LambdaBar_ul[j].i := LambdaBar^i_,j
<?=eqn:makePartial'LambdaBar_U'?>
		//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
		partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);
	}

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK + connHat^I_JK
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

<? if useCalcDeriv_ABar_LL then ?>
	
	//Delta_U.I := e_i^I Delta^i = e_i^I (LambdaBar^i - C^i)
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);

	calcDeriv_ABar_LL(
		solver,
		deriv,
		U,
		x,
		&gammaBar_LL,
		&gammaBar_UU,
		&DBar2_alpha_LL,
		&partial_alpha_L, 
		&partial_phi_L, 
		&connBar_ULL,
		&ABar_UL,
		&partial_beta_UL,
		partial_W_l,
		exp_neg4phi,
		tr_DBar_beta,
		//for RBar_IJ:
		partial_epsilon_LLl,
		&Delta_ULL,
		&Delta_U,
		&connHat_ULL,
		&partial_gammaBar_LLL,
		&partial_LambdaBar_UL
	);

<? end	-- useCalcDeriv_ABar_LL ?>
	
	//////////////////////////////// LambdaBar^i_,t //////////////////////////////// 

<? if useCalcDeriv_LambdaBar_U then ?>

	//partial2_beta_ULL.I.JK = e^i_I (beta^M e^i_M)_,jk e^j_J e^k_K
<?=eqn:makePartial2'beta_U'?>		
	_3sym3 partial2_beta_ULL;
<? for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>	partial2_beta_ULL.<?=xi?>.<?=xjk?> = (
		partial2_beta_Ull[<?=jk-1?>].<?=xi?> + (
			- partial_beta_Ul[<?=j-1?>].<?=xi?> * calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			- partial_beta_Ul[<?=k-1?>].<?=xi?> * calc_partial_len_over_len_<?=xi..xj..xj?>(x)
			+ U->beta_U.<?=xi?> * (
				calc_partial2_len_<?=xi..xj..xk?>(x) / calc_len_<?=xi?>(x)
				+ 2. * calc_partial_len_over_len_<?=xi..xj..xi?>(x) * calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			)
		)
	) / (calc_len_<?=xj?>(x) * calc_len_<?=xk?>(x));
<?	end
end ?>

	//for 1D this is r^2 sin(th) * (either sin(th) or cos(th))
	//in case it can be factored out somewhere
	sym3 partial2_det_gammaHat_over_det_gammaHat_LL = calc_partial2_det_gammaHat_over_det_gammaHat_LL(x);
	
	sym3 partial2_detg_LL = sym3_zero;

	sym3 partial2_det_gammaBar_over_det_gammaHat_LL = sym3_add3(
			sym3_real_mul(partial2_det_gammaHat_over_det_gammaHat_LL, detg),
			sym3_from_real3x3(real3_real3_outer(partial_det_gammaHat_over_det_gammaHat_L, partial_detg_L)),
			partial2_detg_LL);
	
	sym3 partial2_det_gammaBar_over_det_gammaBar_LL = sym3_real_mul(
		partial2_det_gammaBar_over_det_gammaHat_LL,
		det_gammaHat_over_det_gammaBar);
	
	/*
	DHat2_beta_ull.i.j.k = DHat_k DHat_j beta^i
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
	(This is breaking convention.  usually derivative indexes are left-most.)
	*/
	real3x3x3 DHat2_beta_ULL;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
			local jk,xjk = from3x3to6(j,k)
?>	DHat2_beta_ULL.<?=xi?>.<?=xj?>.<?=xk?> = 
		partial2_beta_ULL.<?=xi?>.<?=xjk?>
<?		for l,xl in ipairs(xNames) do
?>		+ calc_partial_connHat_Ulll_<?=xi..sym(l,j)..xk?>(x) * U->beta_U.<?=xl?> 
			/ (calc_len_<?=xj?>(x) * calc_len_<?=xk?>(x) * calc_len_<?=xl?>(x))
		+ connHat_ULL.<?=xi?>.<?=sym(l,j)?> * partial_beta_UL.<?=xl?>.<?=xk?>
		+ connHat_ULL.<?=xi?>.<?=sym(l,k)?> * partial_beta_UL.<?=xl?>.<?=xj?>
		- connHat_ULL.<?=xl?>.<?=sym(j,k)?> * partial_beta_UL.<?=xi?>.<?=xl?>
<?			for m,xm in ipairs(xNames) do
?>		+ connHat_ULL.<?=xi?>.<?=sym(m,k)?> * connHat_ULL.<?=xm?>.<?=sym(l,j)?> * U->beta_U.<?=xl?>
		- connHat_ULL.<?=xi?>.<?=sym(m,l)?> * connHat_ULL.<?=xm?>.<?=sym(k,j)?> * U->beta_U.<?=xl?>
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
	TODO simplify math plz
	*/
	real3 tr12_partial2_beta_L = _3sym3_tr12(partial2_beta_ULL);

	real3 partial2_det_gammaBar_over_det_gammaBar_times_beta_L = sym3_real3_mul(partial2_det_gammaBar_over_det_gammaBar_LL, U->beta_U);

	real partial_det_gammaBar_over_det_gammaBar_dot_beta = real3_dot(partial_det_gammaBar_over_det_gammaBar_L, U->beta_U);

	real3 partial_det_gammaBar_over_det_gammaBar_times_partial_beta_L = real3_real3x3_mul(partial_det_gammaBar_over_det_gammaBar_L, partial_beta_UL);

	/*
	DBar_tr_DBar_beta_l.i = DBar_i DBar_j beta^j
	= DBar_i (beta^j_,j + connBar^j_lj beta^l)
	= (beta^j_,j + connBar^j_lj beta^l)_,i
		+ connBar^j_mi (beta^m_,j + connBar^m_lj beta^l)
		- connBar^m_ji (beta^j_,m + connBar^j_lm beta^l)
	= beta^j_,ji 
		+ connBar^j_lj,i beta^l
		+ connBar^j_lj beta^l_,i

	...using determinant instead...

	= beta^j_,ji 
		+ 1/2 det_gammaBar_,l / det_gammaBar beta^l_,i
		+ (1/2 det_gammaBar_,l / det_gammaBar),i beta^l
	= beta^j_,ji + 1/2 (
		+ det_gammaBar_,l beta^l_,i / det_gammaBar
		+ det_gammaBar_,il beta^l / det_gammaBar
		- det_gammaBar_,l det_gammaBar_,i beta^l / det_gammaBar^2
	)

	Etienne's SENR uses this alternative formulation: 
	= beta^j_,ji
		+ 1/2 (
			gammaBar_,ij beta^j / gammaBar
			+ (gammaBar_,j / gammaBar) beta^j_,i 
			- (gammaBar_,i / gammaBar) (gammaBar_,j / gammaBar) beta^j
		)
	*/
	real3 DBar_tr_DBar_beta_L = real3_add(
		tr12_partial2_beta_L,
		real3_real_mul(
			real3_add3(
				partial2_det_gammaBar_over_det_gammaBar_times_beta_L,
				partial_det_gammaBar_over_det_gammaBar_times_partial_beta_L,
				real3_real_mul(partial_det_gammaBar_over_det_gammaBar_L, -partial_det_gammaBar_over_det_gammaBar_dot_beta)
			), .5)
	);
	
	//DBar_tr_DBar_beta_u.i = DBar^i DBar_k beta^k = gammaBar^ij DBar_j DBar_k beta^k
	real3 DBar_tr_DBar_beta_U = sym3_real3_mul(gammaBar_UU, DBar_tr_DBar_beta_L);

	//tr_gammaBar_DHat2_beta_u.i = gammaBar^jk DHat_j DHat_k beta^i
	real3 tr_gammaBar_DHat2_beta_U = real3x3x3_sym3_dot23(DHat2_beta_ULL, gammaBar_UU);

	real3 ABar_times_partial_alpha_U = sym3_real3_mul(ABar_UU, partial_alpha_L);
	real3 ABar_times_partial_phi_U = sym3_real3_mul(ABar_UU, partial_phi_L);

	real3 Delta_times_ABar_U = _3sym3_sym3_dot23(Delta_ULL, ABar_UU);

	real3 gammaBar_times_partial_K_U = sym3_real3_mul(gammaBar_UU, partial_K_L);

	real3 Lbeta_LambdaBar_U = real3_sub(
		real3x3_real3_mul(partial_LambdaBar_UL, U->beta_U),
		real3x3_real3_mul(partial_beta_UL, U->LambdaBar_U));

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

#define real3_add9(a,b,c,d,e,f,g,h,i) \
	real3_add( \
		real3_add( \
			real3_add( \
				real3_add(a,b), \
				c \
			), \
			real3_add(d,e) \
		), \
		real3_add( \
			real3_add(f,g), \
			real3_add(h,i) \
		) \
	)
	
	real3 dt_LambdaBar_U = real3_add9(
		tr_gammaBar_DHat2_beta_U,
		real3_real_mul(Delta_U, 2. / 3. * tr_DBar_beta),
		real3_real_mul(DBar_tr_DBar_beta_U, 1. / 3.), 
		real3_real_mul(ABar_times_partial_alpha_U, -2.),
		real3_real_mul(ABar_times_partial_phi_U, 12.),
		real3_real_mul(Delta_times_ABar_U, 2. * U->alpha),
		real3_real_mul(gammaBar_times_partial_K_U, -4. / 3. * U->alpha), 
		Lbeta_LambdaBar_U,
		real3_real_mul(
			real3_rescaleFromCoord_u(U->S_u, x),
			- 16. * M_PI * U->alpha / exp_neg4phi)
	);

	deriv->LambdaBar_U = real3_add(deriv->LambdaBar_U, dt_LambdaBar_U);

<? end	-- useCalcDeriv_LambdaBar_U ?>
	//////////////////////////////// beta^i_,t and B^i_,t //////////////////////////////// 

	//partial_B_ul[i] := B^i_,t
<? if eqn.useShift == 'HyperbolicGammaDriver' then ?>
<?=eqn:makePartial'B_U'?>
<? end ?>

<? if useCalcDeriv_beta_U then ?>

<? if eqn.useShift == 'GammaDriver' then ?>
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k (connBar^i_,t + eta connBar^i)
	const real k = 3. / 4.;
	const real eta = 1.;	//1.;	// 1 / (2 M), for total mass M
	deriv->beta_U = real3_add3(
		deriv->beta_U,
		real3_real_mul(dt_LambdaBar_U, k),
		real3_real_mul(U->LambdaBar_U, eta));
<? elseif eqn.useShift == 'HyperbolicGammaDriver' then ?>
	
	//partial_B_UL.I.J := e_i^I (B^M e^i_M)_,j e^j_J
	real3x3 partial_B_UL = real3x3_partial_rescaleFromCoord_Ul(U->B_U, partial_B_Ul, x);

	//SENR's 'shiftadvect':
	real3 partial_beta_times_beta = real3x3_real3_mul(partial_beta_UL, U->beta_U);

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
	deriv->beta_U = real3_add3(
		deriv->beta_U,
		U->B_U,
		partial_beta_times_beta 
	);

	//SENR's 'biadvect':
	real3 partial_beta_times_B = real3x3_real3_mul(partial_beta_UL, U->B_U);
	
	const real eta = 1.;
	deriv->B_U = real3_add4(
		deriv->B_U,
		real3_real_mul(dt_LambdaBar_U, .75),
		real3_real_mul(U->B_U, -eta),
		partial_beta_times_B 
	);

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
	
	sym3 gammaBar_LL = sym3_add(sym3_ident, U->epsilon_LL);
	real det_gammaBarLL = sym3_det(gammaBar_LL);

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
	
	const real det_gammaHatLL = 1.;
	real rescaleMetric = cbrt(det_gammaHatLL/det_gammaBarLL);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_LL.<?=xij?> *= rescaleMetric;
<? 		end ?>
	U->epsilon_LL = sym3_sub(gammaBar_LL, sym3_ident);
	det_gammaBarLL = det_gammaHatLL;

<?	end	-- constrain_det_gammaBar ?>

	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	
	//in Buchman's paper it says he doesn't do this
	//and in the new arbitrary-coord formalism, there is a tr ABar_ij term
<? if eqn.guiVars.constrain_tr_ABar.value then ?>
	U->ABar_LL = tracefree(U->ABar_LL, gammaBar_LL, gammaBar_UU);
<? end	-- constrain_tr_ABar ?>

<? else -- constrain_det_gammaBar or constrain_tr_ABar ?>
	
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	
<? end -- constrain_det_gammaBar or constrain_tr_ABar ?>

	//makes the spinning black hole simulations look ugly
	U->alpha = max(U->alpha, solver->alphaMin);

<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

<?=eqn:makePartial'K'?>				//partial_K_l[i] := K,i
<?=eqn:makePartial'W'?>				//partial_W_l[i] := phi_,i 
<?=eqn:makePartial2'W'?>			//partial2_W_ll.ij := phi_,ij

	real3 partial_phi_l;
<? for i,xi in ipairs(xNames) do
?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>
	
	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);

	sym3 partial2_phi_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	partial2_phi_LL.<?=xij?> = .5 * (
		-partial2_W_ll[<?=ij-1?>] / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
		+ partial_W_l[<?=i-1?>] / calc_len_<?=xi?>(x)
			* partial_W_l[<?=j-1?>] / calc_len_<?=xj?>(x)
			/ U->W
	) / U->W;
<? end ?>

	real exp_neg4phi = calc_exp_neg4phi(U);

	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial'epsilon_LL'?>	
	
	_3sym3 partial_gammaBar_LLL;
	calc_partial_gammaBar_LLL(&partial_gammaBar_LLL, U, x, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar_uu.ij := ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	
	
	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	_3sym3 Delta_LLL = sym3_3sym3_mul(gammaBar_LL, Delta_ULL);
	
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);

	real3x3 partial_LambdaBar_UL;
	{
		//partial_LambdaBar_ul[j].i := connBar^i_,j
<?=eqn:makePartial'LambdaBar_U'?>	
		//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
		partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);
	}

<?=eqn:makePartial2'epsilon_LL'?>

	sym3 trBar_partial2_gammaBar_LL = calc_trBar_partial2_gammaBar_LL(
		U, 
		x, 
		&gammaBar_UU, 
		partial_epsilon_LLl, 
		partial2_epsilon_LLll);
	
	sym3 RBar_LL = calc_RBar_LL(
		U,
		x,
		&gammaBar_LL,
		&gammaBar_UU,
		&connHat_ULL,
		&partial_gammaBar_LLL,
		&trBar_partial2_gammaBar_LL,
		&partial_LambdaBar_UL,
		&Delta_U,
		&Delta_ULL,
		&Delta_LLL);

	//RBar := RBar_ij gammaBar^ij
	real RBar = sym3_dot(gammaBar_UU, RBar_LL);

	//tr_DBar2_phi := gammaBar^ij DBar_i DBar_j phi = gammaBar^ij phi_,ij - connBar^k phi_,k
	real tr_DBar2_phi = sym3_dot(gammaBar_UU, partial2_phi_LL)
		- real3_dot(_3sym3_sym3_dot23(connBar_ULL, gammaBar_UU), partial_phi_L);

	//2017 Ruchlin et al, eqn 46
	//H = 2/3 K^2 - ABar^ij ABar_ij + exp(-4 phi) (RBar - 8 DBar^i phi DBar_i phi - 8 gammaBar^ij DBar_i DBar_j phi)
	U->H = 2. / 3. * U->K * U->K
		- sym3_dot(U->ABar_LL, ABar_UU)
		+ exp_neg4phi * (0.
			+ RBar
			- 8. * real3_weightedLenSq(partial_phi_L, gammaBar_UU) 
			- 8. * tr_DBar2_phi
		)
		- 16. * M_PI * U->rho;

#if 1

<?=eqn:makePartial'ABar_LL'?>		//partial_ABar_lll[k].ij = ABar_ij,k
	_3sym3 partial_ABar_LLL;
	calc_partial_ABar_LLL(&partial_ABar_LLL, U, x, partial_ABar_LLl);

	real3 partial_K_L;
<? for i,xi in ipairs(xNames) do
?>	partial_K_L.<?=xi?> = partial_K_l[<?=i-1?>] / calc_len_<?=xi?>(x);
<? end
?>

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
		- 8. * M_PI * U->S_u.<?=xi?> * calc_len_<?=xi?>(x)
<?	for j,xj in ipairs(xNames) do
?>		+ 6. * ABar_UU.<?=sym(i,j)?> * partial_phi_L.<?=xj?>
		- 2. / 3. * exp_6phi * gammaBar_UU.<?=sym(i,j)?> * partial_K_L.<?=xj?>
<?		for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xi?>.<?=sym(j,k)?> * ABar_UU.<?=sym(j,k)?>
<?			for l,xl in ipairs(xNames) do
?>		- ABar_UU.<?=sym(k,j)?> * gammaBar_UU.<?=sym(l,i)?> * partial_gammaBar_LLL.<?=xj?>.<?=sym(k,l)?>
		- ABar_UU.<?=sym(k,i)?> * gammaBar_UU.<?=sym(l,j)?> * partial_gammaBar_LLL.<?=xj?>.<?=sym(k,l)?>
		+ gammaBar_UU.<?=sym(k,i)?> * gammaBar_UU.<?=sym(l,j)?> * partial_ABar_LLL.<?=xj?>.<?=sym(k,l)?>
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
	
	const global cons_t* U = UBuf + index;
	global cons_t* deriv = derivBuf + index;

	if (solver->diffuseSigma != 0.) { 
		//Kreiss-Oligar dissipation
		//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
		//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
		for (int i = 0; i < numIntStates; ++i) {
<?=eqn:makePartial2('ptr[i]', 'real', 'partial2_Ui_ll', false)?>
			real lap = 0<?
for j,xj in ipairs(xNames) do
	local jj = from3x3to6(j,j)
?> + partial2_Ui_ll[<?=jj-1?>]<?
end ?>;
//TODO should that be a + or a -?
			deriv->ptr[i] += solver->diffuseSigma/16. * lap;
		}
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
end		-- getCommonCode
?>
