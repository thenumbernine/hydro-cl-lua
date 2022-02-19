<? 
-- symmath simplify is going slow ...
-- TODO restore these afterwards in case you want to run multiple solvers
local symmath = require "symmath"
symmath.op.div:pushRule'Prune/conjOfSqrtInDenom'
symmath.op.div:pushRule'Factor/polydiv'
symmath.op.pow.wildcardMatches = nil
symmath.matchMulUnknownSubstitution = false


-- integrates whatsoever.
local useCalcDeriv = true
local useCalcDeriv_alpha = false
local useCalcDeriv_W = false
local useCalcDeriv_K = false
local useCalcDeriv_epsilon_LL = false	-- with epsilon_IJ enabled, after 2 iterations, LambdaBar^I fails
local useCalcDeriv_ABar_LL = false
local useCalcDeriv_LambdaBar_U = false	-- FIXME
local useCalcDeriv_beta_U = false		-- FIXME (B_U specifically)
-- useScalarField components:
local useCalcDeriv_Phi = true
local useCalcDeriv_Psi = true
local useCalcDeriv_Pi = true
-- scalar field source terms?
local coupleTandPhi = true
-- PIRK:
local useLBetaWithPIRK = true

-- constrains det gammaBar_ij = det gammaHat_ij, ABar^i_i = 0, and calculates H and M^i ... if the associated flags are set
local useConstrainU = false

-- Does the scalar field source terms (TODO put them in calcDeriv, since for f.d. calcDeriv() and addSource() are really the same.)
local useAddSource = true

-- Does Kreiss-Oligar dissipation (in the calcDeriv function)
local useKreissOligarDissipation = false


local useSENRShiftAndCoDerivs = true -- insert senr's derivative code in to see how things work
local file = require "ext.file"

-- symmath
local Tensor = require "symmath.Tensor"
local Constant = require "symmath.Constant"
local frac = require "symmath.op.div"


local lenExprs = coord.request"coord_dx"
local gammaHat_ll = coord.request"coord_gHol_ll"

local det_gammaHat = coord.request"coord_det_gHol"
local partial_det_gammaHat_l = coord.request"coord_partial_det_gHol_l"
local partial2_det_gammaHat_ll = coord.request"coord_partial2_det_gHol_ll"
local partial_gammaHat_lll = coord.request"coord_partial_gHol_lll""_kij"():permute"_ijk"
local connHat_lll = coord.request"coord_connHol_lll"
local connHat_ull = coord.request"coord_connHol_ull"
local partial_connHat_ulll = connHat_ull"^i_jk,l"():permute"^i_jkl"

local delta_ul = Tensor("^i_j", function(i,j) return i==j and 1 or 0 end)

local e = coord.request"eToEHol"
local eInv = coord.request"eHolToE"
?>

//// MODULE_NAME: <?=eqn_common?>
// only here to appease eqn/bssnok-fd.lua whose calcDTCell module adds eqn_common, which exists in its other subclasses
//// MODULE_NAME: eqn.macros 

//// MODULE_NAME: <?=calc_partial_det_gammaHat_l?>
//// MODULE_DEPENDS: <?=coord_partial_det_gHol_l?>

#define <?=calc_partial_det_gammaHat_l?> coord_partial_det_gHol_l

//// MODULE_NAME: calc_partial_det_gammaHat_L
//// MODULE_DEPENDS: <?=calc_partial_det_gammaHat_l?> <?=rescaleFromCoord_rescaleToCoord?>

real3 calc_partial_det_gammaHat_L(real3 x) {
	real3 partial_det_gammaHat_l = <?=calc_partial_det_gammaHat_l?>(x);
	real3 partial_det_gammaHat_L = real3_rescaleFromCoord_l(partial_det_gammaHat_l, x);
	return partial_det_gammaHat_L;
}

//// MODULE_NAME: calc_partial2_det_gammaHat_ll
//// MODULE_DEPENDS: <?=coord_partial2_det_gHol_ll?>

#define calc_partial2_det_gammaHat_ll coord_partial2_det_gHol_ll

//// MODULE_NAME: calc_partial2_det_gammaHat_LL
//// MODULE_DEPENDS: calc_partial2_det_gammaHat_ll <?=rescaleFromCoord_rescaleToCoord?>

sym3 calc_partial2_det_gammaHat_LL(real3 x) {
	sym3 partial2_det_gammaHat_ll = calc_partial2_det_gammaHat_ll(x);
	sym3 partial2_det_gammaHat_LL = sym3_rescaleFromCoord_ll(partial2_det_gammaHat_ll, x);
	return partial2_det_gammaHat_LL;
}

//// MODULE_NAME: calc_len_#
//// MODULE_DEPENDS: <?=coord_dx_i?>

<? for i,xi in ipairs(xNames) do
?>#define calc_len_<?=xi?>	coord_dx<?=i-1?>
<? end ?>

//// MODULE_NAME: calc_partial*_len*

/*
e_i^I = delta_i^I f_i is a diagonal matrix with f_i indexed function.  for spherical, f_i = diag(1,r,r sin(theta))
e^i_I = delta^i_I f^i is the inverse, so f^i = 1/(f_i)
I'm trying to keep the differentiations to an absolute minimum in the bssnok-fd-num files
coord_dx_i(x) is the same as f_i 
*/
<?
local partial_len_ll = lenExprs"_i,j"():permute"_ij"
local partial2_len_lll = partial_len_ll"_ij,k"():factorDivision():permute"_ijk"

for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
-- TODO? make this based on (pt) param (just like coord_dx), and get eqn:compile to cooperate
?>#define calc_partial_len_<?=xi..xj?>(pt)	(<?=eqn:compile(partial_len_ll[i][j])?>)
<? 	end
end

-- this shows up often enough in the math, when trying to convert everything to normalized coordinates, but TODO
-- it looks like simplifying this everywhere isn't so numerically stable, but only doing it when necessary
-- starting to look more like my very first implementation of normalized coordinate code...
local partial_len_over_len_lll = Tensor("_ijk", function(i,j,k)
	return (partial_len_ll[i][j] / lenExprs[k])()
end)
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>#define calc_partial_len_over_len_<?=xi..xj..xk?>(pt)	(<?=eqn:compile(partial_len_over_len_lll[i][j][k])?>)
<? 		end
	end
end

for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>#define calc_partial2_len_<?=xi..xj..xk?>(pt)	(<?=eqn:compile(partial2_len_lll[i][j][k])?>)
<? 		end
	end
end


-- shorthand
local len_len_ll = (lenExprs"_i" * lenExprs"_j")()
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>#define calc_len_len_<?=xi..xj?>(pt)	(<?=eqn:compile(len_len_ll[i][j])?>)
<?	end
end

-- this is used for all rank-2 tensors: A_ij,k = (A_IJ e_i^I e_j^J)_,k = A_IJ,k (e_i^I e_j^J) + A_IJ (e_i^I e_j^J)_,k
local partial_len_len_lll = len_len_ll"_ij,k"():permute"_ijk"
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>#define calc_partial_len_len_<?=xi..xj..xk?>(pt)	(<?=eqn:compile(partial_len_len_lll[i][j][k])?>)
<? 		end
	end
end

-- e_deInv_eInv^K_LI = e_k^K e^k_L,i e^i_I
-- or maybe "eInv_deInv_norm" ?
local e_deInv_eInv = (e"_k^K" * eInv"^k_L,i" * eInv"^i_I")():permute"^K_LI"
for K,xK in ipairs(xNames) do
	for L,xL in ipairs(xNames) do
		for I,xI in ipairs(xNames) do
?>#define calc_e_deInv_eInv_<?=xK..xL..xI?>(pt)	(<?=eqn:compile(e_deInv_eInv[K][L][I])?>)
<?		end
	end
end

-- eSq_deSq_norm_IJK^MN = ((e_i^M e_j^N)_,k e^k_K e^i_I e^j_J)
local eSq_deSq_norm = (eInv"^i_I" * eInv"^j_J" * (e"_i^M" * e"_j^N")():permute"_ij^MN""_ij^MN_,k"() * eInv"^k_K")():permute"_IJK^MN"
for I,xI in ipairs(xNames) do
	for J,xJ in ipairs(xNames) do
		for K,xK in ipairs(xNames) do
			for M,xM in ipairs(xNames) do
				for N,xN in ipairs(xNames) do
?>#define calc_eSq_deSq_norm_<?=xI..xJ..xK..xM..xN?>(pt)	(<?=eqn:compile(eSq_deSq_norm[I][J][K][M][N])?>)
<?				end
			end
		end
	end
end
?>

//// MODULE_NAME: cplx3_add5
//// MODULE_DEPENDS: cplx3

#define cplx3_add5(a,b,c,d,e)	cplx3_add(cplx3_add(a,b),cplx3_add3(c,d,e))

//// MODULE_NAME: real3_add5

#define real3_add5(a,b,c,d,e)	real3_add(real3_add(a,b),real3_add3(c,d,e))

//// MODULE_NAME: real3_add6

#define real3_add6(a,b,c,d,e,f)	real3_add(real3_add3(a,b,c),real3_add3(d,e,f))
	
//// MODULE_NAME: sym3_add3

#define sym3_add3(a,b,c)	sym3_add(sym3_add(a,b),c)

//// MODULE_NAME: sym3_add4

#define sym3_add4(a,b,c,d)	sym3_add(sym3_add(a,b),sym3_add(c,d))	

//// MODULE_NAME: real3x3_partial_rescaleFromCoord_Ul
//// MODULE_DEPENDS: calc_len_# calc_partial*_len*

/*
derivative index of the result is last
e_i^I (T^M e^i_M)_,j e^j_J
= e_i^I (T^M_,j e^i_M + T^M e^i_M_,j) e^j_J
= T^I_,j e^j_J + T^M e_i^I e^i_M,j e^j_J
*/
real3x3 real3x3_partial_rescaleFromCoord_Ul(
	real3 const T_U, 
	real3x3 const partial_T_Ul,	// derivative in first index
	real3 const pt
) {
	real3x3 partial_T_UL;

<? 
local bleh = (e"_i^I" * eInv"^i_M,j" * eInv"^j_J")():permute"^I_MJ"

for J,xJ in ipairs(xNames) do
 	for I,xI in ipairs(xNames) do
?>	partial_T_UL.<?=xJ?>.<?=xI?> = 
		partial_T_Ul.<?=xJ?>.<?=xI?> / calc_len_<?=xJ?>(pt)
<?		for M,xM in ipairs(xNames) do
			if not Constant:isa(bleh[I][J][M], 0) then
?>		+ T_U.<?=xM?> * (<?=eqn:compile(bleh[I][J][M])?>)
<? 			end
		end
?>	;
<?	end
end
?>
	
	return partial_T_UL;
}

//// MODULE_NAME: real3x3_partial_rescaleToCoord_Ul
//// MODULE_DEPENDS: calc_len_# calc_partial*_len*

/*
derivative index of result is last 
T^i_,j = (T^I e^i_I)_,j
= T^I_,j e^i_I + T^I e^i_I,j
= (T^I_,j - T^I f_i,j / f_i) / f_i delta^i_I
*/
real3x3 real3x3_partial_rescaleToCoord_Ul(
	real3 T_U,
	real3x3 partial_T_Ul,	//derivative is first index
	real3 x
) {
	real3x3 partial_T_ul;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	partial_T_ul.<?=xi?>.<?=xj?> = (
		partial_T_Ul.<?=xj?>.<?=xi?>
		- T_U.<?=xi?> * calc_partial_len_<?=xi..xj?>(x) / calc_len_<?=xi?>(x)
	) / calc_len_<?=xi?>(x);
<? 	end
end
?>	return partial_T_ul;
}

//// MODULE_NAME: cplx3x3_partial_rescaleFromCoord_Ll
//// MODULE_DEPENDS: calc_len_# calc_partial*_len*

/*
This is breaking my old conventions too
derivative index of result is last 
T_I,j e^j_J = (T_i e^i_I)_,j e^j_J
= (T_i,j e^i_I + T_i e^i_I,j) e^j_J
= (T_i,j f^i + T_i f^i_,j) f^j delta^i_I delta^j_J
= (T_i,j - T_i f_i,j / f_i) / (f_i f_j) delta^i_I delta^j_J
*/
cplx3x3 cplx3x3_partial_rescaleFromCoord_Ll(
	cplx3 T_l,
	cplx3x3 partial_T_ll,	//derivative is first index
	real3 x
) {
	cplx3x3 partial_T_LL;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	partial_T_ll.<?=xi?>.<?=xj?> = cplx_real_mul(
		cplx_sub(
			partial_T_ll.<?=xj?>.<?=xi?>,
			cplx_real_mul(
				T_l.<?=xi?>,
				calc_partial_len_<?=xi..xj?>(x) / calc_len_<?=xi?>(x)
			)
		),
		1. / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
	);
<? 	end
end
?>	return partial_T_LL;
}

//// MODULE_NAME: real3x3_partial_rescaleFromCoord_Ll
//// MODULE_DEPENDS: calc_len_# calc_partial*_len*

/*
derivative index is last
e^i_I (T_M e_i^M)_,j e^j_J
= e^i_I (T_M,j e_i^M + T_M e_i^M_,j) e^j_J
= (T_I,j + T_I f^i f_i,j) e^j_J
= (T_I,j + T_I f_i,j / f_i) / f_j delta^j_J
*/
real3x3 real3x3_partial_rescaleFromCoord_Ll(real3 T_L, real3 const partial_T_Ll[3], real3 x) {
	real3x3 partial_T_LL;
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	partial_T_LL.<?=xi?>.<?=xj?> = (0.
		+ partial_T_Ll[<?=j-1?>].<?=xi?>
		+ T_L.<?=xi?> * calc_partial_len_<?=xi..xj?>(x) / calc_len_<?=xi?>(x) 
	) / calc_len_<?=xj?>(x);
<?	end
end ?>
	return partial_T_LL;
}

//// MODULE_NAME: calc_partial_connHat_Ulll_*

//calc_partial_connHat_Ulll_ijkl := e_i^I connHat^i_jk,l
<? 
local partial_connHat_Ulll = Tensor("^I_jkl", function(i,j,k,l)
	return (partial_connHat_ulll[i][j][k][l] * lenExprs[i])()
end)
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		for l,xl in ipairs(xNames) do
?>#define calc_partial_connHat_Ulll_<?=xi..xjk..xl?>(pt) (<?=eqn:compile(partial_connHat_Ulll[i][j][k][l])?>)
<?		end
	end
end
?>

//// MODULE_NAME: calc_partial*_det_gammaHat_over_det_gammaHat_*

<?
local partial_det_gammaHat_over_det_gammaHat_L = Tensor("_i", function(i)
	return (partial_det_gammaHat_l[i] / (det_gammaHat * lenExprs[i]))()
end)
?>
real3 calc_partial_det_gammaHat_over_det_gammaHat_L(real3 pt) {
	real3 partial_det_gammaHat_over_det_gammaHat_L;
<? 
for i,xi in ipairs(xNames) do
?>	partial_det_gammaHat_over_det_gammaHat_L.<?=xi?> = <?=eqn:compile(partial_det_gammaHat_over_det_gammaHat_L[i])?>;
<? end
?>	return partial_det_gammaHat_over_det_gammaHat_L;
}

<?
local partial2_det_gammaHat_over_det_gammaHat_LL = Tensor("_ij", function(i,j)
	return (partial2_det_gammaHat_ll[i][j] / (det_gammaHat * lenExprs[i] * lenExprs[j]))()
end)
?>
sym3 calc_partial2_det_gammaHat_over_det_gammaHat_LL(real3 pt) {
	sym3 partial2_det_gammaHat_over_det_gammaHat_LL;
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_det_gammaHat_over_det_gammaHat_LL.<?=xij?> = <?=eqn:compile(partial2_det_gammaHat_over_det_gammaHat_LL[i][j])?>;
<?
end
?>	return partial2_det_gammaHat_over_det_gammaHat_LL;
}

//// MODULE_NAME: calc_partial_ABar_LLL
//// MODULE_DEPENDS: calc_len_# calc_partial*_len*

static _3sym3 calc_partial_ABar_LLL(
	real3 const x,
	sym3 const ABar_LL,
	sym3 const partial_ABar_LLl[3]
) {
	_3sym3 partial_ABar_LLL;
	//derivative first
	//partial_ABar_lll.k.ij := ABar_ij,k
	//partial_ABar_LLL.K.IJ := ABar_IJ,k e_i^I e_j^J e_k^K + ABar_IJ (e_i^I e_j^J)_,k e_k^K
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_ABar_LLL.<?=xk?>.<?=xij?> = (0.
		+ partial_ABar_LLl[<?=k-1?>].<?=xij?>
		+ (
			calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			+ calc_partial_len_over_len_<?=xj..xk..xj?>(x)
		) * ABar_LL.<?=xij?>
	) / calc_len_<?=xk?>(x);
<?	end
end
?>	return partial_ABar_LLL;
}

//// MODULE_NAME: calc_partial_gammaBar_LLL
//// MODULE_DEPENDS: calc_partial*_len* calc_len_#

/*
partial_gammaBar_LLL.xK.xIJ = gammaBar_ij,k e^i_I e^j_J e^k_K
= (gammaBar_MN e_i^M e_j^N)_,k e^i_I e^j_J e^k_K
= (gammaBar_MN,k e_i^M e_j^N + gammaBar_MN (e_i^M e_j^N)_,k) e^i_I e^j_J e^k_K
= epsilon_IJ,k e^k_K + gammaBar_MN (e_i^M e_j^N)_,k e^i_I e^j_J e^k_K
*/
static _3sym3 calc_partial_gammaBar_LLL(
	real3 const pt,
	sym3 const epsilon_LL,
	sym3 const partial_epsilon_LLl[3]
) {
	_3sym3 partial_gammaBar_LLL;
	//derivative first
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	//partial_gammaBar_LLL.k.IJ := e_i^I e_j^J e_k^K gammaBar_ij,k
<? 
local bleh = ((e"_i^M" * e"_j^N")():permute"^MN_ij"
	"^MN_ij,k"():permute"^MN_ijk" 
	* eInv"^i_I" * eInv"^j_J" * eInv"^k_K"
)():permute"^MN_IJK"
for IJ,xIJ in ipairs(symNames) do
	local I,J,xI,xJ = from6to3x3(IJ)
	for K,xK in ipairs(xNames) do
?>	partial_gammaBar_LLL.<?=xK?>.<?=xIJ?> = 0.
		+ partial_epsilon_LLl[<?=K-1?>].<?=xIJ?> / calc_len_<?=xK?>(pt)
		+ (epsilon_LL.<?=xIJ?><?=I==J and " + 1." or ""?>) * (0.
<?		for M,xM in ipairs(xNames) do
			for N,xN in ipairs(xNames) do
				if not Constant.isValue(bleh[M][N][I][J][K], 0) then
?>			+ <?=eqn:compile(bleh[M][N][I][J][K])?>
<?				end
			end
		end
?>		);
<?	end
end
?>	return partial_gammaBar_LLL;
}

//// MODULE_NAME: calc_connHat_LLL_and_ULL
//// MODULE_DEPENDS: calc_partial*_len* calc_len_#

static void calc_connHat_LLL_and_ULL(
	_3sym3 * const connHat_LLL,
	_3sym3 * const connHat_ULL,
	global <?=cons_t?> const * const U,
	real3 const pt
) {
<?
local connHat_LLL = (connHat_lll"_ijk" * eInv"^i_I" * eInv"^j_J" * eInv"^k_K")():permute"_IJK"
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
for I,xI in ipairs(xNames) do
	for JK,xJK in ipairs(symNames) do
		local J,K,xJ,xK = from6to3x3(JK)
?>	connHat_LLL-><?=xI?>.<?=xJK?> = <?=eqn:compile(connHat_LLL[I][J][K])?>;
<?	end
end
?>	
	//connHat_ull[i].jk := connHat^i_jk = gammaHat^il connHat_ljk
	// raised by gammaHat_IJ, which equals delta^IJ
	*connHat_ULL = *connHat_LLL;
}

//// MODULE_NAME: calc_connBar_ULL

static _3sym3 calc_connBar_ULL(
	_3sym3 const partial_gammaBar_LLL,
	sym3 const gammaBar_UU
) {
	//connBar_lll.i.jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	//connBar_LLL.I.JK := e^i_I e^j_J e^k_K connBar_ijk
	_3sym3 connBar_LLL;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>	connBar_LLL.<?=xi?>.<?=xjk?> = .5 * (0.
		+ partial_gammaBar_LLL.<?=xk?>.<?=sym(i,j)?>
		+ partial_gammaBar_LLL.<?=xj?>.<?=sym(i,k)?>
		- partial_gammaBar_LLL.<?=xi?>.<?=sym(j,k)?>
	);
<?	end
end
?>	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
	_3sym3 connBar_ULL = sym3_3sym3_mul(gammaBar_UU, connBar_LLL);
	return connBar_ULL;
}

//// MODULE_NAME: calc_trBar_partial2_gammaBar_ll
//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?> from3x3to6 

static sym3 calc_trBar_partial2_gammaBar_ll(
	global <?=cons_t?> const * const U,
	real3 const x,
	sym3 const gammaBar_UU,
	sym3 const partial_epsilon_LLl[3],
	sym3 const partial2_epsilon_LLll[6]
) {
	//gives accuracy errors if I try to inline this
	sym3 const gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x); 
	
	/*
	gammaBar_ij,kl gammaBar^kl
	= (
		(delta_IJ + epsilon_IJ) (
			f_I,kl f_j
			+ 2 f_I,k f_J,l 
			+ f_J,kl f_i 
		)
		+ 2 epsilon_IJ,l (f_I,k f_j + f_J,k f_i)
		+ epsilon_IJ,kl f_i f_j
	) gammaBar^kl delta_i^I delta_j^J
	*/
	sym3 trBar_partial2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	trBar_partial2_gammaBar_ll.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
			local kl,xkl = from3x3to6(k,l)
?>		+ (
			+ (U->epsilon_LL.<?=xij?><?=i==j and " + 1." or ""?>) * (0.
				+ calc_partial2_len_<?=xi..xk..xl?>(x) * calc_len_<?=xj?>(x)
				+ calc_partial2_len_<?=xj..xk..xl?>(x) * calc_len_<?=xi?>(x)
				+ 2. * calc_partial_len_<?=xi..xk?>(x) * calc_partial_len_<?=xj..xl?>(x) 
			)
			+ 2. * partial_epsilon_LLl[<?=l-1?>].<?=xij?> * (0.
				+ calc_partial_len_<?=xi..xk?>(x) * calc_len_<?=xj?>(x)
				+ calc_partial_len_<?=xj..xk?>(x) * calc_len_<?=xi?>(x)
			)
			+ partial2_epsilon_LLll[<?=kl-1?>].<?=xij?> * calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x)
		) * gammaBar_uu.<?=xkl?>
<?		end
	end
?>	;
<? end
?>
	return trBar_partial2_gammaBar_ll;
}

//// MODULE_NAME: <?=applyInitCondCell?>
//// MODULE_DEPENDS: _3sym3 <?=coordMap?> <?=calc_gammaHat_ll?> <?=calc_det_gammaBar?> calc_connHat_LLL_and_ULL

// Should initCond provide a metric in cartesian, or in the background metric?
// I'll say Cartesian for now, and then transform them using the rescaling.

<?
-- look for symmath expressions instead of code
-- also skip the initDerivs finite difference 
if initCond.initAnalytical then
	local partial_gamma0_lll = initCond.gamma0_ll"_ij,k"():permute"_ijk"

	--local symmath = require "symmath"
	--local det_gamma = symmath.Matrix.determinant(initCond.gamma0_ll)
?>
//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?> <?=_3sym3_rescaleFromCoord__3sym3_rescaleToCoord?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	//compile compat. the compile code uses 'pt' as the grid point.
	real3 const pt = x;

	U->alpha = <?=eqn:compile(initCond.alpha0)?>;

	sym3 gamma_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = <?=eqn:compile(initCond.gamma0_ll[i][j])?>,
<? end
?>	};

	sym3 K_ll = (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = <?=eqn:compile(initCond.K0_ll[i][j])?>,
<? end
?>	};

	sym3 gammaHat_ll = <?=calc_gammaHat_ll?>(x);
	real det_gammaHat = sym3_det(gammaHat_ll);
	real det_gammaBar = det_gammaHat;
	real det_gamma = sym3_det(gamma_ll);
	real det_gammaBar_over_det_gamma = det_gammaBar / det_gamma;
	real exp_neg4phi = cbrt(det_gammaBar_over_det_gamma);
	U->W = sqrt(exp_neg4phi);
	sym3 gammaBar_ll = sym3_real_mul(gamma_ll, exp_neg4phi);
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);
	
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);
	U->K = sym3_dot(gamma_uu, K_ll);
	sym3 A_ll = sym3_sub(K_ll, sym3_real_mul(gamma_ll, -U->K/3.));
	sym3 ABar_ll = sym3_real_mul(A_ll, exp_neg4phi);
	U->ABar_LL = sym3_rescaleFromCoord_ll(ABar_ll, x);

	real3 beta_u = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = <?=eqn:compile(initCond.beta0_u[i])?>,
<? end
?>	};
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);
	U->B_U = real3_zero;

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	// gammaBar_ij = exp(-4 phi) gamma_ij
	// gammaBar^ij = exp(4 phi) gamma^ij
	real exp_4phi = 1. / exp_neg4phi;
	sym3 gammaBar_uu = sym3_real_mul(gamma_uu, exp_4phi);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);

	//partial_gamma_lll.k.ij := gamma_ij,k
	_3sym3 partial_gamma_lll;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_gamma_lll.<?=xk?>.<?=xij?> = <?=eqn:compile(partial_gamma0_lll[i][j][k])?>;
<?	end
end ?>

	// (det gamma)_,i = det gamma * gamma^jk gamma_jk,i
	real3 partial_det_gamma_l;
<? for i,xi in ipairs(xNames) do
?>	partial_det_gamma_l.<?=xi?> = det_gamma * (0.
<?	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>		+ gamma_uu.<?=sym(j,k)?> * partial_gamma_lll.<?=xi?>.<?=sym(j,k)?>
<?		end
	end
?>	);
<? end
?>

	real3 partial_det_gammaHat_l;
<? 
for i,xi in ipairs(xNames) do
?>	partial_det_gammaHat_l.<?=xi?> = <?=eqn:compile(partial_det_gammaHat_l[i])?>;
<? end
?>
	//W_,i = exp(2 phi)_,i
	// = ((det gammaBar / det gamma)^(1/6))_,i
	// = 1/6 (det gammaBar / det gamma)^(-5/6) (det gammaBar / det gamma)_,i
	// = 1/6 (det gammaBar / det gamma)^(1/6) / (det gammaBar / det gamma) [(det gammaBar)_,i det gamma - det gammaBar (det gamma)_,i] / (det gamma)^2
	// = 1/6 W / (det gammaBar / det gamma) [(det gammaBar)_,i - (det gamma)_,i (det gammaBar / det gamma)] / (det gamma)
	real3 partial_W_l;
<? for i,xi in ipairs(xNames) do
?>	partial_W_l.<?=xi?> = 1./6. * U->W / det_gammaBar_over_det_gamma * (
		partial_det_gammaHat_l.<?=xi?> - partial_det_gamma_l.<?=xi?> / det_gammaBar_over_det_gamma
	) / det_gamma;
<? end
?>

	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	//= ( W^-2 gamma_ij )_,k
	//= ( W^-2 gamma_ij )_,k
	//= -2 W^-3 W_,k gamma_ij + W^-2 gamma_ij,k
	_3sym3 partial_gammaBar_lll;
<? for ij,xij in ipairs(symNames) do
	for k,xk in ipairs(xNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = (
		-2. * partial_W_l.<?=xk?> / U->W * gamma_ll.<?=xij?> 
		+ partial_gamma_lll.<?=xk?>.<?=xij?>
	) / (U->W * U->W);
<?	end
end ?>
	_3sym3 partial_gammaBar_LLL = _3sym3_rescaleFromCoord_lll(partial_gammaBar_lll, x);

//// MODULE_DEPENDS: calc_connBar_ULL
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

	real3 Delta_U = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);
	real3 LambdaBar_U = real3_add(Delta_U, <?=mystery_C_U?>);

	U->LambdaBar_U = LambdaBar_U;


//TODO initialization of these ...
//how about an initial call to <?=constrainU?> ?	
	U->rho = 0.;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	
	U->H = 0.;
	U->M_U = real3_zero;
}

<? elseif initCond.useBSSNVars then -- and not initCond.initAnalytical ?>

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	<?=SETBOUNDS?>(0,0);
	real3 const x = cell->pos;

	if (<?=OOB?>(solver->numGhost, solver->numGhost)) {
		U->alpha = INFINITY;
		U->W = INFINITY;
		U->K = INFINITY;
		U->beta_U = _real3(INFINITY, INFINITY, INFINITY);
		U->B_U = _real3(INFINITY, INFINITY, INFINITY);
		U->LambdaBar_U = _real3(INFINITY, INFINITY, INFINITY);
		U->epsilon_LL = _sym3(INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY);
		U->ABar_LL = _sym3(INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY);
		U->H = INFINITY;
		U->M_U = _real3(INFINITY, INFINITY, INFINITY);
		U->rho = INFINITY;
		U->S_u = _real3(INFINITY, INFINITY, INFINITY);
		U->S_ll = _sym3(INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY);
		return;
	}


	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	//bssn vars - use these for init.senr:
	real alpha = 1.;
	real W = 1.;
	real K = 0.;
	real3 LambdaBar_U = real3_zero;
	real3 beta_U = real3_zero;
	real3 B_U = real3_zero;
	sym3 epsilon_LL = sym3_zero;
	sym3 ABar_LL = sym3_zero;

	// stress-energy tensor
	real rho = 0.;

<? if eqn.useScalarField then ?>	
	cplx Phi = cplx_zero;
	cplx3 Psi_l = cplx3_zero;
	cplx Pi = cplx_zero;
<? end ?>

	<?=initCode()?>

	//bssn vars - use these for init.senr:
	U->alpha = alpha;
	U->K = K;
	U->W = W;
	U->LambdaBar_U = LambdaBar_U;
	U->beta_U = beta_U;
	U->B_U = B_U;
	U->epsilon_LL = epsilon_LL;
	U->ABar_LL = ABar_LL;

	//stress-energy fields
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0.;
	U->M_U = real3_zero;

<? if eqn.useScalarField then ?>
	U->Phi = Phi;
	U->Psi_l = Psi_l;	//init with a numeric derivative?
	U->Pi = Pi;
<? end ?>
}

<?	else -- not initCond.initAnalytical and not initCond.useBSSNVars 
-- if we're using a SENR init cond then init the components directly
-- TODO port these from sympy into symmath 
?>
//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?> <?=calc_gammaBar_LL?> 

void <?=applyInitCondCell?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const xc = coordMap(x);
	real3 const mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);

	sym3 gammaHat_ll = <?=calc_gammaHat_ll?>(x);

	//initCond will assume it is providing a metric in Cartesian
	real alpha = 1.;
	real3 beta_u = real3_zero;
	real3 B_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;
	real rho = 0.;

<? if eqn.useScalarField then ?>	
	cplx Phi = cplx_zero;
	cplx3 Psi_l = cplx3_zero;
	cplx Pi = cplx_zero;
<? end ?>

	<?=initCode()?>

	//rescale from cartesian to spherical
	//TODO what about rotations for change of coordinates?
	// this is another reason why initial conditions should be symbolic
//	beta_u = real3_rescaleToCoord_U(beta_u, x);
	gamma_ll = sym3_rescaleToCoord_LL(gamma_ll, x);
//	K_ll = sym3_rescaleToCoord_LL(K_ll, x);

	U->alpha = alpha;
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);

	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);
	
	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar = <?=calc_det_gammaBar?>(x); 

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = cbrt(det_gammaBar / det_gamma);

	//W = exp(-2 phi)
	U->W = sqrt(exp_neg4phi);

	sym3 gammaBar_ll = sym3_real_mul(gamma_ll, exp_neg4phi);
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);

	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_real_mul(gamma_ll, 1./3. * U->K));
	sym3 ABar_ll = sym3_real_mul(A_ll, exp_neg4phi);
	U->ABar_LL = sym3_rescaleFromCoord_ll(ABar_ll, x);
	
	U->B_U = real3_rescaleFromCoord_u(B_u, x);

	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0.;
	U->M_U = real3_zero;

<? if eqn.useScalarField then ?>
	U->Phi = Phi;
	U->Psi_l = Psi_l;	//init with a numeric derivative?
	U->Pi = Pi;
<? end ?>

}

//// MODULE_NAME: <?=initDerivs?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS?>

//after popularing gammaBar_ll, use its finite-difference derivative to initialize LambdaBar_u
//TODO do this symbolically.  That's what I originally did, but symbolic calculations were getting complex
// however, with spherical BSSN, you need to 
kernel void <?=initDerivs?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const x = cellBuf[index].pos;
	global <?=cons_t?> * const U = UBuf + index;

	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);

//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

//// MODULE_DEPENDS: calc_partial_gammaBar_LLL 
<?=eqn:makePartial1"epsilon_LL"?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	//Delta^i_jk = connBar^i_jk - connHat^i_jk
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	
	U->LambdaBar_U = real3_add(_3sym3_sym3_dot23(Delta_ULL, gammaBar_UU), <?=mystery_C_U?>);
}

<? end -- initCond.initAnalytical or initCond.useBSSNVars ?>

//// MODULE_NAME: calc_RBar_LL
//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?> calc_partial_connHat_Ulll_* calc_len_# <?=_3sym3_rescaleFromCoord__3sym3_rescaleToCoord?>

static sym3 calc_RBar_LL(
	global <?=cons_t?> const * const U,
	real3 const x,
	sym3 const * const gammaBar_LL,
	sym3 const * const gammaBar_UU,
	_3sym3 const * const connHat_ULL,
	_3sym3 const * const partial_gammaBar_LLL,
	sym3 const * const trBar_partial2_gammaBar_ll,
	real3x3 const * const partial_LambdaBar_UL,
	real3 const * const Delta_U,
	_3sym3 const * const Delta_ULL,
	_3sym3 const * const Delta_LLL
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
	sym3 trBar_DHat2_gammaBar_without_partial2_gammaBar_ll;
<?
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	trBar_DHat2_gammaBar_without_partial2_gammaBar_ll.<?=xij?> = (0.
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
?>	);
<?
end
?>

	//trBar_DHat2_gammaBar_ll.ij := gammaBar^kl DHat_k DHat_l gammaBar_ij
	sym3 trBar_DHat2_gammaBar_ll = sym3_add(
		*trBar_partial2_gammaBar_ll,
		trBar_DHat2_gammaBar_without_partial2_gammaBar_ll);

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
			- .5 * trBar_DHat2_gammaBar_ll.<?=xij?> / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
<?	for k,xk in ipairs(xNames) do
?>			//Less accurate when I defer this
			+ .5 * gammaBar_LL-><?=sym(i,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_LL-><?=sym(j,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xi?>
			
			//Less accurate when deferring these calculations
			+ .5 * Delta_U-><?=xk?> * (0.
				+ Delta_LLL-><?=xi?>.<?=sym(k,j)?> 
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

//// MODULE_NAME: applyKreissOligar
//// MODULE_DEPENDS: <?=coordMapR?> eqn.macros

//////////////////////////////// Kreiss-Oligar dissipation //////////////////////////////// 

static void applyKreissOligar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell,
	global <?=cons_t?> * const deriv,
	real3 const x,
	int const * const fields,
	int const numFields
) {
<? if useKreissOligarDissipation then ?>
	if (solver->dissipationCoeff == 0.) return;
		
	//Kreiss-Oligar dissipation
#if 0	// 2013 Baumgarte et al, section IIIB
	real const coeff = solver->dissipationCoeff / dt;
#elif 1	// 2017 Ruchlin et al eqn 44
	real const r = coordMapR(x);
	// all params > 0
	real const rKO = 2.;		//rKO/M = 2
	real const wKO = .17;		//wKO/M = .17
	//epsKO0 = .99 doubles as my 'dissipationCoeff' var
	real const coeff = .5 * (erf((r - rKO) / wKO) + 1.) * solver->dissipationCoeff;
#endif

	real3 dy = solver->grid_dx;
<? if require "hydro.coord.sphere_sinh_radial":isa(coord) then ?>
	real3 const yR = _real3(cell->r, cell->pos.y, cell->pos.z);
<? for i=1,solver.dim do
	local xi = xNames[i]
?>{
	global <?=cell_t?> const * const cellL = cell - solver->stepsize.<?=xi?>;
	real3 const yL = _real3(cellL->r, cellL->pos.y, cellL->pos.z);
	dy.<?=xi?> = real3_len(real3_sub(yR, yL));
}<? end ?>
<? end ?>
	
	real3 const _1_dy = _real3(1./dy.x, 1./dy.y, 1./dy.z);

	//described in 2008 Babiuc et al as Q = (-1)^r h^(2r-1) (D+)^r rho (D-)^r / 2^(2r)
	//...for r=2... -sigma h^3 (D+)^2 rho (D-)^2 / 16 ... and rho=1, except rho=0 at borders maybe.
	int const * const endOfFields = fields + numFields;
	for (int const * ip = fields; ip < endOfFields; ++ip) {
		int const i = *ip;
		
		deriv->ptr[i] += coeff * (0.
<? for j=1,solver.dim do
	local xj = xNames[j]
?>			+ (
				-20. * U->ptr[i]
				+ 15. * (U[1 * solver->stepsize.<?=xj?>].ptr[i] + U[-1 * solver->stepsize.<?=xj?>].ptr[i])
				+ -6. * (U[2 * solver->stepsize.<?=xj?>].ptr[i] + U[-2 * solver->stepsize.<?=xj?>].ptr[i])
				+ U[3 * solver->stepsize.<?=xj?>].ptr[i] + U[-3 * solver->stepsize.<?=xj?>].ptr[i]
			) * _1_dy.<?=xj?>
<? end
?>		) * (1. / 64.);
	}
<? end	-- useKreissOligarDissipation ?>
}

//// MODULE_NAME: tracefree

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
	return sym3_sub(A_ll, sym3_real_mul(g_ll, tr_A * (1. / 3.)));
}

//// MODULE_NAME: getUpwind

/*
Returns the step coefficients, [+-1, +-1, +-1]
based on vector 'v'
to compute upwind differencing from.
Same thing as (int4)sgn(v)
*/
int4 getUpwind(real3 v) {
	return (int4)(
		v.x >= 0 ? 1 : -1,
		v.y >= 0 ? 1 : -1,
		v.z >= 0 ? 1 : -1,
		0
	);
}

//// MODULE_NAME: from3x3to6

//ok I have this in lua code for inline codegen, but here it is in C code, while debugging:
int from3x3to6(int i, int j) {
	int minij = min(i,j);
	int maxij = max(i,j);
	return maxij + minij + (minij > 1);
}

//// MODULE_NAME: sym3_Lbeta_LL
//// MODULE_DEPENDS: from3x3to6

static sym3 sym3_Lbeta_LL(
	sym3 T_LL,					//T_LL.ij = T_ij
	_3sym3 partial_T_LLL,		//partial_T_LLL.k.ij = T_ij,k
	real3 beta_U,				//beta_U.i = beta^i
	real3x3 partial_beta_UL		//partial_beta_UL.i.j = beta^i_,j
) {
#if 0	
	sym3 beta_times_partial_T_LL = real3_3sym3_dot1(beta_U, partial_T_LLL);
	real3x3 T_times_partial_beta_LL = sym3_real3x3_mul(T_LL, partial_beta_UL);
	sym3 sym_T_times_partial_beta_LL = sym3_real_mul(sym3_from_real3x3(T_times_partial_beta_LL), 2.);
	sym3 Lbeta_T_LL = sym3_add(beta_times_partial_T_LL, sym_T_times_partial_beta_LL);
#else
	
	sym3 Lbeta_T_LL;
	int ij = 0;
	for (int i = 0; i < 3; ++i) {
		for (int j = i; j < 3; ++j, ++ij) {
			Lbeta_T_LL.s[ij] = 0;
			for (int k = 0; k < 3; ++k) {
				//  beta^k ∂up_k(gammaBar_ij)
				//+ gammaBar_kj beta^k_,i
				//+ gammaBar_ki beta^k_,j
				Lbeta_T_LL.s[ij] += partial_T_LLL.v[k].s[ij] * beta_U.s[k]
					+ T_LL.s[from3x3to6(k,j)] * partial_beta_UL.v[k].s[i]
					+ T_LL.s[from3x3to6(k,i)] * partial_beta_UL.v[k].s[j];
			}
		}
	}
#endif
	return Lbeta_T_LL;
}

//// MODULE_NAME: calcDeriv_epsilon_LL
//// MODULE_DEPENDS: sym3_add4 sym3_Lbeta_LL 

//////////////////////////////// epsilon_IJ_,t //////////////////////////////// 

/*
2017 Ruchlin et al, eqn 11a
epsilon_ij,t = 2/3 gammaBar_ij (alpha ABar^k_k - DBar_k beta^k) + DHat_i beta_j + DHat_j beta_i - 2 alpha ABar_ij + epsilon_ij,k beta^k + epsilon_ik beta^k_,j + epsilon_kj beta^k_,i
...using DBar_(i beta_j) = DHat_(i beta_j) + epsilon_k(j beta^k_,i) + 1/2 epsilon_ij,k beta^k
	
	epsilon_ij,t = 	
L1:		+ 2/3 gammaBar_ij (
			alpha ABar^k_k 
			- DBar_k beta^k
		) 
		- 2 alpha ABar_ij 

Lbeta:	+ beta^k ∂up_k(gammaBar_ij)
		+ gammaBar_kj beta^k_,i
		+ gammaBar_ki beta^k_,j

gammaBar_ij,k = connBar_ikj + connBar_jki


	epsilon_IJ,t = 	
L1:		+ 2/3 gammaBar_IJ (
			alpha ABar_KL gammaBar^KL
			- beta^K_,k e^k_K							----+
			- beta^K e^k_K,k								| -DBar_k beta^k
			- beta^K 1/2 e^k_K gammaBar_,l / gammaHat	----+
		) 
		- 2 alpha ABar_IJ 

Lbeta:	+ beta^K ∂up_k(gammaBar_ij) (e^k_K e^i_I e^j_J)	<- most of these will have some r/r's that will cancel)
		+ gammaBar_KJ beta^K_,i e^i_I
		+ gammaBar_KI beta^K_,j e^j_J
		+ gammaBar_KJ beta^L e^k_L_,i e^i_I e_k^K
		+ gammaBar_KI beta^L e^k_L_,j e^j_J e_k^K

senr uses:
∂up_k(gammaBar_ij) = 
	gammaHat_ij,k					<- not upwind
	+ e_i^I e_j^J dup_k epsilon_IJ 	<- no upwind derivative of normalization transforms
	+ epsilon_IJ (e_i^I e_j^J)_,k 	<- all analytically symplified

so dup in normalized coordinates:
beta^K ∂up_k(gammaBar_ij) (e^k_K e^i_I e^j_J) = 
	beta^K [(gammaHat_ij,k e^i_I e^j_J e^k_K)]	<- simplify me plz
	+ beta^K (dup_k epsilon_IJ) e^k_K
	+ beta^K gammaBar_MN (e_i^M e_j^N)_,k (e^k_K e^i_I e^j_J)

and the whole of Lbeta in normalized coordiantes:
	+ beta^K (gammaHat_ij,k e^i_I e^j_J e^k_K)				----\ beta^K e^k_K dup_k (gammaBar_IJ) e_i^I e_j^J
	+ beta^K (dup_k epsilon_IJ) e^k_K						----/
	+ beta^K gammaBar_MN ((e_i^M e_j^N)_,k e^k_K e^i_I e^j_J)
	+ gammaBar_KJ beta^K_,i e^i_I
	+ gammaBar_KI beta^K_,j e^j_J
	+ gammaBar_KJ beta^L (e^k_L,i e^i_I e_k^K)
	+ gammaBar_KI beta^L (e^k_L,j e^j_J e_k^K)
*/
static void calcDeriv_epsilon_LL(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const pt,
	sym3 const * const gammaBar_LL,
	real3x3 const * const partial_beta_Ul,
	real3x3 const * const ABar_UL,
	real const tr_DBar_beta
) {
	real const tr_ABar = real3x3_trace(*ABar_UL);

<?=eqn:makePartialUpwind"epsilon_LL"?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	_3sym3 const partial_gammaBar_LLL_upwind = calc_partial_gammaBar_LLL(
		pt,	//xup?
		U->epsilon_LL,	//upwind avg?
		partial_epsilon_LLl_upwind);

#if 0
	sym3 Lbeta_gammaBar_LL = sym3_Lbeta_LL(
		*gammaBar_LL,
		partial_gammaBar_LLL_upwind,
		U->beta_U,
		*partial_beta_UL);

#else	// doesn't seem to be any better
<?
local partial_gammaHat_LLL = (partial_gammaHat_lll"_ijk" * eInv"^i_I" * eInv"^j_J" * eInv"^k_K")():permute"_IJK"
?>
	_3sym3 partial_gammaHat_LLL;	//TODO use this in the other place
<? 	
for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
		for k,xk in ipairs(xNames) do
?>	partial_gammaHat_LLL.<?=xk?>.<?=xij?> = <?=eqn:compile(partial_gammaHat_LLL[i][j][k])?>;
<?		end
	end
?>

	sym3 Lbeta_gammaBar_LL;
<? 	for IJ,xIJ in ipairs(symNames) do 
		local I,J,xI,xJ = from6to3x3(IJ)
?>	Lbeta_gammaBar_LL.<?=xIJ?> = 0.
<? 		for K,xK in ipairs(xNames) do
?>		+ U->beta_U.<?=xK?> * (0.
			+ partial_gammaHat_LLL.<?=xK?>.<?=xIJ?>												//+ beta^K [(gammaHat_ij,k e^i_I e^j_J e^k_K)]
			+ partial_epsilon_LLl_upwind[<?=K-1?>].<?=xIJ?> / calc_len_<?=xK?>(pt)				//+ beta^K (dup_k epsilon_IJ) e^k_K
<?			for M,xM in ipairs(xNames) do
				for N,xN in ipairs(xNames) do
					if not Constant.isValue(eSq_deSq_norm[I][J][K][M][N], 0) then
?>			+ gammaBar_LL-><?=sym(M,N)?> * calc_eSq_deSq_norm_<?=xI..xJ..xK..xM..xN?>(pt)	//+ beta^K gammaBar_MN (e_i^M e_j^N)_,k (e^k_K e^i_I e^j_J)
<?					end
				end
			end
?>		)
		+ gammaBar_LL-><?=sym(K,J)?> * partial_beta_Ul-><?=xI?>.<?=xK?> / calc_len_<?=xI?>(pt)	//+ gammaBar_KJ beta^K_,i e^i_I
		+ gammaBar_LL-><?=sym(K,I)?> * partial_beta_Ul-><?=xJ?>.<?=xK?> / calc_len_<?=xJ?>(pt)	//+ gammaBar_KI beta^K_,j e^j_J
<?			for L,xL in ipairs(xNames) do 
				if not Constant.isValue(e_deInv_eInv[K][L][I], 0) then
?>		+ gammaBar_LL-><?=sym(K,J)?> * U->beta_U.<?=xL?> * calc_e_deInv_eInv_<?=xK..xL..xI?>(pt)	//+ gammaBar_KJ beta^L (e_k^K e^k_L_,i e^i_I)
<?				end
				if not Constant.isValue(e_deInv_eInv[K][L][J], 0) then
?>		+ gammaBar_LL-><?=sym(K,I)?> * U->beta_U.<?=xL?> * calc_e_deInv_eInv_<?=xK..xL..xJ?>(pt)	//+ gammaBar_KI beta^L (e_k^K e^k_L_,j e^j_J)
<?				end
			end 
		end
?>	;
<? 	end 
?>
#endif

	deriv->epsilon_LL = sym3_add4(
		deriv->epsilon_LL,
		
		//part 2
		sym3_real_mul(*gammaBar_LL, (2. / 3.) * (U->alpha * tr_ABar - tr_DBar_beta)),
	
		//part 3
		sym3_real_mul(U->ABar_LL, -2. * U->alpha),
		
		//Lie derivative
		Lbeta_gammaBar_LL
	);
}



//// MODULE_NAME: calcDeriv_W
//// MODULE_DEPENDS: <?=rescaleFromCoord_rescaleToCoord?>

//////////////////////////////// W_,t //////////////////////////////// 

static void calcDeriv_W(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const x,
	real const tr_DBar_beta
) {
<?=eqn:makePartialUpwind"W"?>
	real3 partial_W_L_upwind = real3_rescaleFromCoord_l(partial_W_l_upwind, x);

	real Lbeta_W = real3_dot(partial_W_L_upwind, U->beta_U);

	real L2_W = (1. / 3.) * U->W * (U->alpha * U->K - tr_DBar_beta);

	/*
	2017 Ruchlin et al eqn 11c
	W_,t =
	L1:		1/3 W (alpha K - beta^k_,k - beta^k connBar^j_kj) 
	Lbeta:	+ beta^k W_,k
	
	connBar^j_kj = log(sqrt(gammaBar))_,k = gammaBar_,k / (2 gammaBar)
	gammaBar = gammaHat
so	
	W_,t =
	L1:		1/3 W (alpha K - beta^I_,i e^i_I - beta^I (1/2 e^i_I gammaHat_,i / gammaHat - e^i_I_,i)) 
	Lbeta:	+ beta^I e^i_I W_,i

	*/
	deriv->W += L2_W + Lbeta_W;
}

//// MODULE_NAME: calc_PIRK_L2_ABar_LL
//// MODULE_DEPENDS: sym3_add4 <?=rescaleFromCoord_rescaleToCoord?> real3x3_partial_rescaleFromCoord_Ul calc_RBar_LL tracefree sym3_add3 calc_trBar_partial2_gammaBar_ll

//////////////////////////////// A_IJ,t //////////////////////////////// 

/*
2017 Ruchlin et al eqn 11b
traceless portion of ...
	- 2 alpha DBar_i DBar_j phi 
	+ 4 alpha DBar_i phi DBar_j phi 
	+ 2 DBar_i phi DBar_j alpha 
	+ 2 DBar_i alpha DBar_j phi 
	- DBar_i DBar_j alpha 
	+ alpha RBar_ij 
	- 8 pi alpha S_ij
*/
static sym3 calc_PIRK_L2_ABar_LL(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const x,
	sym3 const * const gammaBar_LL,
	sym3 const * const gammaBar_UU,
	sym3 const * const DBar2_alpha_LL,
	real3 const * const partial_alpha_L, 
	real3 const * const partial_W_L,
	_3sym3 const * const connBar_ULL,
	//for RBar_IJ:
	sym3 partial_epsilon_LLl[3],
	_3sym3 const * const Delta_ULL,
	real3 const * const Delta_U,
	_3sym3 const * const connHat_ULL,
	_3sym3 const * const partial_gammaBar_LLL,
	sym3 const * const gammaBar_ll,
	sym3 const * const gammaBar_uu
) {
<?=eqn:makePartial2"epsilon_LL"?>

	sym3 trBar_partial2_gammaBar_ll = calc_trBar_partial2_gammaBar_ll(
		U, 
		x, 
		*gammaBar_UU, 
		partial_epsilon_LLl, 
		partial2_epsilon_LLll);

	_3sym3 Delta_LLL = sym3_3sym3_mul(*gammaBar_LL, *Delta_ULL);

	real3x3 partial_LambdaBar_UL;
	{
		//partial_LambdaBar_ul[j].i := LambdaBar^i_,j
<?=eqn:makePartial1"LambdaBar_U"?>
		//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
		partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);
	}

<? if true then -- gets to the 4th step of 2nd iter ?>
	
	sym3 RBar_LL = calc_RBar_LL(
		U,
		x,
		gammaBar_LL,
		gammaBar_UU,
		connHat_ULL,
		partial_gammaBar_LLL,
		&trBar_partial2_gammaBar_ll,
		&partial_LambdaBar_UL,
		Delta_U,
		Delta_ULL,
		&Delta_LLL);

	sym3 TF_RBar_LL = tracefree(RBar_LL, *gammaBar_LL, *gammaBar_UU);

<? elseif false then -- when immediately rescaled to RBar_IJ, fails immediately ?>

//// MODULE_DEPENDS: SENR_calc_RBar_ll
	sym3 RBar_LL = SENR_calc_RBar_ll(solver, U, x);
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	RBar_LL.<?=xij?> /= calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x);
<? end
?>
	sym3 TF_RBar_LL = tracefree(RBar_LL, *gammaBar_LL, *gammaBar_UU);

<? else -- rescale after trace-free, still fails immediately ?>

//// MODULE_DEPENDS: SENR_calc_RBar_ll
	sym3 RBar_ll = SENR_calc_RBar_ll(solver, U, x);
	sym3 TF_RBar_ll = tracefree(RBar_ll, *gammaBar_ll, *gammaBar_uu);
	sym3 TF_RBar_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	TF_RBar_LL.<?=xij?> = TF_RBar_ll.<?=xij?> / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x));
<? end
?>

<? end ?>


<?=eqn:makePartial2"W"?>			//partial2_W_ll.ij := W_,ij
	//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
	sym3 partial2_phi_LL_times_WSq;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	partial2_phi_LL_times_WSq.<?=xij?> = .5 * (
		- partial2_W_ll.<?=xij?> * U->W / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
		+ partial_W_L-><?=xi?> * partial_W_L-><?=xj?>
	);
<? end ?>
	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_LL_times_WSq = sym3_sub(
		partial2_phi_LL_times_WSq,
		sym3_real_mul(
			real3_3sym3_dot1(
				*partial_W_L,
				*connBar_ULL
			),
			-1/2 * U->W
		)
	);

	sym3 TF_DBar2_alpha_LL = tracefree(*DBar2_alpha_LL, *gammaBar_LL, *gammaBar_UU);

	real WSq = U->W * U->W;

	sym3 tracelessPart_LL_times_WSq = tracefree(
		sym3_add4(
			sym3_real_mul(
				sym3_from_real3x3(
					real3_real3_outer(*partial_W_L, *partial_alpha_L)
				),
				-2. * U->W
			),
			sym3_real_mul(
				DBar2_phi_LL_times_WSq,
				-2. * U->alpha
			),
			sym3_real_mul(
				real3_outer(*partial_W_L),
				-2. * U->alpha * U->W
			),
			sym3_real_mul(
				sym3_rescaleFromCoord_ll(U->S_ll, x),
				-8. * M_PI * U->alpha * WSq
			)
		), 
		*gammaBar_LL, 
		*gammaBar_UU);

	return sym3_add3(
		tracelessPart_LL_times_WSq,
		sym3_real_mul(TF_DBar2_alpha_LL, -WSq), 
		sym3_real_mul(TF_RBar_LL, U->alpha * WSq));
}


//// MODULE_NAME: calc_PIRK_L2_K

//////////////////////////////// K_,t //////////////////////////////// 

static inline real calc_PIRK_L2_K(
	global <?=cons_t?> const * const U,
	sym3 const * const gammaBar_uu,
	sym3 const * const DBar2_alpha_ll,
	real3 const * const partial_alpha_l,
	real3 const * const partial_W_l
) {

	//tr_DBar2_alpha := gammaBar^ij DBar_i DBar_j alpha
	real const tr_DBar2_alpha = sym3_dot(*gammaBar_uu, *DBar2_alpha_ll);
	
	//2013 Baumgarte eqn B3: 
	// SENR/NRPy's code has a note that the "+ alpha ABar_ij ABar^ij" term should be in L3
	// L2 K = -exp(-4phi) (DBar^2 alpha + 2 DBar^i alpha DBar_i phi) 
	// 		= -W^2 gammaBar^ij (alpha_,ij - connBar^k_ij alpha_,k + 2 alpha_,i phi_,j) 
	// 		= -W gammaBar^ij (W alpha_,ij - W connBar^k_ij alpha_,k - alpha_,i W_,j) 
	return -U->W * (
		  tr_DBar2_alpha * U->W
		- real3_weightedDot(*partial_W_l, *partial_alpha_l, *gammaBar_uu)
	);
}
// W = exp(-2 phi)
// W_,i = exp(-2 phi) * -2 phi_,i
// -1/2 W_,i = W phi_,i

//// MODULE_NAME: calc_PIRK_L3_ABar_LL

// NOTICE - all my PIRK_L3 functions DO NOT contribute Lbeta
/*
L3:		- 2/3 ABar_IJ DBar_k beta^k
		+ alpha ABar_IJ K
		- 2 alpha ABar_IK ABar_LJ gammaBar^kl e_k^K e_l^L
*/
static sym3 calc_PIRK_L3_ABar_LL(
	global <?=cons_t?> const * const U,
	real const tr_DBar_beta,
	sym3 const * const gammaBar_uu,
	real3 const pt
) {
<? if false then ?>	
	
	sym3 const ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, ABar_UL);
	sym3 const L3_ABar_LL = sym3_add(	
		sym3_real_mul(U->ABar_LL, U->alpha * U->K - (2. / 3.) * tr_DBar_beta),
		sym3_real_mul(ABarSq_LL, -2. * U->alpha)
	);	

<? else ?>
	
	sym3 L3_ABar_LL;
<? for IJ,xIJ in ipairs(symNames) do
	local I,J,xI,xJ = from6to3x3(IJ)
?>	L3_ABar_LL.<?=xIJ?> = 0.
		+ U->ABar_LL.<?=xIJ?> * (0.
			- (2./3.) * tr_DBar_beta
			+ U->alpha * U->K
		)
		- 2. * U->alpha * (0.
<?	for K,xK in ipairs(xNames) do
		for L,xL in ipairs(xNames) do
?>			+ U->ABar_LL.<?=sym(I,K)?> * U->ABar_LL.<?=sym(J,L)?> * gammaBar_uu-><?=sym(K,L)?> * calc_len_<?=xK?>(pt) * calc_len_<?=xL?>(pt)
<?		end
	end
?>		)
	;
<? end
?>

<? end ?>
	
	return L3_ABar_LL;
}

//// MODULE_NAME: calc_PIRK_L3_K

// NOTICE - all my PIRK_L3 functions DO NOT contribute Lbeta
static real calc_PIRK_L3_K(
	global <?=cons_t?> const * const U,
	sym3 const * const ABar_UU
) {
	//tr_ABarSq := ABar_ij ABar^ij = ABar_ij ABar_kl gammaBar^ik gammaBar^jl
	real const tr_ABarSq = sym3_dot(U->ABar_LL, *ABar_UU);

	// SENR/NRPy's code has a note that the "+ alpha ABar_ij ABar^ij" term should be in L3
	return U->alpha * (
		U->K * U->K * (1. / 3.)
		+ tr_ABarSq
	);
}

//// MODULE_NAME: calc_PIRK_L2_LambdaBar_U

//////////////////////////////// LambdaBar^I_,t //////////////////////////////// 

static inline real3 calc_PIRK_L2_LambdaBar_U(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	real3 const pt,
	sym3 const * const gammaBar_uu,
	real3x3 const * const partial_beta_Ul,
	_3sym3 const * const Delta_ULL,
	real3 const * const partial_alpha_L,
	real3 const * const partial_W_L,
	real3 const * const partial_K_l
) {
<? if false then ?>
#if 0	//old way
	real3 const x = pt;

	//partial2_beta_ULL.I.JK = e^i_I (beta^M e^i_M)_,jk e^j_J e^k_K
<?=eqn:makePartial2"beta_U"?>		
#if 1 //removing this for debugging and setting to zero.  no dif.
	_3sym3 partial2_beta_ULL;
<? for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>	partial2_beta_ULL.<?=xi?>.<?=xjk?> = (
		partial2_beta_Ull[<?=jk-1?>].<?=xi?> + (
			- partial_beta_Ul-><?=xj?>.<?=xi?> * calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			- partial_beta_Ul-><?=xk?>.<?=xi?> * calc_partial_len_over_len_<?=xi..xj..xj?>(x)
			+ U->beta_U.<?=xi?> * (
				calc_partial2_len_<?=xi..xj..xk?>(x) / calc_len_<?=xi?>(x)
				+ 2. * calc_partial_len_over_len_<?=xi..xj..xi?>(x) * calc_partial_len_over_len_<?=xi..xk..xi?>(x)
			)
		)
	) / (calc_len_<?=xj?>(x) * calc_len_<?=xk?>(x));
<?	end
end ?>

//// MODULE_DEPENDS: real3x3x3 
//// MODULE_DEPENDS: calc_partial_connHat_Ulll_* 
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
		+ connHat_ULL-><?=xi?>.<?=sym(l,j)?> * partial_beta_UL-><?=xl?>.<?=xk?>
		+ connHat_ULL-><?=xi?>.<?=sym(l,k)?> * partial_beta_UL-><?=xl?>.<?=xj?>
		- connHat_ULL-><?=xl?>.<?=sym(j,k)?> * partial_beta_UL-><?=xi?>.<?=xl?>
<?			for m,xm in ipairs(xNames) do
?>		+ connHat_ULL-><?=xi?>.<?=sym(m,k)?> * connHat_ULL-><?=xm?>.<?=sym(l,j)?> * U->beta_U.<?=xl?>
		- connHat_ULL-><?=xi?>.<?=sym(m,l)?> * connHat_ULL-><?=xm?>.<?=sym(k,j)?> * U->beta_U.<?=xl?>
<?			end
		end
?>	;
<?		end
	end 
end
?>
#else
	_3sym3 partial2_beta_ULL = _3sym3_zero;
	real3x3x3 DHat2_beta_ULL;
	for (int j = 0; j < 3; ++j) {
		for (int k = 0; k < 3; ++k) {
			for (int l = 0; l < 3; ++l) {
				DHat2_beta_ULL.v[l].v[k].s[j] = 0;
			}
		}
	}
#endif

	/*
	tr12_partial2_beta_l.i := beta^j_,ji
	beta^i_,jk is a sym3 of 3's ... so I don't have that struct yet ... 
	 What name could I use? sym3x3?  how about real3s3x3 where we have 's' for symmetric and 'x' for cartesian product.
	 Then sym3 turns into real3s3 and _3sym3 turns into real3x3s3.
	TODO simplify math plz
	*/
	real3 tr12_partial2_beta_L = _3sym3_tr12(partial2_beta_ULL);

	//SENR's formulas
	//detg is the ratio of det(gammaBar_ij)/det(gammaHat_ij)

//// MODULE_DEPENDS: <?=calc_det_gammaHat?>
	real det_gammaHat = <?=calc_det_gammaHat?>(x);
//// MODULE_DEPENDS: calc_partial_det_gammaHat_L
	real3 partial_det_gammaHat_L = calc_partial_det_gammaHat_L(x);
//// MODULE_DEPENDS: calc_partial2_det_gammaHat_LL
	sym3 partial2_det_gammaHat_LL = calc_partial2_det_gammaHat_LL(x);
	
	real det_gammaBar = det_gammaHat * detg;

	real3 partial_det_gammaBar_L = real3_add(
		real3_real_mul(partial_det_gammaHat_L, detg),
		real3_real_mul(*partial_detg_L, det_gammaHat));

//// MODULE_DEPENDS: sym3_add3
	sym3 partial2_det_gammaBar_LL = sym3_add3(
		sym3_real_mul(partial2_det_gammaHat_LL, detg),
		sym3_real_mul(sym3_from_real3x3(real3_real3_outer(partial_det_gammaHat_L, *partial_detg_L)), 2.),
		sym3_real_mul(*partial2_detg_LL, det_gammaHat));
	
	real3 partial_det_gammaBar_times_partial_beta_L = real3_real3x3_mul(partial_det_gammaBar_L, *partial_beta_UL);

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
			gammaBar_,ij beta^j
			- gammaBar_,i gammaBar_,j beta^j / gammaBar
			+ gammaBar_,j beta^j_,i 
		) / gammaBar
	*/
	real3 DBar_tr_DBar_beta_L = real3_add(
		tr12_partial2_beta_L,
		real3_real_mul(
			real3_add3(
				sym3_real3_mul(
					partial2_det_gammaBar_LL,
					U->beta_U),
				real3_real_mul(
					real3_real_mul(
						partial_det_gammaBar_L,
						-real3_dot(partial_det_gammaBar_L, U->beta_U)),
					1. / det_gammaBar),
				partial_det_gammaBar_times_partial_beta_L
			), .5 / det_gammaBar)
	);
	
	//DBar_tr_DBar_beta_u.i = DBar^i DBar_k beta^k = gammaBar^ij DBar_j DBar_k beta^k
	real3 DBar_tr_DBar_beta_U = sym3_real3_mul(*gammaBar_UU, DBar_tr_DBar_beta_L);

	//tr_gammaBar_DHat2_beta_u.i = gammaBar^jk DHat_j DHat_k beta^i
	real3 tr_gammaBar_DHat2_beta_U = real3x3x3_sym3_dot23(DHat2_beta_ULL, *gammaBar_UU);

	real3 ABar_times_partial_alpha_U = sym3_real3_mul(*ABar_UU, *partial_alpha_L);
	real3 ABar_times_partial_W_U = sym3_real3_mul(*ABar_UU, *partial_W_L);

	real3 Delta_times_ABar_U = _3sym3_sym3_dot23(*Delta_ULL, *ABar_UU);

	real3 gammaBar_times_partial_K_U = sym3_real3_mul(*gammaBar_UU, *partial_K_L);

/*
L2 Lambda^I =
		+ gammaBar^JK DHat_J DHat_K beta^I
		+ 1/3 DBar^I DBar_J beta^J
		- 4/3 alpha gammaBar^IJ K_,J
		- 2 ABar^IJ alpha_,J 
		- 6 alpha ABar^IJ W_,J / W
		+ 2 alpha ABar^JK Delta^I_JK
*/
//// MODULE_DEPENDS: real3_add6	
	real3 L2_LambdaBar_U = real3_add6(
		tr_gammaBar_DHat2_beta_U,
		real3_real_mul(DBar_tr_DBar_beta_U, (1. / 3.)),
		real3_real_mul(gammaBar_times_partial_K_U, (-4. / 3.) * U->alpha),
		real3_real_mul(ABar_times_partial_alpha_U, -2.),
		real3_real_mul(ABar_times_partial_W_U, -6. * U->alpha / U->W),
		real3_real_mul(Delta_times_ABar_U, 2. * U->alpha)
	);

	return L2_LambdaBar_U;
#endif	// new way
<? end ?>	

	//ABar^IJ = e_i^I γBar^ik e_k^K ABar_KL e_l^L γBar^lj e_j^J
	sym3 ABar_UU;
<? for IJ,xIJ in ipairs(symNames) do
	local I,J,xI,xJ = from6to3x3(IJ)
?>	ABar_UU.<?=xIJ?> = 0.
<?	for K,xK in ipairs(xNames) do
		for L,xL in ipairs(xNames) do
?>		+ U->ABar_LL.<?=sym(K,L)?> * gammaBar_uu-><?=sym(K,L)?> * calc_len_<?=xI?>(pt) * calc_len_<?=xJ?>(pt) * calc_len_<?=xK?>(pt) * calc_len_<?=xL?>(pt)
<?		end
	end
?>	;
<? end
?>

	//partial2_beta_Ull[jk].I = beta^I_,jk
<?=eqn:makePartial2"beta_U"?>		

	real3 L2_LambdaBar_U;
<? 
local bleh = (0

-- beta^J part of 1/3 Dbar2betacontractionU^i	
	+ frac(1,3) * e"_i^I" * eInv"^k_J,kj"
	+ frac(1,6) * e"_i^I" * eInv"^k_J,j" * partial_det_gammaHat_l"_k" / det_gammaHat
	- frac(1,6) * e"_i^I" * eInv"^k_J" * partial_det_gammaHat_l"_k" * partial_det_gammaHat_l"_j" / det_gammaHat^2
	+ frac(1,6) * e"_i^I" * eInv"^k_J" * partial2_det_gammaHat_ll"_jk" / det_gammaHat

-- beta^J part of gammaBar^jk Dhat2betaUdDD^i_,jk
	+ eInv"^k_J" * e"_n^I" * partial_connHat_ulll"^n_jki"
	+ eInv"^k_J" * e"_n^I" * connHat_ull"^n_mi" * connHat_ull"^m_jk"
	- eInv"^k_J" * e"_n^I" * connHat_ull"^n_mk" * connHat_ull"^m_ji"
	+ 2 * e"_n^I" * connHat_ull"^n_jk" * eInv"^k_J,i"
	- e"_n^I" * connHat_ull"^k_ji" * eInv"^n_J,k"
	+ e"_n^I" * eInv"^n_J,ij"

)():permute"^I_Jij"

local bleh2 = (0

-- beta^J_,k part of 1/3 Dbar2betacontractionU^i	
	+ frac(1,6) * delta_ul"^k_j" * e"_i^I" * eInv"^m_J" * partial_det_gammaHat_l"_m" / det_gammaHat
	+ frac(1,3) * delta_ul"^k_j" * e"_i^I" * eInv"^m_J,m"
	+ frac(1,3) * e"_i^I" * eInv"^k_J,j"

-- beta^J_,k part of gammaBar^jk Dhat2betaUdDD^i_,jk
	+ 2 * delta_ul"^k_i" * e"_m^I" * connHat_ull"^m_jl" * eInv"^l_J"
	+ 2 * delta_ul"^k_i" * e"_m^I" * eInv"^m_J,j"
	- delta_ul"^I_J" * connHat_ull"^k_ij"

)():permute"^I_J^k_ij"

local bleh3 = (0
	+ frac(1,3) * e"_m^I" * eInv"^j_J"
	+ delta_ul"^I_J" * delta_ul"^j_m"
)():permute"^I_J^j_m"

for I,xI in ipairs(xNames) do
?>	L2_LambdaBar_U.<?=xI?> = 0.
<?	for	J,xJ in ipairs(xNames) do
?>	+ U->beta_U.<?=xJ?> * (0.
<?		for i,xi in ipairs(xNames) do
			for j,xj in ipairs(xNames) do
				if not Constant.isValue(bleh[I][J][i][j], 0) then
?>		+ gammaBar_uu-><?=sym(i,j)?> * (<?=eqn:compile(bleh[I][J][i][j])?>)
<?				end
			end
		end
?>	)
<?	end
	for	J,xJ in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>	+ partial_beta_Ul-><?=xk?>.<?=xJ?> * (0.
<?			for i,xi in ipairs(xNames) do
				for j,xj in ipairs(xNames) do
					if not Constant.isValue(bleh2[I][J][k][i][j], 0) then
?>		+ gammaBar_uu-><?=sym(i,j)?> * (<?=eqn:compile(bleh2[I][J][k][i][j])?>)
<?					end
				end
			end
?>	)
<?		end
	end
	for J,xJ in ipairs(xNames) do
		for j,xj in ipairs(xNames) do
			for k,xk in ipairs(xNames) do
				local jk = from3x3to6(j,k)
?>	+ partial2_beta_Ull[<?=jk-1?>].<?=xJ?> * (0.		
<?				for m,xm in ipairs(xNames) do
					if not Constant.isValue(bleh3[I][J][j][m], 0) then
?>		+ gammaBar_uu-><?=sym(m,k)?> * (<?=eqn:compile(bleh3[I][J][j][m])?>)
<?					end
				end
?>	)
<?			end
		end
	end
	for j,xj in ipairs(xNames) do
		local i,xi = I,xI	-- transformed by e_i^I so treat I like i
		local ij,xij = from3x3to6(i,j)
?>	- (4./3.) * U->alpha * gammaBar_uu-><?=xij?> * partial_K_l-><?=xj?> * calc_len_<?=xi?>(pt)
<?	end
	for J,xJ in ipairs(xNames) do
?>	- 2. * ABar_UU.<?=sym(I,J)?> * partial_alpha_L-><?=xJ?>
	- 6. * U->alpha * ABar_UU.<?=sym(I,J)?> * partial_W_L-><?=xJ?>
<?		for K,xK in ipairs(xNames) do
?>	+ 2. * U->alpha * U->ABar_LL.<?=sym(J,K)?> * Delta_ULL-><?=xI?>.<?=sym(J,K)?>
<?		end
	end
?>	;
<? 
end
?>
	return L2_LambdaBar_U;
}

//// MODULE_NAME: calc_PIRK_L3_LambdaBar_U

// NOTICE - all my PIRK_L3 functions DO NOT contribute Lbeta
static real3 calc_PIRK_L3_LambdaBar_U(
	real const tr_DBar_beta,
	real3 const * const Delta_U
) {
	//L3 LambdaBar^I = 2/3 Delta^I DBar_J beta^J
	real3 L3_LambdaBar_U = real3_real_mul(*Delta_U, (2. / 3.) * tr_DBar_beta);
	return L3_LambdaBar_U;
}

//// MODULE_NAME: calc_dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar

//another name for this could be d0_LambdaBar_U (2017 Ruchlin, eqn 15)
/*
	LambdaBar^i_,t = 
		+ gammaBar^jk DHat_j DHat_k beta^i
		+ 2/3 Delta^i DBar_j beta^j
		+ 1/3 DBar^i DBar_j beta^j
		- 2 ABar^ij (alpha_,j - 6 alpha phi_,j)
		+ 2 alpha ABar^jk Delta^i_jk
		- 4/3 alpha gammaBar^ij K_,j
Lbeta:	- LambdaBar^k beta^i_,k
		+ dup_k LambdaBar^i beta^k			<-- omitted from this function
source:	- 16 pi exp(4 phi) alpha S^i
	
	wrt W and IJK:
	LambdaBar^I_,t = 
L2:		+ 1/3 DBar^I DBar_J beta^J
		- 4/3 alpha gammaBar^IJ K_,J
		+ gammaBar^JK DHat_J DHat_K beta^I
		- 2 ABar^IJ (alpha_,J - 6 alpha phi_,J)
		+ 2 alpha ABar^JK Delta^I_JK
L3:		+ 2/3 Delta^I DBar_J beta^J
Lbeta:	- beta^I_,J LambdaBar^J
		+ beta^J dup_J(LambdaBar^I) 		<-- omitted from this function
source:	- 16 pi alpha S^I / W^2



SENR:
	LambdaBar^i_,t =
L2:		+ 1/3 Dbar2betacontractionU^i
		- 4/3 alpha gammaBar^ij K_,j
		+ gammaBar^jk Dhat2betaUdDD^i_,jk
		- 2 ABar^jk (
			+ delta^i_j alpha_,k 
			- 6 alpha delta^i_j phidD_,k 
			- alpha DGammaUDD^i_jk
		)
L3:		+ 2/3 DGammaU^i Dbarbetacontraction
Lbeta:	+ LbetaLambarU^i

detgammahat = defined in hatted quantities in senr.  grid metric
GammahatUDD^i_jk = defined in hatted quantities in senr.  grid connection.
detg = defined as aux var in senr. ratio between grid metric determinant and simulation metric determinant

Dbar2betacontractionD_k = beta^j_,jk
					+ 1/2 (
						det_gammaBar_,jk beta^j
						- det_gammaBar_,j det_gammaBar_,k beta^j / det_gammaBar
						+ det_gammaBar_,j beta^j_,k
					) / det_gammaBar

Dbar2betacontractionU^i = gammaBar^ik Dbar2betacontractionD_k

Dhat2betaUdDD^j_,ik = beta^j_,ik + GammaHat^j_id,k beta^d
					+ GammaHat^j_id beta^d_,k
					+ GammaHat^j_kd beta^d_,i 
					- GammaHat^d_ki beta^j_,d
					+ GammaHat^j_lk GammaHat^l_id beta^d
					- GammaHat^j_ld GammaHat^l_ik beta^d

phidD_,i = -W_,i / (2 W)

GammabarDDD_ljk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
GammabarUDD^i_jk = gammaBar^il GammabarDDD_ljk
DGammaUDD^i_jk = GammabarUDD^i_jk - GammahatUDD^i_jk
DGammaU^i = gammabar^jk DGammaUDD^i_jk
detgammabar = detgammahat * detg
detgammabardD_,i = detgammahatdD_,i * detg + detgammahat * detgdD_,i = (detgammahat * detg)_,i
Dbarbetacontraction = beta^i_,i + beta^i detgammabardD_,i / (2 detgammabar)

LbetaLambarU^i = -Lambda^j beta^i_,j
				+ beta^j dup_j Lambda^i 		<-- omitted from this function


SENR, expanded:
	LambdaBar^i_,t =
L2:		+ 1/3 gammaBar^ik (								----+ 1/3 Dbar2betacontractionU^i
			+ beta^j_,jk									|
			+ 1/2 (											|
				+ gammaBar_,jk beta^j						|
				- gammaBar_,j gammaBar_,k beta^j / gammaBar	|
				+ gammaBar_,j beta^j_,k						|
			) / gammaBar									|
		)												----+
		- 4/3 alpha gammaBar^ij K_,j
		+ gammaBar^jk (									----+ gammaBar^jk Dhat2betaUdDD^i_,jk
			+ beta^i_,jk + connHat^i_jl,k beta^l			|
			+ connHat^i_jl beta^l_,k						|
			+ connHat^i_kl beta^l_,j						|
			- connHat^l_kj beta^i_,l						|
			+ connHat^i_mk connHat^m_jl beta^l				|
			- connHat^i_ml connHat^m_jk beta^l				|
		)												----+
		- 2 ABar^jk (
			+ delta^i_j alpha_,k 
			+ 3 alpha delta^i_j W_,k / W
			- alpha Delta^i_jk
		)
L3:		+ 2/3 Delta^i (									----+ 2/3 DGammaU^i Dbarbetacontraction
			+ beta^j_,j 									|
			+ 1/2 beta^j gammaBar_,j / gammaBar				|
		)												----+
Lbeta:	- Lambda^j beta^i_,j
		+ beta^j dup_j Lambda^i 		<-- omitted from this function since B^i_,t doesn't use it, added later for LambdaBar^i_,t


now with source and normalization:
	ΛBar^I_,t =
L2:		
		+ β^J γBar^ij (												----+ #1 of 1/3 Dbar2betacontractionU^i
			+ 1/3 e_i^I e^k_J,kj										|
			+ 1/6 e_i^I e^k_J,j γHat_,k / γHat							|
			- 1/6 e_i^I e^k_J γHat_,k γHat_,j / γHat^2					|
			+ 1/6 e_i^I e^k_J γHat_,jk / γHat						----+
			
			+ e^k_J e_n^I ΓHat^n_jk,i								----+ #1 of γBar^jk Dhat2betaUdDD^i_,jk
			+ e^k_J e_n^I ΓHat^n_mi ΓHat^m_jk							|
			- e^k_J e_n^I ΓHat^n_mk ΓHat^m_ji							|
			+ 2 e_n^I ΓHat^n_jk e^k_J,i									|
			- e_n^I ΓHat^k_ji e^n_J,k									|
			+ e_n^I e^n_J,ij											|
		)					 										----+
		+ β^J_,k γBar^ij (                                      
			+ 1/6 δ^k_j e_i^I e^m_J γHat_,m / γHat					----+ #2 of 1/3 Dbar2betacontractionU^i
			+ 1/3 δ^k_j e_i^I e^m_J,m 									|
			+ 1/3 e_i^I e^k_J,j										----+
			+ 2 δ^k_i e_m^I ΓHat^m_jl e^l_J							----+ #2 of γBar^jk Dhat2betaUdDD^i_,jk
			+ 2 δ^k_i e_m^I e^m_J,j 									|
			- δ^I_J ΓHat^k_ij										----+
		)                                                       
		+ β^J_,jk γBar^km (                                     
			+ 1/3 e_m^I e^j_J										----+ #3 of 1/3 Dbar2betacontractionU^i
			+ δ^I_J δ^j_m 											----+ #3 of γBar^jk Dhat2betaUdDD^i_,jk
		)                                                       
		- 4/3 α γBar^ij K_,j e_i^I									----- -4/3 α γBar^ij K_,j
		- 2 e_i^I γBar^ij ABar_jk γBar^kl α_,l						----- -2 ABar^ij α_,j
		- 6 α γBar^ik ABar_KL γBar^lj W_,j / W e_i^I e_k^K e_l^L	----- -6 α ABar^ij W_,j / W
		+ 2 α ABar^JK Δ^I_JK										----- 2 α Delta^i_jk ABar^jk

L3:		+ 2/3 Δ^I (													----+ 2/3 DGammaU^i Dbarbetacontraction
			+ β^J_,j e^j_J												|
			+ β^J (e^j_J,j + 1/2 e^j_J γHat_,j / γBar)					|
		)															----+

Lbeta:	- ΛBar^J (
			+ β^I_,j e^j_J
			+ β^K (e_i^I e^i_K,j e^j_J)
		)
		+ β^J dup_j (ΛBar^I) e^j_J 									--+- omitted from this function since B^i_,t doesn't use it, added later for ΛBar^i_,t
		+ β^J ΛBar^K (e_i^I e^i_K,j e^j_J)							--/

source:	- 16 pi α S^I / W^2


*/
static inline real3 calc_dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const x,
	sym3 const * const gammaBar_uu,
	real3 const * const partial_alpha_L,
	real3 const * const partial_W_L,
	real3 const * const partial_K_l,
	real3x3 const * const partial_beta_Ul,
	real3x3 const * const partial_beta_UL,
	_3sym3 const * const Delta_ULL,				//FIXME remove this
	real3 const * const Delta_U,				//FIXME remove this
	real const tr_DBar_beta
) {
//// MODULE_DEPENDS: calc_PIRK_L2_LambdaBar_U	
	real3 const L2_LambdaBar_U = calc_PIRK_L2_LambdaBar_U(
		solver,
		U,
		x,
		gammaBar_uu,
		partial_beta_Ul,
		Delta_ULL,
		partial_alpha_L,
		partial_W_L,
		partial_K_l
	);

//// MODULE_DEPENDS: calc_PIRK_L3_LambdaBar_U
	real3 const L3_LambdaBar_U = calc_PIRK_L3_LambdaBar_U(
		tr_DBar_beta,
		Delta_U
	);

<?=eqn:makePartialUpwind"LambdaBar_U"?>
//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul
	real3x3 const partial_LambdaBar_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(
		U->LambdaBar_U, 	//TODO upwind avg
		partial_LambdaBar_Ul_upwind,
		x					//TODO upwind x
	);

//	real3 const Lbeta_LambdaBar_U = real3_sub(
//		real3x3_real3_mul(partial_LambdaBar_UL_upwind, U->beta_U),
//		real3x3_real3_mul(*partial_beta_UL, U->LambdaBar_U));

	real3 const dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar = real3_add4(
		
		L2_LambdaBar_U,	// L2
		
		L3_LambdaBar_U,	// L3

//no more LBeta
//		Lbeta_LambdaBar_U,
//instead, just the advect term:
		real3_neg(real3x3_real3_mul(*partial_beta_UL, U->LambdaBar_U)),

// source term:
		real3_real_mul(real3_rescaleFromCoord_u(U->S_u, x), -16. * M_PI * U->alpha / (U->W * U->W))
	);

	return dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar;
}

//// MODULE_NAME: calc_PIRK_L2_B_U

//////////////////////////////// B^I_,t //////////////////////////////// 

static real3 calc_PIRK_L2_B_U(
	constant <?=solver_t?> const * const solver,
	real3 const dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar
) {
	return real3_real_mul(
		dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar,
		solver->dt_beta_U_k
	);
}

//// MODULE_NAME: calc_PIRK_L3_B_U

// NOTICE - all my PIRK_L3 functions DO NOT contribute Lbeta
static real3 calc_PIRK_L3_B_U(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U
) {
	return real3_real_mul(U->B_U, -solver->dt_beta_U_eta);
}

//// MODULE_NAME: calcDeriv_K
//// MODULE_DEPENDS: calc_PIRK_L2_K calc_PIRK_L3_K

static inline void calcDeriv_K(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	//used for Lbeta_K
	int4 const updir,	//needed for eqn:makePartialUpwind
	real3 const x,
	//used for PIRK_L2_K
	sym3 const * const gammaBar_uu,
	sym3 const * const DBar2_alpha_ll,
	real3 const * const partial_alpha_l,
	real3 const * const partial_W_l,
	//used for PIRK_L3_K
	sym3 const * const ABar_UU,
	//used for source terms
	real const S
) {
	real const L2_K = calc_PIRK_L2_K(
		U,
		gammaBar_uu,
		DBar2_alpha_ll,
		partial_alpha_l,
		partial_W_l
	);

	real const L3_K = calc_PIRK_L3_K(
		U,
		ABar_UU
	);

<?=eqn:makePartialUpwind"K"?>
	real3 const partial_K_L_upwind = real3_rescaleFromCoord_l(partial_K_l_upwind, x);
	real const Lbeta_K = real3_dot(partial_K_L_upwind, U->beta_U);

	/*
	B&S 11.52
	2008 Alcubierre 2.8.12
	K_,t = -gamma^ij D_i D_j alpha + alpha (ABar_ij ABar^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	2017 Ruchlin et al
	K_,t = 
		1/3 alpha K^2 
		+ alpha ABar_ij ABar^ij 
		- exp(-4 phi) gammaBar^ij (
			DBar_i DBar_j alpha 
			+ 2 alpha_,i phi_,j
		) 
		+ K_,i beta^i
		+ 4 pi alpha (rho + S)
	
	wrt W:
	K_,t = 
L2:		- W^2 gammaBar^ij DBar_i DBar_j alpha 
		+ W gammaBar^ij alpha_,i * W_,j
L3:		+ 1/3 alpha K^2 
		+ alpha ABar_ij ABar^ij 			// this term is in L2 in the paper but L3 in code with note of the error
Lbeta:	+ K_,i beta^i
source:	+ 4 pi alpha (rho + S)
	*/
	deriv->K +=
		L2_K
		+ L3_K
		+ Lbeta_K
		+ 4. * M_PI * U->alpha * (U->rho + S);
}

//// MODULE_NAME: calcDeriv_ABar_LL

/*
2017 Ruchlin et al, eqn. 11b
	ABar_ij,t = 
L2:		exp(-4 phi) (trace-free part above)_ij
L3:		- 2/3 ABar_ij DBar_k beta^k
		+ alpha ABar_ij K
		- 2 alpha ABar_ik ABar^k_j
Lbeta:	+ beta^k_,i ABar_jk
		+ beta^k_,j ABar_ik
		+ ABar_ij,k beta^k

in normalized coordinates, with W:
	ABar_IJ,t = 
L2:		W^2 (trace-free part above)_IJ
L3:		- 2/3 ABar_IJ DBar_k beta^k
		+ alpha ABar_IJ K
		- 2 alpha ABar_IK ABar^K_J
Lbeta:	
		+ beta^K dup_k (ABar_IJ) e^k_K
		+ beta^K ABar_MN ((e_i^M e_j^N)_,k e^k_K e^i_I e^j_J)
		+ ABar_JK beta^K_,i e_j^J
		+ ABar_IK beta^K_,j e_i^I
		+ ABar_JK beta^L (e_k^K e^k_L,i e^j_J)
		+ ABar_IK beta^L (e_k^K e^k_L,j e^i_I)
*/
static void calcDeriv_ABar_LL(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const x,
	sym3 const * const gammaBar_LL,
	sym3 const * const gammaBar_UU,
	sym3 const * const DBar2_alpha_LL,
	real3 const * const partial_alpha_L,
	real3 const * const partial_W_L,
	_3sym3 const * const connBar_ULL,
	real3x3 const * const ABar_UL,
	real3x3 const * const partial_beta_Ul,
	real3 const * const partial_W_l,
	real const tr_DBar_beta,
	//for RBar_IJ:
	sym3 partial_epsilon_LLl[3],
	_3sym3 const * const Delta_ULL,
	real3 const * const Delta_U,
	_3sym3 const * const connHat_ULL,
	_3sym3 const * const partial_gammaBar_LLL,
	sym3 const * const gammaBar_ll,
	sym3 const * const gammaBar_uu
) {
//// MODULE_DEPENDS: calc_PIRK_L2_ABar_LL 
	sym3 L2_ABar_LL = calc_PIRK_L2_ABar_LL(
		solver,
		U,
		x,
		gammaBar_LL,
		gammaBar_UU,
		DBar2_alpha_LL,
		partial_alpha_L,
		partial_W_L,
		connBar_ULL,
		partial_epsilon_LLl,
		Delta_ULL,
		Delta_U,
		connHat_ULL,
		partial_gammaBar_LLL,
		gammaBar_ll,
		gammaBar_uu
	);
	
//// MODULE_DEPENDS: calc_PIRK_L3_ABar_LL
	sym3 L3_ABar_LL = calc_PIRK_L3_ABar_LL(
		U,
		tr_DBar_beta,
		gammaBar_uu,
		x
	);

<?=eqn:makePartialUpwind"ABar_LL"?>		
<? if false then ?>
//// MODULE_DEPENDS: calc_partial_ABar_LLL 
	//partial_ABar_lll[k].ij = ABar_ij,k
	_3sym3 partial_ABar_LLL_upwind = calc_partial_ABar_LLL(
		x,			//xup?
		U->ABar_LL,	//upwind avg?
		partial_ABar_LLl_upwind);

//// MODULE_DEPENDS: sym3_Lbeta_LL
	sym3 Lbeta_ABar_LL = sym3_Lbeta_LL(
		U->ABar_LL,
		partial_ABar_LLL_upwind,
		U->beta_U,
		*partial_beta_UL);
<? else ?>
	real3 const pt = x;
	sym3 Lbeta_ABar_LL;
<? 	for IJ,xIJ in ipairs(symNames) do 
		local I,J,xI,xJ = from6to3x3(IJ)
?>	Lbeta_ABar_LL.<?=xIJ?> = 0.
<? 		for K,xK in ipairs(xNames) do
?>		+ U->beta_U.<?=xK?> * (0.
			+ partial_ABar_LLl_upwind[<?=K-1?>].<?=xIJ?> / calc_len_<?=xK?>(pt)				//+ beta^K (dup_k ABar_IJ) e^k_K
<?			for M,xM in ipairs(xNames) do
				for N,xN in ipairs(xNames) do
					if not Constant.isValue(eSq_deSq_norm[I][J][K][M][N], 0) then
?>			+ U->ABar_LL.<?=sym(M,N)?> * calc_eSq_deSq_norm_<?=xI..xJ..xK..xM..xN?>(pt)	//+ beta^K ABar_MN (e_i^M e_j^N)_,k (e^k_K e^i_I e^j_J)
<?					end
				end
			end
?>		)
		+ U->ABar_LL.<?=sym(K,J)?> * partial_beta_Ul-><?=xI?>.<?=xK?> / calc_len_<?=xI?>(pt)	//+ ABar_KJ beta^K_,i e^i_I
		+ U->ABar_LL.<?=sym(K,I)?> * partial_beta_Ul-><?=xJ?>.<?=xK?> / calc_len_<?=xJ?>(pt)	//+ ABar_KI beta^K_,j e^j_J
<?			for L,xL in ipairs(xNames) do 
				if not Constant.isValue(e_deInv_eInv[K][L][I], 0) then
?>		+ U->ABar_LL.<?=sym(K,J)?> * U->beta_U.<?=xL?> * calc_e_deInv_eInv_<?=xK..xL..xI?>(pt)	//+ ABar_KJ beta^L (e^k_L_,i e^i_I e_k^K)
<?				end
				if not Constant.isValue(e_deInv_eInv[K][L][J], 0) then
?>		+ U->ABar_LL.<?=sym(K,I)?> * U->beta_U.<?=xL?> * calc_e_deInv_eInv_<?=xK..xL..xJ?>(pt)	//+ ABar_KI beta^L (e^k_L_,j e^j_J e_k^K)
<?				end
			end 
		end
?>	;
<? 	end 
?>
<? end ?>

	deriv->ABar_LL = sym3_add4(
		deriv->ABar_LL,
		L2_ABar_LL,
		L3_ABar_LL,
		Lbeta_ABar_LL
	);
}

<? if eqn.useScalarField then ?>

//// MODULE_NAME: calcDeriv_Phi

//////////////////////////////// Phi_,t, Psi_I,t //////////////////////////////// 

static void calcDeriv_Phi(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const x
) {
<?=eqn:makePartialUpwind"Phi"?>
	cplx3 partial_Phi_L_upwind = cplx3_rescaleFromCoord_l(partial_Phi_l_upwind, x);
	
	/*
	Phi_,t = alpha Pi + beta^i Phi_,i
	*/
	deriv->Phi = cplx_add3(
		deriv->Phi, 
		cplx_real_mul(U->Pi, U->alpha),
		cplx_from_real(real3_dot(cplx3_re(partial_Phi_L_upwind), U->beta_U))
	);
}

//// MODULE_NAME: calcDeriv_Psi

static void calcDeriv_Psi(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const x,
	real3 const partial_alpha_l,
	real3x3 const partial_beta_ul
) {
<?=eqn:makePartialUpwind"Psi_l"?>
<?=eqn:makePartial1"Pi"?>

	real3 beta_u = real3_rescaleToCoord_U(U->beta_U, x);

	/*
	Psi_i,t = alpha_,i Pi + alpha Pi_,i + beta^j Psi_i,j + beta^j_,i Phi_,j
	*/
	deriv->Psi_l = cplx3_add5(
		deriv->Psi_l,
		real3_cplx_mul(partial_alpha_l, U->Pi),
		cplx3_real_mul(partial_Pi_l, U->alpha),
		cplx3x3_real3_mul(partial_Psi_ll_upwind, beta_u),
		cplx3_real3x3_mul(U->Psi_l, partial_beta_ul));
}

//// MODULE_NAME: calcDeriv_Pi

static void calcDeriv_Pi(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const deriv,
	global <?=cons_t?> const * const U,
	int4 const updir,
	real3 const x,
	real const dt_alpha,
	real3 const dt_beta_U,
	real3 const partial_alpha_L,
	real3x3 const partial_beta_UL,
	sym3 const gammaBar_uu,
	_3sym3 const connHat_ULL,
	real3 const Delta_U,
	real3 const partial_W_l
) {

	cplx3 Psi_L = cplx3_from_real3_real3(
		real3_rescaleFromCoord_l(cplx3_re(U->Psi_l), x),
		real3_rescaleFromCoord_l(cplx3_im(U->Psi_l), x));

<?=eqn:makePartial1"Psi_l"?>
	cplx3x3 partial_Psi_LL = cplx3x3_partial_rescaleFromCoord_Ll(U->Psi_l, partial_Psi_ll, x);
	
	real exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	sym3 gamma_uu = sym3_real_mul(*gammaBar_uu, exp_neg4phi);
	sym3 gamma_UU = sym3_rescaleFromCoord_uu(gamma_uu, x);

	/*
	4conn^t = -1/alpha^3 (alpha_,t - beta^m alpha_,m + alpha^2 K)
	*/
	real _4conn_t = (
		(dt_alpha - real3_dot(U->beta_U, *partial_alpha_L)) / (U->alpha * U->alpha)
		+ U->K) / U->alpha;
	
	_3sym3 connHat_ull = _3sym3_rescaleToCoord_ULL(connHat_ULL, x);

	/*
	Delta^i_jk = connBar^i_jk - connHat^i_jk
	Delta^i = Delta^i_jk gammaBar^jk = connBar^i - connHat^i_jk gammaBar^jk
	connBar^i = Delta^i + connHat^i_jk gammaBar^jk
	*/
	real3 connHat_u = _3sym3_sym3_dot23(connHat_ull, *gammaBar_uu);
	real3 connHat_U = real3_rescaleFromCoord_u(connHat_u, x); 
	real3 connBar_U = real3_add(*Delta_U, connHat_U);

	/*
	2008 Alcubierre eqn 2.8.15
	conn^i = exp(-4 phi) (connBar^i - 2 gammaBar^ij phi_,j)
	= exp(-4 phi) connBar^i - 2 gamma^ij phi_,j
	= W^2 connBar^i + gamma^ij W_,j / W
	*/
	real3 conn_u = real3_add(
		real3_real_mul(
			connBar_U,
			exp_neg4phi
		),
		real3_real_mul(
			sym3_real3_mul(gamma_uu, *partial_W_l),
			1. / U->W
		)
	);
	real3 conn_U = real3_rescaleFromCoord_u(conn_u, x); 
	
	real3 partial_alpha_U = sym3_real3_mul(gamma_UU, *partial_alpha_L);

	/*
	4conn^i = 3conn^i 
		+ 1/alpha^3 beta^i (alpha_,t - beta^m alpha_,m + alpha^2 K) 
		- 1/alpha^2 (beta^i_,t - beta^m beta^i_,m + alpha alpha_,j gamma^ij)
	*/
	real3 _4conn_i_U = real3_add5(
		conn_U,
		real3_real_mul(U->beta_U, ((dt_alpha - real3_dot(U->beta_U, *partial_alpha_L)) / (U->alpha * U->alpha) + U->K) / U->alpha),
		real3_real_mul(*dt_beta_U, -1. / (U->alpha * U->alpha)),
		real3_real_mul(real3x3_real3_mul(partial_beta_UL, U->beta_U), 1. / (U->alpha * U->alpha)),
		real3_real_mul(partial_alpha_U, -1. / U->alpha)
	);

	/*
	Pi_,t = 
		- alpha_,t Pi / alpha^2
		+ alpha_,i beta^i Pi / alpha
		- conn^t Pi alpha^2
		- beta^i_,t Psi_i / alpha
		+ Psi_j beta^j_,i beta^i / alpha
		- alpha conn^i Psi_i
		- alpha conn^t beta^i Psi_i
		+ alpha gamma^ij Psi_i,j
		- alpha (mu^2 + lambda Phi^2) Phi
	*/
	deriv->Pi = cplx_add3(
		deriv->Pi,
		cplx_real_mul(
			cplx_add4(
				real_cplx_mul(real3_dot(U->beta_U, *partial_alpha_L), U->Pi),
				cplx_neg(real3_cplx3_dot(*dt_beta_U, Psi_L)),
				cplx3_real3_dot(Psi_L, real3x3_real3_mul(partial_beta_UL, U->beta_U)),
				real_cplx_mul(-dt_alpha / U->alpha, U->Pi)
			),
			1. / U->alpha
		),
		cplx_real_mul(
			cplx_add5(
				cplx_neg(real3_cplx3_dot(_4conn_i_U, Psi_L)),
				cplx_real_mul(cplx3_real3_dot(Psi_L, U->beta_U), -_4conn_t),
				cplx3x3_sym3_dot(partial_Psi_LL, gamma_UU),
			
				// dV/d|Phi|^2
				real_cplx_mul(-(solver->scalar_mu * solver->scalar_mu + solver->scalar_lambda * cplx_lenSq(U->Phi)), U->Phi),
				
				real_cplx_mul(-U->alpha * _4conn_t, U->Pi)
			),
			U->alpha
		)
	);
}

<? end	-- eqn.useScalarField ?>


<? if useSENRShiftAndCoDerivs then ?>
<?=eqn:template(file["hydro/eqn/bssnok-fd-num-inject-senr.cl"])?>
<? end -- useSENRShiftAndCoDerivs ?>


//// MODULE_NAME: <?=calcDeriv?>
//// MODULE_DEPENDS: <?=initCond_codeprefix?> <?=solver_macros?> eqn.macros

//TODO if we're calculating the constrains in the derivative
// then we do save calculations / memory on the equations
// but we also, for >FE integrators (which require multiple steps) are duplicating calculations
//TODO, what's the difference between 'calcDeriv' and '<?=addSource?>' in a finite-difference solver?
// or in any other solver for that matter?
kernel void <?=calcDeriv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useCalcDeriv then ?>
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

//// MODULE_DEPENDS: getUpwind 
	int4 const updir = getUpwind(U->beta_U);


	//////////////////////////////// alpha_,t //////////////////////////////// 

<?=eqn:makePartial1"alpha"?>			//partial_alpha_l[i] := alpha_,i
	real3 const partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

<? if useCalcDeriv_alpha then ?>

	real dt_alpha;
	{
<?=eqn:makePartialUpwind"alpha"?>
		real3 partial_alpha_L_upwind = real3_rescaleFromCoord_l(partial_alpha_l_upwind, x);

		/*
		2008 Alcubierre 4.2.52 - Bona-Masso family of slicing
		Q = f(alpha) K
		d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
		alpha_,t = -alpha^2 f(alpha) K + alpha_,i beta^i
		alpha_,t = -alpha^2 f(alpha) K + alpha_,i e^i_I beta^I
		 for f(alpha) = 2/alpha
		alpha_,t = -2 alpha K + alpha_,i e^i_I beta^I
		*/
		dt_alpha = -U->K * calc_f_alphaSq(U->alpha) + real3_dot(partial_alpha_L_upwind, U->beta_U);
	}
	
	deriv->alpha += dt_alpha;
<? end	-- useCalcDeriv_alpha ?>


	//////////////////////////////// W_,t //////////////////////////////// 

<?=eqn:makePartial1"beta_U"?>		//partial_beta_Ul.j..i := beta^I_,j
<?=eqn:makePartial1"W"?>			//partial_W_l.i := W_,i 
	real3 partial_W_L = real3_rescaleFromCoord_l(partial_W_l, x);

//// MODULE_DEPENDS: calc_partial*_det_gammaHat_over_det_gammaHat_* 
	/*
	notice: 
	det(gammaBar_ij)/det(gammaHat_ij) = det(gammaBar_IJ)
	*/
	real3 const partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);

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
	real const detg = 1.;
	real3 const partial_detg_L = real3_zero;
	sym3 const partial2_detg_LL = sym3_zero;
	
	real3 const partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);
	
//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?> 
	real const det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real const det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real const det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 const partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 

//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul 
	//partial_beta_UL.I.J := e_i^I beta^i_,j e^j_J = e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

<? if false then ?>
	
	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = real3x3_trace(partial_beta_UL);
	
	//tr_DBar_beta := DBar_j beta^j = beta^j_,j + connBar^j_kj beta^k	
	//ran sim to t=100, only dif with the det grad was 2e-10 ratio of the alpha value at t=10. so these two match.
	real tr_DBar_beta;
	{	//also in K_,t, so don't calculate this twice

//// MODULE_DEPENDS: <?=calc_gammaBar_LL?> 
		sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
		
		sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k

//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
		_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL 
		_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);
		
		real3 tr_connBar_L = _3sym3_tr12(connBar_ULL);

		tr_DBar_beta = tr_partial_beta + real3_dot(U->beta_U, tr_connBar_L);
	}

<? elseif false then ?>
	
	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = real3x3_trace(partial_beta_UL);
	
	//Etienne's SENR Mathematica notebook uses this instead:
	//tr_DBar_beta := beta^j_,j + beta^j gammaBar_,j / (2 gammaBar)
	real const tr_DBar_beta = tr_partial_beta 
		+ real3_dot(U->beta_U, partial_det_gammaBar_over_det_gammaBar_L) * .5;

<? else ?>
	
	//tr_DBar_beta := beta^J_,j e^j_J + beta^J (e^j_J_,j + 1/2 e^j_J gammaBar_,j / gammaBar)
<?
local bleh = (eInv"^j_J_,j"() + frac(1,2) * (eInv"^j_J" * partial_det_gammaHat_l"_j")() / det_gammaHat)()
?>
	real3 const pt = x;
	real const tr_DBar_beta = real3x3_trace(partial_beta_Ul)
<? 
for J,xJ in ipairs(xNames) do
	if not Constant.isValue(bleh[J], 0) then
?>		+ U->beta_U.<?=xJ?> * <?=eqn:compile(bleh[J])?>
<? 	end
end
?>	;

<? end ?>

<? if useCalcDeriv_W then ?>
//// MODULE_DEPENDS: calcDeriv_W
	calcDeriv_W(
		solver,
		deriv,
		U,
		updir,
		x,
		tr_DBar_beta);
<? end	-- useCalcDeriv_W ?>


	//////////////////////////////// K_,t //////////////////////////////// 

<?=eqn:makePartial1"epsilon_LL"?>	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial2"alpha"?>		//partial2_alpha_ll.ij := alpha_,ij

<?=eqn:makePartial1"K"?>				//partial_K_l.i := K,i
	real3 const partial_K_L = real3_rescaleFromCoord_l(partial_K_l, x);

//// MODULE_DEPENDS: <?=calc_gammaBar_LL?> 
	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

//// MODULE_DEPENDS: <?=calc_exp_neg4phi?> 
	real const exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL 
	_3sym3 const partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
//// MODULE_DEPENDS: calc_connBar_ULL 
	_3sym3 const connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	sym3 const partial2_alpha_LL = sym3_rescaleFromCoord_ll(partial2_alpha_ll, x);
	
	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	//used below for TF_DBar2_alpha_LL
	sym3 const DBar2_alpha_LL = sym3_sub(partial2_alpha_LL, real3_3sym3_dot1(partial_alpha_L, connBar_ULL));

	//ABar^i_j := gammaBar^ik ABar_kj
	real3x3 const ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar^ij := ABar^i_k gammaBar^kj
	sym3 const ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	//ABar^IJ = ABar^I_K gammaBar^KJ


<? if useCalcDeriv_K then ?>

#if 0
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);
#else	
	real const det_gammaHat = <?=calc_det_gammaHat?>(x);
	real const det_gammaBar = det_gammaHat;
	sym3 const gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);
	sym3 const gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
#endif	
	
	/*
	gammaBar_ij = exp(-4 phi) gamma_ij
	gammaBar^ij = exp(4 phi) gamma^ij
	gamma^ij = exp(-4 phi) gammaBar^ij
	S := S_ij gamma^ij = exp(-4 phi) S_ij gammaBar^ij 
	*/
	real const S = exp_neg4phi * sym3_dot(U->S_ll, gammaBar_uu);

#if 0	//BAD.  this does some divisions by r which cause numbers to diverge.
	sym3 const DBar2_alpha_ll = sym3_rescaleToCoord_LL(DBar2_alpha_LL, x);
#else	//GOOD.  no divisions (except maybe gammaBar^ij from gammaBar_ij.  works.
	
//// MODULE_DEPENDS: calc_len_# calc_partial*_len*
	/*
	partial_gammaBar_lll.k.ij := gammaBar_ij,k
	gammaBar_ij,k = ((epsilon_IJ + gammaHat_IJ) e_i^I e_j^J)_,k
	gammaBar_ij,k = epsilon_IJ_,k e_i^I e_j^J + (epsilon_IJ + gammaHat_IJ) (e_i^I_,k e_j^J + e_i^I e_j^J_,k)
	gammaBar_ij,k = epsilon_IJ_,k len_i len_j + gammaBar_IJ (len_i_,k len_j + len_i len_j_,k)
	*/
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j,xi,xj = from6to3x3(ij)
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_LLl[<?=k-1?>].<?=xij?> * calc_len_len_<?=xi..xj?>(x)
		+ gammaBar_LL.<?=xij?> * calc_partial_len_len_<?=xi..xj..xk?>(x);
<?	end
end
?>

	_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>	connBar_lll.<?=xi?>.<?=xjk?> = .5 * (0.
		+ partial_gammaBar_lll.<?=xk?>.<?=sym(i,j)?>
		+ partial_gammaBar_lll.<?=xj?>.<?=sym(i,k)?>
		- partial_gammaBar_lll.<?=xi?>.<?=sym(j,k)?>
	);
<?	end
end
?>
	_3sym3 const connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);
	sym3 const DBar2_alpha_ll = sym3_sub(partial2_alpha_ll, real3_3sym3_dot1(partial_alpha_l, connBar_ull));
#endif

//// MODULE_DEPENDS: calcDeriv_K 
	//TODO I changed these to be coordinate rather than non-coordinate
	// but it made no difference
	// so TODO put it back maybe?
	calcDeriv_K(
		solver,
		deriv,
		U,
		//used for Lbeta_K
		updir,
		x,
		//PIRK_L2_K
		&gammaBar_uu,
		&DBar2_alpha_ll,	//TODO
		&partial_alpha_l,
		&partial_W_l,
		//PIRK_L3_K
		&ABar_UU,			//TODO
		//K source terms
		S
	);

<? end	-- useCalcDeriv_K ?>


	//////////////////////////////// epsilon_ij_,t //////////////////////////////// 

<? if useCalcDeriv_epsilon_LL then ?>
//// MODULE_DEPENDS: calcDeriv_epsilon_LL 
	calcDeriv_epsilon_LL(
		solver,
		deriv,
		U,
		updir,
		x,
		&gammaBar_LL,
		&partial_beta_Ul,
		&ABar_UL,
		tr_DBar_beta
	);
<? end	-- useCalcDeriv_epsilon_LL ?>


	//////////////////////////////// ABar_ij_,t //////////////////////////////// 

//// MODULE_DEPENDS: calc_connHat_LLL_and_ULL 
	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK - connHat^I_JK
	_3sym3 const Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

<? if useCalcDeriv_ABar_LL then ?>
	
//// MODULE_DEPENDS: <?=mystery_C_U?>
	//Delta_U.I := e_i^I Delta^i = e_i^I (LambdaBar^i - C^i)
	real3 const Delta_U = real3_sub(U->LambdaBar_U, <?=mystery_C_U?>);

//// MODULE_DEPENDS: calcDeriv_ABar_LL
	calcDeriv_ABar_LL(
		solver,
		deriv,
		U,
		updir,
		x,
		&gammaBar_LL,
		&gammaBar_UU,
		&DBar2_alpha_LL,
		&partial_alpha_L, 
		&partial_W_L,
		&connBar_ULL,
		&ABar_UL,
		&partial_beta_Ul,
		&partial_W_l,
		tr_DBar_beta,
		//for RBar_IJ:
		partial_epsilon_LLl,
		&Delta_ULL,
		&Delta_U,
		&connHat_ULL,
		&partial_gammaBar_LLL,
		&gammaBar_ll,
		&gammaBar_uu
	);

<? end	-- useCalcDeriv_ABar_LL ?>


	//////////////////////////////// LambdaBar^i_,t //////////////////////////////// 

<? if useCalcDeriv_LambdaBar_U or useCalcDeriv_beta_U then ?>
	//debugging: SENR calculates the Delta^I from Delta^I_JK gammaBar^JK
	real3 const Delta_U_from_Delta_ULL = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);

//// MODULE_DEPENDS: calc_dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar 
	real3 const dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar 
	= calc_dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar(
		solver,
		U,
		updir,
		x,
		&gammaBar_uu,
		&partial_alpha_L,
		&partial_W_L,
		&partial_K_l,
		&partial_beta_Ul,
		&partial_beta_UL,
		&Delta_ULL,
		&Delta_U_from_Delta_ULL,
		tr_DBar_beta
	);

<? end ?>
<? if useCalcDeriv_LambdaBar_U then ?>

<?=eqn:makePartialUpwind"LambdaBar_U"?>

//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul 
	real3x3 partial_LambdaBar_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(
		U->LambdaBar_U, 	//TODO upwind avg
		partial_LambdaBar_Ul_upwind,
		x					//TODO upwind x
	);

#if 0
	real3 dt_LambdaBar_U = real3_add(
		dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar,
		real3x3_real3_mul(partial_LambdaBar_UL_upwind, U->beta_U));
#else	//written out
	real3 dt_LambdaBar_U = real3_zero;
	for (int j = 0; j < 3; ++j) {
		dt_LambdaBar_U.s[j] += dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar.s[j];
		for (int k = 0; k < 3; ++k) {
			dt_LambdaBar_U.s[j] += partial_LambdaBar_UL_upwind.v[k].s[j] * U->beta_U.s[k];
		}
	}
#endif

	deriv->LambdaBar_U = real3_add(deriv->LambdaBar_U, dt_LambdaBar_U);

<? end	-- useCalcDeriv_LambdaBar_U ?>


	//////////////////////////////// beta^i_,t and B^i_,t //////////////////////////////// 

<? if useCalcDeriv_beta_U then ?>
/*
	for k = 3/4
	B^i_,t = 
L2:		k (∂_0 ΛBar^i = ΛBar^i_,t - β^j ∂up_j ΛBar^i)
L3:		- eta B^i
(not a Lie derivative, just an advection)
Lβ:		+ β^j ∂up_j B^i 
	

OK to fully go over this, in 2017 Ruchlin:
eqn 10 defines ∂_⊥ = ∂_t - Lβ
eqn 11e define ∂_⊥ ΛBar^i = ... 
implying (∂_t + Lβ) ΛBar^i = ... 
∂_t ΛBar^i = ... - Lβ ΛBar^i
Lie derivatives of vectors are defined as LA(B) = A^i_,j B^j - A^j B^i_,j
∂_t ΛBar^i = ... - β^i_,j ΛBar^j + β^j ΛBar^i_,j

and that would mean that eqn 14b & 15 would imply ...
∂_0 = ∂_t - β^j ∂_j
so
∂_0 ΛBar^i = ∂_t ΛBar^i - β^j ΛBar^i_,j
so 
∂_0 ΛBar^i = ... - β^i_,j ΛBar^j + β^j ΛBar^i_,j - β^j ΛBar^i_,j
∂_0 ΛBar^i = ... - β^i_,j ΛBar^j

so we can rewrite
∂_t ΛBar^i = ∂_0 ΛBar^i + β^j ΛBar^i_,j

and last, eqn 14b:
∂_0 B^i = k ∂_0 ΛBar^i - η B^i
∂_0 B^i = k (∂_t ΛBar^i - β^j ΛBar^i_,j) - η B^i
∂_0 B^i = k (... - β^i_,j ΛBar^j) - η B^i

so if we store ΛBar^i_,t without Lβ (which is ∂_⊥ ΛBar^i)
or even the ∂_0 ΛBar^i then we can re-add the needed terms later


*/

<? if eqn.useShift == "GammaDriver" then ?>
	
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k LambdaBar^i_,t + eta LambdaBar^i
	//real const k = 3. / 4.;
	//real const eta = 1.;	//1.;	// 1 / (2 M), for total mass M
	real3 const dt_beta_U = real3_add(
		real3_real_mul(dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar, solver->dt_beta_U_k),
		real3_real_mul(U->LambdaBar_U, solver->dt_beta_U_eta));
	
<? elseif eqn.useShift == "HyperbolicGammaDriver" then ?>
	
<?=eqn:makePartialUpwind"beta_U"?>
//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul 
	real3x3 const partial_beta_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul_upwind, x);

	//SENR's 'shiftadvect':
	//advecting the shift vector: beta^i_,j beta^j
	real3 const partial_beta_times_beta_upwind = real3x3_real3_mul(partial_beta_UL_upwind, U->beta_U);

	/*
	hyperbolic Gamma driver 
	2017 Ruchlin et al, eqn 14a, 14b
	β^i_,t = 
L1:					B^i
advect shift field:	+ β^i_,j β^j
	*/
	real3 const dt_beta_U = real3_add(U->B_U, partial_beta_times_beta_upwind);


	//partial_B_ul[i] := B^i_,t
<?=eqn:makePartialUpwind"B_U"?>
//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul 
	//partial_B_UL.I.J := e_i^I (B^M e^i_M)_,j e^j_J
	real3x3 const partial_B_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(U->B_U, partial_B_Ul_upwind, x);

<? if false then ?>
//// MODULE_DEPENDS: calc_PIRK_L2_B_U 
	//L2 B^I = (k=3/4) ΛBar^I_,t
	real3 L2_B_U = calc_PIRK_L2_B_U(solver, dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar);
	
//// MODULE_DEPENDS: calc_PIRK_L3_B_U
	//L3 B^I = -eta B^I 
	real3 L3_B_U = calc_PIRK_L3_B_U(solver, U);


	real3 const partial_B_times_beta_upwind = real3x3_real3_mul(partial_B_UL_upwind, U->beta_U);

	deriv->B_U = real3_add4(
		deriv->B_U,
		L2_B_U,
		L3_B_U,
		partial_B_times_beta_upwind
	);
<? end ?>
<? if true then ?>
	
	// L2_B^i = 3/4 (L2_LambdaBar^i + L3_LambdaBar^i - beta^j dup_j LambdaBar^i)		
	// notice L3_LambdaBar^i includes Lbeta LambdaBar^i = beta^j dup_j LambdaBar^i - beta^i_,j LambdaBar^j
	// L2_B^i = 3/4 (L2_LambdaBar^i + L3_wo_Lbeta_LambdaBar^i - beta^i_,j LambdaBar^j)
	real3 const L2_B_U = real3_real_mul(dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar, solver->dt_beta_U_k);

	// L3_B^i = -eta B^i + beta^j dup_j B^i
	real3 const L3_B_U = real3_real_mul(U->B_U, -solver->dt_beta_U_eta);

	//this is part of L3_B_U
	real3 Ladvect_B_U = real3_zero;
	for (int j = 0; j < 3; ++j) {
		for (int k = 0; k < 3; ++k) {
			Ladvect_B_U.s[j] += U->beta_U.s[k] * partial_B_UL_upwind.v[k].s[j];
		}
	}

	for (int j = 0; j < 3; ++j) {
		deriv->B_U.s[j] += L2_B_U.s[j] + L3_B_U.s[j] + Ladvect_B_U.s[j];
	}
<? end ?>
<? end	-- eqn.useShift ?>

	deriv->beta_U = real3_add(deriv->beta_U, dt_beta_U);
<? end	-- useCalcDeriv_beta_U ?>


	//////////////////////////////// Phi_,t //////////////////////////////// 

<? if eqn.useScalarField then ?>

<? if useCalcDeriv_Phi then ?>
//// MODULE_DEPENDS: calcDeriv_Phi
	calcDeriv_Phi(solver, deriv, U, updir, x);
<? end	-- useCalcDeriv_Phi ?>


	//////////////////////////////// Psi_i,t //////////////////////////////// 

	real3x3 partial_beta_ul = real3x3_partial_rescaleToCoord_Ul(U->beta_U, partial_beta_Ul, x);

<? if useCalcDeriv_Psi then ?>
//// MODULE_DEPENDS: calcDeriv_Psi
	calcDeriv_Psi(solver, deriv, U, updir, x,
		partial_alpha_l,
		partial_beta_ul
	);
<? end	-- useCalcDeriv_Psi ?>


	//////////////////////////////// Pi_,t //////////////////////////////// 

<? if useCalcDeriv_Pi then ?>
//// MODULE_DEPENDS: calcDeriv_Pi
	calcDeriv_Pi(solver, deriv, U, updir, x,
		dt_alpha,
		dt_beta_U,
		partial_alpha_L,
		partial_beta_UL,
		gammaBar_uu,
		connHat_ULL,
		Delta_U,
		partial_W_l
	);
<? end	-- useCalcDeriv_Pi ?>
<? end	--eqn.useScalarField ?>


<? if useSENRShiftAndCoDerivs then ?>{
	<?=eqn:template(file["hydro/eqn/bssnok-fd-num-inject-senr-calcDeriv.cl"], {
		-- TODO how to forward the local env into the inline?
		-- until then, just forward the local vars that I use
		useCalcDeriv_alpha = useCalcDeriv_alpha,
		useCalcDeriv_W = useCalcDeriv_W,
		useCalcDeriv_K = useCalcDeriv_K,
		useCalcDeriv_epsilon_LL = useCalcDeriv_epsilon_LL,
		useCalcDeriv_ABar_LL = useCalcDeriv_ABar_LL,
		useCalcDeriv_LambdaBar_U = useCalcDeriv_LambdaBar_U,
		useCalcDeriv_beta_U = useCalcDeriv_beta_U,
	})?>
}<? end -- useSENRShiftAndCoDerivs ?>


//// MODULE_DEPENDS: applyKreissOligar numberof
	// Kreiss-Oligar dissipation:
	int fields[numIntStates];
	for (int i = 0; i < numberof(fields); ++i) fields[i] = i;
	applyKreissOligar(solver, U, cell, deriv, x, fields, numberof(fields));

<? end 	-- useCalcDeriv ?>
}


//// MODULE_NAME: <?=constrainU?>
//// MODULE_DEPENDS: <?=mystery_C_U?> calc_partial_ABar_LLL tracefree calc_RBar_LL <?=calc_exp_neg4phi?> calc_partial_gammaBar_LLL calc_connHat_LLL_and_ULL

kernel void <?=constrainU?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useConstrainU then ?>	
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const U = UBuf + index;
	
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

#if 1	//rescale based on non-coords
	real const det_gammaHatLL = 1.;
	real const rescaleMetric = cbrt(det_gammaHatLL/det_gammaBarLL);
<? 		for ij,xij in ipairs(symNames) do
?>	gammaBar_LL.<?=xij?> *= rescaleMetric;
<? 		end ?>
	U->epsilon_LL = sym3_sub(gammaBar_LL, sym3_ident);
	det_gammaBarLL = det_gammaHatLL;

#else	//rescale bases on coords
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);
	real det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real rescaleMetric = cbrt(1. / det_gammaBar_over_det_gammaHat);
	gammaBar_ll = sym3_real_mul(gammaBar_ll, rescaleMetric); 
	gammaBar_LL = sym3_rescaleFromCoord_ll(gammaBar_ll, x);
	U->epsilon_LL = sym3_sub(gammaBar_LL, sym3_ident);
	det_gammaBarLL = sym3_det(gammaBar_LL);
#endif

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

	real exp_neg4phi = <?=calc_exp_neg4phi?>(U);	//W^2
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);
	sym3 gamma_uu = sym3_real_mul(gammaBar_uu, exp_neg4phi);
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, 1. / exp_neg4phi);
	
<? if eqn.useScalarField and coupleTandPhi then ?>
	/*
	update T_ab values
	2017 Esorihuela-Tomas et al eqn 2
	T_ab = 1/2 (D_a Phi)* (D_b Phi) + 1/2 (D_a Phi) (D_b Phi)* - 1/2 g_ab g^cd (D_c Phi)* (D_d Phi) - 1/2 mu^2 g_ab |Phi Phi*| - 1/4 lambda g_ab |Phi Phi*|^2
	*/
	
	real Phi_lenSq = cplx_lenSq(U->Phi);
	real Psi_lenSq = cplx3_weightedLenSq(U->Psi_l, gamma_uu);
	real Pi_lenSq = cplx_lenSq(U->Pi);

	// V = mu^2 |Phi|^2 + 1/2 lambda |Phi|^4
	real V = (solver->scalar_mu * solver->scalar_mu + .5 * solver->scalar_lambda * Phi_lenSq) * Phi_lenSq;

	/*
	rho = n^a n^b T_ab = 1/2 (Pi Pi* + (Psi^c)* Psi_c + V) = 1/2 (re(Pi Pi*) + re(Psi^c Psi_c*) + V) = 1/2 (|Pi|^2 + |Psi|^2 + V)
	*/
	U->rho = .5 * (Pi_lenSq + Psi_lenSq + V);
	
	/*
	S^u = -1/2 gamma^ua ((Psi_a)* Pi + Psi_a Pi*)
	
	a + a* = re(a + a*) = 2 re(a)
	re(a*) = re(a)
	re((a b)*) = re(a b)
	so (a* b + a b*) = a b + (a b)* = re(a b + (a b)*) = re(a* b) + re(a b*) = 2 re(a* b) = 2(re(a) re(b) + im(a) im(b))

	a b* = (ar + i ai) (br - i bi) = (ar br + ai bi) + i (ai br - ar bi)
	a* b = (ar - i ai) (br + i bi) = (ar br + ai bi) - i (ai br - ar bi)
	
	a a* = (ar + i ai) (ar - i ai) = ar^2 + ai^2 = re(a a*) = re(a* a) = a* a

	S^u = -1/2 gamma^ua (re(Psi_a) re(Pi) + im(Psi_a) im(Pi))
	*/
	U->S_u = real3_neg(sym3_real3_mul(
		gamma_uu,
		real3_add(
			real3_real_mul(cplx3_re(U->Psi_l), U->Pi.re),
			real3_real_mul(cplx3_im(U->Psi_l), U->Pi.im)
		)
	));

	/*
	S_uv = 1/2 (Psi_u* Psi_v + Psi_u Psi_v*) - 1/2 gamma_uv (gamma^cd Psi_c* Psi_d - Pi* Pi + V)
	= re(Psi_u* Psi_v) + 1/2 gamma_uv (-gamma^cd re(Psi_c* Psi_d) + |Pi|^2 - V)
	*/
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	U->S_ll.<?=xij?> = 
		cplx_dot(U->Psi_l.<?=xi?>, U->Psi_l.<?=xj?>)
		+ .5 * gamma_ll.<?=xij?> * (
			-Psi_lenSq
			+ Pi_lenSq
			- V
		);
<?	end
?>

<? end -- useScalarField and coupleTandPhi ?>

	// update H and M values
<? if eqn.guiVars.calc_H_and_M and eqn.guiVars.calc_H_and_M.value then ?>

<?=eqn:makePartial1"W"?>				//partial_W_l[i] := phi_,i 

<?=eqn:makePartial2"W"?>			//partial2_W_ll.ij := phi_,ij

	sym3 partial2_phi_LL_times_WSq;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	partial2_phi_LL_times_WSq.<?=xij?> = .5 * (
		- partial2_W_ll.<?=xij?> * U->W / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
		+ partial_W_l.<?=xi?> / calc_len_<?=xi?>(x)
			* partial_W_l.<?=xj?> / calc_len_<?=xj?>(x)
	);
<? end ?>

	//partial_epsilon[k].ij := epsilon_ij,k = gammaBar_ij,k
<?=eqn:makePartial1"epsilon_LL"?>	
	
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar_uu.ij := ABar^ij = gammaBar^ik ABar_kl gammaBar^lj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	
	
	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	_3sym3 Delta_LLL = sym3_3sym3_mul(gammaBar_LL, Delta_ULL);
	
	real3 Delta_U = real3_sub(U->LambdaBar_U, <?=mystery_C_U?>);

	real3x3 partial_LambdaBar_UL;
	{
		//partial_LambdaBar_ul[j].i := connBar^i_,j
<?=eqn:makePartial1"LambdaBar_U"?>	
//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul 
		//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
		partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);
	}

<?=eqn:makePartial2"epsilon_LL"?>

	sym3 trBar_partial2_gammaBar_ll = calc_trBar_partial2_gammaBar_ll(
		U, 
		x, 
		gammaBar_UU, 
		partial_epsilon_LLl, 
		partial2_epsilon_LLll);

	sym3 RBar_LL = calc_RBar_LL(
		U,
		x,
		&gammaBar_LL,
		&gammaBar_UU,
		&connHat_ULL,
		&partial_gammaBar_LLL,
		&trBar_partial2_gammaBar_ll,
		&partial_LambdaBar_UL,
		&Delta_U,
		&Delta_ULL,
		&Delta_LLL);

	//RBar := RBar_ij gammaBar^ij
	real RBar = sym3_dot(gammaBar_UU, RBar_LL);
	
	//tr_DBar2_phi := gammaBar^ij DBar_i DBar_j phi = gammaBar^ij phi_,ij - connBar^k phi_,k
	real tr_DBar2_phi_times_WSq = 
		sym3_dot(
			gammaBar_UU,
			partial2_phi_LL_times_WSq
		)
		+ real3_dot(
			_3sym3_sym3_dot23(connBar_ULL, gammaBar_UU),
			partial_W_L
		) * .5 * U->W;
	
	//2017 Ruchlin et al, eqn 46
	//H = 2/3 K^2 - ABar^ij ABar_ij + exp(-4 phi) (RBar - 8 DBar^i phi DBar_i phi - 8 gammaBar^ij DBar_i DBar_j phi)
	// but I will scale down by 1/2 to match other equations
	U->H = (1. / 3.) * U->K * U->K
		- .5 * sym3_dot(U->ABar_LL, ABar_UU)
		+ .5 * (
			  RBar * U->W * U->W
			+ 4. * real3_weightedLenSq(partial_W_L, gammaBar_UU) * U->W
			- 8. * tr_DBar2_phi_times_WSq
		)
		- 8. * M_PI * U->rho;

#if 1

<?=eqn:makePartial1"ABar_LL"?>		//partial_ABar_lll[k].ij = ABar_ij,k
	_3sym3 partial_ABar_LLL = calc_partial_ABar_LLL(x, U->ABar_LL, partial_ABar_LLl);

	real3 partial_K_L;
	{
<?=eqn:makePartial1"K"?>				//partial_K_l.i := K,i
		partial_K_L = real3_rescaleFromCoord_l(partial_K_l, x);
	}

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
	but I will rescale the exp(6 phi) out of it so it matches the other eqns
	*/
	real exp_6phi = 1. / (U->W * U->W * U->W);
<? for i,xi in ipairs(xNames) do
?>	U->M_U.<?=xi?> = 
		- 8. * M_PI * U->S_u.<?=xi?> * calc_len_<?=xi?>(x)
<?	for j,xj in ipairs(xNames) do
?>		- 3. * ABar_UU.<?=sym(i,j)?> * partial_W_L.<?=xj?> / U->W
		- (2. / 3.) * exp_6phi * gammaBar_UU.<?=sym(i,j)?> * partial_K_L.<?=xj?>
<?		for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xi?>.<?=sym(j,k)?> * ABar_UU.<?=sym(j,k)?>
<?			for l,xl in ipairs(xNames) do
?>		- ABar_UU.<?=sym(k,j)?> * gammaBar_UU.<?=sym(l,i)?> * partial_gammaBar_LLL.<?=xj?>.<?=sym(k,l)?>
		- ABar_UU.<?=sym(k,i)?> * gammaBar_UU.<?=sym(l,j)?> * partial_gammaBar_LLL.<?=xj?>.<?=sym(k,l)?>
		+ gammaBar_UU.<?=sym(k,i)?> * gammaBar_UU.<?=sym(l,j)?> * partial_ABar_LLL.<?=xj?>.<?=sym(k,l)?>
<?			end
		end
	end
?>	;
<? end ?>
#else
	//2017 Ruchlin et al, eqn 47
	//M^i = exp(-4 phi) (DHat_j ABar^ij + 2 ABar^k(i Delta^j)_jk + 6 ABar^ij phi_,j - 2/3 gammaBar^ij K_,j)
#endif
<? end	-- calc_H_and_M ?>
<? end	-- useConstrainU ?>
}

//// MODULE_NAME: <?=addSource?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?> <?=cell_t?> <?=SETBOUNDS_NOGHOST?> <?=initCond_codeprefix?>

/*
TODO combine with calcDeriv
all this has in it is Phi and Psi for scalar field

I guess one important distinction between addSource and calcDeriv is
calcDeriv is replaced with the PIRK_L* derivative calculations
while addSource is always added to the RHS of both PIRK and finite-difference solvers
*/
kernel void <?=addSource?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
<? if useAddSource then ?>
	<?=SETBOUNDS_NOGHOST?>();
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;

<? if useScalarField then ?>
<?=eqn:makePartialUpwind"alpha"?>
	real3 partial_alpha_L_upwind = real3_rescaleFromCoord_l(partial_alpha_l_upwind, x);

	//2008 Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 f(alpha) K + alpha,i beta^i
	real dt_alpha = -calc_f_alphaSq(U->alpha) * U->K
		+ real3_dot(partial_alpha_L_upwind, U->beta_U);
	

<?=eqn:makePartial1"alpha"?>			//partial_alpha_l[i] := alpha_,i
	
	real3 partial_alpha_u = sym3_real3_mul(
		gamma_uu,
		*(real3*)partial_alpha_l
	);

	real Phi_lenSq = cplx_lenSq(U->Phi);
	real dV_dPhiSq = solver->scalar_mu * solver->scalar_mu + solver->lambda * Phi_lenSq;
	real phi_source = U->Phi * dV_dPhiSq;

	/*	
	\Pi (
		\frac{1}{\alpha} \alpha_{,t} (1 - \frac{1}{\alpha})
		+ K \alpha 
	)
	
	+ \Psi_i (
		\alpha_{,j} \gamma^{ij}
		- \alpha \cdot {}^{(3)} \Gamma^i 
	)

	- \alpha \frac{dV}{d|\Phi|^2} \Phi
	*/
	deriv->Pi += (
		U->Pi * (
			dt_alpha * (1. - 1. / U->alpha) / U->alpha
			+ U->alpha * U->K
		)
		+ real3_dot(partial_alpha_u, U->Psi_l)
		- U->alpha * (
			real3_dot(conn_u, U->Psi_l)
			+ phi_source
		)
	);

<?=eqn:makePartialUpwind"beta_U"?>
	real3x3 partial_beta_ul = real3x3_partial_rescaleToCoord_Ul(U->beta_U, partial_beta_Ul, x);
	
	//\alpha_{,i} \Pi + {\beta^k}_{,i} \Psi_k
	deriv->Psi_l = real3_add(deriv->Psi_l,
		real3_real_mul(*(real3*)partial_alpha_l, U->Pi),
		real3x3_real3_mul(partial_beta_ul, U->Psi_l)
	);


<? end -- useScalarField ?>
<? end -- useAddSource ?>
}

//// MODULE_NAME: <?=BSSNOK_PIRK?>
//// MODULE_DEPENDS: <?=initCond_codeprefix?>

kernel void copyWAlphaBeta(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcA,
	global <?=cons_t?> const * const srcB
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);

	dstBuf[index] = srcA[index];
	dstBuf[index].alpha = srcB[index].alpha;
	dstBuf[index].W = srcB[index].W;
	dstBuf[index].beta_U = srcB[index].beta_U;
}

kernel void copyLambdaBar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcA,
	global <?=cons_t?> const * const srcB
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);

	dstBuf[index] = srcA[index];
	dstBuf[index].LambdaBar_U = srcB[index].LambdaBar_U;
}

//// MODULE_DEPENDS: getUpwind calc_partial*_det_gammaHat_over_det_gammaHat_* real3x3_partial_rescaleFromCoord_Ul <?=calc_gammaBar_LL?> calcDeriv_epsilon_LL calcDeriv_W 

// epsilon_IJ, W, alpha, beta^I
kernel void calcDeriv_PIRK_L1_EpsilonWAlphaBeta(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	int4 const updir = getUpwind(U->beta_U);
	
<?=eqn:makePartial1"beta_U"?>		//partial_beta_Ul.j.i := beta^I_,j

	real3 const partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);


	//TODO detg ...
	real const detg = 1.;
	real3 const partial_detg_L = real3_zero;
	
	real3 const partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);

//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real const det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real const det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real const det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 const partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 

	//partial_beta_UL.I.J := e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 const partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

	//tr_partial_beta := beta^i_,i
	real const tr_partial_beta = real3x3_trace(partial_beta_UL);

	/*
	tr_DBar_beta := DBar_j beta^j = beta^j_,j + connBar^j_kj beta^k
	Etienne's SENR Mathematica notebook uses this instead:
	tr_DBar_beta := beta^j_,j + beta^j gammaBar_,j / (2 gammaBar)
	*/
	real const tr_DBar_beta = tr_partial_beta + real3_dot(U->beta_U, partial_det_gammaBar_over_det_gammaBar_L) * .5;
	
	
	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

	//////////////////////////////// epsilon_ij_,t //////////////////////////////// 
	
	//ABar^i_j := gammaBar^ik ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	deriv->epsilon_LL = sym3_zero;
	calcDeriv_epsilon_LL(
		solver,
		deriv,
		U,
		updir,
		x,
		&gammaBar_LL,
		&partial_beta_Ul,
		&ABar_UL,
		tr_DBar_beta
	);

	//////////////////////////////// W_,t //////////////////////////////// 

	deriv->W = 0.;
	calcDeriv_W(
		solver,
		deriv,
		U,
		updir,
		x,
		tr_DBar_beta);

	//////////////////////////////// alpha_,t //////////////////////////////// 

<?=eqn:makePartial1"alpha"?>			//partial_alpha_l[i] := alpha_,i
	real3 partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

<?=eqn:makePartialUpwind"alpha"?>
	real3 partial_alpha_L_upwind = real3_rescaleFromCoord_l(partial_alpha_l_upwind, x);

	//2008 Alcubierre 4.2.52 - Bona-Masso family of slicing
	//Q = f(alpha) K
	//d/dt alpha = -alpha^2 Q = alpha,t + alpha,i beta^i
	//alpha,t = -alpha^2 f(alpha) K + alpha,i beta^i
	real dt_alpha = -calc_f_alphaSq(U->alpha) * U->K
		+ real3_dot(partial_alpha_L_upwind, U->beta_U);
	
	deriv->alpha = dt_alpha;
	
	//////////////////////////////// beta^i_,t //////////////////////////////// 

<? if eqn.useShift == "GammaDriver" then ?>
	
	real3 dt_LambdaBar_U = real3_add(
		dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar,
		real3x3_real3_mul(partial_LambdaBar_UL_upwind, U->beta_U));
	
	//Gamma-driver
	//B&S 4.82
	//beta^i_,t = k LambdaBar^i_,t + eta LambdaBar^i
	deriv->beta_U = real3_add(
		real3_real_mul(dt_LambdaBar_U, solver->dt_beta_U_k),
		real3_real_mul(U->LambdaBar_U, solver->dt_beta_U_eta));
			
<? elseif eqn.useShift == "HyperbolicGammaDriver" then ?>
	
<?=eqn:makePartialUpwind"beta_U"?>
	real3x3 partial_beta_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul_upwind, x);

	//SENR's 'shiftadvect':
	//advecting the shift vector: beta^i_,j beta^j
	real3 partial_beta_times_beta_upwind = real3x3_real3_mul(partial_beta_UL_upwind, U->beta_U);

	/*
	hyperbolic Gamma driver 
	2017 Ruchlin et al, eqn 14a, 14b
	beta^i_,t = 
L1:				B^i 
advect shift:	+ beta^i_,j beta^j
	*/
	deriv->beta_U = real3_add(
		U->B_U,
		partial_beta_times_beta_upwind);
<? end	-- eqn.useShift ?>


//// MODULE_DEPENDS: applyKreissOligar numberof
	// Kreiss-Oligar dissipation:
	int fields[] = {0, 1, 3, 4, 5, 12, 13, 14, 15, 16, 17};	//TODO derive this from eqn.consVars, or ptr offsets / sizeof(real), rather than hardcoding here
	applyKreissOligar(solver, U, cell, deriv, x, fields, numberof(fields));
}

//// MODULE_DEPENDS: <?=calc_exp_neg4phi?> calc_connHat_LLL_and_ULL <?=mystery_C_U?> calc_PIRK_L2_ABar_LL calc_PIRK_L2_K

// ABar_IJ, K
kernel void calcDeriv_PIRK_L2_ABarK(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

<?=eqn:makePartial1"alpha"?>			//partial_alpha_l[i] := alpha_,i
	real3 const partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

<?=eqn:makePartial1"W"?>			//partial_W_l.i := W_,i 

	real const exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	
//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real const det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	
	sym3 const gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 const gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	
	real3 const partial_W_L = real3_rescaleFromCoord_l(partial_W_l, x);
	
<?=eqn:makePartial1"epsilon_LL"?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL 
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);
	
	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	
	real3 Delta_U = real3_sub(U->LambdaBar_U, <?=mystery_C_U?>);

	sym3 DBar2_alpha_LL;
	{
<?=eqn:makePartial2"alpha"?>		//partial2_alpha_ll.ij := alpha_,ij
		sym3 partial2_alpha_LL = sym3_rescaleFromCoord_ll(partial2_alpha_ll, x);
		
		//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
		//used below for TF_DBar2_alpha_LL
		DBar2_alpha_LL = sym3_sub(partial2_alpha_LL, real3_3sym3_dot1(partial_alpha_L, connBar_ULL));
	}

	deriv->ABar_LL = calc_PIRK_L2_ABar_LL(
		solver,
		U,
		x,
		&gammaBar_LL,
		&gammaBar_UU,
		&DBar2_alpha_LL,
		&partial_alpha_L,
		&partial_W_L,
		&connBar_ULL,
		partial_epsilon_LLl,
		&Delta_ULL,
		&Delta_U,
		&connHat_ULL,
		&partial_gammaBar_LLL
	);

	deriv->K = calc_PIRK_L2_K(
		U,
		&gammaBar_uu,
		&DBar2_alpha_ll,
		&partial_alpha_l,
		&partial_W_l
	);
}

//// MODULE_DEPENDS: sym3_Lbeta_LL calc_PIRK_L3_K calc_partial_ABar_LLL

// ABar_IJ, K
kernel void calcDeriv_PIRK_L3_ABarK(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);

	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar_uu.i.j := ABar^ij = ABar^i_k gammaBar^kj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);
	
	real3 partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);

	real detg = 1.;
	real3 partial_detg_L = real3_zero;
	
	real3 partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);
	
	real det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 

<?=eqn:makePartial1"beta_U"?>		//partial_beta_Ul.j..i := beta^I_,j
	
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

	//2013 Baumgarte et al eqn B2
//// MODULE_DEPENDS: calc_PIRK_L3_ABar_LL
	deriv->ABar_LL = calc_PIRK_L3_ABar_LL(
		U,
		tr_DBar_beta,
		&gammaBar_uu,
		x
	);

	deriv->K = calc_PIRK_L3_K(
		U,
		&ABar_UU
	);

<? if useLBetaWithPIRK then ?>

	int4 updir = getUpwind(U->beta_U);

<?=eqn:makePartialUpwind"ABar_LL"?>		//partial_ABar_lll[k].ij = ABar_ij,k
	_3sym3 partial_ABar_LLL_upwind = calc_partial_ABar_LLL(
		x,			//xup?
		U->ABar_LL,	//upwind avg?
		partial_ABar_LLl_upwind);

	sym3 Lbeta_ABar_LL = sym3_Lbeta_LL(
		U->ABar_LL,
		partial_ABar_LLL_upwind,
		U->beta_U,
		partial_beta_UL);

	deriv->ABar_LL = sym3_add(deriv->ABar_LL, Lbeta_ABar_LL);


<?=eqn:makePartialUpwind"K"?>
	real3 partial_K_L_upwind = real3_rescaleFromCoord_l(partial_K_l_upwind, x);
	real Lbeta_K = real3_dot(partial_K_L_upwind, U->beta_U);

	deriv->K += Lbeta_K;
<? end	-- useLBetaWithPIRK ?>

	
//// MODULE_DEPENDS: applyKreissOligar numberof
	// Kreiss-Oligar dissipation:
	int fields[] = {2, 18, 19, 20, 21, 22, 23};	//TODO derive this from eqn.consVars, or ptr offsets / sizeof(real), rather than hardcoding here
	applyKreissOligar(solver, U, cell, deriv, x, fields, numberof(fields));
}

// LambdaBar^I
kernel void calcDeriv_PIRK_L2_LambdaBar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

<?=eqn:makePartial1"alpha"?>			//partial_alpha_l[i] := alpha_,i
	real3 partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

<?=eqn:makePartial1"K"?>				//partial_K_l.i := K,i
	real3 partial_K_L = real3_rescaleFromCoord_l(partial_K_l, x);

//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1"beta_U"?>		//partial_beta_Ul.j..i := beta^I_,j
<?=eqn:makePartial1"W"?>			//partial_W_l.i := W_,i 
	
	//partial_beta_UL.I.J := e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

<?=eqn:makePartial1"epsilon_LL"?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL 
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	
	real3 partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);
	
	real detg = 1.;
	real3 partial_detg_L = real3_zero;
	sym3 partial2_detg_LL = sym3_zero;
	
	real3 partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);
	
	real det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 



	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar^ij := ABar^i_k gammaBar^kj
	sym3 const ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	//ABar^IJ = ABar^I_K gammaBar^KJ

//// MODULE_DEPENDS: calc_PIRK_L2_LambdaBar_U
	deriv->LambdaBar_U = calc_PIRK_L2_LambdaBar_U(
		solver,
		U,
		x,
		&gammaBar_uu,
		&partial_beta_Ul,
		&Delta_ULL,
		&partial_alpha_L,
		&partial_W_L,
		&partial_K_l
	);
}

//// MODULE_DEPENDS: calc_PIRK_L3_LambdaBar_U

// LambdaBar^I
kernel void calcDeriv_PIRK_L3_LambdaBar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	
	real3 partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);

	real detg = 1.;
	real3 partial_detg_L = real3_zero;
	
	real3 partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);
	
//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 

<?=eqn:makePartial1"beta_U"?>		//partial_beta_Ul.j..i := beta^I_,j
	//partial_beta_UL.I.J := e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = real3x3_trace(partial_beta_UL);


	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

<?=eqn:makePartial1"epsilon_LL"?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL 
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

	real3 Delta_U = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);

	/*
	tr_DBar_beta := DBar_j beta^j = beta^j_,j + connBar^j_kj beta^k
	Etienne's SENR Mathematica notebook uses this instead:
	tr_DBar_beta := beta^j_,j + beta^j gammaBar_,j / (2 gammaBar)
	*/
	real tr_DBar_beta = tr_partial_beta + real3_dot(U->beta_U, partial_det_gammaBar_over_det_gammaBar_L) * .5;

	deriv->LambdaBar_U = calc_PIRK_L3_LambdaBar_U(
		tr_DBar_beta,
		&Delta_U
	);


<? if useLBetaWithPIRK then ?>
	int4 updir = getUpwind(U->beta_U);

//notice that, despite LambdaBar^I_,t only using one of the two Lie derivative terms, L3 LambdaBar^I uses both
<?=eqn:makePartialUpwind"LambdaBar_U"?>
	real3x3 partial_LambdaBar_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(
		U->LambdaBar_U, 	//TODO upwind avg
		partial_LambdaBar_Ul_upwind,
		x					//TODO upwind x
	);

	real3 Lbeta_LambdaBar_U = real3_sub(
		real3x3_real3_mul(partial_LambdaBar_UL_upwind, U->beta_U),
		real3x3_real3_mul(partial_beta_UL, U->LambdaBar_U));
	
	deriv->LambdaBar_U = real3_add(
		deriv->LambdaBar_U,
		Lbeta_LambdaBar_U);
<? end	-- useLBetaWithPIRK ?>


//// MODULE_DEPENDS: applyKreissOligar numberof
	// Kreiss-Oligar dissipation:
	int fields[] = {9, 10, 11};	//TODO derive this from eqn.consVars, or ptr offsets / sizeof(real), rather than hardcoding here
	applyKreissOligar(solver, U, cell, deriv, x, fields, numberof(fields));
}

//// MODULE_DEPENDS: calc_PIRK_L2_B_U 

// B^I
kernel void calcDeriv_PIRK_L2_B(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	int4 const updir = getUpwind(U->beta_U);

<?=eqn:makePartial1"alpha"?>			//partial_alpha_l[i] := alpha_,i
	real3 partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

<?=eqn:makePartial1"K"?>				//partial_K_l.i := K,i
	real3 partial_K_L = real3_rescaleFromCoord_l(partial_K_l, x);
	
//// MODULE_DEPENDS: <?=calc_det_gammaBarLL?>
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	sym3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1"beta_U"?>		//partial_beta_Ul.j..i := beta^I_,j
<?=eqn:makePartial1"W"?>			//partial_W_l.i := W_,i 
	real3 partial_W_L = real3_rescaleFromCoord_l(partial_W_l, x);
	
	//partial_beta_UL.I.J := e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

<?=eqn:makePartial1"epsilon_LL"?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
//// MODULE_DEPENDS: calc_connBar_ULL 
	_3sym3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

	//debugging: SENR calculates the Delta^I from Delta^I_JK gammaBar^JK
	real3 Delta_U = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);


	real3 partial_det_gammaHat_over_det_gammaHat_L = calc_partial_det_gammaHat_over_det_gammaHat_L(x);

	real detg = 1.;
	real3 partial_detg_L = real3_zero;
	sym3 partial2_detg_LL = sym3_zero;
	
	real3 partial_det_gammaBar_over_det_gammaHat_L = real3_add(
		real3_real_mul(partial_det_gammaHat_over_det_gammaHat_L, detg),
		partial_detg_L);
	
	real det_gammaBar_over_det_gammaHat = det_gammaBarLL;
	real det_gammaHat_over_det_gammaBar = 1. / det_gammaBar_over_det_gammaHat;
	real3 partial_det_gammaBar_over_det_gammaBar_L = real3_real_mul(
			partial_det_gammaBar_over_det_gammaHat_L, 
			det_gammaHat_over_det_gammaBar); 


	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = real3x3_trace(partial_beta_UL);

	/*
	tr_DBar_beta := DBar_j beta^j = beta^j_,j + connBar^j_kj beta^k
	Etienne's SENR Mathematica notebook uses this instead:
	tr_DBar_beta := beta^j_,j + beta^j gammaBar_,j / (2 gammaBar)
	*/
	real tr_DBar_beta = tr_partial_beta + real3_dot(U->beta_U, partial_det_gammaBar_over_det_gammaBar_L) * .5;


	real exp_neg4phi = <?=calc_exp_neg4phi?>(U);


	//ABar_ul.i.j := ABar^i_j = gammaBar^kl ABar_kj
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	
	//ABar^ij := ABar^i_k gammaBar^kj
	sym3 ABar_UU = real3x3_sym3_to_sym3_mul(ABar_UL, gammaBar_UU);	//ABar^IJ = ABar^I_K gammaBar^KJ

//// MODULE_DEPENDS: calc_dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar
	real3 dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar 
	= calc_dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar(
		solver,
		U,
		updir,
		x,
		&gammaBar_uu,
		&partial_alpha_L,
		&partial_W_L,
		&partial_K_l,
		&partial_beta_Ul,
		&partial_beta_UL,
		&Delta_ULL,
		&Delta_U,
		tr_DBar_beta
	);

	deriv->B_U = calc_PIRK_L2_B_U(solver, dt_LambdaBar_U_wo_partial_upwind_beta_LambdaBar);
}

//// MODULE_DEPENDS: calc_PIRK_L3_B_U

// B^I
kernel void calcDeriv_PIRK_L3_B(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	global <?=cell_t?> const * const cell = cellBuf + index;
	real3 const x = cell->pos;
	global <?=cons_t?> * const deriv = derivBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;

	deriv->B_U = calc_PIRK_L3_B_U(solver, U);


<? if useLBetaWithPIRK then ?>
	//partial_beta_Ul.j..i := beta^I_,j
<?=eqn:makePartial1"beta_U"?>
	
	//partial_beta_UL.I.J := e_i^I (beta^M e^i_M)_,j e^j_J
	real3x3 partial_beta_UL = real3x3_partial_rescaleFromCoord_Ul(U->beta_U, partial_beta_Ul, x);

	int4 updir = getUpwind(U->beta_U);

<?=eqn:makePartialUpwind"B_U"?>
	real3x3 partial_B_UL_upwind = real3x3_partial_rescaleFromCoord_Ul(
		U->B_U, 	//TODO upwind avg
		partial_B_Ul_upwind,
		x					//TODO upwind x
	);

	real3 Lbeta_B_U = real3_sub(
		real3x3_real3_mul(partial_B_UL_upwind, U->beta_U),
		real3x3_real3_mul(partial_beta_UL, U->B_U));

	deriv->B_U = real3_add(
		deriv->B_U,
		Lbeta_B_U);
<? end	-- useLBetaWithPIRK ?>


//// MODULE_DEPENDS: applyKreissOligar numberof
	// Kreiss-Oligar dissipation:
	int fields[] = {6, 7, 8};	//TODO derive this from eqn.consVars, or ptr offsets / sizeof(real), rather than hardcoding here
	applyKreissOligar(solver, U, cell, deriv, x, fields, numberof(fields));
}
	
//dst = src + deriv * dt
#define PIRK_EQ1(type, x)	dst->x = type##_add(src->x, type##_real_mul(deriv->x, dt))

//epsilon_IJ, W, alpha, beta^I
kernel void PIRK_Eq1_EpsilonWAlphaBeta(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivBuf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const deriv = derivBuf + index;

	PIRK_EQ1(real, alpha);
	PIRK_EQ1(real, W);
	PIRK_EQ1(real3, beta_U);
	PIRK_EQ1(sym3, epsilon_LL);
}
	
//dst = src + (derivL3_n + (derivL2_n + derivL2_1) * .5) * dt
#define PIRK_EQ2(type, x) dst->x = type##_add(src->x, type##_real_mul(type##_add(derivL3_n->x, type##_real_mul(type##_add(derivL2_n->x, derivL2_1->x), .5)), dt))

//LambdaBar^I
kernel void PIRK_Eq2_LambdaBar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivL2_nBuf,
	global <?=cons_t?> const * const derivL2_1Buf,
	global <?=cons_t?> const * const derivL3_nBuf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const derivL2_n = derivL2_nBuf + index;
	global <?=cons_t?> const * const derivL2_1 = derivL2_1Buf + index;
	global <?=cons_t?> const * const derivL3_n = derivL3_nBuf + index;

	PIRK_EQ2(real3, LambdaBar_U);
}

//ABar_IJ, K
kernel void PIRK_Eq2_ABarK(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivL2_nBuf,
	global <?=cons_t?> const * const derivL2_1Buf,
	global <?=cons_t?> const * const derivL3_nBuf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const derivL2_n = derivL2_nBuf + index;
	global <?=cons_t?> const * const derivL2_1 = derivL2_1Buf + index;
	global <?=cons_t?> const * const derivL3_n = derivL3_nBuf + index;

	PIRK_EQ2(real, K);
	PIRK_EQ2(sym3, ABar_LL);	
}

//B^I
kernel void PIRK_Eq2_B(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivL2_nBuf,
	global <?=cons_t?> const * const derivL2_1Buf,
	global <?=cons_t?> const * const derivL3_nBuf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const derivL2_n = derivL2_nBuf + index;
	global <?=cons_t?> const * const derivL2_1 = derivL2_1Buf + index;
	global <?=cons_t?> const * const derivL3_n = derivL3_nBuf + index;

	PIRK_EQ2(real3, B_U);
}

//dst = .5 * (U + U1 + derivL1_1 * dt)
#define PIRK_EQ3(type, x)	dst->x = type##_real_mul(type##_add3(U->x, U1->x, type##_real_mul(derivL1_1->x, dt)), .5);

//epsilon_IJ, W, alpha, beta^I
kernel void PIRK_Eq3_EpsilonWAlphaBeta(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const UBuf,
	global <?=cons_t?> const * const U1Buf,
	global <?=cons_t?> const * const derivL1_1Buf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const U = UBuf + index;
	global <?=cons_t?> const * const U1 = U1Buf + index;
	global <?=cons_t?> const * const derivL1_1 = derivL1_1Buf + index;

	PIRK_EQ3(real, alpha);
	PIRK_EQ3(real, W);
	PIRK_EQ3(real3, beta_U);
	PIRK_EQ3(sym3, epsilon_LL);
}

//dst = src + (derivL2_n + derivL2_next + derivL3_n + derivL3_1) * .5 * dt
#define PIRK_EQ4(type, x)	dst->x = type##_add(src->x, type##_real_mul(type##_add4(derivL2_n->x, derivL2_next->x, derivL3_n->x, derivL3_1->x), .5 * dt))

//LambdaBar^I
kernel void PIRK_Eq4_LambdaBar(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivL2_nBuf,
	global <?=cons_t?> const * const derivL2_nextBuf,
	global <?=cons_t?> const * const derivL3_nBuf,
	global <?=cons_t?> const * const derivL3_1Buf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const derivL2_n = derivL2_nBuf + index;
	global <?=cons_t?> const * const derivL2_next = derivL2_nextBuf + index;
	global <?=cons_t?> const * const derivL3_n = derivL3_nBuf + index;
	global <?=cons_t?> const * const derivL3_1 = derivL3_1Buf + index;

	PIRK_EQ4(real3, LambdaBar_U);
}

//ABar_IJ, K
kernel void PIRK_Eq4_ABarK(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivL2_nBuf,
	global <?=cons_t?> const * const derivL2_nextBuf,
	global <?=cons_t?> const * const derivL3_nBuf,
	global <?=cons_t?> const * const derivL3_1Buf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const derivL2_n = derivL2_nBuf + index;
	global <?=cons_t?> const * const derivL2_next = derivL2_nextBuf + index;
	global <?=cons_t?> const * const derivL3_n = derivL3_nBuf + index;
	global <?=cons_t?> const * const derivL3_1 = derivL3_1Buf + index;

	PIRK_EQ4(sym3, ABar_LL);
	PIRK_EQ4(real, K);
}

//B^I
kernel void PIRK_Eq4_B(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const dstBuf,
	global <?=cons_t?> const * const srcBuf,
	global <?=cons_t?> const * const derivL2_nBuf,
	global <?=cons_t?> const * const derivL2_nextBuf,
	global <?=cons_t?> const * const derivL3_nBuf,
	global <?=cons_t?> const * const derivL3_1Buf,
	real const dt
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	global <?=cons_t?> * const dst = dstBuf + index;
	global <?=cons_t?> const * const src = srcBuf + index;
	global <?=cons_t?> const * const derivL2_n = derivL2_nBuf + index;
	global <?=cons_t?> const * const derivL2_next = derivL2_nextBuf + index;
	global <?=cons_t?> const * const derivL3_n = derivL3_nBuf + index;
	global <?=cons_t?> const * const derivL3_1 = derivL3_1Buf + index;

	PIRK_EQ4(real3, B_U);
}
