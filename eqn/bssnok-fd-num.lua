--[[
Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer" 2010
Alcubierre "Introduction to Numerical Relativity" 2008

then I'm applying 2017 Ruchlin changes...
*) separate gammaBar_ll - gammaHat_ll = epsilon_ll
*) coordinate-transform beta^i, epsilon_ij, ABar_ij, LambdaBar^i to eliminate singularities from the metric

tensors are denoted with suffixes _u _l etc for upper and lower
rescaled tensors are denoted _U _L etc

Then I'm double checking all against (and borrowing heavily from) Zach Etienne's SENR: https://math.wvu.edu/~zetienne/SENR/
--]]
local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local EinsteinEqn = require 'eqn.einstein'
local makestruct = require 'eqn.makestruct'
local applyCommon = require 'common'
local time, getTime = table.unpack(require 'util.time')

local makePartials = require 'eqn.makepartial'

local BSSNOKFiniteDifferenceEquation = class(EinsteinEqn)
BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 
BSSNOKFiniteDifferenceEquation.hasEigenCode = true
BSSNOKFiniteDifferenceEquation.hasCalcDTCode = true
BSSNOKFiniteDifferenceEquation.hasFluxFromConsCode = true
BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true

-- not used with finite-difference schemes anyways
BSSNOKFiniteDifferenceEquation.weightFluxByGridVolume = false

--[[
args:
	useShift = 'none'
				'GammaDriver'
				'HyperbolicGammaDriver' (default)
--]]
function BSSNOKFiniteDifferenceEquation:init(args)
	-- options:
	-- needs to be defined up front
	-- otherwise rebuild intVars based on it ...
	self.useShift = args.useShift or 'HyperbolicGammaDriver'

	local intVars = table{
		{name='alpha', type='real'},			-- 1
		{name='beta_U', type='real3'},		 	-- 3: beta^i
		{name='epsilon_LL', type='sym3'},		-- 6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='W', type='real'},				-- 1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},				-- 1: K = K^i_i
		{name='ABar_LL', type='sym3'},			-- 6: ABar_ij, only 5 dof since ABar^k_k = 0
		{name='LambdaBar_U', type='real3'},		-- 3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
												-- TODO what is C^i ?
	}
	if self.useShift == 'HyperbolicGammaDriver' then
		intVars:insert{name='B_U', type='real3'}
	end

	self.consVars = table()
	:append(intVars)
	:append{
		--hyperbolic variables:
		--real3 a;				//3: a_i
		--_3sym3 dBar;			//18: dBar_ijk, only 15 dof since dBar_ij^j = 0
		--real3 Phi;			//3: Phi_i

		--stress-energy variables:
		{name='rho', type='real'},				--1: n_a n_b T^ab
		{name='S_u', type='real3'},				--3: -gamma^ij n_a T_aj
		{name='S_ll', type='sym3'},				--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},				--1
		{name='M_u', type='real3'},				--3
	}
	self.numIntStates = makestruct.countScalars(intVars)
	
	-- call construction / build structures	
	BSSNOKFiniteDifferenceEquation.super.init(self, args)
end

function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaBar', value=true, compileTime=true},
		--{name='constrain_det_gammaBar', value=false, compileTime=true},

		--{name='constrain_tr_ABar', value=true, compileTime=true},
		{name='constrain_tr_ABar', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},

		{name='shift_eta', value=1},	--1, or 1 / (2 M), for total mass M
	}
end

function BSSNOKFiniteDifferenceEquation:makePartial(field, fieldType, nameOverride)
	local derivOrder = 2 * self.solver.numGhost
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	return makePartials.makePartial(derivOrder, self.solver, field, fieldType, nameOverride)
end

function BSSNOKFiniteDifferenceEquation:makePartial2(field, fieldType, nameOverride)
	local derivOrder = 2 * self.solver.numGhost
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	return makePartials.makePartial2(derivOrder, self.solver, field, fieldType, nameOverride)
end

function BSSNOKFiniteDifferenceEquation:compile(expr)
	local symmath = require 'symmath'
	local coords = self.solver.coord.coords
	local var = symmath.var
	
	local function isInteger(x) return x == math.floor(x) end
	expr = expr:map(function(x)
		if symmath.op.pow.is(x)
		and symmath.Constant.is(x[2])
		and isInteger(x[2].value)
		and x[2].value > 1
		and x[2].value < 100
		then
			return setmetatable(table.rep({x[1]}, x[2].value), symmath.op.mul)
		end
	end)
	
	return symmath.export.C(expr
		:replace(coords[1], var'x.x')
		:replace(coords[2], var'x.y')
		:replace(coords[3], var'x.z'))
end

function BSSNOKFiniteDifferenceEquation:getEnv()
	return applyCommon{
		eqn = self,
		solver = self.solver,
	}
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template([[
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
local partial_len_ll = Tensor'_ij'
partial_len_ll['_ij'] = solver.coord.lenExprs'_i,j'()				-- derivative is last
local partial2_len_lll = Tensor'_ijk'
partial2_len_lll['_ijk'] = partial_len_ll'_ij,k'():factorDivision()	-- derivative is last
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


//connHat^i_jk,l := partial_connHat_ulll[l].i.jk
void calc_partial_connHat_ulll(_3sym3 partial_connHat_ulll[3], real3 x) {
<? 
local dGamma = solver.coord.Gamma_ull'^i_jk,l'()
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		for l,xl in ipairs(xNames) do
?>	partial_connHat_ulll[<?=l-1?>].<?=xi?>.<?=xjk?> = <?=eqn:compile(dGamma[i][j][k][l])?>;
<?		end
	end
end
?>}

real3 calc_partial_det_gammaHat_l(real3 x) {
	real3 partial_det_gammaHat_l;
<? 
local partial_det_gammaHat_l = Tensor('_i', function(i)
	return solver.coord.det_g:diff(solver.coord.coords[i])()
end)
for i,xi in ipairs(xNames) do
?>	partial_det_gammaHat_l.<?=xi?> = <?=eqn:compile(partial_det_gammaHat_l[i])?>;
<? end
?>	return partial_det_gammaHat_l;
}

sym3 calc_partial2_det_gammaHat_ll(real3 x) {
	sym3 partial2_det_gammaHat_ll;
<?
local partial2_det_gammaHat_ll = partial_det_gammaHat_l'_i,j'():factorDivision()
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_det_gammaHat_ll.<?=xij?> = <?=eqn:compile(partial2_det_gammaHat_ll[i][j])?>;
<?
end
?>	return partial2_det_gammaHat_ll;
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
#define calc_det_gammaBarLL(x) 1

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

]], self:getEnv())
end

-- requires: partial_ABar_LLl
function BSSNOKFiniteDifferenceEquation:getCode_partial_ABar_LLL()
	return template([[
	
	//derivative first
	//partial_ABar_lll.k.ij := ABar_ij,k
	//partial_ABar_LLL.K.IJ := e_i^I e_j^J e_k^K ABar_ij,k
	_3sym3 partial_ABar_LLL;
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_ABar_LLL.<?=xk?>.<?=xij?> = (0.
		+ partial_ABar_LLl[<?=k-1?>].<?=xij?>
		+ (
			calc_partial_len_<?=xi..xk?>(x) / calc_len_<?=xi?>(x)
			+ calc_partial_len_<?=xj..xk?>(x) / calc_len_<?=xj?>(x)
		) * U->ABar_LL.<?=xij?>
	) / calc_len_<?=xk?>(x);
<?	end
end
?>
]], self:getEnv())
end

-- requires: partial_epsilon_LLl
function BSSNOKFiniteDifferenceEquation:getCode_partial_gammaBar_LLL()
	return template([[
	
	//derivative first
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	//partial_gammaBar_LLL.k.IJ := e_i^I e_j^J e_k^K gammaBar_ij,k
	_3sym3 partial_gammaBar_LLL;
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_gammaBar_LLL.<?=xk?>.<?=xij?> = (0.
		+ partial_epsilon_LLl[<?=k-1?>].<?=xij?>
		+ (
			calc_partial_len_<?=xi..xk?>(x) / calc_len_<?=xi?>(x)
			+ calc_partial_len_<?=xj..xk?>(x) / calc_len_<?=xj?>(x)
		) * (U->epsilon_LL.<?=xij?><?=i==j and ' + 1.' or ''?>)
	) / calc_len_<?=xk?>(x);
<?	end
end
?>
]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_connHat_LLL_and_ULL()
 	return template([[
	
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
?>	partial_gammaHat_LLL.<?=xk?>.<?=xij?> = 2. * calc_partial_len_<?=xi..xk?>(x) / (calc_len_<?=xi?>(x) * calc_len_<?=xk?>(x));
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
	_3sym3 connHat_LLL;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>		connHat_LLL.<?=xi?>.<?=xjk?> = .5 * (0.

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
	_3sym3 connHat_ULL = connHat_LLL;
]], self:getEnv())
end


function BSSNOKFiniteDifferenceEquation:getCode_connBar_ULL()
 	return template([[
	_3sym3 connBar_ULL;
	{
		//connBar_lll.i.jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
		//connBar_LLL.I.JK := e^i_I e^j_J e^k_K connBar_ijk
		_3sym3 connBar_LLL;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k,xj,xk = from6to3x3(jk)
?>		connBar_LLL.<?=xi?>.<?=xjk?> = .5 * (0.
			+ partial_gammaBar_LLL.<?=xk?>.<?=sym(i,j)?>
			+ partial_gammaBar_LLL.<?=xj?>.<?=sym(i,k)?>
			- partial_gammaBar_LLL.<?=xi?>.<?=sym(j,k)?>
		);
<?	end
end
?>		//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
		connBar_ULL = sym3_3sym3_mul(gammaBar_UU, connBar_LLL);
	}
]], self:getEnv())
end

-- requires: gammaBar_uu, partial_epsilon_LLl, partial2_epsilon_LLll, 
function BSSNOKFiniteDifferenceEquation:getCode_trBar_partial2_gammaBar_LL()
	return template([[

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
				+ 2. * (calc_partial_len_<?=xi..xk?>(x) / calc_len_<?=xi?>(x)) * (calc_partial_len_<?=xj..xl?>(x) / calc_len_<?=xj?>(x)) 
			)
			+ 2. * partial_epsilon_LLl[<?=l-1?>].<?=xij?> * (0.
				+ calc_partial_len_<?=xi..xk?>(x) / calc_len_<?=xi?>(x)
				+ calc_partial_len_<?=xj..xk?>(x) / calc_len_<?=xj?>(x)
			)
			+ partial2_epsilon_LLll[<?=kl-1?>].<?=xij?>
		) * gammaBar_uu.<?=xkl?>
<?		end
	end
?>	;
<? end
?>

]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_RBar_LL()
	return template([[

	//DHat_gammaBar_lll.k.ij = DHat_k gammaBar_ij 
	// = gammaBar_ij,k - connHat^l_ki gammaBar_lj - connHat^l_kj gammaBar_il
	_3sym3 DHat_gammaBar_LLL;
<?
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i], xNames[j]
?>	DHat_gammaBar_LLL.<?=xk?>.<?=xij?> = 0.
		+ partial_gammaBar_LLL.<?=xk?>.<?=xij?>
<?		for l,xl in ipairs(xNames) do
?>		- connHat_ULL.<?=xl?>.<?=sym(k,i)?> * gammaBar_LL.<?=sym(l,j)?>
		- connHat_ULL.<?=xl?>.<?=sym(k,j)?> * gammaBar_LL.<?=sym(l,i)?>
<?		end
?>	;
<?	end
end
?>
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	_3sym3 partial_gammaBar_lll = _3sym3_rescaleToCoord_LLL(partial_gammaBar_LLL, x);

	_3sym3 connHat_ull = _3sym3_rescaleToCoord_ULL(connHat_ULL, x);

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
//I'm moving this term down below
//		+ partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>
		- gammaBar_ll.<?=sym(m,j)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,i)?>
		- gammaBar_ll.<?=sym(m,i)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,j)?>
		- connHat_ull.<?=xm?>.<?=sym(k,i)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,j)?>
		- connHat_ull.<?=xm?>.<?=sym(k,j)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,i)?>
<?			end
?>	;
<?		end
	end
end
?>

	_3sym3 DHat_gammaBar_lll = _3sym3_rescaleToCoord_LLL(DHat_gammaBar_LLL, x);

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
?>	trBar_DHat2_gammaBar_without_partial2_gammaBar_ll.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * (0.
			+ partial_DHat_gammaBar_without_partial2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>
		)
		+ gammaBar_uu.<?=sym(k,l)?> * (0.
<?			for m,xm in ipairs(xNames) do
?>			- connHat_ull.<?=xm?>.<?=sym(l,k)?> * DHat_gammaBar_lll.<?=xm?>.<?=sym(i,j)?>
			- connHat_ull.<?=xm?>.<?=sym(l,i)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(m,j)?>
			- connHat_ull.<?=xm?>.<?=sym(l,j)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(i,m)?>
<?			end
?>		)
<?		end
	end
?>	;
<?
end
?>

	//trBar_DHat2_gammaBar_ll.ij := gammaBar^kl DHat_k DHat_l gammaBar_ij
	sym3 trBar_DHat2_gammaBar_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	trBar_DHat2_gammaBar_LL.<?=xij?> = 0.
		+ trBar_partial2_gammaBar_LL.<?=xij?>
		+ trBar_DHat2_gammaBar_without_partial2_gammaBar_ll.<?=xij?>
//hmm, moving this division into trBar_DHat2_gammaBar_without_partial2_gammaBar_ll
//is causing it to diverge near r=1
			 / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x));
<? end
?>

	//derivative is the last index, unlike the partial_*'s
	//DHat_LambdaBar_ul.i.j := DHat_j LambdaBar^i = LambdaBar^i_,j + connHat^i_jk LambdaBar^k
	real3x3 DHat_LambdaBar_UL;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	DHat_LambdaBar_UL.<?=xi?>.<?=xj?> = 0.
		+ partial_LambdaBar_UL.<?=xi?>.<?=xj?>
<?		for k,xk in ipairs(xNames) do
?>		+ connHat_ULL.<?=xi?>.<?=sym(j,k)?> * U->LambdaBar_U.<?=xk?>
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
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	RBar_LL.<?=xij?> = 0.
			- .5 * trBar_DHat2_gammaBar_LL.<?=xij?>	
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * gammaBar_LL.<?=sym(i,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_LL.<?=sym(j,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xi?>
			+ .5 * Delta_U.<?=xk?> * (
				Delta_LLL.<?=xi?>.<?=sym(k,j)?> 
				+ Delta_LLL.<?=xj?>.<?=sym(k,i)?>
			)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>			+ gammaBar_UU.<?=sym(k,l)?> * (0.
				+ Delta_ULL.<?=xm?>.<?=sym(k,i)?> * Delta_LLL.<?=xj?>.<?=sym(m,l)?>
				+ Delta_ULL.<?=xm?>.<?=sym(k,j)?> * Delta_LLL.<?=xi?>.<?=sym(m,l)?>
				+ Delta_ULL.<?=xm?>.<?=sym(i,k)?> * Delta_LLL.<?=xm?>.<?=sym(j,l)?>
			)
<?			end
		end
	end
?>
	;
<? end ?>
]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([=[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	sym3 gammaHat_ll = calc_gammaHat_ll(x);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = gammaHat_ll;	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

<? if eqn.initState.name == 'Minkowski' then ?>
	
	setFlatSpace(solver, U, x);

<? else -- not Minkowski ?>

	<?=code?>

	U->alpha = alpha;
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);

	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);
	
	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar = calc_det_gammaBar(x); 

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

<? end -- Minkowski ?>

	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0.;
	U->M_u = real3_zero;
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize LambdaBar_u
//TODO do this symbolically.  That's what I originally did, but symbolic calculations were getting complex
// however, with spherical BSSN, you need to 
kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

<? 
-- TODO get Minkowski working without this
if eqn.initState.name == 'Minkowski' then
?>
	U->LambdaBar_U = real3_zero;

<? else	-- initState == Minkowski 
?>
#if 0	
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		setFlatSpace(solver, U, x);
		return;
	}
#endif
 
 <?=eqn:makePartial'epsilon_LL'?>
	
	//derivative first
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	//partial_gammaBar_LLL.k.IJ := e_i^I e_j^J gammaBar_ij,k
	_3sym3 partial_gammaBar_LLL;
<? 
for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
	for k,xk in ipairs(xNames) do
?>	partial_gammaBar_LLL.<?=xk?>.<?=xij?> = (0.
		+ partial_epsilon_LLl[<?=k-1?>].<?=xij?>
		+ (
			calc_partial_len_<?=xi..xk?>(x) * calc_len_<?=xj?>(x)
			+ calc_partial_len_<?=xj..xk?>(x) * calc_len_<?=xi?>(x)
		) * (U->epsilon_LL.<?=xij?><?=i==j and ' + 1.' or ''?>)
	) / calc_len_<?=xk?>(x);
<?	end
end
?>

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:getCode_connBar_ULL()?>
<?=eqn:getCode_connHat_LLL_and_ULL()?>
	
	//Delta^i_jk = connBar^i_jk - connHat^i_jk
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

	U->LambdaBar_U = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);

<? end -- initState == Minkowski ?>

}
]=], self:getEnv())
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd-num.cl'

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

BSSNOKFiniteDifferenceEquation.predefinedDisplayVars = {
-- [=[
	'U alpha',
--	'U beta_U mag',
	'U beta_U x',
	'U beta_U y',
	'U beta_U z',
--	'U B_U mag',
	'U B_U x',
	'U B_U y',
	'U B_U z',
	--'U epsilon_LL norm',
	'U epsilon_LL xx',
	'U epsilon_LL xy',
	'U epsilon_LL xz',
	'U epsilon_LL yy',
	'U epsilon_LL yz',
	'U epsilon_LL zz',
	'U W',
	'U K',
--	'U ABar_LL tr weighted',
--] =]	
	'U ABar_LL xx',
	'U ABar_LL xy',
	'U ABar_LL xz',
	'U ABar_LL yy',
	'U ABar_LL yz',
	'U ABar_LL zz',
-- [=[	
--	'U LambdaBar_U mag',
	'U LambdaBar_U x',
	'U LambdaBar_U y',
	'U LambdaBar_U z',
	'U H',
--	'U M_u mag',
	'U M_u x',
	'U M_u y',
	'U M_u z',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_LL tr weighted',
--]=]	

-- [[ debugging derivatives
-- [=[	
	'deriv alpha',
	'deriv beta_U x',
	'deriv beta_U y',
	'deriv beta_U z',
	'deriv B_U x',
	'deriv B_U y',
	'deriv B_U z',
	'deriv epsilon_LL xx',
	'deriv epsilon_LL xy',
	'deriv epsilon_LL xz',
	'deriv epsilon_LL yy',
	'deriv epsilon_LL yz',
	'deriv epsilon_LL zz',
	'deriv W',
	'deriv K',
--] =]	
	'deriv ABar_LL xx',
	'deriv ABar_LL xy',
	'deriv ABar_LL xz',
	'deriv ABar_LL yy',
	'deriv ABar_LL yz',
	'deriv ABar_LL zz',
-- [=[	
	'deriv LambdaBar_U x',
	'deriv LambdaBar_U y',
	'deriv LambdaBar_U z',
	'deriv H',
	'deriv M_u x',
	'deriv M_u y',
	'deriv M_u z',
--]=]	
--]]

	--'U tr_DBar2_phi',
	--'U DBar_phi_sq',
	--'U ABarSq tr weighted',

-- [[ should be zero for Minkowski
	'U RBar_LL xx',
	'U RBar_LL xy',
	'U RBar_LL xz',
	'U RBar_LL yy',
	'U RBar_LL yz',
	'U RBar_LL zz',
--]]
--[[ should be zero for Minkowski
	'U Delta^x xx',
	'U Delta^x xy',
	'U Delta^x xz',
	'U Delta^x yy',
	'U Delta^x yz',
	'U Delta^x zz',
	'U Delta^y xx',
	'U Delta^y xy',
	'U Delta^y xz',
	'U Delta^y yy',
	'U Delta^y yz',
	'U Delta^y zz',
	'U Delta^z xx',
	'U Delta^z xy',
	'U Delta^z xz',
	'U Delta^z yy',
	'U Delta^z yz',
	'U Delta^z zz',
--]]
--[[ should be diag(0, 2, 2 cos(theta)^2)
	'U trBar_partial2_gammaBar_LL xx',
	'U trBar_partial2_gammaBar_LL xy',
	'U trBar_partial2_gammaBar_LL xz',
	'U trBar_partial2_gammaBar_LL yy',
	'U trBar_partial2_gammaBar_LL yz',
	'U trBar_partial2_gammaBar_LL zz',
--]]
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)
	local env = self:getEnv()

	vars:append{
		{
			name = 'gamma_ll',
			type = 'sym3',
			code = [[
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	*value_sym3 = gamma_ll;
]], 
		},
		{name='gammaHat_ll', code=[[	*value_sym3 = calc_gammaHat_ll(x);]], type='sym3'},
		{name='gammaBar_ll', code=[[	*value_sym3 = calc_gammaBar_ll(U, x);]], type='sym3'},
		{name='gamma_uu', code=[[	*value_sym3 = calc_gamma_uu(U, x);]], type='sym3'},
		{name='gammaHat_uu', code=[[	*value_sym3 = calc_gammaHat_uu(x);]], type='sym3'},
		{name='gammaBar_uu', code=[[	*value_sym3 = calc_gammaBar_uu(U, x);]], type='sym3'},
		{name='K_ll', code=[[
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	*value_sym3 = sym3_real_mul(
		sym3_add(
			sym3_rescaleToCoord_LL(U->ABar_LL, x),
			sym3_real_mul(gammaBar_ll, U->K / 3.)
		), exp_4phi);
]], type='sym3'},

		{name='det gammaBar - det gammaHat', code=[[
	*value = sym3_det(calc_gammaBar_ll(U, x)) - calc_det_gammaBar(x);
]]},
		{name='det gamma based on phi', code=[[
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = det_gamma;
]]},
		
		{name='S', code='*value = sym3_dot(U->S_ll, calc_gamma_uu(U, x));'},
		
		{
			name='volume', 
			code=[[
	//|g| = exp(12 phi) |g_grid|
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = U->alpha * det_gamma;
]],
		},
		{name='f', code='*value = calc_f(U->alpha);'},
		{name='df/dalpha', code='*value = calc_dalpha_f(U->alpha);'},
	
--[=[	
		{
			name = 'ABarSq',
			type = 'sym3',
			code = [[
	sym3 gammaBar_UU = calc_gammaBar_UU(U, x);
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	sym3 ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, ABar_UL);
	*value_sym3 = ABarSq_LL;
]],
		},
		
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = template([[
	
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>
<?=eqn:getCode_connBar_ULL()?>

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
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

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	sym3 partial2_phi_LL = sym3_rescaleFromCoord_ll(partial2_phi_ll, x);

	//tr_DBar2_phi := gammaBar^ij DBar_i DBar_j phi = gammaBar^ij phi_,ij - connBar^k phi_,k
	real tr_DBar2_phi = 0.
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		local ij,xij = from3x3to6(i,j)
?>		+ gammaBar_UU.<?=sym(i,j)?> * (
			partial2_phi_LL.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?		end
?>		)
<?	end
end
?>	;

	*value = tr_DBar2_phi;
]], env)
		},


		{
			name = 'DBar2_phi_ll',
			code = template([[
<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
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

<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>
<?=eqn:getCode_connBar_ULL()?>

	//DBar2_phi := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_LL.<?=xij?> = partial2_phi_LL.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?		end
?>	;
<?	end
end
?>	

	*value = sym3_dot(gammaBar_UU, DBar2_phi_LL);
]], env),
 		},
 
		{
			name = 'partial_phi_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'W'?>
	real3 partial_phi_l;
<? for i,xi in ipairs(xNames) do
?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>
	
	*value_real3 = partial_phi_l;
]], env),
		},

		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'alpha'?>
	*value_real3 = *(real3*)partial_alpha_l;
]], env),
		},

		{
			name = 'DBar2_alpha_ll',
			code = template([[
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>
<?=eqn:getCode_connBar_ULL()?>

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_LL.<?=xij?> = partial2_alpha_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_alpha_L.<?=xk?>
<?	end
?>	;
<? end
?>

	*value = sym3_dot(gammaBar_UU, DBar2_alpha_LL);
]], env),
		},
	
		{
			name = 'tracelessPart_LL',
			type = 'sym3',
			code = template([[
	
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBar(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);

<?=eqn:getCode_connBar_ULL()?>

	sym3 partial2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	partial2_alpha_LL.<?=xij?> = partial2_alpha_ll[<?=ij-1?>] / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x));
<? end
?>
	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_LL.<?=xij?> = partial2_alpha_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_alpha_L.<?=xk?>
<?	end
?>	;
<? end
?>

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
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

	sym3 partial2_phi_LL = sym3_rescaleFromCoord_ll(partial2_phi_ll, x);

	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_LL.<?=xij?> = partial2_phi_LL.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?	end
?>	;
<? end
?>

	real3 partial_alpha_L;
<? for i,xi in ipairs(xNames) do
?>	partial_alpha_L.<?=xi?> = partial_alpha_l[<?=i-1?>] / calc_len_<?=xi?>(x);
<? end
?>

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	
	sym3 tracelessPart_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	tracelessPart_LL.<?=xij?> = (0.
		+ 2. * partial_phi_L.<?=xi?> * partial_alpha_L.<?=xj?>
			+ 2. * partial_phi_L.<?=xj?> * partial_alpha_L.<?=xi?>
			+ U->alpha * (0.
				- 2. * DBar2_phi_LL.<?=xij?>
				+ 4. * partial_phi_L.<?=xi?> * partial_phi_L.<?=xj?>
				- 8. * M_PI * U->S_ll.<?=xij?> / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
			)
		);
<? end
?>
	tracelessPart_LL = tracefree(tracelessPart_LL, gammaBar_LL, gammaBar_UU);

	*value_sym3 = tracelessPart_LL; 
]], env),
		},
--]=]	
	}

	local derivOrder = 2 * self.solver.numGhost
--[=[
--[[ expansion:
2003 Thornburg:  ... from Wald ...
Theta = n^i_;i + K_ij n^i n^j - K
= n^i_,i + Gamma^i_ji n^j + K_ij (n^i n^j - gamma^ij)
... in ADM: n^i = -beta^i / alpha ...
= (-beta^i / alpha)_,i + Gamma^i_ji (-beta^j / alpha) + K_ij (beta^i beta^j / alpha^2 - gamma^ij)
= -beta^i_,i / alpha
	+ beta^i alpha_,i / alpha^2
	- beta^i (1/2 |g|_,i / |g|) / alpha
	+ K_ij beta^i beta^j / alpha^2
	- K

Gamma^j_ij = (ln sqrt(gamma))_,i 
= .5 (ln gamma)_,i 
= .5 gamma_,i / gamma
using gamma = gammaHat / W^6
= .5 (gammaHat W^-6)_,i / (gammaHat W^-6)
= .5 (gammaHat_,i W^-6 - 3 gammaHat W^-7) / (gammaHat W^-6)
= .5 gammaHat_,i / gammaHat - 3 / W
= GammaHat^j_ij - 3 / W
--]]
	vars:insert{
		name = 'expansion', 
		code = template([[
<?=eqn:makePartial'W'?>
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial'beta_U'?>
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	real exp_4phi = 1. / calc_exp_neg4phi(U);

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma_ll = sym3_real_mul(calc_gammaBar_ll(U, x), exp_4phi);

	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_real_mul(U->ABar_ll, exp_4phi),
		sym3_real_mul(gamma_ll, U->K/3.));

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	*value = -tr_partial_beta / U->alpha
<? 
for i,xi in ipairs(xNames) do
?>		+ U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		- U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		+ 3. * partial_W_l[<?=i-1?>] * U->beta_u.<?=xi?> / (U->W * U->alpha)
<?	for j,xj in ipairs(xNames) do
?>		+ K_ll.<?=sym(i,j)?> * U->beta_u.<?=xi?> * U->beta_u.<?=xj?> / (U->alpha * U->alpha)
<?	end
end
?>		- U->K;
]], env)
	}
--]=]		
--[=[
		--[[ ADM geodesic equation spatial terms:
		-Gamma^i_tt = 
			- gamma^ij alpha_,j

			+ alpha^-1 (
				gamma^ij beta^l gamma_kl beta^k_,j
				+ 1/2 gamma^ij gamma_kl,j beta^k beta^l
				- beta^i_,t
				- gamma^ij beta^k gamma_jk,t

				+ alpha^-1 beta^i (
					alpha_,t
					+ beta^j alpha_,j

					+ alpha^-1 (
						beta^i 1/2 beta^j beta^k gamma_jk,t
						- beta^i 1/2 beta^j beta^k beta^l gamma_kl,j
						- beta^i beta^j beta^l gamma_kl beta^k_,j
					)
				)
			)

		substitute 
		alpha_,t = -alpha^2 f K + beta^j alpha_,j
		beta^k_,t = B^k
		gamma_jk,t = -2 alpha K_jk + gamma_jk,l beta^l + gamma_lj beta^l_,k + gamma_lk beta^l_,j
		--]]
	vars:insert{
		name = 'gravity',
		code = template([[
<?=eqn:makePartial'alpha'?>

	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

<?=eqn:makePartial'beta_U'?>
<?=eqn:makePartial'epsilon_LL'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
<?=eqn:makePartial'W'?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_real_mul(partial_epsilon_lll[<?=i-1?>], _1_W * _1_W),
			sym3_real_mul(gammaBar_ll, 2. * partial_W_l[<?=i-1?>] * _1_W * _1_W * _1_W)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	sym3 dt_gamma_ll = sym3_zero;
	
	real partial_alpha_dot_beta = real3_dot(U->beta_u, *(real3*)partial_alpha_l);	//beta^j alpha_,j

	real3 beta_l = sym3_real3_mul(gamma_ll, U->beta_u);								//beta^j gamma_ij
	real3 beta_dt_gamma_l = sym3_real3_mul(dt_gamma_ll, U->beta_u);					//beta^j gamma_ij,t
	real beta_beta_dt_gamma = real3_dot(U->beta_u, beta_dt_gamma_l);				//beta^i beta^j gamma_ij,t
	
	real3 beta_dt_gamma_u = sym3_real3_mul(gamma_uu, beta_dt_gamma_l);				//gamma^ij gamma_jk,t beta^k

	//beta^i beta^j beta^k gamma_ij,k
	real beta_beta_beta_partial_gamma = 0.<?
for i,xi in ipairs(xNames) do
?> + U->beta_u.<?=xi?> * real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>)<?
end ?>;

	//beta_j beta^j_,i
	real3 beta_dbeta_l;
<? for i,xi in ipairs(xNames) do
?>	beta_dbeta_l.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
?>		+ beta_l.<?=xj?> * partial_beta_ul[<?=i-1?>].<?=xj?>
<?	end
?>	;
<? end
?>

	//beta_j beta^j_,i beta^i
	real beta_beta_dbeta = real3_dot(U->beta_u, beta_dbeta_l);

	//beta_j beta^j_,k gamma^ik
	real3 beta_dbeta_u = sym3_real3_mul(gamma_uu, beta_dbeta_l);

	//gamma_kl,j beta^k beta^l
	real3 beta_beta_dgamma_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>),
<? end
?>	};

	real3 beta_beta_dgamma_u = sym3_real3_mul(gamma_uu, beta_beta_dgamma_l);

<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> +=
		_1_alpha * (
			beta_dbeta_u.<?=xi?>
			+ .5 * beta_beta_dgamma_u.<?=xi?>	
			- U->B_U.<?=xi?>
			- beta_dt_gamma_u.<?=xi?>

			+ _1_alpha * U->beta_u.<?=xi?> * (
				.5 * dt_alpha
				+ partial_alpha_dot_beta

				+ _1_alpha * (
					.5 * beta_beta_dt_gamma
					- .5 * beta_beta_beta_partial_gamma 
					- beta_beta_dbeta
				)
			)
		)
	; 
<? end
?>
<? end	-- eqn.useShift ?>
]],			applyCommon{
				eqn = self,
				solver = self.solver,
			}
		), 
		type = 'real3',
	}
--]=]
-- [=[
	do
		-- should be zero when gammaBar_ij = gammaHat_ij
		vars:insert{
			name = 'RBar_LL',
			type = 'sym3',
			code = template([[
<?=eqn:makePartial'LambdaBar_U'?>

	//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
	real3x3 partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);

<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:getCode_connBar_ULL()?>

	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);
	
	_3sym3 partial_connHat_ulll[3];
	calc_partial_connHat_ulll(partial_connHat_ulll, x);

<?=eqn:getCode_connHat_LLL_and_ULL()?>
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	_3sym3 Delta_LLL = sym3_3sym3_mul(gammaBar_LL, Delta_ULL);

	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);

<?=eqn:makePartial2'epsilon_LL'?>
<?=eqn:getCode_trBar_partial2_gammaBar_LL()?>
<?=eqn:getCode_RBar_LL()?>
	*value_sym3 = sym3_rescaleToCoord_LL(RBar_LL, x);
]], env),
		}
	end
--]=]		
--[=[		
	vars:insert{
		name = 'RPhi_LL',
		type = 'sym3',
		code = template([[

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>

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
	
	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	sym3 partial2_phi_LL = sym3_rescaleFromCoord_ll(partial2_phi_ll, x);
	
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>
<?=eqn:getCode_connBar_ULL()?>

	sym3 DBar2_phi_LL = sym3_sub(
		partial2_phi_LL,
		real3_3sym3_dot1(
			partial_phi_L,
			connBar_ULL
		)
	);

	real tr_DBar2_phi = sym3_dot(gammaBar_UU, DBar2_phi_LL);

	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	sym3 RPhi_LL = {
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 
			- 2. * DBar2_phi_LL.<?=xij?>
			+ 4. * partial_phi_L.<?=xi?> * partial_phi_L.<?=xj?>
			+ gammaBar_LL.<?=xij?> * (
				- 2. * tr_DBar2_phi
				- 4. * real3_weightedLenSq(
					partial_phi_L,
					gammaBar_UU
				)
			),
<? end
?>	};

	*value_sym3 = RPhi_ll;
]], env),
	}

--[[
gammaBar_ij = gammaHat_ij + epsilon_ij
= gammaHat_ij + epsilon_IJ e_i^I e_j^J

gammaBar_ij,kl = gammaHat_ij,kl + (epsilon_IJ e_i^I e_j^J)_,kl
= gammaHat_ij,kl + (epsilon_IJ,k e_ij^IJ + epsilon_IJ e_ij^IJ_,k)_,l
= gammaHat_ij,kl
	+ epsilon_IJ,kl e_ij^IJ 
	+ epsilon_IJ,l e_ij^IJ_,k
	+ epsilon_IJ,k e_ij^IJ_,l
	+ epsilon_IJ e_ij^IJ_,kl

gammaBar^kl gammaBar_ij,kl

gammaBar^kl = inv(gammaBar_kl)
= inv(gammaHat_kl + epsilon_kl)
--]]
--]=]	

-- [=[
	vars:insert{
		name='trBar_partial2_gammaBar_LL',
		type = 'sym3',
		code = template([[
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>
<?=eqn:getCode_trBar_partial2_gammaBar_LL()?>
	*value_sym3 = trBar_partial2_gammaBar_LL;
]], env),
	}
--]=]

-- [=[
	vars:insert{
		name = 'tr34 (gamma*dGamma)',
		type = 'real3x3',
		code = template([[

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);

	_3sym3 partial_connHat_ulll[3];
	calc_partial_connHat_ulll(partial_connHat_ulll, x);

	real3x3 tr34_gamma_dGamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr34_gamma_dGamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * gammaBar_ll.<?=sym(i,m)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(j,k)?>
<?				end
			end
		end
?>	;
<?	end
end
?>
	
	*value_real3x3 = tr34_gamma_dGamma_ll;
]], env),
	}
--]=]

-- [=[
	vars:insert{
		name = 'tr14 (Gamma*dgamma)',
		type = 'real3x3',
		code = template([[
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);

<?=eqn:getCode_connHat_LLL_and_ULL()?>
	_3sym3 connHat_ull = _3sym3_rescaleToCoord_ULL(connHat_ULL, x); 

	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll = _3sym3_rescaleToCoord_LLL(partial_gammaBar_LLL, x);

	real3x3 tr14_Gamma_dgamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr14_Gamma_dgamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,i)?> * connHat_ull.<?=xm?>.<?=sym(k,j)?>
<?				end
			end
		end
?>	;
<?	end
end
?>
	*value_real3x3 = tr14_Gamma_dgamma_ll;
]], env),
	}
--]=]

--[=[
	for i,xi in ipairs(xNames) do
		vars:insert{
			name = 'Delta^'..xi,
			type = 'sym3',
			code = template([[
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:getCode_partial_gammaBar_LLL()?>
<?=eqn:getCode_connBar_ULL()?>

<?=eqn:getCode_connHat_LLL_and_ULL()?>

	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK + connHat^I_JK
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

	*value_sym3 = Delta_ULL.<?=xi?>;
]], table(env, {i=i, xi=xi})),
		}
	end
--]=]

	--[[ hmm? not working.
	vars:insert{name='x', type='real3', code='*value_real3=x;'}
	--]]
	-- [[
	vars:insert{name='x', code='*value=x.x;'}
	vars:insert{name='y', code='*value=x.y;'}
	vars:insert{name='z', code='*value=x.z;'}
	--]]

	return vars
end

function BSSNOKFiniteDifferenceEquation:fillRandom(epsilon)
	local ptr = BSSNOKFiniteDifferenceEquation.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.numCells-1 do
		ptr[i].alpha = ptr[i].alpha + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return BSSNOKFiniteDifferenceEquation
