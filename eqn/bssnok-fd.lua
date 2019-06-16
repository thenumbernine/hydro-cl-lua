--[[
Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer" 2010
Alcubierre "Introduction to Numerical Relativity" 2008

then I'm applying 2018 Ruchlin changes..
separate gammaBar_ll = gammaHat_ll + epsilon_ll
--]]

local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'eqn.einstein'
local makestruct = require 'eqn.makestruct'
local applyCommon = require 'common'

local makePartials = require 'eqn.makepartial'
local makePartial = makePartials.makePartial
local makePartial2 = makePartials.makePartial2

local BSSNOKFiniteDifferenceEquation = class(EinsteinEqn)
BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 
BSSNOKFiniteDifferenceEquation.hasEigenCode = true
BSSNOKFiniteDifferenceEquation.hasCalcDTCode = true
BSSNOKFiniteDifferenceEquation.hasFluxFromConsCode = true
BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true

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
		{name='beta_u', type='real3'},		 	-- 3: beta^i
		{name='epsilon_ll', type='sym3'},		-- 6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='W', type='real'},				-- 1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},				-- 1: K = K^i_i
		{name='ABar_ll', type='sym3'},			-- 6: ABar_ij, only 5 dof since ABar^k_k = 0
		{name='LambdaBar_u', type='real3'},		-- 3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
												-- TODO what is C^i ?
	}
	if self.useShift == 'HyperbolicGammaDriver' then
		intVars:insert{name='B_u', type='real3'}
	end

	self.consVars = table()
	:append(intVars)
	:append{
		--hyperbolic variables:
		--real3 a;			//3: a_i
		--_3sym3 dBar;		//18: dBar_ijk, only 15 dof since dBar_ij^j = 0
		--real3 Phi;			//3: Phi_i

		--stress-energy variables:
		{name='rho', type='real'},				--1: n_a n_b T^ab
		{name='S_u', type='real3'},			--3: -gamma^ij n_a T_aj
		{name='S_ll', type='sym3'},			--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},				--1
		{name='M_u', type='real3'},			--3
	}
	self.numIntStates = makestruct.countScalars(intVars)
	
	-- call construction / build structures	
	BSSNOKFiniteDifferenceEquation.super.init(self, args)
end


function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		--{name='constrain_det_gammaBar_ll', value=true, compileTime=true},
		{name='constrain_det_gammaBar_ll', value=false, compileTime=true},

		--{name='constrain_tr_ABar_ll', value=true, compileTime=true},
		{name='constrain_tr_ABar_ll', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},
	}
end

function BSSNOKFiniteDifferenceEquation:getTemplateEnv()
	local derivOrder = 2 * self.solver.numGhost
	return applyCommon{
		eqn = self,
		solver = self.solver,
		makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
		makePartial2 = function(...) return makePartial2(derivOrder, self.solver, ...) end,
	}
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template([[

//gammaBar_ij = gammaHat_ij + epsilon_ij
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_ll = coord_g_ll(x);
	return sym3_add(gammaHat_ll, U->epsilon_ll);
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2018 Ruchlin
real calc_det_gammaBar_ll(real3 x) {
	return coord_det_g(x);
}

void setFlatSpace(global <?=eqn.cons_t?>* U, real3 x) {
	U->alpha = 1.;
	U->beta_u = real3_zero;
	U->epsilon_ll = sym3_zero;
	U->W = 1;
	U->K = 0;
	U->ABar_ll = sym3_ident;
	U->LambdaBar_u = real3_zero;
<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_u = real3_zero;
<? end
?>

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_u = real3_zero;
}

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

//|g| = exp(12 phi) |g_grid|
real calc_det_gamma_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	real det_gamma_ll = coord_det_g(x) / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	return det_gamma_ll;
}

sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar_ll = calc_det_gammaBar_ll(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar_ll);
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
	real det_gamma_ll = coord_det_g(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma_ll); 
	return gamma_uu;
}

]], {eqn=self})
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	sym3 gammaHat_ll = coord_g_ll(x);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = gammaHat_ll;
	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

	<?=code?>

	U->alpha = alpha;
	U->beta_u = beta_u;

	real det_gamma_ll = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma_ll);
	
	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar_ll = calc_det_gammaBar_ll(x); 

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = cbrt(det_gammaBar_ll / det_gamma_ll);

	//W = exp(-2 phi)
	U->W = sqrt(exp_neg4phi);

	sym3 gammaBar_ll = sym3_real_mul(gamma_ll, exp_neg4phi);
	U->epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);

	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_real_mul(gamma_ll, 1./3. * U->K));
	U->ABar_ll = sym3_real_mul(A_ll, exp_neg4phi);
	
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	
	U->H = 0.;
	U->M_u = real3_zero;
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize LambdaBar_u
kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
<?=makePartial('epsilon_ll', 'sym3')?>

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	
	real det_gammaBar_ll = calc_det_gammaBar_ll(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar_ll);

	_3sym3 connHat_lll = coord_conn_lll(x);
	_3sym3 connHat_ull = coord_conn_ull(x);

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

	//Delta^i_jk = connBar^i_jk - connHat^i_jk
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);

<? for i,xi in ipairs(xNames) do
?>	U->LambdaBar_u.<?=xi?> = sym3_dot(Delta_ull.<?=xi?>, gammaBar_uu);
<? end
?>
}
]], table(self:getTemplateEnv(), {
		code = self.initState:initState(self.solver),
	}))
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd.cl'

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

BSSNOKFiniteDifferenceEquation.predefinedDisplayVars = {
	'U alpha',
	'U beta_u mag',
	--'U epsilon_ll norm',
	'U W',
	'U ABar_ll tr weighted',	-- has problems with spherical vacuum spacetime
	'U ABar_ll x x',
	'U ABar_ll x y',
	'U ABar_ll x z',
	'U ABar_ll y y',
	'U ABar_ll y z',
	'U ABar_ll z z',
	'U K',
	'U LambdaBar_u mag',
	'U H',
	'U M_u mag',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_ll tr weighted',

--[[ debugging derivatives
	'deriv alpha',
	'deriv beta_u mag',
	'deriv W',
	'deriv ABar_ll tr weighted',
	'deriv ABar_ll x x',
	'deriv ABar_ll x y',
	'deriv ABar_ll x z',
	'deriv ABar_ll y y',
	'deriv ABar_ll y z',
	'deriv ABar_ll z z',
	'deriv K',
	'deriv LambdaBar_u mag',
--]]
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:append{
		{gamma_ll = [[	*value_sym3 = calc_gamma_ll(U, x);]], type='sym3'},
		{gammaHat_ll = [[	*value_sym3 = coord_g_ll(x);]], type='sym3'},
		{gammaBar_ll = [[	*value_sym3 = calc_gammaBar_ll(U, x);]], type='sym3'},
		{gamma_uu = [[	*value_sym3 = calc_gamma_uu(U, x);]], type='sym3'},
		{gammaHat_uu = [[	*value_sym3 = coord_g_uu(x);]], type='sym3'},
		{gammaBar_uu = [[	*value_sym3 = calc_gammaBar_uu(U, x);]], type='sym3'},
		{K_ll = [[
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	*value_sym3 = sym3_real_mul(
		sym3_add(
			U->ABar_ll,
			sym3_real_mul(gammaBar_ll, U->K / 3.)
		), exp_4phi);
]], type='sym3'},

		{['det gammaBar - det gammaHat'] = [[
	*value = sym3_det(calc_gammaBar_ll(U, x)) - calc_det_gammaBar_ll(x);
]]},
		{['det gamma based on phi'] = [[
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	*value = exp_12phi * calc_det_gammaBar_ll(x);
]]},
		{S = '*value = sym3_dot(U->S_ll, calc_gamma_uu(U, x));'},
		{volume = '*value = U->alpha * calc_det_gamma_ll(U, x);'},
	}

	local derivOrder = 2 * self.solver.numGhost
	vars:append{
	
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
		{expansion = template([[
	<?=makePartial('W', 'real')?>
	<?=makePartial('alpha', 'real')?>
	<?=makePartial('beta_u', 'real3')?>
	real tr_partial_beta = 0. <?
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
]], 			applyCommon{
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
				}

			)
		},
		
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},

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
		{
			gravity = template([[
	<?=makePartial('alpha', 'real')?>

	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

	<?=makePartial('beta_u', 'real3')?>

	<?=makePartial('epsilon_ll', 'sym3')?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
	<?=makePartial('W', 'real')?>
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
	real3 beta_dbeta_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3_dot(beta_l, partial_beta_ul[<?=i-1?>]),
<? end
?>	};

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
			- U->B_u.<?=xi?>
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
]],				applyCommon{
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
				}
			), 
			type = 'real3',
		},
--[=[	
		{
			Ricci = template([[
	_3sym3 connHat_ull = coord_conn_ull(x);
	
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);
	real3 Delta_u = _3sym3_sym3_dot23(Delta_ull, gammaBar_uu);
	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);
	
	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	sym3 gammaHat_uu = coord_g_uu(x);
	
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	_3sym3 partial_gammaBar_lll = _3sym3_add(*(_3sym3*)partial_epsilon_lll, partial_gammaHat_lll);
	
	_3sym3 partial_connHat_ulll[3];
	calc_partial_conn_ulll(partial_connHat_ulll, connHat_ull, gammaHat_uu, partial_gammaHat_lll, partial2_gammaHat_llll);
	
	sym3sym3 partial2_gammaBar_llll = sym3sym3_add(*(sym3sym3*)partial2_epsilon_llll, partial2_gammaHat_llll);
	sym3 RBar_ll = calc_RBar_ll(U, gammaBar_ll, gammaBar_uu, connHat_ull, partial_connHat_ulll, partial2_gammaBar_llll, partial_gammaBar_lll, partial_LambdaBar_ul, Delta_u, Delta_lll, Delta_ull);
	
	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	sym3 RPhi_ll = {
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
	*value_sym3 = R_ll;
]],				applyCommon{
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
				}
			),
			type = 'sym3',
		},
--]=]	
	}
	
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
