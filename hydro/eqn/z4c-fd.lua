--[[
Alcubierre "Introduction to 3+1 Numerical Relativity"
Baumgarte & Shapiro "Numerical Relativity"
2005 Campanelli
2011 Cao, Hilditch "Numerical stability of the Z4c formulation of general relativity"
2017 Ruchlin

changes I'm making to coincide with 2017 Ruchlin 
1) rename the gammaTilde_ll => gammaBar_ll
2) rename ATilde_ll => ABar_ll
3) separate gammaBar_ll = gammaHat_ll + epsilon_ll

TODO implement these in hydro/eqn/ and hydro/solver/ bssnok-fd.lua
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'
local Struct = require 'hydro.code.struct'

local makePartials = require 'hydro.eqn.makepartial'
local makePartial1 = makePartials.makePartial1
local makePartial2 = makePartials.makePartial2

local Z4cFiniteDifferenceEquation = class(EinsteinEqn)
Z4cFiniteDifferenceEquation.name = 'Z4c finite difference' 

Z4cFiniteDifferenceEquation.weightFluxByGridVolume = false

--[[
args:
	useShift:
		none
		GammaDriver
		HyperbolicGammaDriver
--]]
function Z4cFiniteDifferenceEquation:init(args)	
	self.useShift = args.useShift or 'none'

	local intVars = table{
		{name='alpha', type='real'},			-- 1
		{name='beta_u', type='real3'},         -- 3: beta^i
		{name='epsilon_ll', type='sym3'},		-- 6: epsilon_ij = gammaBar_ij - gammaHat_ij, where gammaHat_ij = grid metric. This has only 5 dof since det gammaBar_ij = 1
		{name='chi', type='real'},				-- 1
		{name='KHat', type='real'},			-- 1
		{name='Theta', type='real'},			-- 1
		{name='ABar_ll', type='sym3'},       	-- 6: ABar_ij, only 5 dof since ABar^k_k = 0
		{name='Delta_u', type='real3'},      	-- 3: Delta^i = gammaBar^jk Delta^i_jk
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
		{name='rho', type='real'},			--1: n_a n_b T^ab
		{name='S_u', type='real3'},			--3: -gamma^ij n_a T_aj
		{name='S_ll', type='sym3'},			--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},			--1
		{name='M_u', type='real3'},			--3

		-- aux variable
		{name='gammaBar_uu', type='sym3'},	--6
	}
	self.numIntStates = Struct.countScalars{vars=intVars}

	-- call construction / build structures	
	Z4cFiniteDifferenceEquation.super.init(self, args)
end

function Z4cFiniteDifferenceEquation:createInitState()
	Z4cFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaBar_ll', value=true, compileTime=true},
		{name='constrain_tr_ABar_ll', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		{name='alphaMin', value=1e-7},
	}
end

-- don't use default
function Z4cFiniteDifferenceEquation:initCodeModule_calcDTCell() end

Z4cFiniteDifferenceEquation.solverCodeFile = 'hydro/eqn/z4c-fd.cl'

Z4cFiniteDifferenceEquation.predefinedDisplayVars = {
	'U alpha',
	'U beta_u mag',
	'U epsilon_ll norm',
	'U chi',
	'U KHat',
	'U Theta',
	'U ABar_ll norm',
	'U ABar_ll tr weighted',
	'U Delta_u mag',
	'U B_u mag',
	'U H',
	--'U M_u mag',	-- not implemented yet ...
	'U det gammaBar - det gammaHat',
	'U det gamma_ij based on phi',
	'U volume',
	'U f',
}

function Z4cFiniteDifferenceEquation:getDisplayVars()	
	local vars = Z4cFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:insert{name='det gammaBar - det gammaHat', code = self:template[[
	value.vreal = sym3_det(<?=calc_gammaBar_ll?>(U, x)) - calc_det_gammaBar_ll(x);
]]}	-- for logarithmic displays
	vars:insert{name='det gamma_ij based on phi', code = self:template[[
	real exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	value.vreal = calc_det_gammaBar_ll(x) / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
]]}
	
	local derivOrder = 2 * self.solver.numGhost
	vars:append{
		{name='S', code='value.vreal = sym3_dot(U->S_ll, <?=calc_gamma_uu?>(U, x));'},
		{name='volume', code='value.vreal = U->alpha * calc_det_gamma_ll(U, x);'},
	
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

Gamma^j_ij = (ln sqrt(g))_,i = .5 (ln g)_,i = .5 g_,i / g

(det g)_,i / (det g)
... using phi ...
=  exp(12 phi)_,i / exp(12 phi)
= 12 exp(12 phi) phi_,i / exp(12 phi)
= 12 phi_,i
... using chi ...
= (chi^-3)_,i / (chi^-3)
= -3 chi_,i / chi^4 / (chi^-3)
= -3 chi_,i / chi
--]]
		{name='expansion', code = self:template([[
	<?=makePartial1('chi', 'real')?>
	<?=makePartial1('alpha', 'real')?>
	<?=makePartial1('beta_u', 'real3')?>
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma_ll = sym3_real_mul(<?=calc_gammaBar_ll?>(U, x), exp_4phi);

	//K = KHat + 2 Theta
	real K = U->KHat + 2. * U->Theta;

	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_real_mul(U->ABar_ll, exp_4phi),
		sym3_real_mul(gamma_ll, K/3.));

	value.vreal = -tr_partial_beta / U->alpha
<? 
for i,xi in ipairs(xNames) do
?>		+ U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		- U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		+ 1.5 * partial_chi_l[<?=i-1?>] / U->chi * U->beta_u.<?=xi?> / U->alpha
<?	for j,xj in ipairs(xNames) do
?>		+ K_ll.<?=sym(i,j)?> * U->beta_u.<?=xi?> * U->beta_u.<?=xj?> / (U->alpha * U->alpha)
<?	end
end
?>		- K;
]], 			{
					makePartial1 = function(...) return makePartial1(derivOrder, self.solver, ...) end,
				}
			)
		},
		
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		{name='gamma_ll', code = self:template[[
	{
		real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);
		sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
		value.vsym3 = sym3_real_mul(gammaBar_ll, exp_4phi);
	}
]], type='sym3'},
	
		-- K_ij = exp(4 phi) ABar_ij + K/3 gamma_ij  
		-- gamma_ij = exp(4 phi) gammaBar_ij
		-- K_ij = exp(4 phi) (ABar_ij + K/3 gammaBar_ij)
		{name='K_ll', code = self:template[[
	real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);
	real K = U->KHat + 2. * U->Theta;
	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	value.vsym3 = sym3_real_mul(
		sym3_add(
			U->ABar_ll,
			sym3_real_mul(gammaBar_ll, K / 3.)
		), exp_4phi);
]], type='sym3'},

--[=[ TODO FIXME
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
			name = 'gravity',
			code = self:template([[
	<?=makePartial1('alpha', 'real')?>
	<?=makePartial1('beta_u', 'real3')?>

	//gammaBar_ij = gammaHat_ij + epsilon_ij
	//gammaBar_ij,k = epsilon_ij,k for static meshes
	<?=makePartial1('epsilon_ll', 'sym3')?>

	//chi = exp(-4 phi)
	real _1_chi = 1. / U->chi;
	
	//gamma_ij = 1/chi gammaBar_ij
	sym3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_chi);
	
	//gamma_ij,k = 1/chi gammaBar_ij,k - chi,k / chi^2 gammaBar_ij
	<?=makePartial1('chi', 'real')?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_real_mul(partial_epsilon_lll[<?=i-1?>], _1_chi),
			sym3_real_mul(gammaBar_ll, partial_chi_l[<?=i-1?>] * _1_chi * _1_chi)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	sym3 dt_gamma_ll = sym3_zero;


	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = <?=calc_gamma_uu?>(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
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
?>	value_real3->s<?=i-1?> = -partial_alpha_u.<?=xi?>

		+ _1_alpha * (
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
]],				{
					makePartial1 = function(...) return makePartial1(derivOrder, self.solver, ...) end,
				}
			), 
			type = 'real3',
		},
--]=]	
	}
	
	return vars
end

return Z4cFiniteDifferenceEquation
