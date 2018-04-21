local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'eqn.einstein'
local makestruct = require 'eqn.makestruct'
require 'common'(_G)

local makePartials = require 'eqn.makepartial'
local makePartial = makePartials.makePartial
local makePartial2 = makePartials.makePartial2

local Z4cFiniteDifferenceEquation = class(EinsteinEqn)
Z4cFiniteDifferenceEquation.name = 'Z4c finite difference' 
Z4cFiniteDifferenceEquation.hasEigenCode = true
Z4cFiniteDifferenceEquation.useConstrainU = true
Z4cFiniteDifferenceEquation.useSourceTerm = true

--[[
args:
	useHypGammaDriver
--]]
function Z4cFiniteDifferenceEquation:init(args)
	-- options:

	-- needs to be defined up front
	-- otherwise rebuild intVars based on it ...
	if args.useHypGammaDriver ~= nil then
		self.useHypGammaDriver = args.useHypGammaDriver
	else
		self.useHypGammaDriver = false
	end

	local intVars = table{
		{alpha = 'real'},			-- 1
		{beta_u = 'real3'},         -- 3: beta^i
		{gammaTilde_ll = 'sym3'},		-- 6: gammaTilde_ij, only 5 dof since det gammaTilde_ij = 1
		{chi = 'real'},				-- 1
		{KHat = 'real'},			-- 1
		{Theta = 'real'},			-- 1
		{ATilde_ll = 'sym3'},       -- 6: ATilde_ij, only 5 dof since ATilde^k_k = 0
		{connBar_u = 'real3'},      -- 3: connBar^i = gammaTilde^jk connBar^i_jk = -partial_j gammaTilde^ij
	}

	if self.useHypGammaDriver then
		intVars:insert{B_u = 'real3'}
	end

	self.consVars = table()
	:append(intVars)
	:append{
		--hyperbolic variables:
		--real3 a;			//3: a_i
		--_3sym3 dTilde;		//18: dTilde_ijk, only 15 dof since dTilde_ij^j = 0
		--real3 Phi;			//3: Phi_i

		--stress-energy variables:
		{rho = 'real'},		--1: n_a n_b T^ab
		{S_u = 'real3'},			--3: -gamma^ij n_a T_aj
		{S_ll = 'sym3'},			--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{H = 'real'},				--1
		{M_u = 'real3'},			--3

		-- aux variable
		{gammaTilde_uu = 'sym3'},		--6
	}
	self.numIntStates = makestruct.countReals(intVars)

	-- call construction / build structures	
	Z4cFiniteDifferenceEquation.super.init(self, args)
end

function Z4cFiniteDifferenceEquation:createInitState()
	Z4cFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaTilde_ll', value=true},
		{name='constrain_tr_ATilde_ll', value=true},
		{name='useGammaDriver', value=true},
		{name='diffuseSigma', value=.01},
	}
end

function Z4cFiniteDifferenceEquation:getTemplateEnv()
	local derivOrder = 2 * self.solver.numGhost
	return {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from3x3to6 = from3x3to6,
		from6to3x3 = from6to3x3,
		sym = sym,
		makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
		makePartial2 = function(...) return makePartial2(derivOrder, self.solver, ...) end,
	}
end

function Z4cFiniteDifferenceEquation:getCommonFuncCode()
	return template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1.;
	U->beta_u = _real3(0,0,0);
	U->gammaTilde_ll = _sym3(1,0,0,1,0,1);
	U->chi = 1;
	U->KHat = 0;
	U->Theta = 0;
	U->ATilde_ll = _sym3(1,0,0,1,0,1);
	U->connBar_u = _real3(0,0,0);
<? if eqn.useHypGammaDriver then
?>	U->B_u = _real3(0,0,0);
<? end
?>	U->gammaTilde_uu = _sym3(1,0,0,1,0,1);

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = _real3(0,0,0);
	U->S_ll = _sym3(0,0,0,0,0,0);
	U->H = 0;
	U->M_u = _real3(0,0,0);
}

#define calc_exp_neg4phi(U) ((U)->chi)

//|g| = exp(12 phi) |g_grid|
real calc_det_gamma(global const <?=eqn.cons_t?>* U, real3 x) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	real det_gamma = sqrt_det_g_grid(x) / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	return det_gamma;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	sym3 gamma_uu = sym3_scale(U->gammaTilde_uu, exp_neg4phi);
	return gamma_uu;
}

]], {eqn=self})
end

function Z4cFiniteDifferenceEquation:getInitStateCode()
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = _real3(0,0,0);
	sym3 gamma_ll = _sym3(1,0,0,1,0,1);
	sym3 K_ll = _sym3(0,0,0,0,0,0);
	real rho = 0.;

	<?=code?>

	U->alpha = alpha;
	U->beta_u = beta_u;

	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);

	//gammaTilde_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real det_gammaGrid = sqrt_det_g_grid(x);
	real exp_neg4phi = cbrt(det_gammaGrid / det_gamma);

	U->chi = exp_neg4phi;
	
	U->gammaTilde_ll = sym3_scale(gamma_ll, exp_neg4phi);
	U->gammaTilde_uu = sym3_inv(U->gammaTilde_ll, 1.);

]]--[[
<? for _,x in ipairs(xNames) do
?>	U->a.<?=x?> = calc_a_<?=x?>(x.x, x.y, x.z);
<? end
?>	
]]..[[	

	U->Theta = 0.;	//TODO ... Theta = -Z^mu n_mu = alpha * Z^t ... which is?

	real K = sym3_dot(K_ll, gamma_uu);
	U->KHat = K - 2. * U->Theta;
	sym3 A_ll = sym3_sub(K_ll, sym3_scale(gamma_ll, 1./3. * K));
	U->ATilde_ll = sym3_scale(A_ll, exp_neg4phi);
	
	U->rho = rho;
	U->S_u = _real3(0,0,0);
	U->S_ll = _sym3(0,0,0,0,0,0);
	
	U->H = 0.;
	U->M_u = _real3(0,0,0);
}

//after popularing gammaTilde_ll, use its finite-difference derivative to initialize connBar_u
kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
<?=makePartial('gammaTilde_uu', 'sym3')?>

	//connBar^i = -gammaTilde^ij_,j + 2 gammaTilde^ij Z_j
	//...= gammaTilde_jk,l gammaTilde^ij gammaTilde^kl + 2 gammaTilde^ij Z_j
<? for i,xi in ipairs(xNames) do
?>	U->connBar_u.<?=xi?> =<?
	for j,xj in ipairs(xNames) do
?> - partial_gammaTilde_uul[<?=j-1?>].<?=sym(i,j)?><?
	end
?>;
<? end
?>
}
]], table(self:getTemplateEnv(), {
	code = self.initState.initState 
		and self.initState:initState(self.solver) 
		or '//no code from InitCond:initState() was provided',
}))
end

function Z4cFiniteDifferenceEquation:getSolverCode()
	return template(file['eqn/z4c-fd.cl'], self:getTemplateEnv())
end

function Z4cFiniteDifferenceEquation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
]], {
		eqn = self,
	})
end

function Z4cFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

function Z4cFiniteDifferenceEquation:getDisplayVars()	
	local vars = Z4cFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:insert{['det gammaTilde-1'] = [[*value = -1. + sym3_det(U->gammaTilde_ll);]]}	-- for logarithmic displays
	vars:insert{['det gamma based on phi'] = [[
	real exp_neg4phi = calc_exp_neg4phi(U);
	*value = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);   
]]}
	
	local derivOrder = 2 * self.solver.numGhost
	vars:append{
		{S = '*value = sym3_dot(U->S_ll, calc_gamma_uu(U));'},
		{volume = '*value = U->alpha * calc_det_gamma(U, x);'},
	
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
		{expansion = template([[
	<?=makePartial('chi', 'real')?>
	<?=makePartial('alpha', 'real')?>
	<?=makePartial('beta_u', 'real3')?>
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	real exp_4phi = 1. / calc_exp_neg4phi(U);

	//gamma_ij = exp(4 phi) gammaTilde_ij
	sym3 gamma_ll = sym3_scale(U->gammaTilde_ll, exp_4phi);

	//K = KHat + 2 Theta
	real K = U->KHat + 2. * U->Theta;

	//K_ij = exp(4 phi) ATilde_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_scale(U->ATilde_ll, exp_4phi),
		sym3_scale(gamma_ll, K/3.));

	*value = -tr_partial_beta / U->alpha
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
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
					xNames = xNames,
					sym = sym,
				}

			)
		},
		
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
		{gamma_x = '*valuevec = real3_scale(sym3_x(U->gammaTilde_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{gamma_y = '*valuevec = real3_scale(sym3_y(U->gammaTilde_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{gamma_z = '*valuevec = real3_scale(sym3_z(U->gammaTilde_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{K_x = '*valuevec = real3_add(sym3_x(U->ATilde_ll), real3_scale(sym3_x(U->gammaTilde_ll), (U->KHat + 2. * U->Theta)/3.));', type='real3'},
		{K_y = '*valuevec = real3_add(sym3_y(U->ATilde_ll), real3_scale(sym3_y(U->gammaTilde_ll), (U->KHat + 2. * U->Theta)/3.));', type='real3'},
		{K_z = '*valuevec = real3_add(sym3_z(U->ATilde_ll), real3_scale(sym3_z(U->gammaTilde_ll), (U->KHat + 2. * U->Theta)/3.));', type='real3'},

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
	<?=makePartial('beta_u', 'real3')?>

	<?=makePartial('gammaTilde_ll', 'sym3')?>
	
	//chi = exp(-4 phi)
	real _1_chi = 1. / U->chi;
	
	//gamma_ij = 1/chi gammaTilde_ij
	sym3 gamma_ll = sym3_scale(U->gammaTilde_ll, _1_chi);
	
	//gamma_ij,k = 1/chi gammaTilde_ij,k - chi,k / chi^2 gammaTilde_ij
	<?=makePartial('chi', 'real')?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_scale(partial_gammaTilde_lll[<?=i-1?>], _1_chi),
			sym3_scale(U->gammaTilde_ll, partial_chi_l[<?=i-1?>] * _1_chi * _1_chi)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	sym3 dt_gamma_ll = _sym3(0,0,0,0,0,0);


	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U);
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
?>	valuevec->s<?=i-1?> = -partial_alpha_u.<?=xi?>

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
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
					xNames = xNames,
					sym = sym,
				}
			), 
			type = 'real3',
		},
	}
	
	return vars
end

function Z4cFiniteDifferenceEquation:fillRandom(epsilon)
	local ptr = Z4cFiniteDifferenceEquation.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gammaTilde_ll.xx = ptr[i].gammaTilde_ll.xx + 1
		ptr[i].gammaTilde_ll.yy = ptr[i].gammaTilde_ll.yy + 1
		ptr[i].gammaTilde_ll.zz = ptr[i].gammaTilde_ll.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

function Z4cFiniteDifferenceEquation:getCalcDTCode() end
function Z4cFiniteDifferenceEquation:getFluxFromConsCode() end

return Z4cFiniteDifferenceEquation
