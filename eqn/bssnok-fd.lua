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


local BSSNOKFiniteDifferenceEquation = class(EinsteinEqn)

BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 

BSSNOKFiniteDifferenceEquation.hasEigenCode = true

-- options:

-- needs to be defined up front
-- otherwise rebuild intVars based on it ...
BSSNOKFiniteDifferenceEquation.useHypGammaDriver = true

-- use chi = 1/psi instead of phi, as described in 2006 Campanelli 
-- it should be used with the hyperbolic gamma driver
BSSNOKFiniteDifferenceEquation.useChi = true

local intVars = table{
	{alpha = 'real'},			-- 1
	{beta_u = 'real3'},         -- 3: beta^i
	{gammaBar_ll = 'sym3'},    -- 6: gammaBar_ij, only 5 dof since det gammaBar_ij = 1
                                                                                                 
	BSSNOKFiniteDifferenceEquation.useChi 
		and {chi = 'real'}		-- 1
		or {phi = 'real'},		-- 1
	
	{K = 'real'},               -- 1
	{ATilde_ll = 'sym3'},       -- 6: ATilde_ij, only 5 dof since ATilde^k_k = 0
	{connBar_u = 'real3'},      -- 3: connBar^i = gammaBar^jk connBar^i_jk = -partial_j gammaBar^ij
}

if BSSNOKFiniteDifferenceEquation.useHypGammaDriver then
	intVars:insert{B_u = 'real3'}
end

local consVars = table()
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
	{gammaBar_uu = 'sym3'},		--6
}

BSSNOKFiniteDifferenceEquation.consVars = consVars
BSSNOKFiniteDifferenceEquation.numIntStates = makestruct.countReals(intVars)

BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true


function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaBar_ll', value=true},
		{name='constrain_tr_ATilde_ll', value=true},
		{name='useGammaDriver', value=false},
		{name='diffuseSigma', value=.01},
	}
end

function BSSNOKFiniteDifferenceEquation:getTemplateEnv()
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

-- this is separated so both BSSNOK:getCodePrefix and GRHDSeparate's HydroSolver:createCodePrefix can get to it
function BSSNOKFiniteDifferenceEquation:getExtraCLFuncs()
	return template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1.;
	U->beta_u = _real3(0,0,0);
	U->gammaBar_ll = _sym3(1,0,0,1,0,1);
<? if eqn.useChi then 
?>	U->chi = 1;
<? else
?>	U->phi = 0;
<? end
?>	U->K = 0;
	U->ATilde_ll = _sym3(1,0,0,1,0,1);
	U->connBar_u = _real3(0,0,0);
<? if eqn.useHypGammaDriver then
?>	U->B_u = _real3(0,0,0);
<? end
?>	U->gammaBar_uu = _sym3(1,0,0,1,0,1);

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = _real3(0,0,0);
	U->S_ll = _sym3(0,0,0,0,0,0);
	U->H = 0;
	U->M_u = _real3(0,0,0);
}

<? if eqn.useChi then
?>#define calc_exp_neg4phi(U) ((U)->chi)
<? else
?>#define calc_exp_neg4phi(U) (exp(-4. * (U)->phi))
<? end
?>

real calc_det_gamma(global const <?=eqn.cons_t?>* U) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	real det_gamma = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	return det_gamma;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	sym3 gamma_uu = sym3_scale(U->gammaBar_uu, exp_neg4phi);
	return gamma_uu;
}

]], {eqn=self})
end

-- should this be getInitStateCode like in eqn/euler?
function BSSNOKFiniteDifferenceEquation:getCodePrefix()
	local lines = table()
	
	lines:insert(BSSNOKFiniteDifferenceEquation.super.getCodePrefix(self))
	
	lines:insert(self:getExtraCLFuncs())
	
	return lines:concat()
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
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

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = 1./cbrt(det_gamma);

<? if eqn.useChi then 
?>	U->chi = exp_neg4phi;
<? else
?>	U->phi = log(det_gamma) / 12.;
<? end 
?>
	U->gammaBar_ll = sym3_scale(gamma_ll, exp_neg4phi);
	U->gammaBar_uu = sym3_inv(U->gammaBar_ll, 1.);

]]--[[
<? for _,x in ipairs(xNames) do
?>	U->a.<?=x?> = calc_a_<?=x?>(x.x, x.y, x.z);
<? end
?>	
]]..[[	

	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_scale(gamma_ll, 1./3. * U->K));
	U->ATilde_ll = sym3_scale(A_ll, exp_neg4phi);
	
	U->rho = rho;
	U->S_u = _real3(0,0,0);
	U->S_ll = _sym3(0,0,0,0,0,0);
	
	U->H = 0.;
	U->M_u = _real3(0,0,0);
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize connBar_u
kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	const global <?=eqn.cons_t?>* Up[dim];
	const global <?=eqn.cons_t?>* Um[dim];
	<? for j=0,solver.dim-1 do ?>{
		Up[<?=j?>] = U + stepsize.s<?=j?>;
		Um[<?=j?>] = U - stepsize.s<?=j?>;
	}<? end ?>

<?=makePartial('gammaBar_uu', 'sym3')?>

	//connBar^i = -gammaBar^ij_,j
<? for i,xi in ipairs(xNames) do
?>	U->connBar_u.<?=xi?> =<?
	for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
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

function BSSNOKFiniteDifferenceEquation:getSolverCode()
	return template(file['eqn/bssnok-fd.cl'], self:getTemplateEnv())
end

function BSSNOKFiniteDifferenceEquation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
]], {
		eqn = self,
	})
end

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:insert{['det gammaBar-1'] = [[*value = -1. + sym3_det(U->gammaBar_ll);]]}	-- for logarithmic displays
	vars:insert{['det gamma based on phi'] = [[
	real exp_neg4phi = calc_exp_neg4phi(U);
	*value = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);   
]]}
	
	local derivOrder = 2 * self.solver.numGhost
	vars:append{
		{S = '*value = sym3_dot(U->S_ll, calc_gamma_uu(U));'},
		{volume = '*value = U->alpha * calc_det_gamma(U);'},
		{expansion = '*value = -U->alpha * U->K;'},
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
		{gamma_x = '*valuevec = real3_scale(sym3_x(U->gammaBar_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{gamma_y = '*valuevec = real3_scale(sym3_y(U->gammaBar_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{gamma_z = '*valuevec = real3_scale(sym3_z(U->gammaBar_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{K_x = '*valuevec = real3_add(sym3_x(U->ATilde_ll), real3_scale(sym3_x(U->gammaBar_ll), U->K/3.));', type='real3'},
		{K_y = '*valuevec = real3_add(sym3_y(U->ATilde_ll), real3_scale(sym3_y(U->gammaBar_ll), U->K/3.));', type='real3'},
		{K_z = '*valuevec = real3_add(sym3_z(U->ATilde_ll), real3_scale(sym3_z(U->gammaBar_ll), U->K/3.));', type='real3'},

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

	<?=makePartial('gammaBar_ll', 'sym3')?>
<? if not eqn.useChi then ?>

	//gamma_ij = exp(4 phi) gammaBar_ij
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_scale(exp_4phi, U->gammaBar_ll);

	//gamma_ij,k = exp(4 phi) gammaBar_ij,k + 4 phi,k exp(4 phi) gammaBar_ij
	<?=makePartial('phi', 'real')?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_add(
			sym3_scale(partial_gammaBar_lll[<?=i-1?>], exp_4phi),
			sym3_scale(U->gammaBar_ll, 4. * exp_4phi * partial_phi_l[<?=i-1?>])),
<? end
?>	};

<? else -- useChi ?>

	//chi = exp(-4 phi)
	real _1_chi = 1. / U->chi;
	
	//gamma_ij = 1/chi gammaBar_ij
	sym3 gamma_ll = sym3_scale(U->gammaBar_ll, _1_chi);
	
	//gamma_ij,k = 1/chi gammaBar_ij,k - chi,k / chi^2 gammaBar_ij
	<?=makePartial('chi', 'real')?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_scale(partial_gammaBar_lll[<?=i-1?>], _1_chi),
			sym3_scale(U->gammaBar_ll, partial_chi_l[<?=i-1?>] * _1_chi * _1_chi)),
<? end
?>	};

<? end -- useChi ?>

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

function BSSNOKFiniteDifferenceEquation:fillRandom(epsilon)
	local ptr = BSSNOKFiniteDifferenceEquation.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gammaBar_ll.xx = ptr[i].gammaBar_ll.xx + 1
		ptr[i].gammaBar_ll.yy = ptr[i].gammaBar_ll.yy + 1
		ptr[i].gammaBar_ll.zz = ptr[i].gammaBar_ll.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return BSSNOKFiniteDifferenceEquation
