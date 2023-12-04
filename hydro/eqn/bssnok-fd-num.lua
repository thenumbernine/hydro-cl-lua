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
local table = require 'ext.table'
local BSSNOKFiniteDifferenceEquationBase = require 'hydro.eqn.bssnok-fd'
local Struct = require 'hydro.code.struct'

local BSSNOKFiniteDifferenceEquation = BSSNOKFiniteDifferenceEquationBase:subclass()

BSSNOKFiniteDifferenceEquation.name = 'bssnok_fd_num'

-- not used with finite-difference schemes anyways
BSSNOKFiniteDifferenceEquation.weightFluxByGridVolume = false

BSSNOKFiniteDifferenceEquation.useScalarField = false

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
	self.useScalarField = args.useScalarField

	local intVars = table{
		{name='alpha', type='real'},						-- 0:	1: alpha
		{name='W', type='real'},							-- 1:	1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},							-- 2:	1: K = K^i_i
		{name='beta_U', type='real3', variance='u'},		-- 3:	3: beta^i
		{name='B_U', type='real3', variance='u'},			-- 6:	3: B^i ... only used with HyperbolicGammaDriver
		{name='LambdaBar_U', type='real3', variance='u'},	-- 9:	3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
		{name='epsilon_LL', type='real3s3', variance='ll'},	-- 12:	6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='ABar_LL', type='real3s3', variance='ll'},		-- 18:	6: ABar_ij, only 5 dof since ABar^k_k = 0
	}														-- 24 = total size so far

	if self.useScalarField then
		intVars:append{
			{name='Phi', type='cplx'},
			{name='Psi_l', type='cplx3'},					-- Psi_i = Phi_,i
			{name='Pi', type='cplx'},						-- Pi = n^a Phi_,a
		}
	end

	self.consVars = table()
	:append(intVars)
	:append{
		--stress-energy variables:
		{name='rho', type='real'},							-- 1: n_a n_b T^ab
		{name='S_u', type='real3'},							-- 3: -gamma^ij n_a T_aj
		{name='S_ll', type='real3s3'},							-- 6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},							-- 1
		{name='M_U', type='real3', variance='u'},			-- 3
	}
	self:cdefAllVarTypes(args.solver, self.consVars)	-- have to call before countScalars in eqn:init
	self.numIntStates = Struct.countScalars{vars=intVars}

	-- call construction / build structures
	BSSNOKFiniteDifferenceEquation.super.init(self, args)
end

function BSSNOKFiniteDifferenceEquation:getSymbolFields()
	return BSSNOKFiniteDifferenceEquation.super.getSymbolFields(self):append{
		'calc_partial_det_gammaHat_l',
	}
end

function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		--{name='constrain_det_gammaBar', value=true, compileTime=true},
		{name='constrain_det_gammaBar', value=false, compileTime=true},

		--{name='constrain_tr_ABar', value=true, compileTime=true},
		{name='constrain_tr_ABar', value=false, compileTime=true},

		{name='calc_H_and_M', value=true, compileTime=true},
		--{name='calc_H_and_M', value=false, compileTime=true},

		-- 2013 Baumgarte et al, section IIIB
		--{name='dissipationCoeff', value=.001/16},
		{name='dissipationCoeff', value=cmdline.dissipationCoeff or .99},	-- value matching SENR.  runs until 12.09.
		--{name='dissipationCoeff', value=0},		-- runs until 12.2275.  This is working better than 0.99.  hmm...
		--{name='dissipationCoeff', value=-.99},	-- nope, this is worse, i must have my sign right

		{name='alphaMin', value=1e-3},
		--{name='alphaMin', value=1e-7},
		--{name='alphaMin', value=-math.huge},

		{name='dt_beta_U_k', value=cmdline.dt_beta_U_k or 3/4},
		{name='dt_beta_U_eta', value=cmdline.dt_beta_U_eta or 2},	--1, or 1 / (2 M), for total mass M, or SENR uses 2
	}

	if self.useScalarField then
		self:addGuiVars{
			{name='scalar_lambda', value=0},
			{name='scalar_mu', value=1},
		}
	end
end

function BSSNOKFiniteDifferenceEquation:compile(expr)
	return self.solver.coord:compile(expr)
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'hydro/eqn/bssnok-fd-num.cl'

BSSNOKFiniteDifferenceEquation.predefinedDisplayVars = {
-- [=[
	'U alpha',
--	'U beta_U mag',
--	'U beta_U x',
--	'U beta_U y',
--	'U beta_U z',
--	'U B_U mag',
--	'U B_U x',
--	'U B_U y',
--	'U B_U z',
	--'U epsilon_LL norm',
--	'U epsilon_LL xx',
--	'U epsilon_LL xy',
--	'U epsilon_LL xz',
--	'U epsilon_LL yy',
--	'U epsilon_LL yz',
--	'U epsilon_LL zz',
	'U W',
	'U K',
--	'U ABar_LL tr weighted',
--] =]
--	'U ABar_LL xx',
--	'U ABar_LL xy',
--	'U ABar_LL xz',
--	'U ABar_LL yy',
--	'U ABar_LL yz',
--	'U ABar_LL zz',
-- [=[
--	'U LambdaBar_U mag',
--	'U LambdaBar_U x',
--	'U LambdaBar_U y',
--	'U LambdaBar_U z',
	'U H',
--	'U M_U mag',
--	'U M_U x',
--	'U M_U y',
--	'U M_U z',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_LL tr weighted',
--]=]

-- [[ debugging derivatives
--[=[
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
	'deriv M_U x',
	'deriv M_U y',
	'deriv M_U z',
--]=]
--]]

	--'U tr_DBar2_phi',
	--'U DBar_phi_sq',
--	'U ABarSq_LL tr weighted gammaBar^IJ',

--[[ should be zero for Minkowski
	'U RBar_LL xx',
	'U RBar_LL xy',
	'U RBar_LL xz',
	'U RBar_LL yy',
	'U RBar_LL yz',
	'U RBar_LL zz',
	'U RBar_LL tr weighted gammaBar^IJ',
--]]
--[[ should be zero for Minkowski
	'U gammaBar_UU xx',
	'U gammaBar_UU xy',
	'U gammaBar_UU xz',
	'U gammaBar_UU yy',
	'U gammaBar_UU yz',
	'U gammaBar_UU zz',
--]]
--[[ should be zero for Minkowski
	'U Delta_ULL x xx',
	'U Delta_ULL x xy',
	'U Delta_ULL x xz',
	'U Delta_ULL x yy',
	'U Delta_ULL x yz',
	'U Delta_ULL x zz',
	'U Delta_ULL y xx',
	'U Delta_ULL y xy',
	'U Delta_ULL y xz',
	'U Delta_ULL y yy',
	'U Delta_ULL y yz',
	'U Delta_ULL y zz',
	'U Delta_ULL z xx',
	'U Delta_ULL z xy',
	'U Delta_ULL z xz',
	'U Delta_ULL z yy',
	'U Delta_ULL z yz',
	'U Delta_ULL z zz',
--]]
--[[ should be zero for Minkowski
	'U Delta_U x',
	'U Delta_U y',
	'U Delta_U z',
--]]
--[[ should be diag(0, 2, 2 cos(theta)^2)
	'U trBar_partial2_gammaBar_ll xx',
	'U trBar_partial2_gammaBar_ll xy',
	'U trBar_partial2_gammaBar_ll xz',
	'U trBar_partial2_gammaBar_ll yy',
	'U trBar_partial2_gammaBar_ll yz',
	'U trBar_partial2_gammaBar_ll zz',
--]]
-- [=[ scalar wave variables
	'U Phi re',
	'U Psi_l x re',
	'U Psi_l y re',
	'U Psi_l z re',
	'U Pi re',
--]=]
--[=[
	'U tr_ABarSq',
	'U tr_DBar2_alpha',
--]=]
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:append{
		{name='S', code = self:template'value.vreal = real3s3_dot(U->S_ll, <?=calc_gamma_uu?>(U, x));'},

		{
			name='volume',
			code = self:template[[
	//|g| = exp(12 phi) |g_grid|
	real exp_neg4phi = <?=calc_exp_neg4phi?>(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
//// MODULE_DEPENDS: <?=calc_det_gammaHat?>
	real det_gamma = exp_12phi * <?=calc_det_gammaHat?>(x);
	value.vreal = U->alpha * det_gamma;
]],
		},
	}

	vars:append{
		{
			name = 'ABarSq_LL',
			type = 'real3s3',
			code = self:template[[
	real3s3 gammaBar_UU = <?=calc_gammaBar_UU?>(U, x);
	real3x3 ABar_UL = real3s3_real3s3_mul(gammaBar_UU, U->ABar_LL);
	real3s3 ABarSq_LL = real3s3_real3x3_to_real3s3_mul(U->ABar_LL, ABar_UL);
	value.vreal3s3 = ABarSq_LL;
]],
		},
	}

--[=[
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = self:template[[

<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>

	//partial_phi_l.i := phi_,i = -W_,i / (2 W)
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);

	real3s3 partial2_phi_ll;
	{
		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll.[<?=xij?>
				+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
			) / U->W;
<? end ?>
	}

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	real3s3 partial2_phi_LL = real3s3_rescaleFromCoord_ll(partial2_phi_ll, x);

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

	value.vreal = tr_DBar2_phi;
]]
		},


		{
			name = 'DBar2_phi_ll',
			code = self:template[[
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>

	//partial_phi_l.i := phi_,i = -W_,i / (2 W)
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);

	real3s3 partial2_phi_ll;
	{
		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll.<?=xij?>
				+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
			) / U->W;
<? end ?>
	}

<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	//DBar2_phi := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	real3s3 DBar2_phi_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_LL.<?=xij?> = partial2_phi_LL.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?		end
?>	;
<?	end
end
?>

	value.vreal = real3s3_dot(gammaBar_UU, DBar2_phi_LL);
]],
 		},
--]=]


	vars:append{
		{
			name = 'partial_phi_l',
			type = 'real3',
			code = self:template[[
<?=eqn:makePartial1'W'?>
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);
	value.vreal3 = partial_phi_l;
]],
		},
	}


	vars:append{
		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = self:template[[
<?=eqn:makePartial1'alpha'?>
	value.vreal3 = partial_alpha_l;
]],
		},
	}

--[=[
		{
			name = 'DBar2_alpha_ll',
			code = self:template[[
<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=eqn:makePartial1'epsilon_LL'?>

//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	real3s3 DBar2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_LL.<?=xij?> = partial2_alpha_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_alpha_L.<?=xk?>
<?	end
?>	;
<? end
?>

	value.vreal = real3s3_dot(gammaBar_UU, DBar2_alpha_LL);
]],
		},

		{
			name = 'tracelessPart_LL',
			type = 'real3s3',
			code = self:template[[

<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=eqn:makePartial1'epsilon_LL'?>

//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = calc_det_gammaBar(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);
	real3s3 gammaBar_ll = real3s3_rescaleToCoord_LL(gammaBar_LL, x);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	real3s3 partial2_alpha_LL = real3s3_rescaleFromCoord_ll(partial2_alpha_ll, x);

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	real3s3 DBar2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_LL.<?=xij?> = partial2_alpha_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_alpha_L.<?=xk?>
<?	end
?>	;
<? end
?>

<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>

	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);

	real3s3 partial2_phi_ll;
	{
		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				- partial2_W_ll.<?=xij?>
				+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
			) / U->W;
<? end ?>
	}

	real3s3 partial2_phi_LL = real3s3_rescaleFromCoord_ll(partial2_phi_ll, x);

	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	real3s3 DBar2_phi_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_LL.<?=xij?> = partial2_phi_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?	end
?>	;
<? end
?>

	real3 partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);

	real3s3 tracelessPart_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
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

	value.vreal3s3 = tracelessPart_LL;
]],
		},
--]=]

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
		code = self:template[[
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial1'beta_U'?>
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);

	//gamma_ij = exp(4 phi) gammaBar_ij
	real3s3 gamma_ll = real3s3_real_mul(<?=calc_gammaBar_ll?>(U, x), exp_4phi);

	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K
	real3s3 K_ll = real3s3_add(
		real3s3_real_mul(U->ABar_ll, exp_4phi),
		real3s3_real_mul(gamma_ll, U->K/3.));

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	value.vreal = -tr_partial_beta / U->alpha
<?
for i,xi in ipairs(xNames) do
?>		+ U->beta_u.<?=xi?> * partial_alpha_l.<?=xi?> / (U->alpha * U->alpha)
		- U->beta_u.<?=xi?> * partial_alpha_l.<?=xi?> / (U->alpha * U->alpha)
		+ 3. * partial_W_l.<?=xi?> * U->beta_u.<?=xi?> / (U->W * U->alpha)
<?	for j,xj in ipairs(xNames) do
?>		+ K_ll.<?=sym(i,j)?> * U->beta_u.<?=xi?> * U->beta_u.<?=xj?> / (U->alpha * U->alpha)
<?	end
end
?>		- U->K;
]],
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
		code = self:template[[
<?=eqn:makePartial1'alpha'?>

	real _1_alpha = 1. / U->alpha;

	real3s3 gamma_uu = <?=calc_gamma_uu?>(U, x);
	real3 partial_alpha_u = real3s3_real3_mul(gamma_uu, partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i

<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

<?=eqn:makePartial1'beta_U'?>
<?=eqn:makePartial1'epsilon_LL'?>

	//W = exp(-2 phi)
	real _1_W = 1. / U->W;

	//gamma_ij = W^-2 gammaBar_ij
	real3s3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
	real3s3 gamma_ll = real3s3_real_mul(gammaBar_ll, _1_W * _1_W);

	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
<?=eqn:makePartial1'W'?>
	real3x3s3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3s3_sub(
			real3s3_real_mul(partial_epsilon_lll[<?=i-1?>], _1_W * _1_W),
			real3s3_real_mul(gammaBar_ll, 2. * partial_W_l.<?=xi?> * _1_W * _1_W * _1_W)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	real3s3 dt_gamma_ll = real3s3_zero;

	real partial_alpha_dot_beta = real3_dot(U->beta_u, partial_alpha_l);	//beta^j alpha_,j

	real3 beta_l = real3s3_real3_mul(gamma_ll, U->beta_u);								//beta^j gamma_ij
	real3 beta_dt_gamma_l = real3s3_real3_mul(dt_gamma_ll, U->beta_u);					//beta^j gamma_ij,t
	real beta_beta_dt_gamma = real3_dot(U->beta_u, beta_dt_gamma_l);				//beta^i beta^j gamma_ij,t

	real3 beta_dt_gamma_u = real3s3_real3_mul(gamma_uu, beta_dt_gamma_l);				//gamma^ij gamma_jk,t beta^k

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
	real3 beta_dbeta_u = real3s3_real3_mul(gamma_uu, beta_dbeta_l);

	//gamma_kl,j beta^k beta^l
	real3 beta_beta_dgamma_l = real3{
<? for i,xi in ipairs(xNames) do
?>		real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>),
<? end
?>	};

	real3 beta_beta_dgamma_u = real3s3_real3_mul(gamma_uu, beta_beta_dgamma_l);

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
]],
		type = 'real3',
	}
--]=]
--[=[	-- TODO think about storing this instead of recalculating it here?
	do
		-- should be zero when gammaBar_ij = gammaHat_ij
		vars:insert{
			name = 'RBar_LL',
			type = 'real3s3',
			code = self:template[[
	//partial_LambdaBar_Ul[j].i := LambdaBar^i_,j
<?=eqn:makePartial1'LambdaBar_U'?>
	//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
//// MODULE_DEPENDS: real3x3_partial_rescaleFromCoord_Ul
	real3x3 partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);

<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

//// MODULE_DEPENDS: <?=mystery_C_U?>
	real3 Delta_U = real3_sub(U->LambdaBar_U, <?=mystery_C_U?>);

	real3x3s3 connHat_LLL, connHat_ULL;
//// MODULE_DEPENDS: calc_connHat_LLL_and_ULL
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	real3x3s3 Delta_ULL = real3x3s3_sub(connBar_ULL, connHat_ULL);
	real3x3s3 Delta_LLL = real3s3_real3x3s3_mul(gammaBar_LL, Delta_ULL);

<?=eqn:makePartial2'epsilon_LL'?>

//// MODULE_DEPENDS: calc_trBar_partial2_gammaBar_ll
	real3s3 trBar_partial2_gammaBar_ll = calc_trBar_partial2_gammaBar_ll(
		U,
		x,
		gammaBar_UU,
		partial_epsilon_LLl,
		partial2_epsilon_LLll);

	real3s3 gammaBar_ll = real3s3_rescaleToCoord_LL(gammaBar_LL, x);
	real3s3 gammaBar_uu = real3s3_rescaleToCoord_UU(gammaBar_UU, x);
//// MODULE_DEPENDS: calc_RBar_LL
	real3s3 RBar_LL = calc_RBar_LL(
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

	value.vreal3s3 = real3s3_rescaleToCoord_LL(RBar_LL, x);
]],
		}
	end

--]=]
--[=[
	vars:insert{
		name = 'RPhi_LL',
		type = 'real3s3',
		code = self:template[[

<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>

	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);

	real3s3 partial2_phi_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_phi_ll.<?=xij?> = .5 * (
			-partial2_W_ll.<?=xij?>
			+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
		) / U->W;
<? end ?>

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	real3s3 partial2_phi_LL = real3s3_rescaleFromCoord_ll(partial2_phi_ll, x);

	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	real3s3 DBar2_phi_LL = real3s3_sub(
		partial2_phi_LL,
		real3_real3x3s3_dot1(
			partial_phi_L,
			connBar_ULL
		)
	);

	real tr_DBar2_phi = real3s3_dot(gammaBar_UU, DBar2_phi_LL);

	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	real3s3 RPhi_LL = {
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

	value.vreal3s3 = RPhi_ll;
]],
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

--[=[
	vars:insert{
		name='trBar_partial2_gammaBar_ll',
		type = 'real3s3',
		code = self:template[[
	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);
<?=eqn:makePartial1'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

//// MODULE_DEPENDS: calc_trBar_partial2_gammaBar_ll
	real3s3 trBar_partial2_gammaBar_ll = calc_trBar_partial2_gammaBar_ll(
		U,
		x,
		gammaBar_UU,
		partial_epsilon_LLl,
		partial2_epsilon_LLll);

	value.vreal3s3 = trBar_partial2_gammaBar_ll;
]],
	}
--]=]

--[=[
	vars:insert{
		name = 'tr34 (gamma*dGamma)',
		type = 'real3x3',
		code = self:template[[

	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);
	real3s3 gammaBar_ll = real3s3_rescaleToCoord_LL(gammaBar_LL, x);

	real3s3 gammaBar_uu = real3s3_rescaleToCoord_UU(gammaBar_UU, x);

	real3x3 tr34_gamma_dGamma_ll;
<?
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr34_gamma_dGamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * gammaBar_LL.<?=sym(i,m)?> * calc_len_<?=xi?>(x) * calc_partial_connHat_Ulll_<?=xm..sym(j,k)..xl?>(x)
<?				end
			end
		end
?>	;
<?	end
end
?>

	value.vreal3x3 = tr34_gamma_dGamma_ll;
]],
	}
--]=]

--[=[
	vars:insert{
		name = 'tr14 (Gamma*dgamma)',
		type = 'real3x3',
		code = self:template[[
<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);

	real3x3s3 connHat_LLL, connHat_ULL;
//// MODULE_DEPENDS: calc_connHat_LLL_and_ULL
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	real3x3s3 connHat_ull = real3x3s3_rescaleToCoord_ULL(connHat_ULL, x);

	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	real3x3s3 partial_gammaBar_lll = real3x3s3_rescaleToCoord_LLL(partial_gammaBar_LLL, x);

	real3s3 gammaBar_uu = real3s3_rescaleToCoord_UU(gammaBar_UU, x);

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
	value.vreal3x3 = tr14_Gamma_dgamma_ll;
]],
	}
--]=]

--[=[
	for i,xi in ipairs(xNames) do
		vars:insert{
			name = 'Delta_ULL '..xi,
			type = 'real3s3',
			code = self:template[[
	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	real3x3s3 connHat_LLL, connHat_ULL;
//// MODULE_DEPENDS: calc_connHat_LLL_and_ULL
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK + connHat^I_JK
	real3x3s3 Delta_ULL = real3x3s3_sub(connBar_ULL, connHat_ULL);

	value.vreal3s3 = Delta_ULL.<?=xi?>;
]], {i=i, xi=xi}),
		}
	end
--]=]

--[=[
	vars:insert{
		name = 'Delta_U',
		type = 'real3',
		code = self:template[[
	real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
	real det_gammaBarLL = <?=calc_det_gammaBarLL?>(x);
	real3s3 gammaBar_UU = gammaBar_LL.inverse(det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);

	real3x3s3 connHat_LLL, connHat_ULL;
//// MODULE_DEPENDS: calc_connHat_LLL_and_ULL
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK + connHat^I_JK
	real3x3s3 Delta_ULL = real3x3s3_sub(connBar_ULL, connHat_ULL);

<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = real3s3_dot(gammaBar_UU, Delta_ULL.<?=xi?>);
<? end
?>
]],
	}
--]=]

-- [=[
	vars:insert{
		name='tr_ABarSq',
		code = self:template[[
	real3s3 gammaBar_UU = <?=calc_gammaBar_UU?>(U, x);
	real3x3 ABar_UL = real3s3_real3s3_mul(gammaBar_UU, U->ABar_LL);
	real3s3 ABar_UU = real3x3_real3s3_to_real3s3_mul(ABar_UL, gammaBar_UU);
	real tr_ABarSq = real3s3_dot(U->ABar_LL, ABar_UU);
	value.vreal = tr_ABarSq;
]],
	}
--]=]

-- [=[
	vars:insert{
		name = 'tr_DBar2_alpha',
		code = self:template[[
<?=eqn:makePartial1'alpha'?>
	real3 partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);
<?=eqn:makePartial2'alpha'?>		//partial2_alpha_ll.ij := alpha_,ij
	real3s3 partial2_alpha_LL = real3s3_rescaleFromCoord_ll(partial2_alpha_ll, x);
	real3s3 gammaBar_UU = <?=calc_gammaBar_UU?>(U, x);
<?=eqn:makePartial1'epsilon_LL'?>
//// MODULE_DEPENDS: calc_partial_gammaBar_LLL
	real3x3s3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
//// MODULE_DEPENDS: calc_connBar_ULL
	real3x3s3 connBar_ULL = calc_connBar_ULL(partial_gammaBar_LLL, gammaBar_UU);
	real3s3 DBar2_alpha_LL = real3s3_sub(partial2_alpha_LL, real3_real3x3s3_dot1(partial_alpha_L, connBar_ULL));
	real tr_DBar2_alpha = real3s3_dot(gammaBar_UU, DBar2_alpha_LL);
	value.vreal = tr_DBar2_alpha;
]],
	}
--]=]

	vars:insert{name='x', type='real3', code='value.vreal3=x;'}

	return vars
end

return BSSNOKFiniteDifferenceEquation
