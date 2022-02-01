--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'
local xNames = common.xNames


local ADM_BonaMasso_3D = class(EinsteinEqn)
ADM_BonaMasso_3D.name = 'ADM_BonaMasso_3D'

--[[
args:
noZeroRowsInFlux = true by default.
	true = use 13 vars of a_x, d_xij, K_ij
	false = use all 30 or so hyperbolic conservation variables
useShift
	
	useShift = 'none'
	
	useShift = 'MinimalDistortionElliptic' -- minimal distortion elliptic via Poisson relaxation.  Alcubierre's book, eqn 4.3.14 and 4.3.15
	
	useShift = 'MinimalDistortionEllipticEvolve' -- minimal distortion elliptic via evolution.  eqn 10 of 1996 Balakrishna et al "Coordinate Conditions and their Implementations in 3D Numerical Relativity" 

	useShift = '2005 Bona / 2008 Yano'
	-- 2008 Yano et al, from 2005 Bona et al "Geometrically Motivated..."
	-- 2005 Bona mentions a few, but 2008 Yano picks the first one from the 2005 Bona paper.

	useShift = 'HarmonicShiftCondition-FiniteDifference'
	-- 2008 Alcubierre 4.3.37
	-- I see some problems in the warp bubble test ...

	useShift = 'LagrangianCoordinates'
	--[=[
	Step backwards along shift vector and advect the state
	Idk how accurate this is ...
	Hmm, even if I implement the Lie derivative as Lagrangian coordinate advection
	I'll still be responsible for setting some beta^i_,t gauge
	so for the L_beta Lie derivative, we have some options:
	1) none (for no-shift)
	2) finite difference
	3) finite volume / absorb into the eigensystem
	4) Lagrangian coordinates
	and this should be a separate variable, separate of the shift gauge

	so 
	one variable for what beta^i_,t is
	another variable for how to 
	--]=]
--]]
function ADM_BonaMasso_3D:init(args)
	local solver = assert(args.solver)

	local fluxVars = table{
		{name='a_l', type='real3'},
		{name='d_lll', type='_3sym3'},
		{name='K_ll', type='sym3'},
		{name='V_l', type='real3'},
	}

	self.consVars = table{
		{name='alpha', type='real'},
		{name='gamma_ll', type='sym3'},
	}:append(fluxVars)


	--[[
	how are shift conditions impmlemented?
	options for determining beta^i:
	1) by solving a constraint equation (minimal distortion elliptic solves it with a numerical Poisson solver)
	2) by solving an initial value problem

	Once beta^i is determined, how are the variables iterated?
	This question is split into a) and b):
	a) How are the source-only variables iterated wrt beta^i?
	options:
	1) put the Lie derivative terms into the source side of these variables 
	2) give them -beta^i eigenvalues?  this is the equivalent of rewriting the hyperbolic vars associated with these (a_i,d_kij) back into first-derivative 0th order vars (alpha, gamma_ij)

	b) How are the flux variables iterated wrt beta^i?
	options:
	1) this would be solved by offsetting the eigenvalues 
		with the offset eigenvalues, 
		the hyperbolic state vars' contributions get incorporated into the flux,
		but there are still some shift-based terms that end up in the source ...
		... so the shift is split between the flux and source ...
	2) ... and adding a few terms to the source
	--]]

	self.useShift = args.useShift or 'none'
	
	-- set to false to disable rho, S_i, S^ij
	self.useStressEnergyTerms = true

	if self.useShift ~= 'none' then
		self.consVars:insert{name='beta_u', type='real3'}

		if self.useShift == 'MinimalDistortionElliptic' 
		or self.useShift == 'MinimalDistortionEllipticEvolve' 
		then
			self.consVars:insert{name='betaLap_u', type='real3'}
		end
	end


	--[[
	solve a smaller eigendecomposition that doesn't include the rows of variables whose d/dt is zero.
	kind of like how ideal MHD is an 8-var system, but the flux jacobian solved is a 7x7 because Bx,t = 0.
	TODO make this a ctor arg - so solvers can run in parallel with and without this
	...or why don't I just scrap the old code, because this runs a lot faster.
	--]]
	self.noZeroRowsInFlux = true
	--self.noZeroRowsInFlux = false
	if args.noZeroRowsInFlux ~= nil then
		self.noZeroRowsInFlux = args.noZeroRowsInFlux
	end

	-- NOTE this doesn't work when using shift ... because then all the eigenvalues are -beta^i, so none of them are zero (except the source-only alpha, beta^i, gamma_ij)
	-- with the exception of the lagrangian shift.  if we split that operator out then we can first solve the non-shifted system, then second advect it by the shift vector ...
	--if self.useShift ~= 'none' then
	--	self.noZeroRowsInFlux = false
	--end
	
	-- only count int vars after the shifts have been added
	self:cdefAllVarTypes(solver, self.consVars)	-- have to call before countScalars in eqn:init

	if not self.noZeroRowsInFlux then
		-- skip alpha and gamma
		self.numWaves = Struct.countScalars{vars=fluxVars}
		assert(self.numWaves == 30)
	else
		-- skip alpha, gamma_ij, a_q, d_qij, V_i for q != the direction of flux
		self.numWaves = 13
	end

	self.numIntStates = Struct.countScalars{vars=self.consVars}

	-- now add in the source terms (if you want them)
	if self.useStressEnergyTerms then
		self.consVars:append{
			--stress-energy variables:
			{name='rho', type='real'},					--1: n_a n_b T^ab
			{name='S_u', type='real3'},				--3: -gamma^ij n_a T_aj
			{name='S_ll', type='sym3'},				--6: gamma_i^c gamma_j^d T_cd
		}								
	end
	self.consVars:append{
		--constraints:              
		{name='H', type='real'},					--1
		{name='M_u', type='real3'},				--3
		-- TODO stress constraint as well?  
	}

	self.eigenVars = table{
		{name='alpha', type='real'},
		{name='alpha_sqrt_f', type='real'},
		{name='gamma_uu', type='sym3'},
		-- sqrt(gamma^jj) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
		{name='sqrt_gammaUjj', type='real3'},
	}

	-- hmm, only certain shift methods actually use beta_u ...
	if self.useShift ~= 'none' then
		self.eigenVars:insert{name='beta_u', type='real3'}
	end

	
	-- build stuff around consVars	
	ADM_BonaMasso_3D.super.init(self, args)


	if self.useShift == 'MinimalDistortionElliptic' then
		local MinimalDistortionEllipticShift = require 'hydro.op.gr-shift-mde'
		self.solver.ops:insert(MinimalDistortionEllipticShift{solver=self.solver})
	elseif self.useShift == 'LagrangianCoordinates' then
		local LagrangianCoordinateShift = require 'hydro.op.gr-shift-lc'
		self.solver.ops:insert(LagrangianCoordinateShift{solver=self.solver})
	end
end

function ADM_BonaMasso_3D:createInitState()
	ADM_BonaMasso_3D.super.createInitState(self)
	
	self:addGuiVars{
		{
			type = 'combo',
			name = 'constrain V',
			options = {
				'none',	-- as long as there is a damping term to the source, direct constraint methods aren't required.
				'replace V',
				'average',	-- TODO add averaging weights, from 100% V (which works) to 100% d (which doesn't yet work)
			},
			compileTime = true,
		},

		-- K_ij source term stress constraint coefficient
		-- 1 of these is occurring in the ADM formalism
		{name='K_ll_srcStressCoeff', value=0},
		
		-- K_ij source term Hamiltonian constraint coefficient
		{name='K_ll_srcHCoeff', value=0},-- -.5},

		-- gamma_ij source term stress constraint coefficient
		{name='gamma_ll_srcStressCoeff', value=0},	--1},
	
		-- gamma_ij source term Hamiltonian constraint coefficient
		{name='gamma_ll_srcHCoeff', value=0}, ---.5},

		-- convergence between finite-difference of alpha,i and alpha a_i
		{name='a_convCoeff', value=0},
		
		-- convergence between finite-difference of 1/2 gamma_ij,k and d_kij
		{name='d_convCoeff', value=0},
		
		-- convergence between V_k and d_kj^j - d^j_jk
		{name='V_convCoeff', value=10},
	}
end

-- don't use default
function ADM_BonaMasso_3D:initCodeModule_fluxFromCons() end

function ADM_BonaMasso_3D:getModuleDepends_waveCode()
	return {
		self.symbols.initCond_codeprefix,	-- calc_f
	}
end

function ADM_BonaMasso_3D:getModuleDepends_displayCode()
	return {
		self.symbols.calc_gamma_ll,
		self.symbols.calc_gamma_uu,
		self.symbols.initCond_codeprefix,	-- calc_f
	}
end

ADM_BonaMasso_3D.solverCodeFile = 'hydro/eqn/adm3d.cl'

ADM_BonaMasso_3D.predefinedDisplayVars = {
	'U alpha',
	'U gamma_ll x x',
	'U a_l x',
	'U d_lll x x x',
	'U K_ll x x',
	'U V_l x',
	'U K_ll tr weighted gamma^ij',	-- same as K = K^i_i = K_ij gamma^ij
	'U H',
	'U volume',
	'U f',
}

function ADM_BonaMasso_3D:getDisplayVars()
	local vars = ADM_BonaMasso_3D.super.getDisplayVars(self)

	vars:append{
		{name='volume', code='value.vreal = U->alpha * sqrt(sym3_det(U->gamma_ll));'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		{name='expansion', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal = -sym3_dot(gamma_uu, U->K_ll);
]]		},
	}:append{
--[=[
--[[
Alcubierre 3.1.1
Baumgarte & Shapiro 2.90
H = R + K^2 - K^ij K_ij - 16 pi rho
for rho = n_a n_b T^ab (B&S eqn 2.89)
and n_a = -alpha t_,a (B&S eqns 2.19, 2.22, 2.24)

momentum constraints
--]]
		{H = [[
	.5 * 
]]		},
--]=]
	}

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
	vars:insert{name='gravity', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal3 = real3_real_mul(sym3_real3_mul(gamma_uu, U->a_l), -U->alpha * U->alpha);
]], type='real3'}

	vars:insert{name='alpha vs a_i', code=self:template[[
	if (<?=OOB?>(1,1)) {
		value.vreal3 = real3_zero;
	} else {
		<? for i=1,solver.dim do
			local xi = xNames[i]
		?>{
			real partial_i_log_alpha = (
				log(U[solver->stepsize.<?=xi?>].alpha) 
				- log(U[-solver->stepsize.<?=xi?>].alpha)
			) / (2. * solver->grid_dx.s<?=i-1?>);
			value.vreal3.<?=xi?> = fabs(partial_i_log_alpha - U->a_l.<?=xi?>);
		}<? end ?>
		<? for i=solver.dim+1,3 do
			local xi = xNames[i]
		?>{
			value.vreal3.<?=xi?> = 0;
		}<? end ?>
	}
]], type='real3'}

	-- d_kij = gamma_ij,k
	for i,xi in ipairs(xNames) do
		vars:insert{name='gamma_ij vs d_'..xi..'ij', code=self:template([[
	if (<?=OOB?>(1,1)) {
		value.vsym3 = sym3_zero;
	} else {
		<? if i <= solver.dim then ?>
		sym3 partial_i_gamma_ll = sym3_real_mul(
			sym3_sub(
				U[solver->stepsize.<?=xi?>].gamma_ll, 
				U[-solver->stepsize.<?=xi?>].gamma_ll
			), 
			1. / (2. * solver->grid_dx.s<?=i-1?>)
		);
		<? else ?>
		sym3 partial_i_gamma_ll = sym3_zero;
		<? end ?>
		value.vsym3 = sym3_sub(sym3_real_mul(partial_i_gamma_ll, .5), U->d_lll.<?=xi?>);
		value.vsym3 = (sym3){<?
	for jk,xjk in ipairs(symNames) do 
?>			.<?=xjk?> = fabs(value.vsym3.<?=xjk?>),
<?	end
?>		};
	}
]], {
	i = i,
	xi = xi,
}), type='sym3'}
	end

	vars:insert{name='V constraint', code=self:template[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d_lll.<?=xi?>, gamma_uu);
		real d2 = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + U->d_lll.<?=xj?>.<?=sym(k,i)?> * gamma_uu.<?=sym(j,k)?><?
		end
	end ?>;
		value.vreal3.<?=xi?> = U->V_l.<?=xi?> - (d1 - d2);
	}<? end ?>
]], type='real3'}

	return vars
end

function ADM_BonaMasso_3D:eigenWaveCodePrefix(args)
	return self:template([[
real sqrt_gammaUjj = 0./0.;
if (<?=n?>.side == 0) {
	sqrt_gammaUjj = (<?=eig?>)->sqrt_gammaUjj.x;
} else if (<?=n?>.side == 1) {
	sqrt_gammaUjj = (<?=eig?>)->sqrt_gammaUjj.y;
} else if (<?=n?>.side == 2) {
	sqrt_gammaUjj = (<?=eig?>)->sqrt_gammaUjj.z;
}
real const eig_lambdaLight = sqrt_gammaUjj * (<?=eig?>)->alpha;
real const eig_lambdaGauge = sqrt_gammaUjj * (<?=eig?>)->alpha_sqrt_f;
]], args)
end

function ADM_BonaMasso_3D:eigenWaveCode(args)
	-- TODO find out if -- if we use the lagrangian coordinate shift operation -- do we still need to offset the eigenvalues by -beta^i?
	local shiftingLambdas = self.useShift ~= 'none'
		--and self.useShift ~= 'LagrangianCoordinates'

	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..args.eig..').beta_u.s['..args.n..'.side]'
	else
		betaUi = '0'
	end

	if not self.noZeroRowsInFlux then
		if args.waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif args.waveIndex >= 1 and args.waveIndex <= 5 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif args.waveIndex >= 6 and args.waveIndex <= 23 then
			return '-'..betaUi
		elseif args.waveIndex >= 24 and args.waveIndex <= 28 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif args.waveIndex == 29 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	else	-- noZeroRowsInFlux 
		-- noZeroRowsInFlux implies useShift == 'none'
		if args.waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif args.waveIndex >= 1 and args.waveIndex <= 5 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif args.waveIndex == 6 then
			return '-'..betaUi
		elseif args.waveIndex >= 7 and args.waveIndex <= 11 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif args.waveIndex == 12 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	end
	error'got a bad waveIndex'
end

function ADM_BonaMasso_3D:consWaveCodePrefix(args)
	return self:template([[
real const det_gamma = sym3_det((<?=U?>)->gamma_ll);
sym3 const gamma_uu = sym3_inv((<?=U?>)->gamma_ll, det_gamma);
real sqrt_gammaUjj = 0./0.;
if (<?=n?>.side == 0) {
	sqrt_gammaUjj = sqrt(gamma_uu.xx);
} else if (<?=n?>.side == 1) {
	sqrt_gammaUjj = sqrt(gamma_uu.yy);
} else if (<?=n?>.side == 2) {
	sqrt_gammaUjj = sqrt(gamma_uu.zz);
}
real const eig_lambdaLight = (<?=U?>)->alpha * sqrt_gammaUjj;
real const f_alphaSq = calc_f_alphaSq((<?=U?>)->alpha);
real const eig_lambdaGauge = sqrt_gammaUjj * sqrt(f_alphaSq);
]], args)
end

function ADM_BonaMasso_3D:consWaveCode(args)
	args = table(args):setmetatable(nil)
	args.U = args.eig
	args.eig = nil
	return self:eigenWaveCode(args)
end

function ADM_BonaMasso_3D:consWaveCodeMinMaxAllSidesPrefix(args)
	return self:template([[
real const f_alphaSq = calc_f_alphaSq((<?=U?>)->alpha);
real const det_gamma = sym3_det((<?=U?>)->gamma_ll);
real const alpha_sqrt_f = sqrt(f_alphaSq);
]], args)
end

--//// MODULE_DEPENDS: <?=SETBOUNDS?> <?=initCond_codeprefix?> <?=eqn_guiVars_compileTime?>
function ADM_BonaMasso_3D:consWaveCodeMinMaxAllSides(args)
	return self:template([[
real gammaUjj = 0./0.;
if ((<?=n?>).side == 0) {
	gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
} else if ((<?=n?>).side == 1) {
	gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
} else if ((<?=n?>).side == 2) {
	gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
}
real const sqrt_gammaUjj = sqrt(gammaUjj);
real const lambdaLight = sqrt_gammaUjj * U->alpha;
real const lambdaGauge = sqrt_gammaUjj * alpha_sqrt_f;

real const lambda = (real)max(lambdaGauge, lambdaLight);

<? if eqn.useShift ~= "none" then ?>
real const betaUi = U->beta_u.s<?=side?>;
<? else ?>
real const betaUi = 0.;
<? end ?>

<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'(real)min((real)0., -betaUi - lambda);',
	'(real)max((real)0., -betaUi + lambda);')?>
]], args)
end

return ADM_BonaMasso_3D
