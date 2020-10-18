--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" 2008 chapter on Toy 1+1 spacetimes.

See comments in my gravitation-waves project adm1d_v1.lua file for the math.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinEqn = require 'hydro.eqn.einstein'

local ADM_BonaMasso_1D_2008Alcubierre = class(EinsteinEqn)

ADM_BonaMasso_1D_2008Alcubierre.name = 'ADM_BonaMasso_1D_2008Alcubierre' 

ADM_BonaMasso_1D_2008Alcubierre.consVars = {
	{name='alpha', type='real'}, 
	{name='gamma_xx', type='real', variance=''}, 
	{name='a_x', type='real', variance=''}, 
	{name='D_g', type='real', variance=''}, 
	{name='KTilde', type='real'},
}
ADM_BonaMasso_1D_2008Alcubierre.numWaves = 3	-- alpha and gamma_xx are source-term only

ADM_BonaMasso_1D_2008Alcubierre.useSourceTerm = true
ADM_BonaMasso_1D_2008Alcubierre.roeUseFluxFromCons = true

ADM_BonaMasso_1D_2008Alcubierre.guiVars = {
	{name='a_x_convCoeff', value=10},
	{name='D_g_convCoeff', value=10},
}

function ADM_BonaMasso_1D_2008Alcubierre:initCodeModules()
	ADM_BonaMasso_1D_2008Alcubierre.super.initCodeModules(self)
	for moduleName, depends in pairs{

		['setFlatSpace'] = {},

--[=[
		-- the PLM version that uses this crashes
		-- so maybe there's something wrong with this
		['fluxFromCons'] = {
			'solver.solver_t',
			'eqn.cons_t',
			'eqn.common',	-- calc_f ... or is it initCond.codeprefix?
		},
--]=]

		['eigen_forInterface'] = {},
		['eigen_forCell'] = {},
		['eigen_left/rightTransform'] = {},
		['eigen_fluxTransform'] = {},
		
		['addSource'] = {
			'initCond.codeprefix',	-- calc_*
		},
	} do
		self:addModuleFromSourceFile{
			name = moduleName,
			depends = depends,
		}
	end
end

--[=[ enable this out if you enable the fluxFromCons above
function ADM_BonaMasso_1D_2008Alcubierre:initCodeModule_fluxFromCons() end
--]=]

-- don't use default
function ADM_BonaMasso_1D_2008Alcubierre:initCodeModule_setFlatSpace() end

function ADM_BonaMasso_1D_2008Alcubierre:getModuleDependsApplyInitCond()
	return table(ADM_BonaMasso_1D_2008Alcubierre.super.getModuleDependsApplyInitCond(self), {
		'sym3',
		'coordMap',
	})
end

-- don't use eqn.einstein, which says calc_gamma_ll and calc_gamma_uu
function ADM_BonaMasso_1D_2008Alcubierre:getModuleDependsSolver() 
	return {}
end

-- don't use eqn.einstein:
function ADM_BonaMasso_1D_2008Alcubierre:createDisplayComponents() end

ADM_BonaMasso_1D_2008Alcubierre.needsInitDerivs = true

ADM_BonaMasso_1D_2008Alcubierre.solverCodeFile = 'hydro/eqn/adm1d_v1.cl'

function ADM_BonaMasso_1D_2008Alcubierre:getDisplayVars()
	return ADM_BonaMasso_1D_2008Alcubierre.super.getDisplayVars(self):append{
		-- adm1d_v2 cons vars:
		{name='d_xxx', code='value.vreal = .5 * U->D_g * U->gamma_xx;'},
		{name='K_xx', code='value.vreal = U->KTilde * sqrt(U->gamma_xx);'},
		-- aux:
		{name='dx_alpha', code='value.vreal = U->alpha * U->a_x;'},
		{name='dx_gamma_xx', code='value.vreal = U->gamma_xx * U->D_g;'},
		{name='volume', code='value.vreal = U->alpha * sqrt(U->gamma_xx);'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		{name='K', code='value.vreal = U->KTilde / sqrt(U->gamma_xx);'},
		{name='expansion', code='value.vreal = -U->KTilde / sqrt(U->gamma_xx);'},
		{name='gravity mag', code='value.vreal = -U->alpha * U->alpha * U->a_x / U->gamma_xx;'},
	
		{name='alpha vs a_x', code=[[
	if (OOB(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_alpha - U->alpha * U->a_x);
	}
]]},

		{name='gamma_xx vs D_g', code=[[
	if (OOB(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_gamma_xx - U->gamma_xx * U->D_g);
	}
]]},
	}
end

ADM_BonaMasso_1D_2008Alcubierre.eigenVars = {
	{name='f', type='real'},
	{name='alpha', type='real'},
	{name='gamma_xx', type='real'},
}

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCodePrefix(n, eig, x, waveIndex)
	return self:template([[
	real eig_lambda = <?=eig?>.alpha * sqrt(<?=eig?>.f / <?=eig?>.gamma_xx);
]], {
		eig = '('..eig..')',
	})
end

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCode(n, eig, x, waveIndex)
	if waveIndex == 0 then
		return '-eig_lambda'
	elseif waveIndex == 1 then
		return '0'
	elseif waveIndex == 2 then
		return 'eig_lambda'
	else
		error'got a bad waveIndex'
	end
end

function ADM_BonaMasso_1D_2008Alcubierre:consWaveCodePrefix(n, U, x, waveIndex)
	return self:template([[
	real f = calc_f(<?=U?>.alpha);
	real eig_lambda = <?=U?>.alpha * sqrt(f / <?=U?>.gamma_xx);
]], {
		U = '('..U..')',
	})
end

ADM_BonaMasso_1D_2008Alcubierre.consWaveCode = ADM_BonaMasso_1D_2008Alcubierre.eigenWaveCode
	
return ADM_BonaMasso_1D_2008Alcubierre
