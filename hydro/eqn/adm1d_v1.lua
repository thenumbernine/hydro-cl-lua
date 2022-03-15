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

ADM_BonaMasso_1D_2008Alcubierre.roeUseFluxFromCons = true

ADM_BonaMasso_1D_2008Alcubierre.guiVars = {
	{name='a_x_convCoeff', value=10},
	{name='D_g_convCoeff', value=10},
}

-- [=[ enable this out if you enable the fluxFromCons above
-- the PLM version that uses this crashes
-- so maybe there's something wrong with this
function ADM_BonaMasso_1D_2008Alcubierre:initCodeModule_fluxFromCons() 
	self.solver.modules:addFromMarkup{
		code = self:template[[
//// MODULE_NAME: <?=fluxFromCons?>
//// MODULE_DEPENDS: <?=solver_t?> <?=cons_t?>
#define <?=fluxFromCons?>(\
	/*<?=cons_t?> * const */resultFlux,\
	/*constant <?=solver_t?> const * const */solver,\
	/*<?=cons_t?> const * const */U,\
	/*<?=cell_t?> const * const */cell,\
	/*<?=normal_t?> const */n\
) {\
//// MODULE_DEPENDS: <?=initCond_codeprefix?>
	real const f = calc_f((U)->alpha);\
	real const alpha_over_sqrt_gamma_xx = (U)->alpha / sqrt((U)->gamma_xx);\
	(resultFlux)->alpha = 0;\
	(resultFlux)->gamma_xx = 0;\
	(resultFlux)->a_x = (U)->KTilde * f * alpha_over_sqrt_gamma_xx;\
	(resultFlux)->D_g = (U)->KTilde * 2. * alpha_over_sqrt_gamma_xx;\
	(resultFlux)->KTilde = (U)->a_x * alpha_over_sqrt_gamma_xx;\
}

]],
	}
end
--]=]

-- don't use eqn.einstein:
function ADM_BonaMasso_1D_2008Alcubierre:createDisplayComponents() end

ADM_BonaMasso_1D_2008Alcubierre.solverCodeFile = 'hydro/eqn/adm1d_v1.cl'

ADM_BonaMasso_1D_2008Alcubierre.predefinedDisplayVars = {
	'U alpha',
	'U gamma_xx',
	'U a_x',
	'U d_xxx',
	'U K_xx',
	'U K',
	-- TODO H constraint
	'U volume',
	'U alpha vs a_x',
	'U gamma_xx vs D_g',
	'U f',
}

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
	
		{name='alpha vs a_x', code=self:template[[
	if (<?=OOB?>(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_alpha - U->alpha * U->a_x);
	}
]]},

		{name='gamma_xx vs D_g', code=self:template[[
	if (<?=OOB?>(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_gamma_xx - U->gamma_xx * U->D_g);
	}
]]},
	}
end

ADM_BonaMasso_1D_2008Alcubierre.eigenVars = {
	{name='alpha', type='real'},
	{name='sqrt_gamma_xx', type='real'},
	{name='sqrt_f', type='real'},
}

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCodePrefix(args)
	return self:template([[
real const eig_lambda = (<?=eig?>)->alpha * (<?=eig?>)->sqrt_f / (<?=eig?>)->sqrt_gamma_xx;
]], args)
end

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCode(args)
	if args.waveIndex == 0 then
		return '-eig_lambda'
	elseif args.waveIndex == 1 then
		return '0'
	elseif args.waveIndex == 2 then
		return 'eig_lambda'
	else
		error'got a bad waveIndex'
	end
end

function ADM_BonaMasso_1D_2008Alcubierre:consWaveCodePrefix(args)
	return self:template([[
real const alphaSq_f = calc_f_alphaSq((<?=U?>)->alpha);
real const eig_lambda = sqrt(alphaSq_f / (<?=U?>)->gamma_xx);
]], args)
end

ADM_BonaMasso_1D_2008Alcubierre.consWaveCode = ADM_BonaMasso_1D_2008Alcubierre.eigenWaveCode
	
return ADM_BonaMasso_1D_2008Alcubierre
