--[[
Based on 2008 Yano et al "Flux Vector Splitting..."
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local Z4_2008Yano = class(EinsteinEqn)
Z4_2008Yano.name = 'Z4 (2008 Yano et al)'

local fluxVars = table{
	{name='a_l', type='real3'},		-- 3:  0-2
	{name='d_lll', type='_3sym3'},	-- 18: 3-20
	{name='K_ll', type='sym3'},		-- 6:  21-26
	{name='Theta', type='real'},	-- 1:  27
	{name='Z_l', type='real3'},		-- 3:  28-30
}

Z4_2008Yano.consVars = table{
	{name='alpha', type='real'},
	{name='gamma_ll', type='sym3'},
}:append(fluxVars)

Z4_2008Yano.numWaves = Struct.countScalars{vars=fluxVars}
assert(Z4_2008Yano.numWaves == 31)
	
Z4_2008Yano.numIntStates = Struct.countScalars{vars=Z4_2008Yano.consVars}

Z4_2008Yano.consVars:append{
	--constraints:              
	{name='H', type='real'},				--1
	{name='M_u', type='real3'},				--3
}

Z4_2008Yano.useSourceTerm = true

function Z4_2008Yano:createInitState()
	Z4_2008Yano.super.createInitState(self)
	self:addGuiVar{name = 'm', value = 2}
end

function Z4_2008Yano:initCodeModules()
	Z4_2008Yano.super.initCodeModules(self)

	for moduleName, depends in pairs{
		['calcDT'] = {},
		['eigen_forCell'] = {},
		['calcCellMinMaxEigenvalues'] = {},
		['eigen_forInterface'] = {},
		['eigen_left/rightTransform'] = {},
		['eigen_fluxTransform'] = {},
		['addSource'] = {},
	} do
		self:addModuleFromSourceFile{
			name = moduleName,
			depends = depends,
		}
	end
end

-- don't use default
function Z4_2008Yano:initCodeModule_calcDT() end

Z4_2008Yano.needsInitDerivs = true

Z4_2008Yano.solverCodeFile = 'hydro/eqn/z4_2008yano.cl'

Z4_2008Yano.predefinedDisplayVars = {
	'U alpha',
	'U gamma_ll x x',
	'U d_lll_x x x',
	'U K_ll x x',
	'U Theta',
	'U Z_l x',
	'U H',
	'U M_u',
	'U volume',
	'U f',
}

function Z4_2008Yano:getDisplayVars()
	local vars = Z4_2008Yano.super.getDisplayVars(self)
	vars:append{
		{name='volume', code='value.vreal = U->alpha * sqrt(sym3_det(U->gamma_ll));'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='expansion', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal = -sym3_dot(gamma_uu, U->K_ll);
]]		},
--[=[
	-- 1998 Bona et al
--[[
H = 1/2 ( R + K^2 - K_ij K^ij ) - alpha^2 8 pi rho
for 8 pi rho = G^00

momentum constraints
--]]
		{H = [[
	.5 * 
]]		},
--]=]

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
		{name='gravity', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal3 = real3_real_mul(sym3_real3_mul(gamma_uu, U->a_l), -U->alpha * U->alpha);
]], type='real3'},
	}
	
	return vars
end

Z4_2008Yano.eigenVars = table{
	{name='alpha', type='real'},
	{name='sqrt_f', type='real'},
	{name='gamma_ll', type='sym3'},
	{name='gamma_uu', type='sym3'},
	{name='sqrt_gammaUjj', type='real3'},
}

function Z4_2008Yano:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	<? if side==0 then ?>
	real eig_lambdaLight = <?=eig?>.alpha * <?=eig?>.sqrt_gammaUjj.x;
	<? elseif side==1 then ?>
	real eig_lambdaLight = <?=eig?>.alpha * <?=eig?>.sqrt_gammaUjj.y;
	<? elseif side==2 then ?>
	real eig_lambdaLight = <?=eig?>.alpha * <?=eig?>.sqrt_gammaUjj.z;
	<? end ?>
	real eig_lambdaGauge = eig_lambdaLight * <?=eig?>.sqrt_f;
]], {
		eig = '('..eig..')',
		side = side,
	})
end

function Z4_2008Yano:eigenWaveCode(side, eig, x, waveIndex)
	local betaUi
	if self.useShift then
		betaUi = '('..eig..').beta_u.'..xNames[side+1]
	else
		betaUi = '0'
	end

	if waveIndex == 0 then
		return '-'..betaUi..' - eig_lambdaGauge'
	elseif waveIndex >= 1 and waveIndex <= 6 then
		return '-'..betaUi..' - eig_lambdaLight'
	elseif waveIndex >= 7 and waveIndex <= 23 then
		return '-'..betaUi
	elseif waveIndex >= 24 and waveIndex <= 29 then
		return '-'..betaUi..' + eig_lambdaLight'
	elseif waveIndex == 30 then
		return '-'..betaUi..' + eig_lambdaGauge'
	end
	error'got a bad waveIndex'
end

function Z4_2008Yano:consWaveCodePrefix(side, U, x, waveIndex)
	return template([[
	real det_gamma = sym3_det(<?=U?>.gamma_ll);
	sym3 gamma_uu = sym3_inv(<?=U?>.gamma_ll, det_gamma);
	<? if side==0 then ?>
	real eig_lambdaLight = <?=U?>.alpha * sqrt(gamma_uu.xx);
	<? elseif side==1 then ?>                          
	real eig_lambdaLight = <?=U?>.alpha * sqrt(gamma_uu.yy);
	<? elseif side==2 then ?>                          
	real eig_lambdaLight = <?=U?>.alpha * sqrt(gamma_uu.zz);
	<? end ?>
	real f = calc_f(<?=U?>.alpha);
	real eig_lambdaGauge = eig_lambdaLight * sqrt(f);
]], {
		U = '('..U..')',
		side = side,
	})
end
Z4_2008Yano.consWaveCode = Z4_2008Yano.eigenWaveCode

return Z4_2008Yano
