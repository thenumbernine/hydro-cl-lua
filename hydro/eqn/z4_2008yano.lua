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
Z4_2008Yano.name = 'z4_2008yano'

-- TODO keep this as 'true' and instead implement 'fluxFromCons'
-- because with this disabled the roe flux assumes dF/dU * U = F (right?)
Z4_2008Yano.roeUseFluxFromCons = false

function Z4_2008Yano:init(args)
	
	self.useShift = args.useShift or 'none'
	
	local fluxVars = table{
		{name='a_l', type='real3'},		-- 3:  0-2
		{name='d_lll', type='_3sym3'},	-- 18: 3-20
		{name='K_ll', type='sym3'},		-- 6:  21-26
		{name='Theta', type='real'},	-- 1:  27
		{name='Z_l', type='real3'},		-- 3:  28-30
	}

	self.consVars = table{
		{name='alpha', type='real'},
		{name='gamma_ll', type='sym3'},
	}:append(fluxVars)
	
	self:cdefAllVarTypes(args.solver, self.consVars)	-- have to call before countScalars in eqn:init

	self.numWaves = Struct.countScalars{vars=fluxVars}
	assert(self.numWaves == 31)
		
	self.numIntStates = Struct.countScalars{vars=self.consVars}

	self.consVars:append{
		--constraints:              
		{name='H', type='real'},				--1
		{name='M_u', type='real3'},				--3
	}

	Z4_2008Yano.super.init(self, args)
end

function Z4_2008Yano:createInitState()
	Z4_2008Yano.super.createInitState(self)
	self:addGuiVar{name = 'm', value = 2}
end

-- don't use default
function Z4_2008Yano:initCodeModule_calcDTCell() end

Z4_2008Yano.solverCodeFile = 'hydro/eqn/z4_2008yano.cl'

Z4_2008Yano.predefinedDisplayVars = {
	'U alpha',
	'U gamma_ll x x',
	'U a_l x',
	'U d_lll_x x x',
	'U K_ll x x',
	'U K_ll tr weighted gamma^ij',
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

function Z4_2008Yano:eigenWaveCodePrefix(args)
	return template([[
real const eig_lambdaLight = (<?=eig?>)->alpha * (<?=eig?>)->sqrt_gammaUjj.s[(<?=n?>).side];
real const eig_lambdaGauge = eig_lambdaLight * (<?=eig?>)->sqrt_f;
]], args)
end

function Z4_2008Yano:eigenWaveCode(args)
	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..args.eig..')->beta_u.s[('..args.n..').side]'
	else
		betaUi = '0'
	end

	local waveIndex = args.waveIndex
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

function Z4_2008Yano:consWaveCodePrefix(args)
	return template([[
real const det_gamma = sym3_det((<?=U?>)->gamma_ll);
sym3 const gamma_uu = sym3_inv((<?=U?>)->gamma_ll, det_gamma);
real eig_lambdaLight = 0./0.;
if ((<?=n?>).side == 0) {
	eig_lambdaLight = (<?=U?>)->alpha * sqrt(gamma_uu.xx);
} else if ((<?=n?>).side == 1) {
	eig_lambdaLight = (<?=U?>)->alpha * sqrt(gamma_uu.yy);
} else if ((<?=n?>).side == 2) {
	eig_lambdaLight = (<?=U?>)->alpha * sqrt(gamma_uu.zz);
}
real const f = calc_f((<?=U?>)->alpha);
real eig_lambdaGauge = eig_lambdaLight * sqrt(f);
]], args)
end
Z4_2008Yano.consWaveCode = Z4_2008Yano.eigenWaveCode

return Z4_2008Yano
