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

Z4_2008Yano.hasCalcDTCode = true
Z4_2008Yano.useSourceTerm = true

function Z4_2008Yano:createInitState()
	Z4_2008Yano.super.createInitState(self)
	self:addGuiVar{name = 'm', value = -1}
end

function Z4_2008Yano:getCommonFuncCode()
	return template([[
void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1;
	U->gamma_ll = sym3_ident;
	U->a_l = real3_zero;
	U->d_lll.x = sym3_zero;
	U->d_lll.y = sym3_zero;
	U->d_lll.z = sym3_zero;
	U->K_ll = sym3_zero;
	U->Theta = 0;
	U->Z_l = real3_zero;
}
]], {
		eqn = self,
		solver = self.solver,
	})
end

Z4_2008Yano.needsInitDerivs = true
Z4_2008Yano.initCondCode = [[
<? 
local common = require 'hydro.common'
local xNames = common.xNames 
local symNames = common.symNames 
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym 
?>
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	setFlatSpace(solver, U, x);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;

	<?=code?>

	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;
	
	//Z_u n^u = 0
	//Theta = alpha n_u Z^u = alpha Z^u
	//for n_a = (-alpha, 0)
	//n^a_l = (1/alpha, -beta^i/alpha)
	//(Z_t - Z_i beta^i) / alpha = Theta ... = ?
	//Z^t n_t + Z^i n_i = -alpha Z^t = Theta
	U->Theta = 0;
	U->Z_l = real3_zero;
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (solver->grid_dx.s<?=i-1?> * U->alpha);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> - U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>) / solver->grid_dx.s<?=i-1?>;
	<? end ?>
<? end
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>
}
]]

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
		{name='det_gamma', code='value.vreal = sym3_det(U->gamma_ll);'},
		{name='volume', code='value.vreal = U->alpha * sqrt(sym3_det(U->gamma_ll));'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='K_ll', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal = sym3_dot(gamma_uu, U->K_ll);
]]		},
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

function Z4_2008Yano:fillRandom(epsilon)
	local ptr = Z4_2008Yano.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma_ll.xx = ptr[i].gamma_ll.xx + 1
		ptr[i].gamma_ll.yy = ptr[i].gamma_ll.yy + 1
		ptr[i].gamma_ll.zz = ptr[i].gamma_ll.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return Z4_2008Yano
