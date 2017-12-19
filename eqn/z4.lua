--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
With some hints from a few of the Z4 papers.
(Which papers again?  The 2008 Alcubierre book doesn't have eigenvectors, and I got them from one of the papers, but the papers don't all have the same eigensystem (my favorite difference being whether sqrt(gamma^xx) shows up in the eigenvalues or not)) 
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local EinsteinEqn = require 'eqn.einstein'
local symmath = require 'symmath'
local makeStruct = require 'eqn.makestruct'
require 'common'(_G)

local Z4 = class(EinsteinEqn)
Z4.name = 'Z4'

local fluxVars = table{
	{a = 'real3'},
	{d = '_3sym3'},
	{K = 'sym3'},
	{Theta = 'real'},
	{Z = 'real3'},
}

Z4.consVars = table{
	{alpha = 'real'},
	{gamma = 'sym3'},
}:append(fluxVars)

Z4.numWaves = makeStruct.countReals(fluxVars)
assert(Z4.numWaves == 31)

Z4.hasCalcDT = true
Z4.hasEigenCode = true
Z4.useSourceTerm = true

function Z4:createInitState()
	Z4.super.createInitState(self)
	self:addGuiVar{name = 'lambda', value = -1}
end

function Z4:getCodePrefix()
	return table{
		Z4.super.getCodePrefix(self),
		template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1;
	U->gamma = _sym3(1,0,0,1,0,1);
	U->a = _real3(0,0,0);
	U->d.x = _sym3(0,0,0,0,0,0);
	U->d.y = _sym3(0,0,0,0,0,0);
	U->d.z = _sym3(0,0,0,0,0,0);
	U->K = _sym3(0,0,0,0,0,0);
	U->Theta = 0;
	U->Z = _real3(0,0,0);
}
]], {eqn=self}),
	}:concat()
end

Z4.initStateCode = [[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	setFlatSpace(U);

	real alpha = 1.;
	real3 beta_u = _real3(0,0,0);
	sym3 gamma_ll = _sym3(1,0,0,1,0,1);
	sym3 K_ll = _sym3(0,0,0,0,0,0);

	<?=code?>

	U->alpha = alpha;
	U->gamma = gamma_ll;
	U->K = K_ll;
	
	//Z_u n^u = 0
	//Theta = alpha n_u Z^u = alpha Z^u
	//for n_a = (-alpha, 0)
	//n^a = (1/alpha, -beta^i/alpha)
	//(Z_t - Z_i beta^i) / alpha = Theta ... = ?
	//Z^t n_t + Z^i n_i = -alpha Z^t = Theta
	U->Theta = 0;
	U->Z = _real3(0,0,0);
}

kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;

<? for i,xi in ipairs(xNames) do ?>
	U->a.<?=xi?> = (U[stepsize.<?=xi?>].alpha - U[-stepsize.<?=xi?>].alpha) / (grid_dx<?=i-1?> * U->alpha);
	<? for j=0,2 do ?>
		<? for k=j,2 do ?>
	U->d.<?=xi?>.s<?=j..k?> = .5 * (U[stepsize.<?=xi?>].gamma.s<?=j..k?> - U[-stepsize.<?=xi?>].gamma.s<?=j..k?>) / grid_dx<?=i-1?>;
		<? end ?>
	<? end ?>
<? end ?>
}
]]

function Z4:getSolverCode()
	return template(file['eqn/z4.cl'], {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
		sym = sym,
	})
end

function Z4:getDisplayVars()
	local vars = Z4.super.getDisplayVars(self)
	vars:append{
		{det_gamma = '*value = sym3_det(U->gamma);'},
		{volume = '*value = U->alpha * sqrt(sym3_det(U->gamma));'},
		{f = '*value = calc_f(U->alpha);'},
		{K = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = sym3_dot(gammaU, U->K);
]]		},
		{expansion = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value = -U->alpha * sym3_dot(gammaU, U->K);
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
		{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*valuevec = real3_scale(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]], type='real3'},
	}
	
	return vars
end

Z4.eigenVars = table{
	{alpha = 'real'},
	{sqrt_f = 'real'},
	{gamma = 'sym3'},
	{gammaU = 'sym3'},
	{sqrt_gammaUjj = 'real3'},
}

function Z4:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	<? if side==0 then ?>
	real eig_lambdaLight = <?=eig?>->alpha * <?=eig?>->sqrt_gammaUjj.x;
	<? elseif side==1 then ?>
	real eig_lambdaLight = <?=eig?>->alpha * <?=eig?>->sqrt_gammaUjj.y;
	<? elseif side==2 then ?>
	real eig_lambdaLight = <?=eig?>->alpha * <?=eig?>->sqrt_gammaUjj.z;
	<? end ?>
	real eig_lambdaGauge = eig_lambdaLight * <?=eig?>->sqrt_f;
]], {
		eig = '('..eig..')',
		side = side,
	})
end

function Z4:eigenWaveCode(side, eig, x, waveIndex)

	local betaUi
	if self.useShift then
		betaUi = eig..'->beta_u.'..xNames[side+1]
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



function Z4:fillRandom(epsilon)
	local ptr = Z4.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma.xx = ptr[i].gamma.xx + 1
		ptr[i].gamma.yy = ptr[i].gamma.yy + 1
		ptr[i].gamma.zz = ptr[i].gamma.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return Z4 
