--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local NumRelEqn = require 'eqn.numrel'
local symmath = require 'symmath'
local makeStruct = require 'eqn.makestruct'

-- TODO assign these as locals instead of globals
require 'common'(_G)

local ADM_BonaMasso_3D = class(NumRelEqn)
ADM_BonaMasso_3D.name = 'ADM_BonaMasso_3D'

local fluxVars = table{
	{a = 'real3'},
	{d = '_3sym3'},
	{K = 'sym3'},
	{V = 'real3'},
}

ADM_BonaMasso_3D.consVars = table{
	{alpha = 'real'},
	{gamma = 'sym3'},
}:append(fluxVars)


--[[
solve a smaller eigendecomposition that doesn't include the rows of variables whose d/dt is zero.
kind of like how ideal MHD is an 8-var system, but the flux jacobian solved is a 7x7 because Bx,t = 0
TODO make this a ctor arg - so solvers can run in parallel with and without this
...or why don't I just scrap the old code, because this runs a lot faster.
--]]
ADM_BonaMasso_3D.noZeroRowsInFlux = true

if not ADM_BonaMasso_3D.noZeroRowsInFlux then
	-- skip alpha and gamma
	ADM_BonaMasso_3D.numWaves = makeStruct.countReals(fluxVars)
	assert(ADM_BonaMasso_3D.numWaves == 30)
else
	-- skip alpha, gamma, a_q, d_qij, V_i for q != the direction of flux
	ADM_BonaMasso_3D.numWaves = 13
end

ADM_BonaMasso_3D.hasCalcDT = true
ADM_BonaMasso_3D.hasEigenCode = true
ADM_BonaMasso_3D.useSourceTerm = true
ADM_BonaMasso_3D.useConstrainU = true

function ADM_BonaMasso_3D:createInitState()
	ADM_BonaMasso_3D.super.createInitState(self)
	self:addGuiVar{
		type = 'combo',
		name = 'constrain V',
		options = {
			'none',
			'replace V',
			'average',	-- TODO add averaging weights, from 100% V (which works) to 100% d (which doesn't yet work)
		}
	}
end

function ADM_BonaMasso_3D:getCodePrefix()
	return table{
		ADM_BonaMasso_3D.super.getCodePrefix(self),
		template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1.;
	U->gamma = _sym3(1,0,0,1,0,1);
	U->a = _real3(0,0,0);
	U->d[0] = _sym3(0,0,0,0,0,0);
	U->d[1] = _sym3(0,0,0,0,0,0);
	U->d[2] = _sym3(0,0,0,0,0,0);
	U->K = _sym3(0,0,0,0,0,0);
	U->V = _real3(0,0,0);
}
]], {eqn=self}),
	}:concat'\n'
end

ADM_BonaMasso_3D.initStateCode = [[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = _real3(0,0,0);
	sym3 gamma_ll = _sym3(1,0,0,1,0,1);
	sym3 K_ll = _sym3(0,0,0,0,0,0);

	<?=code?>

	U->alpha = alpha;
	U->gamma = gamma_ll;
	U->K = K_ll;
	
	//TODO V_i = d_ik^k - d^k_ki 
	U->V = _real3(0,0,0);	
}

kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;

<? for i=0,2 do ?>
	U->a.s<?=i?> = (U[stepsize.s<?=i?>].alpha - U[-stepsize.s<?=i?>].alpha) / (grid_dx<?=i?> * U->alpha);
	<? for j=0,2 do ?>
		<? for k=j,2 do ?>
	U->d[<?=i?>].s<?=j..k?> = .5 * (U[stepsize.s<?=i?>].gamma.s<?=j..k?> - U[-stepsize.s<?=i?>].gamma.s<?=j..k?>) / grid_dx<?=i?>;
		<? end ?>
	<? end ?>
<? end ?>
}
]]

function ADM_BonaMasso_3D:getSolverCode()
	return template(file['eqn/adm3d.cl'], {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
		sym = sym,
	})
end

function ADM_BonaMasso_3D:getDisplayVars()
	local vars = ADM_BonaMasso_3D.super.getDisplayVars(self)

	vars:append{
		{det_gamma = '*value = sym3_det(U->gamma);'},
		{volume = '*value = U->alpha * sqrt(sym3_det(U->gamma));'},
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
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
	}:append{
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
	}

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
	vars:insert{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*valuevec = real3_scale(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]], type='real3'}
	
	vars:insert{constraint_V = template([[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d[<?=i-1?>], gammaU);
		real d2 = 0.<?
	for j=1,3 do
		for k,xk in ipairs(xNames) do
?> + U->d[<?=j-1?>].<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
		end
	end ?>;
		valuevec-><?=xi?> = U->V.<?=xi?> - (d1 - d2);
	}<? end ?>
]], {sym=sym, xNames=xNames}), type='real3'}

	return vars
end

ADM_BonaMasso_3D.eigenVars = table{
	{alpha = 'real'},	--used only by eigen_calcWaves ... makes me think eigen_forCell / eigen_forSide should both calculate waves and basis variables in the same go
	{sqrt_f = 'real'},
	{gammaU = 'sym3'},
	-- sqrt(gamma^jj) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
	{sqrt_gammaUjj = 'real3'},
}

function ADM_BonaMasso_3D:fillRandom(epsilon)
	local ptr = ADM_BonaMasso_3D.super.fillRandom(self, epsilon)
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

return ADM_BonaMasso_3D
