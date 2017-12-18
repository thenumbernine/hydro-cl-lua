--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local EinsteinEqn = require 'eqn.einstein'
local symmath = require 'symmath'
local makeStruct = require 'eqn.makestruct'

-- TODO assign these as locals instead of globals
require 'common'(_G)

local ADM_BonaMasso_3D = class(EinsteinEqn)
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


-- no shift
--ADM_BonaMasso_3D.useShift = false

-- minimal distortion elliptic -- Alcubierre's book, eqn 4.3.14 and 4.3.15
--ADM_BonaMasso_3D.useShift = 'MinimalDistortionElliptic'

-- 2008 Alcubierre 4.3.37
-- I see some problems in the warp bubble test ...
ADM_BonaMasso_3D.useShift = 'HarmonicShiftCondition-FiniteDifference'


if ADM_BonaMasso_3D.useShift == 'MinimalDistortionElliptic' then
	ADM_BonaMasso_3D.consVars:insert{beta_u = 'real3'}
	--[[ and maybe some of these ...
	ADM_BonaMasso_3D.consVars:insert{gamma_uu = 'sym3'}
	ADM_BonaMasso_3D.consVars:insert{conn_ull = '_3sym3'}
	ADM_BonaMasso_3D.consVars:insert{R_ll = 'sym3'}
	--]]

elseif ADM_BonaMasso_3D.useShift == 'HarmonicShiftCondition-FiniteDifference' then
	ADM_BonaMasso_3D.consVars:insert{beta_u = 'real3'}
end


--[[
solve a smaller eigendecomposition that doesn't include the rows of variables whose d/dt is zero.
kind of like how ideal MHD is an 8-var system, but the flux jacobian solved is a 7x7 because Bx,t = 0.
TODO make this a ctor arg - so solvers can run in parallel with and without this
...or why don't I just scrap the old code, because this runs a lot faster.
--]]
ADM_BonaMasso_3D.noZeroRowsInFlux = true

-- NOTE this doesn't work when using shift ... because then all the eigenvalues are -beta^i, so none of them are zero (except the source-only alpha, beta^i, gamma_ij)
if ADM_BonaMasso_3D.useShift then
	ADM_BonaMasso_3D.noZeroRowsInFlux = false
end

if not ADM_BonaMasso_3D.noZeroRowsInFlux then
	-- skip alpha and gamma
	ADM_BonaMasso_3D.numWaves = makeStruct.countReals(fluxVars)
	assert(ADM_BonaMasso_3D.numWaves == 30)
else
	-- skip alpha, gamma, a_q, d_qij, V_i for q != the direction of flux
	ADM_BonaMasso_3D.numWaves = 13
end



-- only count int vars after the shifts have been added
ADM_BonaMasso_3D.numIntStates = makeStruct.countReals(ADM_BonaMasso_3D.consVars)



ADM_BonaMasso_3D.hasCalcDT = true
ADM_BonaMasso_3D.hasEigenCode = true
ADM_BonaMasso_3D.useSourceTerm = true
ADM_BonaMasso_3D.useConstrainU = true

function ADM_BonaMasso_3D:init(solver)
	ADM_BonaMasso_3D.super.init(self, solver)

	if self.useShift == 'MinimalDistortionElliptic' then
		local MinimalDistortionEllipticShift = require 'solver.gr-shift-mde'
		solver.ops:insert(MinimalDistortionEllipticShift{solver=solver})
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
			}
		},
		{name='a_convCoeff', value=10},
		{name='d_convCoeff', value=10},
		{name='V_convCoeff', value=10},
	}
	-- TODO add shift option
	-- but that means moving the consVars construction to the :init()
end

function ADM_BonaMasso_3D:getCodePrefix()
	return table{
		ADM_BonaMasso_3D.super.getCodePrefix(self),
		template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1.;
	U->gamma = _sym3(1,0,0,1,0,1);
	U->a = _real3(0,0,0);
	U->d.x = _sym3(0,0,0,0,0,0);
	U->d.y = _sym3(0,0,0,0,0,0);
	U->d.z = _sym3(0,0,0,0,0,0);
	U->K = _sym3(0,0,0,0,0,0);
	U->V = _real3(0,0,0);
<? if eqn.useShift then 
?>	U->beta_u = _real3(0,0,0);
<? end 
?>
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
	U->V = _real3(0,0,0);
<? if eqn.useShift then
?>	U->beta_u = beta_u;
<? end
?>
}

kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a.<?=xi?> = (U[stepsize.<?=xi?>].alpha - U[-stepsize.<?=xi?>].alpha) / (grid_dx<?=i-1?> * U->alpha);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d.<?=xi?>.<?=xjk?> = .5 * (U[stepsize.<?=xi?>].gamma.<?=xjk?> - U[-stepsize.<?=xi?>].gamma.<?=xjk?>) / grid_dx<?=i-1?>;
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a.<?=xi?> = 0;
	U->d.<?=xi?> = _sym3(0,0,0,0,0,0);
<?
end
?>

//V_i = d_ik^k - d^k_ki 
<? for i,xi in ipairs(xNames) do ?>
	U->V.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + gammaU.<?=sym(j,k)?> * ( U->d.<?=xi?>.<?=sym(j,k)?> - U->d.<?=xj?>.<?=sym(k,i)?> )<?
		end
	end ?>;
<? end ?>
}
]]

function ADM_BonaMasso_3D:getSolverCode()
	local derivOrder = 2 * self.solver.numGhost
	return template(file['eqn/adm3d.cl'], {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
		sym = sym,
		makePartial = function(...) return require 'eqn.makepartial'.makePartial(derivOrder, self.solver, ...) end,
		makePartial2 = function(...) return require 'eqn.makepartial'.makePartial2(derivOrder, self.solver, ...) end,
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

	vars:insert{['alpha vs a_i'] = template([[
	if (OOB(1,1)) {
		*valuevec = _real3(0,0,0);
	} else {
		<? for i=1,solver.dim do
			local xi = xNames[i]
		?>{
			real di_alpha = (U[stepsize.<?=xi?>].alpha - U[-stepsize.<?=xi?>].alpha) / (2. * grid_dx<?=i-1?>);
			valuevec-><?=xi?> = fabs(di_alpha - U->alpha * U->a.<?=xi?>);
		}<? end ?>
		<? for i=solver.dim+1,3 do
			local xi = xNames[i]
		?>{
			valuevec-><?=xi?> = 0;
		}<? end ?>
	}
]], {
	solver = self.solver,
	xNames = xNames,
}), type='real3'}

	-- d_kij = gamma_ij,k
	for i,xi in ipairs(xNames) do
		vars:insert{['gamma_ij vs d_'..xi..'ij'] = template([[
	if (OOB(1,1)) {
		*valuesym3 = (sym3){.s={0,0,0,0,0,0}};
	} else {
		<? if i <= solver.dim then ?>
		sym3 di_gamma_jk = sym3_scale(
			sym3_sub(
				U[stepsize.<?=xi?>].gamma, 
				U[-stepsize.<?=xi?>].gamma
			), 
			1. / (2. * grid_dx<?=i-1?>)
		);
		<? else ?>
		sym3 di_gamma_jk = _sym3(0,0,0,0,0,0);
		<? end ?>
		*valuesym3 = sym3_sub(di_gamma_jk, sym3_scale(U->d.<?=xi?>, 2.));
		*valuesym3 = (sym3){<?
	for jk,xjk in ipairs(symNames) do 
?>			.<?=xjk?> = fabs(valuesym3-><?=xjk?>),
<?	end
?>		};
	}
]], {
	i = i,
	xi = xi,
	xNames = xNames,
	symNames = symNames,
	solver = self.solver,
}), type='sym3'}
	end

	vars:insert{['V constraint'] = template([[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	<? for i,xi in ipairs(xNames) do ?>{
		real d1 = sym3_dot(U->d.<?=xi?>, gammaU);
		real d2 = 0.<?
	for j=1,3 do
		for k,xk in ipairs(xNames) do
?> + U->d.<?=xj?>.<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
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
if ADM_BonaMasso_3D.useShift then
	ADM_BonaMasso_3D.eigenVars:insert{beta_u = 'real3'}
end

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
