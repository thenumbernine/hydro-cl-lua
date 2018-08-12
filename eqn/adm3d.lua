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

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local sym = common.sym


local ADM_BonaMasso_3D = class(EinsteinEqn)
ADM_BonaMasso_3D.name = 'ADM_BonaMasso_3D'
ADM_BonaMasso_3D.hasCalcDTCode = true
ADM_BonaMasso_3D.hasEigenCode = true
ADM_BonaMasso_3D.useSourceTerm = true
ADM_BonaMasso_3D.useConstrainU = true

--[[
args:
useShift
	
	useShift = false	--  no shift
	
	useShift = 'MinimalDistortionElliptic' -- minimal distortion elliptic.  Alcubierre's book, eqn 4.3.14 and 4.3.15

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

	local fluxVars = table{
		{a = 'real3'},
		{d = '_3sym3'},
		{K = 'sym3'},
		{V = 'real3'},
	}

	self.consVars = table{
		{alpha = 'real'},
		{gamma = 'sym3'},
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

	self.useShift = args.useShift

	if self.useShift then
		self.consVars:insert{beta_u = 'real3'}

		--[[ and maybe some of these ...
		if self.useShift == 'MinimalDistortionElliptic' then
			self.consVars:insert{gamma_uu = 'sym3'}
			self.consVars:insert{conn_ull = '_3sym3'}
			self.consVars:insert{R_ll = 'sym3'}
		end
		--]]
	end


	--[[
	solve a smaller eigendecomposition that doesn't include the rows of variables whose d/dt is zero.
	kind of like how ideal MHD is an 8-var system, but the flux jacobian solved is a 7x7 because Bx,t = 0.
	TODO make this a ctor arg - so solvers can run in parallel with and without this
	...or why don't I just scrap the old code, because this runs a lot faster.
	--]]
	self.noZeroRowsInFlux = true

	-- NOTE this doesn't work when using shift ... because then all the eigenvalues are -beta^i, so none of them are zero (except the source-only alpha, beta^i, gamma_ij)
	-- with the exception of the lagrangian shift.  if we split that operator out then we can first solve the non-shifted system, then second advect it by the shift vector ...
	--if self.useShift then
	--	self.noZeroRowsInFlux = false
	--end

	if not self.noZeroRowsInFlux then
		-- skip alpha and gamma
		self.numWaves = makeStruct.countScalars(fluxVars)
		assert(self.numWaves == 30)
	else
		-- skip alpha, gamma, a_q, d_qij, V_i for q != the direction of flux
		self.numWaves = 13
	end

	-- only count int vars after the shifts have been added
	self.numIntStates = makeStruct.countScalars(self.consVars)


	self.eigenVars = table{
		{alpha = 'real'},
		{sqrt_f = 'real'},
		{gammaU = 'sym3'},
		-- sqrt(gamma^jj) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
		{sqrt_gammaUjj = 'real3'},
	}

	-- hmm, only certain shift methods actually use beta_u ...
	if self.useShift then
		self.eigenVars:insert{beta_u = 'real3'}
	end


	
	-- build stuff around consVars	
	ADM_BonaMasso_3D.super.init(self, args)


	if self.useShift == 'MinimalDistortionElliptic' then
		local MinimalDistortionEllipticShift = require 'op.gr-shift-mde'
		self.solver.ops:insert(MinimalDistortionEllipticShift{solver=self.solver})
	elseif self.useShift == 'LagrangianCoordinates' then
		local LagrangianCoordinateShift = require 'op.gr-shift-lc'
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
			-- upon setting this, set useConstrainU accordingly (turn it off if we're not constraining U)
			onChange = function(guivar, value, solver)
				-- disable the solver's constrainU if we're not needing it
				-- note that the kernel is still created, just now won't be called
				solver.useConstrainU = value ~= 0
			end,
		},
		{name='a_convCoeff', value=10},
		{name='d_convCoeff', value=10},
		{name='V_convCoeff', value=10},
	}
	-- TODO add shift option
	-- but that means moving the consVars construction to the :init()
end

function ADM_BonaMasso_3D:getCommonFuncCode()
	return template([[
void setFlatSpace(global <?=eqn.cons_t?>* U, real3 x) {
	U->alpha = 1.;
	U->gamma = sym3_ident;
	U->a = real3_zero;
	U->d.x = sym3_zero;
	U->d.y = sym3_zero;
	U->d.z = sym3_zero;
	U->K = sym3_zero;
	U->V = real3_zero;
<? if eqn.useShift then 
?>	U->beta_u = real3_zero;
<? end 
?>
}
]], {eqn=self})
end

ADM_BonaMasso_3D.initStateCode = [[
<? 
local common = require 'common'()
local xNames = common.xNames 
local symNames = common.symNames 
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym 
?>
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(mins, maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;

	//throw-away for ADM3D ... but not for BSSNOK
	// TODO hold rho somewhere?
	real rho = 0.;

	<?=code?>

	U->alpha = alpha;
	U->gamma = gamma_ll;
	U->K = K_ll;
	U->V = real3_zero;
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
	U->d.<?=xi?> = sym3_zero;
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

ADM_BonaMasso_3D.solverCodeFile = 'eqn/adm3d.cl'

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
	*value = -sym3_dot(gammaU, U->K);
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
	vars:insert{gravity = [[
	real det_gamma = sym3_det(U->gamma);
	sym3 gammaU = sym3_inv(U->gamma, det_gamma);
	*value_real3 = real3_real_mul(sym3_real3_mul(gammaU, U->a), -U->alpha * U->alpha);
]], type='real3'}

	vars:insert{['alpha vs a_i'] = template([[
	if (OOB(1,1)) {
		*value_real3 = real3_zero;
	} else {
		<? for i=1,solver.dim do
			local xi = xNames[i]
		?>{
			real di_alpha = (U[stepsize.<?=xi?>].alpha - U[-stepsize.<?=xi?>].alpha) / (2. * grid_dx<?=i-1?>);
			value_real3-><?=xi?> = fabs(di_alpha - U->alpha * U->a.<?=xi?>);
		}<? end ?>
		<? for i=solver.dim+1,3 do
			local xi = xNames[i]
		?>{
			value_real3-><?=xi?> = 0;
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
		sym3 di_gamma_jk = sym3_real_mul(
			sym3_sub(
				U[stepsize.<?=xi?>].gamma, 
				U[-stepsize.<?=xi?>].gamma
			), 
			1. / (2. * grid_dx<?=i-1?>)
		);
		<? else ?>
		sym3 di_gamma_jk = sym3_zero;
		<? end ?>
		*valuesym3 = sym3_sub(di_gamma_jk, sym3_real_mul(U->d.<?=xi?>, 2.));
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
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + U->d.<?=xj?>.<?=sym(k,i)?> * gammaU.<?=sym(j,k)?><?
		end
	end ?>;
		value_real3-><?=xi?> = U->V.<?=xi?> - (d1 - d2);
	}<? end ?>
]], {sym=sym, xNames=xNames}), type='real3'}

	return vars
end

function ADM_BonaMasso_3D:eigenWaveCodePrefix(side, eig, x, waveIndex)
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

function ADM_BonaMasso_3D:eigenWaveCode(side, eig, x, waveIndex)
	-- TODO find out if -- if we use the lagrangian coordinate shift operation -- do we still need to offset the eigenvalues by -beta^i?
	local shiftingLambdas = self.useShift 
		--and self.useShift ~= 'LagrangianCoordinates'

	local betaUi
	if self.useShift then
		betaUi = eig..'.beta_u.'..xNames[side+1]
	else
		betaUi = '0'
	end


	if not self.noZeroRowsInFlux then
		if waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif waveIndex >= 1 and waveIndex <= 5 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif waveIndex >= 6 and waveIndex <= 23 then
			return '-'..betaUi
		elseif waveIndex >= 24 and waveIndex <= 28 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif waveIndex == 29 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end

	else	-- noZeroRowsInFlux 
		-- noZeroRowsInFlux implies not useShift
		if waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif waveIndex >= 1 and waveIndex <= 5 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif waveIndex == 6 then
			return '-'..betaUi
		elseif waveIndex >= 7 and waveIndex <= 11 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif waveIndex == 12 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	end
	error'got a bad waveIndex'
end

function ADM_BonaMasso_3D:consWaveCodePrefix(side, U, x, waveIndex)
	return template([[
	real det_gamma = sym3_det(<?=U?>.gamma);
	sym3 gammaU = sym3_inv(<?=U?>.gamma, det_gamma);
	<? if side==0 then ?>
	real eig_lambdaLight = <?=U?>.alpha * sqrt(gammaU.xx);
	<? elseif side==1 then ?>                          
	real eig_lambdaLight = <?=U?>.alpha * sqrt(gammaU.yy);
	<? elseif side==2 then ?>                          
	real eig_lambdaLight = <?=U?>.alpha * sqrt(gammaU.zz);
	<? end ?>
	real f = calc_f(<?=U?>.alpha);
	real eig_lambdaGauge = eig_lambdaLight * sqrt(f);
]], {
		U = '('..U..')',
		side = side,
	})
end
ADM_BonaMasso_3D.consWaveCode = ADM_BonaMasso_3D.eigenWaveCode


function ADM_BonaMasso_3D:fillRandom(epsilon)
	local ptr = ADM_BonaMasso_3D.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.numCells-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma.xx = ptr[i].gamma.xx + 1
		ptr[i].gamma.yy = ptr[i].gamma.yy + 1
		ptr[i].gamma.zz = ptr[i].gamma.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return ADM_BonaMasso_3D
