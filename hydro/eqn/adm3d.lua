--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'
local Struct = require 'hydro.struct.struct'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local sym = common.sym


local ADM_BonaMasso_3D = class(EinsteinEqn)
ADM_BonaMasso_3D.name = 'ADM_BonaMasso_3D'
ADM_BonaMasso_3D.hasCalcDTCode = true
ADM_BonaMasso_3D.useSourceTerm = true
ADM_BonaMasso_3D.useConstrainU = true

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

	if not self.noZeroRowsInFlux then
		-- skip alpha and gamma
		self.numWaves = Struct.countScalars{vars=fluxVars}
		assert(self.numWaves == 30)
	else
		-- skip alpha, gamma_ij, a_q, d_qij, V_i for q != the direction of flux
		self.numWaves = 13
	end

	-- only count int vars after the shifts have been added
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
		{name='sqrt_f', type='real'},
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
			-- upon setting this, set useConstrainU accordingly (turn it off if we're not constraining U)
			onChange = function(guivar, value, solver)
				-- disable the solver's constrainU if we're not needing it
				-- note that the kernel is still created, just now won't be called
				solver.useConstrainU = value ~= 0
			end,
			compileTime = true,
		},

		-- K_ij source term stress constraint coefficient
		{name='K_ll_srcStressCoeff', value=1},
		
		-- K_ij source term Hamiltonian constraint coefficient
		{name='K_ll_srcHCoeff', value=-.5},

		-- gamma_ij source term stress constraint coefficient
		{name='gamma_ll_srcStressCoeff', value=1},
	
		-- gamma_ij source term Hamiltonian constraint coefficient
		{name='gamma_ll_srcHCoeff', value=-.5},

		-- convergence between finite-difference of alpha,i and alpha a_i
		{name='a_convCoeff', value=0},
		
		-- convergence between finite-difference of 1/2 gamma_ij,k and d_kij
		{name='d_convCoeff', value=0},
		
		-- convergence between V_k and d_kj^j - d^j_jk
		{name='V_convCoeff', value=10},
	}
	-- TODO add shift option
	-- but that means moving the consVars construction to the :init()
end

function ADM_BonaMasso_3D:getCommonFuncCode()
	return template([[
void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->gamma_ll = sym3_ident;
	U->a_l = real3_zero;
	U->d_lll = _3sym3_zero;
	U->K_ll = sym3_zero;
	U->V_l = real3_zero;
<? if eqn.useShift ~= 'none' then 
?>	U->beta_u = real3_zero;
<? end 
?>	
<? if eqn.useStressEnergyTerms then ?>
	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

]], {
		eqn = self,
		solver = self.solver,
	})
end

ADM_BonaMasso_3D.initStateCode = [[
<? 
local common = require 'hydro.common'
local xNames = common.xNames 
local symNames = common.symNames 
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym 
?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = coord_g_ll(x);
	sym3 K_ll = sym3_zero;

	//throw-away for ADM3D ... but not for BSSNOK
	// TODO hold rho somewhere?
	real rho = 0.;

	<?=code?>

	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;
	U->V_l = real3_zero;
<? if eqn.useShift ~= 'none' then
?>	U->beta_u = beta_u;
<? end
?>
<? if eqn.useStressEnergyTerms then ?>
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (solver->grid_dx.s<?=i-1?> * U->alpha);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> - U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>) / solver->grid_dx.s<?=i-1?>;
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>

//V_i = d_ik^k - d^k_ki 
<? for i,xi in ipairs(xNames) do ?>
	U->V_l.<?=xi?> = 0.<?
	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?> + gamma_uu.<?=sym(j,k)?> * ( U->d_lll.<?=xi?>.<?=sym(j,k)?> - U->d_lll.<?=xj?>.<?=sym(k,i)?> )<?
		end
	end ?>;
<? end ?>
}
]]

ADM_BonaMasso_3D.solverCodeFile = 'hydro/eqn/adm3d.cl'

ADM_BonaMasso_3D.predefinedDisplayVars = {
	'U alpha',
	'U gamma_ll x x',
	'U d_lll x x x',
	'U K_ll x x',
	'U V_l x',
	'U H',
	'U volume',
	'U f',
}

function ADM_BonaMasso_3D:getDisplayVars()
	local vars = ADM_BonaMasso_3D.super.getDisplayVars(self)

	vars:append{
		{name='det_gamma', code='value.vreal = sym3_det(U->gamma_ll);'},
		{name='volume', code='value.vreal = U->alpha * sqrt(sym3_det(U->gamma_ll));'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
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

	vars:insert{name='alpha vs a_i', code=template([[
	if (OOB(1,1)) {
		value.vreal3 = real3_zero;
	} else {
		<? for i=1,solver.dim do
			local xi = xNames[i]
		?>{
			real di_alpha = (U[solver->stepsize.<?=xi?>].alpha - U[-solver->stepsize.<?=xi?>].alpha) / (2. * solver->grid_dx.s<?=i-1?>);
			value.vreal3.<?=xi?> = fabs(di_alpha - U->alpha * U->a_l.<?=xi?>);
		}<? end ?>
		<? for i=solver.dim+1,3 do
			local xi = xNames[i]
		?>{
			value.vreal3.<?=xi?> = 0;
		}<? end ?>
	}
]], {
	solver = self.solver,
	xNames = xNames,
}), type='real3'}

	-- d_kij = gamma_ij,k
	for i,xi in ipairs(xNames) do
		vars:insert{name='gamma_ij vs d_'..xi..'ij', code=template([[
	if (OOB(1,1)) {
		value.vsym3 = sym3_zero;
	} else {
		<? if i <= solver.dim then ?>
		sym3 di_gamma_jk = sym3_real_mul(
			sym3_sub(
				U[solver->stepsize.<?=xi?>].gamma_ll, 
				U[-solver->stepsize.<?=xi?>].gamma_ll
			), 
			1. / (2. * solver->grid_dx.s<?=i-1?>)
		);
		<? else ?>
		sym3 di_gamma_jk = sym3_zero;
		<? end ?>
		value.vsym3 = sym3_sub(di_gamma_jk, sym3_real_mul(U->d_lll.<?=xi?>, 2.));
		value.vsym3 = (sym3){<?
	for jk,xjk in ipairs(symNames) do 
?>			.<?=xjk?> = fabs(value.vsym3.<?=xjk?>),
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

	vars:insert{name='V constraint', code=template([[
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
]], {sym=sym, xNames=xNames}), type='real3'}

	return vars
end

function ADM_BonaMasso_3D:eigenWaveCodePrefix(n, eig, x, waveIndex)
	return template([[
	real eig_lambdaLight;
	if (n.side == 0) {
		eig_lambdaLight = <?=eig?>.alpha * <?=eig?>.sqrt_gammaUjj.x;
	} else if (n.side == 1) {
		eig_lambdaLight = <?=eig?>.alpha * <?=eig?>.sqrt_gammaUjj.y;
	} else if (n.side == 2) {
		eig_lambdaLight = <?=eig?>.alpha * <?=eig?>.sqrt_gammaUjj.z;
	}
	real eig_lambdaGauge = eig_lambdaLight * <?=eig?>.sqrt_f;
]], {
		eig = '('..eig..')',
		n = n,
	})
end

function ADM_BonaMasso_3D:eigenWaveCode(n, eig, x, waveIndex)
	-- TODO find out if -- if we use the lagrangian coordinate shift operation -- do we still need to offset the eigenvalues by -beta^i?
	local shiftingLambdas = self.useShift ~= 'none'
		--and self.useShift ~= 'LagrangianCoordinates'

	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..eig..').beta_u.s[n.side]'
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
		-- noZeroRowsInFlux implies useShift == 'none'
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

function ADM_BonaMasso_3D:consWaveCodePrefix(n, U, x, waveIndex)
	return template([[
	real det_gamma = sym3_det(<?=U?>.gamma_ll);
	sym3 gamma_uu = sym3_inv(<?=U?>.gamma_ll, det_gamma);
	real eig_lambdaLight;
	if (n.side == 0) {
		eig_lambdaLight = <?=U?>.alpha * sqrt(gamma_uu.xx);
	} else if (n.side == 1) {
		eig_lambdaLight = <?=U?>.alpha * sqrt(gamma_uu.yy);
	} else if (n.side == 2) {
		eig_lambdaLight = <?=U?>.alpha * sqrt(gamma_uu.zz);
	}
	real f = calc_f(<?=U?>.alpha);
	real eig_lambdaGauge = eig_lambdaLight * sqrt(f);
]], {
		U = '('..U..')',
		n = n,
	})
end
ADM_BonaMasso_3D.consWaveCode = ADM_BonaMasso_3D.eigenWaveCode


function ADM_BonaMasso_3D:fillRandom(epsilon)
	local ptr = ADM_BonaMasso_3D.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.numCells-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma_ll.xx = ptr[i].gamma_ll.xx + 1
		ptr[i].gamma_ll.yy = ptr[i].gamma_ll.yy + 1
		ptr[i].gamma_ll.zz = ptr[i].gamma_ll.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return ADM_BonaMasso_3D
