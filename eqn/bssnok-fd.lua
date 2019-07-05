--[[
Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer" 2010
Alcubierre "Introduction to Numerical Relativity" 2008

then I'm applying 2017 Ruchlin changes...
*) separate gammaBar_ll - gammaHat_ll = epsilon_ll
*) coordinate-transform beta^i, epsilon_ij, ABar_ij, LambdaBar^i to eliminate singularities from the metric

tensors are denoted with suffixes _u _l etc for upper and lower
rescaled tensors are denoted _U _L etc
--]]

local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'eqn.einstein'
local makestruct = require 'eqn.makestruct'
local applyCommon = require 'common'

local makePartials = require 'eqn.makepartial'
local makePartial = makePartials.makePartial
local makePartial2 = makePartials.makePartial2

local BSSNOKFiniteDifferenceEquation = class(EinsteinEqn)
BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 
BSSNOKFiniteDifferenceEquation.hasEigenCode = true
BSSNOKFiniteDifferenceEquation.hasCalcDTCode = true
BSSNOKFiniteDifferenceEquation.hasFluxFromConsCode = true
BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true

-- not used with finite-difference schemes anyways
BSSNOKFiniteDifferenceEquation.weightFluxByGridVolume = false

--[[
args:
	useShift = 'none'
				'GammaDriver'
				'HyperbolicGammaDriver' (default)
--]]
function BSSNOKFiniteDifferenceEquation:init(args)
	-- options:
	-- needs to be defined up front
	-- otherwise rebuild intVars based on it ...
	self.useShift = args.useShift or 'HyperbolicGammaDriver'

	local intVars = table{
		{name='alpha', type='real'},			-- 1
		{name='beta_U', type='real3'},		 	-- 3: beta^i
		{name='epsilon_LL', type='sym3'},		-- 6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='W', type='real'},				-- 1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},				-- 1: K = K^i_i
		{name='ABar_LL', type='sym3'},			-- 6: ABar_ij, only 5 dof since ABar^k_k = 0
		{name='LambdaBar_U', type='real3'},		-- 3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
												-- TODO what is C^i ?
	}
	if self.useShift == 'HyperbolicGammaDriver' then
		intVars:insert{name='B_U', type='real3'}
	end

	self.consVars = table()
	:append(intVars)
	:append{
		--hyperbolic variables:
		--real3 a;				//3: a_i
		--_3sym3 dBar;			//18: dBar_ijk, only 15 dof since dBar_ij^j = 0
		--real3 Phi;			//3: Phi_i

		--stress-energy variables:
		{name='rho', type='real'},				--1: n_a n_b T^ab
		{name='S_u', type='real3'},				--3: -gamma^ij n_a T_aj
		{name='S_ll', type='sym3'},				--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},				--1
		{name='M_u', type='real3'},				--3
	}
	self.numIntStates = makestruct.countScalars(intVars)
	
	-- call construction / build structures	
	BSSNOKFiniteDifferenceEquation.super.init(self, args)
end

function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaBar', value=true, compileTime=true},
		--{name='constrain_det_gammaBar', value=false, compileTime=true},

		-- turning this off until I can diagnose where the ABar_ij trace derivatives are coming from in my spherical vacuum spacetime
		--{name='constrain_tr_ABar', value=true, compileTime=true},
		{name='constrain_tr_ABar', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},
	}
end

function BSSNOKFiniteDifferenceEquation:makePartial(field, fieldType, nameOverride, rescale)
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	if rescale == nil then rescale = true end
	local derivOrder = 2 * self.solver.numGhost
	return makePartial(derivOrder, self.solver, field, fieldType, nameOverride, rescale)
end

function BSSNOKFiniteDifferenceEquation:makePartial2(field, fieldType, nameOverride, rescale)
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	if rescale == nil then rescale = true end
	local derivOrder = 2 * self.solver.numGhost
	return makePartial2(derivOrder, self.solver, field, fieldType, nameOverride, rescale)
end


function BSSNOKFiniteDifferenceEquation:getTemplateEnv()
	local derivOrder = 2 * self.solver.numGhost
	return applyCommon{
		eqn = self,
		solver = self.solver,
	}
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template([[

//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero

//gammaBar_ij = gammaHat_ij + epsilon_ij
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_ll = coord_g_ll(x);
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);
	return gammaBar_ll;
}

real calc_det_gammaHat(real3 x) {
	real det_gammaHat = coord_det_g(x);
	return det_gammaHat;
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2017 Ruchlin
real calc_det_gammaBar(real3 x) {
	//TODO detg ...
	real detg = 1.;
	real det_gammaHat = calc_det_gammaHat(x);
	real det_gammaBar = det_gammaHat * detg;
	return det_gammaBar;
}

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	real det_gamma = calc_det_gammaBar(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma); 
	return gamma_uu;
}

void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->beta_U = real3_zero;
	U->epsilon_LL = sym3_zero;
	U->W = 1;
	U->K = 0;
	U->ABar_LL = sym3_zero;

	//LambdaBar^i = Delta^i + C^i = Delta^i_jk gammaBar^jk = (connBar^i_jk - connHat^i_jk) gammaBar^jk + C^i
	//but when space is flat we have connBar^i_jk = connHat^i_jk and therefore Delta^i_jk = 0, Delta^i = 0, and LambdaBar^i = 0
	U->LambdaBar_U = mystery_C_U;

<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_U = real3_zero;
<? end
?>

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_u = real3_zero;
}

]], self:getTemplateEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_connBar_ull()
	return template([[
	_3sym3 connBar_ull;
	{
		//connBar_lll.i.jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
		_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
?>		connBar_lll.<?=xi?>.<?=xjk?> = .5 * (0.
			+ partial_gammaBar_lll.<?=xk?>.<?=sym(i,j)?>
			+ partial_gammaBar_lll.<?=xj?>.<?=sym(i,k)?>
			- partial_gammaBar_lll.<?=xi?>.<?=sym(j,k)?>
		);
<?	end
end
?>		//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
		connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);
	}
]], self:getTemplateEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_RBar_ll()
	return template([[

	//partial2_gammaBar_llll.kl.ij := gammaBar_ij,kl
	sym3sym3 partial2_gammaBar_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
?>	partial2_gammaBar_llll.<?=xkl?>.<?=xij?> = partial2_epsilon_llll[<?=kl-1?>].<?=xij?> + partial2_gammaHat_llll.<?=xkl?>.<?=xij?>;
<?	end
end 
?>

	//DHat_gammaBar_lll.k.ij = DHat_k gammaBar_ij 
	// = gammaBar_ij,k - connHat^l_ki gammaBar_lj - connHat^l_kj gammaBar_il
	_3sym3 DHat_gammaBar_lll;
<?
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i], xNames[j]
?>	DHat_gammaBar_lll.<?=xk?>.<?=xij?> = partial_gammaBar_lll.<?=xk?>.<?=xij?>
<?		for l,xl in ipairs(xNames) do
?>		- connHat_ull.<?=xl?>.<?=sym(k,i)?> * gammaBar_ll.<?=sym(l,j)?>
		- connHat_ull.<?=xl?>.<?=sym(k,j)?> * gammaBar_ll.<?=sym(l,i)?>
<?		end
?>	;
<?	end
end
?>

	/*
	partial_DHat_gammaBar_llll[l].k.ij := partial_l DHat_k gammaBar_ij
		= (gammaBar_ij,k - connHat^m_ki gammaBar_mj - connHat^m_kj gammaBar_mi)_,l
		= gammaBar_ij,kl 
			- connHat^m_ki,l gammaBar_mj 
			- connHat^m_kj,l gammaBar_mi
			- connHat^m_ki gammaBar_mj,l
			- connHat^m_kj gammaBar_mi,l
	*/
	_3sym3 partial_DHat_gammaBar_llll[3];
<? 
for k,xk in ipairs(xNames) do
	for l,xl in ipairs(xNames) do
		for ij,xij in ipairs(symNames) do
			local i,j = from6to3x3(ij)
			local xi,xj = xNames[i], xNames[j]
?>	partial_DHat_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
		+ partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>
		- partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,i)?> * gammaBar_ll.<?=sym(m,j)?>		//diverging: RBar_ij = diag(0, -.5, -.75)
		- partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,j)?> * gammaBar_ll.<?=sym(m,i)?>		//diverging: RBar_ij = diag(0, -.5, -.75)
		- connHat_ull.<?=xm?>.<?=sym(k,i)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,j)?>		//diverging: RBar_ij = diag(0, 1, 1)
		- connHat_ull.<?=xm?>.<?=sym(k,j)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,i)?>		//diverging: RBar_ij = diag(0, 1, 1)
<?			end
?>	;
<?		end
	end
end
?>

	/*
	DHat2_gammaBar_llll[l].k.ij = DHat_l DHat_k gammaBar_ij
		= partial_l DHat_k gammaBar_ij
			- connHat^m_lk DHat_m gammaBar_ij
			- connHat^m_li DHat_k gammaBar_mj
			- connHat^m_lj DHat_k gammaBar_im
	*/
	_3sym3 DHat2_gammaBar_llll[3];
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i], xNames[j]
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>	DHat2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
		+ partial_DHat_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>	//diverging, such that in spherical vacuum the trace of these is diag(0, 1, .5)
<?			for m,xm in ipairs(xNames) do
?>		- connHat_ull.<?=xm?>.<?=sym(l,k)?> * DHat_gammaBar_lll.<?=xm?>.<?=sym(i,j)?>
		- connHat_ull.<?=xm?>.<?=sym(l,i)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(m,j)?>
		- connHat_ull.<?=xm?>.<?=sym(l,j)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(i,m)?>
<?			end
?>	;
<?		end
	end
end
?>

	//trBar_DHat2_gammaBar_ll.ij := gammaBar^kl DHat_k DHat_l gammaBar_ij
	sym3 trBar_DHat2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	trBar_DHat2_gammaBar_ll.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * DHat2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>
<?		end
	end
?>	;
<? end
?>

	//derivative is the last index, unlike the partial_*'s
	//DHat_LambdaBar_ul.i.j := DHat_j LambdaBar^i = LambdaBar^i_,j + connHat^i_jk LambdaBar^k
	real3x3 DHat_LambdaBar_ul;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	DHat_LambdaBar_ul.<?=xi?>.<?=xj?> = partial_LambdaBar_ul[<?=j-1?>].<?=xi?> 
<?		for k,xk in ipairs(xNames) do
?>		+ connHat_ull.<?=xi?>.<?=sym(j,k)?> * LambdaBar_u.<?=xk?>
<?		end
?>	;
<?
	end
end
?>

	/*
	2017 Ruchlin eqn 12
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ijk
		+ 1/2 Delta^k Delta_jik
		+ gammaBar^kl (
			Delta^m_ki Delta_jml
			+ Delta^m_kj Delta_iml
			+ Delta^m_ik Delta_mjl
		)
	
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ikj
		+ 1/2 Delta^k Delta_jki
		+ Delta^m_ki Delta_jm^k
		+ Delta^m_kj Delta_im^k
		+ Delta^m_ik Delta_mj^k
	*/
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	RBar_ll.<?=xij?> = 0.			
			- .5 * trBar_DHat2_gammaBar_ll.<?=xij?>	
<?	for k,xk in ipairs(xNames) do
?>
			+ .5 * gammaBar_ll.<?=sym(i,k)?> * DHat_LambdaBar_ul.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_ll.<?=sym(j,k)?> * DHat_LambdaBar_ul.<?=xk?>.<?=xi?> 
			
			+ .5 * Delta_u.<?=xk?> * (Delta_lll.<?=xi?>.<?=sym(k,j)?> + Delta_lll.<?=xj?>.<?=sym(k,i)?>)

<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>			+ gammaBar_uu.<?=sym(k,l)?> * (0.
				+ Delta_ull.<?=xm?>.<?=sym(k,i)?> * Delta_lll.<?=xj?>.<?=sym(m,l)?>
				+ Delta_ull.<?=xm?>.<?=sym(k,j)?> * Delta_lll.<?=xi?>.<?=sym(m,l)?>
				+ Delta_ull.<?=xm?>.<?=sym(i,k)?> * Delta_lll.<?=xm?>.<?=sym(j,l)?>
			)
<?			end
		end
	end
?>;
<? end ?>
]], {
		eqn = self,
		solver = self.solver,
	})
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	sym3 gammaHat_ll = coord_g_ll(x);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = gammaHat_ll;
	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

<? 
-- BIG TODO (think about this one)
-- Provide initial data in algebraic form, and convert to code here.
-- This is what I was doing in HydroGPU, however it was proving slow to invert some matrices
-- (though fwiw I was only using Gauss-Jordan elimination then, now I'm using a simplification for 3x3's)
-- so maybe yes maybe no
-- until then, I'll hard-code the Minkowski option to just keep everthing zero
-- better TODO instead, maybe this is an excuse to move the rescaling code into the initState
-- then just let Minkowski call setFlatSpace() 
-- the downside to this is that whatever code produced by the initState would be specific to BSSN
-- and therefore be incompatible with the adm3d, etc solvers (unless they too adopted identical coordinate rescaling) ?
if eqn.initState.name == 'Minkowski' then ?>
	
	setFlatSpace(solver, U, x);

<? else -- not Minkowski ?>

	<?=code?>

	U->alpha = alpha;
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);

	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);
	
	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar = calc_det_gammaBar(x); 

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = cbrt(det_gammaBar / det_gamma);

	//W = exp(-2 phi)
	U->W = sqrt(exp_neg4phi);

	sym3 gammaBar_ll = sym3_real_mul(gamma_ll, exp_neg4phi);
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);

	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_real_mul(gamma_ll, 1./3. * U->K));
	sym3 ABar_ll = sym3_real_mul(A_ll, exp_neg4phi);
	U->ABar_LL = sym3_rescaleFromCoord_ll(ABar_ll, x);

	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	
	U->H = 0.;
	U->M_u = real3_zero;
	
<? end -- Minkowski ?>
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize LambdaBar_u
kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
#if 0	
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		setFlatSpace(solver, U, x);
		return;
	}
#endif

<?=eqn:makePartial'epsilon_LL'?>

	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>
	
	_3sym3 connHat_lll = coord_conn_lll(x);
	_3sym3 connHat_ull = coord_conn_ull(x);
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);

<?=eqn:getCode_connBar_ull()?>
	
	//Delta^i_jk = connBar^i_jk - connHat^i_jk
#if 0
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);
#else	
	_3sym3 Delta_ull;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
?>	Delta_ull.<?=xi?>.<?=xjk?> = connBar_ull.<?=xi?>.<?=xjk?> - connHat_ull.<?=xi?>.<?=xjk?>;
<?	end
end
?>
#endif

	real3 LambdaBar_u;
#if 0
	LambdaBar_u = _3sym3_sym3_dot23(Delta_ull, gammaBar_uu);
#else
<? for i,xi in ipairs(xNames) do
?>	LambdaBar_u.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
		for k,xk in ipairs(xNames) do
?>		+ Delta_ull.<?=xi?>.<?=sym(j,k)?> * gammaBar_uu.<?=sym(j,k)?>		
<?		end
	end
?>	;
<? end
?>
#endif

	U->LambdaBar_U = real3_rescaleFromCoord_u(LambdaBar_u, x);
}
]], table(self:getTemplateEnv(), {
		code = self.initState:initState(self.solver),
	}))
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd.cl'

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

BSSNOKFiniteDifferenceEquation.predefinedDisplayVars = {
--[=[
	'U alpha',
--	'U beta_U mag',
	'U beta_U x',
	'U beta_U y',
	'U beta_U z',
--	'U B_U mag',
	'U B_U x',
	'U B_U y',
	'U B_U z',
	--'U epsilon_LL norm',
	'U epsilon_LL xx',
	'U epsilon_LL xy',
	'U epsilon_LL xz',
	'U epsilon_LL yy',
	'U epsilon_LL yz',
	'U epsilon_LL zz',
	'U W',
	'U K',
--	'U ABar_LL tr weighted',
--] =]	
	'U ABar_LL xx',
	'U ABar_LL xy',
	'U ABar_LL xz',
	'U ABar_LL yy',
	'U ABar_LL yz',
	'U ABar_LL zz',
--[ =[	
--	'U LambdaBar_U mag',
	'U LambdaBar_U x',
	'U LambdaBar_U y',
	'U LambdaBar_U z',
	'U H',
--	'U M_u mag',
	'U M_u x',
	'U M_u y',
	'U M_u z',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_LL tr weighted',
--]=]	

-- [[ debugging derivatives
--[=[	
	'deriv alpha',
	'deriv beta_U x',
	'deriv beta_U y',
	'deriv beta_U z',
	'deriv B_U x',
	'deriv B_U y',
	'deriv B_U z',
	'deriv epsilon_LL xx',
	'deriv epsilon_LL xy',
	'deriv epsilon_LL xz',
	'deriv epsilon_LL yy',
	'deriv epsilon_LL yz',
	'deriv epsilon_LL zz',
	'deriv W',
	'deriv K',
--] =]	
	'deriv ABar_LL xx',
	'deriv ABar_LL xy',
	'deriv ABar_LL xz',
	'deriv ABar_LL yy',
	'deriv ABar_LL yz',
	'deriv ABar_LL zz',
--[ =[	
	'deriv LambdaBar_U x',
	'deriv LambdaBar_U y',
	'deriv LambdaBar_U z',
	'deriv H',
	'deriv M_u x',
	'deriv M_u y',
	'deriv M_u z',
--]=]	
--]]

	--'U tr_DBar2_phi',
	--'U DBar_phi_sq',
	--'U ABarSq tr weighted',

--[[
diverging terms in ABar_ij,t:
+ tracelessPart_ll.<?=xij?>:
	+ 2. * partial_phi_l.<?=xi?> * partial_alpha_l[<?=j-1?>]
	+ 2. * partial_phi_l.<?=xj?> * partial_alpha_l[<?=i-1?>]
	- 2. * U->alpha * DBar2_phi_ll.<?=xij?>
	+ 4. * U->alpha * partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?>
	- 8. * U->alpha * M_PI * U->S_ll.<?=xij?>

- TF_DBar2_alpha_ll.<?=xij?>:
	DBar2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]
		- partial_alpha_l[<?=k-1?>] * connBar_ull.<?=xk?>.<?=xij?>

+ U->alpha * TF_RBar_ll.<?=xij?>
	... has lots of terms

so what to check?
	RBar_ij terms ...
	especially partial_DHat_gammaBar_llll ... 

*) RBar_ij should be 0, I'm seeing diag(0, 1, cos(θ))
especially its terms ...
*) gammaBar^kl gammaBar_ij,kl should be diag(0,2,2*cos(θ)^2), but I'm seeing 0
*) gammaBar^kl gammaBar_im connHat^m_jk,l should be diag(0,-1,-1) but I'm seeing diag(0,-1,-1.5)
*) gammaBar^kl gammaBar_im,l connHat^m_jk should be diag(0,2,2) ... correct, that works. 

why is del gammaBar_ij wrong?
gamma_ij = diag(1, r^2, r^2 sin(θ)^2), which is right.
gammaHat_ij = diag(1, r^2, r^2 sin(θ)^2), which is right.
gammaBar_ij = diag(1, r^2, r^2 sin(θ)^2), which is right.
...but those partial second derivatives numerically evaluated are going to zero at some point.

replacing it with symmath simplifications gives us ...
*) gammaBar^kl gammaBar_ij,kl gives us diag(0, 2, 2*cos(θ)^2) ... which is CORRECT
--]]
	'U RBar_ll xx',
	'U RBar_ll xy',
	'U RBar_ll xz',
	'U RBar_ll yy',
	'U RBar_ll yz',
	'U RBar_ll zz',
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:append{
		{
			name = 'gamma_ll',
			type = 'sym3',
			code = [[
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	*value_sym3 = gamma_ll;
]], 
		},
		{name='gammaHat_ll', code=[[	*value_sym3 = coord_g_ll(x);]], type='sym3'},
		{name='gammaBar_ll', code=[[	*value_sym3 = calc_gammaBar_ll(U, x);]], type='sym3'},
		{name='gamma_uu', code=[[	*value_sym3 = calc_gamma_uu(U, x);]], type='sym3'},
		{name='gammaHat_uu', code=[[	*value_sym3 = coord_g_uu(x);]], type='sym3'},
		{name='gammaBar_uu', code=[[	*value_sym3 = calc_gammaBar_uu(U, x);]], type='sym3'},
		{name='K_ll', code=[[
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	*value_sym3 = sym3_real_mul(
		sym3_add(
			sym3_rescaleToCoord_LL(U->ABar_LL, x),
			sym3_real_mul(gammaBar_ll, U->K / 3.)
		), exp_4phi);
]], type='sym3'},

		{name='det gammaBar - det gammaHat', code=[[
	*value = sym3_det(calc_gammaBar_ll(U, x)) - calc_det_gammaBar(x);
]]},
		{name='det gamma based on phi', code=[[
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = det_gamma;
]]},
		{name='S', code='*value = sym3_dot(U->S_ll, calc_gamma_uu(U, x));'},
		{
			name='volume', 
			code=[[
	//|g| = exp(12 phi) |g_grid|
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = U->alpha * det_gamma;
]],
		},
		{name='f', code='*value = calc_f(U->alpha);'},
		{name='df/dalpha', code='*value = calc_dalpha_f(U->alpha);'},
	
--[=[	
		{
			name = 'ABarSq',
			type = 'sym3',
			code = [[
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	sym3 ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, ABar_UL);
	*value_sym3 = ABarSq_LL;
]],
		},
		
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = template([[
	
	<?=eqn:makePartial'epsilon_LL'?>

	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>
	
	
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);

<?=eqn:getCode_connBar_ull()?>

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}

	real tr_DBar2_phi = 0.
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		local ij = from3x3to6(i,j)
		local xij = symNames[ij]
?>		+ gammaBar_uu.<?=sym(i,j)?> * (
			partial2_phi_ll.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?		end
?>		)
<?	end
end
?>	;

	*value = tr_DBar2_phi;
]], self:getTemplateEnv())
		},
	
		{
			name = 'partial_phi_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'W'?>
	real3 partial_phi_l;
<? for i,xi in ipairs(xNames) do
?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>
	*value_real3 = partial_phi_l;
]], self:getTemplateEnv()),
		},

		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'alpha'?>
	*value_real3 = *(real3*)partial_alpha_l;
]], self:getTemplateEnv()),
		},

		{
			name = 'DBar2_phi_ll',
			code = template([[
<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}



<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>


	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
<?=eqn:getCode_connBar_ull()?>
	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?	end
?>	;
<? end
?>
	*value = sym3_dot(gammaBar_uu, DBar2_phi_ll);
]], self:getTemplateEnv()),
		},

		{
			name = 'DBar2_alpha_ll',
			code = template([[
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>


<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>


	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
<?=eqn:getCode_connBar_ull()?>

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]
<?	for k,xk in ipairs(xNames) do
?>		- partial_alpha_l[<?=k-1?>] * connBar_ull.<?=xk?>.<?=xij?>
<?	end
?>	;
<? end
?>

	*value = sym3_dot(gammaBar_uu, DBar2_alpha_ll);
]], self:getTemplateEnv()),
		},
	
	
		{
			name = 'tracelessPart_ll',
			type = 'sym3',
			code = template([[
	
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>


<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

<?=eqn:getCode_connBar_ull()?>

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]
<?	for k,xk in ipairs(xNames) do
?>		- partial_alpha_l[<?=k-1?>] * connBar_ull.<?=xk?>.<?=xij?>
<?	end
?>	;
<? end
?>


<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}

	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?	end
?>	;
<? end
?>

	
	sym3 tracelessPart_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	tracelessPart_ll.<?=xij?> = 0.
			+ 2. * partial_phi_l.<?=xi?> * partial_alpha_l[<?=j-1?>]
			+ 2. * partial_phi_l.<?=xj?> * partial_alpha_l[<?=i-1?>]
			+ U->alpha * (0.
				- 2. * DBar2_phi_ll.<?=xij?>
				+ 4. * partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?>
				- 8. * M_PI * U->S_ll.<?=xij?>
			)
		;
<? end
?>
	tracelessPart_ll = tracefree(tracelessPart_ll, gammaBar_ll, gammaBar_uu);

	*value_sym3 = tracelessPart_ll; 
]], self:getTemplateEnv()),
		},
--]=]	
	}

	local derivOrder = 2 * self.solver.numGhost
	vars:append{
--[=[
--[[ expansion:
2003 Thornburg:  ... from Wald ...
Theta = n^i_;i + K_ij n^i n^j - K
= n^i_,i + Gamma^i_ji n^j + K_ij (n^i n^j - gamma^ij)
... in ADM: n^i = -beta^i / alpha ...
= (-beta^i / alpha)_,i + Gamma^i_ji (-beta^j / alpha) + K_ij (beta^i beta^j / alpha^2 - gamma^ij)
= -beta^i_,i / alpha
	+ beta^i alpha_,i / alpha^2
	- beta^i (1/2 |g|_,i / |g|) / alpha
	+ K_ij beta^i beta^j / alpha^2
	- K

Gamma^j_ij = (ln sqrt(gamma))_,i 
= .5 (ln gamma)_,i 
= .5 gamma_,i / gamma
using gamma = gammaHat / W^6
= .5 (gammaHat W^-6)_,i / (gammaHat W^-6)
= .5 (gammaHat_,i W^-6 - 3 gammaHat W^-7) / (gammaHat W^-6)
= .5 gammaHat_,i / gammaHat - 3 / W
= GammaHat^j_ij - 3 / W
--]]
		{name='expansion', code=template([[
	<?=eqn:makePartial'W'?>
	<?=eqn:makePartial'alpha'?>
	<?=eqn:makePartial'beta_U'?>
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	real exp_4phi = 1. / calc_exp_neg4phi(U);

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma_ll = sym3_real_mul(calc_gammaBar_ll(U, x), exp_4phi);

	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_real_mul(U->ABar_ll, exp_4phi),
		sym3_real_mul(gamma_ll, U->K/3.));

	*value = -tr_partial_beta / U->alpha
<? 
for i,xi in ipairs(xNames) do
?>		+ U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		- U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		+ 3. * partial_W_l[<?=i-1?>] * U->beta_u.<?=xi?> / (U->W * U->alpha)
<?	for j,xj in ipairs(xNames) do
?>		+ K_ll.<?=sym(i,j)?> * U->beta_u.<?=xi?> * U->beta_u.<?=xj?> / (U->alpha * U->alpha)
<?	end
end
?>		- U->K;
]], 			applyCommon{
					eqn = self,
					solver = self.solver,
				}

			)
		},
--]=]		
--[=[
		--[[ ADM geodesic equation spatial terms:
		-Gamma^i_tt = 
			- gamma^ij alpha_,j

			+ alpha^-1 (
				gamma^ij beta^l gamma_kl beta^k_,j
				+ 1/2 gamma^ij gamma_kl,j beta^k beta^l
				- beta^i_,t
				- gamma^ij beta^k gamma_jk,t

				+ alpha^-1 beta^i (
					alpha_,t
					+ beta^j alpha_,j

					+ alpha^-1 (
						beta^i 1/2 beta^j beta^k gamma_jk,t
						- beta^i 1/2 beta^j beta^k beta^l gamma_kl,j
						- beta^i beta^j beta^l gamma_kl beta^k_,j
					)
				)
			)

		substitute 
		alpha_,t = -alpha^2 f K + beta^j alpha_,j
		beta^k_,t = B^k
		gamma_jk,t = -2 alpha K_jk + gamma_jk,l beta^l + gamma_lj beta^l_,k + gamma_lk beta^l_,j
		--]]
		{
			name = 'gravity',
			code= template([[
	<?=eqn:makePartial'alpha'?>

	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

	<?=eqn:makePartial'beta_u'?>

	<?=eqn:makePartial'epsilon_LL', 'partial_epsilon_lll'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
	<?=eqn:makePartial'W'?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_real_mul(partial_epsilon_lll[<?=i-1?>], _1_W * _1_W),
			sym3_real_mul(gammaBar_ll, 2. * partial_W_l[<?=i-1?>] * _1_W * _1_W * _1_W)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	sym3 dt_gamma_ll = sym3_zero;
	
	real partial_alpha_dot_beta = real3_dot(U->beta_u, *(real3*)partial_alpha_l);	//beta^j alpha_,j

	real3 beta_l = sym3_real3_mul(gamma_ll, U->beta_u);								//beta^j gamma_ij
	real3 beta_dt_gamma_l = sym3_real3_mul(dt_gamma_ll, U->beta_u);					//beta^j gamma_ij,t
	real beta_beta_dt_gamma = real3_dot(U->beta_u, beta_dt_gamma_l);				//beta^i beta^j gamma_ij,t
	
	real3 beta_dt_gamma_u = sym3_real3_mul(gamma_uu, beta_dt_gamma_l);				//gamma^ij gamma_jk,t beta^k

	//beta^i beta^j beta^k gamma_ij,k
	real beta_beta_beta_partial_gamma = 0.<?
for i,xi in ipairs(xNames) do
?> + U->beta_u.<?=xi?> * real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>)<?
end ?>;

	//beta_j beta^j_,i
	real3 beta_dbeta_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3_dot(beta_l, partial_beta_ul[<?=i-1?>]),
<? end
?>	};

	//beta_j beta^j_,i beta^i
	real beta_beta_dbeta = real3_dot(U->beta_u, beta_dbeta_l);

	//beta_j beta^j_,k gamma^ik
	real3 beta_dbeta_u = sym3_real3_mul(gamma_uu, beta_dbeta_l);

	//gamma_kl,j beta^k beta^l
	real3 beta_beta_dgamma_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>),
<? end
?>	};

	real3 beta_beta_dgamma_u = sym3_real3_mul(gamma_uu, beta_beta_dgamma_l);

<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> +=
		_1_alpha * (
			beta_dbeta_u.<?=xi?>
			+ .5 * beta_beta_dgamma_u.<?=xi?>	
			- U->B_U.<?=xi?>
			- beta_dt_gamma_u.<?=xi?>

			+ _1_alpha * U->beta_u.<?=xi?> * (
				.5 * dt_alpha
				+ partial_alpha_dot_beta

				+ _1_alpha * (
					.5 * beta_beta_dt_gamma
					- .5 * beta_beta_beta_partial_gamma 
					- beta_beta_dbeta
				)
			)
		)
	; 
<? end
?>
<? end	-- eqn.useShift ?>
]],				applyCommon{
					eqn = self,
					solver = self.solver,
				}
			), 
			type = 'real3',
		},
--]=]
-- [=[
		{
			name = 'RBar_ll',
			type = 'sym3',
			code = template([[
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial'LambdaBar_U'?>
<?=eqn:makePartial2'epsilon_LL'?>

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

	_3sym3 connHat_lll = coord_conn_lll(x);
	_3sym3 connHat_ull = coord_conn_ull(x);

	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

<?=eqn:getCode_connBar_ull()?>

	real3 LambdaBar_u = real3_rescaleToCoord_U(U->LambdaBar_U, x);
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);
	real3 Delta_u = real3_rescaleToCoord_U(Delta_U, x);

	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);
	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);

	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	sym3 gammaHat_uu = coord_g_uu(x);
	
	_3sym3 partial_connHat_ulll[3];
	coord_partial_conn_ulll(partial_connHat_ulll, x);

	sym3 RBar_ll;
	<?=eqn:getCode_RBar_ll()?>

// RBar_ij in vacuum:
//	[0,	0,	0]
//	[0,	1,	0]
//	[0,	0,	.5]
	*value_sym3 = RBar_ll;

// ... but (RBar^TF)_ij has a value, (RBar^TF)_xx
//	sym3 TF_RBar_ll = tracefree(RBar_ll, gammaBar_ll, gammaBar_uu);
//	*value_sym3 = TF_RBar_ll;

// so what is its trace?
// very high, since gammaBar_uu diverges near r=0
//	*value = sym3_dot(gammaBar_uu, RBar_ll);

]], applyCommon{
					eqn = self,
					solver = self.solver,
				}
			),
		},
--]=]		
--[=[		
		{
			name = 'RPhi_ll',
			type = 'sym3',
			code = template([[

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>

	real3 partial_phi_l;
<? for i,xi in ipairs(xNames) do
?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

	sym3 partial2_phi_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_phi_ll.<?=xij?> = .5 * (
			-partial2_W_ll[<?=ij-1?>] 
			+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
		) / U->W;
<? end ?>
	
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

<?=eqn:getCode_connBar_ull()?>

	sym3 DBar2_phi_ll = sym3_sub(
		partial2_phi_ll,
		real3_3sym3_dot1(
			partial_phi_l,
			connBar_ull
		)
	);

	real tr_DBar2_phi = sym3_dot(gammaBar_uu, DBar2_phi_ll);

	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	sym3 RPhi_ll = {
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>		.<?=xij?> = 
			- 2. * DBar2_phi_ll.<?=xij?>
			+ 4. * partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?>
			+ gammaBar_ll.<?=xij?> * (
				- 2. * tr_DBar2_phi
				- 4. * real3_weightedLenSq(
					partial_phi_l,
					gammaBar_uu
				)
			),
<? end
?>	};

	*value_sym3 = RPhi_ll;
]],				applyCommon{
					eqn = self,
					solver = self.solver,
				}
			),
		},
--]=]	

		{	-- gammaBar^kl gammaBar_ij,kl
			name = 'del gammaBar',
			type = 'sym3',
			code = template([[
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

<?=eqn:makePartial2'epsilon_LL'?>
	
	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	
	//partial2_gammaBar_llll.kl.ij := gammaBar_ij,kl
	sym3sym3 partial2_gammaBar_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
?>	partial2_gammaBar_llll.<?=xkl?>.<?=xij?> = partial2_epsilon_llll[<?=kl-1?>].<?=xij?> + partial2_gammaHat_llll.<?=xkl?>.<?=xij?>;
<?	end
end 
?>

	sym3 trBar_partial2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
?>	trBar_partial2_gammaBar_ll.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>
<?		end
	end
?>	;
<? end
?>
	*value_sym3 = trBar_partial2_gammaBar_ll;
]], self:getTemplateEnv()),
		},

--[[
gammaBar_ij = gammaHat_ij + epsilon_ij
= gammaHat_ij + epsilon_IJ e_i^I e_j^J

gammaBar_ij,kl = gammaHat_ij,kl + (epsilon_IJ e_i^I e_j^J)_,kl
= gammaHat_ij,kl + (epsilon_IJ,k e_ij^IJ + epsilon_IJ e_ij^IJ_,k)_,l
= gammaHat_ij,kl
	+ epsilon_IJ,kl e_ij^IJ 
	+ epsilon_IJ,l e_ij^IJ_,k
	+ epsilon_IJ,k e_ij^IJ_,l
	+ epsilon_IJ e_ij^IJ_,kl

gammaBar^kl gammaBar_ij,kl

gammaBar^kl = inv(gammaBar_kl)
= inv(gammaHat_kl + epsilon_kl)
--]]
		
		-- symbolically generate ... so far gammaBar^ij	
		-- this fixes it.
		{
			name='del gammaBar sym',
			type = 'sym3',
			code = template([[
<?
do
	local table = require 'ext.table'
	local symmath = require 'symmath'
	local var = symmath.var
	local sin = symmath.sin
	local cos = symmath.cos
	local Tensor = symmath.Tensor
	local Matrix = symmath.Matrix

	local compileVars = table()
	local compileRepls = table()	-- pairs of {from, to}
	local function compile(expr)
		-- TODO What about repeated expressions that could be deferred, like sin(x)?
		-- You could move them to another expression, but that means multiple code expressions.
		-- And outputting multiple lines without a function header would require yet another kind of symmath code output...
	
		local function isInteger(x) return x == math.floor(x) end
		
		-- don't do expand()
		--expr:expand()
		-- instead just expand muls
		expr = expr:map(function(x)
			if symmath.op.pow.is(x)
			and symmath.Constant.is(x[2])
			and isInteger(x[2].value)
			and x[2].value >= 1
			and x[2].value < 100
			then
				return setmetatable(table.rep({x[1]}, x[2].value), symmath.op.mul)
			end
		end)
		for _,repl in ipairs(compileRepls) do
			expr = expr:replace(repl[1], repl[2])
		end
		return expr:compile(compileVars, 'C', {hideHeader=true})
	end
	
	local coords = Tensor.coords()[1].variables
	for i,coord in ipairs(coords) do
		compileVars:insert{['x.'..xNames[i] ] = coord}
	end
	local function compileVar(...)
		local v = var(...)
		compileVars:insert(v)
		return v
	end


	local cos_xs = Tensor('_i', function(i) 
		local v = compileVar('cos_'..i) 
		compileRepls:insert{cos(coords[i]), v}
?>	real <?=v?> = cos(x.<?=xNames[i]?>);
<?		return v
	end)
	local sin_xs = Tensor('_i', function(i) 
		local v = compileVar('sin_'..i) 
		compileRepls:insert{sin(coords[i]), v}
?>	real <?=v?> = sin(x.<?=xNames[i]?>);
<?		return v
	end)


	local gammaHat = Tensor.metric().metric
	
	local det_gammaHat = Matrix.determinant(gammaHat)
?>	real det_gammaHat = <?=compile(det_gammaHat)?>;
<?	local det_gammaHatVar = compileVar('det_gammaHat', coords)
	
	local e = Tensor('_i^I', function(i,j)
		return (i==j and symmath.sqrt(gammaHat[i][i])() or 0)
	end)
	local eu = Tensor('_I^i', function(i,j)
		return i==j and (1/symmath.sqrt(gammaHat[i][i]))() or 0
	end)

	local epsilonVars = Tensor('_IJ', function(I,J)
		return compileVar('U->epsilon_LL.'..sym(I,J), coords)
	end)

	local epsilon = (epsilonVars'_IJ' * e'_i^I' * e'_j^J')()
	local gammaBar = (gammaHat'_ij' + epsilon'_ij')()

	local det_gammaBar = Matrix.determinant(gammaBar)
	det_gammaBar = (det_gammaBar / det_gammaHat * det_gammaHatVar)()
?>	real det_gammaBar = <?=compile(det_gammaBar)?>;
<?	local det_gammaBarVar = compileVar('det_gammaBar', coords)
	
	local gammaBarInv = Tensor('^ij', table.unpack((Matrix.inverse(gammaBar, nil, nil, nil, det_gammaBarVar))))
	-- this isn't always completely effective.  sometimes it introduces sin(theta)'s in the denominator
	--gammaBarInv = (gammaBarInv / det_gammaHat * det_gammaHatVar)()

	local partial_epsilonVars = Tensor('_IJk', function(I,J,k)
		local v = compileVar('partial_epsilon_LLl['..(k-1)..'].'..sym(I,J), coords)
		compileRepls:insert{epsilonVars[I][J]:diff(coords[k]), v}
		return v
	end)

	local partial2_epsilonVars = Tensor('_IJkl', function(I,J,k,l)
		local kl = from3x3to6(k,l)
		local v = compileVar('partial2_epsilon_LLll['..(kl-1)..'].'..sym(I,J), coords)
		compileRepls:insert{epsilonVars[I][J]:diff(coords[k], coords[l]), v}
		return v
	end)

?><?=eqn:makePartial('epsilon_LL', nil, nil, false)?>
<?=eqn:makePartial2('epsilon_LL', nil, nil, false)?>
<?	
	local partial2_gammaBar = gammaBar'_ij,kl'()

	local del_gammaBar = (gammaBarInv'^kl' * partial2_gammaBar'_ijkl')()

?>	sym3 del_gammaBar_ll;
<?	for i=1,3 do
		for j=i,3 do
?>	del_gammaBar_ll.<?=sym(i,j)?> = <?=compile(del_gammaBar[i][j])?>;
<?		end
	end
end
?>
	*value_sym3 = del_gammaBar_ll;
]], self:getTemplateEnv()),
		},
	
		{
			name = 'tr34 (gamma*dGamma)',
			type = 'real3x3',
			code = template([[
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
		
	_3sym3 partial_connHat_ulll[3];
	coord_partial_conn_ulll(partial_connHat_ulll, x);

	real3x3 tr34_gamma_dGamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr34_gamma_dGamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * gammaBar_ll.<?=sym(i,m)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(j,k)?>
<?				end
			end
		end
?>	;
<?	end
end
?>
	
	*value_real3x3 = tr34_gamma_dGamma_ll;
]], self:getTemplateEnv()),
		},

		{
			name = 'tr14 (Gamma*dgamma)',
			type = 'real3x3',
			code = template([[

<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

	_3sym3 connHat_ull = coord_conn_ull(x);
		
	real3x3 tr14_Gamma_dgamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr14_Gamma_dgamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,i)?> * connHat_ull.<?=xm?>.<?=sym(k,j)?>
<?				end
			end
		end
?>	;
<?	end
end
?>
	*value_real3x3 = tr14_Gamma_dgamma_ll;
]], self:getTemplateEnv()),
		},

		--[[ hmm? not working.
		{name='x', type='real3', code='*value_real3=x;'},
		--]]
		-- [[
		{name='x', code='*value=x.x;'},
		{name='y', code='*value=x.y;'},
		{name='z', code='*value=x.z;'},
		--]]
	}

	return vars
end

function BSSNOKFiniteDifferenceEquation:fillRandom(epsilon)
	local ptr = BSSNOKFiniteDifferenceEquation.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.numCells-1 do
		ptr[i].alpha = ptr[i].alpha + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return BSSNOKFiniteDifferenceEquation
