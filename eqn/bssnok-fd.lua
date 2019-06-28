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
		--{name='constrain_det_gammaBar_ll', value=true, compileTime=true},
		{name='constrain_det_gammaBar_ll', value=false, compileTime=true},

		--{name='constrain_tr_ABar_ll', value=true, compileTime=true},
		{name='constrain_tr_ABar_ll', value=false, compileTime=true},
		
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
		makePartial = function(...) return self:makePartial(...) end,
		makePartial2 = function(...) return self:makePartial2(...) end,
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
	return sym3_add(gammaHat_ll, epsilon_ll);
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2018 Ruchlin
real calc_det_gammaBar_ll(real3 x) {
	return coord_det_g(x);
}

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar_ll = calc_det_gammaBar_ll(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar_ll);
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	real det_gamma_ll = coord_det_g(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma_ll); 
	return gamma_uu;
}

void calc_connBar_lll(
	_3sym3* connBar_lll,
	const sym3* partial_epsilon_lll,//[3]
	const _3sym3* connHat_lll,
	real3 x
) {
	/*
	connBar_lll[i].jk := connBar_ijk 
	= 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	= 1/2 (
		epsilon_ij,k + gammaHat_ij,k 
		+ epsilon_ik,j + gammaHat_ik,j 
		- epsilon_jk,i - gammaHat_jk,i
	)
	= 1/2 (epsilon_ij,k + epsilon_ik,j - epsilon_jk,i) + connHat_ijk
	*/
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
?>	connBar_lll-><?=xi?>.<?=xjk?> = .5 * (0.
			+ partial_epsilon_lll[<?=k-1?>].<?=sym(i,j)?>
			+ partial_epsilon_lll[<?=j-1?>].<?=sym(i,k)?> 
			- partial_epsilon_lll[<?=i-1?>].<?=xjk?>
		) + connHat_lll-><?=xi?>.<?=xjk?>;
<?	end
end
?>
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
	U->LambdaBar_U = real3_rescaleFromCoord_u(mystery_C_U, x);

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

function BSSNOKFiniteDifferenceEquation:getCode_RBar_ll()
	return template([[
	/*
	trBar_DHat2_gammaBar_ll.ij =
	gammaBar^kl DHat_k DHat_l gammaBar_ij
expand DHat_l:
	= gammaBar^kl DHat_k (gammaBar_ij,l - connHat^m_il gammaBar_mj - connHat^m_jl gammaBar_im)
expand DHat_k:
	= gammaBar^kl (
		(gammaBar_ij,l - connHat^m_il gammaBar_mj - connHat^m_jl gammaBar_im)_,k
		- connHat^n_ik (gammaBar_nj,l - connHat^m_nl gammaBar_mj - connHat^m_jl gammaBar_nm)
		- connHat^n_jk (gammaBar_in,l - connHat^m_il gammaBar_mn - connHat^m_nl gammaBar_im)
		- connHat^n_lk (gammaBar_ij,n - connHat^m_in gammaBar_mj - connHat^m_jn gammaBar_im)
	)
distribute.  no raises and lowers.
	= gammaBar^kl (
		gammaBar_ij,kl
		- gammaBar_im connHat^m_jk,l
		- gammaBar_jm connHat^m_ik,l
		- 2 gammaBar_im,k connHat^m_jl
		- 2 gammaBar_jm,k connHat^m_il 
		- gammaBar_ij,m connHat^m_kl
		+ connHat^n_ik (gammaBar_mj connHat^m_nl + gammaBar_nm connHat^m_jl)
		+ connHat^n_jk (gammaBar_mn connHat^m_il + gammaBar_im connHat^m_nl)
		+ connHat^n_lk (gammaBar_mj connHat^m_in + gammaBar_im connHat^m_jn)
	)
replace gammaBar_ij,kl with epsilon_ij,kl + gammaHat_ij,kl 
substitute gammaBar_il connHat^l_jk = Q_ijk (symmetric in indexes 2 & 3)
substitute gammaBar_im connHat^m_jk,l = P_ijkl (symmetric in indexes 2 & 3)
substitute gammaBar_im,j connHat^m_kl = S_ijkl (symmetric in 3 & 4)	
	= gammaBar^kl (
		epsilon_ij,kl
		+ gammaHat_ij,kl
		- P_ijkl
		- P_jikl
		- 2 S_iljk
		- 2 S_jlik
		- gammaBar_ij,m connHat^m_kl
		+ connHat^m_ik (Q_jlm + Q_mlj)
		+ connHat^m_jk (Q_ilm + Q_mli)
		+ connHat^m_lk (Q_imj + Q_jmi)
	)
	*/
	sym3 trBar_DHat2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	trBar_DHat2_gammaBar_ll.<?=xij?> = 0.<?
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do 
			local kl = from3x3to6(k,l)
?> 			+ gammaBar_uu.<?=sym(k,l)?> * (0.
				+ partial2_epsilon_llll[<?=kl-1?>].<?=xij?>
				+ partial2_gammaHat_llll.<?=sym(k,l)?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>				- gammaBar_ll.<?=sym(i,m)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(j,k)?>
				- gammaBar_ll.<?=sym(j,m)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(i,k)?>	
				- 2. * partial_gammaBar_lll.<?=xl?>.<?=sym(i,m)?> * connHat_ull.<?=xm?>.<?=sym(j,k)?>
				- 2. * partial_gammaBar_lll.<?=xl?>.<?=sym(j,m)?> * connHat_ull.<?=xm?>.<?=sym(i,k)?>	
				- partial_gammaBar_lll.<?=xm?>.<?=xij?> * connHat_ull.<?=xm?>.<?=sym(k,l)?>
				+ connHat_ull.<?=xm?>.<?=sym(i,k)?> * (0.
<?				for n,xn in ipairs(xNames) do				
?>					+ gammaBar_ll.<?=sym(j,n)?> * connHat_ull.<?=xn?>.<?=sym(l,m)?> 
					+ gammaBar_ll.<?=sym(m,n)?> * connHat_ull.<?=xn?>.<?=sym(l,j)?>
<?				end
?>				)
				+ connHat_ull.<?=xm?>.<?=sym(j,k)?> * (0.
<?				for n,xn in ipairs(xNames) do				
?>					+ gammaBar_ll.<?=sym(i,n)?> * connHat_ull.<?=xn?>.<?=sym(l,m)?> 
					+ gammaBar_ll.<?=sym(m,n)?> * connHat_ull.<?=xn?>.<?=sym(l,i)?>
<?				end
?>				)
				+ connHat_ull.<?=xm?>.<?=sym(l,k)?> * (0.
<?				for n,xn in ipairs(xNames) do				
?>					+ gammaBar_ll.<?=sym(i,n)?> * connHat_ull.<?=xn?>.<?=sym(m,j)?> 
					+ gammaBar_ll.<?=sym(j,n)?> * connHat_ull.<?=xn?>.<?=sym(m,i)?>
<?				end
?>				)
<?			end
?>			)
<?		end
	end	
?>	;
<? end
?>

	/*
	2018 Ruchlin eqn 12
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
#warning the RBar_ij calculation is the only thing not working 
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	RBar_ll.<?=xij?> = 0.			
			- .5 * trBar_DHat2_gammaBar_ll.<?=xij?>		
<?	for k,xk in ipairs(xNames) do
?>
			+ .5 * Delta_u.<?=xk?> * (Delta_lll.<?=xi?>.<?=sym(j,k)?> + Delta_lll.<?=xj?>.<?=sym(i,k)?>)
<?		for l,xl in ipairs(xNames) do
?>

			+ .5 * (
				partial_LambdaBar_ul[<?=i-1?>].<?=xk?> 
				+ connHat_ull.<?=xk?>.<?=sym(i,l)?> * LambdaBar_u.<?=xl?>
			) * gammaBar_ll.<?=sym(k,j)?>
			+ .5 * (
				partial_LambdaBar_ul[<?=j-1?>].<?=xk?> 
				+ connHat_ull.<?=xk?>.<?=sym(j,l)?> * LambdaBar_u.<?=xl?>
			) * gammaBar_ll.<?=sym(k,i)?>

<?			for m,xm in ipairs(xNames) do
?>
			+ gammaBar_uu.<?=sym(k,l)?> * (0.
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

	<?=code?>

	U->alpha = alpha;
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);

	real det_gamma_ll = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma_ll);
	
	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar_ll = calc_det_gammaBar_ll(x); 

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = cbrt(det_gammaBar_ll / det_gamma_ll);

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

<?=makePartial'epsilon_LL'?>

	_3sym3 connHat_lll = coord_conn_lll(x);
	_3sym3 connHat_ull = coord_conn_ull(x);

	// TODO I am suspicious that I should be creating connBar_ull another way, 
	// maybe by using LambdaBar^i and W?	
	//connBar_lll[i].jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
	//Not here, of course.  This is where LambdaBar^i is initialized.
	_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>	connBar_lll.<?=xi?>.<?=xjk?> = .5 * (
//the next three sums are diverging 
// unless I add this '0. +' to it ... hmm ...
		0. +
			partial_epsilon_lll[<?=k-1?>].<?=sym(i,j)?>
			+ partial_epsilon_lll[<?=j-1?>].<?=sym(i,k)?> 
			- partial_epsilon_lll[<?=i-1?>].<?=xjk?>
		)
		+ connHat_lll.<?=xi?>.<?=xjk?>
	;
<?	end
end
?>
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);

	//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
#if 0
	_3sym3 connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);
#else
	_3sym3 connBar_ull;
<? for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
?>	connBar_ull.<?=xi?>.<?=xjk?> = 0.
<?		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(i,l)?> * connBar_lll.<?=xl?>.<?=xjk?>
<?		end
?>	;
<?	end
end
?>
#endif

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
	'U alpha',
	'U beta_U mag',
	--'U epsilon_LL norm',
	'U W',
	'U ABar_LL tr weighted',
	'U ABar_LL x x',
	'U ABar_LL x y',
	'U ABar_LL x z',
	'U ABar_LL y y',
	'U ABar_LL y z',
	'U ABar_LL z z',
	'U K',
	'U LambdaBar_U mag',
	'U H',
	'U M_U mag',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_LL tr weighted',

-- [[ debugging derivatives
	'deriv alpha',
	'deriv beta_U mag',
	'deriv W',
	'deriv ABar_LL tr weighted',
	'deriv ABar_LL x x',
	'deriv ABar_LL x y',
	'deriv ABar_LL x z',
	'deriv ABar_LL y y',
	'deriv ABar_LL y z',
	'deriv ABar_LL z z',
	'deriv K',
	'deriv LambdaBar_U mag',
--]]
	
	-- debugging
	'U ABarSq_ll tr weighted',	-- has problems with spherical vacuum spacetime
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
	*value = sym3_det(calc_gammaBar_ll(U, x)) - calc_det_gammaBar_ll(x);
]]},
		{name='det gamma based on phi', code=[[
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	*value = exp_12phi * calc_det_gammaBar_ll(x);
]]},
		{name='S', code='*value = sym3_dot(U->S_ll, calc_gamma_uu(U, x));'},
		{
			name='volume', 
			code=[[
	//|g| = exp(12 phi) |g_grid|
	real exp_neg4phi = calc_exp_neg4phi(U);
	real det_gamma_ll = coord_det_g(x) / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	*value = U->alpha * det_gamma_ll;
]],
		},
	}

	local derivOrder = 2 * self.solver.numGhost
	vars:append{

		{name='f', code='*value = calc_f(U->alpha);'},
		{name='df/dalpha', code='*value = calc_dalpha_f(U->alpha);'},

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
	<?=makePartial'W'?>
	<?=makePartial'alpha'?>
	<?=makePartial'beta_U'?>
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
					makePartial = function(...) return self:makePartial(...) end,
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
	<?=makePartial'alpha'?>

	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

	<?=makePartial'beta_u'?>

	<?=makePartial'epsilon_LL', 'partial_epsilon_lll'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
	<?=makePartial'W'?>
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
					makePartial = function(...) return self:makePartial(...) end,
				}
			), 
			type = 'real3',
		},
--]=]
--[=[	
		{
			name = 'RBar_ll',
			type = 'sym3',
			code = template([[
	_3sym3 connHat_ull = coord_conn_ull(x);
<?=makePartial'epsilon_LL'?>
	_3sym3 connBar_lll;
	calc_connBar_lll(&connBar_lll, partial_epsilon_lll, &connHat_lll, x);
	_3sym3 connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);
	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);

	sym3sym3 partial2_gammaHat_llll = coord_d2g_llll(x);
	sym3 gammaHat_uu = coord_g_uu(x);
	
	_3sym3 partial_gammaHat_lll = coord_dg_lll(x);
	_3sym3 partial_gammaBar_lll = _3sym3_add(_3sym3_rescaleToCoord(*(_3sym3*)partial_epsilon_lll, x), partial_gammaHat_lll);
	
	_3sym3 partial_connHat_ulll[3];
	calc_partial_conn_ulll(partial_connHat_ulll, &connHat_ull, &gammaHat_uu, &partial_gammaHat_lll, &partial2_gammaHat_llll);
	
	sym3sym3 partial2_gammaBar_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
?>	partial2_gammaBar_llll.<?=xij?>.<?=xkl?> = partial2_epsilon_llll[<?=ij-1?>].<?=xkl?> + partial2_gammaHat_llll.<?=xij?>.<?=xkl?>;
<?	end
end 
?>
	
	sym3 RBar_ll;
	<?=eqn:getCode_RBar_ll()?>
	*value_sym3 = RBar_ll;
]], applyCommon{
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return self:makePartial(...) end,
				}
			),
		},
		{
			name = 'RPhi_ll',
			type = 'sym3',
			code = template([[
	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	sym3 RPhi_ll = {
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = 
			- 2. * DBar2_phi_ll.<?=xij?>
			+ 4. * partial_phi_l[<?=i-1?>] * partial_phi_l[<?=j-1?>]
			+ gammaBar_ll.<?=xij?> * (
				- 2. * tr_DBar2_phi
				- 4. * real3_weightedLenSq(
					*(real3*)partial_phi_l,
					gammaBar_uu
				)
			),
<? end
?>	};

	*value_sym3 = RPhi_ll;
]],				applyCommon{
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return self:makePartial(...) end,
				}
			),
		},
--]=]	
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
