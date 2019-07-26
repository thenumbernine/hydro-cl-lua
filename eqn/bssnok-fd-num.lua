--[[
Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer" 2010
Alcubierre "Introduction to Numerical Relativity" 2008

then I'm applying 2017 Ruchlin changes...
*) separate gammaBar_ll - gammaHat_ll = epsilon_ll
*) coordinate-transform beta^i, epsilon_ij, ABar_ij, LambdaBar^i to eliminate singularities from the metric

tensors are denoted with suffixes _u _l etc for upper and lower
rescaled tensors are denoted _U _L etc

Then I'm double checking all against (and borrowing heavily from) Zach Etienne's SENR: https://math.wvu.edu/~zetienne/SENR/
--]]
local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local EinsteinEqn = require 'eqn.einstein'
local makestruct = require 'eqn.makestruct'
local common = require 'common'
local time, getTime = table.unpack(require 'util.time')
local makePartials = require 'eqn.makepartial'

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
		{name='alpha', type='real'},			-- 0:	1: alpha
		{name='W', type='real'},				-- 1:	1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},				-- 2:	1: K = K^i_i
		{name='beta_U', type='real3'},		 	-- 3:	3: beta^i
		{name='B_U', type='real3'},				-- 6:	3: B^i ... only used with HyperbolicGammaDriver
		{name='LambdaBar_U', type='real3'},		-- 9:	3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
		{name='epsilon_LL', type='sym3'},		-- 12:	6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='ABar_LL', type='sym3'},			-- 18:	6: ABar_ij, only 5 dof since ABar^k_k = 0
	}											-- 24 = total size so far

	self.consVars = table()
	:append(intVars)
	:append{
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

		{name='constrain_tr_ABar', value=true, compileTime=true},
		--{name='constrain_tr_ABar', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		
		-- 2013 Baumgarte et al, section IIIB
		--{name='diffuseCoeff', value=.001/16},
		{name='diffuseCoeff', value=0},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},

		{name='shift_eta', value=1},	--1, or 1 / (2 M), for total mass M
	}
end

function BSSNOKFiniteDifferenceEquation:fieldTypeForVar(varname)
	local _, var = self.consVars:find(nil, function(v) return v.name == varname end)
	return assert(var).type
end

function BSSNOKFiniteDifferenceEquation:makePartial1(field, fieldType, nameOverride)
	local derivOrder = 2 * self.solver.numGhost
	fieldType = fieldType or self:fieldTypeForVar(field)
	return makePartials.makePartial1(derivOrder, self.solver, field, fieldType, nameOverride)
end

function BSSNOKFiniteDifferenceEquation:makePartial2(field, fieldType, nameOverride)
	local derivOrder = 2 * self.solver.numGhost
	fieldType = fieldType or self:fieldTypeForVar(field)
	return makePartials.makePartial2(derivOrder, self.solver, field, fieldType, nameOverride)
end

-- TODO higher orders plz
-- computes the partial between U[stepSize] and U[0]
function BSSNOKFiniteDifferenceEquation:makePartialUpwind(field, fieldType, nameOverride)
	local xNames = common.xNames
	local clnumber = require 'cl.obj.number'
	fieldType = fieldType or self:fieldTypeForVar(field)

	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end

	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	local name = nameOverride or ('partial_'..field..suffix..'_upwind')
	local d1coeffs = {1}	-- TODO more plz 
	
	local lines = table()
	if fieldType == 'real' then
		lines:insert('\treal3 '..name..';')
	else
		lines:insert('\t'..fieldType..' '..name..'[3];')
	end
	for i,xi in ipairs(xNames) do
		local namei
		if fieldType == 'real' then
			namei = name..'.'..xi
		else
			namei = name..'['..(i-1)..']'
		end
		local expr = zero
		if i <= self.solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				local UR = 'Uup['..j..' * solver->stepsize.'..xi..'].'..field
				local UL = 'Uup['..(j-1)..' * solver->stepsize.'..xi..'].'..field
				expr = add(expr, real_mul(sub(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

function BSSNOKFiniteDifferenceEquation:compile(expr)
	return self.solver.coord:compile(expr)
end

function BSSNOKFiniteDifferenceEquation:getEnv()
	return table(common, {
		eqn = self,
		solver = self.solver,
	})
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template(file['eqn/bssnok-fd-num.cl'], table(self:getEnv(), {getCommonCode=true}))
end

--[[
Should initState provide a metric in cartesian, or in the background metric?
I'll say Cartesian for now, and then transform them using the rescaling.
--]]
function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([=[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	sym3 gammaHat_ll = calc_gammaHat_ll(x);

	//initState will assume it is providing a metric in Cartesian
	real alpha = 1.;
	real3 beta_u = real3_zero;
	real3 B_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;
	
	real rho = 0.;

	<?=code?>
	
	//rescale from cartesian to spherical
	//TODO what about rotations for change of coordinates?
	// this is another reason why initial conditions should be symbolic
//	beta_u = real3_rescaleToCoord_U(beta_u, x);
	gamma_ll = sym3_rescaleToCoord_LL(gamma_ll, x);
//	K_ll = sym3_rescaleToCoord_LL(K_ll, x);

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
	
	U->B_U = real3_rescaleFromCoord_u(B_u, x);

	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0.;
	U->M_u = real3_zero;
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize LambdaBar_u
//TODO do this symbolically.  That's what I originally did, but symbolic calculations were getting complex
// however, with spherical BSSN, you need to 
kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	//Delta^i_jk = connBar^i_jk - connHat^i_jk
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	
	U->LambdaBar_U = real3_add(_3sym3_sym3_dot23(Delta_ULL, gammaBar_UU), mystery_C_U);
}
]=], table(self:getEnv(), {
		code = self.initState:initState(self.solver),
	}))
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd-num.cl'

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

BSSNOKFiniteDifferenceEquation.predefinedDisplayVars = {
-- [=[
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
-- [=[	
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
-- [=[	
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
	'U ABarSq_LL tr weighted gammaBar^IJ',

--[[ should be zero for Minkowski
	'U RBar_LL xx',
	'U RBar_LL xy',
	'U RBar_LL xz',
	'U RBar_LL yy',
	'U RBar_LL yz',
	'U RBar_LL zz',
--]]
--[[ should be zero for Minkowski
	'U gammaBar_UU xx',
	'U gammaBar_UU xy',
	'U gammaBar_UU xz',
	'U gammaBar_UU yy',
	'U gammaBar_UU yz',
	'U gammaBar_UU zz',
--]]
--[[ should be zero for Minkowski
	'U Delta_ULL x xx',
	'U Delta_ULL x xy',
	'U Delta_ULL x xz',
	'U Delta_ULL x yy',
	'U Delta_ULL x yz',
	'U Delta_ULL x zz',
	'U Delta_ULL y xx',
	'U Delta_ULL y xy',
	'U Delta_ULL y xz',
	'U Delta_ULL y yy',
	'U Delta_ULL y yz',
	'U Delta_ULL y zz',
	'U Delta_ULL z xx',
	'U Delta_ULL z xy',
	'U Delta_ULL z xz',
	'U Delta_ULL z yy',
	'U Delta_ULL z yz',
	'U Delta_ULL z zz',
--]]
--[[ should be zero for Minkowski
	'U Delta_U x',
	'U Delta_U y',
	'U Delta_U z',
--]]
--[[ should be diag(0, 2, 2 cos(theta)^2)
	'U trBar_partial2_gammaBar_ll xx',
	'U trBar_partial2_gammaBar_ll xy',
	'U trBar_partial2_gammaBar_ll xz',
	'U trBar_partial2_gammaBar_ll yy',
	'U trBar_partial2_gammaBar_ll yz',
	'U trBar_partial2_gammaBar_ll zz',
--]]
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)
	local env = self:getEnv()

	vars:append{
		{name='gamma_ll', code = [[	*value_sym3 = calc_gamma_ll(U, x);]], type='sym3'},
		{name='gamma_uu', code=[[	*value_sym3 = calc_gamma_uu(U, x);]], type='sym3'},
		{name='gammaHat_ll', code=[[	*value_sym3 = calc_gammaHat_ll(x);]], type='sym3'},
		{name='gammaHat_uu', code=[[	*value_sym3 = calc_gammaHat_uu(x);]], type='sym3'},
		{name='gammaBar_ll', code=[[	*value_sym3 = calc_gammaBar_ll(U, x);]], type='sym3'},
		{name='gammaBar_uu', code=[[	*value_sym3 = calc_gammaBar_uu(U, x);]], type='sym3'},
		{name='gammaBar_LL', code=[[	*value_sym3 = calc_gammaBar_LL(U, x);]], type='sym3'},
		{name='gammaBar_UU', code=[[	*value_sym3 = calc_gammaBar_UU(U, x);]], type='sym3'},
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
	
		{
			name = 'ABarSq_LL',
			type = 'sym3',
			code = [[
	sym3 gammaBar_UU = calc_gammaBar_UU(U, x);
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	sym3 ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, ABar_UL);
	*value_sym3 = ABarSq_LL;
]],
		},
		
--[=[	
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = template([[
	
<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>
	
	//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);
	
	sym3 partial2_phi_ll;
	{
		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll.[<?=xij?> 
				+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
			) / U->W;
<? end ?>
	}

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	sym3 partial2_phi_LL = sym3_rescaleFromCoord_ll(partial2_phi_ll, x);

	//tr_DBar2_phi := gammaBar^ij DBar_i DBar_j phi = gammaBar^ij phi_,ij - connBar^k phi_,k
	real tr_DBar2_phi = 0.
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		local ij,xij = from3x3to6(i,j)
?>		+ gammaBar_UU.<?=sym(i,j)?> * (
			partial2_phi_LL.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?		end
?>		)
<?	end
end
?>	;

	*value = tr_DBar2_phi;
]], env)
		},


		{
			name = 'DBar2_phi_ll',
			code = template([[
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>
	
	//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);
	
	sym3 partial2_phi_ll;
	{
		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll.<?=xij?> 
				+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
			) / U->W;
<? end ?>
	}

<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	//DBar2_phi := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_LL.<?=xij?> = partial2_phi_LL.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?		end
?>	;
<?	end
end
?>	

	*value = sym3_dot(gammaBar_UU, DBar2_phi_LL);
]], env),
 		},
--]=]

		{
			name = 'partial_phi_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial1'W'?>
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);
	*value_real3 = partial_phi_l;
]], env),
		},

		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial1'alpha'?>
	*value_real3 = partial_alpha_l;
]], env),
		},

--[=[
		{
			name = 'DBar2_alpha_ll',
			code = template([[
<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=eqn:makePartial1'epsilon_LL'?>
	
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_LL.<?=xij?> = partial2_alpha_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_alpha_L.<?=xk?>
<?	end
?>	;
<? end
?>

	*value = sym3_dot(gammaBar_UU, DBar2_alpha_LL);
]], env),
		},
	
		{
			name = 'tracelessPart_LL',
			type = 'sym3',
			code = template([[
	
<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=eqn:makePartial1'epsilon_LL'?>
	
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBar(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);

	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	sym3 partial2_alpha_LL = sym3_rescaleFromCoord_ll(partial2_alpha_ll, x);
	
	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_LL.<?=xij?> = partial2_alpha_LL.<?=xij?>
<?	for k,xk in ipairs(xNames) do
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_alpha_L.<?=xk?>
<?	end
?>	;
<? end
?>

<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);
	
	sym3 partial2_phi_ll;
	{
		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				- partial2_W_ll.<?=xij?> 
				+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
			) / U->W;
<? end ?>
	}

	sym3 partial2_phi_LL = sym3_rescaleFromCoord_ll(partial2_phi_ll, x);

	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_LL;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_LL.<?=xij?> = partial2_phi_LL.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ULL.<?=xk?>.<?=xij?> * partial_phi_L.<?=xk?>
<?	end
?>	;
<? end
?>

	real3 partial_alpha_L = real3_rescaleFromCoord_l(partial_alpha_l, x);

	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	
	sym3 tracelessPart_LL;
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>	tracelessPart_LL.<?=xij?> = (0.
		+ 2. * partial_phi_L.<?=xi?> * partial_alpha_L.<?=xj?>
			+ 2. * partial_phi_L.<?=xj?> * partial_alpha_L.<?=xi?>
			+ U->alpha * (0.
				- 2. * DBar2_phi_LL.<?=xij?>
				+ 4. * partial_phi_L.<?=xi?> * partial_phi_L.<?=xj?>
				- 8. * M_PI * U->S_ll.<?=xij?> / (calc_len_<?=xi?>(x) * calc_len_<?=xj?>(x))
			)
		);
<? end
?>
	tracelessPart_LL = tracefree(tracelessPart_LL, gammaBar_LL, gammaBar_UU);

	*value_sym3 = tracelessPart_LL; 
]], env),
		},
--]=]	
	}

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
	vars:insert{
		name = 'expansion', 
		code = template([[
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial1'beta_U'?>
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

	//tr_partial_beta := beta^i_,i
	real tr_partial_beta = 0. <?
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	*value = -tr_partial_beta / U->alpha
<? 
for i,xi in ipairs(xNames) do
?>		+ U->beta_u.<?=xi?> * partial_alpha_l.<?=xi?> / (U->alpha * U->alpha) 
		- U->beta_u.<?=xi?> * partial_alpha_l.<?=xi?> / (U->alpha * U->alpha) 
		+ 3. * partial_W_l.<?=xi?> * U->beta_u.<?=xi?> / (U->W * U->alpha)
<?	for j,xj in ipairs(xNames) do
?>		+ K_ll.<?=sym(i,j)?> * U->beta_u.<?=xi?> * U->beta_u.<?=xj?> / (U->alpha * U->alpha)
<?	end
end
?>		- U->K;
]], env)
	}
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
	vars:insert{
		name = 'gravity',
		code = template([[
<?=eqn:makePartial1'alpha'?>

	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

<?=eqn:makePartial1'beta_U'?>
<?=eqn:makePartial1'epsilon_LL'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
<?=eqn:makePartial1'W'?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_real_mul(partial_epsilon_lll[<?=i-1?>], _1_W * _1_W),
			sym3_real_mul(gammaBar_ll, 2. * partial_W_l.<?=xi?> * _1_W * _1_W * _1_W)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	sym3 dt_gamma_ll = sym3_zero;
	
	real partial_alpha_dot_beta = real3_dot(U->beta_u, partial_alpha_l);	//beta^j alpha_,j

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
	real3 beta_dbeta_l;
<? for i,xi in ipairs(xNames) do
?>	beta_dbeta_l.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
?>		+ beta_l.<?=xj?> * partial_beta_ul[<?=i-1?>].<?=xj?>
<?	end
?>	;
<? end
?>

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
]], env), 
		type = 'real3',
	}
--]=]
-- [=[
	do
		-- should be zero when gammaBar_ij = gammaHat_ij
		vars:insert{
			name = 'RBar_LL',
			type = 'sym3',
			code = template([[
	
	//partial_LambdaBar_Ul[j].i := LambdaBar^i_,j
<?=eqn:makePartial1'LambdaBar_U'?>
	//partial_LambdaBar_UL.I.J := e_i^I (Lambda^M e^i_M)_,j e^j_J
	real3x3 partial_LambdaBar_UL = real3x3_partial_rescaleFromCoord_Ul(U->LambdaBar_U, partial_LambdaBar_Ul, x);

<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	_3sym3 Delta_LLL = sym3_3sym3_mul(gammaBar_LL, Delta_ULL);

<?=eqn:makePartial2'epsilon_LL'?>

	sym3 trBar_partial2_gammaBar_ll = calc_trBar_partial2_gammaBar_ll(
		U, 
		x, 
		&gammaBar_UU, 
		partial_epsilon_LLl, 
		partial2_epsilon_LLll);

	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);	
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);	
	sym3 RBar_LL = calc_RBar_LL(
		U,
		x,
		&gammaBar_LL,
		&gammaBar_UU,
		&connHat_ULL,
		&partial_gammaBar_LLL,
		&trBar_partial2_gammaBar_ll,
		&partial_LambdaBar_UL,
		&Delta_U,
		&Delta_ULL,
		&Delta_LLL);

	*value_sym3 = sym3_rescaleToCoord_LL(RBar_LL, x);
]], env),
		}
	end
--]=]		
--[=[		
	vars:insert{
		name = 'RPhi_LL',
		type = 'sym3',
		code = template([[

<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>

	real3 partial_phi_l = real3_real_mul(partial_W_l, -.5 / U->W);

	sym3 partial2_phi_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_phi_ll.<?=xij?> = .5 * (
			-partial2_W_ll.<?=xij?> 
			+ partial_W_l.<?=xi?> * partial_W_l.<?=xj?> / U->W
		) / U->W;
<? end ?>
	
	real3 partial_phi_L = real3_rescaleFromCoord_l(partial_phi_l, x);
	sym3 partial2_phi_LL = sym3_rescaleFromCoord_ll(partial2_phi_ll, x);
	
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	sym3 DBar2_phi_LL = sym3_sub(
		partial2_phi_LL,
		real3_3sym3_dot1(
			partial_phi_L,
			connBar_ULL
		)
	);

	real tr_DBar2_phi = sym3_dot(gammaBar_UU, DBar2_phi_LL);

	//2008 Alcubierre eqn 2.8.18
	//2010 Baumgarte, Shapiro eqn 3.10
	sym3 RPhi_LL = {
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		.<?=xij?> = 
			- 2. * DBar2_phi_LL.<?=xij?>
			+ 4. * partial_phi_L.<?=xi?> * partial_phi_L.<?=xj?>
			+ gammaBar_LL.<?=xij?> * (
				- 2. * tr_DBar2_phi
				- 4. * real3_weightedLenSq(
					partial_phi_L,
					gammaBar_UU
				)
			),
<? end
?>	};

	*value_sym3 = RPhi_ll;
]], env),
	}

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
--]=]	

--[=[
	vars:insert{
		name='trBar_partial2_gammaBar_ll',
		type = 'sym3',
		code = template([[
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
<?=eqn:makePartial1'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>
	
	sym3 trBar_partial2_gammaBar_ll = calc_trBar_partial2_gammaBar_ll(
		U, 
		x, 
		&gammaBar_UU, 
		partial_epsilon_LLl, 
		partial2_epsilon_LLll);

	*value_sym3 = trBar_partial2_gammaBar_ll;
]], env),
	}
--]=]

--[=[
	vars:insert{
		name = 'tr34 (gamma*dGamma)',
		type = 'real3x3',
		code = template([[

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);

	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);
	
	real3x3 tr34_gamma_dGamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr34_gamma_dGamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * gammaBar_LL.<?=sym(i,m)?> * calc_len_<?=xi?>(x) * calc_partial_connHat_Ulll_<?=xm..sym(j,k)..xl?>(x)
<?				end
			end
		end
?>	;
<?	end
end
?>
	
	*value_real3x3 = tr34_gamma_dGamma_ll;
]], env),
	}
--]=]

--[=[
	vars:insert{
		name = 'tr14 (Gamma*dgamma)',
		type = 'real3x3',
		code = template([[
<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);
	
	_3sym3 connHat_ull = _3sym3_rescaleToCoord_ULL(connHat_ULL, x); 

	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll = _3sym3_rescaleToCoord_LLL(partial_gammaBar_LLL, x);

	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);
	
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
]], env),
	}
--]=]

--[=[
	for i,xi in ipairs(xNames) do
		vars:insert{
			name = 'Delta_ULL '..xi,
			type = 'sym3',
			code = template([[
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);

	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK + connHat^I_JK
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

	*value_sym3 = Delta_ULL.<?=xi?>;
]], table(env, {i=i, xi=xi})),
		}
	end
--]=]

--[=[
	vars:insert{
		name = 'Delta_U',
		type = 'real3',
		code = template([[
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);

<?=eqn:makePartial1'epsilon_LL'?>
	_3sym3 partial_gammaBar_LLL = calc_partial_gammaBar_LLL(x, U->epsilon_LL, partial_epsilon_LLl);
	
	_3sym3 connBar_ULL;
	calc_connBar_ULL(&connBar_ULL, &partial_gammaBar_LLL, &gammaBar_UU);
	
	_3sym3 connHat_LLL, connHat_ULL;
	calc_connHat_LLL_and_ULL(&connHat_LLL, &connHat_ULL, U, x);

	//Delta_ULL[I].JK := Delta^I_JK = connBar^I_JK + connHat^I_JK
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);

<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = sym3_dot(gammaBar_UU, Delta_ULL.<?=xi?>);
<? end
?>
]], env),
	}
--]=]

	vars:insert{name='x', type='real3', code='*value_real3=x;'}

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
