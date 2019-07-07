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
local applyCommon = require 'common'
local time, getTime = table.unpack(require 'util.time')

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

		{name='constrain_tr_ABar', value=true, compileTime=true},
		--{name='constrain_tr_ABar', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},
	}
end

function BSSNOKFiniteDifferenceEquation:makePartial(field, fieldType, nameOverride)
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	local derivOrder = 2 * self.solver.numGhost
	return makePartial(derivOrder, self.solver, field, fieldType, nameOverride)
end

function BSSNOKFiniteDifferenceEquation:makePartial2(field, fieldType, nameOverride)
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	local derivOrder = 2 * self.solver.numGhost
	return makePartial2(derivOrder, self.solver, field, fieldType, nameOverride)
end


-- NOTICE you have to use this in the scope of a push/pop of the previous fenv
-- ... and you have to set this's mt's __index to the old fenv
function BSSNOKFiniteDifferenceEquation:getCASEnv()
	if self.casEnv then return self.casEnv end
	local derivOrder = 2 * self.solver.numGhost
	local casEnv = applyCommon{
		eqn = self,
		solver = self.solver,
	}
	setmetatable(casEnv, {})
	self.casEnv = casEnv

	-- NOTICE if there is a local variable in the outer scope
	-- then it will not get set in the nested new env
	-- [[ do this every time you use the casEnv
	local oldEnv = getfenv()
	getmetatable(casEnv).__index = oldEnv
	setfenv(1, casEnv)
	--]]

	-- this is the local env code for casEnv
	symmath = require 'symmath'
symmath.tostring = require 'symmath.tostring.MathJax'
print(symmath.tostring.header)	-- TODO output to a new file?
printbr = print
	var = symmath.var
	sin = symmath.sin
	cos = symmath.cos
	sqrt = symmath.sqrt
	frac = symmath.frac
	Tensor = symmath.Tensor
	Matrix = symmath.Matrix

	-- this is a list of variables passed to compile()
	compileVars = table()
	
	--[[
	this is a list of pairs of {from, to}
	the 'to' is usually a variable
	the 'from' is usually partial derivatives, but can be things such as expensive functions you want to minimize calling (sin, cos, etc)
	it is used to replace expressions with variables before compiling
	--]]
	compileRepls = table()	
	
	function compile(expr)
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
	
	coords = Tensor.coords()[1].variables
	for i,coord in ipairs(coords) do
		compileVars:insert{['x.'..xNames[i] ] = coord}
	end

	-- create a variable to be added to the 'compileVars' list
	function compileVar(...)
		local v = var(...)
		compileVars:insert(v)
		return v
	end
	
	-- compile an expression (typically variable), insert it into compileRepls as being associated with a specific expression
	function compileReplVar(repl, ...)
		local v = compileVar(...)
		compileRepls:insert{repl, v}
		return v
	end

	-- TODO only generate what we need
	cos_xs = Tensor('_i', function(i) 
		return compileReplVar(cos(coords[i]), 'cos_'..i) 
	end)
	sin_xs = Tensor('_i', function(i) 
		return compileReplVar(sin(coords[i]), 'sin_'..i) 
	end)

	-- generates assignment code from specified variables
	function assign(name, with)
		if with == nil then with = casEnv[name] end
		assert(symmath.Expression.is(with), "not an expression")
		return '\treal '..name..' = '..compile(with)..';'
	end
	
	-- generate assignment code for a variable based on its associated expr in the compileRepls table
	function assignFromRepl(v)
		local index = compileRepls:find(nil, function(p)
			if p[2] == v then return p[1] end
		end)
		assert(index, "failed to find var "..v)
		-- hmm, compile() can't handle the repl code, because it is programmed to replace it...
		-- so this rule needs to be pushed/popped ...
		-- seems order of operations is really sensitive here
		local replpair = compileRepls:remove(index)
		local result = assign(v.name, replpair[1])
		compileRepls:insert(index, replpair)
		return result
	end
	
	-- assign any vars in the compileRepls table to their respective compileRepls exprs
	function assignAllFromRepl(vars)
		local lines = table()
		for i,v in ipairs(vars) do
			lines:insert(assignFromRepl(v))
		end
		return lines:concat'\n'
	end

	function assign_real3(name, expr)
		if expr == nil then expr = casEnv[name] end
		local lines = table()
		lines:insert('\treal3 '..name..';')
		for i,xi in ipairs(xNames) do
			assert(expr and expr[i], "failed to find expr for var "..name)
			lines:insert('\t'..name..'.'..xi..' = '..compile(expr[i])..';')
		end
		return lines:concat'\n'
	end

	function assign_sym3(name, expr)
		if expr == nil then expr = casEnv[name] end
		local lines = table()
		lines:insert('\tsym3 '..name..';')
		for ij,xij in ipairs(symNames) do
			local i,j = from6to3x3(ij)
			assert(expr and expr[i] and expr[i][j], "failed to find expr for var "..name)
			lines:insert('\t'..name..'.'..xij..' = '..compile(expr[i][j])..';')
		end
		return lines:concat'\n'
	end

	function assign_real3x3(name, expr)
		if expr == nil then expr = casEnv[name] end
		local lines = table()
		lines:insert('\treal3x3 '..name..';')
		for i,xi in ipairs(xNames) do
			for j,xj in ipairs(xNames) do
				assert(expr and expr[i] and expr[i][j], "failed to find expr for var "..name)
				lines:insert('\t'..name..'.'..xi..'.'..xj..' = '..compile(expr[i][j])..';')
			end
		end
		return lines:concat'\n'
	end

	-- assumes the 2nd and 3rd indexes are symmetric
	function assign_3sym3(name, expr)
		if expr == nil then expr = casEnv[name] end
		local lines = table()
		lines:insert('\t_3sym3 '..name..';')
		for i,xi in ipairs(xNames) do
			for jk,xjk in ipairs(symNames) do
				local j,k = from6to3x3(jk)
				assert(expr and expr[i] and expr[i][j] and expr[i][j][k], "failed to find expr for var "..name)
				lines:insert('\t'..name..'.'..xi..'.'..xjk..' = '..compile(expr[i][j][k])..';')
			end
		end
		return lines:concat'\n'
	end

	-- assumes 2nd and 3rd indexes are symmetric
	-- places the 4th index first, as an array
	function assign_3sym3x3(name, expr)
		if expr == nil then expr = casEnv[name] end
		local lines = table()
		lines:insert('\t_3sym3 '..name..'[3];')
		for i,xi in ipairs(xNames) do
			for jk,xjk in ipairs(symNames) do
				local j,k = from6to3x3(jk)
				for l,xl in ipairs(xNames) do
					assert(expr and expr[i] and expr[i][j] and expr[i][j][k] and expr[i][j][k][l], "failed to find expr for var "..name)
					lines:insert('\t'..name..'['..(l-1)..'].'..xi..'.'..xjk..' = '..compile(expr[i][j][k][l])..';')
				end
			end
		end
		return lines:concat'\n'
	end

	function assign_sym3sym3(name, expr)
		if expr == nil then expr = casEnv[name] end
		local lines = table()
		lines:insert('\tsym3sym3 '..name..';')
		for ij,xij in ipairs(symNames) do
			local i,j = from6to3x3(ij)
			for kl,xkl in ipairs(symNames) do
				local k,l = from6to3x3(kl)
				assert(expr and expr[i] and expr[i][j] and expr[i][j][k] and expr[i][j][k][l], "failed to find expr for var "..name)
				lines:insert('\t'..name..'.'..xij..'.'..xkl..' = '..compile(expr[i][j][k][l])..';')
			end
		end
		return lines:concat'\n'
	end


	-- from here on out is stuff specific to different functions 


	gammaHat_ll = Tensor.metric().metric
	
	det_gammaHat = Matrix.determinant(gammaHat_ll)
	det_gammaHat_var = compileVar('det_gammaHat', coords)

	gammaHat_uu = Tensor('^ij', table.unpack((Matrix.inverse(gammaHat_ll, nil, nil, nil, 
		-- defer det_gammaHat
		--det_gammaHat_var))))
		-- don't bother.  gammaHat^ij is just an inverse scale matrix
		det_gammaHat))))

	partial_gammaHat_lll = gammaHat_ll'_ij,k'()
	partial2_gammaHat_llll = partial_gammaHat_lll'_ijk,l'()
	connHat_lll = ((partial_gammaHat_lll'_ijk' + partial_gammaHat_lll'_ikj' - partial_gammaHat_lll'_jki')/2)()
	connHat_ull = (gammaHat_uu'^im' * connHat_lll'_mjk')()
	partial_connHat_ulll = connHat_ull'^i_jk,l'() 
	partial_det_gammaHat_l = Tensor('_i', function(i) 
		return det_gammaHat:diff(coords[i])() 
	end)
	partial_det_gammaHat_l_vars = Tensor('_i', function(i) 
		return compileReplVar(
			det_gammaHat_var:diff(coords[i])(),
			'partial_det_gammaHat_l.'..xNames[i], 
			coords)
	end)
	partial2_det_gammaHat_ll = partial_det_gammaHat_l'_i,j'()

	-- non-coordinate basis
	-- TODO either remove the default metric, or introduce a separate set of indexes for the non-coordinate basis
	-- either one to prevent symmath from automatically raising or lowering
	-- otherwise using e or eu incorrectly will result in extra operations

	e = Tensor('_i^I', function(i,j)
		return (i==j and sqrt(gammaHat_ll[i][i])() or 0)
	end)
	eu = Tensor('^i_I', function(i,j)
		return i==j and (1/e[i][i])() or 0
	end)

	-- state variables

	alpha = compileVar('U->alpha', coords)
	beta_U_vars = Tensor('^I', function(I)
		return compileVar('U->beta_U.'..xNames[I], coords)
	end)
	epsilon_LL_vars = Tensor('_IJ', function(I,J)
		return compileVar('U->epsilon_LL.'..sym(I,J), coords)
	end)
	W = compileVar('U->W', coords)
	K = compileVar('U->K', coords)
	ABar_LL_vars = Tensor('_IJ', function(I,J)
		return compileVar('U->ABar_LL.'..sym(I,J), coords)
	end)
	LambdaBar_U_vars = Tensor('^I', function(I)
		return compileVar('U->LambdaBar_U.'..xNames[I], coords)
	end)
	B_U_vars = Tensor('^I', function(I)
		return compileVar('U->B_U.'..xNames[I], coords)
	end)
	
	rho = compileVar'U->rho'
	
	pi = compileVar'M_PI'

	-- state variable derivatives (variables created with 'makePartial')

	partial_alpha_l_vars = Tensor('_i', function(i) return compileReplVar(alpha:diff(coords[i]), 'partial_alpha_l['..(i-1)..']', coords) end)
	partial2_alpha_ll_vars = Tensor('_ij', function(i,j)
		local ij = from3x3to6(i,j)	
		return compileReplVar(alpha:diff(coords[i], coords[j]), 'partial2_alpha_ll['..(ij-1)..']', coords) 
	end)
	partial_epsilon_LLl_vars = Tensor('_IJk', function(I,J,k)
		return compileReplVar(epsilon_LL_vars[I][J]:diff(coords[k]), 'partial_epsilon_LLl['..(k-1)..'].'..sym(I,J), coords)
	end)
	partial2_epsilon_LLll_vars = Tensor('_IJkl', function(I,J,k,l)
		local kl = from3x3to6(k,l)
		return compileReplVar(epsilon_LL_vars[I][J]:diff(coords[k], coords[l]), 'partial2_epsilon_LLll['..(kl-1)..'].'..sym(I,J), coords)
	end)
	partial_K_l_vars = Tensor('_i', function(i) return compileReplVar(K:diff(coords[i]), 'partial_K_l['..(i-1)..']', coords) end)
	partial_W_l_vars = Tensor('_i', function(i) return compileReplVar(W:diff(coords[i]), 'partial_W_l['..(i-1)..']', coords) end)
	partial_LambdaBar_Ul_vars = Tensor('^I_j', function(I,j)
		return compileReplVar(LambdaBar_U_vars[I]:diff(coords[j]), 'partial_LambdaBar_Ul['..(j-1)..'].'..xNames[I], coords)
	end)
	partial_ABar_LLl_vars = Tensor('_IJk', function(I,J,k)
		return compileReplVar(ABar_LL_vars[I][J]:diff(coords[k]), 'partial_ABar_LLl['..(k-1)..'].'..sym(I,J), coords)
	end)
	partial_beta_Ul_vars = Tensor('^I_j', function(I,j)
		return compileReplVar(beta_U_vars[I]:diff(coords[j]), 'partial_beta_Ul['..(j-1)..'].'..xNames[I], coords)
	end)
	partial2_beta_Ull_vars = Tensor('^I_jk', function(I,j,k)
		local jk = from3x3to6(j,k)
		return compileReplVar(beta_U_vars[I]:diff(coords[j], coords[k]), 'partial2_beta_Ull['..(jk-1)..'].'..xNames[I], coords)
	end)
	partial_B_Ul_Vars = Tensor('^I_j', function(I,j)
		return compileReplVar(B_U_vars[I]:diff(coords[j]), 'partial_B_Ul['..(j-1)..'].'..xNames[I], coords)
	end)

	-- state variables in coordinate form

	epsilon_ll = (epsilon_LL_vars'_IJ' * e'_i^I' * e'_j^J')()
	LambdaBar_u = (LambdaBar_U_vars'^I' * eu'^i_I')()
	ABar_ll = (ABar_LL_vars'_IJ' * e'_i^I' * e'_j^J')()
	beta_u = (beta_U_vars'^I' * eu'^i_I')()
	B_u = (B_U_vars'^I' * eu'^i_I')()

	-- derivatives in coordinate form
		-- without rescaling
	
	partial_alpha_l = Tensor('_i', function(i) return alpha:diff(coords[i]) end)
	partial2_alpha_ll = partial_alpha_l'_i,j'()
	partial_W_l = Tensor('_i', function(i) return W:diff(coords[i]) end)
	partial_K_l = Tensor('_i', function(i) return K:diff(coords[i]) end)

		-- with rescaling (these will be slower and more complex)

	partial_epsilon_lll = epsilon_ll'_ij,k'()
	partial_ABar_lll = ABar_ll'_ij,k'() 
	partial_LambdaBar_ul = LambdaBar_u'^i_,j'()
	partial_beta_ul = beta_u'^i_,j'()
	partial2_beta_ull = partial_beta_ul'^i_j,k'()
	partial_B_ul = B_u'^i_,j'()

	-- derived values

	S = compileVar'S'
	exp_neg4phi = W * W
	
	gammaBar_ll = (gammaHat_ll'_ij' + epsilon_ll'_ij')()
	gammaBar_LL = (eu'^i_I' * eu'^j_J' * gammaBar_ll'_ij')()

	det_gammaBar = Matrix.determinant(gammaBar_ll)

	-- factor out the coord metric det: r^4 sin(theta)^2 of spherical
	-- but leave the rest as a variable
	det_gammaBar_over_det_gammaHat = (det_gammaBar / det_gammaHat)()
	det_gammaBar_over_det_gammaHat_var = compileVar('det_gammaBar_over_det_gammaHat', coords)  

	-- defer det_gammaHat
	--det_gammaBar = (det_gammaBar_over_det_gammaHat * det_gammaHat_var)()
	-- defer det_gammaBar/det_gammaHat
	det_gammaBar = (det_gammaBar_over_det_gammaHat_var * det_gammaHat)() 
	-- don't
	det_gammaBar = det_gammaBar
	
	det_gammaBar_var = compileVar('det_gammaBar', coords)

	gammaBar_uu_vars = Tensor('^ij', function(i,j)
		return compileVar('gammaBar_uu.'..sym(i,j), coords)
	end)
	gammaBar_uu = Tensor('^ij', table.unpack((Matrix.inverse(gammaBar_ll, nil, nil, nil, 
		-- defer det_gammaBar
		--det_gammaBar_var))))
		-- defer det_gammaBar_over_det_gammaHat (but leave det_gammaHat for further simplification)
		det_gammaBar_over_det_gammaHat_var * det_gammaHat))))
		-- don't
		--det_gammaBar))))
	
	-- this isn't always completely effective.  sometimes it introduces sin(theta)'s in the denominator
	--gammaBar_uu = (gammaBar_uu / det_gammaHat * det_gammaHat_var)()

	partial_gammaBar_lll = gammaBar_ll'_ij,k'()
	partial2_gammaBar_llll = partial_gammaBar_lll'_ijk,l'()
	trBar_partial2_gammaBar_ll = (gammaBar_uu'^kl' * partial2_gammaBar_llll'_ijkl')()

	connBar_lll = ((partial_gammaBar_lll'_ijk' + partial_gammaBar_lll'_ikj' - partial_gammaBar_lll'_jki') / 2)()

	-- defer gammaBar^ij
	--connBar_ull = (gammaBar_uu_vars'^im' * connBar_lll'_mjk')()
	-- don't
	connBar_ull = (gammaBar_uu'^im' * connBar_lll'_mjk')()

	-------------------------------- alpha_,t -------------------------------- 

	-- usually in init/einstein
	local fGuiVar = solver.eqn.guiVars.f_eqn
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alpha, symmath)
	f = symmath.clone(f)
	
	--Alcubierre 4.2.52 - Bona-Masso family of slicing
	Q = f * K
	
	--d/dt alpha alpha,t - alpha_,i beta^i = -alpha^2 Q
	-- already seeing a 1/r and 1/(r sin(θ)) in the denom...
	-- but it is only next to beta^θ and beta^φ
	--simplification is important for the 1/alpha's at least
	dt_alpha = (-alpha^2 * Q + partial_alpha_l'_i' * beta_u'^i')():factorDivision()
printbr(dt_alpha)
--printbr(assign'dt_alpha')

	-------------------------------- W_,t -------------------------------- 

printbr'tr_connBar_l'	
	-- using conn trace
	--tr_connBar_l = connBar_ull'^j_kj'()
	-- using det_gammaBar_,k / det_gammaBar
	-- when deferring det_gammaBar/det_gammaHat we get the terms [4/r, 2 cot(th), 0] out front, and all else has 1/(det gammaBar/det gammaHat)
	tr_connBar_l = Tensor('_i', function(i) 
		local expr = (sqrt(det_gammaBar):diff(coords[i]) / sqrt(det_gammaBar))()
		-- [[ expand (det gammaBar/det gammaHat)_,i
		for j=1,3 do
			expr = expr:replace(
				det_gammaBar_over_det_gammaHat_var:diff(coords[j]),
				det_gammaBar_over_det_gammaHat:diff(coords[j])()
			)
		end
		--]]
		expr = expr:factorDivision() 
		return expr
	end)
	tr_connBar_l_vars = Tensor('_i', function(i)
		return compileVar('tr_connBar_l.'..xNames[i], coords)
	end)
printbr(tr_connBar_l)
--printbr(assign_real3'tr_connBar_l')
printbr'tr_DBar_beta' 
	-- if you do defer gammaBar^ij then you get a 1/(r^6 sin(theta)^6) in the denom
	-- if you don't defer gammaBar^ij but do defer det_gammaBar then you get a 1/(det_gammaBar r^3 sin(theta))^3
	-- if you don't defer anything then you get ... the same thing, but expanded out, of course 
	-- if you defer (det_gammaBar/det_gammaHat) then this still allows the (r^4 sin(theta)^2) of det_gammaHat to simplify out
	--  and you get lots of 1/(det gammaBar/det gammaHat) terms (which are non-singular) 
	--  and a few 1/r terms next to epsilon_IJ (non-coord)
	--  and a few 1/(r sin(theta))'s next to epsilon_IJ square terms
	-- if you defer (det_gammaBar/det_gammaHat) and use sqrt(gbar)_,i/sqrt(gbar)=conn^j_ij for tr_connBar_l then ...
	--  you get a 1/r in front of epsilon_IJ,theta
	--   and a 1/(r sin(theta)) in front of epsilon_IJ,phi
	tr_DBar_beta = (partial_beta_ul'^j_j' + tr_connBar_l_vars'_j' * beta_u'^j')():factorDivision()
	tr_DBar_beta_var = compileVar('tr_DBar_beta', coords)
printbr(tr_DBar_beta)
--printbr(assign'tr_DBar_beta')
printbr(assign'tr_DBar_beta')
	--2017 Ruchlin et al eqn 11c
	--W,t = 1/3 W (alpha K - beta^k connBar^j_kj - beta^k_,k) + beta^k W_,k
	-- if you defer (det_gammaBar/det_gammaHat) then you get ....
	--  beta^X/r * ... + beta^Y/r * ...
	-- so all the 1/r terms are next to beta's
printbr'dt_W'
	-- terms that aren't affected by anything: 
	--  1/3 W (alpha K - beta^R_,r)
	-- terms scaled by 1/(det gammaBar/det gammaHat):
	--	some of -1/3 W beta^I connBar^J_IJ
	dt_W = (frac(1,3) * W * (alpha * K - tr_DBar_beta_var) + beta_u'^k' * partial_W_l'_k')():factorDivision()
printbr(dt_W)
printbr(assign'dt_W')
	-------------------------------- K_,t -------------------------------- 
	
	-- TODO just set the metric to gammaBar
	-- or set some indexes to gammaBar and others to gamma? and yet others to non-coord gammaBar?

printbr('ABar_ul')
	ABar_ul = (gammaBar_uu'^ik' * ABar_ll'_kj')():factorDivision()
	ABar_ul_vars = Tensor('^i_j', function(i,j) return compileVar('ABar_ul.'..xNames[i]..'.'..xNames[j], coords) end)
printbr(ABar_ul)
printbr('ABarSq_ul')
	ABarSq_ul = (ABar_ul_vars'^i_k' * ABar_ul_vars'^k_j')():factorDivision()
	ABarSq_ul_vars = Tensor('^i_j', function(i,j) return compileVar('ABarSq_ul.'..xNames[i]..'.'..xNames[j], coords) end)
printbr(ABarSq_ul)
printbr('tr_ABarSq')
	tr_ABarSq = (ABarSq_ul_vars'^k_k')():factorDivision()
	tr_ABarSq_var = compileVar('tr_ABarSq', coords) 
printbr(tr_ABarSq)

	connBar_ull_vars = Tensor('^i_jk', function(i,j,k)
		return compileVar('connBar_ull.'..xNames[i]..'.'..sym(j,k), coords)
	end)
printbr('DBar2_alpha_ll')
	DBar2_alpha_ll = (partial2_alpha_ll'_ij' - connBar_ull_vars'^k_ij' * partial_alpha_l'_k')():factorDivision()
printbr(DBar2_alpha_ll)
	DBar2_alpha_ll_vars = Tensor('_ij', function(i,j)
		return compileVar('DBar2_alpha_ll.'..sym(i,j), coords)
	end)
printbr('trBar_DBar2_alpha')
	trBar_DBar2_alpha = (gammaBar_uu_vars'^ij' * DBar2_alpha_ll_vars'_ij')():factorDivision()
	trBar_DBar2_alpha_var = compileVar('trBar_DBar2_alpha', coords)
printbr(trBar_DBar2_alpha)

printbr('partial_alpha_u')
	partial_alpha_u = (gammaBar_uu_vars'^ij' * partial_alpha_l'_j')():factorDivision()
printbr(partial_alpha_u)

	-- W = exp(-2 phi)
	-- phi = -log(W)/2
	--W_for_phi = W:eq(exp(-2*phi_var))
	--phi = W_for_phi:solve(phi_var)
	phi = (-symmath.log(W)/2):factorDivision()
	partial_phi_l = Tensor('_i', function(i) 
		return phi:diff(coords[i])():factorDivision() 
	end)
	partial_phi_l_vars = Tensor('_i', function(i) return compileVar('partial_phi_l.'..xNames[i], coords) end)

	--[[
	B&S 11.52
	Alcubierre 2.8.12
	K_,t = -gamma^ij D_i D_j alpha + alpha (ABar_ij ABar^ij + K^2 / 3) + 4 pi alpha (rho + S) + beta^i K_,i
	2017 Ruchlin et al
	K_,t = 
		1/3 alpha K^2 
		+ alpha ABar_ij ABar^ij 
		- exp(-4 phi) (
			DBar^i DBar_i alpha 
			+ 2 gammaBar^ij alpha_,i phi_,j
		) 
		+ K_,i beta^i
		+ 4 pi alpha (rho + S)
	
	what is immediately apparently safe from 1/r?
		1/3 alpha K^2
		+ exp(-4 phi) gammaBar^ij (-alpha_,ij + 2 alpha_,i phi_,j)
		+ 4 pi alpha (rho + S)
	what could be problematic?
		alpha ABar_ij ABar^ij
		+ exp(-4phi) * connections
		+ K_,i beta^i
	--]]
printbr('dt_K')
	dt_K = (frac(1,3) * alpha * K^2
		+ alpha * tr_ABarSq_var
		- exp_neg4phi * (
			trBar_DBar2_alpha_var
			+ 2 * partial_alpha_u'^i' * partial_phi_l'_i'
		)
		+ partial_K_l'_i' * beta_u'^i'
		+ 4 * pi * alpha * (rho + S)
	)():factorDivision()
printbr(dt_K)
	-------------------------------- epsilon_ij,t -------------------------------- 

	-- e^i_I e^j_J (gammaBar_ij,k beta^k + gammaBar_ki * beta^k,_j + gammaBar_kj * beta^k,_i)

printbr'Lbeta_gammaBar_LL'
	Lbeta_gammaBar_LL = (eu'^i_I' * eu'^j_J' * (
			beta_u'^k' * partial_gammaBar_lll'_ijk'
			+ gammaBar_ll'_ki' * partial_beta_ul'^k_j'
			+ gammaBar_ll'_kj' * partial_beta_ul'^k_i'
		))():factorDivision()
printbr(Lbeta_gammaBar_LL)
	LBeta_gammaBar_LL_vars = Tensor('_IJ', function(I,J) return compileVar('LBeta_gammaBar_LL.'..sym(I,J), coords) end)  

	--[[
	2017 Ruchlin et al, eqn 11a
	epsilon_ij,t = 2/3 gammaBar_ij (alpha ABar^k_k - DBar_k beta^k) + DHat_i beta_j + DHat_j beta_i - 2 alpha ABar_ij + epsilon_ij,k beta^k + epsilon_ik beta^k_,j + epsilon_kj beta^k_,i
	...using DBar_(i beta_j) = DHat_(i beta_j) + epsilon_k(j beta^k_,i) + 1/2 epsilon_ij,k beta^k
	= 	
		+ 2/3 gammaBar_ij (
			alpha ABar^k_k 
			- beta^k_,k 
			- connBar^k_lk beta^l
		) 
		- 2 alpha ABar_ij 
		+ .5 DBar_i beta_j
		+ .5 DBar_j beta_i
	

	Etienne's SENR Mathematica notebook:
	= 
		// Lie derivative terms
		beta^k gammaBar_ij,k
		+ gammaBar_ki beta^k_,j
		+ gammaBar_kj beta^k_,i

		- 2/3 gammaBar_ij DBar_k beta^k
		(notice that "2/3 alpha ABar^k_k" is excluded)
		- 2 alpha ABar_ij

	Looks like the paper's notation DBar_i beta_j implies DBar_j (gammaBar_ik beta^k) = gammaBar_ki DBar_j beta^k
	...which is not DBar_j ( gamma_ik beta^k ), which was my interpretation of the definition of beta_i := gamma_ij beta^k
	(denoting indexes for coordinate transforms is easier to interpret than denoting entire tensors for coordinate transforms)
	--]]
printbr'dt_epsilon_LL'
	dt_epsilon_LL = (
		Lbeta_gammaBar_LL'_IJ'
		- frac(2,3) * gammaBar_LL'_IJ' * tr_DBar_beta
		- 2 * alpha * ABar_LL_vars'_IJ'
	)():factorDivision()
printbr(dt_epsilon_LL)
--]=]	


	-- [[ do this every time you stop using the casEnv
	setfenv(1, oldEnv)
	--]]
	
	
	return casEnv
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template([[

//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero

sym3 calc_gammaHat_ll(real3 x) {
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>
	return gammaHat_ll;
}

real calc_det_gammaHat(real3 x) {
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign'det_gammaHat'?>
	return det_gammaHat;
}

sym3 calc_gammaHat_uu(real3 x) {
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign'det_gammaHat'?>
<?=assign_sym3'gammaHat_uu'?>
	return gammaHat_uu;
}

//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign_sym3'gammaBar_ll'?>
	return gammaBar_ll;
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//however det(gammaHat_ij) == det(gammaBar_ij) by the eqn just before (6) in 2017 Ruchlin
real calc_det_gammaBar(real3 x) {
	//TODO detg ...
	real det_gammaHat = calc_det_gammaHat(x);
	//real detg = 1.;
	//real det_gammaBar = det_gammaHat * detg;
	//return det_gammaBar;
	return det_gammaHat;
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

]], self:getCASEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_RBar_ll()
	local casEnv = self:getCASEnv()
--[=[ prereq code for getCode_RBar_ll
-- TODO some kind of request / prefix system so I don't have to gen the cos()'s and sin()'s that I don't need
-- have it work into the coord/coord stuff so I don't gen extra stuff there either
-- (or alternatively just move all the stuff added to coord/coord just for into getCASEnv)
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign'det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
--]=]
	return template([[

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
	partial_DHat_gammaBar_llll_minus_one_term[l].k.ij := partial_l DHat_k gammaBar_ij
		= (gammaBar_ij,k - connHat^m_ki gammaBar_mj - connHat^m_kj gammaBar_mi)_,l
		= gammaBar_ij,kl 
			- connHat^m_ki,l gammaBar_mj 
			- connHat^m_kj,l gammaBar_mi
			- connHat^m_ki gammaBar_mj,l
			- connHat^m_kj gammaBar_mi,l
	*/
	_3sym3 partial_DHat_gammaBar_llll_minus_one_term[3];
<? 
for k,xk in ipairs(xNames) do
	for l,xl in ipairs(xNames) do
		for ij,xij in ipairs(symNames) do
			local i,j = from6to3x3(ij)
			local xi,xj = xNames[i], xNames[j]
?>	partial_DHat_gammaBar_llll_minus_one_term[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
//		+ partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>	// moved into the above gen'd code 
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
	DHat2_gammaBar_llll_minus_one_term[l].k.ij = DHat_l DHat_k gammaBar_ij
		= partial_l DHat_k gammaBar_ij
			- connHat^m_lk DHat_m gammaBar_ij
			- connHat^m_li DHat_k gammaBar_mj
			- connHat^m_lj DHat_k gammaBar_im
	*/
	_3sym3 DHat2_gammaBar_llll_minus_one_term[3];
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i], xNames[j]
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>	DHat2_gammaBar_llll_minus_one_term[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
		+ partial_DHat_gammaBar_llll_minus_one_term[<?=l-1?>].<?=xk?>.<?=xij?>	//diverging, such that in spherical vacuum the trace of these is diag(0, 1, .5)
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
		+ trBar_partial2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * DHat2_gammaBar_llll_minus_one_term[<?=l-1?>].<?=xk?>.<?=xij?>
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
?>	DHat_LambdaBar_ul.<?=xi?>.<?=xj?> = partial_LambdaBar_ul.<?=xi?>.<?=xj?>
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
	return template([=[
<?
-- [[ do this every time you use the casEnv
local casEnv = eqn:getCASEnv()
local oldEnv = getfenv()
getmetatable(casEnv).__index = oldEnv
setfenv(1, casEnv)
--]]
?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>

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
//TODO do this symbolically.  That's what I originally did, but symbolic calculations were getting complex
// however, with spherical BSSN, you need to 
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

<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=eqn:makePartial'epsilon_LL'?>

<?=assign'det_gammaHat'?>
<?=assign_3sym3'connHat_ull'?>

<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_ull'?>
	
	//Delta^i_jk = connBar^i_jk - connHat^i_jk
	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);

	real3 LambdaBar_u;
	LambdaBar_u = _3sym3_sym3_dot23(Delta_ull, gammaBar_uu);

	U->LambdaBar_U = real3_rescaleFromCoord_u(LambdaBar_u, x);
}
<?
-- [[ do this every time you stop using the casEnv
setfenv(1, oldEnv)
--]]
?>
]=], {
		eqn = self,
		code = self.initState:initState(self.solver),
	})
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd.cl'

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
-- ]=]	
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
-- [=[	
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
--]=]	
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

--[[
	'U RBar_ll xx',
	'U RBar_ll xy',
	'U RBar_ll xz',
	'U RBar_ll yy',
	'U RBar_ll yz',
	'U RBar_ll zz',
--]]
--[[
	'U del gammaBar_ll sym xx',
	'U del gammaBar_ll sym xy',
	'U del gammaBar_ll sym xz',
	'U del gammaBar_ll sym yy',
	'U del gammaBar_ll sym yz',
	'U del gammaBar_ll sym zz',
--]]
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
		{name='gammaHat_ll', code=[[	*value_sym3 = calc_gammaHat_ll(x);]], type='sym3'},
		{name='gammaBar_ll', code=[[	*value_sym3 = calc_gammaBar_ll(U, x);]], type='sym3'},
		{name='gamma_uu', code=[[	*value_sym3 = calc_gamma_uu(U, x);]], type='sym3'},
		{name='gammaHat_uu', code=[[	*value_sym3 = calc_gammaHat_uu(x);]], type='sym3'},
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
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	sym3 ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, ABar_UL);
	*value_sym3 = ABarSq_LL;
]],
		},
		
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = template([[

<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_ull'?>

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
]], self:getCASEnv())
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
]], self:getCASEnv()),
		},

		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'alpha'?>
	*value_real3 = *(real3*)partial_alpha_l;
]], self:getCASEnv()),
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

<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_ull'?>
	
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
]], self:getCASEnv()),
		},

		{
			name = 'DBar2_alpha_ll',
			code = template([[
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>
	
<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_ull'?>

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
]], self:getCASEnv()),
		},
	
	
		{
			name = 'tracelessPart_ll',
			type = 'sym3',
			code = template([[
	
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>

<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_ull'?>

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
]], self:getCASEnv()),
		},
--]=]	
	}

	local derivOrder = 2 * self.solver.numGhost
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
<?=eqn:makePartial'W'?>
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial'beta_U'?>
<?=assign_real3x3'partial_beta_ul'?>
	real tr_partial_beta = real3x3_trace(partial_beta_ul);

	real exp_4phi = 1. / calc_exp_neg4phi(U);

<?=assign_sym3'gammaBar_ll'?>

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);

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
]], 		applyCommon{
				eqn = self,
				solver = self.solver,
			}
		)
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
		code= template([[
<?=eqn:makePartial'alpha'?>

	real _1_alpha = 1. / U->alpha;

<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

<?=eqn:makePartial'beta_U'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
<?=assign_sym3'gammaBar_ll'?>
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
<?=eqn:makePartial'W'?>
<?=assign_3sym3('partial_gamma_lll', partial_gamma_lll:permute'_kij'?>

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
	real3 beta_dbeta_l;
<? for i,xi in ipairs(xNames) do
?>	beta_dbeta_l.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
?>		+ beta_l.<?=xj?> * partial_beta_ul.<?=xj?>.<?=xi?>
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
]],			applyCommon{
				eqn = self,
				solver = self.solver,
			}
		), 
		type = 'real3',
	}
--]=]
-- [=[
	do
		-- [[ do this every time you use the casEnv
		local casEnv = self:getCASEnv()
		local oldEnv = getfenv()
		getmetatable(casEnv).__index = oldEnv
		setfenv(1, casEnv)
		--]]	

		vars:insert{
			name = 'RBar_ll',
			type = 'sym3',
			code = template([[
<?=eqn:makePartial'LambdaBar_U'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>

<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>

<?=assign_sym3'gammaBar_ll'?>
<?=assign_sym3'gammaBar_uu'?>

<?=assign'det_gammaHat'?>
<?=assign_3sym3'connHat_ull'?>
<?=assign_3sym3'connBar_ull'?>

<?=assign_real3'LambdaBar_u'?>
<?=assign_real3x3'partial_LambdaBar_ul'?>
	
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);
	real3 Delta_u = real3_rescaleToCoord_U(Delta_U, x);

	_3sym3 Delta_ull = _3sym3_sub(connBar_ull, connHat_ull);
	_3sym3 Delta_lll = sym3_3sym3_mul(gammaBar_ll, Delta_ull);

<?=assign_sym3'gammaHat_uu'?>

<?=assign_real3'partial_det_gammaHat_l'?>
<?=assign_3sym3x3'partial_connHat_ulll'?>

<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>

	sym3 RBar_ll;
<?=eqn:getCode_RBar_ll()?>

	*value_sym3 = RBar_ll;
]], casEnv),
		}
	
		-- [[ do this every time you stop using the casEnv
		setfenv(1, oldEnv)
		--]]
	end
--]=]		
--[=[		
	vars:insert{
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
	
<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>

<?=assign_3sym3'connBar_ull'?>

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
]],			applyCommon{
				eqn = self,
				solver = self.solver,
			}
		),
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

	do
		-- [[ do this every time you use the casEnv
		local casEnv = self:getCASEnv()
		local oldEnv = getfenv()
		getmetatable(casEnv).__index = oldEnv
		setfenv(1, casEnv)
		--]]	

		-- symbolically generate ... so far gammaBar^ij	
		-- this fixes it.
		vars:insert{
			name='del gammaBar_ll sym',
			type = 'sym3',
			code = template([[
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>
<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
	*value_sym3 = trBar_partial2_gammaBar_ll;
]], casEnv),
		}
		-- [[ do this every time you stop using the casEnv
		setfenv(1, oldEnv)
		--]]
	end

	vars:insert{
		name = 'tr34 (gamma*dGamma)',
		type = 'real3x3',
		code = template([[
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>
<?=assign_real3'partial_det_gammaHat_l'?>
<?=assign_3sym3x3'partial_connHat_ulll'?>
	
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
]], self:getCASEnv()),
	}

	vars:insert{
		name = 'tr14 (Gamma*dgamma)',
		type = 'real3x3',
		code = template([[
<?=assignAllFromRepl(cos_xs)?>
<?=assignAllFromRepl(sin_xs)?>
<?=eqn:makePartial'epsilon_LL'?>

<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij') 
-- how does :permute() work? forward or backward? 
?>

<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>

<?=assign'det_gammaHat'?>
<?=assign_3sym3'connHat_ull'?>
		
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
]], self:getCASEnv()),
	}

	--[[ hmm? not working.
	vars:insert{name='x', type='real3', code='*value_real3=x;'}
	--]]
	-- [[
	vars:insert{name='x', code='*value=x.x;'}
	vars:insert{name='y', code='*value=x.y;'}
	vars:insert{name='z', code='*value=x.z;'}
	--]]

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
