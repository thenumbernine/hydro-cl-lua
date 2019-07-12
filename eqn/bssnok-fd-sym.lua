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

		--{name='constrain_tr_ABar', value=true, compileTime=true},
		{name='constrain_tr_ABar', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},

		{name='shift_eta', value=1},	--1, or 1 / (2 M), for total mass M
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
function BSSNOKFiniteDifferenceEquation:getEnv()
	if self.env then return self.env end

	local outfile = assert(io.open('bssn_output.html', 'w'))

time('building symbolic math env...', function()
	
	local derivOrder = 2 * self.solver.numGhost
	local env = applyCommon{
		eqn = self,
		solver = self.solver,
	}
	setmetatable(env, {})
	self.env = env

	-- NOTICE if there is a local variable in the outer scope
	-- then it will not get set in the nested new env
	-- [[ do this every time you use the env
	local oldEnv = getfenv()
	getmetatable(env).__index = oldEnv
	setfenv(1, env)
	--]]

	-- this is the local env code for env
	symmath = require 'symmath'
-- [[ mathjax output
	symmath.tostring = symmath.export.MathJax
	symmath.tostring.useCommaDerivative = true
	outfile:write((tostring(symmath.tostring.header):gsub('tryToFindMathJax.js', 'tryToFindMathJax_andDeferMathRendering.js')))
	local idnum = 1
	local function convertSingle(arg)
		local s = tostring(arg)
		if symmath.Expression.is(arg) then
		-- TODO proper way is to arg:replace() everything
		-- but that is slow
			s = s:gsub('TF_DBar2_alpha_LL%.(.)(.)', '(\\bar{D}_{%1} \\bar{D}_{%2} \\alpha)^{TF}')
			s = s:gsub('tracelessPart_LL%.(.)(.)', '(traceless)_{\\hat{%1} \\hat{%2}}')
			s = s:gsub('U%->alpha', '\\alpha')
			s = s:gsub('partial_alpha_l%[(.)%]', function(i) return '\\alpha_{,'..xNames[i+1]..'}' end)
			s = s:gsub('partial2_alpha_ll%[(.)%]', function(ij) return '\\alpha_{,'..sym(from6to3x3(ij+1))..'}' end)
			s = s:gsub('U%->beta_U%.(.)', '\\beta^{\\hat{%1}}')
			s = s:gsub('partial_beta_Ul%[(.)%]%.(.)', function(j) return '{\\beta^{\\hat{%2}}}_{,'..xNames[j+1]..'}' end)
			s = s:gsub('partial2_beta_Ul%[(.)%]%.(.)', function(jk) return '{\\beta^{\\hat{%2}}}_{,'..sym(from6to3x3(jk+1))..'}' end)
			s = s:gsub('U%->epsilon_LL%.(.)(.)', '\\epsilon_{\\hat{%1}\\hat{%2}}')
			s = s:gsub('partial_epsilon_LLl%[(.)%]%.(.)(.)', function(k) return '\\epsilon_{\\hat{%2}\\hat{%3},'..xNames[i+1]..'}' end)
			s = s:gsub('partial2_epsilon_LLll%[(.)%]%.(.)(.)', function(kl) return '\\epsilon_{\\hat{%2}\\hat{%3},'..sym(from6to3x3(kl+1))..'}' end)
			s = s:gsub('U%->W', 'W')
			s = s:gsub('partial_W_l%[(.)%]', function(i) return 'W_{,'..xNames[i+1]..'}' end)
			s = s:gsub('partial2_W_ll%[(.)%]', function(ij) return 'W_{,'..sym(from6to3x3(ij+1))..'}' end)
			s = s:gsub('U%->K', 'K')
			s = s:gsub('partial_K_l%[(.)%]', function(i) return 'K_{,'..xNames[i+1]..'}' end)
			s = s:gsub('U%->ABar_LL%.(.)(.)', '\\bar{A}_{\\hat{%1}\\hat{%2}}')
			s = s:gsub('partial_ABar_LLl%[(.)%]%.(.)(.)', function(k) return '\\bar{A}_{\\hat{%2}\\hat{%3},'..xNames[k+1]..'}' end)
			s = s:gsub('Lbeta_ABar_LL%.(.)(.)', '\\mathcal{L}_{\\beta} \\bar{A}_{\\hat{%1}\\hat{%2}}')
			s = s:gsub('U%->LambdaBar_U%.(.)', '\\bar{\\Lambda}^{\\hat{%1}}')
			s = s:gsub('partial_LambdaBar_Ul%[(.)%]%.(.)', function(j) return '{\\bar{\\Lambda}^{\\hat{%2}}}_{,'..xNames[j+1]..'}' end)
			s = s:gsub('dt_LambdaBar_U%.(.)', '{\\bar{\\Lambda}^{%1}}_{,t}')
			s = s:gsub('U%->B_U%.(.)', 'B^{\\hat{%1}}')
			s = s:gsub('partial_B_Ul%[(.)%]%.(.)', function(j) return '{B^{\\hat{%2}}}_{,'..xNames[j+1]..'}' end)
			s = s:gsub('U%->rho', '\\rho')
			s = s:gsub('M_PI', '\\pi')
			s = s:gsub('partial2_det_gammaHat_ll%.(..)', '\\hat{\\gamma}_{,%1}')
			s = s:gsub('partial_det_gammaHat_l%.(.)', '\\hat{\\gamma}_{,%1}')
			s = s:gsub('det_gammaBar_over_det_gammaHat', '(\\frac{det(\\bar{\\gamma}_{ij})}{det(\\hat{\\gamma}_{ij})})')
			s = s:gsub('det_gammaHat', 'det(\\hat{\\gamma}_{ij})')
			s = s:gsub('det_gammaBar', 'det(\\bar{\\gamma}_{ij})')
			s = s:gsub('gammaBar_uu%.(..)', '\\bar{\\gamma}^{%1}')
			s = s:gsub('gammaBar_ll%.(..)', '\\bar{\\gamma}_{%1}')
			s = s:gsub('gammaBar_UU%.(.)(.)', '\\bar{\\gamma}^{\\hat{%1}\\hat{%2}}')
			s = s:gsub('gammaBar_LL%.(.)(.)', '\\bar{\\gamma}_{\\hat{%1}\\hat{%2}}')
			s = s:gsub('tr_connBar_l%.(.)', '{\\bar{\\Gamma}^*}_{%1 *}')
			s = s:gsub('ABar_ul%.(.)%.(.)', '{\\bar{A}^%1}_%2')
			s = s:gsub('ABar_uu%.(..)', '\\bar{A}^{%1}')
			s = s:gsub('ABar_UL%.(.)%.(.)', '{\\bar{A}^{\\hat{%1}}}_{\\hat{%2}}')
			s = s:gsub('ABar_UU%.(.)(.)', '\\bar{A}^{\\hat{%1}\\hat{%2}}')
			s = s:gsub('ABarSq_ul%.(.)%.(.)', '{(\\bar{A}^2)^%1}_%2')
			s = s:gsub('ABarSq_ll%.(..)', '\\bar{A}_{%1}')
			s = s:gsub('ABarSq_LL%.(.)(.)', '\\bar{A}_{\\hat{%1}\\hat{%2}}')
			s = s:gsub('tr_ABarSq', '({(\\bar{A}^2)^*}_*)')
			s = s:gsub('TF_DBar2_alpha_ll%.(.)(.)', '(\\bar{D}_%1 \\bar{D}_%2 \\alpha)^{TF}')
			s = s:gsub('trBar_DBar2_alpha', '(\\bar{D}^* \\bar{D}_* \\alpha)')
			s = s:gsub('DBar2_alpha_ll%.(.)(.)', '\\bar{D}_%1 \\bar{D}_%2 \\alpha')
			s = s:gsub('DBar_tr_DBar_beta_u%.(.)', '\\bar{D}^{%1} \\bar{D}_* \\beta^*')
			s = s:gsub('tr_DBar_beta', '(\\bar{D}_* \\beta^*)')
			s = s:gsub('partial_phi_l%.(.)', '\\phi_{,%1}')
			s = s:gsub('partial2_phi_ll%.(..)', '\\phi_{,%1}')
			s = s:gsub('TF_RBar_LL%.(.)(.)', '{(\\bar{R})^{TF}}_{\\hat{%1}\\hat{%2}}')
			s = s:gsub('Delta_ull%.(.)%.(..)', '{\\Delta^%1}_{%2}')
			s = s:gsub('Delta_ULL%.(.)%.(.)(.)', '{\\Delta^{\\hat{%1}}}_{\\hat{%2}\\hat{%3}}')
			s = s:gsub('connBar_ULL%.(.)%.(.)(.)', '{\\bar{\\Gamma}^{\\hat{%1}}}_{\\hat{%2}\\hat{%3}}')
			s = s:gsub('Delta_U%.(.)', '\\Delta^{\\hat{%1}}')
			s = s:gsub('tr_gammaBar_DHat2_beta_u%.(.)', '(\\hat{D}^* \\hat{D}_* \\beta^{%1})')
			s = s:gsub('solver%->shift_eta', '\\eta')

			idnum = idnum + 1
			local idname = 'span'..idnum
			s = template([[
<button onclick='renderMathJaxForDiv("<?=idname?>")'>&gt;</button>
<span id='<?=idname?>'><?=s?></span>]], {
				idname = idname,
				s = s,
			})
		end
		return s
	end
	local function convert(...)
		if select('#', ...) == 0 then return end
		local arg = ...
		return convertSingle(arg), select(2, ...)
	end
	function printbr(...)
		local s = table.concat({convert(...)}, '\t')
		print(s..'<br>')
		outfile:write(s..'<br>\n')
		outfile:flush()
	end
--]]	
	var = symmath.var
	sin = symmath.sin
	cos = symmath.cos
	sqrt = symmath.sqrt
	frac = symmath.frac
	Tensor = symmath.Tensor
	Matrix = symmath.Matrix

	--[[
	this is a list of pairs of {from, to}
	the 'to' is usually a variable
	the 'from' is usually partial derivatives, but can be things such as expensive functions you want to minimize calling (sin, cos, etc)
	it is used to replace expressions with variables before compiling
	--]]
	compileRepls = table()	
	
	coords = Tensor.coords()[1].variables

	local replCoords = table()
	for i=1,#coords do
		replCoords[i] = var('x.'..xNames[i])
	end

	function replaceCompileVars(expr)
		for _,repl in ipairs(compileRepls) do
			expr = expr:replace(repl[1], repl[2])
		end
		for i,coord in ipairs(coords) do
			expr = expr:replace(coord, replCoords[i])
		end
		return expr
	end

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
			and x[2].value > 1
			and x[2].value < 100
			then
				return setmetatable(table.rep({x[1]}, x[2].value), symmath.op.mul)
			end
		end)
		expr = replaceCompileVars(expr)
		return symmath.export.C(expr)
	end
	
	-- compile an expression (typically variable), insert it into compileRepls as being associated with a specific expression
	function compileReplVar(repl, ...)
		local v = var(...)
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
		if with == nil then with = env[name] end
		assert(symmath.Expression.is(with), "not an expression")
		return '\treal '..name..' = '..compile(with)..';'
	end
	
	-- generate assignment code for a variable based on its associated expr in the compileRepls table
	function assignRepl(v)
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
	function assignRepls(vars)
		local lines = table()
		for i,v in ipairs(vars) do
			lines:insert(assignRepl(v))
		end
		return lines:concat'\n'
	end

	function assign_real3(name, expr)
		if expr == nil then expr = env[name] end
		local lines = table()
		lines:insert('\treal3 '..name..';')
		for i,xi in ipairs(xNames) do
			assert(expr and expr[i], "failed to find expr for var "..name)
			lines:insert('\t'..name..'.'..xi..' = '..compile(expr[i])..';')
		end
		return lines:concat'\n'
	end

	function assign_sym3(name, expr)
		if expr == nil then expr = env[name] end
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
		if expr == nil then expr = env[name] end
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
		if expr == nil then expr = env[name] end
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
		if expr == nil then expr = env[name] end
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
		if expr == nil then expr = env[name] end
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

	function makevars_real3(variance, name)
		return Tensor(variance, function(i)
			return var(name..'.'..xNames[i], coords)
		end)
	end

	function makevars_real3x3(variance, name)
		return Tensor(variance, function(i,j)
			return var(name..'.'..xNames[i]..'.'..xNames[j], coords)
		end)
	end
	
	function makevars_sym3(variance, name)
		return Tensor(variance, function(i,j)
			return var(name..'.'..sym(i,j), coords)
		end)
	end
	
	function makevars_3sym3(variance, name)
		return Tensor(variance, function(i,j,k) 
			return var(name..'.'..xNames[i]..'.'..sym(j,k), coords) 
		end)
	end


	-- factorDivision doesn't seem to recurse... 
	-- TODO fix this in symmath
	local oldFactorDivision = symmath.factorDivision
	function symmath.factorDivision(expr, ...)
		if Tensor.is(expr) then
			return Tensor(expr.variance, function(...)
				local x = oldFactorDivision(expr[{...}]())
				if symmath.op.add.is(x) then
					for i=1,#x do
						x[i] = x[i]()
					end
				end
				return x
			end)
		else
			return oldFactorDivision(expr)
		end
	end
	symmath.Expression.factorDivision = symmath.factorDivision


	-- from here on out is stuff specific to different functions 

	-- non-coordinate basis
	-- TODO either remove the default metric, or introduce a separate set of indexes for the non-coordinate basis
	-- either one to prevent symmath from automatically raising or lowering
	-- otherwise using e or eu incorrectly will result in extra operations

printbr'gammaHat_ll'
	gammaHat_ll = Tensor.metric().metric
printbr(gammaHat_ll)

printbr'e'
	e = Tensor('_i^I', function(i,j)
		return (i==j and solver.coord.lenExprs[i] or 0)
	end)
printbr(e)
printbr'eu'
	eu = Tensor('^i_I', function(i,j)
		return i==j and (1/e[i][i])() or 0
	end)
printbr(eu)

printbr'det_gammaHat'
	det_gammaHat = Matrix.determinant(gammaHat_ll)
	det_gammaHat_var = var('det_gammaHat', coords)
printbr(det_gammaHat)
printbr'gammaHat_uu'
	gammaHat_uu = Tensor('^ij', table.unpack((Matrix.inverse(gammaHat_ll, nil, nil, nil, 
		-- defer det_gammaHat
		--det_gammaHat_var))))
		-- don't bother.  gammaHat^ij is just an inverse scale matrix
		det_gammaHat))))
printbr(gammaHat_uu)
printbr'partial_gammaHat_lll'
	partial_gammaHat_lll = gammaHat_ll'_ij,k'()
printbr(partial_gammaHat_lll)
printbr'partial2_gammaHat_llll'
	partial2_gammaHat_llll = partial_gammaHat_lll'_ijk,l'()
printbr(partial2_gammaHat_llll)
printbr'connHat_lll'
	connHat_lll = ((partial_gammaHat_lll'_ijk' + partial_gammaHat_lll'_ikj' - partial_gammaHat_lll'_jki')/2)()
printbr(connHat_lll)

-- for spherical this has 4x 1/r's and 2x 1/sin(theta)'s
printbr'connHat_ull'
	connHat_ull = (gammaHat_uu'^im' * connHat_lll'_mjk')()
printbr(connHat_ull)

-- for spherical this has 9x 1/r's and 3x 1/sin(theta)'s
-- soo ... it's not any better to use this.
printbr'connHat_ULL'
	connHat_ULL = (connHat_ull'^i_jk' * e'_i^I' * eu'^j_J' * eu'^k_K')():factorDivision()
	connHat_ULL = Tensor('^I_JK', function(I,J,K)
		return connHat_ULL[I][J][K]():factorDivision()
	end)
	connHat_ULL_vars = makevars_3sym3('^I_JK', 'connHat_ULL')
printbr(connHat_ULL)

printbr'partial_connHat_ulll'
	partial_connHat_ulll = connHat_ull'^i_jk,l'() 
printbr(partial_connHat_ulll)
printbr'partial_det_gammaHat_l'
	partial_det_gammaHat_l = Tensor('_i', function(i) 
		return det_gammaHat:diff(coords[i])() 
	end)
	partial_det_gammaHat_l_vars = Tensor('_i', function(i) 
		return compileReplVar(
			det_gammaHat_var:diff(coords[i])(),
			'partial_det_gammaHat_l.'..xNames[i], 
			coords)
	end)
printbr(partial_det_gammaHat_l)
printbr'partial2_det_gammaHat_ll'
	partial2_det_gammaHat_ll = partial_det_gammaHat_l'_i,j'()
printbr(partial2_det_gammaHat_ll)


	-- state variables

	alpha = var('U->alpha', coords)
	beta_U_vars = makevars_real3('^I', 'U->beta_U')
	epsilon_LL_vars = makevars_sym3('_IJ', 'U->epsilon_LL')
	W = var('U->W', coords)
	K = var('U->K', coords)
	ABar_LL_vars = makevars_sym3('_IJ', 'U->ABar_LL')
	LambdaBar_U_vars = makevars_real3('^I', 'U->LambdaBar_U')
	B_U_vars = makevars_real3('^I', 'U->B_U')
	
	rho = var'U->rho'
	
	pi = var'M_PI'

	-- state variable derivatives (variables created with 'makePartial')

	partial_alpha_l_vars = Tensor('_i', function(i) 
		return compileReplVar(alpha:diff(coords[i]), 'partial_alpha_l['..(i-1)..']', coords) 
	end)
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
	partial_K_l_vars = Tensor('_i', function(i) 
		return compileReplVar(K:diff(coords[i]), 'partial_K_l['..(i-1)..']', coords) 
	end)
	partial_W_l_vars = Tensor('_i', function(i) 
		return compileReplVar(W:diff(coords[i]), 'partial_W_l['..(i-1)..']', coords) 
	end)
	partial2_W_ll_vars = Tensor('_ij', function(i,j) 
		local ij = from3x3to6(i,j)
		return compileReplVar(W:diff(coords[i], coords[j]), 'partial2_W_ll['..(ij-1)..']', coords) 
	end)
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

		-- some other variables

	S = var'S'
	exp_neg4phi = W * W

		-- derivatives rescaled back to non-coordinates:

printbr'partial_alpha_L'
	partial_alpha_L = (eu'^i_I' * partial_alpha_l_vars'_i')():factorDivision()
printbr(partial_alpha_L)
printbr'partial_K_L'
	partial_K_L = (eu'^i_I' * partial_K_l_vars'_i')():factorDivision()
printbr(partial_K_L)

-- this is identity.  don't bother store it.
printbr'gammaHat_LL'
	gammaHat_LL = (eu'^i_I' * (eu'^j_J' * gammaHat_ll'_ij')())()
printbr(gammaHat_LL)

printbr'gammaBar_LL'
	gammaBar_LL = (gammaHat_LL'_IJ' + epsilon_LL_vars'_IJ')()
	gammaBar_LL_vars = makevars_sym3('_IJ', 'gammaBar_LL')
printbr(gammaBar_LL)

printbr'gammaBar_ll'
	gammaBar_ll = (gammaHat_ll'_ij' + epsilon_ll'_ij')()
printbr(gammaBar_ll)

printbr'det_gammaBar'
	det_gammaBar = Matrix.determinant(gammaBar_ll)
	det_gammaBar_var = var('det_gammaBar', coords)
printbr(det_gammaBar)

	-- factor out the coord metric det: r^4 sin(theta)^2 of spherical
	-- but leave the rest as a variable
printbr'det_gammaBar_over_det_gammaHat'
	det_gammaBar_over_det_gammaHat = (det_gammaBar / det_gammaHat)()
	det_gammaBar_over_det_gammaHat_var = var('det_gammaBar_over_det_gammaHat', coords)  
printbr(det_gammaBar_over_det_gammaHat)

printbr'gammaBar_UU'
	gammaBar_UU = Tensor('^ij', table.unpack((Matrix.inverse(gammaBar_LL, nil, nil, nil,
		-- det(gammaBar_IJ) == det(gammaBar_ij)/det(gammaHat_ij)
		det_gammaBar_over_det_gammaHat_var)))) 
	gammaBar_UU_vars = makevars_sym3('^IJ', 'gammaBar_UU')
printbr(gammaBar_UU)

-- [[ defer det_gammaHat in det_gammaBar's definition
printbr'det_gammaBar'
	det_gammaBar = (det_gammaBar_over_det_gammaHat_var * det_gammaHat)() 
printbr(det_gammaBar)
--]]

-- TODO try to never use this in the code.  
-- it can tend to zero, and and gammaBar^ij can be singular
printbr'gammaBar_uu'
--[[ gammaBar^ij := invert(gammaBar_ij)
	gammaBar_uu = Tensor('^ij', table.unpack((Matrix.inverse(gammaBar_ll, nil, nil, nil, 
		-- defer det_gammaBar
		--det_gammaBar_var))))
		-- defer det_gammaBar_over_det_gammaHat (but leave det_gammaHat for further simplification)
		det_gammaBar_over_det_gammaHat_var * det_gammaHat))))
		-- don't
		--det_gammaBar))))
	-- this isn't always completely effective.  sometimes it introduces sin(theta)'s in the denominator
	--gammaBar_uu = (gammaBar_uu / det_gammaHat * det_gammaHat_var)()
--]]
-- [[ defer gammaBar^IJ
-- use gammaBar^ij = e^i_I e^j_J gammaBar^IJ
	gammaBar_uu = (eu'^i_I' * eu'^j_J' * gammaBar_UU_vars'^IJ')():factorDivision()
--]]
	gammaBar_uu_vars = makevars_sym3('^ij', 'gammaBar_uu')
printbr(gammaBar_uu)

printbr'partial_gammaBar_lll'
	partial_gammaBar_lll = gammaBar_ll'_ij,k'()
printbr(partial_gammaBar_lll)
printbr'partial2_gammaBar_llll'
	partial2_gammaBar_llll = partial_gammaBar_lll'_ijk,l'()
printbr(partial2_gammaBar_llll)
-- slow (but necessary to remove the singularity in RBar_ij:
printbr'trBar_partial2_gammaBar_ll'
	trBar_partial2_gammaBar_ll = (gammaBar_uu'^kl' * partial2_gammaBar_llll'_ijkl')():factorDivision()
printbr(trBar_partial2_gammaBar_ll)
printbr'connBar_lll'
	connBar_lll = ((partial_gammaBar_lll'_ijk' + partial_gammaBar_lll'_ikj' - partial_gammaBar_lll'_jki') / 2)()
	connBar_lll_vars = makevars_3sym3('_ijk', 'connBar_lll')
printbr(connBar_lll)

-- TODO this has 1/r's in it.
-- this can be singular.  don't use it.
printbr'connBar_LLL'
	connBar_LLL = (((connBar_lll'_ijk' * eu'^i_I')() * eu'^j_J')() * eu'^k_K')():factorDivision()
	connBar_LLL = Tensor('_IJK', function(I,J,K)
		return connBar_LLL[I][J][K]():factorDivision()
	end)
	connBar_LLL_vars = makevars_3sym3('_IJK', 'connBar_LLL')
printbr(connBar_LLL)

-- TODO this can have 1/r's too ...
-- NOTICE this is rescaled connBar^i_jk, NOT the connection associated with rescaled metric gammaBar_IJ
-- I want to store these variables instead of connBar^i_jk, because connBar^i_jk has 1/r's
-- I'm putting '_' out there as a warning:
--  DO NOT DIFFERENTIATE THIS AND EXPECT THE DERIVATIVE OF THE NON-RESCALED VERSION
printbr'connBar_ULL'
	connBar_ULL = (gammaBar_UU_vars'^IL' * connBar_LLL'_LJK')():factorDivision()
	connBar_ULL_vars = makevars_3sym3('^I_JK', 'connBar_ULL')
printbr(connBar_ULL)

	local function removeBetas(expr)
		expr = expr:replace(tr_DBar_beta_var, 0)
		for i=1,3 do
			expr = expr:replace(beta_U_vars[i], 0)
			if DBar_tr_DBar_beta_u_vars then
				expr = expr:replace(DBar_tr_DBar_beta_u_vars[i], 0)
			end
			if tr_gammaBar_DHat2_beta_u_vars then
				expr = expr:replace(tr_gammaBar_DHat2_beta_u_vars[i], 0)
			end
			for j=1,3 do
				expr = expr:replace(partial_beta_Ul_vars[i][j], 0)
				for k=1,3 do
					expr = expr:replace(partial2_beta_Ull_vars[i][j][k], 0)
				end
			end
		end
		expr = expr()
		return expr
	end


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
printbr'dt_alpha'
	dt_alpha = (-alpha^2 * Q + partial_alpha_l'_i' * beta_u'^i')():factorDivision()
printbr(dt_alpha)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_alpha)():factorDivision())

--[[
dt_alpha:
	alpha
	partial_alpha_l
	beta_u -> beta_U
--]]
	
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
	tr_connBar_l_vars = makevars_real3('_i', 'tr_connBar_l')
printbr(tr_connBar_l)
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
	tr_DBar_beta_var = var('tr_DBar_beta', coords)
printbr(tr_DBar_beta)
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
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_W)():factorDivision())
	
	-------------------------------- K_,t -------------------------------- 
	
	-- TODO just set the metric to gammaBar
	-- or set some indexes to gammaBar and others to gamma? and yet others to non-coord gammaBar?

printbr'ABar_UL'
	ABar_UL = (gammaBar_UU_vars'^IK' * ABar_LL_vars'_KJ')():factorDivision()
	ABar_UL_vars = makevars_real3x3('^I_J', 'ABar_UL')
printbr(ABar_UL)
	tr_ABar = ABar_UL_vars'^I_I'():factorDivision()
printbr'ABar_UU'
	ABar_UU = (ABar_UL_vars'^I_K' * gammaBar_UU_vars'^KJ')()
	ABar_UU_vars = makevars_sym3('^IJ', 'ABar_UU')
printbr(ABar_UU)
printbr'ABarSq_LL'
	ABarSq_LL = (ABar_LL_vars'_IK' * ABar_UL_vars'^K_J')()
	ABarSq_LL_vars = makevars_sym3('_IJ', 'ABarSq_LL')
printbr(ABarSq_LL)
-- [[
printbr('tr_ABarSq')
	tr_ABarSq = (gammaBar_UU_vars'^IJ' * ABarSq_LL_vars'_IJ')():factorDivision()
printbr(tr_ABarSq)
--]]
--[[
printbr('ABarSq_ll')
	ABarSq_ll = (ABar_ll'_ik' * ABar_ul_vars'^k_j')():factorDivision()
	ABarSq_ll_vars = makevars_sym3('_ij', 'ABarSq_ll')
printbr(ABarSq_ll)
printbr('tr_ABarSq')
	tr_ABarSq = (gammaBar_uu_vars'^ij' * ABarSq_ll'_ij')():factorDivision()
printbr(tr_ABarSq)
--]]
--[[
printbr('ABarSq_ul')
	ABarSq_ul = (ABar_ul_vars'^i_k' * ABar_ul_vars'^k_j')():factorDivision()
	ABarSq_ul_vars = makevars_real3x3('^i_j', 'ABarSq_ul')
printbr(ABarSq_ul)
printbr('tr_ABarSq')
	tr_ABarSq = (ABarSq_ul_vars'^k_k')():factorDivision()
printbr(tr_ABarSq)
--]]
	tr_ABarSq_var = var('tr_ABarSq', coords) 

	-- TODO Don't store this as a local var.  Use it symbolically even? 
printbr'DBar2_alpha_ll'
	DBar2_alpha_ll = (
		partial2_alpha_ll'_ij' 
		- eu'^i_I' * eu'^j_J' * connBar_ULL_vars'^K_IJ' * eu'^k_K' * partial_alpha_l'_k'
	)():factorDivision()
printbr(DBar2_alpha_ll)

printbr'DBar2_alpha_LL'
	DBar2_alpha_LL = (
		eu'^i_I' * eu'^j_J' * partial2_alpha_ll'_ij' 
		- connBar_ULL_vars'^K_IJ' * eu'^k_K' * partial_alpha_l'_k'
	)():factorDivision()
	DBar2_alpha_LL_vars = makevars_sym3('_ij', 'DBar2_alpha_LL')
printbr(DBar2_alpha_LL)

printbr('trBar_DBar2_alpha')
	trBar_DBar2_alpha = (gammaBar_UU_vars'^IJ' * DBar2_alpha_LL_vars'_IJ')():factorDivision()
	trBar_DBar2_alpha_var = var('trBar_DBar2_alpha', coords)
printbr(trBar_DBar2_alpha)

printbr('partial_alpha_u')
	partial_alpha_u = (gammaBar_uu_vars'^ij' * partial_alpha_l'_j')():factorDivision()
printbr(partial_alpha_u)

	-- W = exp(-2 phi)
	-- phi = -log(W)/2
	--W_for_phi = W:eq(exp(-2*phi_var))
	--phi = W_for_phi:solve(phi_var)
printbr'phi'
	phi = (-symmath.log(W)/2):factorDivision()
printbr(phi)
printbr'partial_phi_l'
	partial_phi_l = Tensor('_i', function(i) 
		return phi:diff(coords[i])():factorDivision() 
	end)
	partial_phi_l_vars = makevars_real3('_i', 'partial_phi_l')
printbr(partial_phi_l)
printbr'partial_phi_L'
	partial_phi_L = (eu'^i_I' * partial_phi_l'_i')():factorDivision()
printbr(partial_phi_L)

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
	
in fact, gammaBar^ij alpha_,i phi_,j is problematic.
because gammaBar^ij has inverse scale factors.
so, for alpha_,theta, alpha_,phi, W_,theta, W_,phi nonzero ... how do we avoid diverging?
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
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_K)():factorDivision())
	
	-------------------------------- epsilon_ij,t -------------------------------- 

	-- e^i_I e^j_J (gammaBar_ij,k beta^k + gammaBar_ki * beta^k,_j + gammaBar_kj * beta^k,_i)
printbr'gammaBar_times_partial_beta_ll'
	gammaBar_times_partial_beta_ll = (gammaBar_ll'_ik' * partial_beta_ul'^k_j')()
printbr(gammaBar_times_partial_beta_ll)

-- this is the rescaled version of the Lie derivative of gammaBar_ij
-- slow
printbr'Lbeta_gammaBar_LL'
	Lbeta_gammaBar_LL = (
			eu'^i_I' * (
				eu'^j_J' * (
					(beta_u'^k' * partial_gammaBar_lll'_ijk')()
					+ gammaBar_times_partial_beta_ll'_ij'
					+ gammaBar_times_partial_beta_ll'_ji'
				)()
			)()
		)():factorDivision()
printbr(Lbeta_gammaBar_LL)
	LBeta_gammaBar_LL_vars = makevars_sym3('_IJ', 'LBeta_gammaBar_LL')

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

		+ 2/3 gammaBar_ij (alpha ABar^k_k - DBar_k beta^k)
		- 2 alpha ABar_ij
		(notice in 2017 Ruchlin the 'alpha' on the last term isn't there, but it is in the SENR Mathematica notebook and in the 2012 Baumgarte paper)

	Looks like the paper's notation DBar_i beta_j implies DBar_j (gammaBar_ik beta^k) = gammaBar_ki DBar_j beta^k
	...which is not DBar_j ( gamma_ik beta^k ), which was my interpretation of the definition of beta_i := gamma_ij beta^k
	(denoting indexes for coordinate transforms is easier to interpret than denoting entire tensors for coordinate transforms)
	--]]
-- very slow
printbr'dt_epsilon_LL'
	dt_epsilon_LL = (
		Lbeta_gammaBar_LL'_IJ'
		+ frac(2,3) * gammaBar_LL'_IJ' * (
			alpha * tr_ABar
			- tr_DBar_beta
		)() - 2 * alpha * ABar_LL_vars'_IJ'
	)():factorDivision()
printbr(dt_epsilon_LL)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_epsilon_LL)():factorDivision())
--]=]	

	-------------------------------- ABar_ij_,t -------------------------------- 

printbr'partial2_phi_ll'
	partial2_phi_ll = partial_phi_l'_i,j'():factorDivision()
	partial2_phi_ll_vars = makevars_sym3('_ij', 'partial2_phi_ll') 
printbr(partial2_phi_ll)

printbr'Delta_ULL'
	Delta_ULL = (connBar_ULL_vars'^I_JK' - connHat_ULL'^I_JK')():factorDivision()
	Delta_ULL_vars = makevars_3sym3('^I_JK', 'Delta_ULL')
printbr(Delta_ULL)

printbr'Delta_U'
	Delta_U = (Delta_ULL'^I_JK' * gammaBar_UU_vars'^JK')():factorDivision()
	Delta_U_vars = makevars_real3('^I', 'Delta_U') 
printbr(Delta_U)

printbr'Delta_LLL'
	Delta_LLL = (gammaBar_LL'_IL' * Delta_ULL'^L_JK')():factorDivision()
	Delta_LLL_vars = makevars_3sym3('_IJK', 'Delta_LLL')
printbr(Delta_ULL)

printbr'LambdaBar_u'
	LambdaBar_u = (LambdaBar_U_vars'^I' * eu'^i_I')()
printbr(LambdaBar_u)

printbr'DBar2_phi_LL'
	DBar2_phi_LL = (
		eu'^i_I' * eu'^j_J' * partial2_phi_ll_vars'_ij' 
		- connBar_ULL_vars'^K_IJ' * eu'^k_K' * partial_phi_l_vars'_k'
	)():factorDivision()
	DBar2_phi_LL_vars = makevars_sym3('_IJ', 'DBar2_phi_LL')
printbr(DBar2_phi_LL)

-- not used in ABar_ij,t, but used elsewhere
printbr'tr_DBar2_phi'
	tr_DBar2_phi = (DBar2_phi_LL_vars'_IJ' * gammaBar_UU_vars'^IJ')()
printbr(tr_DBar2_phi)


printbr'Lbeta_ABar_LL'
	Lbeta_ABar_LL = (eu'^i_I' * eu'^j_J' * (
		partial_ABar_lll'_ijk' * beta_u'^k'
		+ ABar_ll'_kj' * partial_beta_ul'^k_i'
		+ ABar_ll'_ik' * partial_beta_ul'^k_j'
	))():factorDivision()
printbr(Lbeta_ABar_LL)
	Lbeta_ABar_LL_vars = makevars_sym3('_IJ', 'Lbeta_ABar_LL')

	S_ll_vars = makevars_sym3('_ij', 'U->S_ll')

-- not used by ABar_IJ,t, but used elsewhere
	--2008 Alcubierre eqn 2.8.18
	--2010 Baumgarte, Shapiro eqn 3.10
printbr'RPhi_LL'
	RPhi_LL = (
		-2 * DBar2_phi_LL'_IJ'
		+ 4 * partial_phi_L'_I' * partial_phi_L'_J'
		+ gammaBar_LL'_IJ' * (
			-2 * tr_DBar2_phi
			- 4 * gammaBar_UU'^IJ' * partial_phi_L'_I' * partial_phi_L'_J'
		)
	)():factorDivision()
printbr(RPhi_LL)

printbr'tracelessPart_LL'
	tracelessPart_LL = (
		2 * partial_phi_L'_I' * partial_alpha_L'_J'
		+ 2 * partial_alpha_L'_I' * partial_phi_L'_J'
		- 2 * alpha * DBar2_phi_LL_vars'_IJ'
		+ 4 * alpha * partial_phi_L'_I' * partial_phi_L'_J'
		- 8 * pi * alpha * eu'^i_I' * eu'^j_J' * S_ll_vars'_ij'
	)()
printbr(tracelessPart_LL)
-- NOTICE after assignment of this variable, I'm expecting tracefree() to be called on it
	tracelessPart_LL_vars = makevars_sym3('_IJ', 'tracelessPart_LL') 
	
	TF_DBar2_alpha_LL = makevars_sym3('_IJ', 'TF_DBar2_alpha_LL') 
	TF_RBar_LL = makevars_sym3('_ij', 'TF_RBar_LL') 

--[[
	2017 Ruchlin et al, eqn. 11b
	ABar_ij,t = 
		+ ABar_ij,k beta^k
		+ beta^k_,i ABar_jk
		+ beta^k_,j ABar_ik
		
		+ ABar_ij (alpha K - 2/3 DBar_k beta^k)
		- 2 alpha ABar_ik ABar^k_j
		+ exp(-4 phi) (trace-free terms)_ij
--]]
printbr'dt_ABar_LL'
	dt_ABar_LL = (
		Lbeta_ABar_LL_vars'_IJ'
		+ ABar_LL_vars'_IJ' * (alpha * K - frac(2,3) * tr_DBar_beta_var)()
		- 2 * alpha * ABarSq_LL_vars'_IJ'
		+ exp_neg4phi * (
			tracelessPart_LL_vars'_IJ'
			- TF_DBar2_alpha_LL'_IJ'
			+ alpha * TF_RBar_LL'_IJ'
		)
	)():factorDivision()
printbr(dt_ABar_LL)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_ABar_LL)():factorDivision())


	-------------------------------- LambdaBar^i_,t -------------------------------- 

printbr'Lbeta_LambaBar_U'
	Lbeta_LambaBar_U = (
		e'_i^I' * (
			partial_LambdaBar_ul'^i_j' * beta_u'^j'
			- LambdaBar_u'^j' * partial_beta_ul'^i_j'
		)
	)():factorDivision()
printbr(Lbeta_LambaBar_U)

	tr_gammaBar_DHat2_beta_u_vars = makevars_real3('^i', 'tr_gammaBar_DHat2_beta_u') 
	DBar_tr_DBar_beta_u_vars = makevars_real3('^i', 'DBar_tr_DBar_beta_u')

	--[[
	LambdaBar^i_,t = 
		gammaBar^jk DHat_j DHat_k beta^i
		+ 2/3 Delta^i DBar_j beta^j
		+ 1/3 DBar^i DBar_j beta^j
		- 2 ABar^ij (alpha_,j - 6 phi_,j)
		+ 2 alpha ABar^jk Delta^i_jk		-- NOTICE the paper doesn't have the extra alpha, but the Etienne SENR Mathematica notebook does have it.  Maybe I'm reading an earlier version of the paper?
		- 4/3 alpha gammaBar^ij K_,j
		- 16 pi exp(4 phi) alpha S^i
		
		+ LambdaBar^i_,j beta^j
		- LambdaBar^j beta^i_,j
	
		seems there is no way to avoid 1/r's in this term. 
	--]]
printbr'dt_LambdaBar_U'
	-- [=[ runs until t=0.046
	dt_LambdaBar_U = (
		Lbeta_LambaBar_U'^I' 
		+ 2 * alpha * Delta_ULL_vars'^I_JK' * ABar_UU_vars'^JK'
		+ frac(2,3) * Delta_U_vars'^I' * tr_DBar_beta_var
		- 2 * ABar_UU_vars'^IJ' * (
			partial_alpha_L'_J'
			- 6 * partial_phi_L'_J'
		)
		- frac(4,3) * alpha * gammaBar_UU_vars'^IJ' * partial_K_L'_J'	
		+ e'_i^I' * (
			tr_gammaBar_DHat2_beta_u_vars'^i'
			+ frac(1,3) * DBar_tr_DBar_beta_u_vars'^i'
		)	
		)():factorDivision()
	--]=]
	--[=[ runs until t=0.2something..
	dt_LambdaBar_U = Tensor'^I'
	--]=]
	dt_LambdaBar_U_vars = makevars_real3('^I', 'dt_LambdaBar_U') 
printbr(dt_LambdaBar_U)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_LambdaBar_U)())
	
	-------------------------------- beta^i_,t and B^i_,t -------------------------------- 

	eta_var = var'solver->shift_eta'

	if eqn.useShift == 'GammaDriver' then
		k_var = var'k'
		--Gamma-driver
		--B&S 4.82
		--beta^i_,t = k (connBar^i_,t + eta connBar^i)
printbr'dt_beta_U'
		dt_beta_U = (k_var * dt_LambdaBar_U_vars'^I' + eta_var * LambdaBar_U_vars'^I')()
printbr(dt_beta_U)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_beta_U)():factorDivision())
	
	elseif eqn.useShift == 'HyperbolicGammaDriver' then
printbr'dt_beta_U'
		dt_beta_U = (
			B_U_vars'^I'
			-- if shiftadvect==true from SENR
			+ e'_i^I' * partial_beta_ul'^i_j' * beta_u'^j'
		)():factorDivision()
printbr(dt_beta_U)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_beta_U)():factorDivision())

		--[[
		hyperbolic Gamma driver 
		2017 Ruchlin et al, eqn 14a, 14b
		beta^i_,t = B^i + beta^i_,j beta^j
		B^i_,t = 3/4 LambdaBar^i_,t - eta B^i 

		Etienne's SENR doesn't do the full Lie derivative though.
		just the B^i_,j beta^j term, and adds the full LambdaBar^i_,t to B^i_,t
		--]]
printbr'dt_B_U'
		dt_B_U = (
			frac(3,4) * dt_LambdaBar_U_vars'^I'
			- eta_var * B_U_vars'^I'
			-- if biadvect==true from SENR
			+ e'_i^I' * partial_beta_ul'^i_j' * B_u'^j'
		)():factorDivision()
printbr(dt_B_U)
printbr'...with $\\beta^i = 0$...'
printbr(removeBetas(dt_B_U)():factorDivision())

	end	-- eqn.useShift

	-- [[ do this every time you stop using the env
	setfenv(1, oldEnv)
	--]]

end)	-- time()

	outfile:close()

	return self.env
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template([[

<? local dim = solver.dim ?>

#if 1	
//rescaling, used for bssn finite-difference, but I am tempted to try it with other coordinate systems with singularities
//TODO for the initial conditions do this symbolically instead of numerically

//apply this to lower indexes to convert from coordinate metric to better metric
//apply this to upper indexes to convert from better metric to coordinate metric
real3 real3_rescaleFromCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x / coord_dx0(x),
		.y = v.y / coord_dx1(x),
		.z = v.z / coord_dx2(x),
	};
}
#define real3_rescaleToCoord_U real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_L(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y * coord_dx1(x),
		.z = v.z * coord_dx2(x),
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_L

sym3 sym3_rescaleFromCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleToCoord_UU sym3_rescaleFromCoord_ll

sym3 sym3_rescaleToCoord_LL(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleFromCoord_uu sym3_rescaleToCoord_LL

_3sym3 _3sym3_rescaleFromCoord_lll(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleToCoord_UUU _3sym3_rescaleFromCoord_lll

_3sym3 _3sym3_rescaleToCoord_LLL(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleFromCoord_uuu _3sym3_rescaleToCoord_LLL

sym3sym3 sym3sym3_rescaleFromCoord_lll(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleToCoord_UUUU sym3sym3_rescaleFromCoord_llll

sym3sym3 sym3sym3_rescaleToCoord_LLLL(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleFromCoord_uuuu sym3sym3_rescaleToCoord_LLLL

#else	//debugging -- turning it off

#define real3_rescaleFromCoord_l(a,x) a
#define real3_rescaleToCoord_U(a,x) a
#define real3_rescaleToCoord_L(a,x) a
#define real3_rescaleFromCoord_u(a,x) a
#define sym3_rescaleFromCoord_ll(a,x) a
#define sym3_rescaleToCoord_UU(a,x) a
#define sym3_rescaleToCoord_LL(a,x) a
#define sym3_rescaleFromCoord_uu(a,x) a
#define _3sym3_rescaleFromCoord_lll(a,x) a
#define _3sym3_rescaleToCoord_UUU(a,x) a
#define _3sym3_rescaleToCoord_LLL(a,x) a
#define _3sym3_rescaleFromCoord_uuu(a,x) a
#define sym3sym3_rescaleFromCoord_lll(a,x) a
#define sym3sym3_rescaleToCoord_UUUU(a,x) a
#define sym3sym3_rescaleToCoord_LLLL(a,x) a
#define sym3sym3_rescaleFromCoord_uuuu (a,x) a

#endif




//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero

sym3 calc_gammaHat_ll(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>
	return gammaHat_ll;
}

real calc_det_gammaHat(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaHat'?>
	return det_gammaHat;
}

sym3 calc_gammaHat_uu(real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaHat'?>
<?=assign_sym3'gammaHat_uu'?>
	return gammaHat_uu;
}

//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
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

sym3 calc_gammaBar_LL(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaBar_LL'?>
	return gammaBar_LL;
}

sym3 calc_gammaBar_UU(global const <?=eqn.cons_t?>* U, real3 x) {
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaHat'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>
	return gammaBar_UU;
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

]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_RBar_LL()
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

<?=assign_sym3'gammaBar_uu'?>
	
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
?>	RBar_LL.<?=xij?> = (0.			
			- .5 * trBar_DHat2_gammaBar_ll.<?=xij?>	
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * gammaBar_ll.<?=sym(i,k)?> * DHat_LambdaBar_ul.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_ll.<?=sym(j,k)?> * DHat_LambdaBar_ul.<?=xk?>.<?=xi?>
<?	end
?>	) / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x))
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * Delta_U.<?=xk?> * (
				Delta_LLL.<?=xi?>.<?=sym(k,j)?> 
				+ Delta_LLL.<?=xj?>.<?=sym(k,i)?>
			)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>			+ gammaBar_UU.<?=sym(k,l)?> * (0.
				+ Delta_ULL.<?=xm?>.<?=sym(k,i)?> * Delta_LLL.<?=xj?>.<?=sym(m,l)?>
				+ Delta_ULL.<?=xm?>.<?=sym(k,j)?> * Delta_LLL.<?=xi?>.<?=sym(m,l)?>
				+ Delta_ULL.<?=xm?>.<?=sym(i,k)?> * Delta_LLL.<?=xm?>.<?=sym(j,l)?>
			)
<?			end
		end
	end
?>
	;
<? end ?>
]], self:getEnv())
end

--[[
Should initState provide a metric in cartesian, or in the background metric?
I'll say Cartesian for now, and then transform them using the rescaling.
--]]
function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([=[
<?
-- [[ do this every time you use the env
local env = eqn:getEnv()
local oldEnv = getfenv()
getmetatable(env).__index = oldEnv
setfenv(1, env)
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

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = gammaHat_ll;
	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

	<?=code?>
//rescale from cartesian to spherical
gamma_ll = sym3_rescaleToCoord_LL(gamma_ll, x);

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
<? if false then ?>

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=eqn:makePartial'epsilon_LL'?>

<?=assign'det_gammaHat'?>
<?=assign_3sym3'connHat_ULL'?>

<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>
<?=assign_3sym3'connBar_lll'?>
<?=assign_3sym3'connBar_ULL'?>
<?=assign_3sym3'Delta_ULL'?>	

	U->LambdaBar_U = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);
<? end ?>
}
<?
-- [[ do this every time you stop using the env
setfenv(1, oldEnv)
--]]
?>
]=], {
		eqn = self,
		code = self.initState:initState(self.solver),
	})
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd-sym.cl'

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
	'U RBar_LL xx',
	'U RBar_LL xy',
	'U RBar_LL xz',
	'U RBar_LL yy',
	'U RBar_LL yz',
	'U RBar_LL zz',
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

<?=assign_sym3'gammaBar_UU'?>
<?=assign_3sym3'connBar_ULL'?>

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
<?=assign_real3'partial_phi_l'?>
<?=assign_sym3'partial2_phi_ll'?>

<?=assign_sym3'DBar2_phi_LL'?>
<?=assign'tr_DBar2_phi'?>

	*value = tr_DBar2_phi;
]], self:getEnv())
		},
	
		{
			name = 'partial_phi_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'W'?>
<?=assign_real3'partial_phi_l'?>
	*value_real3 = partial_phi_l;
]], self:getEnv()),
		},

		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'alpha'?>
	*value_real3 = *(real3*)partial_alpha_l;
]], self:getEnv()),
		},

		{
			name = 'DBar2_alpha_ll',
			code = template([[
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>
	
<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_lll'?>
<?=assign_sym3'DBar2_alpha_ll'?> 

	*value = sym3_dot(gammaBar_uu, DBar2_alpha_ll);
]], self:getEnv()),
		},
	
	
		{
			name = 'tracelessPart_LL',
			type = 'sym3',
			code = template([[
	
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>

<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>
<?=assign_3sym3'connBar_lll'?>
<?=assign_sym3'DBar2_alpha_ll'?> 

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
<?=assign_real3'partial_phi_l'?>
<?=assign_sym3'partial2_phi_ll'?>
<?=assign_sym3'DBar2_phi_LL'?>	
<?=assign_sym3'tracelessPart_LL'?>	
	tracelessPart_LL = tracefree(tracelessPart_LL, gammaBar_LL, gammaBar_UU);
	*value_sym3 = tracelessPart_LL; 
]], self:getEnv()),
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
--[=[
	do
		-- [[ do this every time you use the env
		local env = self:getEnv()
		local oldEnv = getfenv()
		getmetatable(env).__index = oldEnv
		setfenv(1, env)
		--]]	

		vars:insert{
			name = 'RBar_LL',
			type = 'sym3',
			code = template([[
<?=eqn:makePartial'LambdaBar_U'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

<?=assign'det_gammaBar_over_det_gammaHat'?>

<?=assign_sym3'gammaBar_LL'?>
<?=assign_sym3'gammaBar_UU'?>

<?=assign'det_gammaHat'?>
<?=assign_sym3'gammaBar_ll'?>
<?=assign_3sym3'connHat_ull'?>
<?=assign_3sym3'connHat_ULL'?>
<?=assign_3sym3'connBar_lll'?>
<?=assign_3sym3'connBar_ULL'?>

<?=assign_real3'LambdaBar_u'?>
<?=assign_real3x3'partial_LambdaBar_ul'?>
<?=assign_real3'Delta_U'?>	

<?=assign_3sym3'Delta_ULL'?>
<?=assign_3sym3'Delta_LLL'?>	

<?=assign_real3'partial_det_gammaHat_l'?>
<?=assign_3sym3x3'partial_connHat_ulll'?>

<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij')?>

	sym3 RBar_LL;
<?=eqn:getCode_RBar_LL()?>
	*value_sym3 = RBar_LL;
]], env),
		}
	
		-- [[ do this every time you stop using the env
		setfenv(1, oldEnv)
		--]]
	end
--]=]		
--[=[		
	vars:insert{
		name = 'RPhi_LL',
		type = 'sym3',
		code = template([[

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
<?=assign_real3'partial_phi_l'?>
<?=assign_sym3'partial2_phi_ll'?>
	
<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign_sym3'gammaBar_uu'?>

<?=assign_3sym3'connBar_lll'?>

<?=assign_sym3'DBar2_phi_LL'?>
<?=assign'tr_DBar2_phi'?>

<?=assign'RPhi_LL'?>
	*value_sym3 = RPhi_ll;
]], self:getEnv()),
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
	do
		-- [[ do this every time you use the env
		local env = self:getEnv()
		local oldEnv = getfenv()
		getmetatable(env).__index = oldEnv
		setfenv(1, env)
		--]]	

		-- symbolically generate ... so far gammaBar^ij	
		-- this fixes it.
		vars:insert{
			name='del gammaBar_ll sym',
			type = 'sym3',
			code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

<?=assign_sym3'gammaBar_LL'?>
<?=assign_sym3'gammaBar_UU'?>

<?=assign_sym3'trBar_partial2_gammaBar_ll'?>
	*value_sym3 = trBar_partial2_gammaBar_ll;
]], env),
		}
		-- [[ do this every time you stop using the env
		setfenv(1, oldEnv)
		--]]
	end
--]=]

--[=[
	vars:insert{
		name = 'tr34 (gamma*dGamma)',
		type = 'real3x3',
		code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaBar_ll'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>
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
]], self:getEnv()),
	}
--]=]

--[=[
	vars:insert{
		name = 'tr14 (Gamma*dgamma)',
		type = 'real3x3',
		code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=eqn:makePartial'epsilon_LL'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign_sym3'gammaBar_UU'?>

<?=assign_3sym3('partial_gammaBar_lll', partial_gammaBar_lll:permute'_kij') 
-- how does :permute() work? forward or backward? 
?>

<?=assign_sym3'gammaBar_ll'?>
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
]], self:getEnv()),
	}
--]=]

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
