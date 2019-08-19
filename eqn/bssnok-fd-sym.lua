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
local BSSNOKFiniteDifferenceEquationBase = require 'eqn.bssnok-fd'
local makestruct = require 'eqn.makestruct'
local common = require 'common'
local time, getTime = table.unpack(require 'util.time')
local makePartials = require 'eqn.makepartial'

local BSSNOKFiniteDifferenceEquation = class(BSSNOKFiniteDifferenceEquationBase)
BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 
BSSNOKFiniteDifferenceEquation.hasEigenCode = true
BSSNOKFiniteDifferenceEquation.hasCalcDTCode = true
BSSNOKFiniteDifferenceEquation.hasFluxFromConsCode = true
BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true

-- not used with finite-difference schemes anyways
BSSNOKFiniteDifferenceEquation.weightFluxByGridVolume = false

-- seems all the hyperbolic formalisms listed in Alcubierre's book use alpha sqrt(gamma^ii) for the speed-of-light wavespeed
-- however the 2017 Ruchlin paper says to use gamma_ij
BSSNOKFiniteDifferenceEquation.cflMethod = '2008 Alcubierre'
--BSSNOKFiniteDifferenceEquation.cflMethod = '2013 Baumgarte et al, eqn 32'
--BSSNOKFiniteDifferenceEquation.cflMethod = '2017 Ruchlin et al, eqn 53'

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
	self.cflMethod = args.cflMethod

	local intVars = table{
		{name='alpha', type='real'},			-- 0:	1: alpha
		{name='W', type='real'},				-- 1:	1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},				-- 2:	1: K = K^i_i
		{name='beta_U', type='real3'},		 	-- 3:	3: beta^i
		{name='B_U', type='real3'},				-- 6:	3: B^i ... only used with HyperbolicGammaDriver
		{name='LambdaBar_U', type='real3'},		-- 9:	3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
		{name='epsilon_LL', type='sym3'},		-- 12:	6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='ABar_LL', type='sym3'},			-- 18:	6: ABar_ij, only 5 dof since ABar^k_k = 0
	}	

	self.consVars = table()
	:append(intVars)
	:append{
		--stress-energy variables:
		{name='rho', type='real'},				--1: n_a n_b T^ab
		{name='S_u', type='real3'},				--3: -gamma^ij n_a T_aj
		{name='S_ll', type='sym3'},				--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},				--1
		{name='M_U', type='real3'},				--3
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

-- NOTICE you have to use this in the scope of a push/pop of the previous fenv
-- ... and you have to set this's mt's __index to the old fenv
function BSSNOKFiniteDifferenceEquation:getEnv()
	if self.env then return self.env end
	
	--[[ 
	look for a cache of the equations
	how to do this ...
	I just realized I don't have a symmath exporter to symmath
	TODO how to save things easily 
	1) make a 'notebook' function env for exportation
	... which records all your assignments as you make them?
	2) recursive export, which tracks what variables are needed and exports them (though this gets tricky when doing lots of replace()'s)
	3) just manually do it
	
	TODO TODO
	I can't just categorize the cache based on the coordinate system.
	There are a few parameters of solver and eqn referenced by the env, like useShift, which will also determine the uniqueness of the env.
	In the end I'll have to make some kind of unique name generation
	option 1) concat all dependent parameters, like I do in the accuracy test files
	option 2) generate a MD5 hash.  this is what the computer science folks like to do, because it looks technical and is obfuscating.  That's exactly why I don't want to do it.
	--]]
	
	local symdir = 'eqn/bssnok-fd-sym/'..self.solver.coord.name
	local function fixpath(s) return require 'ffi'.os == 'Windows' and s:gsub('/', '\\') or s end
	if require 'ffi'.os == 'Windows' then
		os.execute('mkdir '..fixpath(symdir))
	else
		os.execute('mkdir -p '..fixpath(symdir))
	end

	-- next chore: save everything in env to env.lua
	-- this means saving anything that will be assign'd
	-- and it also means saving all the compileRepls (for assign's sake)
	local oldEnv = getfenv()
	local env
	local envCacheFilename = symdir..'/env.lua'
-- [=[
	local envCacheData = file[envCacheFilename]
	if cmdline.bssnUseCache == false then envCacheData = false end
	local gotCache
	if envCacheData then
		print("found cached equations...")
		env = assert(load(envCacheData))()
		self.env = env
		env.eqn = self
		env.solver = self.solver
		setfenv(1, env)
		-- first setup some of the internal env functions, then return the obj
		gotCache = true
	else
--]=]	
		local symmath = require 'symmath'
	
		-- symmath will mess with the metatable of the env table it is passed
		local symenv = {}
		symmath.setup{env=symenv}
		env = setmetatable({}, {
			__index = function(t,k)
				-- I can't directly reference symmath here because symmath.tostring will override tostring
				-- so instead I'll do so through the symenv I pass to symmath.setup()
				local v = symenv[k] if v ~= nil then return v end
				local v = common[k] if v ~= nil then return v end
				local v = oldEnv[k] if v ~= nil then return v end
			end,
		})
		self.env = env
		env.eqn = self
		env.solver = self.solver
		setfenv(1, env)
		--]]

-- [[ mathjax output
		symmath.tostring = symmath.export.MathJax
		symmath.tostring.useCommaDerivative = true

		outfile = assert(io.open(symdir..'/math.html', 'w'))
		outfile:write((tostring(symmath.tostring.header):gsub('tryToFindMathJax.js', '../tryToFindMathJax_andDeferMathRendering.js')))
		local idnum = 1
		local function convertSingle(arg)
			local s = tostring(arg)
			if symmath.Expression.is(arg) then
			-- TODO proper way is to arg:replace() everything
			-- but that is slow
				s = s:gsub('U%->alpha', '\\alpha')
				s = s:gsub('U%->beta_U%.(.)', '\\beta^{\\hat{%1}}')
				s = s:gsub('U%->epsilon_LL%.(.)(.)', '\\epsilon_{\\hat{%1}\\hat{%2}}')
				s = s:gsub('U%->W', 'W')
				s = s:gsub('U%->K', 'K')
				s = s:gsub('U%->ABar_LL%.(.)(.)', '\\bar{A}_{\\hat{%1}\\hat{%2}}')
				s = s:gsub('U%->LambdaBar_U%.(.)', '\\bar{\\Lambda}^{\\hat{%1}}')
				s = s:gsub('U%->B_U%.(.)', 'B^{\\hat{%1}}')
				s = s:gsub('U%->rho', '\\rho')
				s = s:gsub('U%->S_u%.(.)', 'S^{%1}')
				s = s:gsub('U%->S_ll%.(..)', 'S_{%1}')
				s = s:gsub('partial_alpha_l%.(.)', '\\alpha_{,%1}')
				s = s:gsub('partial_beta_Ul%.(.)%.(.)', '{\\beta^{\\hat{%2}}}_{,%1}')
				s = s:gsub('partial_epsilon_LLl%[(.)%]%.(.)(.)', function(k,xi,xj) return '\\epsilon_{\\hat{'..xi..'}\\hat{'..xj..'},'..xNames[k+1]..'}' end)
				s = s:gsub('partial_W_l%.(.)', 'W_{,%1}')
				s = s:gsub('partial_K_l%.(.)', 'K_{,%1}')
				s = s:gsub('partial_ABar_LLl%[(.)%]%.(.)(.)', function(k,xi,xj) return '\\bar{A}_{\\hat{'..xi..'}\\hat{'..xj..'},'..xNames[k+1]..'}' end)
				s = s:gsub('partial_LambdaBar_Ul%.(.)%.(.)', '{\\bar{\\Lambda}^{\\hat{%2}}}_{,%1}')
				s = s:gsub('partial_B_Ul%.(.)%.(.)', '{B^{\\hat{%2}}}_{,%1}')
				s = s:gsub('partial_alpha_l_upwind%.(.)', '(\\alpha^{up})_{,%1}')
				s = s:gsub('partial_beta_Ul_upwind%.(.)%.(.)', '{(\\beta^{up})^{\\hat{%2}}}_{,%1}')
				s = s:gsub('partial_epsilon_LLl_upwind%[(.)%]%.(.)(.)', function(k,xi,xj) return '(\\epsilon^{up})_{\\hat{'..xi..'}\\hat{'..xj..'},'..xNames[k+1]..'}' end)
				s = s:gsub('partial_W_l_upwind%.(.)', '(W^{up})_{,%1}')
				s = s:gsub('partial_K_l_upwind%.(.)', '(K^{up})_{,%1}')
				s = s:gsub('partial_ABar_LLl_upwind%[(.)%]%.(.)(.)', function(k,xi,xj) return '(\\bar{A}^{up})_{\\hat{'..xi..'}\\hat{'..xj..'},'..xNames[k+1]..'}' end)
				s = s:gsub('partial_LambdaBar_Ul_upwind%.(.)%.(.)', '{(\\bar{\\Lambda}^{up})^{\\hat{%2}}}_{,%1}')
				s = s:gsub('partial_B_Ul_upwind%.(.)%.(.)', '{(B^{up})^{\\hat{%2}}}_{,%1}')
				s = s:gsub('partial2_alpha_ll%.(..)', '\\alpha_{,%1}')
				s = s:gsub('partial2_beta_Ul%[(.)%]%.(.)', function(jk,xi) return '{\\beta^{\\hat{'..xi..'}}}_{,'..sym(from6to3x3(jk+1))..'}' end)
				s = s:gsub('partial2_epsilon_LLll%[(.)%]%.(.)(.)', function(kl,xi,xj) return '\\epsilon_{\\hat{'..xi..'}\\hat{'..xj..'},'..sym(from6to3x3(kl+1))..'}' end)
				s = s:gsub('partial2_W_ll%.(..)', 'W_{,%1}')
				s = s:gsub('M_PI', '\\pi')
				s = s:gsub('det_gammaBar_over_det_gammaHat', '(\\frac{det(\\bar{\\gamma}_{ij})}{det(\\hat{\\gamma}_{ij})})')
				s = s:gsub('det_gammaHat', 'det(\\hat{\\gamma}_{ij})')
				s = s:gsub('det_gammaBar', 'det(\\bar{\\gamma}_{ij})')
				s = s:gsub('ABar_ul%.(.)%.(.)', '{\\bar{A}^%1}_%2')
				s = s:gsub('ABar_uu%.(..)', '\\bar{A}^{%1}')
				s = s:gsub('ABarSq_ul%.(.)%.(.)', '{(\\bar{A}^2)^%1}_%2')
				s = s:gsub('ABarSq_ll%.(..)', '\\bar{A}_{%1}')
				s = s:gsub('solver%->shift_eta', '\\eta')
				s = s:gsub('gammaBar_UU%.(.)(.)', '\\bar{\\gamma}^{\\hat{%1}\\hat{%2}}')
				s = s:gsub('connBar_LLL%.(.)%.(.)(.)', '{\\bar{\\Gamma}^{\\hat{%1}}}_{\\hat{%2}\\hat{%3}}')

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
			local converted = {convert(...)}
			local s = table.concat(converted, '\t')
			print(s..'<br>')
			outfile:write(s..'<br>\n')
			outfile:flush()
		end

		recordedAssignments = table()	-- key, value
		getmetatable(env).__newindex = function(t,k,v)
			-- this will only work for the first assignment, mind you
			-- (unless a subsequent assignment is nil)
			-- but otherwise, subsequent assignments don't call __newindex
			-- unless the vlues are read/written to a separate table
			if symmath.Expression.is(v) then
				recordedAssignments:insert(k)
			end
			rawset(t,k,v)
		end

		--[[
		this is a list of pairs of {from, to}
		the 'to' is usually a variable
		the 'from' is usually partial derivatives, but can be things such as expensive functions you want to minimize calling (sin, cos, etc)
		it is used to replace expressions with variables before compiling
		--]]
		compileRepls = table()	
--]]
	end

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
		expr = replaceCompileVars(expr)
		return symmath.export.C(expr)
	end
	
	-- compile an expression (typically variable), insert it into compileRepls as being associated with a specific expression
	function compileReplVar(repl, ...)
		local v = var(...)
		compileRepls:insert{repl, v}
		return v
	end

	-- generates assignment code from specified variables
	function assign(name, with)
		if with == nil then 
			with = env[name] 
			if with == nil then error("couldn't find "..name) end
		end
		if type(with) == 'number' then with = symmath.Constant(with) end
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
assert(env.assignRepls)

	function assign_real3(name, expr)
		if expr == nil then expr = env[name] end
		local lines = table()
		lines:insert('\treal3 '..name..';')
		for i,xi in ipairs(xNames) do
			assert(expr, "failed to find var "..name)
			assert(expr[i], "failed to find var "..name..'['..i..']')
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
			assert(expr, "failed to find var "..name)
			assert(expr[i], "failed to find var "..name..'['..i..']')
			assert(expr[i][j], "failed to find var "..name..'['..i..']['..j..']')
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
				assert(expr, "failed to find var "..name)
				assert(expr[i], "failed to find var "..name..'['..i..']')
				assert(expr[i][j], "failed to find var "..name..'['..i..']['..j..']')
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
				assert(expr, "failed to find var "..name)
				assert(expr[i], "failed to find var "..name..'['..i..']')
				assert(expr[i][j], "failed to find var "..name..'['..i..']['..j..']')
				assert(expr[i][j][k], "failed to find var "..name..'['..i..']['..j..']['..k..']')
				lines:insert('\t'..name..'.'..xi..'.'..xjk..' = '..compile(expr[i][j][k])..';')
			end
		end
		return lines:concat'\n'
	end

	function assign_real3x3x3(name, expr)
		if expr == nil then expr = env[name] end
		local lines = table()
		lines:insert('\treal3x3x3 '..name..';')
		for i,xi in ipairs(xNames) do
			for j,xj in ipairs(xNames) do
				for k,xk in ipairs(xNames) do
					assert(expr, "failed to find var "..name)
					assert(expr[i], "failed to find var "..name..'['..i..']')
					assert(expr[i][j], "failed to find var "..name..'['..i..']['..j..']')
					assert(expr[i][j][k], "failed to find var "..name..'['..i..']['..j..']['..k..']')
					lines:insert('\t'..name..'.'..xi..'.'..xj..'.'..xk..' = '..compile(expr[i][j][k])..';')
				end
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
					assert(expr, "failed to find var "..name)
					assert(expr[i], "failed to find var "..name..'['..i..']')
					assert(expr[i][j], "failed to find var "..name..'['..i..']['..j..']')
					assert(expr[i][j][k], "failed to find var "..name..'['..i..']['..j..']['..k..']')
					assert(expr[i][j][k][l], "failed to find var "..name..'['..i..']['..j..']['..k..']['..l..']')
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
				assert(expr, "failed to find var "..name)
				assert(expr[i], "failed to find var "..name..'['..i..']')
				assert(expr[i][j], "failed to find var "..name..'['..i..']['..j..']')
				assert(expr[i][j][k], "failed to find var "..name..'['..i..']['..j..']['..k..']')
				assert(expr[i][j][k][l], "failed to find var "..name..'['..i..']['..j..']['..k..']['..l..']')
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

	function makevars_real3x3x3(variance, name)
		return Tensor(variance, function(i,j,k)
			return var(name..'.'..xNames[i]..'.'..xNames[j]..'.'..xNames[k], coords)
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
			local x = oldFactorDivision(expr)
			if symmath.op.add.is(x) then
				for i=1,#x do
					x[i] = x[i]()
				end
			end
			return x
		end
	end
	symmath.Expression.factorDivision = symmath.factorDivision
	symmath.Array.factorDivision = symmath.factorDivision
	symmath.Matrix.factorDivision = symmath.factorDivision
	symmath.Tensor.factorDivision = symmath.factorDivision
	symmath.op.add.factorDivision = symmath.factorDivision
	symmath.op.sub.factorDivision = symmath.factorDivision
	symmath.op.mul.factorDivision = symmath.factorDivision
	symmath.op.div.factorDivision = symmath.factorDivision
	symmath.op.pow.factorDivision = symmath.factorDivision
	symmath.op.unm.factorDivision = symmath.factorDivision

	-- done setting up the env, now we can return the cached copy
	if not gotCache then 

time('building symbolic math env', function()
		-- from here on out is stuff specific to different functions 

		-- TODO only generate what we need
		cos_xs = Tensor('_i', function(i) 
			return compileReplVar(cos(coords[i]), 'cos_'..i) 
		end)
		sin_xs = Tensor('_i', function(i) 
			return compileReplVar(sin(coords[i]), 'sin_'..i) 
		end)

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
	printbr(det_gammaHat)
	printbr'gammaHat_uu'
		gammaHat_uu = Tensor('^ij', table.unpack((Matrix.inverse(gammaHat_ll, nil, nil, nil, det_gammaHat))))
	printbr(gammaHat_uu)
	printbr'partial_gammaHat_lll'
		partial_gammaHat_lll = gammaHat_ll'_ij,k'():permute'_ijk'
	printbr(partial_gammaHat_lll)
	printbr'partial2_gammaHat_llll'
		partial2_gammaHat_llll = partial_gammaHat_lll'_ijk,l'():permute'_ijkl'
	printbr(partial2_gammaHat_llll)
	printbr'connHat_lll'
		connHat_lll = ((partial_gammaHat_lll'_ijk' + partial_gammaHat_lll'_ikj' - partial_gammaHat_lll'_jki')/2)()
	printbr(connHat_lll)

	-- for spherical this has 4x 1/r's and 2x 1/sin(theta)'s
	printbr'connHat_ull'
		connHat_ull = (gammaHat_uu'^im' * connHat_lll'_mjk')()
	printbr(connHat_ull)

	printbr'connHat_ULL'
		connHat_ULL = (connHat_ull'^i_jk' * e'_i^I' * eu'^j_J' * eu'^k_K')():factorDivision()
	printbr(connHat_ULL)

	printbr'partial_connHat_ulll'
		partial_connHat_ulll = connHat_ull'^i_jk,l'():permute'^i_jkl'
	printbr(partial_connHat_ulll)
	printbr'partial_det_gammaHat_l'
		partial_det_gammaHat_l = Tensor('_i', function(i) 
			return det_gammaHat:diff(coords[i])() 
		end)
	printbr(partial_det_gammaHat_l)
	printbr'partial2_det_gammaHat_ll'
		partial2_det_gammaHat_ll = partial_det_gammaHat_l'_i,j'()
	printbr(partial2_det_gammaHat_ll)

		-- state variables

		alpha = var('U->alpha', coords)
		beta_U = makevars_real3('^I', 'U->beta_U')
		epsilon_LL = makevars_sym3('_IJ', 'U->epsilon_LL')
		W = var('U->W', coords)
		K = var('U->K', coords)
		ABar_LL = makevars_sym3('_IJ', 'U->ABar_LL')
		LambdaBar_U = makevars_real3('^I', 'U->LambdaBar_U')
		B_U = makevars_real3('^I', 'U->B_U')
		
		rho = var'U->rho'
		S_u = makevars_real3('^i', 'U->S_u')
		S_ll = makevars_sym3('_ij', 'U->S_ll')
		
		H = var('U->H', coords)
		
		pi = var'M_PI'

		-- state variable derivatives (variables created with 'makePartial1')

		partial_alpha_l_vars = Tensor('_i', function(i) 
			return compileReplVar(alpha:diff(coords[i]), 'partial_alpha_l.'..xNames[i], coords) 
		end)
		partial2_alpha_ll_vars = Tensor('_ij', function(i,j)
			return compileReplVar(alpha:diff(coords[i], coords[j]), 'partial2_alpha_ll.'..sym(i,j), coords) 
		end)
		partial_epsilon_LLl_vars = Tensor('_IJk', function(I,J,k)
			return compileReplVar(epsilon_LL[I][J]:diff(coords[k]), 'partial_epsilon_LLl['..(k-1)..'].'..sym(I,J), coords)
		end)
		partial2_epsilon_LLll_vars = Tensor('_IJkl', function(I,J,k,l)
			local kl = from3x3to6(k,l)
			return compileReplVar(epsilon_LL[I][J]:diff(coords[k], coords[l]), 'partial2_epsilon_LLll['..(kl-1)..'].'..sym(I,J), coords)
		end)
		partial_K_l_vars = Tensor('_i', function(i) 
			return compileReplVar(K:diff(coords[i]), 'partial_K_l.'..xNames[i], coords) 
		end)
		partial_W_l_vars = Tensor('_i', function(i) 
			return compileReplVar(W:diff(coords[i]), 'partial_W_l.'..xNames[i], coords) 
		end)
		partial2_W_ll_vars = Tensor('_ij', function(i,j) 
			return compileReplVar(W:diff(coords[i], coords[j]), 'partial2_W_ll.'..sym(i,j), coords) 
		end)
		partial_LambdaBar_Ul_vars = Tensor('^I_j', function(I,j)
			return compileReplVar(LambdaBar_U[I]:diff(coords[j]), 'partial_LambdaBar_Ul.'..xNames[j]..'.'..xNames[I], coords)
		end)
		partial_ABar_LLl_vars = Tensor('_IJk', function(I,J,k)
			return compileReplVar(ABar_LL[I][J]:diff(coords[k]), 'partial_ABar_LLl['..(k-1)..'].'..sym(I,J), coords)
		end)
		partial_beta_Ul_vars = Tensor('^I_j', function(I,j)
			return compileReplVar(beta_U[I]:diff(coords[j]), 'partial_beta_Ul.'..xNames[j]..'.'..xNames[I], coords)
		end)
		partial2_beta_Ull_vars = Tensor('^I_jk', function(I,j,k)
			local jk = from3x3to6(j,k)
			return compileReplVar(beta_U[I]:diff(coords[j], coords[k]), 'partial2_beta_Ull['..(jk-1)..'].'..xNames[I], coords)
		end)
		partial_B_Ul_vars = Tensor('^I_j', function(I,j)
			return compileReplVar(B_U[I]:diff(coords[j]), 'partial_B_Ul.'..xNames[j]..'.'..xNames[I], coords)
		end)

		-- state variables in coordinate form

		epsilon_ll = (epsilon_LL'_IJ' * e'_i^I' * e'_j^J')()
		LambdaBar_u = (LambdaBar_U'^I' * eu'^i_I')()
		ABar_ll = (ABar_LL'_IJ' * e'_i^I' * e'_j^J')()
		beta_u = (beta_U'^I' * eu'^i_I')()
		B_u = (B_U'^I' * eu'^i_I')()

		-- derivatives in coordinate form
			-- without rescaling
		
		partial_alpha_l = Tensor('_i', function(i) return alpha:diff(coords[i]) end)
		partial2_alpha_ll = partial_alpha_l'_i,j'()
		partial_W_l = Tensor('_i', function(i) return W:diff(coords[i]) end)
		partial_K_l = Tensor('_i', function(i) return K:diff(coords[i]) end)

			-- with rescaling (these will be slower and more complex)

		partial_epsilon_lll = epsilon_ll'_ij,k'():permute'_ijk'
		
		partial_ABar_lll = ABar_ll'_ij,k'():permute'_ijk'
		
		partial_LambdaBar_ul = LambdaBar_u'^i_,j'():permute'^i_j'
		partial_beta_ul = beta_u'^i_,j'():permute'^i_j'
		partial2_beta_ull = partial_beta_ul'^i_j,k'():permute'^i_jk'
		partial_B_ul = B_u'^i_,j'():permute'^i_j'

			-- some other variables

	printbr'exp_neg4phi'
		exp_neg4phi = W * W
	printbr(exp_neg4phi)

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
		gammaBar_LL = (gammaHat_LL'_IJ' + epsilon_LL'_IJ')()
	printbr(gammaBar_LL)

	printbr'gammaBar_ll'
		gammaBar_ll = (gammaHat_ll'_ij' + epsilon_ll'_ij')()
	printbr(gammaBar_ll)

	printbr'det_gammaBar'
		det_gammaBar = Matrix.determinant(gammaBar_ll)
	printbr(det_gammaBar)

		-- factor out the coord metric det: r^4 sin(theta)^2 of spherical
		-- but leave the rest as a variable
	printbr'det_gammaBar_over_det_gammaHat'
		det_gammaBar_over_det_gammaHat = (det_gammaBar / det_gammaHat)()
		det_gammaBar_over_det_gammaHat_var = var('det_gammaBar_over_det_gammaHat', coords)  
	printbr(det_gammaBar_over_det_gammaHat)

	printbr'gammaBar_UU'
		gammaBar_UU = Tensor('^IJ', table.unpack((Matrix.inverse(gammaBar_LL, nil, nil, nil,
			-- det(gammaBar_IJ) == det(gammaBar_ij)/det(gammaHat_ij)
			det_gammaBar_over_det_gammaHat_var)))) 
		-- I'm trying to remove all vars, but this one is pretty essential to lowering the symbolic expression sizes
		-- I picked gammaBar^IJ because it is O(1/det gammaBar_IJ), which is near 1
		gammaBar_UU_vars = makevars_sym3('^IJ', 'gammaBar_UU')
	printbr(gammaBar_UU)

	-- [[ defer det_gammaHat in det_gammaBar's definition
	printbr'det_gammaBar'
		det_gammaBar = (det_gammaBar_over_det_gammaHat_var * det_gammaHat)() 
	printbr(det_gammaBar)
	--]]

	printbr'gammaBar_uu'
		gammaBar_uu = (eu'^i_I' * eu'^j_J' * gammaBar_UU_vars'^IJ')():factorDivision()
	printbr(gammaBar_uu)
	printbr'gamma_uu'
		gamma_uu = (exp_neg4phi * gammaBar_uu'^ij')():factorDivision()
	printbr(gamma_uu)

	printbr'partial_gammaBar_lll'
		partial_gammaBar_lll = gammaBar_ll'_ij,k'():permute'_ijk'
	printbr(partial_gammaBar_lll)
	printbr'partial2_gammaBar_llll'
		partial2_gammaBar_llll = partial_gammaBar_lll'_ijk,l'():permute'_ijkl'
	printbr(partial2_gammaBar_llll)
	-- slow (but necessary to remove the singularity in RBar_ij:
	printbr'trBar_partial2_gammaBar_ll'
		trBar_partial2_gammaBar_ll = (gammaBar_uu'^kl' * partial2_gammaBar_llll'_ijkl')():factorDivision()
	printbr(trBar_partial2_gammaBar_ll)
	printbr'connBar_lll'
		connBar_lll = ((partial_gammaBar_lll'_ijk' + partial_gammaBar_lll'_ikj' - partial_gammaBar_lll'_jki') / 2)()
	printbr(connBar_lll)

	printbr'connBar_LLL'
		connBar_LLL = (((connBar_lll'_ijk' * eu'^i_I')() * eu'^j_J')() * eu'^k_K')():factorDivision()
		connBar_LLL_vars = makevars_3sym3('_IJK', 'connBar_LLL')
	printbr(connBar_LLL)

	printbr'connBar_ULL'
		connBar_ULL = (gammaBar_UU_vars'^IL' * connBar_LLL_vars'_LJK')():factorDivision()
	printbr(connBar_ULL)

		local function removeBetas(expr)
			expr = expr:replace(tr_DBar_beta_var, 0)
			for i=1,3 do
				expr = expr:replace(beta_U[i], 0)
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

	printbr'S'
		S = (gamma_uu'^ij' * S_ll'_ij')():factorDivision()
	printbr(S)
	
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
		
		partial_alpha_l_upwind = Tensor('_i', function(i) return var('partial_alpha_l_upwind.'..xNames[i], coords) end)
		
		--d/dt alpha alpha,t - alpha_,i beta^i = -alpha^2 Q
		-- already seeing a 1/r and 1/(r sin(θ)) in the denom...
		-- but it is only next to beta^θ and beta^φ
		--simplification is important for the 1/alpha's at least
	printbr'dt_alpha'
		dt_alpha = (-alpha^2 * Q + partial_alpha_l_upwind'_i' * beta_u'^i')():factorDivision()
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

	-- [[
	printbr'partial_det_gammaBar_over_det_gammaHat_l'
		partial_det_gammaBar_over_det_gammaHat_l = Tensor('_i', function(i)
			return det_gammaBar_over_det_gammaHat:diff(coords[i])():factorDivision()
		end)
	printbr(partial_det_gammaBar_over_det_gammaHat_l)
		partial_det_gammaBar_over_det_gammaHat_l_vars = Tensor('_i', function(i) 
			return compileReplVar(det_gammaBar_over_det_gammaHat_var:diff(coords[i]), 'partial_det_gammaBar_over_det_gammaHat_l.'..xNames[i], coords) 
		end)
	printbr'partial2_det_gammaBar_over_det_gammaHat_ll'
		partial2_det_gammaBar_over_det_gammaHat_ll = partial_det_gammaBar_over_det_gammaHat_l'_i,j'():factorDivision()
	printbr(partial2_det_gammaBar_over_det_gammaHat_ll)
		partial2_det_gammaBar_over_det_gammaHat_ll_vars = Tensor('_ij', function(i,j) 
			return compileReplVar(det_gammaBar_over_det_gammaHat_var:diff(coords[i], coords[j]), 'partial2_det_gammaBar_over_det_gammaHat_ll.'..sym(i,j), coords) 
		end)
	--]]

	-- [[ connBar^k_jk = sqrt(det gammaBar_mn),j
	printbr'tr_connBar_l'	
		tr_connBar_l  = Tensor('_i', function(i)
			return (sqrt(det_gammaBar):diff(coords[i]) / sqrt(det_gammaBar))():factorDivision()
		end)
	printbr(tr_connBar_l)
	--]]
	--[[ connBar^k_jk as a trace 
	tr_connBar_L = connBar_ULL'^K_JK'()
	tr_connBar_l = (e'_i^I' * tr_connBar_L'_I')()
	--]]

	-- connBar^k_ki,j = (1/2 det gammaBar_mn,i / det gammaBar_mn)_,j
	-- so it is a partial2 of det gammaBar_mn, and is symmetric (even though it's partial1 of connBar^j_ij)
	printbr'partial_tr_connBar_ll'
		partial_tr_connBar_ll = tr_connBar_l'_i,j'():factorDivision()
	printbr(partial_tr_connBar_ll)

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
		tr_DBar_beta = (partial_beta_ul'^j_j' + tr_connBar_l'_j' * beta_u'^j')():factorDivision()
	printbr(tr_DBar_beta)
		partial_W_l_upwind = Tensor('_i', function(i) return var('partial_W_l_upwind.'..xNames[i], coords) end)
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
		dt_W = (frac(1,3) * W * (alpha * K - tr_DBar_beta) + beta_u'^k' * partial_W_l_upwind'_k')():factorDivision()
	printbr(dt_W)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_W)():factorDivision())
		
		-------------------------------- K_,t -------------------------------- 
		
		-- TODO just set the metric to gammaBar
		-- or set some indexes to gammaBar and others to gamma? and yet others to non-coord gammaBar?

	printbr'ABar_UL'
		ABar_UL = (gammaBar_UU_vars'^IK' * ABar_LL'_KJ')():factorDivision()
	printbr(ABar_UL)
		tr_ABar = ABar_UL'^I_I'():factorDivision()
	printbr'ABar_UU'
		ABar_UU = (ABar_UL'^I_K' * gammaBar_UU_vars'^KJ')()
	printbr(ABar_UU)
	printbr'ABarSq_LL'
		ABarSq_LL = (ABar_LL'_IK' * ABar_UL'^K_J')()
	printbr(ABarSq_LL)
	printbr('tr_ABarSq')
		tr_ABarSq = (gammaBar_UU_vars'^IJ' * ABarSq_LL'_IJ')():factorDivision()
	printbr(tr_ABarSq)

	printbr'DBar2_alpha_ll'
		DBar2_alpha_ll = (
			partial2_alpha_ll'_ij' 
			- eu'^i_I' * eu'^j_J' * connBar_ULL'^K_IJ' * eu'^k_K' * partial_alpha_l'_k'
		)():factorDivision()
	printbr(DBar2_alpha_ll)

	printbr'DBar2_alpha_LL'
		DBar2_alpha_LL = (
			eu'^i_I' * eu'^j_J' * partial2_alpha_ll'_ij' 
			- connBar_ULL'^K_IJ' * eu'^k_K' * partial_alpha_l'_k'
		)():factorDivision()
	printbr(DBar2_alpha_LL)

	printbr('tr_DBar2_alpha')
		tr_DBar2_alpha = (gammaBar_UU_vars'^IJ' * DBar2_alpha_LL'_IJ')():factorDivision()
	printbr(tr_DBar2_alpha)

	printbr('partial_alpha_u')
		partial_alpha_u = (gammaBar_uu'^ij' * partial_alpha_l'_j')():factorDivision()
	printbr(partial_alpha_u)

		-- W_def = exp(-2 phi_var)
		-- phi_def = -log(W_def)/2
		--W_for_phi = W_def:eq(exp(-2*phi))
		--phi_for_W = W_for_phi:solve(phi)
	printbr'phi'
		phi = (-symmath.log(W)/2):factorDivision()
	printbr(phi)
	printbr'partial_phi_l'
		partial_phi_l = Tensor('_i', function(i) 
			return phi:diff(coords[i])():factorDivision() 
		end)
	printbr(partial_phi_l)
	printbr'partial_phi_L'
		partial_phi_L = (eu'^i_I' * partial_phi_l'_i')():factorDivision()
	printbr(partial_phi_L)

		partial_K_l_upwind = Tensor('_i', function(i) return var('partial_K_l_upwind.'..xNames[i], coords) end)
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
			+ alpha * tr_ABarSq
			- exp_neg4phi * (
				tr_DBar2_alpha
				+ 2 * partial_alpha_u'^i' * partial_phi_l'_i'
			)
			+ partial_K_l_upwind'_i' * beta_u'^i'
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
		partial_epsilon_LLl_upwind = Tensor('_ijk', function(i,j,k) return var('partial_epsilon_LLl_upwind['..(k-1)..'].'..sym(i,j), coords) end)
	printbr'partial_epsilon_lll_upwind'
		partial_epsilon_lll_upwind = (partial_epsilon_LLl_upwind'_IJk' * e'_i^I' * e'_j^J' + epsilon_LL'_IJ' * (e'_i^I' * e'_j^J')'_ij^IJ_,k')():permute'_ijk'
	printbr(partial_epsilon_lll_upwind)
	printbr'partial_gammaBar_lll_upwind'
		partial_gammaBar_lll_upwind = (partial_epsilon_lll_upwind'_ijk' + partial_gammaHat_lll'_ijk')()
	printbr(partial_gammaBar_lll_upwind)
	
	printbr'Lbeta_gammaBar_LL'
		Lbeta_gammaBar_LL = (
				eu'^i_I' * (
					eu'^j_J' * (
						(beta_u'^k' * partial_gammaBar_lll_upwind'_ijk')()
						+ gammaBar_times_partial_beta_ll'_ij'
						+ gammaBar_times_partial_beta_ll'_ji'
					)()
				)()
			)():factorDivision()
	printbr(Lbeta_gammaBar_LL)

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
			)() - 2 * alpha * ABar_LL'_IJ'
		)():factorDivision()
	printbr(dt_epsilon_LL)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_epsilon_LL)():factorDivision())
	--]=]	
		
		-------------------------------- RBar_ij -------------------------------- 

	printbr'Delta_ULL'
		Delta_ULL = (connBar_ULL'^I_JK' - connHat_ULL'^I_JK')():factorDivision()
	printbr(Delta_ULL)

	printbr'Delta_U'
		--Delta_U = (Delta_ULL'^I_JK' * gammaBar_UU_vars'^JK')():factorDivision()
		Delta_U = LambdaBar_U	-- plus mystery C?
	printbr(Delta_U)

	printbr'Delta_LLL'
		Delta_LLL = (gammaBar_LL'_IL' * Delta_ULL'^L_JK')():factorDivision()
	printbr(Delta_LLL)

	printbr'DHat_gammaBar_lll'
		DHat_gammaBar_lll = (
			partial_gammaBar_lll'_ijk'
			- connHat_ull'^l_ki' * gammaBar_ll'_lj'
			- connHat_ull'^l_kj' * gammaBar_ll'_li'
		)():permute'_ijk'
	printbr(DHat_gammaBar_lll)

		--[[
		partial_DHat_gammaBar_llll.ijkl := 
			gammaBar_ij,kl 
			- connHat^m_ki,l gammaBar_mj 
			- connHat^m_kj,l gammaBar_mi
			- connHat^m_ki gammaBar_mj,l
			- connHat^m_kj gammaBar_mi,l
		--]]
		partial_DHat_gammaBar_llll = (
				partial2_gammaBar_llll'_ijkl'
				- partial_connHat_ulll'^m_kil' * gammaBar_ll'_mj'
				- partial_connHat_ulll'^m_kjl' * gammaBar_ll'_mi'
				- connHat_ull'^m_ki' * partial_gammaBar_lll'_mjl'
				- connHat_ull'^m_kj' * partial_gammaBar_lll'_mil'
			)():permute'_ijkl'

		--[[
		DHat2_gammaBar_llll.ijkl := DHat_l DHat_k gammaBar_ij 
		= partial_l DHat_k gammaBar_ij
			- connHat^m_lk DHat_m gammaBar_ij
			- connHat^m_li DHat_k gammaBar_mj
			- connHat^m_lj DHat_k gammaBar_mi
		--]]
	printbr'trBar_DHat2_gammaBar_ll'
		trBar_DHat2_gammaBar_ll = (gammaBar_uu'^kl' * (
				partial2_gammaBar_llll'_ijkl'
				- connHat_ull'^m_lk' * DHat_gammaBar_lll'_ijm'
				- connHat_ull'^m_li' * DHat_gammaBar_lll'_mjk'
				- connHat_ull'^m_lj' * DHat_gammaBar_lll'_mik'
				)
			)():factorDivision():permute'_ij'
	printbr(trBar_DHat2_gammaBar_ll)

	printbr'DHat_LambdaBar_ul'
		DHat_LambdaBar_ul = (partial_LambdaBar_ul'^i_j' + connHat_ull'^i_kj' * LambdaBar_u'^k')():factorDivision():permute'^i_j'
	printbr(DHat_LambdaBar_ul)

	printbr'DHat_LambdaBar_LL'
		DHat_LambdaBar_LL = (gammaBar_ll'_ik' * DHat_LambdaBar_ul'^k_j' * eu'^i_I' * eu'^j_J')():factorDivision():permute'_IJ'
	printbr(DHat_LambdaBar_LL)

		-- this is (gammaBar^kl DHat_k DHat_l gammaBar_ij) e^i_I e^j_J
		-- it is not gammaBar^kl DHat_k DHat_l (gammaBar_ij e^i_I e^j_J)
	printbr'trBar_DHat2_gammaBar_LL'
		trBar_DHat2_gammaBar_LL = (trBar_DHat2_gammaBar_ll'_ij' * eu'^i_I' * eu'^j_J')():factorDivision():permute'_IJ'
	printbr(trBar_DHat2_gammaBar_LL)

		-- Delta^K Delta_IJK
	printbr'DeltaSq1_LL'
		DeltaSq1_LL = (Delta_U'^K' * Delta_LLL'_IJK')():factorDivision():permute'_IJ'
	printbr(DeltaSq1_LL)

	printbr'Delta_UUL'
		Delta_UUL = (gammaBar_UU_vars'^JL' * Delta_ULL'^I_LK')():factorDivision():permute'^IJ_K'
	printbr(Delta_UUL)

		-- Delta^ML_I * Delta_JML
	printbr'DeltaSq2_LL'
		DeltaSq2_LL = (Delta_UUL'^ML_I' * Delta_LLL'_JML')():factorDivision():permute'_IJ'
	printbr(DeltaSq2_LL)

		-- Delta^ML_I Delta_MJL
	printbr'DeltaSq3_LL'
		DeltaSq3_LL = (Delta_UUL'^ML_I' * Delta_LLL'_MJL')():factorDivision():permute'_IJ'
	printbr(DeltaSq3_LL)

	printbr'RBar_LL'
		--[[
		2017 Ruchlin eqn 12
		RBar_ij = 
			-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
			+ 1/2 gammaBar_ki DHat_j LambdaBar^k
			+ 1/2 gammaBar_kj DHat_i LambdaBar^k
			+ 1/2 Delta^k Delta_ikj
			+ 1/2 Delta^k Delta_jki
			+ Delta^m_ki Delta_jm^k
			+ Delta^m_kj Delta_im^k
			+ Delta^m_ik Delta_mj^k
		--]]	
		RBar_LL = (
			frac(1,2) * (
				- trBar_DHat2_gammaBar_LL'_IJ'
				+ DHat_LambdaBar_LL'_IJ'
				+ DHat_LambdaBar_LL'_JI'
			)
			+ frac(1,2) * (
				DeltaSq1_LL'_IJ'
				+ DeltaSq1_LL'_JI'
			)
			+ DeltaSq2_LL'_IJ'
			+ DeltaSq2_LL'_JI'
			+ DeltaSq3_LL'_IJ'
		)():factorDivision()
	printbr(RBar_LL)

		-------------------------------- ABar_ij_,t -------------------------------- 

	printbr'partial2_phi_ll'
		partial2_phi_ll = partial_phi_l'_i,j'():factorDivision()
	printbr(partial2_phi_ll)

	printbr'LambdaBar_u'
		LambdaBar_u = (LambdaBar_U'^I' * eu'^i_I')()
	printbr(LambdaBar_u)

	printbr'DBar2_phi_LL'
		DBar2_phi_LL = (
			eu'^i_I' * eu'^j_J' * partial2_phi_ll'_ij' 
			- connBar_ULL'^K_IJ' * eu'^k_K' * partial_phi_l'_k'
		)():factorDivision()
	printbr(DBar2_phi_LL)

	-- not used in ABar_ij,t, but used elsewhere
	printbr'tr_DBar2_phi'
		tr_DBar2_phi = (DBar2_phi_LL'_IJ' * gammaBar_UU_vars'^IJ')()
	printbr(tr_DBar2_phi)

	printbr'partial_ABar_lll_upwind'
		partial_ABar_LLl_upwind = Tensor('_ijk', function(i,j,k) return var('partial_ABar_LLl_upwind['..(k-1)..'].'..sym(i,j), coords) end)
		partial_ABar_lll_upwind = (partial_ABar_LLl_upwind'_IJk' * e'_i^I' * e'_j^J' + ABar_LL'_IJ' * (e'_i^I' * e'_j^J')'_ij^IJ_,k')():permute'_ijk'
	printbr(partial_ABar_lll_upwind)

	printbr'Lbeta_ABar_LL'
		Lbeta_ABar_LL = (eu'^i_I' * eu'^j_J' * (
			partial_ABar_lll_upwind'_ijk' * beta_u'^k'
			+ ABar_ll'_kj' * partial_beta_ul'^k_i'
			+ ABar_ll'_ik' * partial_beta_ul'^k_j'
		))():factorDivision()
	printbr(Lbeta_ABar_LL)

		partial_phi_sq = (gammaBar_UU_vars'^KL' * partial_phi_L'_K' * partial_phi_L'_L')():factorDivision()

	-- not used by ABar_IJ,t, but used elsewhere
		--2008 Alcubierre eqn 2.8.18
		--2010 Baumgarte, Shapiro eqn 3.10
	printbr'RPhi_LL'
		RPhi_LL = (
			-2 * DBar2_phi_LL'_IJ'
			+ 4 * partial_phi_L'_I' * partial_phi_L'_J'
			+ gammaBar_LL'_IJ' * (
				-2 * tr_DBar2_phi
				- 4 * partial_phi_sq
			)
		)():factorDivision()
	printbr(RPhi_LL)

		--[[
		2017 Ruchlin et al eqn 11b
		traceless portion of ...
			-2 alpha DBar_i DBar_j phi 
			+ 4 alpha DBar_i phi DBar_j phi 
			+ 2 DBar_i phi DBar_j alpha 
			+ 2 DBar_i alpha DBar_j phi 
			- DBar_i DBar_j alpha 
			+ alpha RBar_ij 
			- 8 pi alpha S_ij
		
		SENR does trace-free on a few terms individually:
		*) -DBar_i DBar_j alpha
		*) alpha RBar_ij
		*) everything else
		--]]
	printbr'tracelessPart_LL'
		tracelessPart_LL = (
			2 * partial_phi_L'_I' * partial_alpha_L'_J'
			+ 2 * partial_alpha_L'_I' * partial_phi_L'_J'
			- 2 * alpha * DBar2_phi_LL'_IJ'
			+ 4 * alpha * partial_phi_L'_I' * partial_phi_L'_J'
			- 8 * pi * alpha * eu'^i_I' * eu'^j_J' * S_ll'_ij'
		)()
	printbr(tracelessPart_LL)
	printbr'TF_tracelessPart_LL'
		TF_tracelessPart_LL = (tracelessPart_LL'_IJ' - frac(1,3) * gammaBar_LL'_IJ' * (gammaBar_UU_vars'^KL' * tracelessPart_LL'_KL'))():factorDivision()
	printbr(TF_tracelessPart_LL)
	
	printbr'TF_DBar2_alpha_LL'
		TF_DBar2_alpha_LL = (DBar2_alpha_LL'_IJ' - frac(1,3) * gammaBar_LL'_IJ' * (gammaBar_UU_vars'^KL' * DBar2_alpha_LL'_KL'))():factorDivision()
	printbr(TF_DBar2_alpha_LL)

	printbr'tr_RBar'
		tr_RBar = (gammaBar_UU_vars'^KL' * RBar_LL'_KL')():factorDivision()
	printbr(tr_RBar)
	printbr'TF_RBar_LL'
		TF_RBar_LL = (RBar_LL'_IJ' - frac(1,3) * gammaBar_LL'_IJ' * tr_RBar)():factorDivision()
	printbr(TF_RBar_LL)

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
			Lbeta_ABar_LL'_IJ'
			+ ABar_LL'_IJ' * (alpha * K - frac(2,3) * tr_DBar_beta)()
			- 2 * alpha * ABarSq_LL'_IJ'
			+ exp_neg4phi * (
				TF_tracelessPart_LL'_IJ'
				- TF_DBar2_alpha_LL'_IJ'
				+ alpha * TF_RBar_LL'_IJ'
			)
		)():factorDivision()
	printbr(dt_ABar_LL)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_ABar_LL)():factorDivision())

		-------------------------------- LambdaBar^i_,t -------------------------------- 

	printbr'partial_LambdaBar_ul_upwind'
		partial_LambdaBar_Ul_upwind = Tensor('^I_j', function(i,j) return var('partial_LambdaBar_Ul_upwind.'..xNames[j]..'.'..xNames[i], coords) end)
		partial_LambdaBar_ul_upwind = (partial_LambdaBar_Ul_upwind'^I_j' * eu'^i_I' + LambdaBar_U'^I' * eu'^i_I,j')():permute'^i_j'
	printbr(partial_LambdaBar_ul_upwind)

	printbr'Lbeta_LambaBar_U'
		Lbeta_LambaBar_U = (
			e'_i^I' * (
				partial_LambdaBar_ul_upwind'^i_j' * beta_u'^j'
				- LambdaBar_u'^j' * partial_beta_ul'^i_j'
			)
		)():factorDivision()
	printbr(Lbeta_LambaBar_U)

	--[[
	DHat2_beta_ull.i.j.k = DHat_k DHat_j beta^i
	= DHat_k (beta^i_,j + connHat^i_lj beta^l)
	= (beta^i_,j + connHat^i_lj beta^l)_,k
		+ connHat^i_mk (beta^m_,j + connHat^m_lj beta^l)
		- connHat^m_jk (beta^i_,m + connHat^i_lm beta^l)
	= beta^i_,jk 
		+ connHat^i_lj,k beta^l
		+ connHat^i_lj beta^l_,k
		+ connHat^i_lk beta^l_,j
		- connHat^l_jk beta^i_,l
		+ connHat^i_mk connHat^m_lj beta^l
		- connHat^m_jk connHat^i_lm beta^l
	--]]
	printbr'DHat_beta_ul'
		DHat_beta_ul = (beta_u'^i_,j' + connHat_ull'^i_kj' * beta_u'^k')():permute'^i_j'
	printbr(DHat_beta_ul)
	printbr'DHat2_beta_ull'
		DHat2_beta_ull = (DHat_beta_ul'^i_j,k' 
			+ connHat_ull'^i_kl' * DHat_beta_ul'^l_j'
			- connHat_ull'^l_kj' * DHat_beta_ul'^i_l'
		)():factorDivision():permute'^i_jk'
	printbr(DHat2_beta_ull)
	printbr'trBar_DHat2_beta_u'
		trBar_DHat2_beta_u = (DHat2_beta_ull'^i_jk' * gammaBar_uu'^jk')():factorDivision()
	printbr(trBar_DHat2_beta_u)
	
	--[[
	tr12_partial2_beta_l.i := beta^j_,ji
	beta^i_,jk is a sym3 of 3's ... so I don't have that struct yet ... 
	 What name could I use? sym3x3?  how about real3s3x3 where we have 's' for symmetric and 'x' for cartesian product.
	 Then sym3 turns into real3s3 and _3sym3 turns into real3x3s3.
	--]]
	printbr'tr12_partial2_beta_l'
		tr12_partial2_beta_l = partial2_beta_ull'^j_ij'():factorDivision()
	printbr(tr12_partial2_beta_l)

	--[[
	DBar_tr_DBar_beta_l.i = DBar_i DBar_j beta^j
	= DBar_i (beta^j_,j + tr_connBar_l beta^l)
	= (beta^j_,j + tr_connBar_l beta^l)_,i
		+ connBar^j_mi (beta^m_,j + connBar^m_lj beta^l)
		- connBar^m_ji (beta^j_,m + connBar^j_lm beta^l)
	= beta^j_,ji
		+ tr_connBar_l,i beta^l
		+ tr_connBar_l beta^l_,i

	...using determinant instead...

	= beta^j_,ji 
		+ 1/2 det_gammaBar_,l / det_gammaBar beta^l_,i
		+ (1/2 det_gammaBar_,l / det_gammaBar),i beta^l
	= beta^j_,ji + 1/2 (
		+ det_gammaBar_,l beta^l_,i / det_gammaBar
		+ det_gammaBar_,il beta^l / det_gammaBar
		- det_gammaBar_,l det_gammaBar_,i beta^l / det_gammaBar^2
	)

	Etienne's SENR uses this alternative formulation: 
	= beta^j_,ji
		+ 1/2 (
			gammaBar_,ij beta^j
			- gammaBar_,i gammaBar_,j beta^j / gammaBar
			+ gammaBar_,j beta^j_,i
		) / gammaBar
	--]]
	printbr'DBar_tr_DBar_beta_l'
		DBar_tr_DBar_beta_l = (tr12_partial2_beta_l'_i'
			+ (tr_connBar_l'_j' * partial_beta_ul'^j_i')()
			+ (partial_tr_connBar_ll'_ij' * beta_u'^j')()
		)():factorDivision()
	printbr(DBar_tr_DBar_beta_l)
	
	printbr'DBar_tr_DBar_beta_u'
		DBar_tr_DBar_beta_u = (gammaBar_uu'^ij' * DBar_tr_DBar_beta_l'_j')():factorDivision()
	printbr(DBar_tr_DBar_beta_u)
		
		--[[
		LambdaBar^i_,t = 
			+ LambdaBar^i_,j beta^j
			- LambdaBar^j beta^i_,j
			
			+ 2/3 Delta^i DBar_j beta^j
			+ 1/3 DBar^i DBar_j beta^j
			- 4/3 alpha gammaBar^ij K_,j
			+ gammaBar^jk DHat_j DHat_k beta^i
			
			- 2 ABar^ij (alpha_,j - 6 alpha phi_,j)
			+ 2 alpha ABar^jk Delta^i_jk
			- 16 pi exp(4 phi) alpha S^i
			
			seems there is no way to avoid 1/r's in this term. 
		--]]
	printbr'dt_LambdaBar_U'
		dt_LambdaBar_U = (
			Lbeta_LambaBar_U'^I' 
			+ 2 * alpha * Delta_ULL'^I_JK' * ABar_UU'^JK'
			+ frac(2,3) * Delta_U'^I' * tr_DBar_beta
			- 2 * ABar_UU'^IJ' * (
				partial_alpha_L'_J'
				- 6 * alpha * partial_phi_L'_J'
			)
			- frac(4,3) * alpha * gammaBar_UU_vars'^IJ' * partial_K_L'_J'	
			+ e'_i^I' * (
				trBar_DHat2_beta_u'^i'
				+ frac(1,3) * DBar_tr_DBar_beta_u'^i'
			)	
			)():factorDivision()
	printbr(dt_LambdaBar_U)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_LambdaBar_U)():factorDivision())
		
		-------------------------------- beta^i_,t and B^i_,t -------------------------------- 
	
		eta_var = var'solver->shift_eta'

		--if eqn.useShift == 'GammaDriver' then
			
			k_var = var'k'
			--Gamma-driver
			--B&S 4.82
			--beta^i_,t = k (connBar^i_,t + eta connBar^i)
	printbr'dt_beta_U_GammaDriver'
			dt_beta_U_GammaDriver = (k_var * dt_LambdaBar_U'^I' + eta_var * LambdaBar_U'^I')()
	printbr(dt_beta_U_GammaDriver)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_beta_U_GammaDriver)():factorDivision())
		
		--elseif eqn.useShift == 'HyperbolicGammaDriver' then
	
	printbr'partial_beta_ul_upwind'
		partial_beta_Ul_upwind = Tensor('^I_j', function(i,j) return var('partial_beta_Ul_upwind.'..xNames[j]..'.'..xNames[i], coords) end)
		partial_beta_ul_upwind = (partial_beta_Ul_upwind'^I_j' * eu'^i_I' + beta_U'^I' * eu'^i_I,j')():permute'^i_j'
	printbr(partial_beta_ul_upwind)
	
	printbr'dt_beta_U_HyperbolicGammaDriver'
			dt_beta_U_HyperbolicGammaDriver = (
				B_U'^I'
				-- if shiftadvect==true from SENR
				+ e'_i^I' * partial_beta_ul_upwind'^i_j' * beta_u'^j'
			)():factorDivision()
	printbr(dt_beta_U_HyperbolicGammaDriver)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_beta_U_HyperbolicGammaDriver)():factorDivision())
	
	printbr'partial_B_ul_upwind'
		partial_B_Ul_upwind = Tensor('^I_j', function(i,j) return var('partial_B_Ul_upwind.'..xNames[j]..'.'..xNames[i], coords) end)
		partial_B_ul_upwind = (partial_B_Ul_upwind'^I_j' * eu'^i_I' + B_U'^I' * eu'^i_I,j')():permute'^i_j'
	printbr(partial_B_ul_upwind)

			--[[
			hyperbolic Gamma driver 
			2017 Ruchlin et al, eqn 14a, 14b
			beta^i_,t = B^i + beta^i_,j beta^j
			B^i_,t = 3/4 LambdaBar^i_,t - eta B^i 

			Etienne's SENR doesn't do the full Lie derivative though.
			just the B^i_,j beta^j term, and adds the full LambdaBar^i_,t to B^i_,t
			--]]
	printbr'dt_B_U_HyperbolicGammaDriver'
			dt_B_U_HyperbolicGammaDriver = (
				frac(3,4) * dt_LambdaBar_U'^I'
				- eta_var * B_U'^I'
				-- if biadvect==true from SENR
				+ e'_i^I' * partial_B_ul_upwind'^i_j' * beta_u'^j'
			)():factorDivision()
	printbr(dt_B_U_HyperbolicGammaDriver)
	printbr'...with $\\beta^i = 0$...'
	printbr(removeBetas(dt_B_U_HyperbolicGammaDriver)():factorDivision())

		--end	-- eqn.useShift
		
		-------------------------------- H & M^i -------------------------------- 
	
		-- 2017 Ruchlin et al, eqn 46
		-- H = 2/3 K^2 - ABar^ij ABar_ij + exp(-4 phi) (RBar - 8 DBar^i phi DBar_i phi - 8 gammaBar^ij DBar_i DBar_j phi)
	printbr('H_def')
		H_def = (
			frac(2,3) * K^2 
			- tr_ABarSq
			+ exp_neg4phi * (
				tr_RBar
				- 8 * partial_phi_sq
				- 8 * tr_DBar2_phi
			)
			- 16 * pi * rho
		)():factorDivision()
	printbr(H_def)
	
	printbr'M_U_def'
		--[[
		
		TODO use this one
		2013 Baumgarte et al, eqn 14
		M^i = exp(-4 phi) (
			1/sqrt(det gammaBar) DHat_j (sqrt(det gammaBar) ABar^ij)
			+ 6 ABar^ij phi_,j
			- 2/3 gammaBar^ij K_,j
			ABar^jk Delta^i_jk
		) - 8 pi S^i
		= exp(-4 phi) (
			connBar^k_jk ABar^ij
			+ DHat_j ABar_kl gammaBar^ik gammaBar^lj
			- (DHat_j gammaBar_mn) (
				gammaBar^im gammaBar^nk ABar_kl gammaBar^lj
				+ gammaBar^lm gammaBar^nj ABar_kl gammaBar^ki
			)
			+ 6 ABar^ij phi_,j
			- 2/3 gammaBar^ij K_,j
			ABar^jk Delta^i_jk
		) - 8 pi S^i
		= exp(-4 phi) (
			+ gammaBar^jk DHat_j ABar_kl gammaBar^li
			
			+ ABar^ij (connBar^k_jk + 6 phi_,j)
			- 2/3 gammaBar^ij K_,j
			
			+ Delta^i_jk ABar^jk
			
			- DHat_j gammaBar_kl (
				gammaBar^ik ABar^lj
				+ gammaBar^kj ABar^li
			)
		) - 8 pi S^i

		--]]
		partial_gammaBar_LLL = (partial_gammaBar_lll'_ijk' * eu'^i_I' * eu'^j_J' * eu'^k_K')():permute'_IJK'
		partial_ABar_LLL = (partial_ABar_lll'_ijk' * eu'^i_I' * eu'^j_J' * eu'^k_K')():permute'_IJK'
		exp_neg6phi = W * W * W
		M_U_def = (
			- frac(2,3) * gammaBar_UU_vars'^IJ' * partial_K_L'_J'
			+ exp_neg6phi * (
				- 8 * pi * S_u'^i' * e'_i^I'
				+ 6 * ABar_UU'^IJ' * partial_phi_L'_J'
				- connBar_ULL'^I_JK' * ABar_UU'^JK'
				- ABar_UU'^KJ' * gammaBar_UU_vars'^LI' * partial_gammaBar_LLL'_KLJ'
				- ABar_UU'^KI' * gammaBar_UU_vars'^LJ' * partial_gammaBar_LLL'_KLJ'
				+ gammaBar_UU_vars'^IK' * gammaBar_UU_vars'^LJ' * partial_ABar_LLL'_KLJ'
			))():factorDivision()
	printbr(M_U_def)

		-------------------------------- init cond stuff -------------------------------- 
		-- the rest if init cond stuff goes in the init state object
		-- this just holds whatever derived values

	printbr'LambdaBar0_U'
		LambdaBar0_U = (Delta_ULL'^I_JK' * gammaBar_UU_vars'^JK')():factorDivision()
	printbr(LambdaBar0_U)


		outfile:close()


		-- TODO add in the solver.coord.coords somewhere
		-- TODO and 'symmath=require'symmath'' and setfenv so we can capture all assignments 

		local function save()
			local lines = table()
			
			lines:insert('--[[ recording these variables:')
			lines:append(recordedAssignments)
			lines:insert('--]]')

			lines:insert[[
-- begin prefix code
local table = require 'ext.table'
local symmath = require 'symmath'
local common = require 'common'
local oldEnv = getfenv()
local symenv = {}
symmath.setup{env=symenv}
local env = setmetatable({}, {
	__index = function(t,k)
		local v = symenv[k] if v ~= nil then return v end
		local v = common[k] if v ~= nil then return v end
		local v = oldEnv[k] if v ~= nil then return v end
	end,
})
setfenv(1, env)
-- end prefix code
]]

			for _,name in ipairs(recordedAssignments) do
				local expr = env[name]
				assert(expr, "couldn't find expression for "..name)
				lines:insert(name..' = '..symmath.export.SymMath(expr))
			end

			lines:insert'compileRepls = table{'
			for _,p in ipairs(compileRepls) do
				local repl, v = table.unpack(p)
				-- TODO don't print 'v', instead print v's name in env[]
				-- buuuuttt what about when it is a sub-assignment?  
				--  like var'cos_1' is equal to cos_xs[1]
				-- I suppose as long as variable equality is by name, no need to worry for now
				lines:insert('\t{'..symmath.export.SymMath(repl)..', '..symmath.export.SymMath(v)..'},')
			end
			lines:insert'}'
			lines:insert[[
-- begin prefix code
return env
-- end prefix code
]]

			file[envCacheFilename] = lines:concat'\n'
		end
		save()
end)	-- time()
	end

	return self.env
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template(file['eqn/bssnok-fd-sym.cl'], setmetatable({getCommonCode=true}, {__index=self:getEnv()}))
end

--[[
Should initState provide a metric in cartesian, or in the background metric?
I'll say Cartesian for now, and then transform them using the rescaling.
--]]
function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	-- do this first to initialize the expression fields
	local env = self:getEnv()
	local initState = self.initState
	
	-- look for symmath expressions instead of code
	-- also skip the initDerivs finite difference 
	if initState.initAnalytical then
		
		local symmath = require 'symmath'
		local Tensor = symmath.Tensor
		local Matrix = symmath.Matrix
		local frac = symmath.frac
		local sqrt = symmath.sqrt

		local eu = env.eu
		local gamma0_ll = initState.gamma0_ll
		local K0_ll = initState.K0_ll

		local gammaHat_ll = Tensor.metric().metric
		local det_gammaHat = Matrix.determinant(gammaHat_ll)
		local det_gammaBar = det_gammaHat 	-- TODO make this constraint a function
		local det_gamma0 = Matrix.determinant(gamma0_ll)
		local exp_neg4phi = ((det_gammaBar / det_gamma0)^frac(1,3))()
		local W0 = sqrt(exp_neg4phi)()
		local gammaBar0_ll = (exp_neg4phi * gamma0_ll'_ij')()
		local epsilon0_ll = (gammaBar0_ll'_ij' - gammaHat_ll'_ij')()
		local epsilon0_LL = (epsilon0_ll'_ij' * eu'^i_I' * eu'^j_J')()
		local gamma0_uu = Tensor('^ij', table.unpack((Matrix.inverse(gamma0_ll, nil, nil, nil, det_gamma0))))
		local K0 = (gamma0_uu'^ij' * K0_ll'_ij')()
		local A0_ll = (K0_ll'_ij' - frac(1,3) * gamma0_ll'_ij' * (gamma0_uu'^kl' * K0_ll'_kl'))()
		local ABar0_ll = (exp_neg4phi * A0_ll'_ij')()
		local ABar0_LL = (ABar0_ll'_ij' * eu'^i_I' * eu'^j_J')()

		env.alpha0 = initState.alpha0
		env.beta0_U = (initState.beta0_u'^i' * eu'^i_I')():factorDivision()
		-- TODO what to initialize B^i to?
		env.B0_U = Tensor('^I', 0,0,0)
		env.W0 = W0
		env.epsilon0_LL = epsilon0_LL
		env.K0 = K0
		env.ABar0_LL = ABar0_LL
-- TODO print these into the math.html file
for k in ('alpha0 W0 K0 beta0_U B0_U epsilon_LL ABar0_LL'):gmatch'%S+' do
	print(k, env[k])
end

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

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

<?=assign'alpha0'?>
	U->alpha = alpha0;
<?=assign_real3'beta0_U'?>	
	U->beta_U = beta0_U;
<?=assign_real3'B0_U'?>	
	U->B_U = B0_U;
<?=assign'W0'?>
	U->W = W0;
<?=assign_sym3'epsilon0_LL'?>
	U->epsilon_LL = epsilon0_LL;
<?=assign'K0'?>
	U->K = K0;
<?=assign_sym3'ABar0_LL'?>
	U->ABar_LL = ABar0_LL;
<?=assign_real3'LambdaBar0_U'?>
	U->LambdaBar_U = LambdaBar0_U;

//TODO initialization of these ...
//how about an initial call to constrainU?	
	U->rho = 0.;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	
	U->H = 0.;
	U->M_U = real3_zero;
}
]], 	setmetatable({
			initState = initState,
		}, {
			__index = env,
		}))
	end
	
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

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign_sym3'gammaHat_ll'?>

	real alpha = 1.;
	real3 beta_u = real3_zero;
	real3 B_u = real3_zero;

	//initState will assume it is providing a metric in Cartesian
	sym3 gamma_ll = sym3_ident;
	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

	<?=code?>
	
	//rescale from cartesian to spherical
	gamma_ll = sym3_rescaleToCoord_LL(gamma_ll, x);

	U->alpha = alpha;
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);
	U->B_U = real3_rescaleFromCoord_u(B_u, x);

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
	U->M_U = real3_zero;
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

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=eqn:makePartial1'epsilon_LL'?>

<?=assign'det_gammaBar_over_det_gammaHat'?>

<?=assign_3sym3'Delta_ULL'?>
<?=assign_real3'LambdaBar0_U'?>
	U->LambdaBar_U = LambdaBar0_U;
}
]=], setmetatable({
		code = initState:initState(self.solver),
	}, {
		__index = env,
	}))
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd-sym.cl'
function BSSNOKFiniteDifferenceEquation:getSolverCode()
	return template(file[self.solverCodeFile], self:getEnv())
end

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
--	'U M_U mag',
	'U M_U x',
	'U M_U y',
	'U M_U z',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_LL tr weighted',
--]=]	

--[[ debugging derivatives
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
	'deriv M_U x',
	'deriv M_U y',
	'deriv M_U z',
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
-- [=[
	'U tr_ABarSq',
	'U tr_DBar2_alpha',
--]=]
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)
	local env = self:getEnv()

	vars:append{
		{name='gamma_ll', type='sym3', code = [[	*value_sym3 = calc_gamma_ll(U, x);]]},
		{name='gamma_uu', type='sym3', code=[[	*value_sym3 = calc_gamma_uu(U, x);]]},
		{name='gammaHat_ll', type='sym3', code=[[	*value_sym3 = calc_gammaHat_ll(x);]]},
		{name='gammaHat_uu', type='sym3', code=[[	*value_sym3 = calc_gammaHat_uu(x);]]},
		{name='gammaBar_ll', type='sym3', code=[[	*value_sym3 = calc_gammaBar_ll(U, x);]]},
		{name='gammaBar_uu', type='sym3', code=[[	*value_sym3 = calc_gammaBar_uu(U, x);]]},
		{name='gammaBar_LL', type='sym3', code=[[	*value_sym3 = calc_gammaBar_LL(U, x);]]},
		{name='gammaBar_UU', type='sym3', code=[[ *value_sym3 = calc_gammaBar_UU(U, x);]]},
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
		{name='det gamma based on phi', code=template([[
<?=assign'exp_neg4phi'?>
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = det_gamma;
]], env)},
		
		{name='S', code='*value = sym3_dot(U->S_ll, calc_gamma_uu(U, x));'},
		
		{
			name='volume', 
			code=template([[
	//|g| = exp(12 phi) |g_grid|
<?=assign'exp_neg4phi'?>
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = U->alpha * det_gamma;
]], env),
		},
		{name='f', code='*value = calc_f(U->alpha);'},
		{name='df/dalpha', code='*value = calc_dalpha_f(U->alpha);'},
	
		{
			name = 'ABarSq_LL',
			type = 'sym3',
			code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>

<?=assign_sym3'ABarSq_LL'?>
	*value_sym3 = ABarSq_LL;
]], env),
		},
		
--[=[	
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>
	*value = tr_DBar2_phi;
]], env)
		},
--]=]	
		
		{
			name = 'partial_phi_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial1'W'?>
<?=assign_real3'partial_phi_l'?>
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
<?=assign_sym3'DBar2_alpha_ll'?> 
	*value = sym3_dot(gammaBar_uu, DBar2_alpha_ll);
]], env),
		},
	
	
		{
			name = 'TF_tracelessPart_LL',
			type = 'sym3',
			code = template([[
<?=eqn:makePartial1'alpha'?>
<?=eqn:makePartial2'alpha'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>
<?=assign_sym3'TF_tracelessPart_LL'?>	
	*value_sym3 = TF_tracelessPart_LL; 
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
		code= template([[
<?=eqn:makePartial1'alpha'?>

	real _1_alpha = 1. / U->alpha;

<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

<?=eqn:makePartial1'beta_U'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
<?=assign_sym3'gammaBar_ll'?>
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
<?=eqn:makePartial1'W'?>
<?=assign_3sym3('partial_gamma_lll', partial_gamma_lll:permute'_kij'?>

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
]],	env), 
		type = 'real3',
	}
--]=]
-- [=[
	vars:insert{
		name = 'RBar_LL',
		type = 'sym3',
		code = template([[
<?=eqn:makePartial1'LambdaBar_U'?>
<?=eqn:makePartial1'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>

<?=assign'det_gammaBar_over_det_gammaHat'?>

<?=assign_sym3'connBar_LLL'?>
<?=assign_sym3'RBar_LL'?>
	*value_sym3 = RBar_LL;
]], env),
	}
--]=]
--[=[
	vars:insert{
		name = 'RPhi_LL',
		type = 'sym3',
		code = template([[
<?=eqn:makePartial1'W'?>
<?=eqn:makePartial2'W'?>
<?=assign'det_gammaBar_over_det_gammaHat'?>
<?=assign'det_gammaBar'?>
<?=assign'RPhi_LL'?>
	*value_sym3 = RPhi_LL;
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
	do
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
<?=eqn:makePartial1'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

	*value_sym3 = trBar_partial2_gammaBar_ll;
]], env),
		}
	end
--]=]

-- [=[
	vars:insert{
		name='tr_ABarSq',
		code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'tr_ABarSq'?>
	*value = tr_ABarSq;
]], env),
	}
--]=]

-- [=[
	vars:insert{
		name='tr_DBar2_alpha',
		code = template([[
<?=assignRepls(cos_xs)?>
<?=assignRepls(sin_xs)?>
<?=assign'tr_DBar2_alpha'?>
	*value = tr_DBar2_alpha;
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
