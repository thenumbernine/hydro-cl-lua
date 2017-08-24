local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local template = require 'template'
local vec3 = require 'vec.vec3'
local clnumber = require 'cl.obj.number'
local InitCond = require 'init.init'

--symmath.tostring = require 'symmath.tostring.SingleLine'		

local xNames = table{'x', 'y', 'z'}
local symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

local from3x3to6_table = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6},}
local function from3x3to6(i,j) return from3x3to6_table[i][j] end

local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
local function from6to3x3(i) return table.unpack(from6to3x3_table[i]) end

local function getTemplateEnv(solver)
	return {
		solver = solver,
		xNames = xNames,
		symNames = symNames,
		from3x3to6 = from3x3to6,
		from6to3x3 = from6to3x3,
		clnumber = clnumber,
	}
end

local function symMat33Det(xx, xy, xz, yy, yz, zz)
	return xx * yy * zz
		+ xy * yz * xz
		+ xz * xy * yz
		- xz * yy * xz
		- yz * yz * xx
		- zz * xy * xy
end

local function symMat33Inv(xx, xy, xz, yy, yz, zz)
	local d = symMat33Det(xx, xy, xz, yy, yz, zz)
	return (yy * zz - yz * yz) / d,		-- xx
			(xz * yz - xy * zz) / d,	-- xy
			(xy * yz - xz * yy) / d,	-- xz
			(xx * zz - xz * xz) / d,	-- yy
			(xz * xy - xx * yz) / d,	-- yz
			(xx * yy - xy * xy) / d		-- zz
end

local function compileC(expr, name, vars)
	assert(type(expr) == 'table', "expected table, found "..type(expr))
	if symmath.Expression.is(expr) then 
		expr = expr()
	
		local range = require 'ext.range'
		-- replace pow(x,2) with x*x
		expr = expr:map(function(x)
			if symmath.op.pow.is(x) 
			and symmath.Constant.is(x[2])
			then
				local value = assert(x[2].value)
				if value > 0 and value == math.floor(value) then
					if value == 1 then
						return x[1]
					else
						return symmath.op.mul(range(value):map(function() 
							return symmath.clone(x[1])
						end):unpack())
					end
				end
			end
		end)		
		
		print('compiling '..name..':')
		print(expr)
		local code = expr:compile(vars, 'C'), name
	
		-- ugh...
		code = code:gsub('sqrt%(', '(real)sqrt((real)')
		
		print(code)
		return code
	end
	return table.map(expr, function(v,k) return compileC(v,k,vars) end), name
end

local NumRelInitCond = class(InitCond)

function NumRelInitCond:getCodePrefix(solver)
	-- looks like all num rel solvers might need this
	-- maybe I should put it in InitCond?
	
	local alphaVar = symmath.var'alpha'
	-- TODO each eqn must be stated as a guiVar of the solver
	-- make a parent class or something for all num rel eqns
	local fGuiVar = solver.eqn.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	local codes = table()
	codes.f = compileC(f, 'f', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})

	return codes:map(function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'		
end

-- I'm working on a unified initial condition code for 1D and 3D NR problems:
--[[
args:
	solver = the solver
	getCodes = function that accepts a table of the expressions alpha, beta, gamma, K
				and returns which expressions are needed by this particular equation

	args = {x,y,z} spatial basis variables
	alpha = lapse expression
	beta = shift expressions 
	gamma = 3-metric
	K = extrinsic curvature

	symmetric 3x3 matrices (particularly gamma & K) are stored as {xx,xy,xz,yy,yz,zz}

	(still working on)
	density, pressure = lua function

	preCompile = function to call on expressions before compile
--]]
local function initNumRel(args)	
	local vars = assert(args.vars)
	
	local exprs = table{
		alpha = assert(args.alpha),
		beta = {table.unpack(args.beta or {0,0,0})},	-- optional
		gamma = {table.unpack(args.gamma)},
		K = {table.unpack(args.K)}
	}
	assert(#exprs.beta == 3, "didn't find 3 beta vars")
	assert(#exprs.gamma == 6, "didn't find 6 gamma vars")
	assert(#exprs.K == 6, "didn't find 6 K vars")

	local function toExpr(expr, name)
		if type(expr) == 'number' then expr = symmath.Constant(expr) end
		if type(expr) == 'table' then
			if not expr.isa then
				expr = table.map(expr, toExpr)
			end
			-- simplify?
			if symmath.Expression.is(expr) then
				expr = expr()
			end
		end
		return expr, name
	end
	print('converting everything to expressions...')
	exprs = table.map(exprs, toExpr)
	print('...done converting everything to expressions')
	
	-- ADM
	
	-- assuming this is one of the ADM 1D solvers
	print('compiling expressions...')
	local codes = table.map(
		assert(args.getCodes(exprs, vars, args), "getCodes needs to return something"),
		function(v,k) 
			if args.preCompile then
				v = args.preCompile(v)
			end
			return compileC(v,k,vars) 
		end)
	print('...done compiling expressions')
	
	-- here's the lapse conditions 

	local alphaVar = symmath.var'alpha'
	-- TODO each eqn must be stated as a guiVar of the solver
	-- make a parent class or something for all num rel eqns
	local fGuiVar = args.solver.eqn.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	codes.f = compileC(f, 'f', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})

	return table.map(codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

local function buildFCCode(solver, diff)
	local alphaVar = symmath.var'alpha'
	local fGuiVar = solver.eqn.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local fCCode = compileC(f, 'f', {alphaVar})

	return fCCode
end

return table{
	-- from 1997 Alcubierre "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
	{
		name = 'gaussian perturbation',
		init = function(self, solver)
			local size = solver.maxs[1] - solver.mins[1]
			
			-- 1997 Alcubierre uses amplitude of 5 on a grid size of 300 with dx=1
			--local H = 5 / 300 * size
			--local sigma = 10 / 300 * size
			-- 1998 Bona et al use amplitude of .05 on a grid size of 100 with dx=.01
			--local H = .05 * size
			--local sigma = math.sqrt(.05) * size
			-- messing with it on my own
			local H = .01 * size
			local sigma = .1 * size

			-- TODO if only OpenCL allowed something like uniforms ...
			solver.eqn:addGuiVar{name = 'H', value = H}
			solver.eqn:addGuiVar{name = 'sigma', value = sigma}
		end,
		getCodePrefix = function(self, solver, getCodes)

			-- here's the coordinates

			local xs = xNames:map(function(x) return symmath.var(x) end)
			symmath.Tensor.coords{{variables=xs}}
			local x,y,z = xs:unpack()

			-- here's the metric

			local H = solver.eqn.guiVars.H.value
			local sigma = solver.eqn.guiVars.sigma.value
		
			local xc = .5 * (solver.mins[1] + solver.maxs[1])
			local yc = .5 * (solver.mins[2] + solver.maxs[2])
			local zc = .5 * (solver.mins[3] + solver.maxs[3])
			local xcs = {xc,yc,zc}
			
			local alpha = 1
			--[=[
			alpha = 1/2 * (
				(1 + symmath.sqrt(1 + kappa))
				* symmath.sqrt((1-h:diff(x))/(1+h:diff(x)))
				- 
				kappa / (1 + symmath.sqrt(1 + kappa))
				* symmath.sqrt((1+h:diff(x))/(1-h:diff(x)))
			),
			--]=]

			--local s = (x - xc)^2 + (y - yc)^2 + (z - zc)^2
			local s = xs:sub(1,solver.dim):map(function(x,i)
				return (x - xcs[i])^2
			end):sum()

			-- to be fair, I don't think Alcubierre ever does this wave in more than 1D 
			local h = H * symmath.exp(-s / sigma^2)
			h = h()
			local dh = symmath.Tensor('_i', function(i) return h:diff(xs[i])() end)
			local d2h = symmath.Tensor('_ij', function(i,j) return h:diff(xs[i],xs[j])() end)

			local delta = symmath.Tensor('_ij', function(i,j) return i == j and 1 or 0 end)
			local gamma = (delta'_ij' - dh'_i' * dh'_j')()
			local K = (-d2h'_ij' / (1 - (dh'_k' * dh'_k')() )^.5 )()
			
			print('...done deriving and compiling.')
			
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = symmath.clone(alpha),
				gamma = symNames:map(function(xij,ij) return gamma[{from6to3x3(ij)}] end),
				K = symNames:map(function(xij,ij) return K[{from6to3x3(ij)}] end),
			}
		end,
	},
	-- from 2012 Alic et. al. "Conformal and covariant formulations of the Z4 system with constraint-violation damping"
	{
		name = 'plane gauge wave',
		getCodePrefix = function(self, solver, getCodes)
			local xs = xNames:map(function(x) return symmath.var(x) end)
			symmath.Tensor.coords{{variables=xs}}
			local x,y,z = xs:unpack()

			local A = .1
			local L = 1
			local h = 1 - A * symmath.sin((2 * math.pi / L) * x)
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = symmath.sqrt(h),
				gamma = {h,0,0,1,0,1},
				K = {0,0,0,0,0,0},
			}
		end,
	},
	{
		name = 'Alcubierre warp bubble',
		getCodePrefix = function(self, solver, getCodes)
			-- [[ safe values
			local R = .5		-- warp bubble radius
			local sigma = 8	-- warp bubble thickness
			local speed = .1	-- warp bubble speed
			--]]
			--[[
			local R = .5
			local sigma = 8
			local speed = 1
			--]]

			local xs = xNames:map(function(x) return symmath.var(x) end)
			symmath.Tensor.coords{{variables=xs}}
			local x,y,z = xs:unpack()
			
			local x_s = 0 -- speed * t
			local v_s = speed -- x_s:diff(t)()
			local r_s = ((x - x_s)^2 + y^2 + z^2)^.5
			local f = (symmath.tanh(sigma * (r_s + R)) - symmath.tanh(sigma * (r_s - R))) / (2 * symmath.tanh(sigma * R))
		
			local betaUx = -v_s * f

			local alpha = 1

			local K_xx = betaUx:diff(x)() / alpha
			local K_xy = betaUx:diff(y)() / (2 * alpha)
			local K_xz = betaUx:diff(z)() / (2 * alpha)

			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				
				vars = xs,
				
				alpha = alpha,
				
				--[[
				interesting note:
				Alcubierre warp bubble drive depends on a beta parameter
				(the local metric is completely flat)
				however so does the toy 1+1 relativity require a beta parameter to produce the stated extrinsic curvature
				yet the toy 1+1 sample set beta=0
				--]]
				beta = {betaUx, 0, 0},

				gamma = {1, 0, 0, 1, 0, 1},		-- identity
				K = {K_xx, K_xy, K_xz, 0, 0, 0},
			}
		end,
	},
	{
		name = 'Schwarzschild black hole',
		getCodePrefix = function(self, solver, getCodes)
			
			local R = .002	-- Schwarzschild radius
			
--[=[ using symbolic calculations
			local x,y,z = symmath.vars('x','y','z')
			local r = symmath.sqrt(x^2 + y^2 + z^2)
		
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = {x,y,z},
				-- 4D metric ADM components:
				alpha = symmath.sqrt(1 - R/r),
				beta = {0,0,0},
				gamma = {
					x^2/((r/R-1)*r^2) + 1,	-- xx
					x*y/((r/R-1)*r^2),		-- xy
					x*z/((r/R-1)*r^2),		-- xz
					y^2/((r/R-1)*r^2) + 1,	-- yy
					y*z/((r/R-1)*r^2),		-- yz
					z^2/((r/R-1)*r^2) + 1,	-- zz
				},
				K = {0,0,0,0,0,0},
			}
--]=]
--[=[ just pass it the cl code.
			local fCCode = buildFCCode(solver)
			
			return template([[
#define calc_f(alpha)			(<?=fCCode?>)
#define rSq(x,y,z) 				(x*x + y*y + z*z)
#define r(x,y,z) 				sqrt(rSq(x,y,z))
#define calc_alpha(x,y,z) 		sqrt(1. - <?=R?>/r(x,y,z))
<? 
for i,xi in ipairs(xNames) do
?>#define calc_beta_<?=xi?>(x,y,z)		0.
<?
end
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>#define calc_gamma_<?=xij?>(x,y,z) 	(<?=xi?>*<?=xj?>/((r(x,y,z) / <?=R?> - 1.) * rSq(x,y,z))<?= i==j and ' + 1.' or ''?>)
<?
end
for ij,xij in ipairs(symNames) do
?>#define calc_K_<?=xij?>(x,y,z)		0.
<?
end
?>]], 		table(getTemplateEnv(solver), {
				fCCode = fCCode:match'{ return (.*); }',
				R = clnumber(R),
			}))
--]=]
		end,
	},
	{	-- Baumgarte & Shapiro, table 2.1, isotropic coordinates
		name = 'black hole - isotropic',
		getCodePrefix = function(self, solver, getCodes)
			local R = .001	-- Schwarzschild radius
		
			local fCCode = buildFCCode(solver)

			solver:setBoundaryMethods'fixed'
--[=[ runs a lot faster than passing r=sqrt(x^2+y^2+z^2) to compile
--	but the equations are too complex -- they tend to crash upon compile
			local xs = xNames:map(function(x) return symmath.var(x) end)
			symmath.Tensor.coords{{variables=xs}}
			local x,y,z = xs:unpack()
			local r = symmath.var('r', xs)
			local r_from_xyz = r:eq(symmath.sqrt(x^2+y^2+z^2))
			local psi = (1 + 1/4 * R / r)
			local alpha = (1 - 1/4 * R / r) / psi
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = alpha,
				beta = {0,0,0},
				gamma = {psi^4, 0, 0, psi^4, 0, psi^4},
				K = {0,0,0,0,0,0},
				preCompile = function(x)
					return (x
						:replace(r:diff(x), x/r)
						:replace(r:diff(y), y/r)
						:replace(r:diff(z), z/r)
						:subst(r_from_xyz))()
				end,
			}
--]=]
-- [=[ the downside to explicitly defining stuff is that it has to be done for all aux vars that a solver might use
			return template([[
inline real calc_f(real alpha) { return <?=fCCode?>; }
inline real bssn_psi(real r) { return 1. + .25 * <?=R?> / r; }	//sum of mass over radius
inline real bssn_psi2(real r) { return 1. - .25 * <?=R?> / r; }

inline real calc_alpha(real x, real y, real z) { 
	real r = sqrt(x*x + y*y + z*z);
	return bssn_psi2(r) / bssn_psi(r);
}

<?
for i,xi in ipairs(xNames) do
?>#define calc_beta_<?=xi?>(x,y,z)		0.
<?
end
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
	if i==j then
?>

inline real calc_gamma_<?=xij?>(real x, real y, real z) {
	real r = sqrt(x*x + y*y + z*z);
	real psi = bssn_psi(r);
	real psiSq = psi*psi;
	return psiSq*psiSq;
}

<?	else
?>#define calc_gamma_<?=xij?>(x,y,z)	0.
<?	end
end
for ij,xij in ipairs(symNames) do
?>#define calc_K_<?=xij?>(x,y,z)		0.
<?
end
?>

//aux vars used by adm3d
<? for i,xi in ipairs(xNames) do
?>

inline real d<?=xi?>_bssn_psi(real x, real y, real z) {
	real r = sqrt(x*x + y*y + z*z);
	real rSq = r * r;
	real rCubed = rSq * r;
	return -.25 * <?=R?> * <?=xi?> / rCubed;
}

inline real d<?=xi?>_bssn_psi2(real x, real y, real z) {
	return -d<?=xi?>_bssn_psi(x,y,z);
}

inline real calc_d<?=xi?>_alpha(real x, real y, real z) {
	real r = sqrt(x*x + y*y + z*z);
	real psi = bssn_psi(r);
	real psi2 = bssn_psi2(r);
	real dpsi = d<?=xi?>_bssn_psi(x,y,z);
	real dpsi2 = d<?=xi?>_bssn_psi2(x,y,z);
	return dpsi2 / psi - psi2 * dpsi2 / (psi * psi);
}

inline real calc_a_<?=xi?>(real x, real y, real z) {
	return calc_d<?=xi?>_alpha(x,y,z) / calc_alpha(x,y,z);
}

<? end
for ij,xij in ipairs(symNames) do
	for k,xk in ipairs(xNames) do
?>
inline real calc_d_<?=xk?><?=xij?>(real x, real y, real z) {
	real r = sqrt(x*x + y*y + z*z);
	real psi = bssn_psi(r);
	real dpsi = d<?=xk?>_bssn_psi(x,y,z);
	return 4. * psi * psi * psi * dpsi;
}
<?	end
end
?>]], 		table(getTemplateEnv(solver), {
				fCCode = fCCode:match'{ return (.*); }',
				R = clnumber(R),
			}))
--]=]
		end,
	},
	{
		name = 'spinning black hole - isotropic',
		getCodePrefix = function(self, solver, getCodes)
			local R = .01
		end,
	},
	{
		name = 'binary black holes - isotropic',
		getCodePrefix = function(self, solver, getCodes)
			local fCCode = buildFCCode(solver)
		
			-- [=[
			return template([[
#define calc_f(alpha)			(<?=fCCode?>)
#define rSq(x,y,z)	 			((x)*(x) + (y)*(y) + (z)*(z))
#define r1(x,y,z) 				sqrt(rSq(x - .25, y, z))
#define r2(x,y,z) 				sqrt(rSq(x + .25, y, z))
#define sq(x)					((x)*(x))
#define bssn_psi(x,y,z)			(1. + .25 * <?=R1?> / r1(x, y, z) + .25 * <?=R2?> / r2(x, y, z))
#define calc_alpha(x,y,z) 		((1. - .25 * <?=R1?> / r1(x, y, z) - .25 * <?=R2?> / r2(x, y, z)) / bssn_psi(x, y, z))
<? 
for i,xi in ipairs(xNames) do
?>#define calc_beta_<?=xi?>(x,y,z)		0.
<?
end
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
	if i==j then
?>#define calc_gamma_<?=xij?>(x,y,z) 	sq(sq(bssn_psi(x,y,z)))
<?	else
?>#define calc_gamma_<?=xij?>(x,y,z)	0.
<?	end
end
for ij,xij in ipairs(symNames) do
?>#define calc_K_<?=xij?>(x,y,z)		0.
<?
end
?>]], 		table(getTemplateEnv(solver), {
				fCCode = fCCode:match'{ return (.*); }',
				R1 = clnumber(.01),
				R2 = clnumber(.01),
			}))	
			--]=]
		end,
	},
	{
		name = 'stellar model',
		getCodePrefix = function(self, solver, getCodes)
			-- hmm, the symbolic stuff seemed to be working so much better in the last project ...
			-- now i'm replacing it with macros for speed's sake ...
			do
				local bodyMass = .001
				local bodyRadius = .1 
				
				local fCCode = buildFCCode(solver)
				
--[[
K_ij = -alpha Gamma^t_ij
...and for Schwarzschild, Gamma^t_ij = 0
...so K_ij = 0

a_i = ln(alpha)_,i = alpha,i / alpha

alpha = (1 - 2 minMass / r)^(1/2)
minMass = bodyMass * min((r/bodyRadius)^3, 1)
d/dr minMass = r >= bodyRadius and 0 or (3 r^2 bodyMass / bodyRadius^3)
d/dr (minMass / r) = (d/dr minMass) / r - minMass / r^2
	= bodyMass ((r >= bodyRadius and 0 or (3 r / bodyRadius^3)) - min((r/bodyRadius)^3, 1) / r^2)
d/dr alpha = 1/2 (1 - 2 minMass / r)^(-1/2) * -2 * d/dr(minMass / r)
	= -1/alpha d/dr(minMass / r) dr/dx
dr/dx = x / r

g_ij = x_i x_j / (r^2 (r / minRadius - 1)) + delta_ij
d_kij = 1/2 g_ij,k  = 1/2 (
	x_i,k x_j / (r^2 (r / minRadius - 1))
	+ x_i x_j,k / (r^2 (r / minRadius - 1))
	- x_i x_j / (r^3 (r / minRadius - 1)^2 ) * (
		2 (r / minRadius - 1) + r (
			(minRadius - r d/dr minRadius) / minRadius^2
		)
	)

--]]
				return template([[
#define calc_f(alpha)			(<?=fCCode?>)
#define rSq(x,y,z) 				(x*x + y*y + z*z)
#define r(x,y,z) 				sqrt(rSq(x,y,z))
#define rCubed(x,y,z)			(r(x,y,z) * rSq(x,y,z))
#define bodyMass				<?=bodyMass?>
#define bodyRadius				<?=bodyRadius?>
#define cubed(x)				((x)*(x)*(x))
#define minMass(x,y,z)			(bodyMass * min(rCubed(x,y,z)/cubed(bodyRadius), 1.))
#define minRadius(x,y,z)		(2.*minMass(x,y,z))
#define calc_alpha(x,y,z)		sqrt(1. - 2*minMass(x,y,z)/r(x,y,z))
<?  for i,xi in ipairs(xNames) do
?>#define calc_beta_<?=xi?>(x,y,z)		0.
<? end
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi, xj = xNames[i], xNames[j]
?>#define calc_gamma_<?=xij?>(x,y,z)	(<?=i==j and '1.+' or ''?> + <?=xi?> * <?=xj?> / ((r(x,y,z) / minRadius(x,y,z) - 1.) * rSq(x,y,z)))
<? end
for ij,xij in ipairs(symNames) do
?>#define calc_K_<?=xij?>(x,y,z)		0.
<? end
?>#define dr_minMass_over_r(x,y,z)	(bodyMass * ((r >= bodyRadius ? 0. : 3. * r(x,y,z) / cubed(bodyRadius))- min(rCubed(x,y,z)/cubed(bodyRadius), 1.) / rSq(x,y,z)))
<?
for i,xi in ipairs(xNames) do
?>#define calc_a_<?=xi?>(x,y,z)	(-dr_minMass_over_r(x,y,z) / (calc_alpha(x,y,z) * calc_alpha(x,y,z)) * x / r(x,y,z))
<? end
?>

]], 			table(getTemplateEnv(solver), {
					fCCode = fCCode:match'{ return (.*); }',
					bodyRadius  = clnumber(bodyRadius),
					bodyMass = clnumber(bodyMass),
				}))
			end


			--[[ this is technically correct, but gets some artifacts with the radial boundary on the cartesian grid
			local H = symmath.Heaviside
			--]]
			-- [[ this looks a bit smoother.  sharper edges mean greater artifacts.
			local function H(u) return symmath.tanh((u) * 10) * .5 + .5 end
			--]]
				
			local min = class(require 'symmath.Function')
			min.name = 'min'
			min.func = math.min
			-- derivative wrt 1st param... 
			function min:evaluateDerivative(deriv, ...)
				local a = self[1]
				local b = self[2]
				return H(b - a) * deriv(a, ...) + H(a - b) * deriv(b, ...)
			end
			
			local bodies = args and args.bodies or {{
				pos = {0,0,0},
				mass = .001,
				radius = .1,
			}}

			local x,y,z = symmath.vars('x','y','z')
			
			print('building variables ...')
			
			local alpha = 1
			local gamma = table{1,0,0,1,0,1}
			for _,body in ipairs(bodies) do
				local M = body.mass 
				local R = body.radius
				
				local x_ = x - body.pos[1] 
				local y_ = y - body.pos[2] 
				local z_ = z - body.pos[3] 
				local rSq = x_^2 + y_^2 + z_^2
				local r = rSq^.5
				local m = M * min(r/R, 1)^3
				local R = 2*m
				
				alpha = alpha - 2*m/r
				gamma[1] = gamma[1] + x_^2/((r/R-1)*rSq)
				gamma[2] = gamma[2] + x_*y_/((r/R-1)*rSq)
				gamma[3] = gamma[3] + x_*z_/((r/R-1)*rSq)
				gamma[4] = gamma[4] + y_^2/((r/R-1)*rSq)
				gamma[5] = gamma[5] + y_*z_/((r/R-1)*rSq)
				gamma[6] = gamma[6] + z_^2/((r/R-1)*rSq)
			end
			alpha = alpha^.5
			
			print('...done building variables')
			
			print('initializing numerical relativity variables ...')
			local codes = initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = {x,y,z},
				-- 4D metric ADM components:
				alpha = alpha,
				beta = {0,0,0},
				gamma = gamma,
				K = {0,0,0,0,0,0},
--				useNumericInverse = true,	-- if gamma gets too complex ...
				-- hmm, why do I need gammaU again?  just for initialization of V^i I think ...
				-- hmm would be nice if any field could be a function, algebra, or constant ...
--[[				
				density = function(x,y,z)
					local density = 0
					for _,body in ipairs(bodies) do
						local x_ = x - body.pos[1] 
						local y_ = y - body.pos[2] 
						local z_ = z - body.pos[3]
						local rSq = x_*x_ + y_*y_ + z_*z_
						if rSq < body.radius * body.radius then
							density = density + (body.density or (body.mass / (4/3 * math.pi * body.radius * body.radius * body.radius)))
						end
					end
					return density
				end,
				pressure = function(x,y,z)
					local pressure = 0
					for _,body in ipairs(bodies) do
						local M = body.mass
						local R = body.radius
						local x_ = x - body.pos[1] 
						local y_ = y - body.pos[2] 
						local z_ = z - body.pos[3]
						local rSq = x_*x_ + y_*y_ + z_*z_
						if rSq < body.radius * body.radius then
							local r = math.sqrt(rSq)
							local rho0 = body.pressure or (body.mass / (4/3 * math.pi * body.radius * body.radius * body.radius))
							-- TOV pressure solution: http://physics.stackexchange.com/questions/69953/solving-the-tolman-oppenheimer-volkoff-tov-equation
							-- ... solution to constant pressure?  doesn't density decrease with radius? time to find a better source.
							pressure = pressure + (body.pressure or (rho0 * (
								(
									math.sqrt(1 - 2 * M / R) - math.sqrt(1 - 2 * M * r * r / (R * R * R))
								) / (
									math.sqrt(1 - 2 * M * r * r / (R * R * R)) - 3 * math.sqrt(1 - 2 * M / R)
								)
							)))
						end
					end
					return pressure
				end,
--]]
			}
			print('...done initializing numerical relativity variables') 
			return codes	
		end,
	},
	{
		name = 'stellar model 2',
		getCodePrefix = function(self, solver)
			-- planet plucked out of existence
			initNumRel{
				bodies={
					{pos = {0,0,0}, radius = .1, mass = .001, density=0, pressure=0}
				}
			}
		end,
	},
	{
		name = 'stellar model 3',
		getCodePrefix = function(self, solver)
			--[[
			earth radius = 6.37101e+6 m
			domain: 10x radius = 6.37101e+7 m
			earth mass = 5.9736e+24 kg = 5.9736e+24 * 6.6738480e-11 / 299792458^2 m
			earth mass = Em * G / c^2 in meters
			--]]
			local G = 6.6738480e-11	-- kg m^3/s^2
			local c = 299792458	-- m/s
			-- massInRadii is the order of 1e-9.  much more subtle than the default 1e-3 demo
			local earth = {radiusInM = 6.37101e+6, massInKg = 5.9736e+24}
			-- massInRadii is on the order of 1e-6
			local sun = {radiusInM = 6.960e+8, massInKg = 1.9891e+30}

			local planet = sun
			planet.massInM = planet.massInKg * G / c^2
			planet.massInRadii = planet.massInM / planet.radiusInM
			planet.radiusInCoords = .1
			planet.massInCoords = planet.massInRadii * planet.radiusInCoords
			for k,v in pairs(planet) do print(k,v) end

			local gridUnitsInM = planet.radiusInM / planet.radiusInCoords
			initNumRel{
				bodies={
					{pos = {0,0,0}, radius = planet.radiusInCoords, mass = planet.massInCoords},
				}
			}
		end,
	},
	--[[
	2007 Alic et al "Efficient Implementation of finite volume methods in Numerical Relativity"
	1D Black Hole in wormhole form
	the paper says for >1D use the eqn 24 and do the 'stuffed black hole' fix
	the paper describes the metric in terms of eta, and then gives a relationship between r and eta ...
	but never labels their graphs, whether they are in terms of eta or r ...
	I tried eta and got an initial value of a lapse that looked, well, like a tanh would look (that is the lapse, according to the metric)
	 which is not what their graphs look like.  their graphs look like a sigmoid/tanh that goes from .9 to 1 ... not the range of tanh(eta) which goes from 0 to 1 on a domain of 0 to 8
	--]]
	{
		name = '1D black hole - wormhole form',
		init = function(self, solver)
			solver.eqn:addGuiVar{name='m', value=1}
		end,
		getCodePrefix = function(self, solver, getCodes)
			local m = solver.eqn.guiVars.m.value
		
			-- TODO refreshCodePrefix ... except this is called by that ...
			-- hmm, account for initState's defining the bounds
			solver.mins = require 'vec.vec3'(0,0,0) 
			solver.maxs = require 'vec.vec3'(8,8,8) 
			
			-- coords: t, eta, Omega
			-- g_tt = -(tanh eta)^2 = -alpha^2 <=> alpha = tanh eta
			-- g_eta_eta = 4 m^2 cosh(eta/2)^4
			-- g_Omega_Omega = 4 m^2 cosh(eta/2)^4
			-- using r = m/2 exp(Eta)
			-- m = mass
			-- TODO how should the finite volume scheme be modified to incorporate the volume element -- especially a dynamic volume element
			--local xNames = table{'eta', 'theta', 'phi'}
			-- TODO mapping to isotropic coordinates?
			-- or separate models or automatic transforms between isotropic and spherical?
			local xs = xNames:map(function(x) return symmath.var(x) end)
			
			-- [[ using eta as a coord
			local eta, theta, phi = xs:unpack()
			local alpha = symmath.tanh(eta)
			local g_eta_eta = 4 * m^2 * symmath.cosh(eta/2)^4
			local gamma = {g_eta_eta, 0, 0, g_eta_eta, 0, g_eta_eta}  
			--]]
			--[[
			local x,y,z = xs:unpack()
			local r = (x^2+y^2+z^2)^.5
			-- r = m/2 exp(eta)
			-- eta = log(2r/m)
			-- alpha = tanh(eta) = tanh(log(2r/m)) = sinh(eta) / cosh(eta) = (exp(eta) - exp(-eta)) / (exp(eta) + exp(-eta))
			-- = (2r/m - m/2r) / (2r/m + m/2r)
			-- = (4r^2 - m^2) / (4r^2 + m^2)
			local tanh_eta = (4 * r^2 - m^2) / (4 * r^2 + m^2)
			local alpha = tanh_eta
			-- or maybe the paper just initializes a lapse of 1? 
			--local alpha = 1
			-- 4 m^2 cosh(eta/2)^4 deta^2 = 4 m^2 cosh(log(2r/m)/2)^4 /r^2 dr^2
			-- = 4 (m/r)^2 (sqrt(2r/m) + sqrt(m/2r) )^4 dr^2
			-- but now we have an infinity as r->0 ...
			local g_rr = 4 * (m/r)^2 * (symmath.sqrt(2*r/m) + symmath.sqrt(.5*m/r))^4
			-- ... unless that "cosh eta / 2" means "cosh(eta)/2" .... but then why not just write 4/16 = 1/4 out front?  grrr this paper  ....
			-- 1/4 (m/r)^2 (2r/m + m/2r)^4
			--local g_rr = (m/r)^2 * (2*r/m + m/(2*r))^4 / 4
			local gamma = {g_rr, 0, 0, g_rr, 0, g_rr}
			--]]
			
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = alpha,
				gamma = gamma,
				K = {0,0,0,0,0,0},
			}
		end,
	},
--[[
2003 Alcubierre et al "Toward standard testbeds for numerical relativity"
2004 Bona et al "A symmetry-breaking mechanism..."
2009 Bona, Bona-Casas "Gowdy waves..." https://arxiv.org/pdf/0911.1208v1.pdf
Gaudy wave tests
1/sqrt(t) exp(Q/2) (-dt^2 + dz^2) + t (exp(P) dx^2 + exp(-P) dy^2)
Q & P are functions of t & z, only periodic in z
t = t0 exp(-tau/tau0)
dt = -t0/tau0 exp(-tau/tau0) dtau

specifically for this problem:
P = J0(2 pi t) cos(2 pi z)
Q = pi J0(2 pi) J1(2 pi) - 2 pi t J0(2 pi t) J1(2 pi t) cos(2 pi z)^2
		+ 2 pi^2 t^2 (J0(2 pi t)^2 + J1(2 pi t)^2 - J0(2 pi)^2 - J1(2 pi)^2)
2 pi t0 is the 20th root of the Bessel function J0 <=> t0 ~ 9.88

at time tau=0 <=> t=t0, and picking 2 pi t0 to be a root of J0 ...
Q = pi J0(2 pi) J1(2 pi) - 2 pi^2 t0^2 (J0(2 pi)^2 + J1(2 pi)^2)
--]]
	{
		name = 'Gowdy waves',
		getCodePrefix = function(self, solver)
			--solver.mins = require 'vec.vec3'(0,0,0)
			--solver.maxs = require 'vec.vec3'(10,10,10)
		end,
	},

-- 2003 Alcubierre et al "Toward standard testbeds for numerical relativity"
-- grids are from [-.5, .5] with 50 rho cells, rho = 1,2,4, and timesteps of .01/rho
	{
		name = 'testbed - robust',
		-- pick epsilon so epsilon^2 = 0
		-- eqn 4.1: pick epsilon from -1e-10 / rho^2 to 1e=10 / rho^2
		init = function(self, solver)
			solver.eqn:addGuiVar{name='epsilon', value=1e-10}
		end,
		resetState = function(self, solver)
			local epsilon = solver.eqn.guiVars.epsilon.value
			solver.eqn:fillRandom(epsilon)
		end,
		refreshInitStateProgram = function(self, solver) end,
		getInitStateCode = function(self, solver) end,
	},
	{
		name = 'testbed - gauge wave',
		init = function(self, solver)
--			solver.mins = vec3(-.5, -.5, -.5)
--			solver.maxs = vec3(-.5, -.5, -.5)
			solver:setBoundaryMethods'periodic'
			solver.eqn:addGuiVar{name='A', value=.1}	-- .1, .01
			solver.eqn:addGuiVar{name='d', value=1}
--			solver.eqn.guiVars.f.value = solver.eqn.guiVars.f.options:find'1'	-- set f=1
		end,
		getCodePrefix = function(self, solver, getCodes)
			local A = solver.eqn.guiVars.A.value
			local d = solver.eqn.guiVars.d.value
			local xs = xNames:map(function(x) return symmath.var(x) end)
			local x,y,z = xs:unpack()
			local t = 0
			local theta = (2 * math.pi / d) * (x - t)
			local H = 1 + A * symmath.sin(theta)
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = symmath.sqrt(H),
				gamma = {H, 0, 0, 1, 0, 1},
				K = {-math.pi * A / d * symmath.cos(theta) / symmath.sqrt(H), 0,0,0,0,0},
			}
		end,
	},
	{
		name = 'testbed - gauge wave - diagonal',
		init = function(self, solver)
			error"finishme"
		end,
	},
	{
		name = 'testbed - linear wave',
		init = function(self, solver)
-- TODO changing the range upon init causes something to freeze up ...
--			solver.mins = vec3(-.5, -.5, -.5)
--			solver.maxs = vec3(-.5, -.5, -.5)
			solver:setBoundaryMethods'periodic'
			solver.eqn:addGuiVar{name='A', value=1e-8}
			solver.eqn:addGuiVar{name='d', value=1}
		end,
		getCodePrefix = function(self, solver, getCodes)
			local A = solver.eqn.guiVars.A.value
			local d = solver.eqn.guiVars.d.value
			local xs = xNames:map(function(x) return symmath.var(x) end)
			local x,y,z = xs:unpack()
			local t = 0
			local theta = (2 * math.pi / d) * (x - t)
			local b = A * symmath.sin(theta)
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = 1,
				gamma = {1,0,0,1+b,0,1-b},
				K = {0,0,0,b:diff(t)/2,0,-b:diff(t)/2},
			}
		end,
	},
	{
		name = 'testbed - linear wave - diagonal',
		init = function(self, solver)
			error"finishme"
		end,
	},
	{
		name = 'testbed - Gowdy',
		init = function(self, solver)
		
		end,
		getCodePrefix = function(self, solver, getCodes)
			local frac = symmath.frac
			local xs = xNames:map(function(x) return symmath.var(x) end)
			local x,y,z = xs:unpack()
			local t = 0
			--[[ general solution
			local lambda = -2 * pi * t * J0(2 * pi * t) * J1(2 * pi * t) * symmath.cos(2 * pi * z)^2
				+ 2 * pi^2 * t^2 * (J0(2 * pi * t)^2 + J1(2 * pi * t)^2)
				- 1/2 * ( (2 * pi)^2 * (J0(2 * pi)^2 + J1(2 * pi)^2) - 2 * pi * J0(2 * pi) * J1(2 * pi))
			--]]
			local lambda = -.5 * ( (2 * math.pi)^2 * (J0(2 * math.pi)^2 + J1(2 * math.pi)^2) - 2 * math.pi * J0(2 * math.pi) * J1(2 * math.pi))
			local alpha = t^symmath.frac(-1,4) * symmath.exp(lambda/4)
			local P = J0(2 * pi * t) * symmath.cos(2 * pi * z)
			local gamma = {t * symmath.exp(P), 0, 0, t * symmath.exp(-P), 0, 1/symmath.sqrt(t) * symmath.exp(lambda/2)}
			local K = {-1/2 * t^(1/4) * symmath.exp(P-lambda/4) * (1 + t * dP_dt), 0, 0, -1/2 * t^(1/4) * symmath.exp(-P-lambda/4) * (1 - t * dP_dt), 0, 1/4*t^(-1/4) * symmath.exp(lambda/4) * (1/t - dlambda_dt)}
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = alpha,
				gamma = gamma,
				K = K,
			}
		end,
	},
}:map(function(cl)
	return class(NumRelInitCond, cl)
end)
