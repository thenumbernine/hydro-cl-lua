-- based on the book "Introduction to 3+1 Numerical Relativity" and on the paper "Introduction to Numerical Relativity", both by Alcubierre

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'
local symmath = require 'symmath'
symmath.tostring = require 'symmath.tostring.SingleLine'		

local ADM_BonaMasso_1D_3Var = class(Equation)
ADM_BonaMasso_1D_3Var.name = 'ADM_BonaMasso_1D_3Var' 

ADM_BonaMasso_1D_3Var.numStates = 5
ADM_BonaMasso_1D_3Var.numWaves = 3	-- alpha and gamma_xx are source-term only

ADM_BonaMasso_1D_3Var.consVars = {'alpha', 'gamma_xx', 'a_x', 'd_xxx', 'KTilde_xx'}
ADM_BonaMasso_1D_3Var.mirrorVars = {{'gamma_xx', 'a_x', 'd_xxx', 'KTilde_xx'}}

ADM_BonaMasso_1D_3Var.useSourceTerm = true

ADM_BonaMasso_1D_3Var.initStateNames = {'gaussianWave'}

ADM_BonaMasso_1D_3Var.displayVars = table()
	:append(ADM_BonaMasso_1D_3Var.consVars)
	:append{'K_xx', 'volume'}

-- I'm working on a unified initial condition code for 1D and 3D NR problems:
--[[
args:
	args = {x,y,z} spatial basis variables
	alpha = lapse expression
	(beta isn't used yet)
	gamma = 3-metric
	K = extrinsic curvature

	gamma & K are stored as {xx,xy,xz,yy,yz,zz}

	density = lua function
--]]
local function initNumRel(args)	
	local vars = assert(args.vars)
	
	local exprs = table{
		alpha = assert(args.alpha),
		gamma = {table.unpack(args.gamma)},
		K = {table.unpack(args.K)}
	}
	assert(#exprs.gamma == 6)
	assert(#exprs.K == 6)

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
	
	local function compileC(expr, name)
		assert(type(expr) == 'table')
		if symmath.Expression.is(expr) then 
			expr = expr()
			print('compiling '..expr)
			return require 'symmath.tostring.C':compile(expr, vars), name
		end
		return table.map(expr, compileC), name
	end
	
	-- ADM
--	if solverName == 'ADM1DRoe' then
		exprs.gamma = exprs.gamma[1]	-- only need g_xx
		exprs.a = (exprs.alpha:diff(vars[1]) / exprs.alpha)()	-- only need a_x
		exprs.d = (exprs.gamma:diff(vars[1])/2)()	-- only need D_xxx
		exprs.K = exprs.K[1]	-- only need K_xx
		print('compiling expressions...')
		local codes = table.map(exprs, compileC)
		print('...done compiling expressions')
		
		local xs = table{'x', 'y', 'z'}
		local sym = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}
	
		return {
			alpha = codes.alpha,
			gamma_xx = codes.gamma,
			a_x = codes.a,
			d_xxx = codes.d,
			K_xx = codes.K
		}
--[[	
	elseif solverName == 'ADM3DRoe' then
		
		-- for complex computations it might be handy to extract the determinant first ...
		-- or even just perform a numerical inverse ...
		if not args.useNumericInverse then
			print('inverting spatial metric...')
			exprs.gammaU = {mat33.inv(exprs.gamma:unpack())}
			print('...done inverting spatial metric')
		end

		-- this takes forever.  why is that?  differentiation?
		print('building metric partials...')
		exprs.D = table.map(vars, function(x_k)
			return table.map(exprs.gamma, function(g_ij)
				print('differentiating '..g_ij)
				return (g_ij:diff(x_k)/2)()
			end)
		end)
		print('...done building metric partials')
	
		print('building lapse partials...')
		exprs.A = table.map(vars, function(var)
			return (exprs.alpha:diff(var) / exprs.alpha)()
		end)
		print('...done building lapse partials')
		
		print('compiling expressions...')
		local calc = table.map(exprs, buildCalc)
		print('...done compiling expressions')
	
		local densityFunc = args.density or 0
		if symmath.Expression.is(densityFunc) then
			densityFunc = densityFunc():compile(vars)
		end

		local pressureFunc = args.pressure or 0
		if symmath.Expression.is(pressureFunc) then
			pressureFunc = pressureFunc():compile(vars)
		end

		initState = function(x,y,z)
		
			local alpha = calc.alpha(x,y,z)
			local gamma = calc.gamma:map(function(g_ij) return g_ij(x,y,z) end)
			local A = calc.A:map(function(A_i) return A_i(x,y,z) end)
			local D = calc.D:map(function(D_i) return D_i:map(function(D_ijk) return D_ijk(x,y,z) end) end)
			local gammaU = args.useNumericInverse and table{mat33.inv(gamma:unpack())} or calc.gammaU:map(function(gammaUij) return gammaUij(x,y,z) end) 
			
			local function sym3x3(m,i,j)
				local m_xx, m_xy, m_xz, m_yy, m_yz, m_zz = m:unpack()
				if i==1 then
					if j==1 then return m_xx end
					if j==2 then return m_xy end
					if j==3 then return m_xz end
				elseif i==2 then
					if j==1 then return m_xy end
					if j==2 then return m_yy end
					if j==3 then return m_yz end
				elseif i==3 then
					if j==1 then return m_xz end
					if j==2 then return m_yz end
					if j==3 then return m_zz end
				end
				error'here'
			end
			local V = range(3):map(function(i)
				local s = 0
				for j=1,3 do
					for k=1,3 do
						local D_ijk = sym3x3(D[i],j,k)
						local D_kji = sym3x3(D[k],j,i)
						local gammaUjk = sym3x3(gammaU,j,k)
						local dg = (D_ijk - D_kji) * gammaUjk
						s = s + dg
					end
				end
				return s
			end)
			local K = {}
			for i=1,6 do
				K[i] = calc.K[i](x,y,z)
			end

			local density = 0
			if type(densityFunc) == 'number' then
				density = densityFunc
			elseif type(densityFunc) == 'function' then
				density = densityFunc(x,y,z)
			elseif densityFunc ~= nil then
				error("don't know how to handle density")
			end
			
			local pressure = 0
			if type(pressureFunc) == 'number' then
				pressure = pressureFunc
			elseif type(pressureFunc) == 'function' then
				pressure = pressureFunc(x,y,z)
			elseif pressureFunc ~= nil then
				error("don't know how to handle pressure")
			end
		
			local velocityX, velocityY, velocityZ
			if args.velocity then velocityX, velocityY, velocityZ = args.velocity(x,y,z) end
			velocityX = velocityX or 0
			velocityY = velocityY or 0
			velocityZ = velocityZ or 0

			return
				alpha,
				gamma[1], gamma[2], gamma[3], gamma[4], gamma[5], gamma[6],
				A[1], A[2], A[3],
				D[1][1], D[1][2], D[1][3], D[1][4], D[1][5], D[1][6],
				D[2][1], D[2][2], D[2][3], D[2][4], D[2][5], D[2][6],
				D[3][1], D[3][2], D[3][3], D[3][4], D[3][5], D[3][6],
				K[1], K[2], K[3], K[4], K[5], K[6],
				V[1], V[2], V[3],
				density,
				velocityX, velocityY, velocityZ,
				pressure
		end
	end
--]]
end

ADM_BonaMasso_1D_3Var.initStates = {
	{
		name = 'gauge shock wave',
		init = function(solver)
			local alphaVar = symmath.var'alpha'

			--local f = 1
			--local f = 1.69
			--local f = .49
			local f = 1 + 1 / alphaVar^2
			local dalpha_f = f:diff(alphaVar)()
		
			local Tensor = require 'symmath.Tensor'

			local H, sigma
			-- [[ not unit domain
			H = 5
			sigma = 10
			--]]
			--[[ unit domain
			H = 1/1000
			sigma = 1/10
			--]]
		
			local xc = .5 * (solver.mins[1] + solver.maxs[1])
			local yc = .5 * (solver.mins[2] + solver.maxs[2])
			local zc = .5 * (solver.mins[3] + solver.maxs[3])
			local x, y, z = symmath.vars('x', 'y', 'z')
			Tensor.coords{{variables={x,y,z}}}
			
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

			local h = H * symmath.exp(-((x - xc)^2 + (y - yc)^2 + (z - zc)^2) / sigma^2)
			gamma = {
				1 - h:diff(x)^2,
				-h:diff(x) * h:diff(y),
				-h:diff(x) * h:diff(z),
				1 - h:diff(y)^2,
				-h:diff(y) * h:diff(z),
				1 - h:diff(z)^2,
			}
			local div_h = h:diff(x)^2 + h:diff(y)^2 + h:diff(z)^2
			local K_denom = (1 - div_h)^.5
			K = {
				-h:diff(x,x) / K_denom,
				-h:diff(x,y) / K_denom,
				-h:diff(x,z) / K_denom,
				-h:diff(y,y) / K_denom,
				-h:diff(y,z) / K_denom,
				-h:diff(z,z) / K_denom,
			}

			print('...done deriving and compiling.')
			
			return initNumRel{vars={x,y,z}, alpha=alpha, gamma=gamma, K=K}
		end,
	},
	{
		name = 'Alcubierre warp bubble',
		init = function(solver)
		end,
	},
	{
		name = 'Schwarzschild black hole',
		init = function(solver)
		end,
	},
	{
		name = 'stellar model',
		init = function(solver)
			setup()
		end,
	},
	{
		name = 'stellar model 2',
		init = function(solver)
			-- planet plucked out of existence
			setup{
				bodies={
					{pos = {0,0,0}, radius = .1, mass = .001, density=0, pressure=0}
				}
			}
		end,
	},
	{
		name = 'stellar model 3',
		init = function(solver)
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
			setup{
				bodies={
					{pos = {0,0,0}, radius = planet.radiusInCoords, mass = planet.massInCoords},
				}
			}
		end,
	},
}
ADM_BonaMasso_1D_3Var.initStateNames = table.map(ADM_BonaMasso_1D_3Var.initStates, function(initState) return initState.name end)

function ADM_BonaMasso_1D_3Var:getInitStateCode(solver)

--[[ moving it all here, and making it 3d-capable
	local initState = self.initStates[solver.initStatePtr[0]+1]
	local codes = initState.init(solver)
--]]
-- [[
		
	local x,y,z = symmath.vars('x','y','z')
	local alphaVar = symmath.var'alpha'
	
	local xc = 150
	local H = 5
	local sigma = 10
	
	--if self.initState[0] == 0 then ...
	local h = H * symmath.exp(-(x - xc)^2 / sigma^2)
	
	local kappa = 1
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
	local gamma_xx = 1 - h:diff(x)^2
	local K_xx = -h:diff(x,x) / gamma_xx^.5

	--local f = 1
	--local f = 1.69
	--local f = .49
	local f = 1 + kappa / alphaVar^2

	-- above this line is the initial condition expressions
	-- below is the symmath codegen

	local function makesym(expr, k)
		return symmath.clone(expr):simplify(), k
	end
	
	local exprs = table{
		alpha = alpha,
		gamma_xx = gamma_xx,
		K_xx = K_xx,
	}:map(makesym)
	exprs.a_x = (exprs.alpha:diff(x) / exprs.alpha)()
	exprs.d_xxx = (exprs.gamma_xx:diff(x) / 2)()

	local function compileC(expr, args)
print('compile',expr,table.map(args,tostring):concat', ')
		local code = require 'symmath.tostring.C':compile(expr, args)
print('...to',code)
		return code
	end

	local codes = exprs:map(function(expr, name)
print(name)		
		return compileC(expr, {x,y,z}), name
	end)
	
	f = makesym(f)
print'f'
	codes.f = compileC(f, {alphaVar})

print'dalpha_f'
	local dalpha_f = f:diff(alphaVar):simplify()
	codes.dalpha_f = compileC(dalpha_f, {alphaVar})
--]]

	self.codes = codes

	return table{
		self:codePrefix(),
		[[
__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real4 xs = CELL_X(i);
	real x = xs[0];
	real y = 0;
	real z = 0;
	__global cons_t* U = UBuf + index;
	
	U->alpha = init_calc_alpha(x,y,z);
	U->gamma_xx = init_calc_gamma_xx(x,y,z);
	U->a_x = init_calc_a_x(x,y,z);
	U->d_xxx = init_calc_d_xxx(x,y,z);
	real K_xx = init_calc_K_xx(x,y,z);
	U->KTilde_xx = K_xx / sqrt(U->gamma_xx);
}
]],
	}:concat'\n'
end

--[[
this is called by getInitStateCode and by solverCode
getInitStateCode is called when initState changes, which rebuilds this, since it is dependent on initState (even though I only have one right now)
solverCode also uses this, but doesn't calculate it
so this will all only work right so long as solverCode() is only ever called after getInitStateCode()
--]]
function ADM_BonaMasso_1D_3Var:codePrefix()
	return table.map(self.codes, function(code,name,t)
		return 'real init_calc_'..name..code, #t+1
	end):concat'\n'
end

function ADM_BonaMasso_1D_3Var:solverCode()
	return table{
		self:codePrefix(),
		'#include "adm1d3to5var.cl"',
	}:concat'\n'
end

ADM_BonaMasso_1D_3Var.eigenVars = {'f'}
function ADM_BonaMasso_1D_3Var:getEigenInfo()
	local makeStruct = require 'makestruct'
	return {
		typeCode =
			makeStruct('eigen_t', self.eigenVars) .. '\n' ..
			makeStruct('fluxXform_t', {'alpha', 'gamma_xx', 'f'}),
		-- don't use the default matrix stuff. instead it'll be provided in adm1d3to5var.cl
		code = nil,
		displayVars = self.eigenVars,
	}
end

return ADM_BonaMasso_1D_3Var 
