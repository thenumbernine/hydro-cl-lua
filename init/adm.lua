local symmath = require 'symmath'
local table = require 'ext.table'

--symmath.tostring = require 'symmath.tostring.SingleLine'		

local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
local function from6to3x3(i)
	return table.unpack(from6to3x3_table[i])
end

local xNames = table{'x', 'y', 'z'}

-- symmetric indexes: xx xy xz yy yz zz
local symNames = table()
for i,xi in ipairs(xNames) do
	for j=i,3 do
		local xj = xNames[j]
		symNames:insert(xi..xj)
	end
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
		
code = code:gsub('sqrt%(', '(real)sqrt((real)')
		
		print(code)
		return code
	end
	return table.map(expr, function(v,k) return compileC(v,k,vars) end), name
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
	
	-- ADM
	
	-- assuming this is one of the ADM 1D solvers
	print('compiling expressions...')
	local codes = table.map(
		assert(args.getCodes(exprs, vars, args), "getCodes needs to return something"),
		function(v,k) return compileC(v,k,vars) end)
	print('...done compiling expressions')
	
	-- here's the lapse conditions 

	local alphaVar = symmath.var'alpha'
	-- TODO each eqn must be stated as a guiVar of the solver
	-- make a parent class or something for all num rel eqns
	local fGuiVar = args.solver.eqn.guiVarsForName.f
	local fLuaCode = fGuiVar.options[fGuiVar.value[0]+1]
	
	local f = assert(loadstring('local alpha = ... return '..fLuaCode))(alphaVar)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	codes.f = compileC(f, 'f', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})

	return table.map(codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

return {
	{
		name = 'gauge shock wave',
		init = function(solver, getCodes)

			-- here's the coordinates

			local xs = xNames:map(function(x) return symmath.var(x) end)
			symmath.Tensor.coords{{variables=xs}}
			local x,y,z = xs:unpack()

			-- here's the metric

			local size = solver.maxs[1] - solver.mins[1]
			local H = 5 * size / 300
			local sigma = 10 * size / 300
		
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
	{
		name = 'Alcubierre warp bubble',
		init = function(solver, getCodes)
			local R = .5		-- warp bubble radius
			local sigma = 8	-- warp bubble thickness
			local speed = .1	-- warp bubble speed

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
		init = function(solver, getCodes)
			local R = .002	-- Schwarzschild radius
			
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
		end,
	},
	{
		name = 'stellar model',
		init = function(solver, getCodes)
				
			--[[ this is technically correct, but gets some artifacts with the radial boundary on the cartesian grid
			local H = symmath.Heaviside
			--]]
			-- [[ this looks a bit smoother.  sharper edges mean greater artifacts.
			local function H(u) return symmath.tanh((u) * 10) * .5 + .5 end
			--]]
				
			local class = require 'ext.class'
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
