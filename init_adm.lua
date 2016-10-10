local symmath = require 'symmath'
local table = require 'ext.table'

--symmath.tostring = require 'symmath.tostring.SingleLine'		

local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
local function from6to3x3(i)
	return table.unpack(from6to3x3_table[i])
end


-- I'm working on a unified initial condition code for 1D and 3D NR problems:
--[[
args:
	args = {x,y,z} spatial basis variables
	alpha = lapse expression
	(beta isn't used yet)
	gammaLL = 3-metric
	K = extrinsic curvature

	gammaLL & K are stored as {xx,xy,xz,yy,yz,zz}

	density = lua function
--]]
local function initNumRel(args)	
	local vars = assert(args.vars)
	
	local exprs = table{
		alpha = assert(args.alpha),
		gammaLL = {table.unpack(args.gammaLL)},
		K = {table.unpack(args.K)}
	}
	assert(#exprs.gammaLL == 6)
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
		exprs.gammaLL = exprs.gammaLL[1]	-- only need g_xx
		exprs.a = (exprs.alpha:diff(vars[1]) / exprs.alpha)()	-- only need a_x
		exprs.d = (exprs.gammaLL:diff(vars[1])/2)()	-- only need D_xxx
		exprs.K = exprs.K[1]	-- only need K_xx
		print('compiling expressions...')
		local codes = table.map(exprs, compileC)
		print('...done compiling expressions')
		
		return {
			alpha = codes.alpha,
			gamma_xx = codes.gammaLL,
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
			exprs.gammaU = {mat33.inv(exprs.gammaLL:unpack())}
			print('...done inverting spatial metric')
		end

		-- this takes forever.  why is that?  differentiation?
		print('building metric partials...')
		exprs.D = table.map(vars, function(x_k)
			return table.map(exprs.gammaLL, function(g_ij)
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
			local gammaLL = calc.gammaLL:map(function(g_ij) return g_ij(x,y,z) end)
			local A = calc.A:map(function(A_i) return A_i(x,y,z) end)
			local D = calc.D:map(function(D_i) return D_i:map(function(D_ijk) return D_ijk(x,y,z) end) end)
			local gammaU = args.useNumericInverse and table{mat33.inv(gammaLL:unpack())} or calc.gammaU:map(function(gammaUij) return gammaUij(x,y,z) end) 
			
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
				gammaLL[1], gammaLL[2], gammaLL[3], gammaLL[4], gammaLL[5], gammaLL[6],
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

return {
	{
		name = 'gauge shock wave',
		init = function(solver, args)
			local var = symmath.var
			
			local xNames = table{'x', 'y', 'z'}
	
			-- symmetric indexes: xx xy xz yy yz zz
			local symNames = table()
			for i,xi in ipairs(xNames) do
				for j=i,3 do
					local xj = xNames[j]
					symNames:insert(xi..xj)
				end
			end


			-- lapse function

			local alphaVar = args.alphaVar or var'alpha'
			--local kappa = 1 -- TODO specify this? or should I bother?
			local f = symmath.clone(args.f or 1)
			
			local dalpha_f = f:diff(alphaVar)()
	
			-- the rest

			local Tensor = require 'symmath.Tensor'
			local xs = xNames:map(function(x) return var(x) end)
			Tensor.coords{{variables=xs}}
			local x,y,z = xs:unpack()

			local size = solver.maxs[1] - solver.mins[1]
			local H = 5 * size / 300
			local sigma = 10 * size / 300
		
			local xc = .5 * (solver.mins[1] + solver.maxs[1])
			local yc = .5 * (solver.mins[2] + solver.maxs[2])
			local zc = .5 * (solver.mins[3] + solver.maxs[3])
			
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
			h = h()
local dh = Tensor('_i', function(i) return h:diff(xs[i])() end)
local d2h = Tensor('_ij', function(i,j) return h:diff(xs[i],xs[j])() end)
-- hmm why isn't comma derivative working?
--print('h',h)
--print('h_,i', h'_,i')
--print('h_,i', h'_,i'())
--os.exit()

			local deltaLL = Tensor('_ij', function(i,j) return i == j and 1 or 0 end)
			local gammaLL = (deltaLL'_ij' - dh'_i' * dh'_j')()
			local KLL = (-d2h'_ij' / (1 - (dh'_k' * dh'_k')() )^.5 )()
			
			print('...done deriving and compiling.')
			
			local codes = initNumRel{
				vars = xs,
				alpha = alpha,
				gammaLL = symNames:map(function(xij,ij) return gammaLL[{from6to3x3(ij)}] end),
				K = symNames:map(function(xij,ij) return KLL[{from6to3x3(ij)}] end),
			}
	
			codes.f = require 'symmath.tostring.C':compile(f, {alphaVar})
			codes.dalpha_f = require 'symmath.tostring.C':compile(dalpha_f, {alphaVar})

			return codes
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
