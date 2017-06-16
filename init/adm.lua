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


-- I'm working on a unified initial condition code for 1D and 3D NR problems:
--[[
args:
	solver = the solver

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
		assert(type(expr) == 'table', "expected table, found "..type(expr))
		if symmath.Expression.is(expr) then 
			expr = expr()
			print('compiling '..expr)
			return expr:compile(vars, 'C'), name
		end
		return table.map(expr, compileC), name
	end
	
	-- ADM
	
	-- assuming this is one of the ADM 1D solvers
	if args.solver.eqn.numWaves == 3 then
		print('compiling expressions...')
		local codes = table{
			alpha  = exprs.alpha,
			gamma_xx = exprs.gammaLL[1],	-- only need g_xx
			a_x = (exprs.alpha:diff(vars[1]) / exprs.alpha)(),	-- only need a_x
			d_xxx = (exprs.gammaLL[1]:diff(vars[1])/2)(),	-- only need D_xxx
			K_xx = exprs.K[1],	-- only need K_xx
		}:map(compileC)
		print('...done compiling expressions')
		return codes

	-- assuming this is an ADM 3D solver
	elseif args.solver.eqn.numWaves == 30
	or args.solver.eqn.numWaves == 45
	then 

		-- for complex computations it might be handy to extract the determinant first ...
		-- or even just perform a numerical inverse ...
		if not args.useNumericInverse then
			print('inverting spatial metric...')
			exprs.gammaUU = {symMat33Inv(exprs.gammaLL:unpack())}
			print('...done inverting spatial metric')
		end

		-- this takes forever.  why is that?  differentiation?
		print('building metric partials...')
		exprs.d = table.map(vars, function(xk)
			return table.map(exprs.gammaLL, function(gamma_ij)
				print('differentiating '..gamma_ij)
				return (gamma_ij:diff(xk)/2)()
			end)
		end)
		print('...done building metric partials')
	
		print('building lapse partials...')
		exprs.a = table.map(vars, function(var)
			return (exprs.alpha:diff(var) / exprs.alpha)()
		end)
		print('...done building lapse partials')
	
		-- this requires turning d from sym into a Tensor
		--  and requires gammaUU being stored as the Tensor metric inverse
		--exprs.V = (d'_ik^k' - d'^k_ki')()

		print('compiling expressions...')
		local codes = {}
		local codes = table(
			--f = exprs.f,
			--dalpha_f = exprs.dalpha_f,
			{alpha = exprs.alpha},
			symNames:map(function(xij,ij)
				return exprs.gammaLL[ij], 'gamma_'..xij
			end),
			xNames:map(function(xi,i)
				return exprs.a[i], 'a_'..xi
			end),
			table(xNames:map(function(xk,k,t)
				return symNames:map(function(xij,ij)
					return exprs.d[k][ij], 'd_'..xk..xij
				end), #t+1
			end):unpack()),
			symNames:map(function(xij,ij)
				return exprs.K[ij], 'K_'..xij
			end)
		):map(compileC)

		print('...done compiling expressions')

		-- TODO if useNumericInverse is true then add a gammaUij that's based on Cramer's rule

		return codes
	else
		error "don't have support for this ADM solver yet"
	end
end

return {
	{
		name = 'gauge shock wave',
		init = function(solver, args)
			local var = symmath.var
			
			-- lapse function

			local alphaVar = args.alphaVar or var'alpha'
			--local kappa = 1 -- TODO specify this? or should I bother?
			local f = symmath.clone(args.f or 1)
			
			local dalpha_f = f:diff(alphaVar)()
	
			-- the rest

			local Tensor = symmath.Tensor
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
				solver = solver,
				vars = xs,
				alpha = alpha,
				gammaLL = symNames:map(function(xij,ij) return gammaLL[{from6to3x3(ij)}] end),
				K = symNames:map(function(xij,ij) return KLL[{from6to3x3(ij)}] end),
			}
	
			codes.f = f:compile({alphaVar}, 'C')
			codes.dalpha_f = dalpha_f:compile({alphaVar}, 'C')

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
