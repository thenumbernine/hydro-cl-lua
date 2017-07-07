--[[

This tells us the coordinate chart for our embedding in nD Cartesian (Euclidian?) geometry.
There are a few options on how to do this.

1) The Engineer way:
	Represent our eqn vector coordinates in Cartesian coordinates,
	use our Geometry to determine the cell volumes, areas, centers, normals, etc.,
	compute flux in Cartesian coordinates.
	
	This seems like the least change from a Cartesian grid, and should be easy to implement:
	All you have to change is the face normals and the dx values of the finite volume update 
	and change dx's in any operators of divergence (div B=0 magnetic monopole constraint) and Laplacian (del rho = 4 pi G for gravitation).
	
	This treats the problem as if it were made up a rigid grid, i.e. flattening off the curves in cylindrical and treating each cell like a polygon.
	For that reason, face center positions should probably be calculated as the average of vertices rather than the coordinate intermediate position,
	otherwise they won't represent the flat geometry and will introduce errors.
	This misses out on the perk of simulating problems whose components are purely rotational about an origin
	with much greater accuracy than a Cartesian grid.

2) The Physicist way:
	Represent our eqn vector coordinates in anholonomic coordinates, normalizing the vector components.
	Ex. for cylindrical the holonomic basis is e_r = [cos phi, sin phi], e_phi = [-r sin phi, r cos phi], 
	while the anholonomic normalized basis is only different by e_phiHat = 1/r e_phi = [-sin phi, cos phi].
	This creates orthonormal basis vectors and creates an identity metric (which means the equations don't have to be changed)
	but requires keeping track of both the holonomic and anholonomic coordinates,
	and requires lots of extra changes:
	Finite volume update requires scaling the flux up by the volume *and* back down by the anholonomic normalization.
	The dx's should be coordinate based, and therefore not metric-adjusted?
	The coordinate conversion needs the holonomic basis information still.
	The centrifugal forces need the anholonomic connections.
	The grid length and distance calculations need the holonomic connections.
	It looks like a lot of extra bookkeeping.
	Buuut the Euler fluid equations at least stay untouched (except the added Coriolis source term)
		
	Most sources take that a step further for the Euler fluid equations and say to represent the velocity in covariant form
	so that the pressure term contribution to the flux needs no rescaling by the holonomic cell volume (while everything else does rescale ... how does the eigen-decomposition represent that?)
	by making the tensor next to that pressure term a delta^i_j instead of a g^ij.
	This seems great for a single equation, but might make a Roe implementation problematic if I want to make it flexible. 
	
	This also means representing covariant derivatives (divergence and Laplace-Beltrami)
	with connections / holonomic volume derivatives *with* extra rescaling to the anholonomic basis.

3) The Mathematician way:
	Represent our eqn vector coordinates as holonomic coordinates.  No normalization.
	This would mean keeping that g^ij in front of the pressure term, and therefore adjusting the eigen-decomposition of the Euler fluid equations.
	It would also mean performing our dx's in coordinate space.  No anholonomic readjustments.

	This is like the Physicist way except *without* extra normalization for the holonomic-anholonomic basis exchanges.
	It is also made to work with a metric, and therefore works easiest when the metric is a dynamic variable (right?  where will the change-in-metric stuff have to go?)
	Turns out this is similar to the conservative form for relativistic solvers mentioned in papers by Font and company.

TODO
charts per dimension
that way every class isn't repeating 1D
and that way we can have multiple 2D options per-permutation of 3D variables

places where geometry is used

- the initial conditions ... should params be specified in coordinate chart (r,theta,phi) or basis (x,y,z) ?   ... coordinate chart for now
- the vector structures -- same question?  coordinate chart
- the flux
- the calcDT
- anything dealing with cell volume or surface
- anything dealing with lengths
- anything dealing with coordinates




here's another rehashing of my options
- contra- vs co-variant form
	this is equation specific
	but certain eqns like Euler have the benefit of turning g's into deltas

- holonomic vs anholonomic
	This is whether vector components are normalized as they are transported through normals.
	Conceptually easier, but mathematically more difficult when working with covariant derivatives in curved space.

- coordinate components vs curvilinear components
	One easy fix would be to just use geom for computing normals and volumes and surface areas,
	and just keep all coordinates in Euclidian space.
	This would certainly be easiest when adopting the equations to unstructured meshes.
	This does throw out any benefits of when velocity streamlines coincide with coordinate derivatives. 

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'

local Geometry = class()

function Geometry:init(args)
	local symmath = require 'symmath'
	local const = symmath.Constant
	
	local dim = args.solver.dim
	local var = symmath.var
	local vars = symmath.vars
	local Matrix = symmath.Matrix
	local Tensor = symmath.Tensor
	
	local flatMetric = Matrix:lambda({dim, dim}, function(i,j) return i==j and 1 or 0 end)
	local coords = args.coords
	local embedded = args.embedded
	
	Tensor.coords{
		{variables=coords},
		{variables=embedded, symbols='IJKLMN', metric=flatMetric},
	}

	local baseCoords = table.map(coords, function(coord)
		return coord.base or coord
	end)

--print('coordinates:', table.unpack(coords))
--print('base coords:', table.unpack(baseCoords))
--print('embedding:', table.unpack(embedded))
	
	local eta = Tensor('_IJ', table.unpack(flatMetric)) 
--print'flat metric:'
--print(var'\\eta''_IJ':eq(eta'_IJ'()))
--print()

	local u = args.chart()
print'coordinate chart:'
print(var'u''^I':eq(u'^I'()))
print()

	local e = Tensor'_u^I'
	e['_u^I'] = u'^I_,u'()
print'embedded:'
print(var'e''_u^I':eq(var'u''^I_,u'):eq(e'_u^I'()))
print()

	local anholonomic
	for i=1,#coords do
		if coords[i].base then
			anholonomic = true
			break
		end
	end
print('is anholonomic? '..tostring(anholonomic))

	-- for the sake of grid lengths, 
	-- I will need the basis and metric of the holonomic version as well
	local eHol
	if not anholonomic then
		eHol = e
	else
		eHol = Tensor('_u^I', function(a,I)
			return u[I]:diff(baseCoords[a])()
		end)
print'holonomic embedded:'
print(var'e''_u^I':eq(var'u''^I_,u'):eq(eHol'_u^I'()))
	end

	-- commutation coefficients
	local c = Tensor'_ab^c'
	if anholonomic then
--print'connection coefficients:'
--print(var'c''_uv^w' * var'e''_w','$=[ e_u, e_v ]$')
		for i,ui in ipairs(coords) do
			for j,uj in ipairs(coords) do
				local psi = var('\\psi', baseCoords)
				local diff = ui:applyDiff(uj:applyDiff(psi)) - uj:applyDiff(ui:applyDiff(psi))
				local diffEval = diff()
				if diffEval ~= const(0) then
--print('$[',ui.name,',',uj.name,'] =$',diff:eq(diffEval))
					diff = diff()
--print('factor division',diff)
					local dpsi = table.map(baseCoords, function(uk) return psi:diff(uk) end)
--print('dpsi', dpsi:unpack())
					local A,b = symmath.factorLinearSystem({diff}, dpsi)
					-- now extract psi:diff(uk)
					-- and divide by e_k to get the correct coefficient
					-- TODO this assumes that e_a is only a function of partial_a
					-- if e_a is a linear combination of e_a^b partial_b then you can work it out to find
					-- c_ab^d = (e^-1)_c^d (e_a^r e_b^c_,r - e_b^r e_a^c_,r)
					-- TODO put this somewhere else so everyone can use it
					assert(b[1][1] == const(0))
					for k,uk in ipairs(coords) do
						local coeff = (A[1][k] * dpsi[k] / uk:applyDiff(psi))()
						-- assert dphi is nowhere in coeff ...
						c[i][j][k] = coeff 
					end
				end
			end
		end
	end
print'commutation:'
print(var'c''_uv^w':eq(c'_uv^w'()))

	local g = (e'_u^I' * e'_v^J' * eta'_IJ')()
	--[[ TODO automatically do this ...
	g = g:map(function(expr)
		if symmath.op.pow.is(expr)
		and expr[2] == const(2)
		and symmath.cos.is(expr[1])
		then
			return 1 - symmath.sin(expr[1][1]:clone())^2
		end
	end)()
	--]]
print'metric:'
print(var'g''_uv':eq(var'e''_u^I' * var'e''_v^J' * var'\\eta''_IJ'):eq(g'_uv'()))

	local gHol
	if not anholonomic then
		gHol = g
	else
		gHol = (eHol'_u^I' * eHol'_v^J' * eta'_IJ')()
print'holonomic metric:'
print(var'g''_uv':eq(var'e''_u^I' * var'e''_v^J' * var'\\eta''_IJ'):eq(gHol'_uv'()))
	end	

	Tensor.metric(g)

	local GammaL = Tensor'_abc'
	GammaL['_abc'] = ((g'_ab,c' + g'_ac,b' - g'_bc,a' + c'_abc' + c'_acb' - c'_bca') / 2)()
print'1st kind Christoffel:'
print(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(GammaL'_abc'()))

	local Gamma = Tensor'^a_bc'
	Gamma['^a_bc'] = GammaL'^a_bc'()
print'connection:'
print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma'^a_bc'()))


	-- code generation

	local paramU = Tensor('^a', function(a) 
		return var('{v^'..a..'}')
	end)
	
	local toC = require 'symmath.tostring.C'
	local toC_coordArgs = table.map(baseCoords, function(coord, i)
		return {['{x^'..i..'}'] = coord}	-- 1-based
	end):append(range(dim):map(function(a)
		return {[paramU[a].name] = paramU[a]}
	end))
	local function compile(expr)
		local orig = expr	
		-- replace pow(x,2) with x*x
		expr = expr:map(function(x)
			if symmath.op.pow.is(x) 
			and const.is(x[2])
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
		
		local code = toC:compile(
			expr,
			toC_coordArgs
		):match'return (.*);'
	
--print'compiling...'
--print(orig)
--print'...to...'
--print(code)

		return code
	end

	-- uCode is used to project the grid for displaying
	self.uCode = range(dim):map(function(i) 
		local uCode = compile(u[i])
print('uCode['..i..'] = '..uCode)
		return uCode
	end)

	-- just giving up and manually writing this out
	
	local function cross(a,b)
		return table{
			(a[2]*b[3]-b[2]*a[3])(),
			(a[3]*b[1]-b[3]*a[1])(),
			(a[1]*b[2]-b[1]*a[2])(),
		}
	end

	-- extend 'e' to full R3 
	-- TODO should I do this from the start?
	-- just provide the full R3 coordinates, and make no 'eExt' struct
	local eExt = table()
	eExt[1] = range(3):map(function(i) return e[1][i] or const(0) end)
	if dim >= 2 then
		eExt[2] = range(3):map(function(i) return e[2][i] or const(0) end)
	else
		-- e2 is the maximal length orthogonal to e1
		local found = false
		for i =1,3 do
			local n = table{0,0,0}:map(function(x) return const(x) end)
			n[i] = 1
			eExt[2] = cross(n,eExt[1])
			local e2len = symmath.sqrt(eExt[2]:sum())()
			if e2len ~= const(0) then
				found = true
			end
		end
		if not found then
			error("how can you not find a perpendicular basis for e1="..eExt[1]:map(tostring):concat', ')
		end
	end
	if dim >= 3 then
		eExt[3] = range(3):map(function(i) return e[3][i] or const(0) end)
	else
		eExt[3] = cross(eExt[1],eExt[2])
	end

	self.eCode = eExt:map(function(ei,i) 
		return ei:map(function(eij,j)
			local eijCode = compile(eij) 
print('eCode['..i..']['..j..'] = ' .. tostring(eijCode))
			return eijCode 
		end)
	end)
	
--[=[ not being used
	local eHolLen = range(#eHol):map(function(i)
		return symmath.sqrt(
			range(#eHol):map(function(j)
				return eHol[i][j]^2
			end):sum()
		)()
	end)

	self.eHolLenCode = eHolLen:map(function(eiHolLen, i)
		local eiHolLenCode = compile(eiHolLen)
print('eHolLen['..i..'] = '..eiHolLenCode)
		return eiHolLenCode
	end)

	local eExtLen = eExt:map(function(ei,i)
		return symmath.sqrt(ei:map(function(x) return x^2 end):sum())()
	end)
	local eExtUnit = eExt:map(function(ei,i)
		return ei:map(function(eij) return (eij/eExtLen[i])() end)
	end)
	self.eUnitCode = eExtUnit:map(function(ei_unit,i) return ei_unit:map(compile) end)
print('eUnitCode = ', tolua(self.eUnitCode, {indent=true}))
--]=]

	local lowerExpr = paramU'_a'()
	self.lowerCodes = range(dim):map(function(i)
		local lowerCode = compile(lowerExpr[i])
print('lowerCode['..i..'] = '..lowerCode)
		return lowerCode
	end)

	local lenSqExpr = (paramU'^a' * paramU'_a')()
	self.uLenSqCode = compile(lenSqExpr)
print('uLenSqCodes = '..self.uLenSqCode)

	-- Conn^i_jk(x) v^j v^k
	local connExpr = (Gamma'^a_bc' * paramU'^b' * paramU'^c')()
	self.connCodes = range(dim):map(function(i)
		local conniCode = compile(connExpr[i])
print('connCode['..i..'] = '..conniCode)
		return conniCode
	end)

	-- dx is the change across the grid
	-- therefore it is based on the holonomic metric
	self.dxCodes = range(dim):map(function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'^b' * gHol'_ab')()
		local lenCode = compile((symmath.sqrt(lenSqExpr))())
print('dxCode['..i..'] = '..lenCode)
		return lenCode
	end)

	self.g = g
	self.gCode = range(dim):map(function(i)
		return range(i,dim):map(function(j)
			local gijCode = compile(self.g[i][j])
print('g['..i..']['..j..'] = '..gijCode)
			return gijCode, j
		end)
	end)

	local gU = Tensor('^ab', table.unpack((Matrix.inverse(g))))
	self.gU = gU
	self.gUCode = range(dim):map(function(i)
		return range(i,dim):map(function(j)
			local gUijCode = compile(self.gU[i][j])
print('gU['..i..']['..j..'] = '..gUijCode)
			return gUijCode, j
		end)
	end)
	
	local sqrt_gU = Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
	self.sqrt_gUCode = range(dim):map(function(i)
		return range(i,dim):map(function(j)
			local sqrt_gUijCode = compile(sqrt_gU[i][j])
print('sqrt(gU['..i..']['..j..']) = '..sqrt_gUijCode)
			return sqrt_gUijCode, j
		end)
	end)

	local volumeExpr = symmath.sqrt(symmath.Matrix.determinant(g))()
	self.volumeCode = compile(volumeExpr)
print('volumeCode = '..self.volumeCode)

end

return Geometry
