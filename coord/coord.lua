--[[

TODO rename all these to 'Coord'
because they are coordinate system properties, not geometry properties



This tells us the coordinate chart for our embedding in nD Cartesian (Euclidian?) geometry.
There are a few options on how to do this.

1) The Engineer way:
	Represent our eqn vector coordinates in Cartesian coordinates,
	use our CoordinateSystem to determine the cell volumes, areas, centers, normals, etc.,
	compute flux in Cartesian coordinates.
	
	This seems like the least change from a Cartesian grid, and should be easy to implement:
	All you have to change is the face normals and the dx values of the finite volume update 
	and change dx's in any operators of divergence (div B=0 magnetic monopole constraint) and Laplacian (del rho = 4 pi G for gravitation).
	
	This treats the problem as if it were made up a rigid grid, i.e. flattening off the curves in cylindrical and treating each cell like a polygon.
	For that reason, face center positions should probably be calculated as the average of vertices rather than the coordinate intermediate position,
	otherwise they won't represent the flat geometry and will introduce errors.
	This misses out on the perk of simulating problems whose components are purely rotational about an origin
	with much greater accuracy than a Cartesian grid.

	...in fact, this can be implemented as a 'grid' class on top of whatever curvilinear coordinates I choose --
	I should be able to use a cylindrical grid and cartesian coordaintes with no extra modifications to the 
	connections, or raising/lowering of coordinates -- only modifying my volume and side computations.

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

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local CoordinateSystem = class()

function CoordinateSystem:init(args)
	local symmath = require 'symmath'
	local const = symmath.Constant

	self.verbose = cmdline.coordVerbose or cmdline.verbose

	local dim = 3
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

	local baseCoords = table.mapi(coords, function(coord)
		return coord.base or coord
	end)
	--[[
	if self.verbose then
		print('coordinates:', table.unpack(coords))
		print('base coords:', table.unpack(baseCoords))
		print('embedding:', table.unpack(embedded))
	end
	--]]

	local eta = Tensor('_IJ', table.unpack(flatMetric)) 
	--[[
	if self.verbose then
		print'flat metric:'
		print(var'\\eta''_IJ':eq(eta'_IJ'()))
		print()
	end
	--]]

	local u = args.chart()
	if self.verbose then
		print'coordinate chart:'
		print(var'u''^I':eq(u'^I'()))
		print()
	end

	local e = Tensor'_u^I'
	e['_u^I'] = u'^I_,u'()
	if self.verbose then
		print'embedded:'
		print(var'e''_u^I':eq(var'u''^I_,u'):eq(e'_u^I'()))
		print()
	end

	local anholonomic
	for i=1,#coords do
		if coords[i].base then
			anholonomic = true
			break
		end
	end
	self.anholonomic = anholonomic
	if self.verbose then
		print('is anholonomic? '..tostring(anholonomic))
	end

	-- for the sake of grid lengths, 
	-- I will need the basis and metric of the holonomic version as well
	local eHol
	if not anholonomic then
		eHol = e
	else
		eHol = Tensor('_u^I', function(a,I)
			return u[I]:diff(baseCoords[a])()
		end)
		if self.verbose then
			print'holonomic embedded:'
			print(var'eHol''_u^I':eq(var'u''^I_,u'):eq(eHol'_u^I'()))
		end
	end

	-- commutation coefficients
	local c = Tensor'_ab^c'
	if anholonomic then
		if self.verbose then
			print'connection coefficients:'
			print(var'c''_uv^w' * var'e''_w','$=[ e_u, e_v ]$')
		end
		for i,ui in ipairs(coords) do
			for j,uj in ipairs(coords) do
				local psi = var('\\psi', baseCoords)
				local diff = ui:applyDiff(uj:applyDiff(psi)) - uj:applyDiff(ui:applyDiff(psi))
				local diffEval = diff()
				if diffEval ~= const(0) then
					if self.verbose then
						print('$[',ui.name,',',uj.name,'] =$',diff:eq(diffEval))
					end
					diff = diff()
					if self.verbose then
						print('factor division',diff)
					end
					local dpsi = table.mapi(baseCoords, function(uk) return psi:diff(uk) end)
					if self.verbose then
						print('dpsi', dpsi:unpack())
					end
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
	if self.verbose then
		print'commutation:'
		print(var'c''_uv^w':eq(c'_uv^w'()))
	end

	local g = (e'_u^I' * e'_v^J' * eta'_IJ')()

	if self.verbose then
		print'metric:'
		print(var'g''_uv':eq(var'e''_u^I' * var'e''_v^J' * var'\\eta''_IJ'):eq(g'_uv'()))
	end

	local gHol
	if not anholonomic then
		gHol = g
	else
		gHol = (eHol'_u^I' * eHol'_v^J' * eta'_IJ')()
		if self.verbose then
			print'holonomic metric:'
			print(var'gHol''_uv':eq(var'eHol''_u^I' * var'eHol''_v^J' * var'\\eta''_IJ'):eq(gHol'_uv'()))
		end
	end	

	Tensor.metric(g)

	local dg = Tensor'_cab'
	dg['_cab'] = g'_ab,c'()
self.dg = dg
	if self.verbose then
		print'metric partial:'
		print(var'g''_ab,c':eq(dg'_cab'()))
	end
	
	local Gamma_lll = Tensor'_abc'
	Gamma_lll['_abc'] = ((dg'_cab' + dg'_bac' - dg'_abc' + c'_abc' + c'_acb' - c'_bca') / 2)()
self.Gamma_lll = Gamma_lll	
	if self.verbose then
		print'1st kind Christoffel:'
		print(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(Gamma_lll'_abc'()))
	end

	local Gamma_ull = Tensor'^a_bc'
	Gamma_ull['^a_bc'] = Gamma_lll'^a_bc'()
self.Gamma_ull = Gamma_ull	
	if self.verbose then
		print'connection:'
		print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma_ull'^a_bc'()))
	end

	-- for anholonomic coordinates, we also need the holonomic connections
	--  for calculating the sqrt(det(g)) = grid volume
	if self.anholonomic then
		Tensor.metric(gHol)
		
		local GammaHol_lll = Tensor'_abc'
		GammaHol_lll['_abc'] = ((gHol'_ab,c' + gHol'_ac,b' - gHol'_bc,a') / 2)()
		if self.verbose then
			print'1st kind Christoffel of holonoic basis / Levi-Civita connection:'
			print(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'gHol''_ab,c' + var'gHol''_ac,b' - var'gHol''_bc,a')):eq(GammaHol_lll'_abc'()))
		end

		local GammaHol_ull = Tensor'^a_bc'
		GammaHol_ull['^a_bc'] = GammaHol_lll'^a_bc'()
		if self.verbose then
			print'holonomic / Levi-Civita connection:'
			print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(GammaHol_ull'^a_bc'()))
		end

		Tensor.metric(g)
	end

	local d2g = Tensor'_abcd'
	d2g['_cdab'] = dg'_cab,d'()

	-- code generation

	local function substCoords(code)
		code = code:gsub('{pt^(%d)}', function(i)
			return self.coords[i+0]
		end)
		code = code:gsub('{u^(%d)}', function(i)
			return 'u'..self.coords[i+0]
		end)
		code = code:gsub('{v^(%d)}', function(i)
			return 'v'..self.coords[i+0]
		end)
		return code
	end

	local paramU = Tensor('^a', function(a)
		return var('{u^'..a..'}')
	end)
	
	local paramV = Tensor('^a', function(a)
		return var('{v^'..a..'}')
	end)
	
	local toC = require 'symmath.tostring.C'
	local toC_coordArgs = table.mapi(baseCoords, function(coord, i)
		return {['{pt^'..i..'}'] = coord}	-- 1-based
	end):append(range(dim):mapi(function(a)
		return {[paramU[a].name] = paramU[a]}
	end)):append(range(dim):mapi(function(a)
		return {[paramV[a].name] = paramV[a]}
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
						return symmath.op.mul(range(value):mapi(function() 
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

		--[[
		if self.verbose then
			print('compiling...'..orig..'...to...'..code)
		end
		--]]

		return code
	end

	-- uCode is used to project the grid for displaying
	self.uCode = range(dim):mapi(function(i) 
		local uCode = compile(u[i])
		if cmdline.coordVerbose then
			if uCode ~= '0.' then
				print('uCode['..i..'] = '..substCoords(uCode))
			end
		end
		return uCode
	end)
	
	-- extend 'e' to full R3 
	-- TODO should I do this from the start?
	-- just provide the full R3 coordinates, and make no 'eExt' struct
	local eExt = range(dim):mapi(function(j)
		return range(dim):mapi(function(i) return e[j][i] or const(0) end)
	end)

	self.eCode = eExt:mapi(function(ei,i) 
		return ei:mapi(function(eij,j)
			local eijCode = compile(eij) 
			if cmdline.coordVerbose then
				if eijCode ~= '0.' then
					print('eCode['..i..']['..j..'] = '..substCoords(eijCode))
				end			
			end			
			return eijCode 
		end)
	end)
	
--[=[ not being used
	local eHolLen = range(#eHol):mapi(function(i)
		return symmath.sqrt(
			range(#eHol):mapi(function(j)
				return eHol[i][j]^2
			end):sum()
		)()
	end)
	
	self.eHolLenCode = eHolLen:mapi(function(eiHolLen, i)
		local eiHolLenCode = compile(eiHolLen)
		if cmdline.coordVerbose then
			print('eHolLen['..i..'] = '..substCoords(eiHolLenCode))
		end
		return eiHolLenCode
	end)
	
	local eExtLen = eExt:mapi(function(ei,i)
		return symmath.sqrt(ei:mapi(function(x) return x^2 end):sum())()
	end)
	local eExtUnit = eExt:mapi(function(ei,i)
		return ei:mapi(function(eij) return (eij/eExtLen[i])() end)
	end)
	self.eUnitCode = eExtUnit:mapi(function(ei_unit,i) return ei_unit:mapi(compile) end)
	if cmdline.coordVerbose then
		print('eUnitCode = ', tolua(self.eUnitCode, {indent=true}))
	end
--]=]

	-- v^k -> v_k
	local lowerExpr = paramU'_a'()
	self.lowerCodes = range(dim):mapi(function(i)
		local lowerCode = compile(lowerExpr[i])
		if cmdline.coordVerbose then
			print('lowerCode['..i..'] = '..substCoords(lowerCode))
		end
		return lowerCode
	end)

	-- v^k -> v_k
	local raiseExpr = paramU'_a'()
	self.raiseCodes = range(dim):mapi(function(i)
		local raiseCode = compile(raiseExpr[i])
		if cmdline.coordVerbose then
			print('raiseCode['..i..'] = '..substCoords(raiseCode))
		end
		return raiseCode
	end)

	-- v^k v_k
	local lenSqExpr = (paramU'^a' * paramU'_a')()
	self.uLenSqCode = compile(lenSqExpr)
	if cmdline.coordVerbose then
		print('uLenSqCodes = '..substCoords(self.uLenSqCode))
	end

	self.dg_lll_codes = range(dim):mapi(function(k)
		return range(dim):mapi(function(i)
			return range(dim):map(function(j)
				local code = compile(dg[k][i][j])
				if cmdline.coordVerbose then
					if code ~= '0.' then
						print('dg['..k..i..j..'] = '..substCoords(code))
					end
				end
				return code
			end)
		end)
	end)

	self.d2g_llll_codes = range(dim):mapi(function(k)
		return range(dim):mapi(function(l)
			return range(dim):mapi(function(i)
				return range(dim):map(function(j)
					local code = compile(d2g[k][l][i][j])
					if cmdline.coordVerbose then
						if code ~= '0.' then
							print('d2g['..k..l..i..j..'] = '..substCoords(code))
						end
					end
					return code
				end)
			end)
		end)
	end)

	self.conn_lll_codes = range(dim):mapi(function(i)
		return range(dim):mapi(function(j)
			return range(dim):mapi(function(k)
				local code = compile(Gamma_lll[i][j][k])
				if cmdline.coordVerbose then
					if code ~= '0.' then
						print('connCode_lll['..i..j..k..'] = '..substCoords(code))
					end
				end
				return code
			end)
		end)
	end)

	self.conn_ull_codes = range(dim):mapi(function(i)
		return range(dim):mapi(function(j)
			return range(dim):mapi(function(k)
				local code = compile(Gamma_ull[i][j][k])
				if cmdline.coordVerbose then
					if code ~= '0.' then
						print('connCode_ull['..i..j..k..'] = '..substCoords(code))
					end
				end
				return code
			end)
		end)
	end)
	
	-- Conn^i_jk(x) u^j v^k
	local connExpr = (Gamma_ull'^a_bc' * paramU'^b' * paramV'^c')()
	self.connApply23Codes = range(dim):mapi(function(i)
		local conniCode = compile(connExpr[i])
		if cmdline.coordVerbose then
			if conniCode ~= '0.' then
				print('connApply23Code['..i..'] = '..substCoords(conniCode))
			end		
		end		
		return conniCode
	end)

	-- u^j v^k Conn_jk^i(x)
	local connLastExpr = (paramU'^b' * paramV'^c' * Gamma_ull'_bc^a')()
	self.connApply12Codes = range(dim):mapi(function(i)
		local code = compile(connLastExpr[i])
		if cmdline.coordVerbose then
			if code ~= '0.' then
				print('connApply12Code['..i..'] = '..substCoords(code))
			end		
		end		
		return code
	end)
	
	-- Conn^i = Conn^i_jk g^jk
	self.tr23_conn_u_expr = (Gamma_ull'^a_b^b')()
	self.tr23_conn_u_codes = range(dim):mapi(function(i)
		local code = compile(self.tr23_conn_u_expr[i])
		if cmdline.coordVerbose then
			if code ~= '0.' then
				print('tr23_conn_u_code['..i..'] = '..substCoords(code))
			end		
		end		
		return code
	end)

	-- sqrt(g)_,i / sqrt(g) = Conn^j_ij
	local tr13_conn_l_codes = (Gamma_ull'^b_ab')()
	self.tr13_conn_l_codes = range(dim):mapi(function(i)
		local code = compile(tr13_conn_l_codes[i])
		if cmdline.coordVerbose then
			if code ~= '0.' then
				print('tr13_conn_l_code['..i..'] = '..substCoords(code))
			end		
		end		
		return code
	end)

	self.lenExprs = range(dim):mapi(function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'^b' * gHol'_ab')()
		local lenExpr = symmath.sqrt(lenSqExpr)()
		return lenExpr
	end)

	-- dx is the change across the grid
	-- therefore it is based on the holonomic metric
	self.dxCodes = range(dim):mapi(function(i)
		local lenCode = compile(self.lenExprs[i])
		if cmdline.coordVerbose then
			print('dxCode['..i..'] = '..substCoords(lenCode))
		end
		return lenCode
	end)

	local areaExprs = range(dim):mapi(function(i)
		local area = const(1)
		for j=1,dim do
			if j ~= i then
				area = area * self.lenExprs[j]
			end
		end
		area = area()
		return area
	end)

	-- area of the side in each direction
	self.areaCodes = range(dim):mapi(function(i)
		local areaCode = compile(areaExprs[i])
		if cmdline.coordVerbose then
			print('areaCode['..i..'] = '..substCoords(areaCode))
		end
		return areaCode
	end)

	self.g = g
	self.gCode = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local gijCode = compile(self.g[i][j])
			if cmdline.coordVerbose then
				if gijCode ~= '0.' then
					print('g['..i..']['..j..'] = '..substCoords(gijCode))
				end			
			end			
			return gijCode, j
		end)
	end)

	local gU = Tensor('^ab', table.unpack((Matrix.inverse(g))))
	self.gU = gU
	self.gUCode = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local gUijCode = compile(self.gU[i][j])
			if cmdline.coordVerbose then
				if gUijCode ~= '0.' then
					print('gU['..i..']['..j..'] = '..substCoords(gUijCode))
				end			
			end			
			return gUijCode, j
		end)
	end)
	
	local sqrt_gU = Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
	self.sqrt_gUCode = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local sqrt_gUijCode = compile(sqrt_gU[i][j])
			if cmdline.coordVerbose then
				if sqrt_gUijCode ~= '0.' then
					print('sqrt(gU['..i..']['..j..']) = '..substCoords(sqrt_gUijCode))
				end			
			end			
			return sqrt_gUijCode, j
		end)
	end)

	local det_g_expr = symmath.Matrix.determinant(gHol)
	self.det_g_code = compile(det_g_expr)
	if cmdline.coordVerbose then
		print('det_g_code = '..substCoords(self.det_g_code))
	end

	local volumeExpr = symmath.sqrt(det_g_expr)()
	self.volumeCode = compile(volumeExpr)
	if cmdline.coordVerbose then
		print('volumeCode = '..substCoords(self.volumeCode))
	end




	-- Now, for num rel in spherical, I need to transform the metric to remove the singularities ...
	-- It is funny to include this alongside the anholonomic stuff.  It seems you would need one or the other but not both.
	
	local J = Tensor('^a_b', function(i,j)
		return i == j and (1 / self.lenExprs[i]) or 0
	end)
	if self.verbose then
		print'rescaling transform:'
		print(J)
	end
	local Tensor = symmath.Tensor
	local JInv = Tensor('_a^b', table.unpack((Matrix.inverse(J))))
	if self.verbose then
		print'rescaling inverse transform:'
		print(JInv)
	end
	local gJ_ll = (self.g'_ab' * J'^a_u' * J'^b_v')()
	if self.verbose then
		print'g_ij, in rescaled coordinates:'
		print(gJ_ll)
	end
	local gJ_uu = (self.gU'^ab' * JInv'_a^u' * JInv'_b^v')()
	if self.verbose then
		print'g^ij, in rescaled coordinates:'
		print(gJ_uu)
	end
	self.gJ_uu_code = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local gJ_uu_ij_code = compile(gJ_uu[i][j])
			if cmdline.coordVerbose then
				if gJ_uu_ij_code ~= '0.' then
					print('gJ_uu['..i..']['..j..'] = '..substCoords(gJ_uu_ij_code))
				end			
			end			
			return gJ_uu_ij_code, j
		end)
	end)

	local sqrt_gJ_uu_expr = Tensor('^ab', function(a,b) return symmath.sqrt(gJ_uu[a][b])() end)
	self.sqrt_gJ_uu_code = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local sqrt_gJ_uu_ij_code = compile(sqrt_gJ_uu_expr[i][j])
			if cmdline.coordVerbose then
				if sqrt_gJ_uu_ij_code ~= '0.' then
					print('sqrt(gJ_uu['..i..']['..j..']) = '..substCoords(sqrt_gJ_uu_ij_code))
				end			
			end			
			return sqrt_gJ_uu_ij_code, j
		end)
	end)
	
	local dgJ_lll = (self.dg'_abc' * J'^a_u' * J'^b_v' * J'^c_w')()
	if self.verbose then
		print'g_ij,k in rescaled coordinates:'
		print(dgJ_lll)
	end
	local connJ_lll = (self.Gamma_lll'_abc' * J'^a_u' * J'^b_v' * J'^c_w')()
	if self.verbose then
		print'Gamma_ijk in rescaled coordinates:'
		print(connJ_lll)
	end
	local connJ_ull = (self.Gamma_ull'^a_bc' * JInv'_a^u' * J'^b_v' * J'^c_w')()
	if self.verbose then
		print'Gamma^i_jk in rescaled coordinates:'
		print(connJ_ull)
	end
end



local template = require 'template'
local clnumber = require 'cl.obj.number'

local xs = table{'x', 'y', 'z'}

local function convertParams(code)
	code = code:gsub('{pt^(%d)}', function(i)
		return 'pt.'..xs[i+0]
	end)
	code = code:gsub('{u^(%d)}', function(i)
		return 'u.'..xs[i+0]
	end)
	code = code:gsub('{v^(%d)}', function(i)
		return 'v.'..xs[i+0]
	end)
	return code
end

local function getCode_real3_to_real(name, code)
	return template([[
real <?=name?>(real3 pt) {
	return <?=code?>;
}]], {
		name = name,
		code = convertParams(code),
	})
end

-- f(x) where x is a point in the coordinate chart
local function getCode_real3_to_real3(name, exprs)
	return template([[
real3 <?=name?>(real3 pt) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] and convertParams(exprs[i]) or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end

-- f(v,x) where x is a point on the coordinate chart and v is most likely a tensor
local function getCode_real3_real3_to_real(name, expr)
	return template([[
real <?=name?>(real3 u, real3 pt) {
	return <?=convertParams(expr)?>;
}]], {
		name = name,
		expr = expr,
		convertParams = convertParams,
	})
end

local function getCode_real3_real3_to_real3(name, exprs)
	return template([[
real3 <?=name?>(real3 u, real3 pt) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] and convertParams(exprs[i]) or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end

local function getCode_real3_real3_real3_to_real3(name, exprs)
	return template([[
real3 <?=name?>(real3 u, real3 v, real3 pt) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] and convertParams(exprs[i]) or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end


local function getCode_real3_to_sym3(name, exprs)
	return template([[
sym3 <?=name?>(real3 pt) {
	return (sym3){
<? for i=1,3 do
	for j=i,3 do
?>		.<?=xs[i]..xs[j]?> = <?=exprs[i] and exprs[i][j] 
			and convertParams(exprs[i][j]) or '0.'?>,
<?	end
end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end

-- symmetric on 2nd & 3rd indexes
local function getCode_real3_to_3sym3(name, exprs)
	return template([[
_3sym3 <?=name?>(real3 pt) {
	return (_3sym3){
<? 
for i=1,3 do
?>	.<?=xs[i]?> = {
<?	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>		.<?=xjk?> = <?=exprs[i] and exprs[i][j] and exprs[i][j][k]
			and convertParams(exprs[i][j][k]) or '0.'?>,
<?	end	
?>	},
<?
end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
		convertParams = convertParams,
		symNames = symNames,
		from6to3x3 = from6to3x3,
	})
end

-- symmetric on 1st & 2nd and on 3rd & 4th
local function getCode_real3_to_sym3sym3(name, exprs)
	return template([[
sym3sym3 <?=name?>(real3 pt) {
	return (sym3sym3){
<? for kl,xkl in ipairs(symNames) do
?>		.<?=xkl?> = (sym3){
<?	for ij,xij in ipairs(symNames) do
?>			.<?=xij?> = <?=exprs[i] and exprs[i][j] and exprs[i][j][k] and exprs[i][j][k][l]
				and convertParams(exprs[i][j][k][l]) or '0.'?>,
<?	end
?>		},
<? end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
		convertParams = convertParams,
		symNames = symNames,
	})
end

--[[
ok standardizing these macros ...
we have a few manifolds to deal with ...
1) grid coordiantes, I = 0..size-1
2) manifold coordinates, X = grid coordinates rescaled to fit between solver->mins and solver->maxs  
3) result coordinates, Y = cartesian output of the mapping of the coordinate system

coord_dx[_for_coord]<?=side?>(pt) gives the length of the dx of the coordinate in 'side' direction, for 'pt' in coordinates
cell_dx[_for_coord]<?=side?>(pt) gives the length of the dx of the cell, which is just coord_dx times the grid_dx
solver->grid_dx is the size of each solver cell, in coordintes 

but then we have these cell_x[_for_grid]<?=side?> functions which take an input of 'i'

should I add these _for_coord _for_grid suffixes to specify what manfiold system the input parameter is? 

--]]
function CoordinateSystem:getCode(solver)
	self.solver = solver
	local dim = solver.dim
	local lines = table()

	local function convertInputFromGrid(code, inputVarName)
		for j=1,3 do
			code = code:gsub(
				'{pt^'..j..'}',
				'cell_x'..(j-1)..'('..inputVarName..'.'..xs[j]..')')
		end
		return code
	end

	local function convertInputFromCoord(code, inputVarName)
		for j=1,3 do
			code = code:gsub('{pt^'..j..'}', inputVarName..'.'..xs[j])
		end
		return code
	end
	
	-- dx0, ...
	-- this is the change in cartesian wrt the change in grid
	lines:append(range(dim):mapi(function(i)
		local code = convertInputFromCoord(self.dxCodes[i], 'pt')
		return '#define coord_dx'..(i-1)..'(pt) ('..code..')'
	end))

	-- area0, ...
	lines:append(range(dim):mapi(function(i)
		local code = convertInputFromCoord(self.areaCodes[i], 'pt')
		return '#define coord_area'..(i-1)..'(pt) ('..code..')'
	end))

	lines:insert'\n'

	lines:append(range(dim):mapi(function(i)
		return '#define cell_dx'..(i-1)..'(pt) (coord_dx'..(i-1)..'(pt) * solver->grid_dx.s'..(i-1)..')'
	end))
	
	lines:append(range(dim):mapi(function(i)
		local prod
		local code = '#define cell_area'..(i-1)..'(pt) (coord_area'..(i-1)..'(pt)'
		for j=1,dim do
			if j ~= i then
				code = code .. ' * solver->grid_dx.s'..(j-1)
			end
		end
		code = code .. ')'
		return code
	end))
	
	lines:insert'\n'

	-- metric determinant = volume^2
	local det_g_code = '(' .. self.det_g_code .. ')'
	lines:insert(getCode_real3_to_real('coord_det_g', det_g_code))

	-- volume
	local volumeCode = '(' .. self.volumeCode .. ')'
	lines:insert(getCode_real3_to_real('coord_volume', volumeCode))
	
	-- coord len code: l(v) = v^i v^j g_ij
	lines:append{
		getCode_real3_real3_to_real('coordLenSq', self.uLenSqCode),
		[[
real coordLen(real3 r, real3 pt) {
	return sqrt(coordLenSq(r, pt));
}]],
	}

	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply23', self.connApply23Codes))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply12', self.connApply12Codes))
	lines:insert(getCode_real3_to_real3('coord_conn_trace23', self.tr23_conn_u_codes))
	lines:insert(getCode_real3_to_real3('coord_conn_trace13', self.tr13_conn_l_codes))
	lines:insert(getCode_real3_to_3sym3('coord_dg_lll', self.dg_lll_codes))
	lines:insert(getCode_real3_to_sym3sym3('coord_d2g_llll', self.d2g_llll_codes))
	lines:insert(getCode_real3_to_3sym3('coord_conn_lll', self.conn_lll_codes))
	lines:insert(getCode_real3_to_3sym3('coord_conn_ull', self.conn_ull_codes))

	--[[
	for i=0,dim-1 do
		lines:insert(getCode_real3_to_real('coordHolBasisLen'..i, self.eHolLenCode[i+1]))
	end
	--]]

	for i,eiCode in ipairs(self.eCode) do
		lines:insert(getCode_real3_to_real3('coordBasis'..(i-1), eiCode))
	end

	lines:insert(getCode_real3_real3_to_real3('coord_lower', self.lowerCodes))
	lines:insert(getCode_real3_real3_to_real3('coord_raise', self.raiseCodes))

	do
		local function addSym3Components(name, codes)
			for i=1,3 do
				for j=i,3 do
					local code = (codes[i] and codes[i][j] and convertParams(codes[i][j]) or clnumber(i==j and 1 or 0))
					lines:insert('#define '..name..(i-1)..(j-1)..'(pt) '..code)
					if i ~= j then
						lines:insert('#define '..name..(j-1)..(i-1)..'(pt) '..code)
					end
				end
			end
		end
		
		addSym3Components('coord_g_ll', self.gCode)
		lines:insert(getCode_real3_to_sym3('coord_g_ll', self.gCode))
		
		addSym3Components('coord_g_uu', self.gUCode)
		lines:insert(getCode_real3_to_sym3('coord_g_uu', self.gUCode))
		
		addSym3Components('coord_sqrt_g_uu', self.sqrt_gUCode)
	
		
		-- gJ = g times jacobian = g times rescaling
		-- TODO better naming for this
		
		addSym3Components('coordRescaled_g_uu', self.gJ_uu_code)
		
		addSym3Components('coordRescaled_sqrt_g_uu', self.sqrt_gJ_uu_code)
	
		-- coord rescaling
		lines:insert(template([[
<? local dim = solver.dim ?>

#if 0	//rescaling, used for bssn finite-difference, but I am tempted to try it with other coordinate systems with singularities

//apply this to lower indexes to convert from coordinate metric to better metric
//apply this to upper indexes to convert from better metric to coordinate metric
real3 real3_rescaleFromCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x / coord_dx0(x),
		.y = v.y<? if dim > 1 then ?> / coord_dx1(x)<? end ?>,
		.z = v.z<? if dim > 2 then ?> / coord_dx2(x)<? end ?>,
	};
}
#define real3_rescaleToCoord_u real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y<? if dim > 1 then ?> * coord_dx1(x)<? end ?>,
		.z = v.z<? if dim > 2 then ?> * coord_dx2(x)<? end ?>,
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_l

sym3 sym3_rescaleFromCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (<?
	if i <= dim then ?>coord_dx<?=i-1?>(x)<? else ?>1.<? end
			?> * <?
	if j <= dim then ?>coord_dx<?=j-1?>(x)<? else ?>1.<? end
			?>),
<? end
?>	};
}
#define sym3_rescaleToCoord_uu sym3_rescaleFromCoord_ll

sym3 sym3_rescaleToCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (<?
	if i <= dim then ?>coord_dx<?=i-1?>(x)<? else ?>1.<? end
			?> * <?
	if j <= dim then ?>coord_dx<?=j-1?>(x)<? else ?>1.<? end
			?>),
<? end
?>	};
}
#define sym3_rescaleFromCoord_uu sym3_rescaleToCoord_ll

_3sym3 _3sym3_rescaleFromCoord_lll(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_rescaleFromCoord_ll(a.<?=xi?>, x),
<? end
?>	};
}
#define _3sym3_rescaleToCoord_uuu _3sym3_rescaleFromCoord_lll

_3sym3 _3sym3_rescaleToCoord_lll(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_rescaleToCoord_ll(a.<?=xi?>, x),
<? end
?>	};
}
#define _3sym3_rescaleFromCoord_uuu _3sym3_rescaleToCoord_lll

sym3sym3 sym3sym3_rescaleFromCoord_lll(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = sym3_rescaleFromCoord_ll(a.<?=xij?>, x),
<? end
?>	};
}
#define sym3sym3_rescaleToCoord_uuuu sym3sym3_rescaleFromCoord_llll

sym3sym3 sym3sym3_rescaleToCoord_llll(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
?>		.<?=xij?> = sym3_rescaleToCoord_ll(a.<?=xij?>, x),
<? end
?>	};
}
#define sym3sym3_rescaleFromCoord_uuuu sym3sym3_rescaleToCoord_llll

#else	//debugging -- turning it off

#define real3_rescaleFromCoord_l(a,x) a
#define real3_rescaleToCoord_u(a,x) a
#define real3_rescaleToCoord_l(a,x) a
#define real3_rescaleFromCoord_u(a,x) a
#define sym3_rescaleFromCoord_ll(a,x) a
#define sym3_rescaleToCoord_uu(a,x) a
#define sym3_rescaleToCoord_ll(a,x) a
#define sym3_rescaleFromCoord_uu(a,x) a
#define _3sym3_rescaleFromCoord_lll(a,x) a
#define _3sym3_rescaleToCoord_uuu(a,x) a
#define _3sym3_rescaleToCoord_lll(a,x) a
#define _3sym3_rescaleFromCoord_uuu(a,x) a
#define sym3sym3_rescaleFromCoord_lll(a,x) a
#define sym3sym3_rescaleToCoord_uuuu(a,x) a
#define sym3sym3_rescaleToCoord_llll(a,x) a
#define sym3sym3_rescaleFromCoord_uuuu (a,x) a

#endif

]], 	{
			solver = solver,
			xNames = xNames,
			symNames = symNames,
			from6to3x3 = from6to3x3 ,
		}))
	end

	lines:insert(template([[

//converts a vector from cartesian coordinates to grid coordinates
//by projecting the vector into the grid basis vectors 
//at x, which is in grid coordinates
real3 cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3_zero;
	<? for i=0,solver.dim-1 do ?>{
		real3 e = coordBasis<?=i?>(pt);
<? if coord.anholonomic then	-- anholonomic normalized
?>		uCoord.s<?=i?> = real3_dot(e, u); // / real3_len(e);
<? else		-- holonomic
?>		real uei = real3_dot(e, u) / real3_lenSq(e);
<? end		
?>		uCoord.s<?=i?> = uei;
		//subtract off this basis component from u
		u = real3_sub(u, real3_real_mul(e, uei));
	}<? end ?>
	//add whatever's left of u
	uCoord = real3_add(uCoord, u);
	return uCoord;
}

//converts a vector from cartesian to grid
//by projecting it onto the basis ... ?
real3 cartesianFromCoord(real3 u, real3 pt) {
	real3 uGrid = real3_zero;
	<? for i=0,solver.dim-1 do ?>{
		real3 e = coordBasis<?=i?>(pt);
		uGrid = real3_add(uGrid, real3_real_mul(e, u.s<?=i?>));
	}<? end ?>
	return uGrid;
}

]], {
		solver = solver,
		coord = self,
	}))

	lines:insert(self:getCoordMapCode())

--print(require 'template.showcode'(lines:concat'\n'))

	return lines:concat'\n'
end

function CoordinateSystem:getCoordMapCode()
	return table{
		getCode_real3_to_real3('coordMap', range(3):mapi(function(i)
			return self.uCode[i] or '{pt^'..i..'}'
		end)),
	}:concat'\n'
end

function CoordinateSystem:getCoordMapGLSLCode()
	return (self:getCoordMapCode()
		:gsub('inline%S*', '')
		:gsub('_real3', 'vec3')	-- real3 constructor
		:gsub('real3', 'vec3')	-- real3 type
	)
end

-- until I get inverses on trig functions working better, I'll have this manually specified
function CoordinateSystem:getCoordMapInvGLSLCode()
	return [[
vec3 coordMapInv(vec3 x) { return x; }
]]
end

return CoordinateSystem
