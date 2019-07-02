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
	Represent our eqn vector components in a holonomic basis.  No normalization.
	This would mean keeping that g^ij in front of the pressure term, and therefore adjusting the eigen-decomposition of the Euler fluid equations.
	It would also mean performing our dx's in coordinate space.  No anholonomic readjustments.

	This is like the Physicist way except *without* extra normalization for the holonomic-anholonomic basis exchanges.
	It is also made to work with a metric, and therefore works easiest when the metric is a dynamic variable (right?  where will the change-in-metric stuff have to go?)
	Turns out this is similar to the conservative form for relativistic solvers mentioned in papers by Font and company.

(and after reading the coordinate-invariant BSSN papers)
4) the Mathematician+Physicist way:
	Represent all equations as tensors.  Use a holonomic basis.  Use a final change of coordinates that normalizes the basis vector lengths.
	You get the best of both worlds: coordinates invariant to transform, and unit length basis vectors.
	But how about the partial spatial derivatives -- wouldn't those also incur an extra term of the partial of the extra change-of-coordinate system?

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
local template = require 'template'
local clnumber = require 'cl.obj.number'

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local CoordinateSystem = class()

--[[
args:
	anholonomic = use anholonomic coordinate system (normalized basis vectors)
		this causes the Coord object to use commutation coefficients, and therefore will have asymmetric connections
--]]
function CoordinateSystem:init(args)
	local symmath = require 'symmath'
	local const = symmath.Constant

	self.verbose = cmdline.coordVerbose or cmdline.verbose
	self.anholonomic = args.anholonomic

	local dim = 3
	local var = symmath.var
	local vars = symmath.vars
	local Matrix = symmath.Matrix
	local Tensor = symmath.Tensor

	local eHolToE = self.eHolToE
	if not eHolToE then
		eHolToE = Matrix.identity(3)
	end

	local baseCoords = self.baseCoords
	if not self.anholonomic then
		self.coords = table(baseCoords)
	else
		
		-- TODO this is why symmath needs CAS function objects
		-- and instead of overriding :applyDiff, just make operators a CAS function object
		local nonCoords = table()
		local nonCoordLinExpr = (eHolToE * baseCoords)()
		for i=1,3 do
			local baseCoord = baseCoords[i]
			-- the non-coordinate = the coordinate, so use the original variable 
			if nonCoordLinExpr[i] == baseCoord then
			-- the non-coordinate ~= the coordinate, so make a new non-coord var and modify its 'applyDiff' function
				nonCoords[i] = baseCoord
			else
				local nonCoord = symmath.var(baseCoord.name..'Hat')
				nonCoord.base = baseCoord
				function nonCoord:applyDiff(x)
					local xPartial = symmath.Matrix:lambda({dim, 1}, function(j,_)
						return x:diff(baseCoords[j])
					end)
					local result = Matrix:lambda({1,dim}, function(_,j) 
						assert(not symmath.Array.is(eHolToE[i][j]))
						return eHolToE[i][j] 
					end) 
						* xPartial
					result = result()
					assert(symmath.Matrix.is(result))
					result = result[1][1]
					assert(symmath.Expression.is(result) and not symmath.Array.is(result))
					return result
				end
				nonCoords[i] = nonCoord
			end
		end
		self.coords = nonCoords
	end
	local coords = self.coords


	local flatMetric = Matrix:lambda({dim, dim}, function(i,j) return i==j and 1 or 0 end)
	local embedded = self.embedded

	Tensor.coords{
		{variables=coords},
		{variables=embedded, symbols='IJKLMN', metric=flatMetric},
	}

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

	local u = self.chart()
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

	local isItReallyAnholonomic
	for i=1,#coords do
		if coords[i].base then
			isItReallyAnholonomic = true
			break
		end
	end
	assert(isItReallyAnholonomic == self.anholonomic)
	local anholonomic = self.anholonomic
	
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
			print(var'{e_{iHol}}^i':eq(eHolToE))
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
	Tensor.metric(g)

	
	-- code generation

	-- for output
	local function substCoords(code)
		code = code:gsub('{pt^(%d)}', function(i)
			return self.baseCoords[i+0].name
		end)
		code = code:gsub('{u^(%d)}', function(i)
			return 'u'..self.baseCoords[i+0].name
		end)
		code = code:gsub('{v^(%d)}', function(i)
			return 'v'..self.baseCoords[i+0].name
		end)
		code = code:gsub('{w^(%d)}', function(i)
			return 'w'..self.baseCoords[i+0].name
		end)
		return code
	end

	local paramU = Tensor('^a', function(a)
		return var('{u^'..a..'}')
	end)
	
	local paramV = Tensor('^a', function(a)
		return var('{v^'..a..'}')
	end)

	local paramW = Tensor('^a', function(a)
		return var('{w^'..a..'}')
	end)
	
	local toC = require 'symmath.tostring.C'
	local toC_coordArgs = table.mapi(baseCoords, function(coord, i)
		return {['{pt^'..i..'}'] = coord}	-- 1-based
	end):append(range(dim):mapi(function(a)
		return {[paramU[a].name] = paramU[a]}
	end)):append(range(dim):mapi(function(a)
		return {[paramV[a].name] = paramV[a]}
	end)):append(range(dim):mapi(function(a)
		return {[paramW[a].name] = paramW[a]}
	end))
	local function compile(expr, extraArgs)
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
	
		local args = toC_coordArgs
		if extraArgs then
			args = table(args):append(extraArgs)
		end

		local code = toC:compile(expr, args):match'return (.*);'

		--[[
		if self.verbose then
			print('compiling...'..orig..'...to...'..code)
		end
		--]]

		return code
	end

	local function printNonZero(name, code)
		if not self.verbose then return end
		local codetype = type(code)
		if codetype == 'string' then
			if code ~= '0.' then
				print(name..' = '..substCoords(code))
			end
		elseif codetype == 'table' then
			for i=1,#code do
				printNonZero(name..'['..i..']', code[i])
			end
		else
			error("can't print code "..('%q'):format(name).." of type "..codetype)
		end
	end
	local function printNonZeroField(field)
		printNonZero(field, self[field])
	end

	--compile a tensor of expressions to a nested table of codes
	local function compileTensor(expr, extraArgs)
		if symmath.Array.is(expr) then
			return table.mapi(expr, function(expri) 
				return compileTensor(expri, extraArgs)
			end)
		elseif symmath.Expression.is(expr) then
			return compile(expr, extraArgs)
		elseif type(expr) == 'number' then
			return clnumber(expr)
		else
			error("don't know how to compile "
				..symmath.Verbose(expr))
				--..require 'ext.tolua'(expr))
		end
	end
	local function compileTensorField(field, expr, extraArgs)
		self[field] = compileTensor(expr, extraArgs)
		printNonZeroField(field)
	end

	-- uCode is used to project the grid for displaying
	compileTensorField('uCode', u)

	-- extend 'e' to full R3 
	-- TODO should I do this from the start?
	-- just provide the full R3 coordinates, and make no 'eExt' struct?
	local eExt = symmath.Array:lambda({dim, dim}, function(i,j)
		return e[i][j] or const(0)
	end)
	compileTensorField('eCode', eExt)
	
--[=[ not being used
	local eHolLen = range(#eHol):mapi(function(i)
		return symmath.sqrt(
			range(#eHol):mapi(function(j)
				return eHol[i][j]^2
			end):sum()
		)()
	end)
	compileTensorField('eHolLenCode', eHolLen)

	local eExtLen = eExt:mapi(function(ei,i)
		return symmath.sqrt(ei:mapi(function(x) return x^2 end):sum())()
	end)
	local eExtUnit = eExt:mapi(function(ei,i)
		return ei:mapi(function(eij) return (eij/eExtLen[i])() end)
	end)
	compileTensorField('eUnitCode', eExtUnit)
--]=]

	-- v^k -> v_k
	local lowerExpr = paramU'_a'()
	compileTensorField('lowerCodes', lowerExpr)

	-- v^k -> v_k
	local raiseExpr = paramU'_a'()
	compileTensorField('raiseCodes', raiseExpr)

	-- v^k v_k
	local lenSqExpr = (paramU'^a' * paramU'_a')()
	compileTensorField('uLenSqCode', lenSqExpr)

	-- c_ab^b
	local tr23_c = c'_ab^b'()
	if self.verbose then
		print(var'c''_ab^b':eq(tr23_c))
	end
	compileTensorField('tr23_cCode', tr23_c)

	local dg = Tensor'_cab'
	dg['_cab'] = g'_ab,c'()
self.dg = dg
	if self.verbose then
		print'metric partial:'
		print(var'g''_ab,c':eq(dg'_cab'()))
	end
	compileTensorField('dg_lll_codes', dg)

	local d2g = Tensor'_abcd'
	d2g['_cdab'] = dg'_cab,d'()
	compileTensorField('d2g_llll_codes', d2g)

	local Gamma_lll = Tensor'_abc'
	Gamma_lll['_abc'] = ((dg'_cab' + dg'_bac' - dg'_abc' + c'_abc' + c'_acb' - c'_bca') / 2)()
	if self.verbose then
		print'1st kind Christoffel:'
		print(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(Gamma_lll'_abc'()))
	end
	compileTensorField('conn_lll_codes', Gamma_lll)

	local Gamma_ull = Tensor'^a_bc'
	Gamma_ull['^a_bc'] = Gamma_lll'^a_bc'()
	if self.verbose then
		print'connection:'
		print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma_ull'^a_bc'()))
	end
	compileTensorField('conn_ull_codes', Gamma_ull)

	local Gamma_ulll = Tensor'^a_bcd'
	Gamma_ulll['^a_bcd'] = Gamma_ull'^a_bc,d'()
	if self.verbose then
		print'connection:'
		print(var'\\Gamma''^a_bc,d':eq(Gamma_ulll'^a_bcd'()))
	end
	compileTensorField('partial_conn_ulll_codes', Gamma_ulll)

	-- u^i v^j b^k Conn_ijk(x)
	local connExpr = (paramU'^a' * paramV'^b' * paramW'^c' * Gamma_lll'_abc')()
	compileTensorField('connApply123Code', connExpr)

	-- u^j v^k Conn_jk^i(x)
	local connExpr = (paramU'^b' * paramV'^c' * Gamma_lll'_bc^a')()
	compileTensorField('connApply12Codes', connExpr)

	-- u^j v^k Conn_j^i_k(x)
	local connExpr = (paramU'^b' * paramV'^c' * Gamma_lll'_b^a_c')()
	compileTensorField('connApply13Codes', connExpr)
	
	-- Conn^i_jk(x) u^j v^k
	local connExpr = (Gamma_ull'^a_bc' * paramU'^b' * paramV'^c')()
	compileTensorField('connApply23Codes', connExpr)
	
	-- Conn^i = Conn^i_jk g^jk
	local connExpr = (Gamma_ull'^a_b^b')()
	compileTensorField('tr23_conn_u_codes', connExpr)

	-- sqrt(g)_,i / sqrt(g) - c_ij^j = Conn^j_ji
	local connExpr = (Gamma_ull'^b_ba')()
	compileTensorField('tr12_conn_l_codes', connExpr)

	-- sqrt(g)_,i / sqrt(g) = Conn^j_ij
	local connExpr = (Gamma_ull'^b_ab')()
	compileTensorField('tr13_conn_l_codes', connExpr)

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

	self.lenExprs = symmath.Array:lambda({dim}, function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'^b' * gHol'_ab')()
		local lenExpr = symmath.sqrt(lenSqExpr)()
		return lenExpr
	end)
	printNonZeroField'lenExprs'

	-- dx is the change across the grid
	-- therefore it is based on the holonomic metric
	compileTensorField('dxCodes', self.lenExprs)


	local integralArgs = table()
	for i=1,dim do
		integralArgs:insert(symmath.var('u'..i..'L'))
		integralArgs:insert(symmath.var('u'..i..'R'))
	end
	
	-- mapping 
	-- u#L -> pt.s# - .5 * grid_dx#
	-- u#R -> pt.s# + .5 * grid_dx#
	local mappedIntegralArgs = integralArgs:map(function(var)
		local i, LR = var.name:match'^u(%d)([LR])$'
		assert(i and LR)
		local addsub = assert(({L='-', R='+'})[LR])
		return {['({pt^'..i..'} '..addsub..' .5 * solver->grid_dx.s'..(i-1)..')'] = var}
	end)
	
	local coord_area_exprs = symmath.Array:lambda({dim}, function(i)
		local area = const(1)
		for j=1,dim do
			if j ~= i then
				area = area * self.lenExprs[j]
			end
		end
		area = area()

		local params = table()
		for j=1,dim do
			params:append{uL, uR}
			if j ~= i then
				local u = self.baseCoords[j]
				local uL, uR = integralArgs[2*j-1], integralArgs[2*j]
				area = area:integrate(u, uL, uR)()
			end
		end

		if self.verbose then
			print(var'area'('_'..i):eq(area))
		end

		-- TODO add in extra code function parameters
		return area
	end)

	-- area of the side in each direction
	compileTensorField('cell_area_codes', coord_area_exprs, mappedIntegralArgs)

	do
		local volume = const(1)
		for j=1,dim do
			volume = volume * self.lenExprs[j]
		end
		local volumeSq = (volume^2)()
		local gHolDet = Matrix.determinant(gHol)()
		if volumeSq ~= gHolDet then
			print('gHolDet')
			print(gHolDet)
			print('volumeSq')
			print(volumeSq)
			error'these should be the same'
		end
		for j=1,dim do
			local u = self.baseCoords[j]
			local uL, uR = integralArgs[2*j-1], integralArgs[2*j]
			volume = volume:integrate(u, uL, uR)()
		end
		if self.verbose then
			print(var'volume':eq(volume))
		end
		compileTensorField('cell_volume_code', volume, mappedIntegralArgs)
	end

	self.g = g
	compileTensorField('gCode', self.g)

	local gU = Tensor('^ab', table.unpack((Matrix.inverse(g))))
	self.gU = gU
	compileTensorField('gUCode', self.gU)
	
	local sqrt_gU = Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
	compileTensorField('sqrt_gUCode', sqrt_gU)

	local det_g_expr = symmath.Matrix.determinant(g)
	compileTensorField('det_g_code', det_g_expr)

	local partial_det_g_expr = Tensor('_i', function(i)
		return det_g_expr:diff(coords[i])()
	end)
	compileTensorField('partial_det_g_code', partial_det_g_expr)

	local partial2_det_g_expr = partial_det_g_expr'_i,j'()
	compileTensorField('partial2_det_g_code', partial2_det_g_expr)

	local sqrt_det_gExpr = symmath.sqrt(det_g_expr)()
	compileTensorField('sqrt_det_gCode', sqrt_det_gExpr)


	-- Now, for num rel in spherical, I need to transform the metric to remove the singularities ...
	-- TODO this is redundant, it is identical to anholonomic normalized 'eHolToE'

	local J = Tensor('^a_b', function(a,b)
		return a ~= b
			and 0
			or (self.anholonomic 
					and 1 
					or (1 / self.lenExprs[a])
				)
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
	compileTensorField('gJ_uu_code', gJ_uu)

	local sqrt_gJ_uu_expr = Tensor('^ab', function(a,b) return symmath.sqrt(gJ_uu[a][b])() end)
	compileTensorField('sqrt_gJ_uu_code', sqrt_gJ_uu_expr)

	local dgJ_lll = (self.dg'_abc' * J'^a_u' * J'^b_v' * J'^c_w')()
	if self.verbose then
		print'g_ij,k in rescaled coordinates:'
		print(dgJ_lll)
	end
end


local xs = table{'x', 'y', 'z'}

-- for code generation
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
	code = code:gsub('{w^(%d)}', function(i)
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

local function getCode_real3_real3_real3_to_real(name, expr)
	return template([[
real <?=name?>(real3 u, real3 v, real3 w, real3 pt) {
	return <?=convertParams(expr)?>;
}]], {
		name = name,
		expr = expr,
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

-- this is my exception to the rule, which accepts a pointer
local function getCode_real3_to_3sym3x3(name, exprs)
	return template([[
void <?=name?>(_3sym3 a[3], real3 pt) {
<?
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		for l,xl in ipairs(xNames) do
?>	a[<?=l-1?>].<?=xi?>.<?=xjk?> = <?

if exprs[i] and exprs[i][j] and exprs[i][j][k] and exprs[i][j][k][l] then
	?><?=convertParams(exprs[i][j][k][l])?><?
else
	?>0.<?
end
?>;
<?		end
	end
end
?>}
]], {
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
	if self.code then return self.code end
	
	self.solver = solver
	-- 3 since all our base types are in 'real3', 'sym3', etc
	-- what about removing this restriction?
	local dim = 3
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
	
	-- dx0, dx1, ...
	-- this is the change in cartesian wrt the change in grid
	-- this is also the normalization factor for the anholonomic ( ... is it?)
	lines:append(range(dim):mapi(function(i)
		local code = convertInputFromCoord(self.dxCodes[i], 'pt')
		return '#define coord_dx'..(i-1)..'(pt) ('..code..')'
	end))

	-- area0, area1, ...
	-- area_i = integral of u_j, j!=i of product of dx_j, j!=i
	lines:append(range(dim):mapi(function(i)
		local code = convertInputFromCoord(self.cell_area_codes[i], 'pt')
		return '#define cell_area'..(i-1)..'(pt) ('..code..')'
	end))

	lines:insert('#define cell_volume(pt) ('..convertInputFromCoord(self.cell_volume_code, 'pt')..')')

	lines:insert'\n'

	lines:append(range(dim):mapi(function(i)
		return '#define cell_dx'..(i-1)..'(pt) (coord_dx'..(i-1)..'(pt) * solver->grid_dx.s'..(i-1)..')'
	end))
	
	lines:insert'\n'

	-- metric determinant ...  det_g = volume^2 for holonomic basis
	local det_g_code = '(' .. self.det_g_code .. ')'
	lines:insert(getCode_real3_to_real('coord_det_g', det_g_code))

	lines:insert(getCode_real3_to_real3('coord_partial_det_g', self.partial_det_g_code))
	lines:insert(getCode_real3_to_sym3('coord_partial2_det_g', self.partial2_det_g_code))

	-- sqrt_det_g ... volume for holonomic basis
	local sqrt_det_gCode = '(' .. self.sqrt_det_gCode .. ')'
	lines:insert(getCode_real3_to_real('coord_sqrt_det_g', sqrt_det_gCode))
	
	-- coord len code: l(v) = v^i v^j g_ij
	lines:insert(getCode_real3_real3_to_real('coordLenSq', self.uLenSqCode))
	lines:insert[[
real coordLen(real3 r, real3 pt) {
	return sqrt(coordLenSq(r, pt));
}]]
	lines:insert(getCode_real3_to_real3('coord_tr23_c', self.tr23_cCode))
	lines:insert(getCode_real3_real3_real3_to_real('coord_conn_apply123', self.connApply123Code))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply12', self.connApply12Codes))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply13', self.connApply13Codes))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply23', self.connApply23Codes))
	lines:insert(getCode_real3_to_real3('coord_conn_trace12', self.tr12_conn_l_codes))
	lines:insert(getCode_real3_to_real3('coord_conn_trace13', self.tr13_conn_l_codes))
	lines:insert(getCode_real3_to_real3('coord_conn_trace23', self.tr23_conn_u_codes))
	lines:insert(getCode_real3_to_3sym3('coord_dg_lll', self.dg_lll_codes))
	lines:insert(getCode_real3_to_sym3sym3('coord_d2g_llll', self.d2g_llll_codes))
	lines:insert(getCode_real3_to_3sym3('coord_conn_lll', self.conn_lll_codes))
	lines:insert(getCode_real3_to_3sym3('coord_conn_ull', self.conn_ull_codes))
	-- TODO why compute all of these for each eqn that doesn't use them?
	-- how about, instead, compute these upon request ...
	-- and while you're at it you can inline them everywhere as well (or is that too much?)
	lines:insert(getCode_real3_to_3sym3x3('coord_partial_conn_ulll', self.partial_conn_ulll_codes))


	--[[
	for i=0,dim-1 do
		lines:insert(getCode_real3_to_real('coordHolBasisLen'..i, self.eHolLenCode[i+1]))
	end
	--]]

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

]], 	{
			solver = solver,
			xNames = xNames,
			symNames = symNames,
			from6to3x3 = from6to3x3 ,
		}))
	end

	lines:insert(self:getCoordMapCode())

--print(require 'template.showcode'(lines:concat'\n'))

	self.code = lines:concat'\n'
	return self.code
end

function CoordinateSystem:getCoordMapCode()
	local lines = table()
		
	lines:insert(getCode_real3_to_real3(
		'coordMap',
		range(3):mapi(function(i)
			return self.uCode[i] or '{pt^'..i..'}'
		end)))
	
	for i,eiCode in ipairs(self.eCode) do
		lines:insert(getCode_real3_to_real3('coordBasis'..(i-1), eiCode))
	end

	lines:insert(template([[
//converts a vector from cartesian coordinates to grid coordinates
//by projecting the vector into the grid basis vectors 
//at x, which is in grid coordinates
real3 cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3_zero;
	<? for i=0,solver.dim-1 do 
	local xi = xNames[i+1]
	?>{
		real3 e = coordBasis<?=i?>(pt);
<? if coord.anholonomic then	-- anholonomic normalized
?>		real uei = real3_dot(e, u);
<? else		-- holonomic
?>		real uei = real3_dot(e, u) / real3_lenSq(e);
<? end		
?>		uCoord.<?=xi?> = uei;
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
	<? for i=0,solver.dim-1 do 
	local xi = xNames[i+1]
	?>{
		real3 e = coordBasis<?=i?>(pt);
		uGrid = real3_add(uGrid, real3_real_mul(e, u.<?=xi?>));
	}<? end ?>
	return uGrid;
}
]], {
		solver = self.solver,
		coord = self,
		xNames = xNames,
	}))
	
	return lines:concat'\n'
end

function CoordinateSystem:getCoordMapGLSLCode()
	return table{
		'#define inline',
		'#define real float',
		'#define _real3 vec3',
		'#define real3 vec3',
		'#define real3_zero vec3(0.,0.,0.)',
		'#define real3_add(a,b) ((a)+(b))',
		'#define real3_sub(a,b) ((a)-(b))',
		'#define real3_real_mul(a,b) ((a)*(b))',
		'#define real3_dot dot',
		'real real3_lenSq(real3 a) { return dot(a,a); }',
		self:getCoordMapCode(),
	}:concat'\n'
end

-- until I get inverses on trig functions working better, I'll have this manually specified
function CoordinateSystem:getCoordMapInvGLSLCode()
	return [[
vec3 coordMapInv(vec3 x) { return x; }
]]
end

return CoordinateSystem
