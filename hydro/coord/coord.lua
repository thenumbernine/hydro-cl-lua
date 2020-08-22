--[[

for bssnok-fd
5495646 bytes with static inline 
5529070 bytes without it
dif of 30k

TODO get rid of all the bloat that I added for bssn
because I'm moving it directly into BSSN



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


Also I need a better naming system for all the expressions, and what they pertain to
we have 
- components based on the coordinate metric (metric is non-identity)
- components based on the non-coordinate orthonormal metric (commutation is non-zero)
- translations between coord, orthonormal, and a cartesian global basis
Also, instead of doing these calcuations over and over again, especially for things like parallel propagators, or more-complicated expressions, just run a script once to perform all calculations, update the chart properties, and store them in one place.
--]]


--[[
NEW IDEA

store the coordinte chart as a buffer.
dynamically generate the structure.
create as many fields as necessary.
let the lua object designate which 3 coordinates are x1 x2 x3 (chart parameters) and which are u1 u2 u3 (chart mapping in embedded space)
...and any other aux vars as well (i.e. radial coord, x,y,z, etc)
	(this would be more interchangeable with the unstructured-mesh solver)
]]

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local clnumber = require 'cl.obj.number'
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local CoordinateSystem = class()

local print = print

--[[
args:
	anholonomic = use anholonomic coordinate system (normalized basis vectors)
		this causes the Coord object to use commutation coefficients, and therefore will have asymmetric connections

	TODO replace that with
	vectorComponent
		holonomic = default, e_i = ∂_i = ∂/∂x^i, no guarantees of orthogonality or normalization
		anholonomic = same as above but orthonormalized
		cartesian = cartesian components even in the presence of curvilinear coordinates 
--]]
function CoordinateSystem:init(args)
	self.solver = assert(args.solver)
	self.replvars = self.replvars or table()

	local symmath = require 'symmath'
	local const = symmath.Constant
		
	self.verbose = cmdline.coordVerbose

	self.vectorComponent = args.vectorComponent or 'holonomic'
assert(args.anholonomic == nil, "coord.anholonomic is deprecated.  instead you should use coord.vectorComponent=='anholonomic'")

	if self.verbose then
		symmath.tostring = symmath.export.MathJax
		symmath.tostring.useCommaDerivative = true
		print(symmath.tostring.header)
		print = symmath.tostring.print
	end

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
	if self.vectorComponent == 'holonomic' then
		self.coords = table(baseCoords)
	elseif self.vectorComponent == 'anholonomic' then
		-- TODO this is why symmath needs CAS function objects
		-- and instead of overriding :applyDiff, just make operators a CAS function object
		local nonCoords = table()
		local nonCoordLinExpr = (eHolToE * Matrix(baseCoords):T())()
		for i=1,3 do
			local baseCoord = baseCoords[i]
			-- the non-coordinate = the coordinate, so use the original variable 
			if nonCoordLinExpr[i] == baseCoord then
			-- the non-coordinate ~= the coordinate, so make a new non-coord var and modify its 'applyDiff' function
				nonCoords[i] = baseCoord
			else
				local nonCoord = symmath.var('\\hat{'..baseCoord.name..'}')
				nonCoord.base = baseCoord
				function nonCoord:applyDiff(x)
					local xPartial = symmath.Matrix:lambda({dim, 1}, function(j,_)
						return x:diff(baseCoords[j])
					end)
					local result = Matrix:lambda({1,dim}, function(_,j) 
						if symmath.Array.is(eHolToE[i][j]) then
							io.stderr:write('eHolToE:\n'..eHolToE..'\n')
							io.stderr:write('eHolToE['..i..']['..j..']:\n'..eHolToE[i][j]..'\n')
							error'invalid eHolToE'
						end
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
	elseif self.vectorComponent == 'cartesian' then
		--start with baseCoords
		-- look at chart derivs
		-- set 'e' to the inverse
		-- TODO or don't
		-- if you are going to use a cartesian normal of the cell as your boundary then the embedded should be the basis
		-- in fact, what is 'e' used for?
		-- e is copied to eHol (only for 'holonomic'), which is used for gHol, which is used for the volume element of the grid (which cartesian will use)
-- TODO but technically, if our manifold coordinate system is non-cartesian, 
--  then our anholonomic transform should be the inverse of the chart
		self.coords = self.embedded
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

	-- comma derivatives are non-coordinates
	local e
	if self.vectorComponent == 'cartesian' then
		-- one option: use the inverse of the coordinate derivatives of the chart
		-- another option: use cartesian
		e = Tensor('_u^I', function(a,I) return a==I and 1 or 0 end)
	else
		e = Tensor'_u^I'
		e['_u^I'] = u'^I_,u'()
	end
	if self.verbose then
		print'embedded:'
		print(var'e''_u^I':eq(var'u''^I_,u'):eq(e'_u^I'()))
		print()
	end


	-- for the sake of grid lengths, 
	-- I will need the basis and metric of the holonomic version as well
	local eHol
	if self.vectorComponent == 'holonomic' then
		eHol = e
	else
		eHol = Tensor('_u^I', function(a,I)
			return u[I]:diff(baseCoords[a])()
		end)
		if self.verbose then
			print'holonomic embedded:'
			print(var'e''_u^I':eq(var'u''^I_,u'):eq(eHol'_u^I'()))
			print(var'eHol''_i^j':eq(eHolToE))
		end
	end

	-- commutation coefficients
	local c = Tensor'_ab^c'
	if self.vectorComponent == 'anholonomic' then
		if self.verbose then
			print'connection coefficients:'
			print(var'c''_uv^w' * var'e''_w','$=[ e_u, e_v ]$')
		end
		for i,ui in ipairs(coords) do
			for j,uj in ipairs(coords) do
				local zeta = var('\\zeta', baseCoords)
				local diff = ui:applyDiff(uj:applyDiff(zeta)) - uj:applyDiff(ui:applyDiff(zeta))
				local diffEval = diff()
				if diffEval ~= const(0) then
					if self.verbose then
						print('$[',ui.name,',',uj.name,'] =$',diff:eq(diffEval))
					end
					diff = diff()
					if self.verbose then
						print('factor division',diff)
					end
					local dpsi = table.mapi(baseCoords, function(uk) return zeta:diff(uk) end)
					if self.verbose then
						print('dpsi', dpsi:unpack())
					end
					local A,b = symmath.factorLinearSystem({diff}, dpsi)
					-- now extract zeta:diff(uk)
					-- and divide by e_k to get the correct coefficient
					-- TODO this assumes that e_a is only a function of partial_a
					-- if e_a is a linear combination of e_a^b partial_b then you can work it out to find
					-- c_ab^d = (e^-1)_c^d (e_a^r e_b^c_,r - e_b^r e_a^c_,r)
					-- TODO put this somewhere else so everyone can use it
					assert(b[1][1] == const(0))
					for k,uk in ipairs(coords) do
						local coeff = (A[1][k] * dpsi[k] / uk:applyDiff(zeta))()
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

	local paramU = Tensor('^i', function(i) return var('u.'..xNames[i]) end)
	local paramV = Tensor('^i', function(i) return var('v.'..xNames[i]) end)
	local paramW = Tensor('^i', function(i) return var('w.'..xNames[i]) end)

	local function printNonZero(name, code)
		if not self.verbose then return end
		local codetype = type(code)
		if codetype == 'string' then
			if code ~= '0.' then
				print(name..' = '..code)
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
	local function compileTensor(expr)
		if symmath.Array.is(expr) then
			return table.mapi(expr, function(expri) 
				return compileTensor(expri)
			end)
		elseif symmath.Expression.is(expr) then
			return self:compile(expr)
		elseif type(expr) == 'number' then
			return clnumber(expr)
		else
			error("don't know how to compile "
				..symmath.Verbose(expr))
				--..require 'ext.tolua'(expr))
		end
	end
	local function compileTensorField(field, expr)
		self[field] = compileTensor(expr)
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


	if self.vectorComponent == 'cartesian' then
		local eHolLen = range(#eHol):mapi(function(i)
			return symmath.sqrt(
				range(#eHol):mapi(function(j)
					return eHol[i][j]^2
				end):sum()
			)()
		end)

		local eHolUnitExt = symmath.Array:lambda({dim, dim}, function(i,j)
			return (eHol[i][j] / eHolLen[i])()
		end)
		compileTensorField('eHolUnitCode', eHolUnitExt)
	end

--[=[ not being used
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

	local Gamma_lll = Tensor'_abc'
	Gamma_lll['_abc'] = ((dg'_cab' + dg'_bac' - dg'_abc' + c'_abc' + c'_acb' - c'_bca') / 2)()
	if self.verbose then
		print'1st kind Christoffel:'
		print(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(Gamma_lll'_abc'()))
	end
self.Gamma_lll = Gamma_lll
	compileTensorField('conn_lll_codes', Gamma_lll)

	local Gamma_ull = Tensor'^a_bc'
	Gamma_ull['^a_bc'] = Gamma_lll'^a_bc'()
self.Gamma_ull = Gamma_ull
	if self.verbose then
		print'connection:'
		print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma_ull'^a_bc'()))
	end
	compileTensorField('conn_ull_codes', Gamma_ull)

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
	if self.vectorComponent == 'holonomic' then
		gHol = g
	else
		gHol = (eHol'_u^I' * eHol'_v^J' * eta'_IJ')()
		if self.verbose then
			print'holonomic metric:'
			print(var'gHol''_uv':eq(var'eHol''_u^I' * var'eHol''_v^J' * var'\\eta''_IJ'):eq(gHol'_uv'()))
		end
	end

	-- not really a tensor.
	self.lenExprs = Tensor('_i', function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'^b' * gHol'_ab')()
		local lenExpr = symmath.sqrt(lenSqExpr)()
		return lenExpr
	end)
	printNonZeroField'lenExprs'

	-- dx is the change across the grid
	-- therefore it is based on the holonomic metric
	compileTensorField('dxCodes', self.lenExprs)


	local integralGridDx = range(dim):mapi(function(i)
		return symmath.var('solver->grid_dx.'..xNames[i])
	end)
	local integralArgs = table()
	for i=1,dim do
		local u = self.baseCoords[i]
		integralArgs:insert(u - .5 * integralGridDx[i])
		integralArgs:insert(u + .5 * integralGridDx[i])
	end

-- [[ no one is using this yet
	local coord_area_exprs = symmath.Array:lambda({dim}, function(i)
		local area = const(1)
		for j=1,dim do
			if j ~= i then
				area = area * self.lenExprs[j]
			end
		end
		area = area()

		for j=1,dim do
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
	-- TODO only by request -- for finite-volume solvers 
	if require 'hydro.solver.fvsolver'.is(assert(args.solver)) then
		compileTensorField('cell_area_codes', coord_area_exprs)
	end

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
--print('volume was', volume)
--print('integrating', u, 'from', uL, 'to', uR)
			volume = volume:integrate(u, uL, uR)()
--print('volume is now', volume)
		end
		if self.verbose then
			print(var'volume':eq(volume))
			print(var'gHolDet':eq(gHolDet))
		end
		if require 'hydro.solver.fvsolver'.is(assert(args.solver)) then
			compileTensorField('cell_volume_code', volume)
		end
	end
--]]

	self.g = g
	compileTensorField('gCode', self.g)

	local gU = Tensor('^ab', table.unpack((Matrix.inverse(g))))
	self.gU = gU
	compileTensorField('gUCode', self.gU)
	
	local sqrt_gU = Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
	compileTensorField('sqrt_gUCode', sqrt_gU)
	
	local sqrt_g = Tensor('_ab', function(a,b) return symmath.sqrt(g[a][b])() end)
	compileTensorField('sqrt_gCode', sqrt_g)

	local det_g_expr = symmath.Matrix.determinant(g)
self.det_g = det_g_expr
	compileTensorField('det_g_code', det_g_expr)

	local sqrt_det_gExpr = symmath.sqrt(det_g_expr)()
	compileTensorField('sqrt_det_gCode', sqrt_det_gExpr)


	self:createCellStruct()
	self.cellStruct:makeType()
	self.cell_t = self.cellStruct.typename
end

function CoordinateSystem:createCellStruct()
--[[
ok here's a dilemma ...
gridSolver has cellBuf that holds cell pos and any other aux vars used for cell calculations
meshsolver has cellBuf that holds cell pos and mesh info
meshsolver needs to pass 'cellBuf'

--]]
	self.cellStruct = Struct{
		solver = self.solver,
		name = 'cell_t',
		dontUnion = true,
		vars = {
			{name='pos', type='real3'},	-- x1 x2 x3 input coordinates to the chart
		},
	}
	
	-- here's the mesh-specific stuff
	if require 'hydro.solver.meshsolver'.is(assert(self.solver)) then
		self.cellStruct.vars:append{
			{type='real', name='volume'},	--volume of the cell
			{type='int', name='faceOffset'},
			{type='int', name='faceCount'},
			{type='int', name='vtxOffset'},
			{type='int', name='vtxCount'},
		}
	end


	-- cartesian code takes ~20 seconds to compile
	-- cylindrical code takes ~60 seconds.
	-- the only main difference in the code is the # of normal computations ... and at that, directly calling cos() and sin()
	-- I'm going to see if reducing the trig calls helps but giving gridsolvers their own faceBuf
	self.faceStruct = Struct{
		solver = self.solver,
		name = 'face_t',
		dontUnion = true,
		vars = {
			{type='real3', name='pos'},		--center.  realN.
			{type='real3', name='normal'},	--normal pointing from first to second
			{type='real3', name='normal2'},	--orthonormal basis around normal
			{type='real3', name='normal3'},
			{type='real', name='area'},		--edge length / surface area
			{type='real', name='cellDist'},	--dist between cell centers along 'normal'
		},
	}

	-- here is the mesh-specific face_t fields
	if require 'hydro.solver.meshsolver'.is(assert(self.solver)) then
		self.faceStruct.vars:append{
			{type='vec2i_t', name='cells'},	--indexes of cells
			{type='int', name='vtxOffset'},
			{type='int', name='vtxCount'},
		}
	end
end

function CoordinateSystem:fillGridCellBuf(cellsCPU)
	local solver = self.solver
	local index = 0
	for k=0,tonumber(solver.gridSize.z)-1 do
		local w = solver.dim >= 3 
			and ((k + .5 - solver.numGhost) / (tonumber(solver.gridSize.z) - 2 * solver.numGhost) * (solver.maxs.z - solver.mins.z) + solver.mins.z)
			or (.5 * (solver.maxs.z + solver.mins.z))
		for j=0,tonumber(solver.gridSize.y)-1 do
			local v = solver.dim >= 2
				and ((j + .5 - solver.numGhost) / (tonumber(solver.gridSize.y) - 2 * solver.numGhost) * (solver.maxs.y - solver.mins.y) + solver.mins.y)
				or (.5 * (solver.maxs.y + solver.mins.y))
			for i=0,tonumber(solver.gridSize.x)-1 do
				local u = solver.dim >= 1
					and ((i + .5 - solver.numGhost) / (tonumber(solver.gridSize.x) - 2 * solver.numGhost) * (solver.maxs.x - solver.mins.x) + solver.mins.x)
					or (.5 * (solver.maxs.x + solver.mins.x))
				cellsCPU[index].pos.x = u
				cellsCPU[index].pos.y = v
				cellsCPU[index].pos.z = w
				index = index + 1
			end
		end
	end
end

function CoordinateSystem:compile(expr)
	if type(expr) == 'number' then return clnumber(expr) end
	
	local symmath = require 'symmath'
	local const = symmath.Constant

	for _,repl in ipairs(self.replvars) do
		expr = expr:replace(repl[1], repl[2])
	end

	-- replace pow(x, .5) with sqrt(x)
	expr = expr:map(function(x)
		if symmath.op.pow.is(x)
		and const.is(x[2])
		and x[2].value == .5
		then
			return symmath.sqrt(x[1])
		end
	end)

	for i,coord in ipairs(self.baseCoords) do
		expr = expr:replace(coord, symmath.var('pt.'..xNames[i]))
	end
print('compiling\n', expr..'\n')
	local code = symmath.export.C(expr)

	return code
end


local xs = table{'x', 'y', 'z'}

local function getCode_real3_to_real(name, code)
	return template([[
static inline real <?=name?>(real3 pt) {
	return <?=code?>;
}]], {
		name = name,
		code = code,
	})
end

-- f(x) where x is a point in the coordinate chart
local function getCode_real3_to_real3(name, exprs)
	return template([[
static inline real3 <?=name?>(real3 pt) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
	})
end

-- f(v,x) where x is a point on the coordinate chart and v is most likely a tensor
local function getCode_real3_real3_to_real(name, expr)
	return template([[
static inline real <?=name?>(real3 u, real3 pt) {
	return <?=expr?>;
}]], {
		name = name,
		expr = expr,
	})
end

local function getCode_real3_real3_to_real3(name, exprs)
	return template([[
static inline real3 <?=name?>(real3 u, real3 pt) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
	})
end

local function getCode_real3_real3_real3_to_real(name, expr)
	return template([[
static inline real <?=name?>(real3 u, real3 v, real3 w, real3 pt) {
	return <?=expr?>;
}]], {
		name = name,
		expr = expr,
	})
end

local function getCode_real3_real3_real3_to_real3(name, exprs)
	return template([[
static inline real3 <?=name?>(real3 u, real3 v, real3 pt) {
	return _real3(
<? for i=1,3 do
?>		<?=exprs[i] or '0.'
		?><?=i==3 and '' or ','?>
<? end
?>	);
}]], {
		name = name,
		exprs = exprs,
	})
end


local function getCode_real3_to_sym3(name, exprs)
	return template([[
sym3 <?=name?>(real3 pt) {
	return (sym3){
<? for i=1,3 do
	for j=i,3 do
?>		.<?=xs[i]..xs[j]?> = <?=exprs[i] and exprs[i][j] 
			and exprs[i][j] or '0.'?>,
<?	end
end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
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
			and exprs[i][j][k] or '0.'?>,
<?	end	
?>	},
<?
end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
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
	?><?=exprs[i][j][k][l]?><?
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
?>			.<?=xij?> = <?=exprs[i] and exprs[i][j] and exprs[i][j][k] and exprs[i][j][k][l] or '0.'?>,
<?	end
?>		},
<? end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
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
function CoordinateSystem:initCodeModules()
	-- 3 since all our base types are in 'real3', 'sym3', etc
	-- what about removing this restriction?
	local dim = 3
	local lines = table()

	-- dx0, dx1, ...
	-- this is the change in cartesian wrt the change in grid
	-- this is also the normalization factor for the anholonomic ( ... is it?)
	lines:append(range(dim):mapi(function(i)
		local code = self.dxCodes[i]
		return '#define coord_dx'..(i-1)..'(pt) ('..code..')'
	end))

-- [[
	-- area0, area1, ...
	-- area_i = integral of u_j, j!=i of product of dx_j, j!=i
	if require 'hydro.solver.fvsolver'.is(self.solver) then
		lines:append(range(dim):mapi(function(i)
			local code = self.cell_area_codes[i]
			return '#define cell_area'..(i-1)..'(pt) ('..code..')'
		end))

		lines:insert('#define cell_volume(pt) ('..self.cell_volume_code..')')
	end
--]]

	lines:insert'\n'

	lines:append(range(dim):mapi(function(i)
		return '#define cell_dx'..(i-1)..'(pt) (coord_dx'..(i-1)..'(pt) * solver->grid_dx.s'..(i-1)..')'
	end))
	
	lines:insert'\n'

	-- metric determinant ... det_g = volume^2 for holonomic basis
	local det_g_code = '(' .. self.det_g_code .. ')'
	lines:insert(getCode_real3_to_real('coord_det_g', det_g_code))

	-- sqrt_det_g ... volume for holonomic basis
	local sqrt_det_gCode = '(' .. self.sqrt_det_gCode .. ')'
	lines:insert(getCode_real3_to_real('coord_sqrt_det_g', sqrt_det_gCode))

	-- coord len code: l(v) = v^i v^j g_ij
	lines:insert(getCode_real3_real3_to_real('coordLenSq', self.uLenSqCode))
	lines:insert[[
static inline real coordLen(real3 r, real3 pt) {
	return sqrt(coordLenSq(r, pt));
}]]

	do
		local lines = table()
		lines:insert(getCode_real3_to_real3('coord_tr23_c', self.tr23_cCode))
		lines:insert(getCode_real3_real3_real3_to_real('coord_conn_apply123', self.connApply123Code))
		lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply12', self.connApply12Codes))
		lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply13', self.connApply13Codes))
		lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply23', self.connApply23Codes))
		lines:insert(getCode_real3_to_real3('coord_conn_trace12', self.tr12_conn_l_codes))
		lines:insert(getCode_real3_to_real3('coord_conn_trace13', self.tr13_conn_l_codes))
		lines:insert(getCode_real3_to_real3('coord_conn_trace23', self.tr23_conn_u_codes))
		lines:insert(getCode_real3_to_3sym3('coord_conn_lll', self.conn_lll_codes))
		lines:insert(getCode_real3_to_3sym3('coord_conn_ull', self.conn_ull_codes))
		self.solver.modules:add{
			name = 'conn',
			depends = {'real3', 'sym3', '_3sym3'},
			code = lines:concat'\n',
		}
	end

	--[[
	for i=0,dim-1 do
		lines:insert(getCode_real3_to_real('coordHolBasisLen'..i, self.eHolLenCode[i+1]))
	end
	--]]

	self.solver.modules:add{
		name = 'coord_lower',
		code = getCode_real3_real3_to_real3('coord_lower', self.lowerCodes),
	}

	self.solver.modules:add{
		name = 'coord_raise',
		code = getCode_real3_real3_to_real3('coord_raise', self.raiseCodes),
	}

	do
		local function addSym3Components(name, codes)
			local lines = table()
			for i=1,3 do
				for j=i,3 do
					local code = (codes[i] and codes[i][j] or clnumber(i==j and 1 or 0))
					lines:insert('#define '..name..(i-1)..(j-1)..'(pt) '..code)
					if i ~= j then
						lines:insert('#define '..name..(j-1)..(i-1)..'(pt) '..code)
					end
				end
			end
			return lines:concat'\n'
		end

		self.solver.modules:add{
			name = 'coord_g_ll##',
			code = addSym3Components('coord_g_ll', self.gCode),
		}

		self.solver.modules:add{
			name = 'coord_g_uu##',
			code = addSym3Components('coord_g_uu', self.gUCode),
		}
		
		-- curvilinear grid normals use sqrt(g^ii), as it is the metric-weighted coordinate normal magnitude
		self.solver.modules:add{
			name = 'coord_sqrt_g_uu##',
			code = addSym3Components('coord_sqrt_g_uu', self.sqrt_gUCode),
		}
	
		self.solver.modules:add{
			name = 'coord_sqrt_g_ll##',
			code = addSym3Components('coord_sqrt_g_ll', self.sqrt_gCode),
		}
	end

	self.solver.modules:add{
		name = 'coord_g_ll',
		depends = {'sym3'},
		code = getCode_real3_to_sym3('coord_g_ll', self.gCode),
	}

	self.solver.modules:add{
		name = 'coord_g_uu',
		depends = {'sym3'},
		code = getCode_real3_to_sym3('coord_g_uu', self.gUCode),
	}

	self:initCodeModule_coordMap()

	-- parallel propagate code
	if require 'hydro.solver.fvsolver'.is(self.solver) then
		local lines = table()
		
		-- parallel-propagate a vector from point 'x' along coordinate 'k' (suffix of func name) by amount 'dx'
		-- TODO derive this from the metric
		-- upper parallel propagator = exp(-int_x^(x+dx) of conn_side dx)
		-- lower parallel propagator = exp(int_x^(x+dx) of conn_side^T dx) = upper^-T
		--
		-- if we're using a cartesian basis then no need to transport anything
		-- ... (? except maybe the flux differential across the cell, done in curvilinear coordinates ?)
		if self.vectorComponent == 'cartesian' 
		or require 'hydro.coord.cartesian'.is(self)
		then
			-- general case for a fixed global orthonormal basi:
			lines:insert(template([[

// Propagate an upper index / propagate the components of a vector.
<? for side=0,solver.dim-1 do
?>#define coord_parallelPropagateU<?=side?>(v, x, dx) (v)
<? end ?>

// To propagate a one-form, apply the inverse-transpose.
// In the case of anholonomic orthonormal basis, the transpose equals the inverse so propagateL == propagateR.
//  This fits with the fact that the metric of an orthonormal basis is identity, so the upper and lower components are equal.
<? for side=0,solver.dim-1 do 
?>#define coord_parallelPropagateL<?=side?>(v, x, dx) (v)
<? end ?>

]], 		{
				solver = self.solver,
			}))
		
		-- anholonomic orthonormal = metric is identity, so inverse = transpose, so upper propagate = lower propagate ... so we only need to define the upper 
		elseif self.vectorComponent == 'anholonomic' then	-- TODO rename this to 'orthonormal' (since anholonomic doesn't imply orthonormal ... and this assumption is only true for orthonormal basis)
			lines:insert(self:getParallelPropagatorCode())
			lines:insert(template([[
<? for side=0,solver.dim-1 do 
?>#define coord_parallelPropagateL<?=side?>(v, x, dx) coord_parallelPropagateU<?=side?>(v, x, dx)
<? end ?>
]], 		{
				solver = self.solver,
			}))
		elseif self.vectorComponent == 'holonomic' then	-- TODO rename this to 'coordinate'
			lines:insert(self:getParallelPropagatorCode())
		else
			error'here'
		end
	
		self.solver.modules:add{
			name = 'coord.coord_parallelPropagate',
			code = lines:concat'\n',
		}
	end

--print(require 'template.showcode'(lines:concat'\n'))

	self:getNormalCodeGenerator()
	self.solver.modules:add{
		name = 'coord.normal',
		typecode = self.normalTypeCode,
		code = self.normalInfoCode,
	}

	self.solver.modules:add{
		name = 'coord',
		code = lines:concat'\n',
	}

	-- store all input coordinates for our cells
	-- for holonomic/anholonomic this is just the linearly interpolated
	self.solver.modules:add{
		name = 'coord.cell_t',
		structs = {self.cellStruct},
		headercode = [[
#define cell_x(i)	cellBuf[INDEXV(i)].pos
]],
	}
end

function CoordinateSystem:initCodeModule_coordMap()
	local solver = self.solver

	-- get back the Cartesian coordinate for some provided chart coordinates
	-- TODO make this a macro based on cellBuf[index]
	-- and make it custom per coord system (just like the cellBuf fields are)
	solver.modules:add{
		name = 'coordMap',
		code = getCode_real3_to_real3('coordMap', range(3):mapi(function(i) return self.uCode[i] or 'pt.'..xNames[i] end)),
	}

	-- get back the radial distance for some provided chart coordinates
	solver.modules:add{
		name = 'coordMapR',
		code = getCode_real3_to_real('coordMapR', self:compile(self.vars.r)),
	}

	do
		local lines = table()
		for i,eiCode in ipairs(self.eCode) do
			lines:insert(getCode_real3_to_real3('coordBasis'..(i-1), eiCode))
		end
		solver.modules:add{
			name = 'coordBasis#',
			code = lines:concat'\n',
		}
	end

	if self.eHolUnitCode then
		local lines = table()
		for i,eHolUnitiCode in ipairs(self.eHolUnitCode) do
			lines:insert(getCode_real3_to_real3('holunit_coordBasis'..(i-1), eHolUnitiCode))
		end
		solver.modules:add{
			name = 'holunit_coordBasis#',
			code = lines:concat'\n',
		}
	end

	do
		local fromlines = table()
		local tolines = table()
		local depends = table()
		local env = {
			solver = self.solver,
			coord = self,
			xNames = xNames,
		}
		if self.vectorComponent == 'cartesian' 
		or require 'hydro.coord.cartesian'.is(coord)
		then
			if not require 'hydro.coord.cartesian'.is(coord) then
				
				depends:insert'holunit_coordBasis#'
				tolines:insert(template([[
//converts a vector from cartesian coordinates to grid curvilinear coordinates
//by projecting the vector into the grid basis vectors
//at x, which is in grid curvilinear coordinates
real3 coord_cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3_zero;
	<? for i=0,solver.dim-1 do 
	local xi = xNames[i+1]
	?>{
		real3 e = holunit_coordBasis<?=i?>(pt);
<? if coord.vectorComponent == 'anholonomic' then	-- anholonomic normalized
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
]], env))
				fromlines:insert(template([[
//converts a vector from cartesian to grid curvilinear coordinates
//by projecting it onto the basis ... ?
real3 coord_cartesianFromCoord(real3 u, real3 pt) {
	real3 uGrid = real3_zero;
	<? for i=0,solver.dim-1 do 
	local xi = xNames[i+1]
	?>{
		real3 e = holunit_coordBasis<?=i?>(pt);
		uGrid = real3_add(uGrid, real3_real_mul(e, u.<?=xi?>));
	}<? end ?>
	return uGrid;
}
]], env))

			else	-- cartesian.is(coord)
				fromlines:insert[[
#define coord_cartesianFromCoord(u, pt) (u)
]]
				tolines:insert[[
#define coord_cartesianToCoord(u, pt) 	(u)
]]
			end 

--TODO change names
--'coord' is ambiguous
-- this is relative to teh vector component basis
			fromlines:insert[[
#define cartesianFromCoord(u, pt) 	(u)
]]
			tolines:insert[[
#define cartesianToCoord(u, pt) 	(u)
]]

		else	-- coord.vectorComponent

			depends:insert'coordBasis#'
			tolines:insert(template([[
//converts a vector from cartesian coordinates to grid curvilinear coordinates
//by projecting the vector into the grid basis vectors
//at x, which is in grid curvilinear coordinates
real3 coord_cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3_zero;
	<? for i=0,solver.dim-1 do 
	local xi = xNames[i+1]
	?>{
		real3 e = coordBasis<?=i?>(pt);
<? if coord.vectorComponent == 'anholonomic' then	-- anholonomic normalized
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

#define cartesianToCoord(u, pt) 	coord_cartesianToCoord(u, pt)
]], env))
			fromlines:insert(template([[
//converts a vector from cartesian to grid curvilinear coordinates
//by projecting it onto the basis ... ?
real3 coord_cartesianFromCoord(real3 u, real3 pt) {
	real3 uGrid = real3_zero;
	<? for i=0,solver.dim-1 do 
	local xi = xNames[i+1]
	?>{
		real3 e = coordBasis<?=i?>(pt);
		uGrid = real3_add(uGrid, real3_real_mul(e, u.<?=xi?>));
	}<? end ?>
	return uGrid;
}

#define cartesianFromCoord(u, pt) 	coord_cartesianFromCoord(u, pt)
]], env))

		end -- coord.vectorComponent
		
		solver.modules:add{
			name = 'cartesianFromCoord',
			depends = depends,
			code = fromlines:concat'\n',
		}
		solver.modules:add{
			name = 'cartesianToCoord',
			depends = depends,
			code = tolines:concat'\n',
		}
	end
end

local function makeGLSL(code)
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
		'',
		'real real3_lenSq(real3 a) { return dot(a,a); }',
		'',
		(code:gsub('static inline ', '')),
	}:concat'\n'
end

function CoordinateSystem:getCoordMapGLSLCode()
	return makeGLSL(self.solver.modules:getCodeAndHeader('coordMap', 'cartesianFromCoord'))
end

-- until I get inverses on systems of equations working, I'll have this manually specified
function CoordinateSystem:getCoordMapInvGLSLCode()
	return [[
vec3 coordMapInv(vec3 x) { return x; }
]]
end

function CoordinateSystem:getGeodesicGLSLCode()
	local lines = table()
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply23', self.connApply23Codes))
	return makeGLSL(lines:concat'\n')
end




--[[
TODO put this somewhere above where display functions can use it
TODO have this return a structure based on the solver and coordinate system. 
This can store all the values that would otherwise be recalculated again and again - like normal lower, upper, and magnitudes. 
In cartesian coords, and anholonomic vector components, (anything with a unit metric and coordinate-aligned normals) this can just return nothing
In non-cartesian coord, holonomic vector components this needs to hold the upper normals (the lower normal is identity) and the normal length (equal to sqrt(g^jj) for side j)
In non-cartesian coord, cartesian components, this can hold the lower normals (they are non coordinate aligned) but the upper normals match since the metric is identity ... and the normal length is 1

Coordinate System			Vector Components:		nL								nU							nLen
Cartesian					Cartesian				identity						identity					1
Cartesian					Holonomic				identity						identity					1			<- The chart is U(x,y,z) = (x,y,z), so the holonomic and anholonomic version is the same as the Cartesian-component version.
Cartesian					Orthonormal				identity						identity					1			_/
curvilinear					Cartesian				normalized chart gradient		normalized chart gradient	1
curvilinear					Holonomic				identity (coordinate gradient)	inverse metric				sqrt(g^jj)
curvilinear					Orthonormal				identity						identity					1

Random note, the curvilinear holonomic covariant normals are defined as the gradient along the coordinates: n_i = partial_i ..  for Cartesian partial_i = e_i this works out to the identity matrix.
But the curvilinear anholonomic orthonormal metric normals have to be defined as n_i = e_i = e_i^coord{i} e_coord{i}, which is a linear transformation of the gradient along the coordinates.

so for cartesian coord, or curvilinear orthonormal, we need to store nothing.
for curvilinear cartesian we need to store the chart gradient (lower normals).
for curvilinear holonomic we need to store the inverse metric and the normal length.

How to organize this?
1) put it in coord
2) have a C struct for holding all our required info.
3) have code gen for (1) calculating the struct and (2) accessing the struct.
--]]
function CoordinateSystem:getNormalCodeGenerator()
	if require 'hydro.solver.meshsolver'.is(self.solver) then
--[[
mesh vertexes are provided in Cartesian coordinates
so their normals are as well
 so face->normal will be in Cartesian components

however their use with the cell vector components requires them to be converted to whatever coord.vectorComponents specifies
maybe I will just do that up front?  save a frame basis & dual (that is orthonormal to the basis)
though for now I'll just support Cartesian / identity metric
--]]
		self.normalTypeCode = template([[
typedef struct {
	real3x3 n;
} normalInfo_t;
]])

		self.normalInfoCode = template([[
#define normalInfo_forFace(face) \
	((normalInfo_t){ \
		.n = (real3x3){ \
			.x = face->normal, \
			.y = face->normal2, \
			.z = face->normal3, \
		}, \
	})

#define normalInfo_len(n)	1.
#define normalInfo_lenSq(n)	1.

<? 
for j,xj in ipairs(xNames) do
	for i,xi in ipairs(xNames) do
?>
#define normalInfo_l<?=j?><?=xi?>(normal)		(normal.n.<?=xj?>.<?=xi?>)
#define normalInfo_u<?=j?><?=xi?>(n)			normalInfo_l<?=j?><?=xi?>(n)
#define normalInfo_l<?=j?><?=xi?>_over_len(n)	normalInfo_l<?=j?><?=xi?>(n)
#define normalInfo_u<?=j?><?=xi?>_over_len(n)	normalInfo_u<?=j?><?=xi?>(n)
<?
	end
end
?>

//v^i (nj)_i
#define normalInfo_vecDotNs(normal, v) (real3x3_real3_mul(normal.n, v))

//v^i (n1)_i
#define normalInfo_vecDotN1(normal, v) (real3_dot(normal.n.x, v))

// w_j <=> (n_j)^i w_i = v_j
#define normalInfo_vecFromNs(normal, v) \
	real3_add3( \
		real3_real_mul(normal.n.x, v.x), \
		real3_real_mul(normal.n.y, v.y), \
		real3_real_mul(normal.n.z, v.z))

]],		{
			xNames = xNames,
		})	
	else	-- not meshsolver

		if require 'hydro.coord.cartesian'.is(self)
		or self.vectorComponent == 'anholonomic'
		then
			--[[
			n_i = n^i = delta_ij for side j
			|n| = 1
			--]]
			self.normalTypeCode = template([[
typedef struct {
	int side;		//0, 1, 2
} normalInfo_t;		//nL = nU = normalBasisForSide (permutation of I), nLen = 1
]])

			-- this is overwhelmingly interface-based
			-- with two exceptions: calcDT and some displayVars
			self.normalInfoCode = template([[
<? for side=0,solver.dim-1 do ?>
#define normalInfo_forSide<?=side?>(x) \
	((normalInfo_t){ \
		.side = <?=side?>, \
	})
<? end ?>

//|n|
#define normalInfo_len(n)	1.
#define normalInfo_lenSq(n)	1.

//(nj)_i, (nj)^i, (nj)_i/|nj|, (nj)^i/|nj|
<? 
for j=1,3 do
	for i,xi in ipairs(xNames) do
?>
#define normalInfo_l<?=j?><?=xi?>(n)			(n.side == <?=(i-j)%3?> ? 1. : 0.)
#define normalInfo_u<?=j?><?=xi?>(n)			normalInfo_l<?=j?><?=xi?>(n)
#define normalInfo_l<?=j?><?=xi?>_over_len(n)	normalInfo_l<?=j?><?=xi?>(n)
#define normalInfo_u<?=j?><?=xi?>_over_len(n)	normalInfo_u<?=j?><?=xi?>(n)
<? 
	end 
end
?>

// this is the same as converting 'v' in global cartesian to 'v' in the basis of nj
// v^i (nj)_i for side j 
#define normalInfo_vecDotNs(n, v) \
	(_real3( \
		v.s[n.side], \
		v.s[(n.side+1)%3], \
		v.s[(n.side+2)%3]))

//v^i (n1)_i
#define normalInfo_vecDotN1(n, v)	(v.s[n.side])

// ...and this is the same as converting v in the basis of nj to v in global cartesian
// v.x * e[side] + v.y * e[side+1] + v.z * e[side+2]
#define normalInfo_vecFromNs(n, v) \
	(_real3( \
		v.s[(3-n.side)%3], \
		v.s[(3-n.side+1)%3], \
		v.s[(3-n.side+2)%3]))

]],			{
				eqn = self,
				solver = self.solver,
				xNames = xNames,
			})

		elseif self.vectorComponent == 'cartesian' then

			--[[
			n_i = n^i = unit(u^i_,j) for side j
			|n| = sqrt(n^i n_i) = 1 (since g_ij = g^ij = delta_ij)
			--]]
			self.normalTypeCode = template([[
typedef struct {
	real3x3 n;		// nL = nU, both are orthonormal so nLen = 1
	real len;
} normalInfo_t;
]])

			-- this would call coord_cartesianFromCoord
			-- which itself aligns with the holunit_coordBasis
			self.normalInfoCode = template([[
<? for side=0,solver.dim-1 do ?>
#define normalInfo_forSide<?=side?>(pt) \
	((normalInfo_t){ \
		.n = (real3x3){ \
			.x = holunit_coordBasis<?=side?>(pt), \
			.y = holunit_coordBasis<?=(side+1)%3?>(pt), \
			.z = holunit_coordBasis<?=(side+2)%3?>(pt), \
		}, \
		.len = 1., \
	})
<? end ?>

//|n1|
#define normalInfo_len(normal)		(normal.len)
#define normalInfo_lenSq(normal)	(normal.len * normal.len)

//(nj)_i, (nj)^i, (nj)_i / |nj|, (nj)^i / |nj|
<? 
for j,xj in ipairs(xNames) do
	for i,xi in ipairs(xNames) do
?>
#define normalInfo_l<?=j?><?=xi?>(normal)			(normal.n.<?=xj?>.<?=xi?>)
#define normalInfo_u<?=j?><?=xi?>(normal)			normalInfo_l<?=j?><?=xi?>(normal)
#define normalInfo_l<?=j?><?=xi?>_over_len(normal)	(normalInfo_l<?=j?><?=xi?>(normal) / normal.len)
#define normalInfo_u<?=j?><?=xi?>_over_len(normal)	normalInfo_l<?=j?><?=xi?>_over_len(normal)
<?
	end
end
?>

//v^i (nj)_i for side j 
#define normalInfo_vecDotNs(normal, v) (real3x3_real3_mul(normal.n, v))

//v^i (n1)_i
#define normalInfo_vecDotN1(normal, v) (real3_dot(normal.n.x, v))


// ...and this is the same as converting v in the basis of nj to v in global cartesian
// v.x * e[side] + v.y * e[side+1] + v.z * e[side+2]
#define normalInfo_vecFromNs(normal, v) \
	real3_add3( \
		real3_real_mul(normal.n.x, v.x), \
		real3_real_mul(normal.n.y, v.y), \
		real3_real_mul(normal.n.z, v.z))


]],			{
				eqn = self,
				solver = self.solver,
				xNames = xNames,
			})	
		
		elseif self.vectorComponent == 'holonomic' then

			--[[
			n_i = delta_ij for side j
			n^i = g^ik delta_kj
			|n| = sqrt(n^i n_i) = sqrt(g^jj)
			--]]
			self.normalTypeCode = template([[
typedef struct {
	int side;
	real3x3 U;
	real len;
} normalInfo_t;
		]])

			self.normalInfoCode = template([[
<? for side=0,solver.dim-1 do ?>
#define normalInfo_forSide<?=side?>(x) \
	((normalInfo_t){ \
		.side = <?=side?>, \
		.U = _real3x3( \
<? 
for j=0,2 do
	for i=0,2 do 
?>			coord_g_uu<?=(side+j)%3?><?=i?>(x)<?=i+3*j < 8 and ',' or ''?> \
<? 
	end
end 
?>		), \
		.len = coord_sqrt_g_uu<?=side..side?>(x), \
	})
<? end ?>

//|n1|
#define normalInfo_len(n)	(n.len)
#define normalInfo_lenSq(n)	(n.len * n.len)

//(nj)_i, (nj_i / |nj|
<?
for j=1,3 do
	for i,xi in ipairs(xNames) do
?>
#define normalInfo_l<?=j?><?=xi?>(n)			(n.side == <?=(i-j)%3?> ? 1. : 0.)
#define normalInfo_l<?=j?><?=xi?>_over_len(n)	(n.side == <?=(i-j)%3?> ? (1./n.len) : 0.)
<?
	end
end
?>

//(nj)^i, (nj)^i / |nj|
<?
for j,xj in ipairs(xNames) do
	for i,xi in ipairs(xNames) do
?>
#define normalInfo_u<?=j?><?=xi?>(n)			(n.U.<?=xj?>.<?=xi?>)
#define normalInfo_u<?=j?><?=xi?>_over_len(n)	(normalInfo_u<?=j?><?=xi?>(n) / n.len)
<?
	end
end
?>

//v^i (nj)_i for side j
#define normalInfo_vecDotNs(n, v) \
	(_real3( \
		v.s[n.side], \
		v.s[(n.side+1)%3], \
		v.s[(n.side+2)%3]))

//v^i (n1)_i
#define normalInfo_vecDotN1(n, v) 	(v.s[n.side])

// ...and this is the same as converting v in the basis of nj to v in global cartesian
// v.x * e[side] + v.y * e[side+1] + v.z * e[side+2]
#define normalInfo_vecFromNs(n, v) \
	(_real3( \
		v.s[(3-n.side)%3], \
		v.s[(3-n.side+1)%3], \
		v.s[(3-n.side+2)%3]))


]],			{
				eqn = self,
				solver = self.solver,
				xNames = xNames
			})

		else
			error'here'
		end


	end	-- meshsolver

	self.normalInfoCode = self.normalInfoCode .. template[[

<? for side=1,3 do ?>
#define normalInfo_l<?=side?>(n) \
	(_real3( \
		normalInfo_l<?=side?>x(n), \
		normalInfo_l<?=side?>y(n), \
		normalInfo_l<?=side?>z(n)))

#define normalInfo_u<?=side?>(n) \
	(_real3( \
		normalInfo_u<?=side?>x(n), \
		normalInfo_u<?=side?>y(n), \
		normalInfo_u<?=side?>z(n)))
<? end ?>

]]

end


return CoordinateSystem
