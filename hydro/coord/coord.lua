--[[

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
local path = require 'ext.path'
local tolua = require 'ext.tolua'
local template = require 'template'
local clnumber = require 'cl.obj.number'
local Struct = require 'hydro.code.struct'
local math = require 'ext.math'		--isfinite

local half = require 'cl.obj.half'
local fromreal, toreal = half.fromreal, half.toreal

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from6to3x3 = common.from6to3x3


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
	-- put all unique code module names here
	require 'hydro.code.symbols'(self, {
		'coord_dx_i',
		'coord_lower',
		'coord_raise',
		'coordLenSq',
		'coordLen',
		'coord_holBasisLen_i',
		'coordMap',
		'coordMapR',
		'coordMapInv',
		'coordMapGLSL',
		'coordMapInvGLSL',
		'coord_basis_i',
		'coord_basisHolUnit_i',
		'cartesianFromCoord',
		'cartesianToCoord',
		'coord_parallelPropagate',
		'normal_t',

		'coord_tr23_c',

		'coord_g_ll',
		'coord_g_ll_ij',
		'coord_g_uu',
		'coord_g_uu_ij',
		'coord_sqrt_g_uu',
		'coord_sqrt_g_uu_ij',
		'coord_sqrt_g_ll',
		'coord_sqrt_g_ll_ij',
		'coord_det_g',
		'coord_sqrt_det_g',
		'coord_partial_det_g',
		'coord_partial2_det_g',
		'coord_partial_g_lll',
		'coord_conn_lll',
		'coord_conn_ull',
		'coord_conn_apply12',
		'coord_conn_apply13',
		'coord_conn_apply23',
		'coord_conn_apply123',
		'coord_conn_trace12',
		'coord_conn_trace13',
		'coord_conn_trace23',

		'coord_gHol_ll',
		'coord_gHol_uu',
		'coord_det_gHol',
		'coord_sqrt_gHol_ll',
		'coord_sqrt_gHol_ll_ij',
		'coord_partial_det_gHol_l',
		'coord_partial2_det_gHol_ll',
		'coord_partial_gHol_lll',
		'coord_connHol_lll',
		'coord_connHol_ull',
		'coord_connHol_trace23',

		'cell_area_i',
		'cell_dx_i',
		'cell_volume',
		'cell_sqrt_det_g',	-- seems that, when this is used, it would most often be used with gHol...
		'cell_calcAvg_withPt',

		'cell_dxs',
		'cell_areas',
		'coord_basisHolUnits',
	})

	local solver = assert(args.solver)
	self.solver = solver

	-- these are for replacing one expression with another
	-- it's useful for simplifying calculations, especially complex ones involving derivatives.  just perform the derivatives separately and replace them later.
	self.repls = self.repls or table()

	-- these are for replacing values, especially dynamic values.
	-- compile your expressions with variables matching #defines in CL code.
	-- then, if you want to evaluate them, you can use this, but it isn't done as often as using 'repls'.
	self.replDefines = self.replDefines or table()

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

	-- 3 since all our base types are in 'real3', 'real3s3', etc
	-- what about removing this restriction?
	local dim = 3

	local var = symmath.var
	local Matrix = symmath.Matrix
	local Tensor = symmath.Tensor
	local frac = symmath.frac

	local eHolToE = self.eHolToE
	if not eHolToE then
		eHolToE = Matrix.identity(3)
	end

	local tangentSpaceOperators
	local baseCoords = self.baseCoords
	if self.vectorComponent == 'holonomic' then
		self.coords = table(baseCoords)
	elseif self.vectorComponent == 'anholonomic' then
		tangentSpaceOperators = table()
		local nonCoords = table()
		local nonCoordLinExpr = (eHolToE * Matrix(baseCoords):T())()
		for i=1,3 do
			local baseCoord = baseCoords[i]
			-- the non-coordinate = the coordinate, so use the original variable
			if nonCoordLinExpr[i] == baseCoord then
				nonCoords[i] = baseCoord
				tangentSpaceOperators[i] = function(x) return x:diff(baseCoord)() end

			-- the non-coordinate ~= the coordinate, so make a new non-coord var
			else
				nonCoords[i] = symmath.var('\\hat{'..baseCoord.name..'}')

				tangentSpaceOperators[i] = function(x)
					local xPartial = symmath.Matrix:lambda({dim, 1}, function(j,_)
						return x:diff(baseCoords[j])
					end)
					local result = Matrix:lambda({1,dim}, function(_,j)
						if symmath.Array:isa(eHolToE[i][j]) then
							io.stderr:write('eHolToE:\n'..eHolToE..'\n')
							io.stderr:write('eHolToE['..i..']['..j..']:\n'..eHolToE[i][j]..'\n')
							error'invalid eHolToE'
						end
						return eHolToE[i][j]
					end)
						* xPartial
					result = result()
					assert(symmath.Matrix:isa(result))
					result = result[1][1]
					assert(symmath.Expression:isa(result) and not symmath.Array:isa(result))
					return result
				end
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

	if self.verbose then
		print('flatMetric:')
		print(flatMetric)
		print('embedded:', table.mapi(embedded, tostring):concat', ')
	end

	self.manifold = Tensor.Manifold()
	self.symchart = self.manifold:Chart{
		coords = coords,
		tangentSpaceOperators = tangentSpaceOperators,
	}
	self.symEmbeddedChart = self.manifold:Chart{
		coords = embedded,
		symbols = 'IJKLMN',
		metric = function() return flatMetric end,
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
		print()
		print'point on chart, in embedded coordinates:'
		print(
			var'P':eq(
				var'P''^I' * var'e''_I'
			):eq(
				-- hmm, Tensor tostring puts one forms as cols regardless of indexing (should it? maybe I should do rows for lower indexes?)
				-- while Matrix doesn't show indexes
				Matrix(
					u'^I'()
				) * Matrix{
					var'e''_x',
					var'e''_y',
					var'e''_z'
				}:T()
			)
		)
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
		print(var'e''_I'..'= embedded basis')
		print(var'e'' _\\tilde{u}'..'$= \\partial_{\\tilde{u}} =$ chart holonomic coordinate basis')
		print(var'e'' _\\hat{u}'..'= chart anholonomic orthonormal basis')
		print()
		print'chart basis, in terms of embedded basis:'
		print()
		print(
			var'e'' _\\hat{u}'
			:eq(
				var'e'' _\\hat{u} ^I'
				* var'e''_I'
			):eq(
				--var'P'' ^I _,\\hat{u}'
				var[[e_{\hat{u}}( P^I )]]
				* var'e''_I'
			):eq(
				e'_u^I'()
				* Tensor('_I', function(I)
					return var'e'(' _'..embedded[I].name)
				end)
			):eq(
				Tensor('_u', function(u)
					return var'e'(' _'..coords[u].name)
				end)
			)
		)
		print()
	end


	-- for the sake of grid lengths,
	-- I will need the basis and metric of the holonomic version as well
	local eHol
	if self.vectorComponent == 'holonomic' then
		if self.verbose then
			print'using holonomic coordinates, so the chart basis operator is equal to the partial derivative operator'
		end
		eHol = e
	else
		if self.verbose then
			print'using non-holonomic, so separately evaluating our chart holonomic basis (for stuff like the volume element etc)'
			print()
		end
		eHol = Tensor('_u^I', function(a,I)
			return u[I]:diff(baseCoords[a])()
		end)
		if self.verbose then
			print()
		end
	end
	if self.verbose then
		print'chart holonomic basis, in terms of embedded basis:'
		print()
		print(
			var'e'' _\\tilde{u}'
			:eq(
				var'e'' _\\tilde{u} ^I'
				* var'e''_I'
			):eq(
				--var'P'' ^I _,\\tilde{u}'
				var[[e_{\tilde{u}}( P^I )]]
				* var'e''_I'
			):eq(
				eHol' _\\tilde{u} ^I'()
				* Tensor('_I', function(I)
					return var'e'(' _'..embedded[I].name)
				end)
			):eq(
				Tensor(' _\\tilde{u}', function(u)
					return var'e'(' _'..baseCoords[u].name)
				end)
			)
		)
		print()

		print'transform from chart holonomic basis to chart anholonomic orthonormal basis:'
		print()
		print(var'e'' _\\hat{u} ^\\tilde{v}':eq(eHolToE))
		print()

		print'such that'
		print()
		print(
			var'e'' _\\hat{u}'
			:eq(
				var'e'' _\\hat{u} ^\\tilde{v}'
				* var'e'' _\\tilde{v} ^I'
				* var'e'' _I'
			):eq(
				-- eHolToE is a matrix, so ..
				Tensor(' _\\hat{u} ^\\tilde{v}', table.unpack(eHolToE))
				* eHol' _\\tilde{v} ^I'()
				* Tensor('_I', function(I)
					return var'e'(' _'..embedded[I].name)
				end)
			)
		)
		print()
	end

	-- commutation coefficients
	local c = self.symchart.commutation
	if self.vectorComponent == 'anholonomic' then
		if self.verbose then
			print'connection coefficients:'
			print(var'c''_uv^w' * var'e''_w','$=[ e_u, e_v ]$')
		end
	end
	if self.verbose then
		print'commutation:'
		print()
		print(var'c'' _\\hat{u} _\\hat{v} ^\\hat{w}':eq(c' _\\hat{u} _\\hat{v} ^\\hat{w}'()))
		print()
	end

	local g = (e'_u^I' * e'_v^J' * eta'_IJ')()
	if self.verbose then
		print'metric:'	-- anholonomic if requested
		print()
		print(
			var'g'' _\\hat{u} _\\hat{v}':eq(
				var'e'' _\\hat{u} ^I'
				* var'e'' _\\hat{v} ^J'
				* var'\\eta''_IJ'
			):eq(
				g' _\\hat{u} _\\hat{v}'()
			)
		)
		print()
	end
	self.symchart:setMetric(g)


	-- code generation

	local paramU = Tensor('^i', function(i) return var('u.'..xNames[i]) end)
	local paramV = Tensor('^i', function(i) return var('v.'..xNames[i]) end)
	local paramW = Tensor('^i', function(i) return var('w.'..xNames[i]) end)

	local function printNonZero(name, code)
		if not self.verbose then return end
		local codetype = type(code)
		if codetype == 'string' then
			if code ~= '0.' then
				print(name..' = <pre style="display:inline">'..code..'</pre>')
			end
		elseif codetype == 'table' then
			for i=1,#code do
				printNonZero(name..'['..i..']', code[i])
			end
		else
			error("can't print code "..('%q'):format(name).." of type "..codetype)
		end
	end
self.printNonZero = printNonZero

	--compile a tensor of expressions to a nested table of codes
	local function compileTensor(expr)
		if symmath.Array:isa(expr) then
			return table.mapi(expr, function(expri)
				return compileTensor(expri)
			end)
		elseif symmath.Expression:isa(expr) then
			return self:compile(expr)
		elseif type(expr) == 'number' then
			return clnumber(expr)
		else
			error("don't know how to compile "
				..symmath.Verbose(expr))
				--..tolua(expr))
		end
	end
self.compileTensor = compileTensor
	local function compilePrintRequestTensor(field)
		local expr = self.request(field)
		local result = compileTensor(expr)
		printNonZero(field, result)
		return result
	end
self.compilePrintRequestTensor = compilePrintRequestTensor




	self.cached = {}
	self.calc = {}
	self.request = function(name)
		if not self.cached[name] then
			local build = assert(self.calc[name], "couldn't find coord.calc["..('%q'):format(name).."]").build
			if not build then
				error("requested calculation of '"..name.."' but couldn't find it")
			end
			self.cached[name] = build() or true
		end
		return self.cached[name]
	end

	--[[
	self.calc[moduleName] = {
		many = true,
		{
			field = (optional) key to self.request(), and also the name of the generated function. default is moduleName.
			build = function that generates / returns the Expression
			result = result type / used for depends.
			args = (optional) input args of generated function.
			define = (optional) use #define in code def instead of function
		},
	}
	or don't set many=true and just have a single entry
	--]]

	-- u is used to project the grid for displaying
	-- TODO rename to P, right, that's a convention of sorts, right?
	self.calc.u = {
		build = function()
			return u
		end,
	}

	-- extend 'e' to full R3
	-- TODO should I do this from the start?
	-- just provide the full R3 coordinates, and make no 'eExt' struct?
	self.calc.eExt = {
		build = function()
			return symmath.Array:lambda({dim, dim}, function(i,j)
				return e[i][j] or const(0)
			end)
		end,
	}

	self.calc.coord_holBasisLen = {
		build = function()
			return range(#eHol):mapi(function(i)
				return symmath.sqrt(
					range(#eHol):mapi(function(j)
						return eHol[i][j]^2
					end):sum()
				)()
			end)
		end,
	}

	self.calc.eHolUnitExt = {
		build = function()
			local eHolLen = self.request'coord_holBasisLen'
			return symmath.Array:lambda({dim, dim}, function(i,j)
				return (eHol[i][j] / eHolLen[i])()
			end)
		end,
	}

	self.calc.eExtLen = {
		build = function()
			return self.request'eExt':mapi(function(ei,i)
				return symmath.sqrt(ei:mapi(function(x) return x^2 end):sum())()
			end)
		end,
	}

	self.calc.eExitUnit = {
		build = function()
			local eExtLen = self.request'eExtLen'
			return self.request'eExt':mapi(function(ei,i)
				return ei:mapi(function(eij) return (eij/eExtLen[i])() end)
			end)
		end,
	}

	-- g_ij
	self.calc.coord_g_ll = {
		build = function()
			return g
		end,
		result = 'real3s3',
	}

	-- g^ij
	self.calc.coord_g_uu = {
		build = function()
			local g = self.request'coord_g_ll'
			return Tensor('^ab', table.unpack((Matrix.inverse(g))))
		end,
		result = 'real3s3',
	}

	-- sqrt(g^ij)
	self.calc.coord_sqrt_g_uu = {
		build = function()
			local gU = self.request'coord_g_uu'
			return Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
		end,
		result = 'real3s3',
	}

	-- sqrt(g_ij)
	self.calc.coord_sqrt_g_ll = {
		build = function()
			local g = self.request'coord_g_ll'
			return Tensor('_ab', function(a,b) return symmath.sqrt(g[a][b])() end)
		end,
		result = 'real3s3',
	}

	-- det(g_ij)
	self.calc.coord_det_g = {
		build = function()
			local g = self.request'coord_g_ll'
			return symmath.Matrix.determinant(g)
		end,
		result = 'real',
	}

	-- sqrt(det(g_ij))
	self.calc.coord_sqrt_det_g = {
		build = function()
			return symmath.sqrt(self.request'coord_det_g')()
		end,
		result = 'real',
	}

	-- det(g)_,i
	self.calc.coord_partial_det_g = {
		build = function()
			return Tensor('_a', function(a)
				return self.symchart.tangentSpaceOperators[a](self.request'coord_det_g')()
			end)
		end,
		result = 'real3',
	}

	-- det(g)_,ij
	self.calc.coord_partial2_det_g = {
		build = function()
			return self.request'coord_partial_det_g''_a,b'()
		end,
		result = 'real3s3',
	}

	-- v^k v_k
	self.calc.coordLenSq = {
		build = function()
			return (paramU'^a' * paramU'_a')()
		end,
		args = 'real3',
		result = 'real',
	}

	-- v^k v_k
	self.calc.coordLen = {
		build = function()
			return symmath.sqrt(self.request'coordLenSq')()
		end,
		args = 'real3',
		result = 'real',
	}

	-- v^k -> v_k
	self.calc.coord_lower = {
		build = function()
			return paramU'_a'()
		end,
		args = 'real3',
		result = 'real3',
	}

	-- v^k -> v_k
	self.calc.coord_raise = {
		build = function()
			return paramU'_a'()
		end,
		args = 'real3',
		result = 'real3',
	}

	-- c_ab^b
	self.calc.coord_tr23_c = {
		build = function()
			local tr23_c = c'_ab^b'()
			if self.verbose then
				print(var'c''_ab^b':eq(tr23_c))
			end
			return tr23_c
		end,
		result = 'real3',
	}

	self.calc.coord_partial_g_lll = {
		build = function()
			local g = self.request'coord_g_ll'
			local dg = g'_ab,c'():permute'_cab'
			if self.verbose then
				print()
				print'metric partial:'	-- anholonomic if requested
				print(var'g'' _\\hat{a} _\\hat{b} _,\\hat{c}':eq(dg' _\\hat{c} _\\hat{a} _\\hat{b}'()))
				print()
			end
			return dg
		end,
		result = 'real3x3s3',
	}

	-- Levi-Civita unique metric-cancelling torsion-free connection for a basis that is a linear transform of a coordinate basis
	self.calc.coord_conn_lll = {
		build = function()
			local dg = self.request'coord_partial_g_lll'
			local Gamma_lll = ((dg'_cab' + dg'_bac' - dg'_abc' + c'_abc' + c'_acb' - c'_bca') / 2)():permute'_abc'
			if self.verbose then
				print'1st kind Christoffel:'	-- anholonomic if requested
				print(
					var'\\Gamma'' _\\hat{a} _\\hat{b} _\\hat{c}':eq(
						frac(1,2)*(
							var'g'' _\\hat{a} _\\hat{b} _,\\hat{c}'
							+ var'g'' _\\hat{a} _\\hat{c} _,\\hat{b}'
							- var'g'' _\\hat{b} _\\hat{c} _,\\hat{a}'
							+ var'c'' _\\hat{a} _\\hat{b} _\\hat{c}'
							+ var'c'' _\\hat{a} _\\hat{c} _\\hat{b}'
							- var'c'' _\\hat{b} _\\hat{c} _\\hat{a}'
						)
					):eq(Gamma_lll' _\\hat{a} _\\hat{b} _\\hat{c}'())
				)
			end
			return Gamma_lll
		end,
		result = 'real3x3s3',
	}

	self.calc.coord_conn_ull = {
		build = function()
			--local g = self.request'coord_g_ll'
			local Gamma_lll = self.request'coord_conn_lll'
			local Gamma_ull = (g'^ad' * Gamma_lll'_dbc')():permute'^a_bc'
			if self.verbose then
				print'connection:'	-- anholonomic if requested
				print(var'\\Gamma'' ^\\hat{a} _\\hat{b} _\\hat{c}':eq(
					var'g'' ^\\hat{a} ^\\hat{d}'
					* var'\\Gamma'' _\\hat{d} _\\hat{b} _\\hat{c}'
				):eq(Gamma_ull' ^\\hat{a} _\\hat{b} _\\hat{c}'()))
			end
			return Gamma_ull
		end,
		result = 'real3x3s3',
	}

	-- u^j v^k Conn_jk^i(x)
	self.calc.coord_conn_apply12 = {
		build = function()
			local Gamma_lll = self.request'coord_conn_lll'
			return (paramU'^b' * paramV'^c' * Gamma_lll'_bc^a')()
		end,
		args = 'real3_real3',
		result = 'real3',
	}

	-- u^j v^k Conn_j^i_k(x)
	self.calc.coord_conn_apply13 = {
		build = function()
			local Gamma_lll = self.request'coord_conn_lll'
			return (paramU'^b' * paramV'^c' * Gamma_lll'_b^a_c')()
		end,
		args = 'real3_real3',
		result = 'real3',
	}

	-- Conn^i_jk(x) u^j v^k
	self.calc.coord_conn_apply23 = {
		build = function()
			local Gamma_ull = self.request'coord_conn_ull'
			return (Gamma_ull'^a_bc' * paramU'^b' * paramV'^c')()
		end,
		args = 'real3_real3',
		result = 'real3',
	}

	-- u^i v^j b^k Conn_ijk(x)
	self.calc.coord_conn_apply123 = {
		build = function()
			local Gamma_lll = self.request'coord_conn_lll'
			return (paramU'^a' * paramV'^b' * paramW'^c' * Gamma_lll'_abc')()
		end,
		args = 'real3_real3_real3',
		result = 'real',
	}

	-- sqrt(g)_,i / sqrt(g) - c_ij^j = Conn^j_ij - c_ij^j = Conn^j_ji
	self.calc.coord_conn_trace12 = {
		build = function()
			local Gamma_ull = self.request'coord_conn_ull'
			return (Gamma_ull'^b_ba')()
		end,
		result = 'real3',
	}

	-- sqrt(g)_,i / sqrt(g) = Conn^j_ij
	self.calc.coord_conn_trace13 = {
		build = function()
			local Gamma_ull = self.request'coord_conn_ull'
			return Gamma_ull'^b_ab'()
		end,
		result = 'real3',
	}

	-- Conn^i = Conn^i_jk g^jk
	self.calc.coord_conn_trace23 = {
		build = function()
			local gU = self.request'coord_g_uu'
			local Gamma_ull = self.request'coord_conn_ull'
			return (Gamma_ull'^a_bc' * gU'^bc')()
		end,
		result = 'real3',
	}

	self.calc.coord_gHol_ll = {
		build = function()
			if self.vectorComponent == 'holonomic' then
				return g
			else
				local gHol = (eHol'_u^I' * eHol'_v^J' * eta'_IJ')()
				if self.verbose then
					print()
					print'holonomic metric:'
					print()
					print(
						var'g'' _\\tilde{u} _\\tilde{v}':eq(
							var'e'' _\\tilde{u} ^I'
							* var'e'' _\\tilde{v} ^J'
							* var'\\eta''_IJ'
						):eq(
							gHol'_uv'()
						)
					)
					print()
				end
				return gHol
			end
		end,
		result = 'real3s3',
	}

	self.calc.coord_gHol_uu = {
		build = function()
			local gHol = self.request'coord_gHol_ll'
			return Tensor('^ab', table.unpack((Matrix.inverse(gHol))))
		end,
		result = 'real3s3',
	}

	self.calc.coord_det_gHol = {
		build = function()
			local gHol = self.request'coord_gHol_ll'
			return symmath.Matrix.determinant(gHol)
		end,
		result = 'real',
	}

	-- not really a tensor.
	-- dx is the change across the grid
	-- therefore it is based on the holonomic metric
	self.calc.coord_dx = {
		build = function()
			local gHol = self.request'coord_gHol_ll'
			return Tensor('_i', function(i)
				local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
				local lenSqExpr = (dir'^a' * dir'^b' * gHol'_ab')()
				local lenExpr = symmath.sqrt(lenSqExpr)()
				return lenExpr
			end)
		end,
	}

	-- this is a 2-point tensor though
	-- TODO just use self.eToEHol?
	self.calc.eToEHol = {
		build = function()
			local lenExprs = self.request"coord_dx"
			local e = Tensor("_i^I",
				{lenExprs[1], 0, 0},
				{0, lenExprs[2], 0},
				{0, 0, lenExprs[3]})
			return e
		end,
	}

	self.calc.eHolToE = {
		build = function()
			local lenExprs = self.request"coord_dx"
			local eInv = Tensor("^i_I",
				{1/lenExprs[1], 0, 0},
				{0, 1/lenExprs[2], 0},
				{0, 0, 1/lenExprs[3]})
			return eInv
		end,
	}

	local integralGridDx = range(dim):mapi(function(i)
		return symmath.var('solver.grid_dx.'..xNames[i])
	end)
	local integralArgs = table()
	for i=1,dim do
		local u = self.baseCoords[i]
		integralArgs:insert(u - frac(1,2) * integralGridDx[i])
		integralArgs:insert(u + frac(1,2) * integralGridDx[i])
	end

	-- replace the long variable names with math symbols
	if self.verbose then
		function self.fixVerbose(expr)
			return expr
				:replace(integralGridDx[1], var'\\Delta x_1')
				:replace(integralGridDx[2], var'\\Delta x_2')
				:replace(integralGridDx[3], var'\\Delta x_3')
		end
	else
		function self.fixVerbose(...) return ... end
	end

	-- area of the side in each direction
	self.calc.cell_area = {
		build = function()
			local lenExprs = self.request'coord_dx'
			return symmath.Array:lambda({dim}, function(i)
				local area = const(1)
				for j=1,dim do
					if j ~= i then
						area = area * lenExprs[j]
					end
				end
				area = area()

				for j=1,dim do
					if j ~= i then
						local u = self.baseCoords[j]
						local uL, uR = integralArgs[2*j-1], integralArgs[2*j]
						area = self:applyReplVars(area)	-- just because of sphere_sinh_radial, insert repls beforehand
						area = area:integrate(u, uL, uR)()
					end
				end

				if self.verbose then
					print(var'area'('_'..i):eq(self.fixVerbose(area)))
				end

				-- TODO add in extra code function parameters
				return area
			end)
		end,
	}

	self.calc.cell_volume = {
		build = function()
			local lenExprs = self.request'coord_dx'
			local volume = const(1)
			for j=1,dim do
				volume = volume * lenExprs[j]
			end
			local volumeSq = (volume^2)()
			local gHolDet = self.request'coord_det_gHol'
			if volumeSq ~= gHolDet then
				print('$det(g_{\\tilde{u}\\tilde{v}} =$', gHolDet)
				print()
				print('$vol^2 =$', volumeSq)
				print()
				error'these should be the same'
			end

			for j=1,dim do
				local u = self.baseCoords[j]
				local uL, uR = integralArgs[2*j-1], integralArgs[2*j]
				if self.verbose then
					print('volume was', self.fixVerbose(volume))
					print()
					print('integrating', u, 'from', self.fixVerbose(uL), 'to', self.fixVerbose(uR))
					print()
				end
				volume = self:applyReplVars(volume)	-- just because of sphere_sinh_radial, insert repls beforehand
				volume = volume:integrate(u, uL, uR)()
				if self.verbose then
					print('volume is now', self.fixVerbose(volume))
				print()
				end
			end
			if self.verbose then
				print()
				print(var'vol':eq(self.fixVerbose(volume)))
				print()
				print(var'det(g_{\\tilde{u}\\tilde{v}})':eq(gHolDet))
				print()
			end
			return volume
		end,
		result = 'real',
		define = 'with solver arg',	-- anything that references integraGridDx needs to be #define, or needs to add a solver param
	}

	self.calc.coord_sqrt_gHol_ll = {
		build = function()
			local gHol = self.request'coord_gHol_ll'
			return Tensor('_ab', function(a,b) return symmath.sqrt(gHol[a][b])() end)
		end,
		result = 'real3s3',
	}

	self.calc.coord_partial_gHol_lll = {
		build = function()
			local gHol = self.request'coord_gHol_ll'
			local dgHol = gHol'_ab,c'():permute'_cab'
			if self.verbose then
				print()
				print'holonomic metric partial:'
				print()
				print(
					var'g'' _\\tilde{a} _\\tilde{b} _,\\tilde{c}':eq(
						dgHol' _\\tilde{c} _\\tilde{a} _\\tilde{b}'()
					)
				)
				print()
			end
			return dgHol
		end,
		result = 'real3x3s3',
	}

	self.calc.coord_connHol_lll = {
		build = function()
			local dgHol = self.request'coord_partial_gHol_lll'
			local GammaHol_lll = (frac(1,2) * (dgHol'_cab' + dgHol'_bac' - dgHol'_abc'))():permute'_abc'
			if self.verbose then
				print()
				print'1st kind Christoffel of holonomic basis:'
				print()
				print(var'\\Gamma'' _\\tilde{a} _\\tilde{b} _\\tilde{c}':eq(
					frac(1,2)*(var''' _\\tilde{a} _\\tilde{b} _,\\tilde{c}'
					+ var'g'' _\\tilde{a} _\\tilde{c} _,\\tilde{b}'
					- var'g'' _\\tilde{b} _\\tilde{c} _,\\tilde{a}')
				):eq(GammaHol_lll' _\\tilde{a} _\\tilde{b} _\\tilde{c}'()))
				print()
			end
			return GammaHol_lll
		end,
		result = 'real3x3s3',
	}

	self.calc.coord_connHol_ull = {
		build = function()
			local gHolU = self.request'coord_gHol_uu'
			local GammaHol_lll = self.request'coord_connHol_lll'
			local GammaHol_ull = (gHolU'^ad' * GammaHol_lll'_dbc')():permute'^a_bc'
			if self.verbose then
				print()
				print'connection:'
				print()
				print(
					var'\\Gamma'' ^\\tilde{a} _\\tilde{b} _\\tilde{c}':eq(
						var'g'' ^\\tilde{a} ^\\tilde{d}'
						* var'\\Gamma'' _\\tilde{d} _\\tilde{b} _\\tilde{c}'
					):eq(
						GammaHol_ull' ^\\tilde{a} _\\tilde{b} _\\tilde{c}'()
					)
				)
				print()
			end
			return GammaHol_ull
		end,
		result = 'real3x3s3',
	}

	-- ConnHol^i = ConnHol^i_jk gHol^jk
	self.calc.coord_connHol_trace23 = {
		build = function()
			local gHolU = self.request'coord_gHol_uu'
			local GammaHol_ull = self.request'coord_connHol_ull'
			return (GammaHol_ull'^a_bc' * gHolU'^bc')()
		end,
		result = 'real3',
	}

	-- gHol_,i
	self.calc.coord_partial_det_gHol_l = {
		build = function()
			return Tensor('_a', function(a)
				return self.symchart.tangentSpaceOperators[a](self.request'coord_det_gHol')()
			end)
		end,
		result = 'real3',
	}

	-- gHol_,ij
	self.calc.coord_partial2_det_gHol_ll = {
		build = function()
			return self.request'coord_partial_det_gHol_l''_a,b'()
		end,
		result = 'real3s3',
	}


	--[[
	args:
		srcname = model to start with, whose build() produces a real3 object
		subprefix = modules to make for each individual component of it
			default = srcname (it appends 0,1,2 onto the submodules)
		dstname = module to make to quick-include all the submodules
			default = srcname..'_i'
		define = true/false whether to use #define or functions
	--]]
	local function addReal3Components(args)
		local srcname = assert(args.srcname)
		local prefix = args.subprefix or srcname
		local dstname = args.dstname or srcname..'_i'

		-- put the individual elements into requests of their own, for codegen request's sake
		for i=1,dim do
			self.calc[prefix..(i-1)] = {
				build = function()
					local t = self.request(srcname)
					local elem = t[i]
					-- TODO assert elem is Expression but not Array?
					return elem
				end,
			}
		end

		self.calc[dstname] = table(
			{many = true},
			range(dim):mapi(function(i)
				local field = prefix..(i-1)
				return {
					field = field,	-- function name & request tensor name
					build = function()
						return self.request(field)
					end,
					result = 'real',
					define = args.define,
					-- TODO also the defines should go in headercode, not code
				}
			end)
		)
	end

	-- put all 'coord_dx#'s into one module called 'coord_dx_i'
	-- put the individual elements of 'coord_dx' request into requests called 'coord_dx#'
	-- dx0, dx1, ...
	-- this is the change in cartesian wrt the change in grid
	-- this is also the normalization factor for the anholonomic ( ... is it?)
	addReal3Components{
		srcname = 'coord_dx',
		define = true,
	}

	-- put individual 'cell_area' into their distinct requests
	-- area0, area1, ...
	-- area_i = integral of u_j, j!=i of product of dx_j, j!=i
	addReal3Components{
		srcname = 'cell_area',
		define = 'with solver arg',
	}

	addReal3Components{
		srcname = 'coord_holBasisLen',
	}

-- [[
	-- scale coord_dx by the solver.grid_dx var to get cell_dx:
	-- TODO use the 'e' and 'eHol' tensors?  then just make this a tensor product?
	--  and then we could use the 'addReal3Components' to build these all at once?
	for i=1,dim do
		self.calc['cell_dx'..(i-1)] = {
			build = function()
				return self.request('coord_dx'..(i-1)) * integralGridDx[i]
			end,
		}
	end

	self.calc.cell_dx_i = table(
		{many=true},
		range(dim):mapi(function(i)
			return {
				field = 'cell_dx'..(i-1),
				build = function()
					return self.request('cell_dx'..(i-1))
				end,
				result = 'real',
				define = 'with solver arg',
			}
		end)
	)
--]]

-- [[
	-- volume of a cell = volume element times grid dx's
	self.calc.cell_sqrt_det_g = {
		build = function()
			local coord_sqrt_det_g = self.request'coord_sqrt_det_g'
			for i=1,dim do
				coord_sqrt_det_g = coord_sqrt_det_g * integralGridDx[i]
			end
			return coord_sqrt_det_g
		end,
		result = 'real',

		-- TODO this  macro had a solver arg, but the other #define's (like cell_dx#) didn't
		-- sooo ... how to specify when to use each?
		define = 'with solver arg',	-- because it uses solver.grid_dx vars

		depends = {solver.solver_t},	-- if you use integralGridDx
	}
--]]


	local function addSym3Components(args)
		local srcname = assert(args.srcname)
		local prefix = args.subprefix or srcname
		local dstname = args.dstname or srcname..'_ij'

		for ij,xij in ipairs(symNames) do -- dim is fixed at 3 so just use symNames?
			local i,j = from6to3x3(ij)
			self.calc[prefix..(i-1)..(j-1)] = {
				build = function()
					local t = self.request(srcname)
					local elem = t[i][j]
					-- TODO assert elem is Expression but not Array?
					return elem
				end,
			}
		end

		self.calc[dstname] = table(
			{many = true},
			range(dim*(dim+1)/2):mapi(function(ij)
				local i,j = from6to3x3(ij)
				local field = prefix..(i-1)..(j-1)
				return {
					field = field,
					build = function()
						return self.request(field)
					end,
					result = 'real',
					define = args.define,
					-- TODO also the defines should go in headercode, not code
				}
			end)
		)
	end

	addSym3Components{
		srcname = 'coord_g_ll',
		define = true,
	}

	addSym3Components{
		srcname = 'coord_g_uu',
		define = true,
	}

	addSym3Components{
		srcname = 'coord_sqrt_g_ll',
		define = true,
	}

	-- curvilinear grid normals use sqrt(g^ii), as it is the metric-weighted coordinate normal magnitude
	addSym3Components{
		srcname = 'coord_sqrt_g_uu',
		define = true,
	}

	addSym3Components{
		srcname = 'coord_sqrt_gHol_ll',
		define = true,
	}
end

-- called after coord creation, before finalize creates the type.  between, other objects can modify it.
function CoordinateSystem:createCellStruct()
	local solver = self.solver

	--[[
	ok here's a dilemma ...
	gridSolver has cellBuf that holds cell pos and any other aux vars used for cell calculations
	meshsolver has cellBuf that holds cell pos and mesh info
	meshsolver needs to pass 'cellBuf'
	--]]
	self.cellStruct = Struct{
		solver = solver,
		name = 'cell_t',
		dontUnion = true,
		vars = {
			{name='pos', type='real3'},		-- x1 x2 x3 input coordinates to the chart
			{name='volume', type='real'},	-- volume of the cell
		},
	}

	-- MeshSolver mesh generation statistics:
	-- cartesian code takes ~20 seconds to compile
	-- cylindrical code takes ~60 seconds.
	-- the only main difference in the code is the # of normal computations ... and at that, directly calling cos() and sin()
	-- I'm going to see if reducing the trig calls helps but giving gridsolvers their own faceBuf
	self.faceStruct = Struct{
		solver = solver,
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
end

function CoordinateSystem:finalizeCellStruct()
	self.cellStruct:makeType()
	self.cell_t = self.cellStruct.typename

	self.faceStruct:makeType()
	self.face_t = self.faceStruct.typename
end

-- TODO this as a CL kernel?
function CoordinateSystem:fillGridCellBuf(cellCpuBuf)
	local solver = self.solver

-- [[ here replace 'solver.' with 'solver.solverPtr.'
	local symmath = require 'symmath'
	local var = symmath.var
	local u, v, w = self.baseCoords:unpack()
	-- TODO instead of passing each function as an arg,
	--  how about adding stuff to the loadstring env in symmath.export.Lua:toFunc ?
	local calcVolume = assert(symmath.export.Lua:toFunc{
		output = {
			self.request'cell_volume'
				:replace(var'solver.grid_dx.x', var'fromreal(solver.solverPtr.grid_dx.x)')
				:replace(var'solver.grid_dx.y', var'fromreal(solver.solverPtr.grid_dx.y)')
				:replace(var'solver.grid_dx.z', var'fromreal(solver.solverPtr.grid_dx.z)')
			,
		},
		input = {{u=u}, {v=v}, {w=w}, var'solver', var'fromreal'},
	})
--]]

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
				cellCpuBuf[index].pos.x = toreal(u)
				cellCpuBuf[index].pos.y = toreal(v)
				cellCpuBuf[index].pos.z = toreal(w)
				cellCpuBuf[index].volume = toreal(calcVolume(u,v,w, solver, fromreal))
				--[[ TODO this is all to evaluate the half-precision inverse of 65520, which I'm failing to do here in luajit (without implementing the bitwise half16 function myself)
				if not math.isfinite(fromreal(cellCpuBuf[index].pos.x)) then
					print("!!! WARNING !!! cell pos.x is non-finite.  Something is bound to go wrong.")
				end
				if not math.isfinite(fromreal(cellCpuBuf[index].pos.y)) then
					print("!!! WARNING !!! cell pos.y is non-finite.  Something is bound to go wrong.")
				end
				if not math.isfinite(fromreal(cellCpuBuf[index].pos.z)) then
					print("!!! WARNING !!! cell pos.z is non-finite.  Something is bound to go wrong.")
				end
				if not math.isfinite(fromreal(cellCpuBuf[index].volume)) then
					print("!!! WARNING !!! cell volume is non-finite.  Something is bound to go wrong.")
				end
				if not math.isfinite(fromreal(cellCpuBuf[index].volume)) then
					print("!!! WARNING !!! cell volume is zero.  Something is bound to go wrong.")
				end
				--]]
				index = index + 1
			end
		end
	end
end

function CoordinateSystem:applyReplVars(expr)
	for _,repl in ipairs(self.repls) do
		expr = expr:subst(repl)
	end
	return expr
end

function CoordinateSystem:applyReplDefines(expr)
	for _,kv in ipairs(self.replDefines) do
		local find, getter = table.unpack(kv)
		local repl = assert(getter(), "getter didn't produce a value")
		expr = expr:replace(find, repl)
	end
	return expr
end

function CoordinateSystem:compile(expr)
	if type(expr) == 'number' then return clnumber(expr) end

	local symmath = require 'symmath'
	local const = symmath.Constant

	expr = self:applyReplVars(expr)

	-- replace pow(x, .5) with sqrt(x)
	expr = expr:map(function(x)
		if symmath.op.pow:isa(x)
		and const:isa(x[2])
		and x[2].value == .5
		then
			return symmath.sqrt(x[1])
		end
	end)

	-- print before replacing base-coordinates
	if self.verbose then
		print('compiling\n', self.fixVerbose(expr)..'\n')
	end
	for i,coord in ipairs(self.baseCoords) do
		expr = expr:replace(coord, symmath.var('pt.'..xNames[i]))
	end
	local code = symmath.export.C(expr)

	return code
end

-- code building functions
local getCode = {}

-- [=====[ as defines.  warning, some of the code calling these functions has some algebra in the arguments (looking at you GLSL code)
-- so in that case, defines can't work until you rewrite the calling code.
-- but lets keep this around as an option for now ...

getCode.real3_to_real_define = function(coord, name, code)
	return template([[
static inline real <?=name?>(real3 pt) {
	return <?=code?>;
}
]], {
		name = name,
		code = code,
	})
end

-- f(x) where x is a point in the coordinate chart
getCode.real3_to_real3_define = function(coord, name, exprs)
	return template([[
static inline real3 <?=name?>(pt) {
	return real3(
<? for i=1,3 do
?>		<?=exprs[i] or '0.'?><?=i==3 and '' or ','?>
<? end
?>	);
}
]], {
		name = name,
		exprs = exprs,
	})
end

-- f(v,x) where x is a point on the coordinate chart and v is most likely a tensor
getCode.real3_real3_to_real_define = function(coord, name, expr)
	return template([[
static inline real3 <?=name?>(real3 u, real3 pt) {
	return <?=expr?>;
}
]], {
		name = name,
		expr = expr,
	})
end

getCode.real3_real3_to_real3_define = function(coord, name, exprs)
	return template([[
static inline real3 <?=name?>(real3 u, real3 pt) {
	return real3( \
<? for i=1,3 do
?>		<?=exprs[i] or '0.'?><?=i==3 and '' or ','?> \
<? end
?>	);
}
]], {
		name = name,
		exprs = exprs,
	})
end

getCode.real3_real3_real3_real3_to_real_define = function(coord, name, expr)
	return template([[
static inline real <?=name?>(real3 u, real3 v, real3 w, real3 pt) {
	return <?=expr?>;
}
]], {
		name = name,
		expr = expr,
	})
end

getCode.real3_real3_real3_to_real3_define = function(coord, name, exprs)
	return template([[
static inline real3 <?=name?>(real3 u, real3 v, real3 pt) {
	return real3( \
<? for i=1,3 do
?>		<?=exprs[i] or '0.'?><?=i==3 and '' or ','?> \
<? end
?>	);
}
]], {
		name = name,
		exprs = exprs,
	})
end


getCode.real3_to_real3s3_define = function(coord, name, exprs)
	return template([[
static inline real3s3 <?=name?>(real3 pt) {
	return real3s3{
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		<?=exprs[i] and exprs[i][j] and exprs[i][j] or '0.'?>,
<? end
?>	};
}
]], {
		symNames = symNames,
		from6to3x3 = from6to3x3,
		name = name,
		exprs = exprs,
	})
end

-- symmetric on 2nd & 3rd indexes
getCode.real3_to_real3x3s3_define = function(coord, name, exprs)
	return template([[
static inline real3x3s3 <?=name?>(real3 pt) {
	return real3x3s3{
<?
for i,xi in ipairs(xNames) do
?>	.<?=xi?> = {
<?	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>		.<?=xjk?> = <?=exprs[i] and exprs[i][j] and exprs[i][j][k]
			and exprs[i][j][k] or '0.'?>,
<?	end
?>	},
<?
end
?>	}
}
]], {
		name = name,
		exprs = exprs,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
	})
end

-- ugly hack, TODO, just use one interface for all
getCode.real3_to_real_defineWithSolver = function(coord, name, code)
	return template([[
static inline real <?=name?>(
	constant <?=solver_t?> const & solver,
	real3 pt
) {
	return <?=code?>;
}
]], {
		name = name,
		code = code,
		solver_t = coord.solver.solver_t,
	})
end

--]=====]
-- [=====[ as functions:

getCode.real3_to_real = function(coord, name, code)
	return template([[
static inline real <?=name?>(real3 pt) {
	return <?=code?>;
}]], {
		name = name,
		code = code,
	})
end

-- f(x) where x is a point in the coordinate chart
getCode.real3_to_real3 = function(coord, name, exprs)
	return template([[
static inline real3 <?=name?>(real3 pt) {
	return real3(
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
getCode.real3_real3_to_real = function(coord, name, expr)
	return template([[
static inline real <?=name?>(real3 u, real3 pt) {
	return <?=expr?>;
}]], {
		name = name,
		expr = expr,
	})
end

getCode.real3_real3_to_real3 = function(coord, name, exprs)
	return template([[
static inline real3 <?=name?>(real3 u, real3 pt) {
	return real3(
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

getCode.real3_real3_real3_real3_to_real = function(coord, name, expr)
	return template([[
static inline real <?=name?>(real3 u, real3 v, real3 w, real3 pt) {
	return <?=expr?>;
}]], {
		name = name,
		expr = expr,
	})
end

getCode.real3_real3_real3_to_real3 = function(coord, name, exprs)
	return template([[
static inline real3 <?=name?>(real3 u, real3 v, real3 pt) {
	return real3(
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


getCode.real3_to_real3s3 = function(coord, name, exprs)
	return template([[
real3s3 <?=name?>(real3 pt) {
	return real3s3{
<? for ij,xij in ipairs(symNames) do
	local i,j,xi,xj = from6to3x3(ij)
?>		<?=exprs[i] and exprs[i][j] and exprs[i][j] or '0.'?>,
<? end
?>	};
}]], {
		symNames = symNames,
		from6to3x3 = from6to3x3,
		name = name,
		exprs = exprs,
	})
end

-- symmetric on 2nd & 3rd indexes
getCode.real3_to_real3x3s3 = function(coord, name, exprs)
	return template([[
real3x3s3 <?=name?>(real3 pt) {
	return real3x3s3{
<?
for i,xi in ipairs(xNames) do
?>	{
<?	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
?>		<?=exprs[i] and exprs[i][j] and exprs[i][j][k]
			and exprs[i][j][k] or '0.'?>,
<?	end
?>	},
<?
end
?>	};
}]], {
		name = name,
		exprs = exprs,
		xNames = xNames,
		symNames = symNames,
		from6to3x3 = from6to3x3,
	})
end

-- this is my exception to the rule, which accepts a pointer
getCode.real3_toreal3x3s3x3 = function(coord, name, exprs)
	return template([[
void <?=name?>(real3x3s3 a[3], real3 pt) {
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
		name = name,
		exprs = exprs,
		symNames = symNames,
		from6to3x3 = from6to3x3,
	})
end

--]=====]


--[[
ok standardizing these macros ...
we have a few manifolds to deal with ...
1) grid coordiantes, I = 0..size-1
2) manifold coordinates, X = grid coordinates rescaled to fit between solver.mins and solver.maxs
3) result coordinates, Y = cartesian output of the mapping of the coordinate system

coord_dx[_for_coord]<?=side?>(pt) gives the length of the dx of the coordinate in 'side' direction, for 'pt' in coordinates
cell_dx[_for_coord]<?=side?>(pt) gives the length of the dx of the cell, which is just coord_dx times the grid_dx
solver.grid_dx is the size of each solver cell, in coordintes

should I add these _for_coord _for_grid suffixes to specify what manfiold system the input parameter is?

--]]
function CoordinateSystem:initCodeModules()
	local solver = self.solver

	for moduleName,infos in pairs(self.calc) do
		if not infos.many then infos = {infos} end
		local depends = table()
		local buildsAndExprNames = table()
		for _,info in ipairs(infos) do
			if not info.result then
				--print("can't add module for "..moduleName.." which is missing its calc result")
				-- should I warn here?  clutters the output
				-- or should I warn if someone tries to add the module when it's not there?  I guess that's already happening as an error.
				-- or should I just make 'real' the default return type and generate all modules?
			else
				local name = info.field or moduleName
				local buildNameParts = table{'real3','to', info.result}
				if info.args then
					buildNameParts:insert(2, info.args)
				end
				if info.define == 'with solver arg' then
					buildNameParts:insert'defineWithSolver'
				elseif info.define then
					buildNameParts:insert'define'
				end
				local buildName = buildNameParts:concat'_'
				local build = assert(getCode[buildName])
				-- don't include anything that's typedef'd by app
				-- otherwise it will screw up the separation between C, CL, and GL
				if info.result ~= 'real'
				and info.result ~= 'real3'
				then
					depends:insert(info.result)
				end
				if info.depends then
					depends:append(info.depends)
				end
				buildsAndExprNames:insert{build=build, name=name}
			end
		end
		if #buildsAndExprNames > 0 then
			xpcall(function()
				solver.modules:add{
					name = self.symbols[moduleName] or error("failed to find symbol for coord-depend "..moduleName),
					depends = depends,
					code = function()
						return buildsAndExprNames:mapi(function(buildAndName)
							local build = buildAndName.build
							local name = buildAndName.name
							return build(self, name, self.compilePrintRequestTensor(name))
						end):concat'\n'
					end,
				}
			end, function(err)
				io.stderr:write('failed for module: '..moduleName..'\n'..err..'\n'..debug.traceback())
				os.exit(1)
			end)
		end
	end


	self:initCodeModule_coordMap()

	-- parallel propagate code
	if require 'hydro.solver.fvsolver':isa(solver)
	-- TODO only if it's a mesh solver using a flux integrator ... which is currently all mesh solvers
	or require 'hydro.solver.meshsolver':isa(solver)
	then
		local lines = table()

		-- parallel-propagate a vector from point 'x' along coordinate 'k' (suffix of func name) by amount 'dx'
		-- TODO derive this from the metric
		-- upper parallel propagator = exp(-int_x^(x+dx) of conn_side dx)
		-- lower parallel propagator = exp(int_x^(x+dx) of conn_side^T dx) = upper^-T
		--
		-- if we're using a cartesian basis then no need to transport anything
		-- ... (? except maybe the flux differential across the cell, done in curvilinear coordinates ?)
		if self.vectorComponent == 'cartesian'
		or require 'hydro.coord.cartesian':isa(self)
		then
			-- general case for a fixed global orthonormal basi:
			lines:insert(template([[

// Propagate an upper index / propagate the components of a vector.
<? for side=0,solver.dim-1 do ?>
real3 coord_parallelPropagateU<?=side?>(real3 v, real3 x, real dx) {
	return v;
}
<? end ?>

// To propagate a one-form, apply the inverse-transpose.
// In the case of anholonomic orthonormal basis, the transpose equals the inverse so propagateL == propagateR.
//  This fits with the fact that the metric of an orthonormal basis is identity, so the upper and lower components are equal.
<? for side=0,solver.dim-1 do ?>
real3 coord_parallelPropagateL<?=side?>(real3 v, real3 x, real dx) {
	return v;
}
<? end ?>

]], 		{
				solver = solver,
			}))

		-- anholonomic orthonormal = metric is identity, so inverse = transpose, so upper propagate = lower propagate ... so we only need to define the upper
		elseif self.vectorComponent == 'anholonomic' then	-- TODO rename this to 'orthonormal' (since anholonomic doesn't imply orthonormal ... and this assumption is only true for orthonormal basis)
			lines:insert(self:getParallelPropagatorCode())
			lines:insert(template([[
<? for side=0,solver.dim-1 do?>
real3 coord_parallelPropagateL<?=side?>(real3 v, real3 x, real dx) {
	return coord_parallelPropagateU<?=side?>(v, x, dx);
}
<? end ?>
]], 		{
				solver = solver,
			}))
		elseif self.vectorComponent == 'holonomic' then	-- TODO rename this to 'coordinate'
			lines:insert(self:getParallelPropagatorCode())
		else
			error'here'
		end

		solver.modules:add{
			name = self.symbols.coord_parallelPropagate,
			depends = self:getModuleDepends_coord_parallelPropagate(),
			code = lines:concat'\n',
		}
	end

	self:initCodeModule_normal()

	-- store all input coordinates for our cells
	-- for holonomic/anholonomic this is just the linearly interpolated
	solver.modules:add{
		name = self.cell_t,
		structs = {self.cellStruct:getForModules()},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.cell_t..' cell_t;',
	}

	solver.modules:add{
		name = self.face_t,
		structs = {self.faceStruct:getForModules()},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.face_t..' face_t;',
	}

	solver.modules:add{
		name = self.symbols.cell_calcAvg_withPt,
		depends = {self.cell_t},
		code = function()
print("WARNING - haven't finished implementing cell_calcAvg_withPt")
			return self.solver.eqn:template[[
#define cell_calcAvg_withPt(\
	/*<?=cell_t?> * const */resultCell,\
	/*<?=cell_t?> const * const */cellL,\
	/*<?=cell_t?> const * const */cellR,\
	/*real3 const */pt\
) {\
	/* TODO average any other fields here .... */\
	(resultCell)->pos = pt;\
}
]]
		end,
	}

	solver.modules:add{
		name = self.symbols.cell_dxs,
		depends = {self.symbols.cell_dx_i},
		code = solver.eqn:template[[
// TODO std::function and std::tuple ...
/*
// error: 'cell_areas' declared as array of 'auto'
// error: taking address of function is not allowed
auto cell_dxs[3] = {
	cell_dx0,
	cell_dx1,
	cell_dx2,
};
*/

// until then, here's templated:
template<int side> real cell_dxs(constant solver_t_1 const & solver, real3 pt);
template<> real cell_dxs<0>(constant solver_t_1 const & solver, real3 pt) { return cell_dx0(solver, pt); }
template<> real cell_dxs<1>(constant solver_t_1 const & solver, real3 pt) { return cell_dx1(solver, pt); }
template<> real cell_dxs<2>(constant solver_t_1 const & solver, real3 pt) { return cell_dx2(solver, pt); }

real3 cell_dx_vec(
	constant solver_t_1 const & solver,
	real3 const pt
) {
	return real3(
		cell_dxs<0>(solver, pt),
		cell_dxs<1>(solver, pt),
		cell_dxs<2>(solver, pt)
	);
}
]],
	}

	solver.modules:add{
		name = self.symbols.cell_areas,
		depends = {self.symbols.cell_area_i},
		code = solver.eqn:template[[
/*
auto cell_areas[3] = {
	cell_area0,
	cell_area1,
	cell_area2,
};
*/

template<int side> real cell_areas(constant solver_t_1 const & solver, real3 pt);
template<> real cell_areas<0>(constant solver_t_1 const & solver, real3 pt) { return cell_area0(solver, pt); }
template<> real cell_areas<1>(constant solver_t_1 const & solver, real3 pt) { return cell_area1(solver, pt); }
template<> real cell_areas<2>(constant solver_t_1 const & solver, real3 pt) { return cell_area2(solver, pt); }
]],
	}

	solver.modules:add{
		name = self.symbols.coord_basisHolUnits,
		depends = {self.symbols.coord_basisHolUnit_i},
		code = solver.eqn:template[[
/*
auto coord_basisHolUnits[3] = {
	coord_basisHolUnit0,
	coord_basisHolUnit1,
	coord_basisHolUnit2,
};
*/

template<int side> real3 coord_basisHolUnits(real3 pt);
template<> real3 coord_basisHolUnits<0>(real3 pt) { return coord_basisHolUnit0(pt); }
template<> real3 coord_basisHolUnits<1>(real3 pt) { return coord_basisHolUnit1(pt); }
template<> real3 coord_basisHolUnits<2>(real3 pt) { return coord_basisHolUnit2(pt); }
]],
	}



end

function CoordinateSystem:getModuleDepends_coordMap()
end
function CoordinateSystem:getModuleDepends_coordMapInv()
end

function CoordinateSystem:getModuleDepends_coordMapGLSL()
	local glslDeps = table()
--[[
	-- search for any functions in the expression
	-- and auto-insert them into getModuleDepends_coordMap
	-- then for those functions in cl but not glsl,
	-- insert those there

	do
		-- this is the expression initCodeModule_coordMap uses
		local symmath = require 'symmath'
		local expr = self.request'u'
		for _,f in ipairs{
			-- here are functions builtin for cl but not glsl
			...
		} do
			if expr:findLambda(function(x)
				return f:isa(x)
			end) then
				glslDeps:insert(f.name)
			end
		end
	end
--]]
	return table(self:getModuleDepends_coordMap()):append(glslDeps)
end
function CoordinateSystem:getModuleDepends_coordMapInvGLSL()
	return self:getModuleDepends_coordMapInv()
end

function CoordinateSystem:getModuleDepends_coord_parallelPropagate()
end

function CoordinateSystem:initCodeModule_coordMap()
	local solver = self.solver

	-- get back the Cartesian coordinate for some provided chart coordinates
	-- TODO make this a macro based on cellBuf[index]
	-- and make it custom per coord system (just like the cellBuf fields are)
	local uCode = self.compilePrintRequestTensor'u'
	local code_coordMap = getCode.real3_to_real3(self, 'coordMap', range(3):mapi(function(i)
		return uCode[i] or 'pt.'..xNames[i]
	end))

	solver.modules:add{
		name = self.symbols.coordMap,
		depends = self:getModuleDepends_coordMap(),
		code = function()
			return code_coordMap
		end,
	}

	-- get back the radial distance for some provided chart coordinates
	solver.modules:add{
		name = self.symbols.coordMapR,
		code = function()
			return getCode.real3_to_real(self, 'coordMapR', self:compile(self.vars.r))
		end,
	}

	local code_coordMapInv = self:getCoordMapInvModuleCode()	-- until i can autogen this ...

	solver.modules:add{
		name = self.symbols.coordMapInv,
		depends = self:getModuleDepends_coordMapInv(),
		code = function()
			return code_coordMapInv
		end,
	}


	-- GLSL needs some extra depends
	--  and I can't think of how to add them in except by doing this ...
	-- see my rant in SphereLogRadial:getModuleDepends_coordMap
	solver.modules:add{
		name = self.symbols.coordMapGLSL,
		depends = self:getModuleDepends_coordMapGLSL(),
		code = function()
			return code_coordMap
		end,
	}

	solver.modules:add{
		name = self.symbols.coordMapInvGLSL,
		depends = self:getModuleDepends_coordMapInvGLSL(),
		code = function()
			return code_coordMapInv
		end,
	}


	solver.modules:add{
		name = self.symbols.coord_basis_i,
--		depends = {'real3'},	-- modules can't explicitly include 'real3' and be used with GLSL at the same time ...
		code = function()
			local eExt = self.compilePrintRequestTensor'eExt'
			return eExt:mapi(function(eiCode,i)
				return getCode.real3_to_real3(self, 'coordBasis'..(i-1), eiCode)
			end):concat'\n'
		end,
	}

	solver.modules:add{
		name = self.symbols.coord_basisHolUnit_i,
		code = function()
			local eHolUnitCode = self.compilePrintRequestTensor'eHolUnitExt'
			return eHolUnitCode:mapi(function(eHolUnitiCode,i)
				return getCode.real3_to_real3(self, 'coord_basisHolUnit'..(i-1), eHolUnitiCode)
			end):concat'\n'
		end,
	}

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
		or require 'hydro.coord.cartesian':isa(self)
		then
			if not require 'hydro.coord.cartesian':isa(self) then
				depends:insert(self.symbols.coord_basisHolUnit_i)
				tolines:insert(template([[
//converts a vector from cartesian coordinates to grid curvilinear coordinates
//by projecting the vector into the grid basis vectors
//at x, which is in grid curvilinear coordinates
real3 coord_cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3(0,0,0);	//for gl compat
	<? for i=0,solver.dim-1 do
	local xi = xNames[i+1]
	?>{
		real3 e = coord_basisHolUnit<?=i?>(pt);
<? if coord.vectorComponent == 'anholonomic' then	-- anholonomic normalized
?>		real uei = e.dot(u);
<? else		-- holonomic
?>		real uei = e.dot(u) / e.lenSq();
<? end
?>		uCoord.<?=xi?> = uei;
		//subtract off this basis component from u
		u -= e * uei;
	}<? end ?>
	//add whatever's left of u
	return uCoord + u;
}
]], env))
				fromlines:insert(template([[
//converts a vector from cartesian to grid curvilinear coordinates
//by projecting it onto the basis ... ?
real3 coord_cartesianFromCoord(real3 u, real3 pt) {
	real3 uGrid = real3(0,0,0);	//for gl compat
	<? for i=0,solver.dim-1 do
	local xi = xNames[i+1]
	?>{
		real3 e = coord_basisHolUnit<?=i?>(pt);
		uGrid += e * u.<?=xi?>;
	}<? end ?>
	return uGrid;
}
]], env))

			else	-- cartesian:isa(coord)
				fromlines:insert[[
// TODO for when I uncouple this from GLSL
// ... or make the shader code SPIR-V from OpenCL-C++ also
real3 coord_cartesianFromCoord(real3 u, real3 pt) { return u; }
]]
				tolines:insert[[
// TODO for when I uncouple this from GLSL
// ... or make the shader code SPIR-V from OpenCL-C++ also
real3 coord_cartesianToCoord(real3 u, real3 pt) { return u; }
]]
			end

--TODO change names
--'coord' is ambiguous
-- this is relative to teh vector component basis
			fromlines:insert[[
// TODO for when I uncouple this from GLSL
// ... or make the shader code SPIR-V from OpenCL-C++ also
real3 cartesianFromCoord(real3 u, real3 private const pt) { return u; }
]]
			tolines:insert[[
// TODO for when I uncouple this from GLSL
// ... or make the shader code SPIR-V from OpenCL-C++ also
real3 cartesianToCoord(real3 u, real3 pt) { return u; }
]]

		else	-- coord.vectorComponent

			depends:insert(self.symbols.coord_basis_i)
			tolines:insert(template([[
//converts a vector from cartesian coordinates to grid curvilinear coordinates
//by projecting the vector into the grid basis vectors
//at x, which is in grid curvilinear coordinates
real3 coord_cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3(0,0,0);	//for gl compat
	<? for i=0,solver.dim-1 do
	local xi = xNames[i+1]
	?>{
		real3 e = coordBasis<?=i?>(pt);
<? if coord.vectorComponent == 'anholonomic' then	-- anholonomic normalized
?>		real uei = e.dot(u);
<? else		-- holonomic
?>		real uei = e.dot(u) / e.lenSq();
<? end
?>		uCoord.<?=xi?> = uei;
		//subtract off this basis component from u
		u -= e * uei;
	}<? end ?>
	//add whatever's left of u
	return uCoord + u;
}

real3 cartesianToCoord(real3 u, real3 pt) { return coord_cartesianToCoord(u, pt); }
]], env))
			fromlines:insert(template([[
//converts a vector from cartesian to grid curvilinear coordinates
//by projecting it onto the basis ... ?
real3 coord_cartesianFromCoord(real3 u, real3 pt) {
	real3 uGrid = real3(0,0,0);	//for gl compat
	<? for i=0,solver.dim-1 do
	local xi = xNames[i+1]
	?>{
		uGrid += coordBasis<?=i?>(pt) * u.<?=xi?>;
	}<? end ?>
	return uGrid;
}

real3 cartesianFromCoord(real3 u, real3 pt) { return coord_cartesianFromCoord(u, pt); }
]], env))

		end -- coord.vectorComponent

		solver.modules:add{
			name = self.symbols.cartesianFromCoord,
			depends = depends,
			code = fromlines:concat'\n',
		}
		solver.modules:add{
			name = self.symbols.cartesianToCoord,
			depends = depends,
			code = tolines:concat'\n',
		}
	end
end

-- until I get inverses on systems of equations working, I'll have this manually specified
function CoordinateSystem:getCoordMapInvModuleCode()
	return [[
real3 coordMapInv(real3 x) { return x; }
]]
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
function CoordinateSystem:initCodeModule_normal()
	local typecode	-- only for lua ffi, inserted into c++ so i can static_assert sizeof ==
	local code = self.solver.eqn:template(path'hydro/coord/coord.clcpp':read())

	if self.verbose then
		print[[
normals:<br>
$n^i =$ i'th normal, along the i'th coordinate basis.<br>
$(n^i)_\\hat{j}$ = i'th basis direction, j'th component.<br>
]]
	end
	if require 'hydro.solver.meshsolver':isa(self.solver) then
		if self.verbose then
			print[[
mesh-based normals are assuming a cartesian coordinate system, and assumed to be normalized:<br>
$(n^i)_k (n^j)^k = g^{ij} = \delta^{ij}$<br>
]]
		end

		typecode = self.solver.eqn:template[[
typedef struct {
	real3x3 n;
} <?=normal_t?>;
]]

		code = code .. self.solver.eqn:template[[
namespace <?=Solver?> {
using Normal = MeshNormal;
}	//namespace <?=Solver?>
]]
	else	-- not meshsolver

		if require 'hydro.coord.cartesian':isa(self)
		or self.vectorComponent == 'anholonomic'
		then
			if self.verbose then
				print[[
anholonomic normals are coordinate-aligned but orthonormalized by the locally-Cartesian basis:<br>
$(n^i)^k (n^j)^l g_{kl} = g^{ij} = \delta^{ij}$, for anholonomic (orthonormal) basis.<br>
$(n^i)_j = (n^i)^j = \delta^{ij}$.<br>
]]
			end

			--[[
			n_i = n^i = delta_ij for side j
			|n| = 1
			--]]
			-- typecode is now only used for lua ffi cdef and tested to sizeof the c++ struct with static_assert
			typecode = self.solver.eqn:template[[
typedef struct {
	int side;			//0, 1, 2
} <?=normal_t?>;		//nL = nU = normalBasisForSide (permutation of I), nLen = 1
]]

			-- this is overwhelmingly interface-based
			-- with two exceptions: calcDT and some displayVars
			code = code .. self.solver.eqn:template[[
namespace <?=Solver?> {
using Normal = AnholonomicNormal;
}	// namespace <?=Solver?>
]]
		elseif self.vectorComponent == 'cartesian' then
			if self.verbose then
				print[[
cartesian-component normals:<br>
$n_i = e_i$<br>
$(n_i)_j = e_{ij}$<br>
$(n_i)^j = {e_i}^j$<br>
$n_i \cdot n_j = (n_i)_k (n_j)^l = \delta_{ij}$<br>
]]
			end

			--[[
			n_i = n^i = unit(u^i_,j) for side j
			|n| = sqrt(n^i n_i) = 1 (since g_ij = g^ij = delta_ij)
			--]]
			typecode = self.solver.eqn:template[[
typedef struct {
	real3x3 n;		// nL = nU, both are orthonormal so nLen = 1
	real len;
} <?=normal_t?>;
]]

			-- this would call coord_cartesianFromCoord
			-- which itself aligns with the coord_basisHolUnit
			code = code .. self.solver.eqn:template[[
namespace <?=Solver?> {
using Normal = CartesianNormal;
}	// namespace <?=Solver?>
]]
		elseif self.vectorComponent == 'holonomic' then
			if self.verbose then
				print[[
cartesian normals are aligned to the embedding cartesian coordinate space:<br>
$n_i = e_i$<br>
$(n_i)_j = \delta_{ij}$<br>
$(n_i)^j = g^{jk} \delta_{ik}$<br>
$(n_i)^k (n_j)^l = g^{km} \delta_{mi} g^{ln} \delta_{nj}$<br>
]]
			end


			--[[
			n_i = delta_ij for side j
			n^i = g^ik delta_kj
			|n| = sqrt(n^i n_i) = sqrt(g^jj)
			--]]
			typecode = self.solver.eqn:template[[
typedef struct {
	int side;
	real3x3 U;
	real len;
} <?=normal_t?>;
]]

			code = code .. self.solver.eqn:template[[
namespace <?=Solver?> {
using Normal = HolonomicNormal;
}	// namespace <?=Solver?>
]]
		else
			error("unknown vectorComponent == "..tolua(self.vectorComponent))
		end
	end	-- meshsolver

	code = code .. self.solver.eqn:template[[
static_assert(sizeof(<?=Solver?>::Normal) == sizeof(<?=normal_t?>));
]]

	-- TODO if you use multiple solvers that have differing vectorComponents
	--  then this will cause a silent ffi error.  only the first normal_t will be defined.
	-- solution: rename all the normal_t C types to <?=normal_t?>
	self.solver.modules:add{
		name = self.symbols.normal_t,
		depends = {
			'real3x3',
			self.symbols.coord_basisHolUnits,
			self.symbols.coord_g_uu_ij,
			self.symbols.coord_sqrt_g_uu_ij,
		},
		typecode = typecode,
		code = code,
	}
end


return CoordinateSystem
