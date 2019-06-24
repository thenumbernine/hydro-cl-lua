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

(and after reading the coordinate-invariant BSSN papers0
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
	-- TODO do we?  finite volume uses the holonomic metric determinant derivative, which is ident -> 1 -> zero
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

	-- for output
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
		code = code:gsub('{w^(%d)}', function(i)
			return 'w'..self.coords[i+0]
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
	local function compileTensor(expr)
		if symmath.Array.is(expr) then
			return table.mapi(expr, compileTensor)
		elseif symmath.Expression.is(expr) then
			return compile(expr)
		elseif type(expr) == 'number' then
			return clnumber(expr)
		else
			error("don't know how to compile "
				..symmath.Verbose(expr))
				--..require 'ext.tolua'(expr))
		end
	end
	--[[
	local function compileTensorField(field)
		local codeField = field..'Code'
		self[codeField] = compileTensor(self[field])
		printNonZeroField(codeField)
	end
	--]]
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
		return e[j][i] or const(0)
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

	compileTensorField('dg_lll_codes', dg)

	compileTensorField('d2g_llll_codes', d2g)
	
	compileTensorField('conn_lll_codes', Gamma_lll)

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

	-- sqrt(g)_,i / sqrt(g) = Conn^j_ij
	local connExpr = (Gamma_ull'^b_ab')()
	compileTensorField('tr13_conn_l_codes', connExpr)

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

	local areaExprs = symmath.Array:lambda({dim}, function(i)
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
	compileTensorField('areaCodes', areaExprs)

	self.g = g
	compileTensorField('gCode', self.g)

	local gU = Tensor('^ab', table.unpack((Matrix.inverse(g))))
	self.gU = gU
	compileTensorField('gUCode', self.gU)
	
	local sqrt_gU = Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
	compileTensorField('sqrt_gUCode', sqrt_gU)

	local det_g_expr = symmath.Matrix.determinant(g)
	compileTensorField('det_g_code', det_g_expr)

	local sqrt_det_gExpr = symmath.sqrt(det_g_expr)()
	compileTensorField('sqrt_det_gCode', sqrt_det_gExpr)


	-- Now, for num rel in spherical, I need to transform the metric to remove the singularities ...
	-- It is funny to include this alongside the anholonomic stuff.  It seems you would need one or the other but not both.
	
	local J = Tensor('^a_b', function(a,b)
		return a == b and (1 / self.lenExprs[a]) or 0
--[[
print('i',i,'j',j,self.lenExprs[i])		
		return i ~= j 
			and 0
			or (self.anholonomic 
					and 1 
					or (1 / self.lenExprs[i])
				)
--]]	
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

	-- metric determinant ...  det_g = volume^2 for holonomic basis
	local det_g_code = '(' .. self.det_g_code .. ')'
	lines:insert(getCode_real3_to_real('coord_det_g', det_g_code))

	-- sqrt_det_g ... volume for holonomic basis
	local sqrt_det_gCode = '(' .. self.sqrt_det_gCode .. ')'
	lines:insert(getCode_real3_to_real('coord_sqrt_det_g', sqrt_det_gCode))
	
	-- coord len code: l(v) = v^i v^j g_ij
	lines:append{
		getCode_real3_real3_to_real('coordLenSq', self.uLenSqCode),
		[[
real coordLen(real3 r, real3 pt) {
	return sqrt(coordLenSq(r, pt));
}]],
	}

	lines:insert(getCode_real3_real3_real3_to_real('coord_conn_apply123', self.connApply123Code))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply12', self.connApply12Codes))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply13', self.connApply13Codes))
	lines:insert(getCode_real3_real3_real3_to_real3('coord_conn_apply23', self.connApply23Codes))
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
		.y = v.y / coord_dx1(x),
		.z = v.z / coord_dx2(x),
	};
}
#define real3_rescaleToCoord_u real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y * coord_dx1(x),
		.z = v.z * coord_dx2(x),
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_l

sym3 sym3_rescaleFromCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleToCoord_uu sym3_rescaleFromCoord_ll

sym3 sym3_rescaleToCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
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
?>		real uei = real3_dot(e, u); // / real3_len(e);
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
