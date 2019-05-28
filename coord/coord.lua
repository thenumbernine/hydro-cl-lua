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

local CoordinateSystem = class()

function CoordinateSystem:init(args)
	local symmath = require 'symmath'
	local const = symmath.Constant


	-- debugging
	local dprint = args.solver.app.verbose and print or function() end



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

--dprint('coordinates:', table.unpack(coords))
--dprint('base coords:', table.unpack(baseCoords))
--dprint('embedding:', table.unpack(embedded))
	
	local eta = Tensor('_IJ', table.unpack(flatMetric)) 
--dprint'flat metric:'
--dprint(var'\\eta''_IJ':eq(eta'_IJ'()))
--dprint()

	local u = args.chart()
dprint'coordinate chart:'
dprint(var'u''^I':eq(u'^I'()))
dprint()

	local e = Tensor'_u^I'
	e['_u^I'] = u'^I_,u'()
dprint'embedded:'
dprint(var'e''_u^I':eq(var'u''^I_,u'):eq(e'_u^I'()))
dprint()

	local anholonomic
	for i=1,#coords do
		if coords[i].base then
			anholonomic = true
			break
		end
	end
	self.anholonomic = anholonomic
dprint('is anholonomic? '..tostring(anholonomic))

	-- for the sake of grid lengths, 
	-- I will need the basis and metric of the holonomic version as well
	local eHol
	if not anholonomic then
		eHol = e
	else
		eHol = Tensor('_u^I', function(a,I)
			return u[I]:diff(baseCoords[a])()
		end)
dprint'holonomic embedded:'
dprint(var'eHol''_u^I':eq(var'u''^I_,u'):eq(eHol'_u^I'()))
	end

	-- commutation coefficients
	local c = Tensor'_ab^c'
	if anholonomic then
--dprint'connection coefficients:'
--dprint(var'c''_uv^w' * var'e''_w','$=[ e_u, e_v ]$')
		for i,ui in ipairs(coords) do
			for j,uj in ipairs(coords) do
				local psi = var('\\psi', baseCoords)
				local diff = ui:applyDiff(uj:applyDiff(psi)) - uj:applyDiff(ui:applyDiff(psi))
				local diffEval = diff()
				if diffEval ~= const(0) then
--dprint('$[',ui.name,',',uj.name,'] =$',diff:eq(diffEval))
					diff = diff()
--dprint('factor division',diff)
					local dpsi = table.mapi(baseCoords, function(uk) return psi:diff(uk) end)
--dprint('dpsi', dpsi:unpack())
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
dprint'commutation:'
dprint(var'c''_uv^w':eq(c'_uv^w'()))

	local g = (e'_u^I' * e'_v^J' * eta'_IJ')()

dprint'metric:'
dprint(var'g''_uv':eq(var'e''_u^I' * var'e''_v^J' * var'\\eta''_IJ'):eq(g'_uv'()))

	local gHol
	if not anholonomic then
		gHol = g
	else
		gHol = (eHol'_u^I' * eHol'_v^J' * eta'_IJ')()
dprint'holonomic metric:'
dprint(var'gHol''_uv':eq(var'eHol''_u^I' * var'eHol''_v^J' * var'\\eta''_IJ'):eq(gHol'_uv'()))
	end	

	Tensor.metric(g)

	local Gamma_lll = Tensor'_abc'
	Gamma_lll['_abc'] = ((g'_ab,c' + g'_ac,b' - g'_bc,a' + c'_abc' + c'_acb' - c'_bca') / 2)()
dprint'1st kind Christoffel:'
dprint(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(Gamma_lll'_abc'()))

	local Gamma = Tensor'^a_bc'
	Gamma['^a_bc'] = Gamma_lll'^a_bc'()
dprint'connection:'
dprint(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma'^a_bc'()))


	-- for anholonomic coordinates, we also need the holonomic connections
	--  for calculating the sqrt(det(g)) = grid volume
	if self.anholonomic then
		Tensor.metric(gHol)
		
		local GammaHol_lll = Tensor'_abc'
		GammaHol_lll['_abc'] = ((gHol'_ab,c' + gHol'_ac,b' - gHol'_bc,a') / 2)()
dprint'1st kind Christoffel of holonoic basis / Levi-Civita connection:'
dprint(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'gHol''_ab,c' + var'gHol''_ac,b' - var'gHol''_bc,a')):eq(GammaHol_lll'_abc'()))
	
		local GammaHol_ull = Tensor'^a_bc'
		GammaHol_ull['^a_bc'] = GammaHol_lll'^a_bc'()
dprint'holonomic / Levi-Civita connection:'
dprint(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(GammaHol_ull'^a_bc'()))
		
		Tensor.metric(g)
	end


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
	
--dprint'compiling...'
--dprint(orig)
--dprint'...to...'
--dprint(code)

		return code
	end

	-- uCode is used to project the grid for displaying
	self.uCode = range(dim):mapi(function(i) 
		local uCode = compile(u[i])
if uCode ~= '0.' then
	dprint('uCode['..i..'] = '..substCoords(uCode))
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
if eijCode ~= '0.' then
	dprint('eCode['..i..']['..j..'] = '..substCoords(eijCode))
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
dprint('eHolLen['..i..'] = '..substCoords(eiHolLenCode))
		return eiHolLenCode
	end)
	
	local eExtLen = eExt:mapi(function(ei,i)
		return symmath.sqrt(ei:mapi(function(x) return x^2 end):sum())()
	end)
	local eExtUnit = eExt:mapi(function(ei,i)
		return ei:mapi(function(eij) return (eij/eExtLen[i])() end)
	end)
	self.eUnitCode = eExtUnit:mapi(function(ei_unit,i) return ei_unit:mapi(compile) end)
dprint('eUnitCode = ', tolua(self.eUnitCode, {indent=true}))
--]=]

	-- v^k -> v_k
	local lowerExpr = paramU'_a'()
	self.lowerCodes = range(dim):mapi(function(i)
		local lowerCode = compile(lowerExpr[i])
dprint('lowerCode['..i..'] = '..substCoords(lowerCode))
		return lowerCode
	end)

	-- v^k -> v_k
	local raiseExpr = paramU'_a'()
	self.raiseCodes = range(dim):mapi(function(i)
		local raiseCode = compile(raiseExpr[i])
dprint('raiseCode['..i..'] = '..substCoords(raiseCode))
		return raiseCode
	end)

	-- v^k v_k
	local lenSqExpr = (paramU'^a' * paramU'_a')()
	self.uLenSqCode = compile(lenSqExpr)
dprint('uLenSqCodes = '..substCoords(self.uLenSqCode))
	
	self.connCodes = range(dim):mapi(function(i)
		return range(dim):mapi(function(j)
			return range(dim):mapi(function(k)
				local code = compile(Gamma[i][j][k])
if code ~= '0.' then
	dprint('connCode['..i..j..k..'] = '..substCoords(code))
end
				return code
			end)
		end)
	end)
	
	-- Conn^i_jk(x) u^j v^k
	local connExpr = (Gamma'^a_bc' * paramU'^b' * paramV'^c')()
	self.connApply23Codes = range(dim):mapi(function(i)
		local conniCode = compile(connExpr[i])
if conniCode ~= '0.' then
	dprint('connApply23Code['..i..'] = '..substCoords(conniCode))
end		
		return conniCode
	end)

	-- u^j v^k Conn_jk^i(x)
	local connLastExpr = (paramU'^b' * paramV'^c' * Gamma'_bc^a')()
	self.connApply12Codes = range(dim):mapi(function(i)
		local connLastiCode = compile(connLastExpr[i])
if connLastiCode ~= '0.' then
	dprint('connApply12Code['..i..'] = '..substCoords(connLastiCode))
end		
		return connLastiCode
	end)
	
	-- Conn^i = Conn^i_jk g^jk
	local connTrace23Expr = (Gamma'^a_b^b')()
	self.connTrace23Codes = range(dim):mapi(function(i)
		local connTraceiCode = compile(connTrace23Expr[i])
if connTraceiCode ~= '0.' then
	dprint('connTrace23Code['..i..'] = '..substCoords(connTraceiCode))
end		
		return connTraceiCode
	end)

	-- sqrt(g)_,i / sqrt(g) = Conn^j_ij
	local connTrace13Expr = (Gamma'^b_ab')()
	self.connTrace13Codes = range(dim):mapi(function(i)
		local connTraceiCode = compile(connTrace13Expr[i])
if connTraceiCode ~= '0.' then
	dprint('connTrace13Code['..i..'] = '..substCoords(connTraceiCode))
end		
		return connTraceiCode
	end)


	-- dx is the change across the grid
	-- therefore it is based on the holonomic metric
	self.dxCodes = range(dim):mapi(function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'^b' * gHol'_ab')()
		local lenCode = compile((symmath.sqrt(lenSqExpr))())
if lenCode ~= '0.' then
	dprint('dxCode['..i..'] = '..substCoords(lenCode))
end	
		return lenCode
	end)

	self.g = g
	self.gCode = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local gijCode = compile(self.g[i][j])
if gijCode ~= '0.' then
	dprint('g['..i..']['..j..'] = '..substCoords(gijCode))
end			
			return gijCode, j
		end)
	end)

	local gU = Tensor('^ab', table.unpack((Matrix.inverse(g))))
	self.gU = gU
	self.gUCode = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local gUijCode = compile(self.gU[i][j])
if gUijCode ~= '0.' then
	dprint('gU['..i..']['..j..'] = '..substCoords(gUijCode))
end			
			return gUijCode, j
		end)
	end)
	
	local sqrt_gU = Tensor('^ab', function(a,b) return symmath.sqrt(gU[a][b])() end)
	self.sqrt_gUCode = range(dim):mapi(function(i)
		return range(i,dim):mapi(function(j)
			local sqrt_gUijCode = compile(sqrt_gU[i][j])
if sqrt_gUijCode ~= '0.' then
	dprint('sqrt(gU['..i..']['..j..']) = '..substCoords(sqrt_gUijCode))
end			
			return sqrt_gUijCode, j
		end)
	end)

	local volumeExpr = symmath.sqrt(symmath.Matrix.determinant(gHol))()
	self.volumeCode = compile(volumeExpr)
dprint('volumeCode = '..substCoords(self.volumeCode))
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

local function getCode_real3_to_3sym3(name, exprs)
	return template([[
_3sym3 <?=name?>(real3 pt) {
	return (_3sym3){
<? for i=1,3 do
?>	.<?=xs[i]?> = {
<?	for j=1,3 do
		for k=j,3 do
?>		.<?=xs[j]..xs[k]?> = <?=exprs[i] and exprs[i][j] and exprs[i][j][k]
			and convertParams(exprs[i][j][k]) or '0.'?>,
<?		end	
	end
?>	},
<?
end
?>	};
}]], {
		xs = xs,
		name = name,
		exprs = exprs,
		convertParams = convertParams,
	})
end


function CoordinateSystem:getCode(solver)
	self.solver = solver
	local dim = solver.dim
	local lines = table()
	
	-- dx0, ...
	-- this is the change in cartesian wrt the change in grid
	lines:append(range(dim):mapi(function(i)
		local code = self.dxCodes[i]
		for j=1,3 do
			code = code:gsub(
				'{pt^'..j..'}',
				'cell_x'..(j-1)..'(i.'..xs[j]..')')
		end
		return '#define dx'..(i-1)..'_at(i) (solver->grid_dx.s'..(i-1)..' * ('..code..'))'
	end))
	
	-- volume
	local volumeCode = '(' .. self.volumeCode .. ')'
	lines:insert('static inline '..getCode_real3_to_real('sqrt_det_g_grid', volumeCode))
	
	-- coord len code: l(v) = v^i v^j g_ij
	lines:append{
		'static inline '..getCode_real3_real3_to_real('coordLenSq', self.uLenSqCode),
		[[
static inline real coordLen(real3 r, real3 pt) {
	return sqrt(coordLenSq(r, pt));
}]],
	}

	lines:insert('static inline '..getCode_real3_real3_real3_to_real3('coord_conn_apply23', self.connApply23Codes))
	lines:insert('static inline '..getCode_real3_real3_real3_to_real3('coord_conn_apply12', self.connApply12Codes))
	lines:insert('static inline '..getCode_real3_to_real3('coord_conn_trace23', self.connTrace23Codes))
	lines:insert('static inline '..getCode_real3_to_real3('coord_conn_trace13', self.connTrace13Codes))
	lines:insert('static inline '..getCode_real3_to_3sym3('coord_conn', self.connCodes))

	--[[
	for i=0,dim-1 do
		lines:insert(getCode_real3_to_real('coordHolBasisLen'..i, self.eHolLenCode[i+1]))
	end
	--]]

	for i,eiCode in ipairs(self.eCode) do
		lines:insert('static inline '..getCode_real3_to_real3('coordBasis'..(i-1), eiCode))
	end

	lines:insert('static inline '..getCode_real3_real3_to_real3('coord_lower', self.lowerCodes))
	lines:insert('static inline '..getCode_real3_real3_to_real3('coord_raise', self.raiseCodes))

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
		
		addSym3Components('coord_g', self.gCode)
		addSym3Components('coord_gU', self.gUCode)
		addSym3Components('coord_sqrt_gU', self.sqrt_gUCode)
		lines:insert('static inline '..getCode_real3_to_sym3('coord_g', self.gCode))
		lines:insert('static inline '..getCode_real3_to_sym3('coord_gU', self.gUCode))
	end

	lines:insert(template([[

//converts a vector from cartesian coordinates to grid coordinates
//by projecting the vector into the grid basis vectors 
//at x, which is in grid coordinates
static inline real3 cartesianToCoord(real3 u, real3 pt) {
	real3 uCoord = real3_zero;
	<? for i=0,solver.dim-1 do ?>{
		real3 e = coordBasis<?=i?>(pt);
		//anholonomic normalized
		//uCoord.s<?=i?> = real3_dot(e, u); // / real3_len(e);
		//holonomic
		real uei = real3_dot(e, u) / real3_lenSq(e);
		uCoord.s<?=i?> = uei;
		//subtract off this basis component from u
		u = real3_sub(u, real3_real_mul(e, uei));
	}<? end ?>
	//add whatever's left of u
	uCoord = real3_add(uCoord, u);
	return uCoord;
}

//converts a vector from cartesian to grid
//by projecting it onto the basis ... ?
static inline real3 cartesianFromCoord(real3 u, real3 pt) {
	real3 uGrid = real3_zero;
	<? for i=0,solver.dim-1 do ?>{
		real3 e = coordBasis<?=i?>(pt);
		uGrid = real3_add(uGrid, real3_real_mul(e, u.s<?=i?>));
	}<? end ?>
	return uGrid;
}

]], {
		solver = solver,
	}))

	lines:insert('static inline '..self:getCoordMapCode())

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
