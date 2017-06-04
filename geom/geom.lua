--[[

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

	-- commutation coefficients
	local c = Tensor'_ab^c'
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

	local g = (e'_u^I' * e'_v^J' * eta'_IJ')()
	
	-- TODO automatically do this ...
	g = g:map(function(expr)
		if symmath.op.pow.is(expr)
		and expr[2] == const(2)
		and symmath.cos.is(expr[1])
		then
			return 1 - symmath.sin(expr[1][1]:clone())^2
		end
	end)()
print'metric:'
print(var'g''_uv':eq(var'e''_u^I' * var'e''_v^J' * var'\\eta''_IJ'):eq(g'_uv'()))
	Tensor.metric(g)

print'commutation:'
print(var'c''_uv^w':eq(c'_uv^w'()))

	local GammaL = Tensor'_abc'
	GammaL['_abc'] = ((g'_ab,c' + g'_ac,b' - g'_bc,a' + c'_abc' + c'_acb' - c'_bca') / 2)()
print'1st kind Christoffel:'
print(var'\\Gamma''_abc':eq(symmath.op.div(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(GammaL'_abc'()))

	local Gamma = Tensor'^a_bc'
	Gamma['^a_bc'] = GammaL'^a_bc'()
print'connection:'
print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma'^a_bc'()))


	-- code generation

	
	local toC = require 'symmath.tostring.C'
	local toC_coordArgs = table.map(baseCoords, function(coord, i)
		return {[coord] = '{x'..i..'}'}	-- 1-based
	end)	
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


	self.uCode = range(dim):map(function(i) return compile(u[i]) end)

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

	local eExtLen = eExt:map(function(ei,i)
		return symmath.sqrt(ei:map(function(x) return x^2 end):sum())
	end)
	local eExtUnit = eExt:map(function(ei,i)
		return ei:map(function(eij) return (eij/eExtLen[i])() end)
	end)

	self.eCode = eExt:map(function(ei,i)
		return ei:map(compile)
	end)

	self.eUnitCode = eExtUnit:map(function(ei_unit,i)
		return ei_unit:map(compile)
	end)

	local coordU = Tensor('^a', function(a) return baseCoords[a] end)
	
	local lenSqExpr = (coordU'^a' * coordU'_a')()
	self.uLenSqCode = compile(lenSqExpr)
	self.uLenCode = compile((symmath.sqrt(lenSqExpr))())
	
	self.dxCodes = range(dim):map(function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'_a')()
		local lenCode = compile((symmath.sqrt(lenSqExpr))())
		return lenCode
	end)

	local volumeExpr = symmath.sqrt(symmath.Matrix.determinant(g))()
	self.volumeCode = compile(volumeExpr)
print('volume code',self.volumeCode)
end

return Geometry
