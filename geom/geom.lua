--[[

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
	local dim = args.solver.dim
	local symmath = require 'symmath'
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

	print('coordinates:', table.unpack(coords))
	print('embedding:', table.unpack(embedded))
	
	local eta = Tensor('_IJ', table.unpack(flatMetric)) 
	print'flat metric:'
	print(var'\\eta''_IJ':eq(eta'_IJ'()))
	print()

	local u = args.chart()
	print'coordinate chart:'
	print(var'u''^I':eq(u'^I'()))
	print()
	
	local e = Tensor'_u^I'
	e['_u^I'] = u'^I_,u'()
	print'embedded:'
	print(var'e''_u^I':eq(var'u''^I_,u'):eq(e'_u^I'()))
	print()

	local g = (e'_u^I' * e'_v^J' * eta'_IJ')()
	-- TODO automatically do this ...
	g = g:map(function(expr)
		if symmath.powOp.is(expr)
		and expr[2] == symmath.Constant(2)
		and symmath.cos.is(expr[1])
		then
			return 1 - symmath.sin(expr[1][1]:clone())^2
		end
	end)()
	print'metric:'
	print(var'g''_uv':eq(var'e''_u^I' * var'e''_v^J' * var'\\eta''_IJ'):eq(g'_uv'()))
	Tensor.metric(g)

	local GammaL = Tensor'_abc'
	GammaL['_abc'] = ((g'_ab,c' + g'_ac,b' - g'_bc,a') / 2)()
	print'1st kind Christoffel:'
	print(var'\\Gamma''_abc':eq(symmath.divOp(1,2)*(var'g''_ab,c' + var'g''_ac,b' - var'g''_bc,a' + var'c''_abc' + var'c''_acb' - var'c''_bca')):eq(GammaL'_abc'()))

	local Gamma = Tensor'^a_bc'
	Gamma['^a_bc'] = GammaL'^a_bc'()
	print'connection:'
	print(var'\\Gamma''^a_bc':eq(var'g''^ad' * var'\\Gamma''_dbc'):eq(Gamma'^a_bc'()))

	local toC = require 'symmath.tostring.C'
	local toC_coordArgs = table.map(coords, function(coord, i)
		return {[coord] = '{x'..i..'}'}	-- 1-based
	end)	
	
	self.uCode = range(dim):map(function(i)
		return toC:compile(u[i], toC_coordArgs)
			:match'return (.*);'
	end)

	-- just giving up and manually writing this out

	local const = symmath.Constant
	
	local function cross(a,b)
		return table{
			(a[2]*b[3]-b[2]*a[3])(),
			(a[3]*b[1]-b[3]*a[1])(),
			(a[1]*b[2]-b[1]*a[2])(),
		}
	end

	-- extend 'e' to full R3 ... should I do this from the start?
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
		return ei:map(function(eij,j)
			return toC:compile(eij, toC_coordArgs):match'return (.*);'
		end)
	end)

	self.eUnitCode = eExtUnit:map(function(ei_unit,i)
		return ei_unit:map(function(eij_unit,j)
			return toC:compile(eij_unit, toC_coordArgs):match'return (.*);'
		end)
	end)

	local coordU = Tensor('^a', function(a) return coords[a] end)
	
	local lenSqExpr = (coordU'^a' * coordU'_a')()
	self.uLenSqCode = toC:compile(lenSqExpr, toC_coordArgs):match'return (.*);'
	self.uLenCode = toC:compile((symmath.sqrt(lenSqExpr))(), toC_coordArgs):match'return (.*);'
	
	self.dxCodes = range(dim):map(function(i)
		local dir = Tensor('^a', function(a) return a==i and 1 or 0 end)
		local lenSqExpr = (dir'^a' * dir'_a')()
		local lenSqCode = toC:compile(lenSqExpr, toC_coordArgs):match'return (.*);'
		local lenCode = toC:compile((symmath.sqrt(lenSqExpr))(), toC_coordArgs):match'return (.*);'
		return lenCode
	end)

end

return Geometry
