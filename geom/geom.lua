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
	
	local e1, e2, e3
	e1 = table{e[1][1] or const(0), e[1][2] or const(0), e[1][3] or const(0)}
	if dim >= 2 then
		e2 = table{e[2][1] or const(0), e[2][2] or const(0), e[2][3] or const(0)}
	else
		-- e2 is the maximal length orthogonal to e1
		local found = false
		for i =1,3 do
			local n = table{0,0,0}:map(function(x) return const(x) end)
			n[i] = 1
			e2 = cross(n,e1)
			local e2len = symmath.sqrt(e2:sum())()
			if e2len ~= const(0) then
				found = true
			end
		end
		if not found then
			error("how can you not find a perpendicular basis for e1="..e1:map(tostring):concat', ')
		end
	end
	if dim >= 3 then
		e3 = table{e[3][1] or const(0), e[3][2] or const(0), e[3][3] or const(0)}
	else
		e3 = cross(e1,e2)
	end

	local e1len = symmath.sqrt(e1[1]^2 + e1[2]^2 + e1[3]^2)()
	local e2len = symmath.sqrt(e2[1]^2 + e2[2]^2 + e2[3]^2)()
	local e3len = symmath.sqrt(e3[1]^2 + e3[2]^2 + e3[3]^2)()

	local e1unit = e1:map(function(e1i) return (e1i/e1len)() end)
	local e2unit = e2:map(function(e2i) return (e2i/e2len)() end)
	local e3unit = e3:map(function(e3i) return (e3i/e3len)() end)
	
	local es = table{e1,e2,e3}
	local eUnits = table{e1unit,e2unit,e3unit}

	self.eCode = es:map(function(ei,i)
		return ei:map(function(eij,j)
			return toC:compile(eij, toC_coordArgs):match'return (.*);'
		end)
	end)

	self.eUnitCode = eUnits:map(function(ei_unit,i)
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
