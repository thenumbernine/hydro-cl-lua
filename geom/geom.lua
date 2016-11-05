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
