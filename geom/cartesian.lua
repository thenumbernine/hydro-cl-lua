local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local Geometry = require 'geom.geom'

local Cartesian = class(Geometry)

Cartesian.name = 'cartesian' 
Cartesian.coords = {'x', 'y', 'z'}

function Cartesian:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}:sub(1, args.solver.dim)
	args.coords = args.embedded
	args.chart = function() 
		return symmath.Tensor('^I', table.unpack(args.coords)) 
	end
	Cartesian.super.init(self, args)
end

return Cartesian
