local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local Cartesian = class(CoordinateSystem)

Cartesian.name = 'cartesian' 
Cartesian.coords = {'x', 'y', 'z'}

function Cartesian:init(args)
	local x,y,z = symmath.vars('x', 'y', 'z')
	args.embedded = table{x,y,z}
	args.coords = args.embedded
	args.chart = function() return symmath.Tensor('^I', x,y,z) end
	Cartesian.super.init(self, args)
end

return Cartesian
