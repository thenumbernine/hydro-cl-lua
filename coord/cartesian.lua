local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local Cartesian = class(CoordinateSystem)
Cartesian.name = 'cartesian' 

function Cartesian:init(args)
	local x,y,z = symmath.vars('x', 'y', 'z')
	self.embedded = table{x,y,z}
	self.baseCoords = self.embedded
	self.chart = function() return symmath.Tensor('^I', x,y,z) end
	Cartesian.super.init(self, args)
end

return Cartesian
