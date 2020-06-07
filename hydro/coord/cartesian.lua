local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'hydro.coord.coord'

local Cartesian = class(CoordinateSystem)
Cartesian.name = 'cartesian' 

function Cartesian:init(args)
	local x,y,z = symmath.vars('x', 'y', 'z')
	self.embedded = table{x,y,z}
	self.baseCoords = self.embedded
	self.chart = function() return symmath.Tensor('^I', x,y,z) end
	Cartesian.super.init(self, args)

	self.vars = {
		x = x,
		y = y,
		z = z,
		r = (x^2 + y^2 + z^2)^.5,
	}
end

return Cartesian
