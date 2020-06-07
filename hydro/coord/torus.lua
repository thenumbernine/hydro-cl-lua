local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'hydro.coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Torus = class(CoordinateSystem)

Torus.name = 'torus' 

function Torus:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	self.embedded = table{x,y,z}
	local R = 1	-- TODO make this a param
	local u, v, r = symmath.vars('u', 'v', 'r')
	self.baseCoords = table{r, u, v}
	local a = (R - r) / 2
	local c = (R + r) / 2
	self.chart = function() 
		return Tensor('^I', 
			(c + a * cos(v)) * cos(u),
			(c + a * cos(v)) * sin(u),
			a * sin(v)) 
	end
	Torus.super.init(self, args)
end

return Torus
