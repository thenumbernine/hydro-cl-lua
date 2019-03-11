--[[
This is just spherical
I also used inclination instead of declination (lambda instead of theta)
so I could keep coord1 at 0 and still be fine.
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(CoordinateSystem)
Sphere.name = 'sphere1d'
Sphere.coords = table{'r', 'λ', 'φ'}
function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, lambda, phi = symmath.vars('r', 'λ', 'φ')
	args.coords = table{r, lambda, phi}
	args.chart = function() 
		return Tensor('^I', 
			r * cos(lambda) * cos(phi), 
			r * cos(lambda) * sin(phi), 
			r * sin(lambda)
		)
	end
	Sphere.super.init(self, args)
end

return Sphere
