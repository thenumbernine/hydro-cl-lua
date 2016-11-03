local class = require 'ext.class'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(geometry)

function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, theta, phi = symmath.vars('r', 'theta', 'phi')
	args.coords = table{r, theta, phi}
	args.chart = function() 
		return ({
			Tensor('^I', r),
			Tensor('^I', r * cos(theta), r * sin(theta)),
			Tensor('^I', r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)),
		})[args.solver.dim]
	end
	Sphere.super.init(self, args)
end

return Sphere
