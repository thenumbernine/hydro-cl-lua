local class = require 'ext.class'
local symmath = require 'symmath'
local Geometry = require 'geom.geom'

local Cylinder = class(Geometry)

function Cylinder:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, theta, z = symmath.vars('r', 'theta', 'z')
	args.coords = table{r, theta, z}
	args.chart = function() 
		return ({
			Tensor('^I', r),
			Tensor('^I', r * symmath.cos(theta), r * symmath.sin(theta)),
			Tensor('^I', r * symmath.cos(theta), r * symmath.sin(theta), z),
		})[args.solver.dim]
	end
	Cylinder.super.init(self, args)
end

return Cylinder
