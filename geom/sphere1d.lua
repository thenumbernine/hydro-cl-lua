local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(geometry)
Sphere.name = 'sphere1d'
function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}:sub(1,args.solver.dim)
	local r, theta, phi = symmath.vars('r', 'theta', 'phi')
	
	args.coords = table{r, theta, phi}:sub(1, args.solver.dim)
	args.chart = function() 
		return assert(({
			function() return Tensor('^I', r * r) end,
		})[args.solver.dim]())
	end
	Sphere.super.init(self, args)
end

return Sphere
