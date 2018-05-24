local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(geometry)
Sphere.name = 'sphere1d'
Sphere.coords = table{'r', 'θ', 'φ'}
function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, theta, phi = symmath.vars('r', 'θ', 'φ')
	args.coords = table{r, theta, phi}
	args.chart = function() return Tensor('^I', r * r, 0, 0) end
	Sphere.super.init(self, args)
end

return Sphere
