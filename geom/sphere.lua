local class = require 'ext.class'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(geometry)
Sphere.name = 'sphere' 
function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}:sub(1,args.solver.dim)
	local r, theta, phi = symmath.vars('r', 'theta', 'phi')
	
	local thetaHat = symmath.var'thetaHat'
	thetaHat.base = theta
	function thetaHat:applyDiff(expr) return expr:diff(theta) / r end

	local phiHat = symmath.var'phiHat'
	phiHat.base = phi
	function phiHat:applyDiff(expr) return expr:diff(phi) / (r * sin(theta)) end
	
	args.coords = table{thetaHat, phiHat, r}:sub(1, args.solver.dim)
	args.chart = function() 
		return ({
			function() return Tensor('^I', theta) end,
			function() return Tensor('^I', sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)) end,
			function() return Tensor('^I', r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)) end,
		})[args.solver.dim]()
	end
	Sphere.super.init(self, args)
end

return Sphere
