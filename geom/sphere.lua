local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local Geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(Geometry)

Sphere.name = 'sphere' 
Sphere.coords = {'r', 'θ', 'φ'}

function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}:sub(1,args.solver.dim)
	local r, theta, phi = symmath.vars('r', 'θ', 'φ')

	-- [[ holonomic
	args.coords = ({
		function() return table{theta} end,
		function() return table{theta, phi} end,
		function() return table{r, theta, phi} end,
	})[args.solver.dim]()
	--]]
	--[[ anholonomic
	local thetaHat = symmath.var'thetaHat'
	thetaHat.base = theta
	function thetaHat:applyDiff(expr) return expr:diff(theta) / r end

	local phiHat = symmath.var'phiHat'
	phiHat.base = phi
	function phiHat:applyDiff(expr) return expr:diff(phi) / (r * sin(theta)) end
	
	args.coords = table{thetaHat, phiHat, r}:sub(1, args.solver.dim)
	--]]
	
	args.chart = ({
		function() return Tensor('^I', theta) end,
		function() return Tensor('^I', sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)) end,
		function() 
			return Tensor('^I', 
				r * sin(theta) * cos(phi), 
				r * sin(theta) * sin(phi), 
				r * cos(theta)
			) 
		end,
	})[args.solver.dim]
	
	Sphere.super.init(self, args)
end

return Sphere
