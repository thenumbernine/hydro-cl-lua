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
	local theta, phi, r = symmath.vars('θ', 'φ', 'r')

	-- [[ holonomic
	args.coords = table{theta, phi, r}	--:sub(1, args.solver.dim)
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
	
	args.chart = function() 
		return Tensor('^I', 
			r * sin(theta) * cos(phi), 
			r * sin(theta) * sin(phi), 
			r * cos(theta)
		) 
	end
	
	Sphere.super.init(self, args)
end

return Sphere
