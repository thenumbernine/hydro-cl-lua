local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Cylinder = class(geometry)
Cylinder.name = 'cylinder' 
function Cylinder:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	args.embedded = table{x,y,z}:sub(1, args.solver.dim)
	
	local r, theta = symmath.vars('r', 'theta')

	--[[ holonomic (explodes, probably because it still need to include g^ij into P and v_i)
	args.coords = table{r, theta, z}:sub(1, args.solver.dim)
	--]]
	-- [[ anholonomic
	local thetaHat = symmath.var'thetaHat'
	thetaHat.base = theta
	function thetaHat:applyDiff(x) return x:diff(theta) / r end
	args.coords = table{r, thetaHat, z}:sub(1, args.solver.dim)
	--]]

	args.chart = function() 
		return ({
			function() return Tensor('^I', r) end,
			function() return Tensor('^I', r * cos(theta), r * sin(theta)) end,
			function() return Tensor('^I', r * cos(theta), r * sin(theta), z) end,
		})[args.solver.dim]()
	end
	
	Cylinder.super.init(self, args)
end

return Cylinder
