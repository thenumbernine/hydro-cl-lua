local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Cylinder = class(geometry)

function Cylinder:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}:sub(1, args.solver.dim)
	local r, theta, z = symmath.vars('r', 'theta', 'z')
	args.coords = table{r, theta, z}:sub(1, args.solver.dim)
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
