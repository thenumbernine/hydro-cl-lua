local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Torus = class(geometry)

Torus.name = 'torus' 
Torus.coords = {'u', 'v', 'r'}

function Torus:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	args.embedded = table{x,y,z}
	local R = 1	-- TODO make this a geom param
	local u, v, r = symmath.vars('u', 'v', 'r')
	args.coords = table{r, u, v}
	local a = (R - r) / 2
	local c = (R + r) / 2
	args.chart = function() 
		return Tensor('^I', 
			(c + a * cos(v)) * cos(u),
			(c + a * cos(v)) * sin(u),
			a * sin(v)) 
	end
	Torus.super.init(self, args)
end

return Torus
