local class = require 'ext.class'
local table = require 'ext.table'
local Cube3DMeshFactory = require 'hydro.mesh.cube3d'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'

local Cylinder3DMeshFactory = class(Cube3DMeshFactory)

Cylinder3DMeshFactory.name = 'Cylinder3DMesh'

function Cylinder3DMeshFactory:init(args) 
	args = table(args)
	args.mins = vec3d(args.mins or {.5, .5, 0})
	args.maxs = vec3d(args.mins or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {0, 1, 0})
	args.capmin = vec3i(args.capmin or {1, 0, 0})
	Cylinder3DMeshFactory.super.init(self, args)
end

function Cylinder3DMeshFactory:coordChart(x) 
	local r = x.x
	local theta = x.y * 2 * math.pi
	local z = x.z
	return self.mesh.real3(
		r * math.cos(theta),
		r * math.sin(theta),
		z)
end

return Cylinder3DMeshFactory 
