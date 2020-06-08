local class = require 'ext.class'
local table = require 'ext.table'
local Cube3DMeshFactory = require 'hydro.mesh.cube3d'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'

local Sphere3DMeshFactory = class(Cube3DMeshFactory)

Sphere3DMeshFactory.name = 'Sphere3DMesh'

function Sphere3DMeshFactory:init(args) 
	args = table(args)
	args.mins = vec3d(args.mins or {.5, .5, 0})
	args.maxs = vec3d(args.mins or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {0, 0, 1})
	args.capmin = vec3i(args.capmin or {1, 0, 0})
	Sphere3DMeshFactory.super.init(self, args)
end

function Sphere3DMeshFactory:coordChart(x) 
	local r = x.x
	local theta = x.y * math.pi 
	local phi = x.z * 2 * math.pi
	local sinth = sin(theta)
	return mesh.real3(
		r * math.cos(phi) * sinth,
		r * math.sin(phi) * sinth,
		r * math.cos(theta))
end

return Sphere3DMeshFactory
