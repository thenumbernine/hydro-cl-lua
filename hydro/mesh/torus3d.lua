local class = require 'ext.class'
local table = require 'ext.table'
local Cube3DMeshFactory = require 'hydro.mesh.cube3d'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'

local Torus3DMeshFactory = class(Cube3DMeshFactory)

Torus3DMeshFactory.name = 'Torus3DMesh'

Torus3DMeshFactory.R = 2

function Torus3DMeshFactory:init(args) 
	args = table(args)
	args.mins = vec3d(args.mins or {0, 0, 0})
	args.maxs = vec3d(args.mins or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {0, 1, 1})
	args.capmin = vec3i(args.capmin or {0, 0, 0})
	Torus3DMeshFactory.super.init(self, args)
end

function Torus3DMeshFactory:coordChart(x) 
	local r = x.x
	local theta = x.y * 2 * math.pi
	local phi = x.z * 2 * math.pi
	return mesh.real3(
		(r * math.cos(theta) + self.R) * math.cos(phi), 
		(r * math.cos(theta) + self.R) * math.sin(phi), 
		-r * math.sin(theta))
end

return Torus3DMeshFactory 