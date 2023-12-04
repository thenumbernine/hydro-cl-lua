local table = require 'ext.table'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'
local Cube3DMeshFactory = require 'hydro.mesh.cube3d'

local Cylinder3DMeshFactory = Cube3DMeshFactory:subclass()

Cylinder3DMeshFactory.name = 'cylinder3d'

function Cylinder3DMeshFactory:init(args) 
	args = table(args)
	args.mins = vec3d(args.mins or {0, 0, -.5})
	args.maxs = vec3d(args.maxs or {1, 1, .5})
	args.wrap = vec3i(args.wrap or {0, 1, 0})
-- TODO capmin.x==1 isn't working	
	args.capmin = vec3i(args.capmin or {0, 0, 0})
	Cylinder3DMeshFactory.super.init(self, args)
end

function Cylinder3DMeshFactory:coordChart(r, theta, z)
	theta = theta * 2 * math.pi
	return 
		r * math.cos(theta),
		r * math.sin(theta),
		z
end

return Cylinder3DMeshFactory 
