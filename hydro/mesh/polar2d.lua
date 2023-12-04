local table = require 'ext.table'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Polar2DMeshFactory = Quad2DMeshFactory:subclass()

Polar2DMeshFactory.name = 'polar2d'

function Polar2DMeshFactory:init(args)
	args = table(args)
	args.mins = vec3d(args.mins or {0, 0, -1})
	args.maxs = vec3d(args.maxs or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {0, 1, 0})
	args.capmin = vec3i(args.capmin or {0, 0, 0})
	Polar2DMeshFactory.super.init(self, args)
end

function Polar2DMeshFactory:coordChart(r,phi,z)
	phi = phi * 2 * math.pi
	return r * math.cos(phi), r * math.sin(phi), z
end

return Polar2DMeshFactory 
