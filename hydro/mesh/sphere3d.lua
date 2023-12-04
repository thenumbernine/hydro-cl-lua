local table = require 'ext.table'
local Cube3DMeshFactory = require 'hydro.mesh.cube3d'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'

local Sphere3DMeshFactory = Cube3DMeshFactory:subclass()

Sphere3DMeshFactory.name = 'sphere3d'

function Sphere3DMeshFactory:init(args) 
	args = table(args)
-- TODO rmin=0 fails without any sort of caps	
	args.mins = vec3d(args.mins or {.1, 0, 0})
	args.maxs = vec3d(args.maxs or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {0, 0, 1})
--TODO cap for 3D not working yet, and we can't do thetaMin=0 without capMin
--	args.capmin = vec3i(args.capmin or {1, 1, 0})
--	args.capmax = vec3i(args.capmin or {0, 1, 0})
	Sphere3DMeshFactory.super.init(self, args)
end

function Sphere3DMeshFactory:coordChart(r, theta, phi)
	theta = theta * math.pi 
	phi = phi * 2 * math.pi
	local sinth = math.sin(theta)
	return 
		r * math.cos(phi) * sinth,
		r * math.sin(phi) * sinth,
		r * math.cos(theta)
end

return Sphere3DMeshFactory
