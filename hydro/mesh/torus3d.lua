local class = require 'ext.class'
local table = require 'ext.table'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'
local Cube3DMeshFactory = require 'hydro.mesh.cube3d'

local Torus3DMeshFactory = class(Cube3DMeshFactory)

Torus3DMeshFactory.name = 'torus3d'

Torus3DMeshFactory.R = 2

function Torus3DMeshFactory:init(args) 
	args = table(args)
	args.mins = vec3d(args.mins or {0, 0, 0})
	args.maxs = vec3d(args.mins or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {0, 1, 1})
	args.capmin = vec3i(args.capmin or {0, 0, 0})
	Torus3DMeshFactory.super.init(self, args)
end

function Torus3DMeshFactory:coordChart(r, theta, phi)
	theta = theta * 2 * math.pi
	phi = phi * 2 * math.pi
	return 
		(r * math.cos(theta) + self.R) * math.cos(phi), 
		(r * math.cos(theta) + self.R) * math.sin(phi), 
		-r * math.sin(theta)
end

return Torus3DMeshFactory 
