local class = require 'ext.class'
local table = require 'ext.table'
local vec3d = require 'vec-ffi.vec3d'
local vec3i = require 'vec-ffi.vec3i'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Cylinder2DMeshFactory = class(Quad2DMeshFactory)

Cylinder2DMeshFactory.name = 'cylinder2d' 
Cylinder2DMeshFactory.R = 1

function Cylinder2DMeshFactory:init(args)
	args = table(args)
	args.mins = vec3d(args.mins or {0, 0, -1})
	args.maxs = vec3d(args.maxs or {1, 1, 1})
	args.wrap = vec3i(args.wrap or {1, 0, 0})
	args.capmin = vec3i(args.capmin or {0, 0, 0})
	Cylinder2DMeshFactory.super.init(self, args)
end

-- hmm, same as cylinder3d ...
function Cylinder2DMeshFactory:coordChart(phi,z,r)
	r = r + self.R
	phi = phi * 2 * math.pi
	return r * math.cos(phi), r * math.sin(phi), z
end

return Cylinder2DMeshFactory 
