local class = require 'ext.class'
local math = require 'ext.math'
local vec3i = require 'vec-ffi.vec3i'
local vec3d = require 'vec-ffi.vec3d'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Polar2DMeshFactory = class(Quad2DMeshFactory)

function Polar2DMeshFactory:init(args)
	Polar2DMeshFactory.super.init(self, args)
	args.mins = vec3d(.1, 0, -1)
	args.maxs = vec3d(1, 1, 1)
	args.wrap = vec3i(0, 1, 0)
	args.capmin = vec3i(0, 0, 0)
end

-- hmm, same as cylinder3d ...
function Polar2DMeshFactory:coordChart(r,theta,z)
	theta = theta * 2 * math.pi
	return r * math.cos(theta), r * math.sin(theta), z
end

return Polar2DMeshFactory 
