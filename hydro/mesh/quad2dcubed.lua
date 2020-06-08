local class = require 'ext.class'
local math = require 'ext.math'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DCubedMeshFactory = class(Quad2DMeshFactory)

local function cubed(x)
	return x * x * x
end

function Quad2DCubedMeshFactory:coordChart(x)
	return self.mesh.real3(
		cubed(x.x),
		cubed(x.y),
		cubed(x.z)
	)
end

return Quad2DCubedMeshFactory 
