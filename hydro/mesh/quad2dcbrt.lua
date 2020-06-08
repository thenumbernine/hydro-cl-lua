local class = require 'ext.class'
local math = require 'ext.math'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DCbrtMeshFactory = class(Quad2DMeshFactory)

function Quad2DCbrtMeshFactory:coordChart(x)
	return self.mesh.real3(
		math.cbrt(x.x),
		math.cbrt(x.y),
		math.cbrt(x.z)
	)
end

return Quad2DCbrtMeshFactory 
