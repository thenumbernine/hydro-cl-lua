local class = require 'ext.class'
local math = require 'ext.math'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DCbrtMeshFactory = class(Quad2DMeshFactory)

Quad2DCbrtMeshFactory.name = 'quad2dcbrt' 

function Quad2DCbrtMeshFactory:coordChart(x,y,z)
	return 
		math.cbrt(x),
		math.cbrt(y),
		math.cbrt(z)
end

return Quad2DCbrtMeshFactory 
