local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DCubedMeshFactory = Quad2DMeshFactory:subclass()

Quad2DCubedMeshFactory.name = 'quad2dcubed'

local function cubed(x)
	return x * x * x
end

function Quad2DCubedMeshFactory:coordChart(x,y,z)
	return 
		cubed(x.x),
		cubed(x.y),
		cubed(x.z)
end

return Quad2DCubedMeshFactory 
