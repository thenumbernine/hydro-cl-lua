local math = require 'ext.math'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DTwistMeshFactory = Quad2DMeshFactory:subclass()

Quad2DTwistMeshFactory.name = 'quad2dtwist' 

function Quad2DTwistMeshFactory:coordChart(x,y,z)
	local r = math.sqrt(x*x + y*y)
	--local theta = math.max(0, 1 - r)
	local sigma = 3		-- almost 0 at r=1
	local rotationAmplitude = 3
	local theta = rotationAmplitude * sigma * r * math.exp(-sigma * sigma * r * r)
	local costh = math.cos(theta)
	local sinth = math.sin(theta)
	return 
		costh * x - sinth * y,
		sinth * x + costh * y,
		0

end

return Quad2DTwistMeshFactory 
