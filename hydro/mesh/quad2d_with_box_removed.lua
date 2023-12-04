local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DWithBoxRemovedMeshFactory = Quad2DMeshFactory:subclass()

function Quad2DWithBoxRemovedMeshFactory:testMakeCell(i,j)
	local sx = tonumber(self.size.x)
	local sy = tonumber(self.size.y)
	local u = (i + .5) / sx
	local v = (j + .5) / sy
	return u <= .4 or u >= .6
		or v <= .4 or v >= .6
end

return Quad2DWithBoxRemovedMeshFactory 
