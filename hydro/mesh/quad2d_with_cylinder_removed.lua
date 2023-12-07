local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DWithCylinderRemovedMeshFactory = Quad2DMeshFactory:subclass()

Quad2DWithCylinderRemovedMeshFactory.cylinderRadius = .2

function Quad2DWithCylinderRemovedMeshFactory:coordChart(x, y, z)
	-- morph it from a square to a circle to a square
	local ax = math.abs(x)
	local ay = math.abs(y)
	local az = math.abs(z)
	local linf = math.max(ax,ay,az)
	-- linf=1 => f=1, linf=r => f=0, linf=0 => f=0
	local r = self.cylinderRadius
	local f = math.max(linf - .5, 0) / .5
	local l2 = math.sqrt(x*x + y*y + z*z)
	if l2 == 0 then return 0,0,0 end
	--[[
	at linf=1 we want xyz unchanged <=> xyz' = xyz
	at linf=.5 we want |xyz|=.5 <=> xyz' = xyz*linf/l2
	so let's do xyz' = xyz * l
	and linf=1 => l=1=l2/l2, linf=.5 => l=linf/l2, linf=0 => l=0
	so f=1 => l=1, f=0 => l=|xyz|
	--]]
	local l = f + (1 - f) * linf/l2 * 2 * self.cylinderRadius
	x = x * l
	y = y * l
	z = z * l
	return x, y, z
end

function Quad2DWithCylinderRemovedMeshFactory:testMakeCell(i,j)
	local sx = tonumber(self.size.x)
	local sy = tonumber(self.size.y)
	local u = (i + .5) / sx
	local v = (j + .5) / sy
	return u <= .25 or u >= .75
		or v <= .25 or v >= .75
end

return Quad2DWithCylinderRemovedMeshFactory 
