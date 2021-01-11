local class = require 'ext.class'
local Quad2DMeshFactory = require 'hydro.mesh.quad2d'

local Quad2DWithCylinderRemovedMeshFactory = class(Quad2DMeshFactory)

Quad2DWithCylinderRemovedMeshFactory.cylinderRadius = .1

function Quad2DWithCylinderRemovedMeshFactory:coordChart(x, y, z)
	-- morph it from a square to a circle to a square
	local ax = math.abs(x)
	local ay = math.abs(y)
	local az = math.abs(z)
	local linf = math.max(ax,ay,az)
	-- linf=1 => f=1, linf=r => f=0, linf=0 => f=0
	local r = self.cylinderRadius 
	local f = math.max(linf - r, 0) / (1 - r)
	local l2 = math.sqrt(x*x + y*y + z*z)
	--[[
	at linf=1 we want xyz unchanged <=> xyz' = xyz
	at linf=.5 we want |xyz|=.5 <=> xyz' = xyz*linf/l2
	so let's do xyz' = xyz * l
	and linf=1 => l=1=l2/l2, linf=.5 => l=linf/l2, linf=0 => l=0
	so f=1 => l=1, f=0 => l=|xyz|
	--]]
	local l = (f*l2 + (1-f)*linf)/l2
	x = x * l
	y = y * l
	z = z * l
	return x, y, z
end

function Quad2DWithCylinderRemovedMeshFactory:testMakeCell(i,j)
	local r = self.cylinderRadius 
	local sx = tonumber(self.size.x)
	local sy = tonumber(self.size.y)
	local u = (i + .5) / sx
	local v = (j + .5) / sy
	return u <= (.5-r) or u >= (.5+r)
		or v <= (.5-r) or v >= (.5+r)
end

return Quad2DWithCylinderRemovedMeshFactory 
