local class = require 'ext.class'
local vec2i = require 'vec-ffi.vec2i'
local vec2d = require 'vec-ffi.vec2d'
local vector = require 'hydro.util.vector'
local Mesh = require 'hydro.mesh.mesh'
local Chart2DMeshFactory = require 'hydro.mesh.chart2d'

local Quad2DMeshFactory = class(Chart2DMeshFactory)

Quad2DMeshFactory.name = 'Quad2DMesh'

function Quad2DMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)

	local n = self.size + 1
	local step = vec2i(1, n.x)
	local vtxsize = n:volume()
	if self.capmin.x ~= 0 then vtxsize = vtxsize + 1 end
	mesh.vtxs:resize(vtxsize)

	local coordRangeMax = vec2i(self.size:unpack())
	if self.wrap.x ~= 0 or self.capmin.x ~= 0 then coordRangeMax.x = coordRangeMax.x + 1 end
	if self.wrap.y ~= 0 or self.capmin.y ~= 0 then coordRangeMax.y = coordRangeMax.y + 1 end

	local iofs = vec2i()
	if self.capmin.x ~= 0 then iofs.x = 1 end
	if self.capmin.y ~= 0 then iofs.y = 1 end

	local i = vec2i()
	for iy=0,tonumber(n.y)-1 do
		i.y = iy
		for ix=0,tonumber(n.x)-1 do
			i.x = ix
			local x = vec2d((i + iofs):unpack()) / vec2d(coordRangeMax:unpack()) * (self.maxs - self.mins) + self.mins
			local u = self:coordChart(x)
			mesh.vtxs.v[tonumber(i:dot(step))] = mesh.real3(u:unpack())
		end
	end
	
	local capindex = n:volume()
	if self.capmin.x ~= 0 then
		local sum = mesh.real3()
		for j=0,n.y-1 do
			sum = sum + mesh.vtxs.v[tonumber(0 + n.x * j)]
		end
		mesh.vtxs.v[tonumber(capindex)] = sum / tonumber(n.y)
	end

	local imax = vec2i()
	for j=0,1 do
		imax.s[j] = self.wrap.s[j] ~= 0 and n.s[j] or n.s[j]-1
	end

	local ni = vec2i()
	for iy=0,tonumber(imax.y-1) do
		i.y = iy
		ni.y = (i.y + 1) % n.y
		for ix=0,tonumber(imax.x-1) do
			i.x = ix
			ni.x = (i.x + 1) % n.x
			mesh:addCell(vector('int',{
				tonumber(i.x + n.x * i.y),
				tonumber(ni.x + n.x * i.y),
				tonumber(ni.x + n.x * ni.y),
				tonumber(i.x + n.x * ni.y),
			}))
		end
	end

	if self.capmin.x ~= 0 then
		for j=0,imax.y-1 do
			local jn = (j + 1) % n.y
			mesh.addCell(vector('int',{
				tonumber(0 + n.x * j), 
				tonumber(0 + n.x * jn), 
				tonumber(capindex),
			}))
		end
	end

	mesh:calcAux()
	return mesh

end

return Quad2DMeshFactory
