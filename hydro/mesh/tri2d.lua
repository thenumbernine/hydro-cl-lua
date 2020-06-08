local class = require 'ext.class'
local vec2d = require 'vec-ffi.vec2d'
local vector = require 'hydro.util.vector'
local Mesh = require 'hydro.mesh.mesh'
local Chart2DMeshFactory = require 'hydro.mesh.chart2d'

local Tri2DMeshFactory = class(Chart2DMeshFactory)

Tri2DMeshFactory.name = 'Tri2DMesh'

function Tri2DMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)

	local nx = tonumber(self.size.x) + 1
	local ny = tonumber(self.size.y) + 1
	
	local stepx = 1
	local stepy = nx
	
	local vtxsize = nx * ny
	if self.capmin.x ~= 0 then vtxsize = vtxsize + 1 end

	mesh.vtxs:resize(vtxsize)

	local coordRangeMaxX = tonumber(self.size.x)
	local coordRangeMaxY = tonumber(self.size.y)
	if self.wrap.x ~= 0 or self.capmin.x ~= 0 then coordRangeMaxX = coordRangeMaxX + 1 end
	if self.wrap.y ~= 0 or self.capmin.y ~= 0 then coordRangeMaxY = coordRangeMaxY + 1 end
	
	local iofsx = self.capmin.x ~= 0 and 1 or 0
	local iofsy = self.capmin.y ~= 0 and 1 or 0

	for iy=0,ny-1 do
		for ix=0,nx-1 do
			local x = vec2d(
				(ix + iofsx) / coordRangeMaxX * (self.maxs.x - self.mins.x) + self.mins.x,
				(iy + iofsy) / coordRangeMaxY * (self.maxs.y - self.mins.y) + self.mins.y)
			local u = self:coordChart(x)
			mesh.vtxs.v[ix * stepx + iy * stepy] = mesh.real3(u:unpack())
		end
	end
	
	local capindex = nx * ny
	if self.capmin.x then
		local sum = mesh.real3()
		for j=0,ny-1 do
			sum = sum + mesh.vtxs.v[0 + nx * j]
		end
		mesh.vtxs.v[capindex] = sum / ny
	end
	
	local imaxx = self.wrap.x ~= 0 and nx or nx - 1
	local imaxy = self.wrap.y ~= 0 and ny or ny - 1
	
	for iy=0,imaxy-1 do
		local niy = (iy + 1) % ny
		for ix=0,imaxx-1 do
			local nix = (ix + 1) % nx
-- TODO just combine quad2d and tri2d, and only vary the 'addTri' function
-- [[ tris
			mesh:addCell(vector('int',{
				ix + nx * iy,
				nix + nx * iy,
				nix + nx * niy,
			}))
			mesh:addCell(vector('int',{
				nix + nx * niy,
				ix + nx * niy,
				ix + nx * iy,
			}))
--]]
--[[ quads
			mesh:addCell(vector('int',{
				ix + nx * iy,
				nix + nx * iy,
				nix + nx * niy,
				ix + nx * niy,
			}))
--]]
		end
	end

	if self.capmin.x ~= 0 then
		for j=0,imaxy-1 do
			local jn = (j + 1) % ny
			mesh.addCell(vector('int',{
				0 + nx * j,
				0 + nx * jn,
				capindex,
			}))
		end
	end

	mesh:calcAux()
	return mesh
end

return Tri2DMeshFactory
