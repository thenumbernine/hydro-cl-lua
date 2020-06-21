local class = require 'ext.class'
local vec3sz = require 'vec-ffi.vec3sz'
local vec3i = require 'vec-ffi.vec3i'
local vec3d = require 'vec-ffi.vec3d'
local vector = require 'hydro.util.vector'
local Mesh = require 'hydro.mesh.mesh'
local MeshSolver = require 'hydro.mesh.factory'

local Quad2DMeshFactory = class(MeshFactory)

Quad2DMeshFactory.name = 'quad2d'
Quad2DMeshFactory.dim = 2

function Quad2DMeshFactory:init(args)
	args = args or {}
	self.size = vec3sz(args.size or {20, 20, 1})
	self.mins = vec3d(args.mins or {-1, -1, -1})
	self.maxs = vec3d(args.maxs or {1, 1, 1})
	self.wrap = vec3i(args.wrap or {0, 0, 0})
	self.capmin = vec3i(args.capmin or {0, 0, 0})
	self.triangulate = args.triangulate
end

function Quad2DMeshFactory:coordChart(x, y, z)
	return x, y, z
end

function Quad2DMeshFactory:testMakeCell(i, j) 
	return true
end

function Quad2DMeshFactory:addPoly(mesh, ...)
	if not self.triangulate then
		mesh:addCell(vector('int',{...}))
	else
		local va = select(1, ...)
		local vb = select(2, ...)
		for i=3,select('#', ...) do
			local vc = select(i, ...)
			mesh:addCell(vector('int', {va,vb,vc}))
			vb = vc
		end
	end
end

function Quad2DMeshFactory:createMesh(solver)
assert(solver.dim == 2)	-- or maybe manually set it to 2?
	local mesh = Mesh(solver)
	self.mesh = mesh

	local nx = tonumber(self.size.x) + 1
	local ny = tonumber(self.size.y) + 1
	if self.wrap.x ~= 0 or self.capmin.x ~= 0 then nx = nx - 1 end
	if self.wrap.y ~= 0 or self.capmin.y ~= 0 then ny = ny - 1 end
	
	local stepx = 1
	local stepy = nx
	
	local vtxsize = nx * ny
	if self.capmin.x ~= 0 then vtxsize = vtxsize + 1 end
	
	mesh.vtxs:resize(vtxsize)

	local vtxmaxX = tonumber(self.size.x)
	local vtxmaxY = tonumber(self.size.y)
	if self.capmin.x ~= 0 then vtxmaxX = vtxmaxX + 1 end
	if self.capmin.y ~= 0 then vtxmaxY = vtxmaxY + 1 end

	local iofsx = self.capmin.x ~= 0 and 1 or 0
	local iofsy = self.capmin.y ~= 0 and 1 or 0

	for iy=0,ny-1 do
		for ix=0,nx-1 do
			local x = (ix + iofsx) / vtxmaxX * (self.maxs.x - self.mins.x) + self.mins.x
			local y = (iy + iofsy) / vtxmaxY * (self.maxs.y - self.mins.y) + self.mins.y
			local z = 0
			mesh.vtxs.v[ix * stepx + iy * stepy] = mesh.real3(self:coordChart(x,y,z))
		end
	end
	
	local capindex = nx * ny
	if self.capmin.x ~= 0 then
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
			if self:testMakeCell(ix,iy) then
				self:addPoly(mesh, 
					ix + nx * iy,
					nix + nx * iy,
					nix + nx * niy,
					ix + nx * niy)
			end
		end
	end

	if self.capmin.x ~= 0 then
		for j=0,imaxy-1 do
			local jn = (j + 1) % ny
			self:addPoly(mesh, 
				0 + nx * j,
				0 + nx * jn,
				capindex)
		end
	end

	mesh:calcAux()
	return mesh
end

return Quad2DMeshFactory
