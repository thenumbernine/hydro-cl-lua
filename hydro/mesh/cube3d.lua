local class = require 'ext.class'
local vec3sz = require 'vec-ffi.vec3sz'
local vec3i = require 'vec-ffi.vec3i'
local vec3d = require 'vec-ffi.vec3d'
local vector = require 'hydro.util.vector'
local Mesh = require 'hydro.mesh.mesh'
local MeshFactory = require 'hydro.mesh.factory'

local Cube3DMeshFactory = class(MeshFactory)

Cube3DMeshFactory.name = 'cube3d'
Cube3DMeshFactory.dim = 3

function Cube3DMeshFactory:init(args)
	args = args or {}
	self.size = vec3sz(args.size or {10,10,10})
	self.mins = vec3d(args.mins or {-1, -1, -1})
	self.maxs = vec3d(args.maxs or {1, 1, 1})
	self.wrap = vec3i(args.wrap or {0, 0, 0})
	self.capmin = vec3i(args.capmin or {0,0,0})
end

function Cube3DMeshFactory:coordChart(x,y,z) 
	return x,y,z
end

function Cube3DMeshFactory:createMesh(solver)
	assert(solver.dim == 3)
	local mesh = Mesh(solver)
	self.mesh = mesh

	local nx = tonumber(self.size.x) + 1
	local ny = tonumber(self.size.y) + 1
	local nz = tonumber(self.size.z) + 1
	local stepx = 1
	local stepy = nx
	local stepz = nx * ny
	
	local vtxsize = nx * ny * nz
	mesh.vtxs:resize(vtxsize)

	local vtxmaxX = tonumber(self.size.x)
	local vtxmaxY = tonumber(self.size.y)
	local vtxmaxZ = tonumber(self.size.z)
	if self.wrap.x ~= 0 --[[or self.capmin.x ~= 0]] then vtxmaxX = vtxmaxX + 1 end
	if self.wrap.y ~= 0 --[[or self.capmin.y ~= 0]] then vtxmaxY = vtxmaxY + 1 end
	if self.wrap.z ~= 0 --[[or self.capmin.z ~= 0]] then vtxmaxZ = vtxmaxZ + 1 end

	local iofsx = --[[self.capmin.x ~= 0 and 1 or]] 0
	local iofsy = --[[self.capmin.y ~= 0 and 1 or]] 0
	local iofsz = --[[self.capmin.z ~= 0 and 1 or]] 0

	for iz=0,nz-1 do
		for iy=0,ny-1 do
			for ix=0,nx-1 do
				local x = (ix + iofsx) / vtxmaxX * (self.maxs.x - self.mins.x) + self.mins.x
				local y = (iy + iofsy) / vtxmaxY * (self.maxs.y - self.mins.y) + self.mins.y
				local z = (iz + iofsz) / vtxmaxZ * (self.maxs.z - self.mins.z) + self.mins.z
				mesh.vtxs.v[ix * stepx + iy * stepy + iz * stepz] = mesh.real3(self:coordChart(x,y,z))
			end
		end
	end
	
	local imaxx = self.wrap.x ~= 0 and nx or nx - 1
	local imaxy = self.wrap.y ~= 0 and ny or ny - 1
	local imaxz = self.wrap.z ~= 0 and nz or nz - 1
	
	for iz=0,imaxz-1 do
		local niz = (iz + 1) % nz
		for iy=0,imaxy-1 do
			local niy = (iy + 1) % ny
			for ix=0,imaxx-1 do
				local nix = (ix + 1) % nx
				mesh:addCell(vector('int',{
					--using z-order
					ix + nx * (iy + ny * iz),
					nix + nx * (iy + ny * iz),
					ix + nx * (niy + ny * iz),
					nix + nx * (niy + ny * iz),
					
					ix + nx * (iy + ny * niz),
					nix + nx * (iy + ny * niz),
					ix + nx * (niy + ny * niz),
					nix + nx * (niy + ny * niz),
				}))
			end
		end
	end

	mesh:calcAux()
	return mesh
end

return Cube3DMeshFactory
