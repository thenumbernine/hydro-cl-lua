local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local path = require 'ext.path'
local vector = require 'ffi.cpp.vector'
local Mesh = require 'hydro.mesh.mesh'
local MeshFactory = require 'hydro.mesh.factory'

local P2DFMTMeshFactory = class(MeshFactory)

P2DFMTMeshFactory.name = 'p2dfmt'
P2DFMTMeshFactory.dim = 2

function P2DFMTMeshFactory:init(args)
	self.meshfile = assert(args.meshfile)
end

function P2DFMTMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)
	
	local fn = 'grids/'..self.meshfile
	local ls = string.split(assert(string.trim(path(fn):read(), "failed to open "..fn)), '\n')
	local first = ls:remove(1)
	local m, n = string.split(string.trim(ls:remove(1)), '%s+'):mapi(function(l) return tonumber(l) end):unpack()
	local x = string.split(string.trim(ls:concat()), '%s+'):mapi(function(l) return tonumber(l) end)
	assert(#x == 2*m*n)
	-- [[
	local us = x:sub(1,m*n)
	local vs = x:sub(m*n+1)
	local numVtxs = #us
	assert(numVtxs == #vs)

	mesh.vtxs:resize(numVtxs)
	for i=1,numVtxs do
		mesh.vtxs.v[i-1] = mesh.real3(us[i],vs[i],0)
	end
	
	for j=0,n-2 do
		for i=0,m-2 do
			mesh:addCell(vector('int',{
				i + m * j,
				i+1 + m * j,
				i+1 + m * (j+1),
				i + m * (j+1),
			}))
		end
	end

	mesh:calcAux()
	return mesh
end

return P2DFMTMeshFactory 
