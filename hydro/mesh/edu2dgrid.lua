local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local path = require 'ext.path'
local vector = require 'ffi.cpp.vector'
local Mesh = require 'hydro.mesh.mesh'
local MeshFactory = require 'hydro.mesh.factory'

local Edu2DGridMeshFactory = class(MeshFactory)

-- .grid files from ossanworld.com / "I Do Like CFD"
Edu2DGridMeshFactory.name = 'edu2dgrid'
Edu2DGridMeshFactory.dim = 2

function Edu2DGridMeshFactory:init(args)
	self.meshfile = assert(args.meshfile)
end

local function lineToNums(l, n)
	assert(l, "expected newline")
	local t = {string.split(string.trim(l), '%s+'):mapi(function(l) return tonumber(l) end):unpack()}
	assert(#t == n, "expected "..n.." cols")
	return table.unpack(t)
end

function Edu2DGridMeshFactory:createMesh(solver)
	local mesh = Mesh(solver)

	-- 1) read the .grid file

	local fn = 'grids/'..self.meshfile
	local ls = string.split(string.trim(assert(path(fn):read(), "failed to open "..fn)), '\n')
	local n = #ls
	local numVtxs, numTris, numQuads = lineToNums(ls:remove(1), 3)

	mesh.vtxs:resize(numVtxs)
	for i=1,numVtxs do
		-- is the vtx always 2D? can it be 3D?
		local x, y = lineToNums(ls:remove(1), 2)
		mesh.vtxs.v[i-1] = mesh.real3(x,y,0)
	end
	for i=1,numTris do
		local a,b,c = lineToNums(ls:remove(1), 3)
		mesh:addCell(vector('int',{a-1,b-1,c-1}))
	end
	for i=1,numQuads do
		local a,b,c,d = lineToNums(ls:remove(1), 4)
		mesh:addCell(vector('int',{a-1,b-1,c-1,d-1}))
	end

	local numBoundaryMethods = lineToNums(ls:remove(1), 1)

	local numBoundVtxs = {}
	for i=1,numBoundaryMethods do
		numBoundVtxs[i] = lineToNums(ls:remove(1), 1)
	end

	-- hmm, ossanworld's edu2d .grid format stores boundary info as lists of vtxs
	for boundaryMethodIndex=1,numBoundaryMethods do
		local prevvi
		for j=1,numBoundVtxs[boundaryMethodIndex] do
			local vi = lineToNums(ls:remove(1), 1)
			if not(1 <= vi and vi <= #mesh.vtxs) then
				error("got an oob index for the vertex of our boundary class "..tostring(vi).." on line "..(n-#ls))
			end
			vi = vi - 1	-- make it 0-based
			if j > 1 then
				local fi = mesh:findFaceForVtxs(prevvi, vi)
				assert(fi, "failed to find a face between vertexes "..prevvi.." and "..vi)
				-- find the face from prevvi to vi
				mesh.faces.v[fi].boundaryMethodIndex = boundaryMethodIndex
			end
			prevvi = vi
		end
	end

	assert(#ls == 0, "found remaining lines")

	mesh:calcAux()
	
	return mesh
end

return Edu2DGridMeshFactory 
