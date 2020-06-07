local class = require 'ext.class'
local string = require 'ext.string'
local file = require 'ext.file'
local MeshFactory = require 'hydro.mesh.factory'

local P2DFMTMeshFactory = class(MeshFactory)

function P2DFMTMeshFactory:createMesh(args)
	-- TODO
	local meshfilename = assert(args.meshfile)
	local ls = string.split(assert(string.trim(file['grids/'..meshfilename..'.p2dfmt'])), '\n')
	local first = ls:remove(1)
	local m, n = string.split(string.trim(ls:remove(1)), '%s+'):map(function(l) return tonumber(l) end):unpack()
	local x = string.split(string.trim(ls:concat()), '%s+'):map(function(l) return tonumber(l) end)
	assert(#x == 2*m*n)
	print(m, n, m*n)
	-- [[
	local us = x:sub(1,m*n)
	local vs = x:sub(m*n+1)
	assert(#us == #vs)
	local vtxs = table()
	for i=1,#us do
		local u,v = us[i], vs[i]
		addVtx(vtxs, {u,v})
	end

	assert(#vtxs == m*n, "expected "..#vtxs.." to equal "..(m*n))
end

return P2DFMTMeshFactory 
