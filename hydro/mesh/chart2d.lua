local class = require 'ext.class'
local vec2sz = require 'vec-ffi.create_vec2'{ctype='size_t'}
local vec2i = require 'vec-ffi.vec2i'
local vec2d = require 'vec-ffi.vec2d'
local MeshFactory = require 'hydro.mesh.factory'

local Chart2DMeshFactory = class(MeshFactory)

function Chart2DMeshFactory:init(args)
	args = args or {}
	self.size = vec2sz(args.size or {20,20})
	self.mins = vec2d(args.mins or {-1, -1})
	self.maxs = vec2d(args.maxs or {1, 1})
	self.wrap = vec2i(args.wrap or {0, 0})
	self.capmin = vec2i(args.capmin or {0, 0})
end

function Chart2DMeshFactory:coordChart(x) return x end

return Chart2DMeshFactory 
