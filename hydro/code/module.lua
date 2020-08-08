--[[
since compile times on my intel card are so slow
and because I am fearing that a mystery error that I seem to get on both my Intel HD 520 and on the NVidia 1080 Ti's could be due to oversized binaries
I'm going to make a system for only generating what code is required

modules will have:
- name
- other required module's names
- type code
- Struct objects (which will be generated & added to type & cl code)
- header code / #defines
- cl code

so the solver base will have to first build allll possible modules, before typecode-gen or code-gen
and then during codegen it will have to only pick out the modules that it needs
--]]

local class = require 'ext.class'
local table = require 'ext.table'

local Module = class()

function Module:init(args)
	self.name = assert(args.name)
	self.depends = table(args.depends)
	-- structs & typedefs
	self.typecode = args.typecode or ''
	-- structs
	self.structs = table(args.structs)
	-- #defines 
	self.headercode = args.headercode or ''
	-- functions
	self.code = args.code or ''
end

return Module
