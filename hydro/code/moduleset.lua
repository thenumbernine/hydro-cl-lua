local class = require 'ext.class'
local table = require 'ext.table'
local Module = require 'hydro.code.module'

--[[
in some ways I could reproduce this behavior with a bunch of #includes
but the .cl and .h code is templated with lua code
which means I'd have to put the templated .cl/.h code in one dir tree,
then generate it (per-solver at that) and put it in another dir tree,
and then use the #includes (templated per-solver #include locations) to handle dependencies.

so i could either make my dir structure incredibly complex
or i could just write this module system.
--]]

local ModuleSet = class()

function ModuleSet:init(...)
	self.set = {}
	-- copy any into it on ctor
	for i=1,select('#', ...) do
		local set = select(i, ...)
		assert(ModuleSet.is(set))
		for k,v in pairs(set.set) do
			self.set[k] = v
		end
	end
end

-- during init, building modules, do this:

function ModuleSet:add(args)
	local module = Module(args)
	if self.set[module.name] then
		error("trying to re-add module "..tostring(module.name))
	end
	self.set[module.name] = module
end

-- during init, generating code, do this:

function ModuleSet:getCodeForGetter(getter, ...)
	local addedkeys = {}
	local added = table()
	local function add(name)
		if addedkeys[name] then return end	-- don't include twice
		local module = self.set[name]
		if not module then
			error("failed to find module "..name)
		end
		addedkeys[name] = true	
		for _,dep in ipairs(module.depends) do
			add(dep)
		end
		added:insert(module)
	end
	for i=1,select('#', ...) do
		local name = select(i, ...)
		add(name)
	end
	return added:mapi(function(module)
		local code, desc = getter(module)
		if code == '' then return '' end
		return '\n////////////// '..module.name..' '..desc..' //////////////\n\n'
			..code
	end):concat'\n'
end

function ModuleSet:getTypeHeader(...)
	return table{
		self:getCodeForGetter(function(module)
			return module.typecode, 'typecode'
		end, ...),
		self:getCodeForGetter(function(module) 
			return module.structs:mapi(function(struct)
				-- this needs makeType() called first, which generates the .typecode
				-- but it also calls the ffi.metatype (which can only be done once)
				-- and also the ffi.cdef (which is only effective the first time it's done)
				return struct.typecode, 'typecode'
			end):concat'\n'
		end, ...),
	}:concat'\n'
end

function ModuleSet:getHeader(...)
	return table{
		self:getCodeForGetter(function(module) 
			return module.typecode, 'typecode' 
		end, ...),
		self:getCodeForGetter(function(module) 
			return module.structs:mapi(function(struct)
				-- this needs makeType() called first, which generates the .typecode
				-- but it also calls the ffi.metatype (which can only be done once)
				-- and also the ffi.cdef (which is only effective the first time it's done)
				return struct.typecode
			end):concat'\n', 'typecode'
		end, ...),	
		self:getCodeForGetter(function(module) 
			return module.headercode, 'headercode' 
		end, ...),
	}:concat'\n'
end

function ModuleSet:getCode(...)
	return self:getCodeForGetter(function(module) 
		return module.code, 'code' 
	end, ...)
end

return ModuleSet
