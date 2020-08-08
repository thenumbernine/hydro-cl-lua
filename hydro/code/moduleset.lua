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

function ModuleSet:getCodeForField(field, name)
	local addedkeys = {}
	local added = table()
	local function add(name)
		if addedkeys[name] then return end	-- don't include twice
		local module = assert(self.set[name])
		addedkeys[name] = true	
		for _,dep in ipairs(module.depends) do
			add(dep)
		end
		added:insert(module)
	end
	add(name)
	
	return added:mapi(function(module)
		return '\n////////////// '..name..' '..field..' //////////////\n\n'
			..module[field]
	end):concat'\n'
end

function ModuleSet:getTypeHeader(name)
	return self:getCodeForField('typecode', name)
end

function ModuleSet:getHeader(name)
	return self:getCodeForField('typecode', name)..'\n'
		.. self:getCodeForField('headercode', name)
end

function ModuleSet:getCode(name)
	return self:getCodeForField('code', name)
end

return ModuleSet
