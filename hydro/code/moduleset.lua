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
print('building:')	
	local function add(name, from, indent, last)
		local range = require 'ext.range'
		indent = indent or ''
		
		local module, deps
		if not addedkeys[name] then 
			module = self.set[name]
			if not module then
				error("failed to find module "..name
					..(from and (" when requesting from module "..from)
						or " when requesting from root function call")
				)
			end
			addedkeys[name] = true	
			deps = module.depends:filter(function(dep)
				return not addedkeys[dep]
			end)
		end	
		local numdeps = deps and #deps or 0
		local indentlen = #indent
		local str = range(indentlen):mapi(function(i)
			if i < indentlen then 
				return '│' 
			end	
			return last and '┌' or '├'	-- └
		end):concat()
			..(numdeps > 0 and '┴' or '─')	-- ┬
			..'─> '..name
		-- don't include twice
		if module then
			assert(addedkeys[name])
			for i=1,numdeps do
				local dep = deps[i]
				add(dep, name, indent..' ', i == 1)--numdeps)
				--pm1 = pm1:sub(1, -2):gsub('.', '│')  .. (i < n and '├' or '└')
			end
			added:insert(module)
		end
		print(str)
	end
	local numModules = select('#', ...)
	for i=1,numModules do
		local name = select(i, ...)
		add(name, '', ' ', i == 1)--numModules)
	end
	print()
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
