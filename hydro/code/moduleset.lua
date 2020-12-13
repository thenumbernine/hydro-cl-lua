local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
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
	self.verbose = cmdline.moduleVerbose
end

-- during init, building modules, do this:

function ModuleSet:add(args)
	local module = Module(args)
	module.traceback = debug.traceback()
	if self.set[module.name] then
		error("trying to re-add module "..tostring(module.name).."\n"
			.."originally added from "..self.set[module.name].traceback)
	end
	self.set[module.name] = module
end

-- during init, generating code, do this:

function ModuleSet:getDependentModules(...)
	local addedkeys = {}
	local added = table()
if self.verbose then
	print('building:')
end	
	local function add(name, from, indent, last)
		local range = require 'ext.range'
		indent = indent or ''
		local indentlen = #indent
		
		local module, deps
		if not addedkeys[name] then
			module = self.set[name]
			if not module then
				error(
					(from and ("module %q"):format(from) or "function call").." requested to include module "..("%q"):format(name)
					.." but I failed to find it.\n"
					.."modules defined so far: "..table.keys(self.set):concat' '
				)
			end
			addedkeys[name] = true
			deps = module.depends:filter(function(dep)
				return not addedkeys[dep]
			end)
		end
		local numdeps = deps and #deps or 0
		
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
				add(dep, name, indent..' ', i == 1)
			end
			added:insert(module)
		end
if self.verbose then	
	print(str)
end	
	end
	local numModules = select('#', ...)
	for i=1,numModules do
		local name = select(i, ...)
		add(name, nil, ' ', i == 1)
	end
if self.verbose	then
	print()
end
	return added
end

local function comment(code, name, desc)
	code = string.trim(code)
	if code == '' then return end
	return '\n////////////// '..name..' '..desc..' //////////////\n\n'..code..'\n'
end

-- typecode + structs
function ModuleSet:getTypeHeader(...)
	local deps = self:getDependentModules(...)
	local lines = table()
	-- typecode & structs
	for _,module in ipairs(deps) do
		lines:insert(comment(table{
				module.typecode
			}:append(module.structs:mapi(function(struct)
				-- this needs makeType() called first, which generates the .typecode
				-- but it also calls the ffi.metatype (which can only be done once)
				-- and also the ffi.cdef (which is only effective the first time it's done)
				return struct.typecode
			end)):concat'\n', module.name, 'typecode & structs'
		) or nil)
	end
	return lines:concat'\n'
end

-- typecode + structs + header
function ModuleSet:getHeader(...)
	local deps = self:getDependentModules(...)
	local lines = table()
	-- typecode & structs
	for _,module in ipairs(deps) do
		lines:insert(comment(table{
				module.typecode
			}:append(module.structs:mapi(function(struct)
				-- this needs makeType() called first, which generates the .typecode
				-- but it also calls the ffi.metatype (which can only be done once)
				-- and also the ffi.cdef (which is only effective the first time it's done)
				return struct.typecode
			end)):concat'\n', module.name, 'typecode & structs'
		) or nil)
	end
	-- headercode
	for _,module in ipairs(deps) do
		-- allow lazy code generation
		if type(module.headercode) == 'function' then
			module.headercode = module.headercode()
		end
		
		lines:insert(comment(module.headercode, module.name, 'headercode') or nil)
	end
	return lines:concat'\n'
end

-- typecode + structs + header + code
function ModuleSet:getCodeAndHeader(...)
	local deps = self:getDependentModules(...)
	local lines = table()
	-- typecode & structs
	for _,module in ipairs(deps) do
		lines:insert(comment(table{
				module.typecode
			}:append(module.structs:mapi(function(struct)
				-- this needs makeType() called first, which generates the .typecode
				-- but it also calls the ffi.metatype (which can only be done once)
				-- and also the ffi.cdef (which is only effective the first time it's done)
				return struct.typecode
			end)):concat'\n', module.name, 'typecode & structs'
		) or nil)
	end
	-- headercode
	for _,module in ipairs(deps) do
		-- allow lazy code generation
		if type(module.headercode) == 'function' then
			module.headercode = module.headercode()
		end

		lines:insert(comment(module.headercode, module.name, 'headercode') or nil)
	end
	-- code
	for _,module in ipairs(deps) do
		-- allow lazy code generation
		if type(module.code) == 'function' then
			module.code = module.code()
		end
		
		lines:insert(comment(module.code, module.name, 'code') or nil)
	end
	return lines:concat'\n'
end

--[[
add using the following markup:
//// MODULE_NAME: (name)
//// MODULE_DEPENDS: space-separated-dependencies

//// MODULE_CODE:
(code)

//// MODULE_HEADER: 
(header)

//// MODULE_TYPE: 
(type)

args:
	code
	onAdd = function(moduleArgs)
--]]
function ModuleSet:addFromMarkup(args)
	if type(args) == 'string' then args = {code = args} end

	local srcLines = string.split(args.code, '\n')

	local name, dstLines, lineTarget, depends
	
	local function resetState()
		name = nil
		dstLines = {
			typecode = table(),
			headercode = table(),
			code = table(),
		}
		lineTarget = 'code'
		depends = table()
	end
	resetState()

	local function makeModule()
		if name then
			local moduleArgs = {
				name = name,
				depends = depends,
				typecode = dstLines.typecode:concat'\n',
				headercode = dstLines.headercode:concat'\n',
				code = dstLines.code:concat'\n',
			}
			if args.onAdd then
				args.onAdd(moduleArgs)
			end
			if string.trim(moduleArgs.headercode) == '' then moduleArgs.headercode = nil end
			if string.trim(moduleArgs.code) == '' then moduleArgs.code = nil end
			self:add(moduleArgs)
		else
			if dstLines.headercode then
				local code = dstLines.headercode:concat'\n'
				if code:match'%S' then
					print('!!! throwing away headercode !!!!:\n'..code)
				end
			end
			if dstLines.code then
				local code = dstLines.code:concat'\n'
				if code:match'%S' then
					print('!!! throwing away code !!!!:\n'..code)
				end
			end
		end
	
		resetState()
	end
	for _,line in ipairs(srcLines) do
		local readname = line:match'^//// MODULE_NAME: (.*)'
		if readname then
			makeModule()
			name = string.trim(readname)
			assert(#name > 0, "got an empty module name")
		else
			local deps = line:match'^//// MODULE_DEPENDS: (.*)'
			if deps then
				depends:append(string.split(string.trim(deps), ' '))
			else
				if string.trim(line) == '//// MODULE_CODE:' then
					lineTarget = 'code'
				elseif string.trim(line) == '//// MODULE_HEADER:' then
					lineTarget = 'headercode'
				elseif string.trim(line) == '//// MODULE_TYPE:' then
					lineTarget = 'typecode'
				else
					table.insert(dstLines[lineTarget], line)
				end
			end
		end
	end
	makeModule()
end

return ModuleSet
