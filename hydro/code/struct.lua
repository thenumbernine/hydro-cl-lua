local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
--local CLBuffer = require 'cl.obj.buffer'
local safecdef = require 'hydro.code.safecdef'

--[[
We have a few different variable sets used for structs, gui, display, etc
Here's the kinds of things different ones hold:
*) name
*) type
*) units
*) code (for calculating display var values)
*) compile-time or not
*) class body (now that we're in C++)

This holds ...
*) a typename for the struct in OpenCL code
*) a table of vars, each with 'name', and 'type', and other things

Usage:
initialization:
1) struct = Struct() to create Struct object
2) struct:addVar[s]() to add necessary vars
3) struct:alloc() to allocate a CLBuffer of this type
	(and TODO let the cl.obj.buffer keep its cpu ptr around)

then CLBuffer update fromCPU:
4) update struct.var.value's
5) struct:update() to update buffers
--]]
local Struct = class()

function Struct:init(args)
	self.name = assert(args.name)
	self.vars = table(args.vars)
	self.dontUnion = args.dontUnion
	self.dontMakeUniqueName = args.dontMakeUniqueName
	if not self.dontMakeUniqueName then
		self.app = assert(assert(args.solver).app)
	end
	self.body = args.body
end

-- TODO make this static
function Struct:countScalars(scalar)
	scalar = scalar or 'real'
	local structSize = 0
	for _,var in ipairs(self.vars) do
		local res, err = xpcall(function()
			structSize = structSize + ffi.sizeof(var.type)
		end, function(err)
			return 'ffi.sizeof('..var.type..') failed for var '..require 'ext.tolua'(var)..'\n'
				..tostring(err)..'\n'
				..debug.traceback()
		end)
		if not res then error(err) end
	end
	local numScalars = structSize / ffi.sizeof(scalar)
	return numScalars
end

function Struct:makeType()
	assert(not self.typename, "don't call makeType() twice")

	-- TODO no more uniqueName and typename ~= name .. instead force names to be unique
	-- and make them unique before passing them in by appending the lua object uid or something
	if self.dontMakeUniqueName then
		self.typename = self.name
	else
		self.typename = assert(self.app):uniqueName(self.name)
	end
	
	do
		local typecode = self:getTypeCode(self.typename)
		if self.typecode then
			assert(typecode == self.typecode)
		else
			self.typecode = typecode
		end
	end
	safecdef(assert(self.typecode))

	-- make the metatype here
	local struct = self
	local metatable = {
		toLua = function(self)
			local result = {}
			for _,field in ipairs(struct.vars) do
				local name = field.name
				local ctype = field.type
				local value = self[name]
				if ctype.toLua then
					value = value:toLua()
				end
				result[name] = value
			end
			return result
		end,
		__tostring = function(ptr)
			local t = table()
			for _,field in ipairs(struct.vars) do
				local name = field.name
				local s = tostring(ptr[name])
				t:insert(name..'='..s)
			end
			return struct.typename..'{'..t:concat', '..'}'
		end,
		__concat = function(a,b)
			return tostring(a) .. tostring(b)
		end,
		__eq = function(a,b)
			for _,field in ipairs(struct.vars) do
				local name = field.name
				if a[name] ~= b[name] then return false end
			end
			return true
		end,
	}
	metatable.__index = metatable

	local metatype
	local status, err = xpcall(function()
		metatype = ffi.metatype(self.typename, metatable)
	end, function(err)
		return err..' typename='..self.typename..'\n'..debug.traceback()
	end)
	if not status then error(err) end

	local sizeOfFields = table.mapi(struct.vars, function(field)
		local ctype = field.type
		return ffi.sizeof(ctype)
	end):sum() or 0
	if ffi.sizeof(self.typename) ~= sizeOfFields then
		print("ffi.sizeof("..self.typename..") = "..ffi.sizeof(self.typename))
		print('sizeOfFields = '..sizeOfFields)
		for _,field in ipairs(struct.vars) do
			print("ffi.sizeof("..field.type..' '..self.typename..'.'..field.name..') = '..ffi.sizeof(field.type))
		end
		print('typecode:\n'..self.typecode)
		error("struct "..self.typename.." isn't packed!")
	end

	self.metatype = metatype
end

-- this gets the type code
function Struct:getTypeCode(typename)
	local lines = table()
	local scalar = 'real'
--print('getTypeCode '..typename)
	local tab
	local classType
	if self.dontUnion then
		classType = 'struct'
		lines:insert(classType..' '..typename..' {')
		tab = '\t'
	else
		classType = 'union'
		lines:insert(classType..' '..typename..' {')
		local numScalars = self:countScalars(scalar)
		lines:insert('	'..scalar..' ptr['..math.max(1, math.floor(numScalars))..'];')
		lines:insert('	struct {')
		tab = '\t\t'
	end
	for _,var in ipairs(self.vars) do
		lines:insert(
			tab
			..var.type

			-- fixing 'half' and 'double' alignment in solver_t
			-- dontUnion is only used by solver_t
			-- and solver_t is the only one with this C/CL alignment problem
			..(self.dontUnion and ' __attribute__ ((packed))' or '')

			..' '..var.name..';')
	end
	if not self.dontUnion then
		lines:insert('	};')
	end
	-- for opencl-cpp, not for ffi.cef
	lines:insert'//// BEGIN EXCLUDE FOR FFI_CDEF'
	lines:insert('\t'..typename..'() {}')
	for _,var in ipairs(self.vars) do
		lines:insert('\tconstexpr '..typename..' & set_'..var.name..'('..var.type..' const & value_) { '..var.name..' = value_; return *this; }')
	end
	if not self.dontUnion then
		lines:insert(template([[
	constexpr <?=scalar?> & operator[](int i) { return ptr[i]; }
	constexpr <?=scalar?> const & operator[](int i) const { return ptr[i]; }
	constexpr constant <?=scalar?> & operator[](int i) constant { return ptr[i]; }
	constexpr constant <?=scalar?> const & operator[](int i) constant const { return ptr[i]; }
]], 	{
			scalar = scalar,
		}))
	end
	if self.body then lines:insert(self.body) end
	lines:insert'//// END EXCLUDE FOR FFI_CDEF'
	lines:insert('};')
	lines:insert('typedef '..classType..' '..typename..' '..typename..';')
	return lines:concat'\n'
end

return Struct
