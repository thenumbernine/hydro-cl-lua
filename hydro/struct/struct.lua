local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local CLBuffer = require 'cl.obj.buffer'

--[[
We have a few different variable sets used for structs, gui, display, etc
Here's the kinds of things different ones hold:
*) name
*) type
*) units
*) code (for calculating display var values)
*) compile-time or not

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
	self.solver = args.solver
	self.name = assert(args.name)
	self.vars = table(args.vars)
	self.dontUnion = args.dontUnion
end

-- static
function Struct.makeStruct(name, vars, scalar, dontUnion)
	scalar = scalar or 'real'
	local numScalars = require 'hydro.struct.struct'.countScalars({vars=vars}, scalar)

	local lines = table()
	local tab
	if dontUnion then
		lines:insert'typedef struct {'
		tab = '\t'
	else
		lines:insert'typedef union {'
		lines:insert('	'..scalar..' ptr['..numScalars..'];')
		lines:insert('	struct {')
		tab = '\t\t'
	end	
	for _,var in ipairs(vars) do
		
		lines:insert(
			tab
			..var.type
			
			-- fixing 'half' and 'double' alignment in solver_t
			-- dontUnion is only used by solver_t
			-- and solver_t is the only one with this C/CL alignment problem
			..(dontUnion and ' __attribute__ ((packed))' or '')
			
			..' '..var.name..';')
	end
	if not dontUnion then
		lines:insert('	};')
	end
	lines:insert('} '..name..';')
	return lines:concat'\n'
end

-- static
function Struct.safeFFICDef(code)
	xpcall(function()
		ffi.cdef(code)
	end, function(msg)
		print(require 'template.showcode'(code))
		io.stderr:write(msg..'\n'..debug.traceback())
		os.exit(1)
	end)
end

function Struct:makeType()
	assert(not self.typename, "don't call makeType() twice")
	if self.solver then
		self.typename = self.solver.app:uniqueName(self.name)
	else
		self.typename = self.name
	end
	Struct.safeFFICDef(self:getTypeCode())

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
				local ctype = field.type
				
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
				local ctype = field.type
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

	local sizeOfFields = table.map(struct.vars, function(field)
		local name = field.name
		local ctype = field.type
		return ffi.sizeof(ctype)
	end):sum()
	assert(ffi.sizeof(self.typename) == sizeOfFields, "struct "..self.typename.." isn't packed!")
	
	self.metatype = metatype
end

function Struct:getTypeCode()
	if not self.typename then
		self:makeType()
	end
	local code = Struct.makeStruct(self.typename, self.vars, nil, self.dontUnion)
	if self.typecode then
		assert(code == self.typecode)
	else
		self.typecode = code
	end
	return code
end

function Struct:countScalars(scalar)
	scalar = scalar or 'real'
	local structSize = 0
	for _,var in ipairs(self.vars) do
		structSize = structSize + ffi.sizeof(var.type)
	end
	local numScalars = structSize / ffi.sizeof(scalar)
	return numScalars
end

return Struct
