local class = require 'ext.class'
local table = require 'ext.table'
local CLBuffer = require 'cl.obj.buffer'
local makestruct = require'eqn.makestruct'
local ffi = require 'ffi'

--[[
TODO combine this with the solver_t stuff in SolverBase

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
	self.solver = assert(args.solver)
	self.name = assert(args.name)
	self.vars = table()
	if args.vars then
		self:addVars(args.vars)
	end
end

function Struct:makeType()
	assert(not self.typename)
	if not self.typename then
		self.typename = self.solver.app:uniqueName(self.name)
	end
	makestruct.safeFFICDef(self:getTypeCode())
end

function Struct:getTypeCode()
	local code = makestruct.makeStruct(self.typename, self.vars, nil, true)
	if self.typecode then
		assert(code == self.typecode)
	else
		self.typecode = code
	end
	return code
end

return Struct
