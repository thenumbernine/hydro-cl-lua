--[[
TODO use struct-lua here

	Struct:
hydro/eqn/grhd.lua:25:	self.consOnlyStruct = Struct{
hydro/eqn/grhd.lua:35:	self.primOnlyStruct = Struct{
hydro/eqn/glm-mhd.lua:109:	self.roeStruct = Struct{solver=solver, name='roe_t', vars=self.roeVars}
hydro/eqn/srhd.lua:45:	self.consOnlyStruct = Struct{
hydro/eqn/srhd.lua:55:	self.primOnlyStruct = Struct{
hydro/eqn/mhd.lua:94:	self.roeStruct = Struct{solver=solver, name='roe_t', vars=self.roeVars}
hydro/eqn/eqn.lua:145:		self.consStruct = Struct{
hydro/eqn/eqn.lua:205:		self.primStruct = Struct{
hydro/eqn/eqn.lua:247:		self.eigenStruct = Struct{
hydro/eqn/eqn.lua:256:		self.eigenStruct = Struct{solver=solver, name='eigen_t', vars=self.eigenVars}
hydro/eqn/eqn.lua:300:	self.wavesStruct = Struct{
hydro/solver/solverbase.lua:520:	self.solverStruct = Struct{
hydro/solver/solverbase.lua:2702:		local _3sym3Struct = Struct{
hydro/solver/solverbase.lua:3823:		if Struct:isa(typeinfo) then

	consStruct:
hydro/eqn/z4.lua:817:local has_b_ul = eqn.consStruct.vars:find(nil, function(var) return var.name == "b_ul" end)
hydro/eqn/z4.lua:825:local has_B_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end)
hydro/eqn/z4.lua:914:	if self.consStruct.vars:find(nil, function(var) return var.name == 'b_ul' end) then
hydro/eqn/z4.cl:6:local has_beta_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "beta_u" end)
hydro/eqn/z4.cl:7:local has_betaLap_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "betaLap_u" end)
hydro/eqn/z4.cl:8:local has_b_ul = eqn.consStruct.vars:find(nil, function(var) return var.name == "b_ul" end)
hydro/eqn/z4.cl:9:local has_B_u = eqn.consStruct.vars:find(nil, function(var) return var.name == "B_u" end)
hydro/eqn/eqn.lua:118:	getParityVars(self.consStruct.vars, sign, parityVars)
hydro/eqn/eqn.lua:127:	make sure self.consStruct and self.primStruct is defined beforehand
hydro/eqn/eqn.lua:143:	if not self.consStruct then
hydro/eqn/eqn.lua:145:		self.consStruct = Struct{
hydro/eqn/eqn.lua:179:	when we makeType the consStruct
hydro/eqn/eqn.lua:181:	(but eqn:initCodeModule is called after the consStruct type is defined)
hydro/eqn/eqn.lua:189:	self:cdefAllVarTypes(solver, self.consStruct.vars)
hydro/eqn/eqn.lua:191:	self.consStruct:makeType()
hydro/eqn/eqn.lua:193:	self.symbols.cons_t = self.consStruct.typename
hydro/eqn/eqn.lua:200:	self.consStruct.eqn = self
hydro/eqn/eqn.lua:201:	solver.structForType[self.consStruct.typename] = self.consStruct
hydro/eqn/eqn.lua:269:	if self.consStruct.vars then
hydro/eqn/eqn.lua:270:		numReals = self.consStruct:countScalars()
hydro/eqn/eqn.lua:280:			error("you either need to define numStates or consVars or consStruct")
hydro/eqn/eqn.lua:354:		getParityVars(self.consStruct.vars, sign, parityVars, field)
hydro/eqn/eqn.lua:645:		for _,var in ipairs(self.consStruct.vars) do
hydro/eqn/eqn.lua:674:					self.consStruct.vars:find(nil, function(var)
hydro/eqn/eqn.lua:679:					self.consStruct.vars:find(nil, function(var)
hydro/eqn/eqn.lua:698:<?		for _,var in ipairs(eqn.consStruct.vars) do
hydro/eqn/eqn.lua:793:	for _,var in ipairs(self.consStruct.vars) do
hydro/eqn/eqn.lua:799:	assert(self.consStruct)
hydro/eqn/eqn.lua:802:		structs = {self.consStruct:getForModules()},
hydro/eqn/eqn.lua:893:	return self.solver:createDisplayVarArgsForStructVars(self.consStruct.vars)
hydro/eqn/eqn.lua:1317:	local _, var = self.consStruct.vars:find(nil, function(v) return v.name == varname end)
hydro/solver/solverbase.lua:2685:		args.vars = self:createDisplayVarArgsForStructVars(self.eqn.consStruct.vars)
hydro/solver/solverbase.lua:3254:				local vars = self.eqn.consStruct.vars
hydro/solver/solverbase.lua:3255:				local numScalars = self.eqn.consStruct:countScalars()
hydro/solver/gridsolver.lua:697:		for _,var in ipairs(eqn.consStruct.vars) do
hydro/solver/meshsolver.lua:267:	for _,var in ipairs(eqn.consStruct.vars) do

	primStruct:
hydro/eqn/eqn.lua:127:	make sure self.consStruct and self.primStruct is defined beforehand
hydro/eqn/eqn.lua:204:	if not self.primStruct and self.primVars then
hydro/eqn/eqn.lua:211:	if self.primStruct then
hydro/eqn/eqn.lua:214:		self:cdefAllVarTypes(solver, self.primStruct.vars)
hydro/eqn/eqn.lua:217:			self.primStruct:makeType()
hydro/eqn/eqn.lua:219:			return "eqn "..self.name.." primStruct:makeType() failed for type "..self.primStruct.name..'\n'
hydro/eqn/eqn.lua:221:					self.primStruct,
hydro/eqn/eqn.lua:229:		self.symbols.prim_t = self.primStruct.typename
hydro/eqn/eqn.lua:236:		self.primStruct.eqn = self
hydro/eqn/eqn.lua:237:		solver.structForType[self.primStruct.typename] = self.primStruct
hydro/eqn/eqn.lua:271:		if self.primStruct then
hydro/eqn/eqn.lua:272:			local numPrimReals = self.primStruct:countScalars()
hydro/eqn/eqn.lua:357:		getParityVars(self.primStruct.vars, sign, parityVars, field)
hydro/eqn/eqn.lua:808:	if self.primStruct then
hydro/eqn/eqn.lua:811:			structs = {self.primStruct:getForModules()},
hydro/eqn/eqn.lua:1218:	assert(not self.primStruct, "if you're using the default prim<->cons code then you shouldn't have any primStruct")
hydro/init/euler.lua:1041:<? if eqn.primStruct.vars:find(nil, function(var)
hydro/init/euler.lua:1061:<? if eqn.primStruct.vars:find(nil, function(var)

	eigenStruct:

	wavesStruct:

	roeStruct:

	consOnlyStruct:

	primOnlyStruct:

	solverStruct:

--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
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

	-- generate the typecode *except* the typename
	-- then compare it to a map from typecode => typename
	-- if it matches any, use the old typecode, typename, and metatype
	-- otherwise generate a new one

	local codeWithoutTypename = self:getTypeCodeWithoutTypeName()

-- new issue with multi-solvers and modules
-- setting the type code to this typedef also means the module using this code will need to depend on the previous def
-- which I can't communicate just yet
-- so in the mean time, I'll just disable this for now
--[=[
	local app = assert(self.app)
	app.typeInfoForCode = app.typeInfoForCode or {}
	local info = app.typeInfoForCode[codeWithoutTypename]
	if info then
		print('reusing matching struct '..info.typename..' for '..self.name)

		--[[ use cached version
		self.typename = info.typename
		self.typecode = info.typecode
		self.metatype = info.metatype
		return
		--]]
		-- but our module system now expects typenames to be unique names
		-- which means that we can't use a matching typename but instead must do a unique name + a typedef
		-- [[
		self.typename = app:uniqueName(self.name)
		codeWithoutTypename = 'typedef '..info.typename
		self.typecode = codeWithoutTypename .. ' ' .. self.typename .. ';'
		safecdef(assert(self.typecode))
		self.metatype = info.metatype
		return
		--]]
	end
--]=]

	-- TODO no more uniqueName and typename ~= name .. instead force names to be unique
	-- and make them unique before passing them in by appending the lua object uid or something
	if self.dontMakeUniqueName then
		self.typename = self.name
	else
		self.typename = assert(self.app):uniqueName(self.name)
	end
	do
		local typecode = codeWithoutTypename .. ' ' .. self.typename .. ';'
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
		__concat = string.concat,
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

--[=[
	-- store it
	app.typeInfoForCode[codeWithoutTypename] = {
		typename = self.typename,
		typecode = self.typecode,
		metatype = self.metatype,
	}
--]=]
end

-- this gets the type code *except* the typename
-- this way I can compare it to other type codes previously generated,
-- so I only need to generate struct codes once
-- I can use this to reduce the number of ffi.cdef's that are called.
function Struct:getTypeCodeWithoutTypeName()
	local lines = table()
	local scalar = 'real'

	local tab
	if self.dontUnion then
		lines:insert'typedef struct {'
		tab = '\t'
	else
		lines:insert'typedef union {'
		local numScalars = self:countScalars(scalar)
		lines:insert('	'..scalar..' ptr['..math.max(1, math.floor(numScalars))..'];')
		lines:insert('	struct {')
		tab = '\t\t'
	end
	for _,var in ipairs(self.vars) do
		local artype, arsize = var.type:match'(.-)%[(.-)%]'
		if not artype then
			artype = var.type
		end
		lines:insert(
			tab
			..artype

			-- fixing 'half' and 'double' alignment in solver_t
			-- dontUnion is only used by solver_t
			-- and solver_t is the only one with this C/CL alignment problem
			..(self.dontUnion and ' __attribute__ ((packed))' or '')

			..' '..var.name
			..(arsize and '['..arsize..']' or '')
			..';')
	end
	if not self.dontUnion then
		lines:insert('	};')
	end
	lines:insert('}')
	return lines:concat'\n'
end

-- glue function between old hydro struct and new struct-lua which the modules system now works with
function Struct:getForModules()
	return {
		code = self.typecode,
		fielditer = function()
			return coroutine.wrap(function()
				for _,field in ipairs(self.vars) do
					coroutine.yield(field.name, field.type, field)
				end
			end)
		end,
	}
end

return Struct
