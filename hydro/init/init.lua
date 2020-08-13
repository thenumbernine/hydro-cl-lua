local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local time = table.unpack(require 'hydro.util.time')
local Struct = require 'hydro.code.struct'

--[[
name = name of the initial condition

guiVars = any gui variables that the initial conditions wants to expose
TODO call it 'initCondVars'

solverVars = any of the solver/equation's solver_t vars that the initial conditions wants to override

getInitCondCode = function(self, solver) returns the OpenCL code for the initial conditions
--]]

local InitCond = class()

function InitCond:createInitStruct(solver)
	self.initStruct = Struct{
		solver = solver,
		name = 'initCond_t',
		dontUnion = true,
	}
	
	-- then setup the gui vars
	-- this assumes guiVars are index-keyed
	local initGuiVars = self.guiVars
	-- it in turn also indexes them by their name in self.guiVars
	self.guiVars = table()
	if initGuiVars then
		self:addGuiVars(initGuiVars)
	end
end

function InitCond:finalizeInitStruct(solver)
	-- if initCond_t is empty then allocating the initCondBuf will fail
	-- solve this by:
	-- a) not allocating it if initCond_t is empty
	-- b) allocating it at a minimium of size 1
	-- c) putting a temp field in empty initCond_t structures
	-- I'll pick c) for now
	if #self.initStruct.vars == 0 then
		self.initStruct.vars:insert{name='tmp', type='real'}
	end
	self.initStruct:makeType()
	self.initCond_t = self.initStruct.typename
end

-- [[ TODO this is similar to what's in Equation
function InitCond:addGuiVars(args)
	for _,arg in ipairs(args) do
		self:addGuiVar(arg)
	end
end

function InitCond:addGuiVar(args)
	local vartype = args.type
	if not vartype then
		vartype = type(args.value)
	end
	local cl = require('hydro.guivar.'..vartype)

	-- no non-strings allowed as names, especially not numbers,
	-- because I'm going to map names into guiVars' keys, and I don't want to overwrite any integer-indexed guiVars
	assert(args.name and type(args.name) == 'string')
	local var = cl(args)
	
	-- assumes self.guiVars[i] already points to this var
	self.guiVars:insert(var)
	self.guiVars[var.name] = var

	if not args.compileTime then 
		self.initStruct.vars:insert{
			name = var.name,
			type = var.ctype,
		}
	end
end
--]]

function InitCond:initCodeModules(solver)
	solver.modules:add{
		name = 'initCond-codeprefix',
		code = self.getCodePrefix and self:getCodePrefix(solver) or '',
	}
end

function InitCond:refreshInitStateProgram(solver)
	local initCondCode 
	time('generating init state code', function()
		local codePrefix = table{
			solver.modules:getHeader(solver.solverModuleNames:unpack()),
			solver.modules:getCode(solver.solverModuleNames:unpack()),
		}:concat'\n'

		initCondCode = table{
			codePrefix,
			
			self.initStruct.typecode,
		
			self.guiVars:mapi(function(var,i,t) 
				return (var.compileTime and var:getCode() or nil), #t+1
			end):concat'\n',

			self.header and self:header(solver) or '',
			
			-- this in turn calls self:getInitCondCode() but with proper template args applied
			solver.eqn:getInitCondCode(),
		}:concat'\n'
	end)

	local file = require 'ext.file'
	time('building init cond program', function()
		solver.initCondProgramObj = solver.Program{name='initCond', code=initCondCode}
		solver.initCondProgramObj:compile()
	end)

	-- how to give initCond access to rand()?
	-- fill UBuf with random numbers before calling it
	time('randomizing UBuf...', function()
		local ptr = solver.UBufObj:toCPU()
		for i=0,solver.numCells-1 do
			for j=0,solver.eqn.numStates-1 do
				ptr[i].ptr[j] = math.random()
			end
		end
		solver.UBufObj:fromCPU(ptr)
	end)

	solver.applyInitCondKernelObj = solver.initCondProgramObj:kernel('applyInitCond', solver.solverBuf, solver.initCondBuf, solver.UBuf, solver.cellBuf)

	-- here's an ugly hack ...
	-- I need a custom init state kernel for the GLM_MHD only
	-- and it shares init conditions with a lot of other solvers
	-- so ...
	-- (don't the Einstein solvers also use initDerivs?)
	if require 'hydro.eqn.glm-mhd'.is(solver.eqn) then
		solver.initDerivsKernelObj = solver.initCondProgramObj:kernel('initDerivs', solver.solverBuf, solver.UBuf, solver.cellBuf)
	end
end

function InitCond:getInitCondCode(solver)
	return '//no code from InitCond:getInitCondCode() was provided'
end

-- called when the solver resets
function InitCond:resetState(solver)
	solver.applyInitCondKernelObj()

	if cmdline.printBufs then
		print('init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
	
	if require 'hydro.eqn.glm-mhd'.is(solver.eqn) then
		solver:boundary()
		solver.initDerivsKernelObj()
	end
	solver:boundary()

	if cmdline.printBufs then
		print('post-boundary init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
end

return InitCond
