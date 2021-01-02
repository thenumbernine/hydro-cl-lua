local class = require 'ext.class'
local table = require 'ext.table'
local time = table.unpack(require 'hydro.util.time')
local Struct = require 'hydro.code.struct'

--[[
name = name of the initial condition

guiVars = any gui variables that the initial conditions wants to expose
TODO call it 'initCondVars'

solverVars = any of the solver/equation's solver_t vars that the initial conditions wants to override

getInitCondCode = function(self) returns the OpenCL code for the initial conditions
--]]

local InitCond = class()

function InitCond:init(args)
	self.solver = assert(args.solver, "expected solver")
end

function InitCond:createInitStruct()
	self.initStruct = Struct{
		solver = self.solver,
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

function InitCond:finalizeInitStruct()
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

-- these depends are added to the applyInitCond module in the eqn.initCodeModules function
-- so are any in the initCond.depends table
function InitCond:getBaseDepends()
	local solver = assert(self.solver)
	return {
		-- if an InitCond provides codeprefix, it is for code it expects to reference from within 'applyInitCond()'
		solver.eqn.symbols.initCond_codeprefix,
		-- applyInitCond uses these:
		solver.solver_t,
		self.initCond_t,
		solver.coord.cell_t,
		solver.eqn.symbols.initCond_guiVars_compileTime,
		'INDEX',
		'INDEXV',
		solver.symbols.OOB,
		solver.symbols.SETBOUNDS,
		-- enough use #if dim that i'll put this here:
		solver.symbols.solver_macros,
		-- initCond code is specified in terms of primitives, so if the eqn has prim<->cons then it will be needed
		solver.eqn.symbols.consFromPrim,
	}
end

function InitCond:initCodeModules()
	local solver = assert(self.solver)
	solver.modules:add{
		name = self.initCond_t,
		structs = {self.initStruct},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.initCond_t..' initCond_t;',
	}

	solver.modules:add{
		name = solver.eqn.symbols.initCond_guiVars_compileTime,
		headercode = table.mapi(self.guiVars or {}, function(var,i,t) 
			return (var.compileTime and var:getCode() or nil), #t+1
		end):concat'\n',
	}

	solver.modules:add{
		name = solver.eqn.symbols.initCond_codeprefix,
		depends = {
			self.initCond_t,
		},
		code = self.getCodePrefix and self:getCodePrefix() or nil,
	}
end

function InitCond:refreshInitStateProgram()
	local solver = assert(self.solver)
	local eqn = solver.eqn
	
	solver.initModulesEnabled[eqn.symbols.applyInitCond] = true
	if solver:hasModule(eqn.symbols.initDerivs) then
		solver.initModulesEnabled[eqn.symbols.initDerivs] = true
	end

	local initCondCode 
	time('generating init state code', function()
		local moduleNames = table(solver.sharedModulesEnabled, solver.initModulesEnabled):keys()
print('initCond modules: '..moduleNames:sort():concat', ')
		initCondCode = solver.modules:getCodeAndHeader(moduleNames:unpack())
	end)
	
	time('building initCond program', function()
		solver.initCondProgramObj = solver.Program{name='initCond', code=initCondCode}
		solver.initCondProgramObj:compile()
	end)
	
	solver.applyInitCondKernelObj = solver.initCondProgramObj:kernel(eqn.symbols.applyInitCond, solver.solverBuf, solver.initCondBuf, solver.UBuf, solver.cellBuf)
	
	if solver:hasModule(eqn.symbols.initDerivs) then
		solver.initDerivsKernelObj = solver.initCondProgramObj:kernel(eqn.symbols.initDerivs, solver.solverBuf, solver.UBuf, solver.cellBuf)
	end
end

function InitCond:getInitCondCode()
	return '//no code from InitCond:getInitCondCode() was provided'
end

-- called when the solver resets
function InitCond:resetState()
	local solver = assert(self.solver)

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

	solver.applyInitCondKernelObj()

	if cmdline.printBufs then
		print('init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
	
	if solver:hasModule(solver.eqn.symbols.initDerivs) then
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
