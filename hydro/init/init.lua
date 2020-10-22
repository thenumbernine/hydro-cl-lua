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
		name = 'initCond.initCond_t',
		structs = {self.initStruct},
	}

	solver.modules:add{
		name = 'initCond.guiVars.compileTime',
		headercode = table.mapi(self.guiVars or {}, function(var,i,t) 
			return (var.compileTime and var:getCode() or nil), #t+1
		end):concat'\n',
	}

	solver.modules:add{
		name = 'initCond.codeprefix',
		depends = {'initCond.initCond_t'},
		code = self.getCodePrefix and self:getCodePrefix(solver) or nil,
	}

	local eqn = solver.eqn
	
	-- depends of applyInitCond:
	eqn.codeDepends = table{
		-- if an InitCond provides codeprefix, it is for code it expects to reference from within 'applyInitCond()'
		'initCond.codeprefix',
		-- applyInitCond uses these:
		'solver_t',
		'initCond.initCond_t',
		'cons_t',
		'coord.cell_t',
		'initCond.guiVars.compileTime',
		'INDEX', 'INDEXV', 'OOB', 'SETBOUNDS',
		-- initCond code is specified in terms of primitives, so if the eqn has prim<->cons then it will be needed
		'eqn.prim-cons',
	}
	:append(self.depends)
	
	-- TODO get rid of this and switch all over to depmod ... maybe?
	if eqn.getModuleDependsApplyInitCond then
		eqn.codeDepends:append(eqn:getModuleDependsApplyInitCond())
	end

	-- maybe it's wrong to put this code into the solverCodeFile because technically it goes into initCond.cl
	local code = eqn:template(
		require 'ext.file'[eqn.solverCodeFile],
		{
			moduleName = 'applyInitCond',
		}
	)

	solver.modules:add{
		name = 'applyInitCond',
		depends = eqn.codeDepends,
		code = code,
	}
	
	eqn.codeDepends = nil
end

function InitCond:refreshInitStateProgram(solver)

	solver.initModulesEnabled['applyInitCond'] = true

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

	solver.applyInitCondKernelObj = solver.initCondProgramObj:kernel('applyInitCond', solver.solverBuf, solver.initCondBuf, solver.UBuf, solver.cellBuf)

	if solver.eqn.needsInitDerivs then
		solver.initDerivsKernelObj = solver.initCondProgramObj:kernel('initDerivs', solver.solverBuf, solver.UBuf, solver.cellBuf)
	end
end

function InitCond:getInitCondCode(solver)
	return '//no code from InitCond:getInitCondCode() was provided'
end

-- called when the solver resets
function InitCond:resetState(solver)

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
	
	if solver.eqn.needsInitDerivs then
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
