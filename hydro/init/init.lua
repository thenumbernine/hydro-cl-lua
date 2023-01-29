local class = require 'ext.class'
local table = require 'ext.table'
local time = table.unpack(require 'hydro.util.time')
local Struct = require 'hydro.code.struct'
local half = require 'cl.obj.half'

--[[
name = name of the initial condition

guiVars = any gui variables that the initial conditions wants to expose
TODO call it 'initCondVars'

solverVars = any of the solver/equation's solver_t vars that the initial conditions wants to override

getInitCondCode = function(self) returns the OpenCL code for the initial conditions
--]]

local InitCond = class()

--[[
for now it depends on the initCond
but how about I make some stanard overrides?

args:
	solverVars = solverVars to override, overriding initCond and of course overriding solver
	... etc = any args that match any initCond guiVars will automatically override the initCond guiVar value
--]]
function InitCond:init(args)
	self.solver = assert(args.solver, "expected solver")
	self.args = args

	self.solverVars = self.solverVars or {}
	for k,v in pairs(args and args.solverVars or {}) do
		self.solverVars[k] = v
	end
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
	args = table(args)
	-- ok 'args' is the guivar/initcondvar args ...
	-- and 'self.args' is the InitCond ctor args ...
	-- HERE: allow the InitCond ctor args to override any guiVars
	if self.args
	and self.args
	and self.args[args.name] ~= nil
	then
		args.value = self.args[args.name]
	end

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

	-- if it's not a compile-time var then it is for manipulating the initCond struct
	if not args.compileTime then
		-- ... and that means we need to tell the var what struct/field to write when it is modified ...
		function var:refreshInStruct(solver)
			solver.initCondPtr[self.name] = self.value
			solver:refreshInitCondBuf()
		end
		-- ... and add it to the initCond_t
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

	solver.initModulesEnabled.units = true		-- I think this is safe to assume
	solver.initModulesEnabled[eqn.symbols.Equation] = true
	if solver:hasModule(eqn.symbols.initDerivs) then
		solver.initModulesEnabled[eqn.symbols.initDerivs] = true
	end

	local initCondCode
	time('generating init state code', function()
		local moduleNames = table(solver.sharedModulesEnabled, solver.initModulesEnabled):keys()
		if solver.app.verbose then
			print('initCond modules: '..moduleNames:concat', ')
		end
		initCondCode = solver.modules:getCodeAndHeader(moduleNames:unpack())
	end)
	
	time('building program cache/'..solver:getIdent()..'/src/initCond.clcpp ', function()
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
				ptr[i].ptr[j] = half.toreal(math.random())
			end
		end
		solver.UBufObj:fromCPU(ptr)
	end)

	solver.applyInitCondKernelObj()

	if cmdline.printBufs then
		print('init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
	
	if solver.initDerivsKernelObj then
		solver:boundary()
		solver.initDerivsKernelObj()
	end
	solver:boundary()

	if cmdline.printBufs then
		print('post-boundary init UBuf:')
		solver:printBuf(solver.UBufObj)
	end
	
	if solver.constrainUKernelObj then
		-- this calls constrainUKernelObj
		-- and then calls :boundary()
		solver:constrainU()
		if cmdline.printBufs then
			print('post-constrainU init UBuf:')
			solver:printBuf(solver.UBufObj)
		end
	end
end

return InitCond
