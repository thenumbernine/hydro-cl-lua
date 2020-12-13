local table = require 'ext.table'
local class = require 'ext.class'
local file = require 'ext.file'
local Struct = require 'hydro.code.struct'
local Equation = require 'hydro.eqn.eqn'

local CompositeEquation = class(Equation)

CompositeEquation.name = 'composite'

CompositeEquation.initConds = require 'hydro.init.euler':getList()

--[[
args:
	subeqns = equations (instances of 'hydro.eqn.eqn') to combine into one
--]]
function CompositeEquation:init(args)
	local solver = assert(args.solver)
	self.eqns = table(assert(args.subeqns))
assert(#self.eqns > 0, "you need at least one entry in args.subeqns")
	for i=1,#self.eqns do
		if type(self.eqns[i]) == 'string' then
			self.eqns[i] = require('hydro.eqn.'..self.eqns[i]){
				solver = assert(solver, "expected solver"),
			}
		-- else see if it's a class
		-- else, object? how do you pass objects, when objects need solver as an arg, and solver needs eqn as an arg
		else
			error("couldn't create sub-eqn "..i.." from "..require 'ext.tolua'(self.eqns[i]))
		end
	end

	-- TODO
	-- another pain point: how to handle ops?
	-- since ops reference fields inside the eqn structs
	-- and now all our eqn structs are merged / require prefixes
	-- and in some cases (like ePot), we don't want to share them necessarily
	-- or if we do then we want the ops to add them all together
	-- so for now:
	solver.ops = table()

--[[
	self.submodules = self.eqns:mapi(function()
		return require 'hydro.code.moduleset'()
	end)
--]]

	--[[
	initCodeModules of subeqns needs to be called here before super.init
	 which calls cdefAllVarTypes, and needs the sub-eqns' cons_t's defined
	actually in this case we need to add it to the app's modules, not the solver's modules
	and this begs the question, why separate the two?
	why not just replace all solver-specific modules with unique names?
	well, that is tedious.
	in the mean time, just add cons_t's to the app
	
	ok ... 
	I have to call initCodeModules here, to get the sub-eqn code modules into hydro.app
	because hydro.app is where cdefAllVarTypes gets its module code from, in order to count reals used in the cons_t
	but initCodeModules depends on createInitState being called
	createInitState isn't called until after the init()

	so if I really want to call sub-eqn initCodeModules here then I also need to give them an initCond somehow 
	... or convince them to ingore their own initCond initCodeModule call
	... and what are the implications of that?  nothing ... since the only initCond we need is from the CompositeEquation, right?  right?
	
	in fact, the initCond stuff is one of a few points of contention of the whole CompositeEquation.  here's the list:
	- initCond
	- addSource
	- constrainU
	
	so lets change the subeqns to have a solver whose .modules really points to app
	but solver can't be changed, it's the same as this solver
	so replace with the a forwrarding meta
	--]]
	for i,eqn in ipairs(self.eqns) do
		-- alright I don't just want to override modules:add => app.modules:add
		-- I also want to do so to solver.modules:add as well
		-- I can either manually copy the modules after-the-fact (modify callback)
		-- or if I want to modify the mt then I'll have to make a new solver.modules (via override field)
		-- which itself forwards everything to one of the modules, except :add, which adds to both
		eqn.solver = setmetatable({
			modules = solver.app.modules,
--[[
			modules = {
				add = function(fakesubmodules, args)
					solver.app.modules:add(args)
				end,
				addFromMarkup = function(fakesubmodules, args)
					solver.app.modules:addFromMarkup(args)
				end,
			},
--]]
		}, {
			__index = solver,
		})
		eqn.initCond = {
			initCodeModules = function(fakeInitCond, solver)
			end,
			getInitCondCode = function(fakeInitCond, solver)
				return ''
			end,
		}
		
		-- this is to prevent initCond_t from being re-added
		eqn.createInitState_createInitState = function() end
		
		--eqn:initCodeModules()
		-- assign the cons_t, prim_t, and eigen_t to app
		-- so that, when querying modules in the composite's cdefAllVarTypes, it looks in app's modules and finds them
		--rawset(eqn.solver, 'modules', solver.app.modules)
		eqn:initCodeModule_cons_prim_eigen()
		-- hack here: remove the typedefs from the cons_t and prim_t
		solver.app.modules.set[eqn.symbols.cons_t].headercode = ''
		solver.app.modules.set[eqn.symbols.prim_t].headercode = ''
		solver.app.modules.set[eqn.symbols.eigen_t].headercode = ''
	end

	-- now set the subeqns' modules to our local stored copies
	-- so we can store them for later and use them as needed by the composite 
--[[
	for i,eqn in ipairs(self.eqns) do
		rawset(eqn.solver, 'modules', self.submodules[i])
	end
--]]
-- [[
	for i,eqn in ipairs(self.eqns) do
		rawset(eqn.solver, 'modules', nil)
	end
--]]

	self.consVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.symbols.cons_t), name='eqn'..i}
	end)

	self.primVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.symbols.prim_t), name='eqn'..i}
	end)
	
	self.eigenVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.symbols.eigen_t), name='eqn'..i}
	end)

	self.numWaves = self.eqns:mapi(function(eqn)
		return eqn.numWaves
	end):sum()

	CompositeEquation.super.init(self, args)
end

function CompositeEquation:getModuleDepends_cons_t()
	return self.eqns:mapi(function(eqn)
		return eqn.symbols.cons_t
	end)
end

function CompositeEquation:createInitState()
	CompositeEquation.super.createInitState(self)
--[[ but we need the fake initCond during initCodeModules to prevent initCond_t from being re-added	
	for _,eqn in ipairs(self.eqns) do
		eqn.initCond = self.initCond
	end
--]]

-- [[ now we have to deal with sub-eqn solver_t fields
-- but we don't want to re-add the units
	for _,eqn in ipairs(self.eqns) do
		eqn:createInitState()
--[=[ but what about overlapping names?		
		self.guiVars:append(eqn.guiVars)
--]=]	
-- [=[
		for _,var in ipairs(eqn.guiVars) do
			if not self.guiVars:find(nil, function(var2) return var.name == var2.name end) then
				self.guiVars:insert(var)
			end
		end
--]=]	
	end
--]]
end

CompositeEquation.solverCodeFile = 'hydro/eqn/composite.cl'

--[=[
function CompositeEquation:initCodeModules()
	local solver = self.solver
	
	self:initCodeModule_cons_prim_eigen()

	solver.modules:add{
		name = self.symbols.waves_t,
		depends = {'real'},
		typecode = self:template[[
typedef union { 
	real ptr[<?=numWaves?>]; 
} <?=waves_t?>;
]],
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.waves_t..' waves_t;',
	}

-- TODO from here on down, it is a lot in common with eqn/eqn.lua

	self:initCodeModule_calcDT()
	self:initCodeModule_fluxFromCons()
	self:initCodeModule_waveCode()

	solver.modules:addFromMarkup{
		code = self:template(file[self.solverCodeFile]),
		onAdd = function(args)
			-- special case for applyInitCond ...
			if args.name == self.symbols.applyInitCond then
				args.depends:append(self.initCond:getBaseDepends(solver))
				args.depends:append(self.initCond.depends)
				-- only used by hydro/eqn/bssnok-fd.lua:
				if self.getModuleDependsApplyInitCond then
					args.depends:append(self:getModuleDependsApplyInitCond())
				end
			end
		end,
	}

	solver.modules:add{
		name = self.symbols.applyInitCond,
		depends = {
			self.symbols.cell_t,
			self.symbols.initCond_t,
			self.symbols.applyInitCondCell,
			'SETBOUNDS',
		},
		code = self:template[[
kernel void <?=applyInitCond?>(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?>* UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	SETBOUNDS(0,0);
	global <?=cons_t?> * const U = UBuf + index;
	global <?=cons_t?> const * const cell = cellBuf + index;
	<?=applyInitCondCell?>(solver, initCond, U, cell);
}
]],
	}
end
--]=]
function CompositeEquation:initCodeModules()
	-- build the submodules' code
	-- store them locally in composite's submodules[]
	-- but don't initCodeModule_cons_prim_eigen because we already did that (for app.modules for cdefAll...)
	for _,eqn in ipairs(self.eqns) do
		local old = eqn.initCodeModule_cons_prim_eigen
		eqn.initCodeModule_cons_prim_eigen = function() end
		eqn:initCodeModules()
		eqn.initCodeModule_cons_prim_eigen = old
		self.solver.modules.set[eqn.symbols.waves_t].headercode = ''
	end

	CompositeEquation.super.initCodeModules(self)
end

function CompositeEquation:initCodeModule_cons_parallelPropagate()
end

function CompositeEquation:initCodeModule_fluxFromCons()
end


--[[
next issue: implementing functions/modules
I will want to just copy the sub-equations' module code into here
but I also don't want namespace collisions
which means ... suffix on the code?  or grep it away?

next: useAddSource?  TODO remove this and replace it with detection of the 'addSource' module
and for composite equations, just check all equation modules. 
--]]
function CompositeEquation:initCodeModule_consFromPrim_primFromCons()
end

function CompositeEquation:addDisplayVarInfosForType(args)
	local eqnAndStructForType = {}
	for _,eqn in ipairs(self.eqns) do
		eqnAndStructForType[eqn.symbols.cons_t] = {eqn = eqn, struct = eqn.consStruct}
		eqnAndStructForType[eqn.symbols.eigen_t] = {eqn = eqn, struct = eqn.eigenStruct}
	end
	local eqnAndStruct = eqnAndStructForType[args.type]
	if eqnAndStruct then
		local eqn = eqnAndStruct.eqn
		local struct = eqnAndStruct.struct
		local vars = eqn:getDisplayVarsForStructVars(struct.vars)
		return table.unpack(vars)
	else
		return CompositeEquation.super.addDisplayVarInfosForType(self, args)
	end
end

-- calcDT ... this should be the min of all sub-calcDT's
-- next big issue is that we need a unique name for each function 
-- this is where classes come in handy
-- why can't OpenCL use C++?

function CompositeEquation:getModuleDepends_waveCode()
	return table():append(
		self.eqns:mapi(function(eqn)
			return eqn:getModuleDepends_waveCode()
		end):unpack()
	)
end

function CompositeEquation:getEigenDisplayVars()
	return table():append(
		self.eqns:mapi(function(eqn)
			return eqn:getEigenDisplayVars() or {}
		end):unpack()
	)
end

-- TODO - prevent variable collisions - especially from multiple matching subeqns
-- this might require some kind of namespace
function CompositeEquation:eigenWaveCodePrefix(n, eig, x)
	return self.eqns:mapi(function(eqn,i)
		return eqn:eigenWaveCodePrefix(n, '&('..eig..')->eqn'..i, x)
	end):concat'\n'
end

function CompositeEquation:eigenWaveCode(n, eig, x, waveIndex)
	local origWaveIndex = waveIndex
	for i,eqn in ipairs(self.eqns) do
		if waveIndex >= 0 and waveIndex < eqn.numWaves then
			return eqn:eigenWaveCode(n, '&('..eig..')->eqn'..i, x, waveIndex)
		end
		waveIndex = waveIndex - eqn.numWaves
	end
	error("couldn't find waveIndex "..origWaveIndex.." in any sub-eqns")
end

-- TODO same as eigenWaveCodePrefix
function CompositeEquation:consWaveCodePrefix(n, U, x)
	return self.eqns:mapi(function(eqn,i)
		return eqn:consWaveCodePrefix(n, '&('..U..')->eqn'..i, x)
	end):concat'\n'
end

function CompositeEquation:consWaveCode(n, U, x, waveIndex)
	local origWaveIndex = waveIndex
	for i,eqn in ipairs(self.eqns) do
		if waveIndex >= 0 and waveIndex < eqn.numWaves then
			return eqn:consWaveCode(n, '&('..U..')->eqn'..i, x, waveIndex)
		end
		waveIndex = waveIndex - eqn.numWaves
	end
	error("couldn't find waveIndex "..origWaveIndex.." in any sub-eqns")
end

-- TODO eigenMinWaveCode
-- TODO eigenMaxWaveCode
-- TODO consMinWaveCode
-- TODO consMaxWaveCode

return CompositeEquation
