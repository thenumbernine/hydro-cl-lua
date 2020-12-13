local table = require 'ext.table'
local class = require 'ext.class'
local file = require 'ext.file'
local Struct = require 'hydro.code.struct'
local Equation = require 'hydro.eqn.eqn'

local Composite = class(Equation)

Composite.name = 'composite'

-- TODO ... ? list from sub-eqns? unique?
Composite.initConds = require 'hydro.init.euler':getList()

Composite.solverCodeFile = 'hydro/eqn/composite.cl'

--[[
args:
	subeqns = equations (instances of 'hydro.eqn.eqn') to combine into one
--]]
function Composite:init(args)
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
	... and what are the implications of that?  nothing ... since the only initCond we need is from the Composite, right?  right?
	
	in fact, the initCond stuff is one of a few points of contention of the whole Composite.  here's the list:
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
		}
		eqn.initCond.getInitCondCode = function(fakeInitCond, solver)
			return self.initCond:getInitCondCode(solver)
		end
		
		-- this is to prevent initCond_t from being re-added
		eqn.createInitState_createInitState = function() end
		
		--eqn:initCodeModules()
		-- assign the cons_t, prim_t, and eigen_t to app
		-- so that, when querying modules in the composite's cdefAllVarTypes, it looks in app's modules and finds them
		--rawset(eqn.solver, 'modules', solver.app.modules)
		eqn:initCodeModule_cons_prim_eigen_waves()
		-- hack here: remove the typedefs from the cons_t and prim_t
		solver.app.modules.set[eqn.symbols.cons_t].headercode = ''
		solver.app.modules.set[eqn.symbols.prim_t].headercode = ''
		solver.app.modules.set[eqn.symbols.eigen_t].headercode = ''
		solver.app.modules.set[eqn.symbols.waves_t].headercode = ''
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
	
	self.wavesVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.symbols.waves_t), name='eqn'..i}
	end)

	Composite.super.init(self, args)
end

function Composite:getModuleDepends_cons_t()
	return self.eqns:mapi(function(eqn)
		return eqn.symbols.cons_t
	end)
end

function Composite:createInitState()
	Composite.super.createInitState(self)
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

function Composite:initCodeModules()
	-- build the submodules' code
	-- store them locally in composite's submodules[]
	-- but don't initCodeModule_cons_prim_eigen_waves because we already did that (for app.modules for cdefAll...)
	for _,eqn in ipairs(self.eqns) do
		local old = eqn.initCodeModule_cons_prim_eigen_waves
		eqn.initCodeModule_cons_prim_eigen_waves = function() end
		eqn:initCodeModules()
		eqn.initCodeModule_cons_prim_eigen_waves = old
	end

	Composite.super.initCodeModules(self)
end

-- in the solverCodeFile
function Composite:initCodeModule_cons_parallelPropagate() end
function Composite:initCodeModule_fluxFromCons() end
function Composite:initCodeModule_consFromPrim_primFromCons() end
function Composite:initCodeModule_calcDTCell() end

function Composite:getModuleDepends_waveCode()
	return table():append(
		self.eqns:mapi(function(eqn)
			return eqn:getModuleDepends_waveCode()
		end):unpack()
	)
end

-- TODO - prevent variable collisions - especially from multiple matching subeqns
-- this might require some kind of namespace
function Composite:eigenWaveCodePrefix(n, eig, x)
	return self.eqns:mapi(function(eqn,i)
		return eqn:eigenWaveCodePrefix(n, '&('..eig..')->eqn'..i, x)
	end):concat'\n'
end

function Composite:eigenWaveCode(n, eig, x, waveIndex)
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
function Composite:consWaveCodePrefix(n, U, x)
	return self.eqns:mapi(function(eqn,i)
		return eqn:consWaveCodePrefix(n, '&('..U..')->eqn'..i, x)
	end):concat'\n'
end

function Composite:consWaveCode(n, U, x, waveIndex)
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

return Composite
