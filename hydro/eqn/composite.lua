local table = require 'ext.table'
local class = require 'ext.class'
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

	self.submodules = self.eqns:mapi(function()
		return require 'hydro.code.moduleset'()
	end)

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
	for _,eqn in ipairs(self.eqns) do
		eqn.solver = setmetatable({}, {
			__index = solver,
		})
		eqn.initCond = {
			initCodeModules = function(eqn, solver)
			end,
			getInitCondCode = function(eqn, solver)
				return ''
			end,
		}
		--eqn:initCodeModules()
		-- assign the cons_t, prim_t, and eigen_t to app
		-- so that, when querying modules in the composite's cdefAllVarTypes, it looks in app's modules and finds them
		rawset(eqn.solver, 'modules', solver.app.modules)
		eqn:initCodeModule_cons_prim_eigen()
	end

	-- now set the subeqns' modules to our local stored copies
	-- so we can store them for later and use them as needed by the composite 
	for i,eqn in ipairs(self.eqns) do
		rawset(eqn.solver, 'modules', self.submodules[i])
	end

	self.consVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.cons_t), name='eqn'..i}
	end)

	self.primVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.prim_t), name='eqn'..i}
	end)
	
	self.eigenVars = self.eqns:mapi(function(eqn, i)
		return {type=assert(eqn.eigen_t), name='eqn'..i}
	end)

	self.numWaves = self.eqns:mapi(function(eqn)
		return eqn.numWaves
	end):sum()

	CompositeEquation.super.init(self, args)
end

function CompositeEquation:getModuleDepends_cons_t()
	return self.eqns:mapi(function(eqn)
		return eqn.cons_t
	end)
end

function CompositeEquation:createInitState()
	CompositeEquation.super.createInitState(self)
	for _,eqn in ipairs(self.eqns) do
		eqn.initCond = self.initCond
	end
end

function CompositeEquation:initCodeModules()
	-- first build the submodules' code
	-- store them locally in composite's submodules[]
	for _,eqn in ipairs(self.eqns) do
		eqn:initCodeModules()
	end
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

function CompositeEquation:getDisplayVars()
	return table():append(
		self.eqns:mapi(function(eqn)
			return eqn:getDisplayVars()
		end):unpack()
	)
end

return CompositeEquation
