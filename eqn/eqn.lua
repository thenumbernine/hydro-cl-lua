local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local makestruct = require 'eqn.makestruct'

local Equation = class()

Equation.hasEigenCode = nil
Equation.hasCalcDT = nil
Equation.hasFluxFromCons = nil
Equation.useSourceTerm = nil

-- static so that no names overlap across all equations
local uid = 0
function Equation:unique(name)
	uid = uid + 1
	return name..'_'..uid
end

function Equation:init(solver)
	self.solver = assert(solver)

	self.prim_t = self:unique'prim_t'
	self.cons_t = self:unique'cons_t'
	self.consLR_t = self:unique'consLR_t'
	self.eigen_t = self:unique'eigen_t'

	local numReals
	if self.consVars then
		numReals = makestruct.countReals(self.consVars)
	end

	if not self.numStates then 
		self.numStates = numReals
		if not self.numStates then
			error("you either need to define numStates or consVars")
		end
	else
		if numReals then
			assert(self.numStates == numReals)
		end
	end
	-- default # waves is the # of states
	if not self.numWaves then self.numWaves = self.numStates end 

	-- how many states are integratable
	-- (put static states at the end of your cons_t structures)
	if not self.numIntStates then self.numIntStates = self.numStates end

	self.initStateNames = table.map(self.initStates, function(info) return info.name end)
end

function Equation:addGuiVars(args)
	for _,arg in ipairs(args) do
		self:addGuiVar(arg)
	end
end

function Equation:addGuiVar(args)
	local vartype = args.type
	if not vartype then
		vartype = type(args.value)
	end
	local cl = require('guivar.'..vartype)

	-- no non-strings allowed as names, especially not numbers,
	-- because I'm going to map names into guiVars' keys, and I don't want to overwrite any integer-indexed guiVars
	assert(args.name and type(args.name) == 'string')

	local var = cl(args)
	self.guiVars:insert(var)
	self.guiVars[var.name] = var
end

-- always call super first
function Equation:createInitState()
	self.guiVars = table()
	assert(self.initStates, "expected Eqn.initStates")
	self.initState = self.initStates[self.solver.initStateIndex](self.solver)
	assert(self.initState, "couldn't find initState "..self.solver.initStateIndex)	
end

function Equation:getCodePrefix()
	return (self.guiVars and table.map(self.guiVars, function(var) 
		return var:getCode()
	end) or table()):concat'\n'
end

function Equation:getTypeCode()
	assert(self.consVars)
	return makestruct.makeStruct(self.cons_t, self.consVars)
end

function Equation:getSolverCode()
	return template(file[self.solverCodeFile], {eqn=self, solver=self.solver})
end

function Equation:getInitStateCode()
	return self.initState:getInitStateCode(self.solver)
end

function Equation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
]], {
	eqn = self,
})
end

-- TODO autogen the name so multiple solvers don't collide
function Equation:getEigenTypeCode()
	return template([[
typedef struct {
	real evL[<?=numIntStates * numWaves?>];
	real evR[<?=numIntStates * numWaves?>];
<? if solver.checkFluxError then ?>
	real A[<?=numIntStates * numIntStates?>];
<? end ?>
} <?=eqn.eigen_t?>;
]], {
		numIntStates = self.numIntStates,
		numWaves = self.numWaves,
		solver = self.solver,
		eqn = self,
	})
end

function Equation:getEigenCode()
	if self.hasEigenCode then return end
	return template(file['solver/eigen.cl'], {
		solver = self.solver,
		eqn = self,
	})
end

function Equation:getEigenDisplayVars()
	return range(self.numIntStates * self.numWaves):map(function(i)
		local row = (i-1)%self.numWaves
		local col = (i-1-row)/self.numWaves
		return {['evL_'..row..'_'..col] = '*value = eigen->evL['..i..'];'}
	end):append(range(self.numIntStates * self.numWaves):map(function(i)
		local row = (i-1)%self.numIntStates
		local col = (i-1-row)/self.numIntStates
		return {['evR_'..row..'_'..col] = '*value = eigen->evR['..i..'];'}
	end)):append(self.solver.checkFluxError and range(self.numIntStates * self.numIntStates):map(function(i)
		local row = (i-1)%self.numIntStates
		local col = (i-1-row)/self.numIntStates
		return {['A_'..row..'_'..col] = '*value = eigen->A['..i..'];'}
	end) or nil)
end

function Equation:resetState()
	self.initState:resetState(self.solver)
end

return Equation
