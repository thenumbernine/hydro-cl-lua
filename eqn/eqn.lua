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

	self.waves_t = self:unique'waves_t'

	local numReals
	if self.consVars then
		numReals = makestruct.countReals(self.consVars)
		if self.primVars then
			local numPrimReals = makestruct.countReals(self.primVars)
			assert(numPrimReals <= numReals, "hmm, this is awkward")
		end
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
	local mt = getmetatable(self)
	if mt.guiVars then
		self:addGuiVars(mt.guiVars)
	end
	
	-- first create the init state
	assert(self.initStates, "expected Eqn.initStates")
	self.initState = self.initStates[self.solver.initStateIndex](self.solver, self.solver.initStateArgs)
	assert(self.initState, "couldn't find initState "..self.solver.initStateIndex)

	-- then setup the gui vars
	if self.initState.guiVars then
		self:addGuiVars(self.initState.guiVars)
	end
end

function Equation:getCodePrefix()
	return (self.guiVars and table.map(self.guiVars, function(var) 
		return var:getCode()
	end) or table()):concat'\n'
end

function Equation:getExtraTypeCode()
	return template([[
typedef struct { 
	real ptr[<?=eqn.numWaves?>]; 
} <?=eqn.waves_t?>;
]],		{
			eqn = self,
		})
end

function Equation:getTypeCode()
	assert(self.consVars)
	local lines = table{
		makestruct.makeStruct(self.cons_t, self.consVars),
	}
	if self.primVars then
		lines:insert(makestruct.makeStruct(self.prim_t, self.primVars))
	else
		lines:insert('typedef '..self.cons_t..' '..self.prim_t..';')
	end
	return lines:concat'\n'
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

-- accepts a list of struct var info {name=type}
-- returns a list of display var construction info
function Equation:getDisplayVarsForStructVars(structVarInfos, ptrName)
	ptrName = ptrName or 'U'
	ptrName = ptrName .. '->'
	local displayVarInfos = table()
	for _,structVarInfo in ipairs(structVarInfos) do
		local varname, vartype = next(structVarInfo)
		if vartype == 'real' then
			displayVarInfos:insert{[varname] = '*value = '..ptrName..varname..';'}
		elseif vartype == 'real3' then
			displayVarInfos:insert{[varname] = '*valuevec = '..ptrName..varname..';', type='real3'}
		elseif vartype == 'sym3' then
			displayVarInfos:insert{[varname] = '*valuesym3 = '..ptrName..varname..';', type='sym3'}
		elseif vartype == '_3sym3' then
			for i,xi in ipairs(xNames) do
				displayVarInfos:insert{[varname..'_'..xi] = '*valuesym3 = '..ptrName..varname..'.'..xi..';', type='sym3'}
			end
		end
	end
	return displayVarInfos	
end

function Equation:getDisplayVars()
	return self:getDisplayVarsForStructVars(self.consVars)
end

-- TODO autogen the name so multiple solvers don't collide
function Equation:getEigenTypeCode()
	if self.eigenVars then
		return makestruct.makeStruct(self.eigen_t, self.eigenVars)
	
	-- use the default matrix structures	
	-- whose code is in solver/eigen.cl (included below)
	else
		return template([[
typedef struct {
	real evL[<?=numIntStates * numWaves?>];
	real evR[<?=numIntStates * numWaves?>];
<? if solver.checkFluxError then ?>
	real A[<?=numIntStates * numIntStates?>];
<? end ?>
} <?=eqn.eigen_t?>;
]], 	{
			numIntStates = self.numIntStates,
			numWaves = self.numWaves,
			solver = self.solver,
			eqn = self,
		})
	end
end

function Equation:getEigenCode()
	if self.hasEigenCode then return end
	return template(file['solver/eigen.cl'], {
		solver = self.solver,
		eqn = self,
	})
end

function Equation:getEigenDisplayVars()
	-- use the automatic codegen for display vars
	if self.eigenVars then
		return self:getDisplayVarsForStructVars(self.eigenVars, 'eigen')

	-- use the autogen left & right eigenvector matrices
	else
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
end

function Equation:eigenWaveCodePrefix(side, eig, x)
	return ''
end

function Equation:resetState()
	self.initState:resetState(self.solver)
end

return Equation
