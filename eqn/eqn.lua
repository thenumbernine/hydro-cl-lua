local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'

local Equation = class()

Equation.hasEigenCode = nil
Equation.hasCalcDT = nil
Equation.hasFluxFromCons = nil
Equation.useSourceTerm = nil

-- TODO pass this to Equation from some way other than solver
Equation.prim_t = 'prim_t'
Equation.cons_t = 'cons_t'
Equation.consLR_t = 'consLR_t'

function Equation:init(solver)
	self.solver = assert(solver)

	-- TODO get rid of consVars and primVars
	-- they're still used by ADM1D
	-- default # states is # of conservative variables
	if not self.numStates then 
		self.numStates = #self.consVars 
	else
		if self.consVars then
			assert(self.numStates == #self.consVars)
		end
	end
	-- default # waves is the # of states
	if not self.numWaves then self.numWaves = self.numStates end 

	self.initStateNames = table.map(self.initStates, function(info) return info.name end)
	self.guiVarsForName = table.map(self.guiVars, function(var) return var, var.name end)
end

function Equation:getCodePrefix()
	return (self.guiVars and table.map(self.guiVars, function(var) 
		return var:getCode()
	end) or table()):concat'\n'
end

function Equation:getTypeCode()
	assert(self.consVars)
	return table{
		'typedef union {',
		'	real ptr['..self.numStates..'];',
		'	struct {',
		'		real '..table.concat(self.consVars, ', ')..';',
		'	};',
		'} '..self.cons_t..';',
	}:concat'\n'
end

function Equation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
]], {
	eqn = self,
})
end

-- TODO autogen the name so multiple solvers don't collide
function Equation:getEigenTypeCode(solver)
	return template([[
typedef struct {
	real evL[<?=numStates * numWaves?>];
	real evR[<?=numStates * numWaves?>];
<? if solver.checkFluxError then ?>
	real A[<?=numStates * numStates?>];
<? end ?>
} eigen_t;
]], {
		numStates = self.numStates,
		numWaves = self.numWaves,
		solver = solver,
	})
end

function Equation:getEigenCode(solver)
	if self.hasEigenCode then return end
	return template(file['solver/eigen.cl'], {solver=solver})
end

function Equation:getEigenDisplayVars(solver)
	return range(self.numStates * self.numWaves):map(function(i)
		local row = (i-1)%self.numWaves
		local col = (i-1-row)/self.numWaves
		return {['evL_'..row..'_'..col] = 'value = eigen->evL['..i..'];'}
	end):append(range(self.numStates * self.numWaves):map(function(i)
		local row = (i-1)%self.numStates
		local col = (i-1-row)/self.numStates
		return {['evR_'..row..'_'..col] = 'value = eigen->evR['..i..'];'}
	end)):append(solver.checkFluxError and range(self.numStates * self.numStates):map(function(i)
		local row = (i-1)%self.numStates
		local col = (i-1-row)/self.numStates
		return {['A_'..row..'_'..col] = 'value = eigen->A['..i..'];'}
	end) or nil)
end

return Equation
