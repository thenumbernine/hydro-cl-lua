local class = require 'ext.class'
local range = require 'ext.range'
local file = require 'ext.file'
local processcl = require 'processcl'

local Equation = class()

Equation.hasEigenCode = nil
Equation.hasCalcDT = nil
Equation.useSourceTerm = nil

function Equation:init(solver)
	self.solver = assert(solver)

	-- default # states is # of conservative variables
	if not self.numStates then 
		self.numStates = #self.consVars 
	else
		assert(self.numStates == #self.consVars)
	end
	-- default # waves is the # of states
	if not self.numWaves then self.numWaves = self.numStates end 
end

function Equation:getCodePrefix()
	return (self.guiVars and self.guiVars:map(function(var) 
		return var:getCode()
	end) or table()):concat'\n'
end

function Equation:getTypeCode()
	return require 'eqn.makestruct'('cons_t', self.consVars)
end

Equation.getCalcDisplayVarCode = nil

-- TODO autogen the name so multiple solvers don't collide
function Equation:getEigenTypeCode(solver)
	return processcl([[
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
	return processcl(file['solver/eigen.cl'], {solver=solver})
end

function Equation:getEigenDisplayVars(solver)
	return  range(self.numStates * self.numWaves):map(function(i)
		local row = (i-1)%self.numWaves
		local col = (i-1-row)/self.numWaves
		return 'evL_'..row..'_'..col
	end):append(range(self.numStates * self.numWaves):map(function(i)
		local row = (i-1)%self.numStates
		local col = (i-1-row)/self.numStates
		return 'evR_'..row..'_'..col
	end)):append(solver.checkFluxError and range(self.numStates * self.numStates):map(function(i)
		local row = (i-1)%self.numStates
		local col = (i-1-row)/self.numStates
		return 'A_'..row..'_'..col
	end) or nil)
end

function Equation:getEigenCalcDisplayVarCode()
	return [[
	int k = displayVar - displayFirst_eigen;
	if (k < numStates * numWaves) {
		value = eigen->evL[k];
	} else {
		k -= numStates * numWaves;
		if (k < numStates * numWaves) {
			value = eigen->evR[k];
		} else {
			value = eigen->A[k];
		}
	}
]]
end

return Equation
