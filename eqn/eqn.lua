local class = require 'ext.class'
local range = require 'ext.range'

local Equation = class()

function Equation:init()
	-- default # states is # of conservative variables
	if not self.numStates then 
		self.numStates = #self.consVars 
	else
		assert(self.numStates == #self.consVars)
	end
	-- default # waves is the # of states
	if not self.numWaves then self.numWaves = self.numStates end 
end

function Equation:getTypeCode()
	return require 'eqn.makestruct'('cons_t', self.consVars)
end

function Equation:getEigenInfo()
	-- TODO autogen the name so multiple solvers don't collide
	return {
		typeCode =
			'typedef struct { real evL[' .. (self.numStates * self.numWaves) .. '], evR[' .. (self.numStates * self.numWaves) .. ']; } eigen_t;\n'..
			'typedef struct { real A[' .. (self.numStates * self.numStates) .. ']; } fluxXform_t;',
		code = '#include "solver/eigen.cl"',
		displayVars = range(self.numStates * self.numWaves):map(function(i)
			local row = (i-1)%self.numWaves
			local col = (i-1-row)/self.numWaves
			return 'evL_'..row..'_'..col
		end):append(range(self.numStates * self.numWaves):map(function(i)
			local row = (i-1)%self.numStates
			local col = (i-1-row)/self.numStates
			return 'evR_'..row..'_'..col
		end)),
	}
end

Equation.getCalcDisplayVarCode = nil

function Equation:getCalcDisplayVarEigenCode()
	return [[
	int k = displayVar - displayFirst_eigen;
	if (k < numStates * numWaves) {
		value = eigen->evL[k];
	} else {
		k -= numStates * numWaves;
		if (k < numStates * numWaves) {
			value = eigen->evR[k];
		}
	}
]]
end

return Equation
