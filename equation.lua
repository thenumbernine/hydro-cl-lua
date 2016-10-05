local class = require 'ext.class'

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
	-- default size of eigen transform struct is two matrices, one from state->wave, one from wave->state
	if not self.numEigen then
		self.numEigen = 2 * self.numWaves * self.numStates
	end
end

function Equation:getTypeCode()
	return require 'makestruct'('cons_t', self.consVars)
end

Equation.eigenType = 'eigen_t'
function Equation:getEigenTypeCode()
	return 'typedef struct { real evL[' .. (self.numStates * self.numWaves) .. '], evR[' .. (self.numStates * self.numWaves) .. ']; } ' .. self.eigenType .. ';'
end

function Equation:getEigenCode()
	return '#include "eigen.cl"'
end


return Equation
