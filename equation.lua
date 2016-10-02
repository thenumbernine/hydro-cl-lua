local class = require 'ext.class'

local Equation = class()

function Equation:init()
	-- default # states is # of conservative variables
	self.numStates = #self.consVars
	-- default # waves is the # of states
	if not self.numWaves then self.numWaves = self.numStates end 
	-- default size of eigen transform struct is two matrices, one from state->wave, one from wave->state
	if not self.numEigen then
		self.numEigen = 2 * self.numWaves * self.numStates
	end
end

Equation.consType = 'cons_t'
-- separate of header() so it can be executed by ffi.cdef and by OpenCL
function Equation:getTypeCode()
	return require 'makestruct'(self.consType, self.consVars)
end

return Equation
