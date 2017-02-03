local class = require 'ext.class'
local SelfGrav = require 'solver.selfgrav'
local EulerEqn = require 'eqn.euler'

-- Euler equation behavior
-- this isn't quite an Equation because self-gravitation is involved
-- apply to a super-class (roe or roe_implicit_linear for now) to make a Euler solver

local function EulerBehavior(parent)
	local template = class(SelfGrav(parent))

	function template:createEqn()
		self.eqn = EulerEqn(self)
	end

	return template
end

return EulerBehavior
