local class = require 'ext.class'
local NoDiv = require 'solver.nodiv'

local function MaxwellBehavior(parent)
	local template = class(NoDiv(parent))

	function template:createEqn()
		self.eqn = require 'eqn.maxwell'(self)
	end

	return template
end

return MaxwellBehavior
