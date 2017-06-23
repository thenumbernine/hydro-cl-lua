local class = require 'ext.class'
local GRHDEqn = require 'eqn.grhd'
local RHDBehavior = require 'solver/rhd-behavior'

return function(parent)
	local template = class(RHDBehavior(parent))

	function template:createEqn()
		self.eqn = GRHDEqn(self)
	end

	return template
end

