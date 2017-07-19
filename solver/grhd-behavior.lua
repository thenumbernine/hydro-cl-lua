local class = require 'ext.class'
local GRHDEqn = require 'eqn.grhd'
local RHDBehavior = require 'solver/rhd-behavior'

return function(parent)
	local templateClass = class(RHDBehavior(parent))
	
	function templateClass:createEqn()
		self.eqn = GRHDEqn(self)
	end
	
	return templateClass
end
