local class = require 'ext.class'
local SRHDEqn = require 'eqn.srhd'
local RHDBehavior = require 'solver/rhd-behavior'

return function(parent)
	local template = class(RHDBehavior(parent))

	function template:createEqn()
		self.eqn = SRHDEqn(self)
	end

	return template
end
