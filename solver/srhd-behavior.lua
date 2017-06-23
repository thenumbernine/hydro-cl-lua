local class = require 'ext.class'
local SRHDSelfGrav = require 'solver.srhd-selfgrav'
local SRHDEqn = require 'eqn.srhd'
local RHDBehavior = require 'solver/rhd-behavior'

return function(parent)
	local template = class(RHDBehavior(SRHDSelfGrav(parent)))

	function template:createEqn()
		self.eqn = SRHDEqn(self)
	end

	return template
end
