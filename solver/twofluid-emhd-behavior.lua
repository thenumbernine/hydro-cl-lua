local class = require 'ext.class'
local TwoFluidEMHDEqn = require 'eqn.twofluid-emhd'

local function TwoFluidEMHDBehavior(parent)
	local templateClass = class(parent)

	function templateClass:createEqn()
		self.eqn = TwoFluidEMHDEqn(self)
	end

	return templateClass
end

return TwoFluidEMHDBehavior
