local class = require 'ext.class'
local SelfGrav = require 'solver.selfgrav'
local NoDiv = require 'solver.nodiv'
local TwoFluidEMHDEqn = require 'eqn.twofluid-emhd'

local function TwoFluidEMHDBehavior(parent)
	local templateClass = class(
		--SelfGrav of ion_rho
		--SelfGrav of elec_rho
		NoDiv(parent))

	function templateClass:createEqn()
		self.eqn = TwoFluidEMHDEqn(self)
	end

	return templateClass
end

return TwoFluidEMHDBehavior
