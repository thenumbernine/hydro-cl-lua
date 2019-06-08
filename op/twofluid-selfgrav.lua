local class = require 'ext.class'
local SelfGrav = require 'op.selfgrav'
local TwoFluidSelfGrav = class(SelfGrav)

function TwoFluidSelfGrav:getPoissonDivCode()
	return [[
	source = solver->gravitationalConstant
		* unit_m3_per_kg_s2
		* (U->ion_rho + U->elec_rho);
]]
end

return TwoFluidSelfGrav 
