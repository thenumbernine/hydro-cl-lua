local class = require 'ext.class'
local table = require 'ext.table'
local InitState = require 'init.init'
return table{
	{
		name = 'Brio-Wu',
		initState = function(self, solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 2
			end
			return [[
	ion_rho = lhs ? 1 : .125;
	ion_P = lhs ? 1 : .1;
	elec_rho = lhs ? 1 : .125;
	elec_P = lhs ? 1 : .1;
	B.x = .75;
	B.y = lhs ? 1 : -1;
	B.z = 0;
]]
		end,
	}
}:map(function(cl)
	return class(InitState, cl)
end)
