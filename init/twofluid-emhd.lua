local InitState = require 'init.init'
local table = require 'ext.table'
return table{
	{
		name = 'Brio-Wu',
		initState = function(solver)
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value = 2
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
}:map(InitState)
