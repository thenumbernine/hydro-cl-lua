local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local HLL = Flux:subclass()

HLL.name = 'hll'
HLL.solverCodeFile = 'hydro/flux/hll.clcpp'

--[[
hllCalcWaveMethod:
0 = Davis direct
1 = Davis direct bounded
--]]
function HLL:init(args)
	self.hllCalcWaveMethod = args.hllCalcWaveMethod or 1

	HLL.super.init(self, args)
	
	self.solver.solverStruct.vars:append{
		{name='flux_hllCalcWaveMethod', type='int'},
	}
end

function HLL:initCodeModules()
	HLL.super.initCodeModules(self)
	self.solver.solverPtr.flux_hllCalcWaveMethod = self.hllCalcWaveMethod
end


return HLL
