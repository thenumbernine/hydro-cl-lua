local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local HLL = Flux:subclass()

HLL.name = 'hll'
HLL.solverCodeFile = 'hydro/flux/hll.cl'

--HLL.hllCalcWaveMethod = 'Davis direct'
HLL.hllCalcWaveMethod = 'Davis direct bounded'

function HLL:init(args)
	HLL.super.init(self, args)

	self.hllCalcWaveMethod = args.hllCalcWaveMethod
end

return HLL
