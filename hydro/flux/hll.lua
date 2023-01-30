local class = require 'ext.class'
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local HLL = class(Flux)

HLL.name = 'hll'
HLL.solverCodeFile = 'hydro/flux/hll.clcpp'

--HLL.hllCalcWaveMethod = 'Davis direct'
HLL.hllCalcWaveMethod = 'Davis direct bounded'

function HLL:init(args)
	HLL.super.init(self, args)

	self.hllCalcWaveMethod = args.hllCalcWaveMethod
end

return HLL
