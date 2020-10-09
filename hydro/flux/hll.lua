local class = require 'ext.class'
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local HLL = class(Flux)

HLL.name = 'hll'
HLL.solverCodeFile = 'hydro/flux/hll.cl'

--HLL.hllCalcWaveMethod = 'Davis direct'
HLL.hllCalcWaveMethod = 'Davis direct bounded'

function HLL:init(args)
	HLL.super.init(self, args)

	self.hllCalcWaveMethod = args.hllCalcWaveMethod
end

function HLL:getModuleDepends_calcFlux()
	local depends = table(HLL.super.getModuleDepends_calcFlux(self))
	:append{
		'fluxFromCons',
		'eigen_forInterface',
	}
	if not require 'hydro.solver.meshsolver'.is(solver) then
		depends:insert'cell_area#'
	end
	return depends
end

return HLL
