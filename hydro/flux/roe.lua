local class = require 'ext.class'
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local Roe = class(Flux)
Roe.name = 'roe'
Roe.solverCodeFile = 'hydro/flux/roe.cl'

function Roe:getModuleDepends_calcFlux()
	local depends = table(Roe.super.getModuleDepends_calcFlux(self))
	depends:insert'fluxLimiter'
	depends:insert'fluxFromCons'
	if self.solver.eqn.roeUseFluxFromCons then
		depends:insert'fluxFromCons'
	end
	if not require 'hydro.solver.meshsolver'.is(solver) then
		depends:insert'cell_area#'
	end
	return depends
end

return Roe
