--[[
Marquina flux from Appendix A of 1997 Marti, Muller, Font, Ibanez, Marquina "Morphology and Dynamics of Relativistic Jets"
--]]
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local Marquina = Flux:subclass()

Marquina.name = 'hll'
Marquina.solverCodeFile = 'hydro/flux/marquina.cl'

return Marquina
