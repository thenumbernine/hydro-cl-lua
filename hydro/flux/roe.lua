local class = require 'ext.class'
local Flux = require 'hydro.flux.flux'

local Roe = class(Flux)
Roe.name = 'roe'
Roe.solverCodeFile = 'hydro/flux/roe.cl'

return Roe
