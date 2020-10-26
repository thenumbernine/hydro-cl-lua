local class = require 'ext.class'
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local Roe = class(Flux)
Roe.name = 'roe'
Roe.solverCodeFile = 'hydro/flux/roe.cl'
Roe.usesFluxLimiter = true

return Roe
