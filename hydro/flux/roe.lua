local class = require 'ext.class'
local table = require 'ext.table'
local Flux = require 'hydro.flux.flux'

local Roe = class(Flux)
Roe.name = 'roe'
Roe.solverCodeFile = 'hydro/flux/roe.clcpp'
Roe.usesFluxLimiter = true

-- specifically for Euler, but I guess anyone can enable this
Roe.useEntropyFluxFix = false
Roe.entropyFluxFixLimitingFactor = 0.2

function Roe:init(args)
	Roe.super.init(self, args)
	self.useEntropyFluxFix = args.useEntropyFluxFix
	self.entropyFluxFixLimitingFactor = args.entropyFluxFixLimitingFactor
end

return Roe
