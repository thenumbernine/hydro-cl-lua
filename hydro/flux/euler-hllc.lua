local class = require 'ext.class'
local table = require 'ext.table'
local HLL = require 'hydro.flux.hll'

local EulerHLLC = class(HLL)

EulerHLLC.name = 'euler-hllc'
EulerHLLC.solverCodeFile = 'hydro/flux/euler-hllc.clcpp'

--[[
args:
	hllcMethod = hllcMethod option
		options from 2012 Toro "The HLLC Riemann Solver"
		hllcMethod == 0 <=> eqns 38-39
		hllcMethod == 0 <=> 'variation 1' of the paper: eqns 40-41
		hllcMethod == 1 <=> 'variation 2' of the paper: eqns 42-44
-]]
function EulerHLLC:init(args)
	self.hllcMethod = args.hllcMethod or 2
	EulerHLLC.super.init(self, args)
end

return EulerHLLC
