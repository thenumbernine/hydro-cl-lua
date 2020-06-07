local class = require 'ext.class'
local HLL = require 'hydro.solver.hll'


local EulerHLLC = class(HLL)

EulerHLLC.solverCodeFile = 'hydro/solver/euler-hllc.cl'
EulerHLLC.name = 'EulerHLLC'
EulerHLLC.eqnName = 'euler'

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
