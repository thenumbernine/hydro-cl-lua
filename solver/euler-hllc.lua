local class = require 'ext.class'
local HLL = require 'solver.hll'

local EulerHLLC = class(HLL)

EulerHLLC.solverCodeFile = 'solver/euler-hllc.cl'
EulerHLLC.name = 'EulerHLLC'
EulerHLLC.eqnName = 'euler'

return EulerHLLC
