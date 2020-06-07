local class = require 'ext.class'
local EinsteinFiniteDifferenceSolver = require 'hydro.solver.einstein-fd' 

local Z4cFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
Z4cFiniteDifferenceSolver.name = 'Z4c_FiniteDifference'
Z4cFiniteDifferenceSolver.eqnName = 'z4c-fd'
return Z4cFiniteDifferenceSolver
