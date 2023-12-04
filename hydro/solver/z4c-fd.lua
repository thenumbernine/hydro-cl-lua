local EinsteinFiniteDifferenceSolver = require 'hydro.solver.einstein-fd' 

local Z4cFiniteDifferenceSolver = EinsteinFiniteDifferenceSolver:subclass()
Z4cFiniteDifferenceSolver.name = 'Z4c_FiniteDifference'
Z4cFiniteDifferenceSolver.eqnName = 'z4c-fd'
return Z4cFiniteDifferenceSolver
