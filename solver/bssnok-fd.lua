local class = require 'ext.class'
local EinsteinFiniteDifferenceSolver = require 'solver.einstein-fd' 
local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'
BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'
return BSSNOKFiniteDifferenceSolver
