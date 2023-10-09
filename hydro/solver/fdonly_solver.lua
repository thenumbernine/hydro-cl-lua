-- finite-difference solver, driven fully by rhs in eqn's 'addSource'
-- TODO tbf if I want an eqn to be compatible with both this fd and with fv, then it should have a *new* routine for calculating the rhs which is turned into the flux
-- or idk some way of setting part of the eqn as optionally hyperbolic-conservation-law *or* finite-difference,
--  and other parts solely rhs source terms.
-- TODO name this 'fdsolver' and name 'fdsolver' => 'fd_from_fv_solver'
--
-- honestly, this is just GridSolver ...
-- ...until I straighten out how to modularize the equations
local GridSolver = require 'hydro.solver.gridsolver'

-- TODO TODO this should be fv_to_fd_solver because it takes the flux and turns it into finite-difference
local FiniteDifferenceSolver = GridSolver:subclass()
FiniteDifferenceSolver.name = 'FiniteDifference'
FiniteDifferenceSolver.numGhost = 2
return FiniteDifferenceSolver 
