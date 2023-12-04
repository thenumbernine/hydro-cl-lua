local table = require 'ext.table'
local EinsteinFiniteDifferenceSolver = require 'hydro.solver.einstein-fd' 

local BSSNOKFiniteDifferenceSolver = EinsteinFiniteDifferenceSolver:subclass()
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'

-- split into bssnok-fd-num and bssnok-fd-sym
--BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

-- right now certain constraints like mirror are designed to work with numGhost=2
-- use numGhost=3 for compat with SENR
-- 3 = 4th order
BSSNOKFiniteDifferenceSolver.numGhost = 3

-- for certain hydro/eqn/bssnok-fd calculations, dt is based on grid only and no state vars
-- so we only need to calculate it once
function BSSNOKFiniteDifferenceSolver:calcDT()
	local dt = BSSNOKFiniteDifferenceSolver.super.calcDT(self)
	if not self.useFixedDT then
		if self.eqn.cflMethod == '2013 Baumgarte et al, eqn 32' 
		or self.eqn.cflMethod == '2017 Ruchlin et al, eqn 53'
		then
			-- fixedDT is already updated to the latest dt by SolverBase.calcDT
			self.useFixedDT = true
		end
	end
	return dt
end

return BSSNOKFiniteDifferenceSolver
