local class = require 'ext.class'
local EinsteinFiniteDifferenceSolver = require 'solver.einstein-fd' 
local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'
BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

function BSSNOKFiniteDifferenceSolver:getDisplayInfosForType()
	
	-- skip EinsteinFiniteDifferenceSolver
	local t = EinsteinFiniteDifferenceSolver.super.getDisplayInfosForType(self)

	-- hmm, only works with L ... so it only applies to L ...
	table.insert(t.real3, {
		name = ' norm weighted',
		code = '	*value = real3_weightedLen(*value_real3, calc_gammaBar_LL(U, x));',
	})

	-- hmm, how to do the weighting stuff with gammaBar_ll ... 
	-- also, how to determine which metric to raise by ... gamma vs gammaBar
	table.insert(t.sym3, {
		name = ' tr weighted',
		code = '	*value = sym3_dot(*value_sym3, calc_gamma_UU(U, x));',
	})

	return t
end


return BSSNOKFiniteDifferenceSolver
