local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinFiniteDifferenceSolver = require 'solver.einstein-fd' 
local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'
BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

function BSSNOKFiniteDifferenceSolver:createDisplayComponents()
	-- skip EinsteinFiniteDifferenceSolver create components
	EinsteinFiniteDifferenceSolver.super.createDisplayComponents(self)
	
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted',
		code = [[
		sym3 gammaBar_LL = sym3_rescaleFromCoord_ll(calc_gammaBar_ll(U, x), x);
		*value = real3_weightedLen(*value_real3, gammaBar_LL);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		norm = 'tr weighted',
		code = [[
		sym3 gamma_UU = sym3_rescaleFromCoord_uu(calc_gamma_uu(U, x), x);
		*value = sym3_dot(*value_sym3, gamma_UU);]],
	})
end

return BSSNOKFiniteDifferenceSolver
