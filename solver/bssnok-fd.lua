local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinFiniteDifferenceSolver = require 'solver.einstein-fd' 
local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'

-- split into bssnok-fd-num and bssnok-fd-sym
--BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

BSSNOKFiniteDifferenceSolver.numGhost = 4

function BSSNOKFiniteDifferenceSolver:createDisplayComponents()
	-- skip EinsteinFiniteDifferenceSolver create components
	EinsteinFiniteDifferenceSolver.super.createDisplayComponents(self)
	
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	*value = real3_weightedLen(*value_real3, gammaBar_LL);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gamma_UU = sym3_rescaleFromCoord_uu(calc_gamma_uu(U, x), x);
	*value = sym3_dot(*value_sym3, gamma_UU);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'trBar weighted',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_UU = calc_gammaBar_UU(U, x);
	*value = sym3_dot(*value_sym3, gammaBar_UU);]],
	})
end

return BSSNOKFiniteDifferenceSolver
