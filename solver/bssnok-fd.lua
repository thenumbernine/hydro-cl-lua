local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinFiniteDifferenceSolver = require 'solver.einstein-fd' 
local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'

-- split into bssnok-fd-num and bssnok-fd-sym
--BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

-- right now certain constraints like mirror are designed to work with numGhost=2
BSSNOKFiniteDifferenceSolver.numGhost = 2

function BSSNOKFiniteDifferenceSolver:createDisplayComponents()
	-- skip EinsteinFiniteDifferenceSolver create components
	EinsteinFiniteDifferenceSolver.super.createDisplayComponents(self)
	
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_IJ',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	*value = real3_weightedLen(*value_real3, gammaBar_LL);]],
	})
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_ij',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	*value = real3_weightedLen(*value_real3, gammaBar_ll);]],
	})
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma^ij',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gamma_uu = calc_gamma_uu(U, x);
	*value = real3_weightedLen(*value_real3, gamma_uu);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^IJ',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gamma_UU = sym3_rescaleFromCoord_uu(calc_gamma_uu(U, x), x);
	*value = sym3_dot(*value_sym3, gamma_UU);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gammaBar^IJ',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_UU = calc_gammaBar_UU(U, x);
	*value = sym3_dot(*value_sym3, gammaBar_UU);]],
	})
end

return BSSNOKFiniteDifferenceSolver
