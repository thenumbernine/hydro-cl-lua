local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinFiniteDifferenceSolver = require 'hydro.solver.einstein-fd' 

local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'

-- split into bssnok-fd-num and bssnok-fd-sym
--BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

-- 3 = 4th order
BSSNOKFiniteDifferenceSolver.numGhost = 3

-- [=[ commenting out to try to reduce build time
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
	value->vreal = real3_weightedLen(value->vreal3, gammaBar_LL);]],
	})
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_ij',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	value->vreal = real3_weightedLen(value->vreal3, gammaBar_ll);]],
	})
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma^ij',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gamma_uu = calc_gamma_uu(U, x);
	value->vreal = real3_weightedLen(value->vreal3, gamma_uu);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^IJ',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gamma_UU = sym3_rescaleFromCoord_uu(calc_gamma_uu(U, x), x);
	value->vreal = sym3_dot(value->vsym3, gamma_UU);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gammaBar^IJ',
		code = [[
	int index = INDEXV(i);
	real3 x = cell_x(i);
	const global <?=eqn.cons_t?>* U = buf + index;
	sym3 gammaBar_UU = calc_gammaBar_UU(U, x);
	value->vreal = sym3_dot(value->vsym3, gammaBar_UU);]],
	})
end
--]=]

return BSSNOKFiniteDifferenceSolver
