local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinFiniteDifferenceSolver = require 'hydro.solver.einstein-fd' 

local BSSNOKFiniteDifferenceSolver = class(EinsteinFiniteDifferenceSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'

-- split into bssnok-fd-num and bssnok-fd-sym
--BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

-- right now certain constraints like mirror are designed to work with numGhost=2
-- use numGhost=3 for compat with SENR
-- 3 = 4th order
BSSNOKFiniteDifferenceSolver.numGhost = 3

function BSSNOKFiniteDifferenceSolver:createDisplayComponents()
	BSSNOKFiniteDifferenceSolver.super.createDisplayComponents(self)
	
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_IJ',
		code = self.eqn:template[[
	int index = INDEXV(i);
	real3 x = cellBuf[index].pos;
	const global <?=cons_t?>* U = buf + index;
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	value->vreal = real3_weightedLen(value->vreal3, gammaBar_LL);]],
	})
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_ij',
		code = self.eqn:template[[
	int index = INDEXV(i);
	real3 x = cellBuf[index].pos;
	const global <?=cons_t?>* U = buf + index;
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	value->vreal = real3_weightedLen(value->vreal3, gammaBar_ll);]],
	})
--[=[ now in eqn.einstein
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma^ij',
		code = self.eqn:template[[
	int index = INDEXV(i);
	real3 x = cellBuf[index].pos;
	const global <?=cons_t?>* U = buf + index;
	sym3 gamma_uu = <?=calc_gamma_uu?>(U, x);
	value->vreal = real3_weightedLen(value->vreal3, gamma_uu);]],
	})
--]=]	
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^IJ',
		code = self.eqn:template[[
	int index = INDEXV(i);
	real3 x = cellBuf[index].pos;
	const global <?=cons_t?>* U = buf + index;
	sym3 gamma_UU = sym3_rescaleFromCoord_uu(<?=calc_gamma_uu?>(U, x), x);
	value->vreal = sym3_dot(value->vsym3, gamma_UU);]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gammaBar^IJ',
		code = self.eqn:template[[
	int index = INDEXV(i);
	real3 x = cellBuf[index].pos;
	const global <?=cons_t?>* U = buf + index;
	sym3 gammaBar_UU = calc_gammaBar_UU(U, x);
	value->vreal = sym3_dot(value->vsym3, gammaBar_UU);]],
	})
end

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
