-- based on 2014 Oliveira et al - Ergoregion instability- The hydrodynamic vortex
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'

local common = require 'common'
local xNames = common.xNames

local WaveFDEqn = class(Equation)
WaveFDEqn.name = 'wave-fd'

WaveFDEqn.hasCalcDTCode = true			-- don't use default
WaveFDEqn.hasEigenCode = true			-- don't use default
WaveFDEqn.hasFluxFromConsCode = true	-- don't use default
WaveFDEqn.useSourceTerm = true

WaveFDEqn.consVars = table{
	-- typically phi
	{name='psi', type='cplx', units='kg/(m^3)'},
	-- typically pi
	{name='zeta', type='cplx', units='kg/(m^3*s)'},
}

function WaveFDEqn:createBoundaryOptions()
	local BoundaryFixed = class(self.solver.BoundaryFixed)
	BoundaryFixed.fixedCode = '('..self.cons_t..'){.psi=cplx_zero, .zeta=cplx_zero}'
	self.solver:addBoundaryOption(BoundaryFixed)
end

function WaveFDEqn:createInitState()
	WaveFDEqn.super.createInitState(self)
	self:addGuiVars{
		{name='m', value=2},	-- mode of fourier transform of wave eqn
		{name='C', value=.5},
	}
end

-- the default display-all is broken since i switched to the pick-component option
WaveFDEqn.predefinedDisplayVars = {
	'U psi re',
	'U psi im',
	'U psi abs',
	--'U psi arg',
	'U zeta re',
	'U zeta im',
	'U zeta abs',
	--'U zeta arg',
}

WaveFDEqn.initStates = require 'init.nls'

WaveFDEqn.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	
	real r = fabs(x.x);
	cplx q = cplx_zero;
	<?=code?>
	UBuf[index] = (<?=eqn.cons_t?>){
		.psi = q,
		.zeta = cplx_zero,
	};
}
]]

WaveFDEqn.solverCodeFile = 'eqn/wave-fd.cl'

-- SolverBase adds eqn:getEigenTypeCode(), and eqn/eqn.lua provides a default, so this overrides it 
function WaveFDEqn:getEigenTypeCode() end

function WaveFDEqn:getDisplayVars()
	local vars = WaveFDEqn.super.getDisplayVars(self)
	vars:append{
		{name='J0', code='value.vreal = BESSJ0(x.x) * cos(t);', units='kg/(m^3)'},
	}
	return vars
end

return WaveFDEqn
