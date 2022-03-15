-- based on 2014 Oliveira et al - Ergoregion instability- The hydrodynamic vortex
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'

local WaveFDEqn = class(Equation)
WaveFDEqn.name = 'wave_fd'

WaveFDEqn.consVars = table{
	-- typically phi
	{name='psi', type='cplx', units='kg/(m^3)'},
	-- typically pi
	{name='zeta', type='cplx', units='kg/(m^3*s)'},
}

function WaveFDEqn:createBoundaryOptions()
	local BoundaryFixed = class(self.solver.BoundaryFixed)
	BoundaryFixed.fixedCode = '('..self.symbols.cons_t..'){.psi=cplx_zero, .zeta=cplx_zero}'
	self.solver:addBoundaryOption(BoundaryFixed)
end

function WaveFDEqn:createInitState()
	WaveFDEqn.super.createInitState(self)
	self:addGuiVars{
		{name='m', value=2},	-- mode of fourier transform of wave eqn
		{name='C', value=.5},
	}
end

WaveFDEqn.initConds = require 'hydro.init.nls':getList()

WaveFDEqn.solverCodeFile = 'hydro/eqn/wave-fd.cl'

-- don't use default
function WaveFDEqn:initCodeModule_calcDT() end
function WaveFDEqn:initCodeModule_calcDTCell() end

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

function WaveFDEqn:getDisplayVars()
	local vars = WaveFDEqn.super.getDisplayVars(self)
	vars:append{
		{name='J0', code=[[
//// MODULE_DEPENDS: Bessel
value.vreal = BESSJ0(x.x) * cos(t);
]], units='kg/(m^3)'},
	}
	return vars
end

return WaveFDEqn
