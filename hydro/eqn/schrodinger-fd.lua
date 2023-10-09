--[[
Schrodinger equation
because it's parabolic, no hyperbolic conservation law.
For that use eqn/wave with scalar=cplx
or use one of the eqn/maxwell's

i ℏ Ψ_,t = - ℏ^2/(2 m) Ψ_;j^j + V Ψ
[ℏ] = kg m^2 / s
[m] = kg
[V] = kg m^2 / s^2 = energy
[ℏ^2 / (2m)] = kg m^4 / s^2
[ℏ^2 / (2m)] / m^2 = kg m^2 / s^2 = energy
[ℏ] / s = kg m^2 / s^2 = energy
[ℏ]/s [Ψ] = [ℏ^2 / (2m)] [Ψ] / m^2 = [V] [Ψ]
[Ψ] is arbitrary

Ψ_,t = i ℏ/(2 m) Ψ_;j^j - i V / ℏ Ψ
--]]
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'

local SchrodingerFDEqn = Equation:subclass()
SchrodingerFDEqn.name = 'schrodinger_fd'

SchrodingerFDEqn.consVars = table{
	{name='psi', type='cplx', units='kg*m^2/(s^2)'},
}

function SchrodingerFDEqn:createInitState()
	SchrodingerFDEqn.super.createInitState(self)
	self:addGuiVars{
		{name='hBar', value=1, units='kg*m^2/s'},	-- TODO I could just use the real constant and let units determine this' magnitude ...
		{name='m', value=1, units='kg'},			-- mass ... of what? of the field? but the field is *everywhere* ... mass of *everything* ?
	}
end

-- TODO what kind of init cond for this?
SchrodingerFDEqn.initConds = require 'hydro.init.nls':getList()

SchrodingerFDEqn.solverCodeFile = 'hydro/eqn/schrodinger-fd.cl'

-- don't use default
function SchrodingerFDEqn:initCodeModule_calcDT() end
function SchrodingerFDEqn:initCodeModule_calcDTCell() end

-- the default display-all is broken since i switched to the pick-component option
SchrodingerFDEqn.predefinedDisplayVars = {
	'U psi re',
	'U psi im',
	'U psi abs',
	--'U psi arg',
}

return SchrodingerFDEqn
