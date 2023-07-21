--[[
here's a copy of euler (which is based on conserved quantities)
except that this is using primitives as state variables
and that means the flux-jacobian is now the acoustic matrix + normal vel component (times I)
and that means ... i don't think the d/dW (A + I v) * dW is now integratable into any sort of flux vector ...
so now HLL won't work with this
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local Euler = class(Equation)
Euler.name = 'euler_prim'

Euler.initConds = require 'hydro.init.euler':getList()

function Euler:init(args)
	
	self.consVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m*s)', variance='u'},	-- contravariant
		{name='P', type='real', units='kg/(m*s^2)'},	-- per volume.  energy units are kg m^2 / s^2, but energy per volume is kg / (m s^2)
	}

	Euler.super.init(self, args)
end


function Euler:createInitState()
	Euler.super.createInitState(self)
	self:addGuiVars{	
		{name='heatCapacityRatio', value=7/5},				-- unitless
	}
end

-- don't use default
-- in fact, idk / there is no flux vector solvable from the (A + I v) dW/dx system ...
function Euler:initCodeModule_fluxFromCons() end

function Euler:initCodeModule_calcDTCell()
	-- ugly hack to insert more dependent modules
	local path = require 'ext.path'
	self.solver.modules:addFromMarkup(self:template(path'hydro/eqn/cl/calcDT.cl':read()))
end

Euler.solverCodeFile = 'hydro/eqn/euler_prim.clcpp'

-- [=[
Euler.predefinedDisplayVars = {
-- [[
	'U rho',
	'U v',
	'U P',
--]]
--[[
	-- now that I've switched to components, I can't display these
	-- TODO ...unless I allow for multiple displays of the same displayVar...
	-- TODO TODO after building kernels, duplicate all display vars for each component type
	-- then things will look just like they did before, except with less kernels
	'U v x',
	'U v y',
	'U v z',
--]]
	'U div v',
	'U curl v',
}
--]=]

function Euler:getDisplayVars()
	local vars = Euler.super.getDisplayVars(self)
	vars:append{
		-- TODO should the default display generated of variables be in solver units or SI units?
		-- should SI unit displays be auto generated as well?
		{name='v', code='value.vreal3 = U.v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = U.P;', units='kg/(m*s^2)'},
		{name='eInt', code=self:template'value.vreal = Eqn::calc_eInt(solver, U);', units='m^2/s^2'},
		{name='eKin', code=self:template'value.vreal = Eqn::calc_eKin(U, x);', units='m^2/s^2'},
		{name='eTotal', code=self:template'value.vreal = Eqn::calc_ETotal(solver, U, x) / U.rho;', units='m^2/s^2'},
		{name='EInt', code=self:template'value.vreal = Eqn::calc_EInt(solver, U);', units='kg/(m*s^2)'},
		{name='EKin', code=self:template'value.vreal = Eqn::calc_EKin(U, x);', units='kg/(m*s^2)'},
		{name='S', code='value.vreal = U.P / pow(U.rho, (real)solver.heatCapacityRatio);'},
		{name='H', code=self:template'value.vreal = Eqn::calc_H(solver, U.P);', units='kg/(m*s^2)'},
		{name='h', code=self:template'value.vreal = Eqn::calc_h(solver, U.rho, U.P);', units='m^2/s^2'},
		{name='HTotal', code=self:template[[
real const ETotal = Eqn::calc_ETotal(solver, U, x);
value.vreal = Eqn::calc_HTotal(U.P, ETotal);
]], units='kg/(m*s^2)'},
		{name='hTotal', code=self:template[[
real const ETotal = Eqn::calc_ETotal(solver, U, x);
value.vreal = Eqn::calc_hTotal(U.rho, U.P, ETotal);
]], units='m^2/s^2'},
		{name='speed of sound', code=self:template'value.vreal = Eqn::calc_Cs(solver, U);', units='m/s'},
		{name='Mach number', code=self:template'value.vreal = coordLen(U.v, x) / Eqn::calc_Cs(solver, U);'},
		{name='temperature', code=self:template'value.vreal = Eqn::calc_T(solver, U);', units='K'},
	}
	
	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		units = '1/s',
	} or nil)

	-- special for 1d_state_line
	vars:insert{
		name = 'state line',
		type = 'real3',
		units = '1',
		code = 'value.vreal3 = real3(U.rho, coordLen(U.v, x) * sign(U.v.x), U.P);',
	}

	return vars
end

Euler.eigenVars = table{
	-- Roe-averaged vars
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	{name='hTotal', type='real', units='m^2/s^2'},
	-- derived vars
	{name='Cs', type='real', units='m/s'},
	{name='vSq', type='real', units='m^2/s^2'},
	{name='vL', type='real3', units='m/s'},
}

return Euler

--[[
scratch work:

P V = n R T = m R_spec T
ideal gas EOS: P V_m = R T
T = T_C + 273.15' C
R = gas constant = 8.314459848 J / (mol K) 
R_spec = 287.058 J / (kg K) for dry air 
R_spec = R / M = k_B / m
M = molar mass
k_B = Boltzmann constant = 1.3806485279e-23 J / K
m = mass within volume V	<-> ρ = m / V
R_spec = k_B ρ V

caloric perfect:
P = (gamma - 1) ρ e_int
gamma = C_p / C_v = adiabatic index / ratio of specific heats
e_int = C_v T = internal energy per unit mass / specific internal energy
C_v = specific heat at constant volume
C_p = specific heat at constant pressure
P = ρ R_spec T
e_int = C_v T
P = ρ (C_p - C_v) e_int / C_v = ρ (C_p / C_v - 1) e_int

0 C = 273.15 K

using some real-world numbers ...
P = (gamma - 1) ρ e_int
101325 kg / (m s^2) = (C_p - C_v) / C_v (1.2754 kg / m^3) C_v T 
101325 kg / (m s^2) = (1006 - 717.1) J / (kg K) (1.2754 kg / m^3) T 
T = 274.99364522457 K
T = 1.8436452245715 C ... should be
--]]
