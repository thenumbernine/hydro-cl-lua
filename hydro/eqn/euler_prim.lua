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

-- true by default, disable since there is no flux vector of the primitive-state-variable form
Euler.roeUseFluxFromCons = false

Euler.initConds = require 'hydro.init.euler':getList()

function Euler:init(args)
	
	self.consVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m*s)', variance='u'},	-- contravariant
		{name='P', type='real', units='kg/(m*s^2)'},	-- per volume.  energy units are kg m^2 / s^2, but energy per volume is kg / (m s^2)
	}

	Euler.super.init(self, args)
end


function Euler:getSymbolFields()
	return table(Euler.super.getSymbolFields(self)):append{
		'calc_H',
		'calc_h',
		'calc_HTotal',
		'calc_hTotal',
		'calc_eKin',
		'calc_EKin',
		'calc_eInt',
		'calc_EInt',
		'calc_ETotal',
		'calc_Cs',
		'calc_T',
	}
end

function Euler:createInitState()
	Euler.super.createInitState(self)
	self:addGuiVars{	
		{name='heatCapacityRatio', value=7/5},				-- unitless
	}
end

-- don't use default
-- in fact, idk / there is no flux vector solvable from the (A + I v) dW/dx system ...
function Euler:initCodeModule_fluxFromCons() 
end

function Euler:initCodeModule_calcDTCell()
	-- ugly hack to insert more dependent modules
	local path = require 'ext.path'
	self.solver.modules:addFromMarkup(self:template(path'hydro/eqn/cl/calcDT.cl':read()
	..self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
]]))
end

Euler.solverCodeFile = 'hydro/eqn/euler_prim.cl'

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
		{name='v', code='value.vreal3 = U->v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = U->P;', units='kg/(m*s^2)'},
		{name='eInt', code=self:template'value.vreal = <?=calc_eInt?>(solver, U);', units='m^2/s^2'},
		{name='eKin', code=self:template'value.vreal = <?=calc_eKin?>(U, x);', units='m^2/s^2'},
		{name='eTotal', code=self:template'value.vreal = <?=calc_ETotal?>(solver, U, x) / U->rho;', units='m^2/s^2'},
		{name='EInt', code=self:template'value.vreal = <?=calc_EInt?>(solver, U);', units='kg/(m*s^2)'},
		{name='EKin', code=self:template'value.vreal = <?=calc_EKin?>(U, x);', units='kg/(m*s^2)'},
		{name='S', code='value.vreal = U->P / pow(U->rho, (real)solver->heatCapacityRatio);'},
		{name='H', code=self:template'value.vreal = <?=calc_H?>(solver, U->P);', units='kg/(m*s^2)'},
		{name='h', code=self:template'value.vreal = <?=calc_h?>(solver, U->rho, U->P);', units='m^2/s^2'},
		{name='HTotal', code=self:template[[
real const ETotal = <?=calc_ETotal?>(solver, U, x);
value.vreal = <?=calc_HTotal?>(U->P, ETotal);
]], units='kg/(m*s^2)'},
		{name='hTotal', code=self:template[[
real const ETotal = <?=calc_ETotal?>(solver, U, x);
value.vreal = <?=calc_hTotal?>(U->rho, U->P, ETotal);
]], units='m^2/s^2'},
		{name='speed of sound', code=self:template'value.vreal = <?=calc_Cs?>(solver, U);', units='m/s'},
		{name='Mach number', code=self:template'value.vreal = coordLen(U->v, x) / <?=calc_Cs?>(solver, U);'},
		{name='temperature', code=self:template'value.vreal = <?=calc_T?>(solver, U);', units='K'},
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
		code = 'value.vreal3 = real3(U->rho, coordLen(U->v, x) * sign(U->v.x), U->P);',
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

function Euler:eigenWaveCodePrefix(args)
	return self:template([[
real const <?=eqn.symbolPrefix?>Cs_nLen = normal_len(<?=n?>) * (<?=eig?>)->Cs;
real const <?=eqn.symbolPrefix?>v_n = normal_vecDotN1(<?=n?>, (<?=eig?>)->v);
]],	args)
end

function Euler:eigenWaveCode(args)
	if args.waveIndex == 0 then
		return self:template'<?=eqn.symbolPrefix?>v_n - <?=eqn.symbolPrefix?>Cs_nLen'
	elseif args.waveIndex >= 1 and args.waveIndex <= 3 then
		return self:template'<?=eqn.symbolPrefix?>v_n'
	elseif args.waveIndex == 4 then
		return self:template'<?=eqn.symbolPrefix?>v_n + <?=eqn.symbolPrefix?>Cs_nLen'
	end
	error'got a bad waveIndex'
end

-- W is an extra param specific to Euler's calcDT in this case
-- but then I just explicitly wrote out the calcDT, so the extra parameters just aren't used anymore.
function Euler:consWaveCodePrefix(args)
	return self:template([[
real <?=eqn.symbolPrefix?>Cs_nLen = <?=calc_Cs?>(solver, <?=U?>) * normal_len(<?=n?>);
real const <?=eqn.symbolPrefix?>v_n = normal_vecDotN1(<?=n?>, (<?=U?>)->v);
]], args)
end

-- as long as U or eig isn't used, we can use this for both implementations
Euler.consWaveCode = Euler.eigenWaveCode

--Euler.eigenWaveCodeMinMax uses default
--Euler.consWaveCodeMinMax uses default

-- alright, thanks to eqn/euler you now have to always call this before the for-loop along sides of cells in calcDT
-- ok so this goes before consMin/MaxWaveCode
function Euler:consWaveCodeMinMaxAllSidesPrefix(args)
	return self:template([[
real <?=eqn.symbolPrefix?>Cs = <?=calc_Cs?>(solver, <?=U?>);\
]],	args)
end

--[[
set resultMin or resultMax to store the resulting min / max in either
you have to call 'consWaveCodeMinMaxAllSidesPrefix' before you call this
but you don't have to call 'consWaveCodePrefix'
should this function create the vars or assign the vars?
I'll go with create so it can create them const.

so in calcDT this is used with the allsides prefix.
but in hll flux this is used with one specific side.
--]]
function Euler:consWaveCodeMinMaxAllSides(args)
	return self:template([[
real const <?=eqn.symbolPrefix?>Cs_nLen = <?=eqn.symbolPrefix?>Cs * normal_len(<?=n?>);
real const <?=eqn.symbolPrefix?>v_n = normal_vecDotN1(<?=n?>, (<?=U?>)->v);
<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	eqn.symbolPrefix..'v_n - '..eqn.symbolPrefix..'Cs_nLen',
	eqn.symbolPrefix..'v_n + '..eqn.symbolPrefix..'Cs_nLen'
)?>
]], args)
end

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
