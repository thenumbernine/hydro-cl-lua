-- check out my 'wave equation hyperbolic form' worksheet
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local Wave = class(Equation)
Wave.name = 'wave'

-- depending on which you use, the source terms change
--Wave.weightFluxByGridVolume = true
Wave.weightFluxByGridVolume = false

Wave.useSourceTerm = true

-- TODO count this from consVars?
Wave.numStates = 4

Wave.mirrorVars = {{'phi_i.x'}, {'phi_i.y'}, {'phi_i.z'}}

Wave.hasEigenCode = true
Wave.hasFluxFromConsCode = true
Wave.roeUseFluxFromCons = true

Wave.initStates = require 'init.euler'	 -- use rho as our initial condition

-- Sedov initial conditions uses pressure
Wave.usePressure = true
-- gaussian, etc use density
--Wave.usePressure = false


if Wave.usePressure then
	Wave.consVars = table{
		{name='phi_t', type='real', units='kg/(m*s^3)'},
		{name='phi_i', type='real3', units='kg/(m^2*s^2)'},
	}
else
	Wave.consVars = table{
		{name='phi_t', type='real', units='kg/(m^3*s)'},
		{name='phi_i', type='real3', units='kg/(m^4)'},
	}
end

function Wave:createInitState()
	Wave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},
	}
end

Wave.initStateCode = [[
<?
local common = require 'common'()
local xNames = common.xNames
?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;
	
	<?=code?>

	UBuf[index] = (<?=eqn.cons_t?>){
<? if eqn.usePressure then
?>		.phi_t = P,
<? else		
?>		.phi_t = rho,
<? end		
?>		.phi_i = cartesianToCoord(v, x),
	};
}
]]

Wave.solverCodeFile = 'eqn/wave.cl'

Wave.predefinedDisplayVars = {
	'U phi_t (kg/(m s^3))',
	'U phi_i (kg/(m^2 s^2))',
	'U phi_i x (kg/(m^2 s^2))',
	'U phi_i y (kg/(m^2 s^2))',
	'U phi_i z (kg/(m^2 s^2))',
	'U phi_i mag (kg/(m^2 s^2))',
}

Wave.eigenVars = {
	{name='unused', type='real'},
}

function Wave:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real c_sqrt_gU = solver->wavespeed / unit_m_per_s * coord_sqrt_g_uu<?=side..side?>(<?=x?>);
]], {
		side = side,
		x = x,
	})
end

function Wave:consWaveCodePrefix(side, U, x, W)
	return self:eigenWaveCodePrefix(side, nil, x)
end

function Wave:consWaveCode(side, U, x, waveIndex)
	if waveIndex == 0 then
		return '-c_sqrt_gU' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return '0'
	elseif waveIndex == 3 then
		return 'c_sqrt_gU' 
	end
	error'got a bad waveIndex'
end

Wave.eigenWaveCode = Wave.consWaveCode

return Wave
