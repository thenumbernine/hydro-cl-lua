-- check out my 'wave equation hyperbolic form' worksheet
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local common = require 'common'
local xNames = common.xNames

local Wave = class(Equation)
Wave.name = 'wave'

-- depending on which you use, the source terms change
--Wave.weightFluxByGridVolume = true
Wave.weightFluxByGridVolume = false

Wave.useSourceTerm = true

-- TODO count this from consVars?
Wave.numStates = 4

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
local common = require 'common'
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

function Wave:getCommonFuncCode()
	return template([[
/*
background metric ADM decomposition
this assumes the ADM spatial metric gamma_ij is equal to the grid metric

for the wave equation d'Lambertian phi = 0
i.e. g^tt phi_,tt + 2 g^ti phi_;ti + g^ij phi_;ij = 0
ADM is defined such that 
 	g^tt = -alpha^-2
	g^ti = alpha^-2 beta^i 
	g^ij = gamma^ij - alpha^-2 beta^i beta^j
TODO make this configurable somehow
or make it modular enough to merge with BSSNOK
*/

real metric_alpha(real3 x) { return 1.; }
real metric_dalpha_t(real3 x) { return 0.; }
real3 metric_dalpha_l(real3 x) { return real3_zero; }

real3 metric_beta_u(real3 x) { return real3_zero; }
real3x3 metric_dbeta_ul(real3 x) { return _real3x3(0,0,0,0,0,0,0,0,0); }

real metric_K(real3 x) { return 0.; }

real eqn_source(real3 x) { return 0.; }

]], {
	})
end

Wave.predefinedDisplayVars = {
	'U phi_t',
	'U phi_i',
	'U phi_i x',
	'U phi_i y',
	'U phi_i z',
	'U phi_i mag',
}

Wave.eigenVars = {
	{name='unused', type='real'},
}

function Wave:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real wavespeed = solver->wavespeed / unit_m_per_s;
	real alpha_sqrt_gU = metric_alpha(<?=x?>) * coord_sqrt_g_uu<?=side..side?>(<?=x?>);
	real3 beta_u = metric_beta_u(<?=x?>);
]], {
		side = side,
		x = x,
	})
end

function Wave:consWaveCodePrefix(side, U, x, W)
	return self:eigenWaveCodePrefix(side, nil, x)
end

function Wave:consWaveCode(side, U, x, waveIndex)
	local xside = xNames[side+1]
	if waveIndex == 0 then
		return 'wavespeed * (-beta_u.'..xside..' - alpha_sqrt_gU)' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return 'wavespeed * -beta_u.'..xside
	elseif waveIndex == 3 then
		return 'wavespeed * (-beta_u.'..xside..' + alpha_sqrt_gU)' 
	end
	error'got a bad waveIndex'
end

Wave.eigenWaveCode = Wave.consWaveCode

return Wave
