-- check out my 'wave equation hyperbolic form' worksheet
local ffi = require 'ffi'
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

Wave.hasEigenCode = true
Wave.hasFluxFromConsCode = true
Wave.roeUseFluxFromCons = true

Wave.initStates = require 'init.euler'	 -- use rho as our initial condition

-- Sedov initial conditions uses pressure
Wave.usePressure = true
-- gaussian, etc use density
--Wave.usePressure = false

function Wave:init(args)
	self.scalar = 'real'
	--self.scalar = 'cplx'
	
	self.vec3 = self.scalar..'3'
	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'
	
	self.numStates = 4 * self.numRealsInScalar

	if self.usePressure then
		self.consVars = table{
			{name='Pi', type=self.scalar, units='kg/(m*s^3)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^2*s^2)'},
		}
	else
		self.consVars = table{
			{name='Pi', type=self.scalar, units='kg/(m^3*s)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^4)'},
		}
	end

	local suffix = self.scalar == 'real' and '' or ' re'
	self.predefinedDisplayVars = {
		'U Pi'..suffix,
		'U Psi_l'..suffix,
		'U Psi_l x'..suffix,
		'U Psi_l y'..suffix,
		'U Psi_l z'..suffix,
		'U Psi_l mag',
	}

	Wave.super.init(self, args)
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

local scalar = eqn.scalar
local vec3 = eqn.vec3
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
?>		.Pi = <?=scalar?>_from_real(P),
<? else		
?>		.Pi = <?=scalar?>_from_real(rho),
<? end		
?>		.Psi_l = <?=vec3?>_from_real3(cartesianToCoord(v, x)),
	};
}
]]

Wave.solverCodeFile = 'eqn/wave.cl'

function Wave:getCommonFuncCode()
	return template([[
<?
local scalar = eqn.scalar
local vec3 = eqn.vec3
?>
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
real3 metric_beta_u(real3 x) { return real3_zero; }
real metric_K(real3 x) { return 0.; }

real metric_dt_alpha(real3 x) { return 0.; }
real3 metric_partial_alpha_l(real3 x) { return real3_zero; }
real3x3 metric_partial_beta_ul(real3 x) { return _real3x3(0,0,0,0,0,0,0,0,0); }

<?=scalar?> eqn_source(real3 x) { return <?=scalar?>_zero; }

]], {
		eqn = self,
	})
end

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
	waveIndex = math.floor(waveIndex / self.numRealsInScalar)
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
