-- check out my 'wave equation hyperbolic form' worksheet
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local Wave = class(Equation)
Wave.name = 'wave'

Wave.weightFluxByGridVolume = true

Wave.useSourceTerm = false
--Wave.useSourceTerm = true

Wave.numStates = 4

Wave.mirrorVars = {{'phi_i.x'}, {'phi_i.y'}, {'phi_i.z'}}

Wave.hasEigenCode = true
Wave.hasFluxFromConsCode = true
Wave.roeUseFluxFromCons = true

Wave.initStates = require 'init.euler'	 -- use rho as our initial condition

Wave.consVars = table{
	{phi_t = 'real'},
	{phi_i = 'real3'},
}

function Wave:createInitState()
	Wave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1},
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
		//hmm, this might deserve its own initial conditions ...
		// for the initialization of these variables:
		//.phi_t = rho,
		.phi_t = P,
		.phi_i = real3_zero,
	};
}
]]

Wave.solverCodeFile = 'eqn/wave.cl'

Wave.eigenVars = {{unused = 'real'}}

function Wave:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real c_sqrt_gU = solver->wavespeed * coord_sqrt_gU<?=side..side?>(<?=x?>);
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
