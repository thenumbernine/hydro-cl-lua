local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local GuiFloat = require 'guivar.float'
local clnumber = require 'clnumber'
local template = require 'template'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.numStates = 6
Maxwell.consVars = {'epsEx', 'epsEx', 'epsEz', 'Bx', 'By', 'Bz'}
Maxwell.mirrorVars = {{'epsE.x', 'B.x'}, {'epsE.y', 'B.y'}, {'epsE.z', 'B.z'}}

Maxwell.hasEigenCode = true
Maxwell.useSourceTerm = true

Maxwell.initStates = require 'init.euler'

Maxwell.guiVars = table{
	GuiFloat{name='eps0', value=1},	-- permittivity
	GuiFloat{name='mu0', value=1},	-- permeability
	GuiFloat{name='sigma', value=1},-- conductivity
}

function Maxwell:getTypeCode()
	return [[
typedef union {
	real ptr[6];
	struct {
		real3 epsE;
		real3 B;
	};
} cons_t;
]]
end

function Maxwell:getCodePrefix()
	return table{
		Maxwell.super.getCodePrefix(self),
		'#define sqrt_eps0 '..clnumber(math.sqrt(self.guiVarsForName.eps0.value[0])),
		'#define sqrt_mu0 '..clnumber(math.sqrt(self.guiVarsForName.mu0.value[0])),
		[[
real ESq(cons_t U) { 
	return coordLenSq(U.epsE) / (eps0 * eps0);
}

real BSq(cons_t U) {
	return coordLenSq(U.B);
}
]]
	}:concat'\n'
end

function Maxwell:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local code = initState.init(solver)	
	return [[
kernel void initState(
	global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	global cons_t* U = UBuf + index;

	//used
	real3 E = _real3(0,0,0);
	real3 B = _real3(0,0,0);
	//throw-away
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	real ePot = 0;

]]..code..[[
	
	U->epsE = real3_scale(E, eps0);
	U->B = B;
}
]]
end

function Maxwell:getSolverCode(solver)
	return template(file['eqn/maxwell.cl'], {solver=solver})
end

function Maxwell:getDisplayVars(solver)
	return table{
		{Ex = 'value = U->epsE.x / eps0;'},
		{Ey = 'value = U->epsE.y / eps0;'},
		{Ez = 'value = U->epsE.z / eps0;'},
		{E = 'value = sqrt(ESq(*U));'},
		{Bx = 'value = U->B.x;'},
		{By = 'value = U->B.y;'},
		{Bz = 'value = U->B.z;'},
		{B = 'value = sqrt(BSq(*U));'},
		{energy = 'value = .5 * (coordLen(U->epsE) + coordLen(U->B) / mu0);'},
	}:append(table{'E','B'}:map(function(var,i)
		local field = assert( ({E='epsE', B='B'})[var] )
		return {['div '..var] = template([[
	value = 0;
	<? for j=0,solver.dim-1 do ?>{
		value += (U[stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- buf[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / grid_dx<?=j?>;
	}<? end ?>
	value *= .5;
	<? if field == 'epsE' then ?>
	value /= eps0;
	<? end ?>
]], {solver=solver, field=field})}
	end))
end

-- can it be zero sized?
function Maxwell:getEigenTypeCode(solver)
	return 'typedef struct { char mustbesomething; } eigen_t;'
end

function Maxwell:getEigenDisplayVars(solver)
	return {}
end

return Maxwell
