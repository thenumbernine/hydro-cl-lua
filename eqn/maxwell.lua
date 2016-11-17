local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local GuiFloat = require 'guivar.float'
local clnumber = require 'clnumber'
local processcl = require 'processcl'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.numStates = 6
Maxwell.consVars = {'epsEx', 'epsEx', 'epsEz', 'Bx', 'By', 'Bz'}
Maxwell.mirrorVars = {{'epsE.x', 'B.x'}, {'epsE.y', 'B.y'}, {'epsE.z', 'B.z'}}

Maxwell.hasEigenCode = true
Maxwell.useSourceTerm = true

Maxwell.initStates= {{name='default'}}

Maxwell.guiVars = table{
	GuiFloat{name='eps0', value=1},	-- permittivity
	GuiFloat{name='mu0', value=1},	-- permeability
	GuiFloat{name='sigma', value=1},-- conductivity
}

function Maxwell:getTypeCode()
	return [[
typedef struct {
	real3 epsE;
	real3 B;
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
	return table{
		processcl(
		[[
__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
<? if solver.dim > 1 then ?>
		&& x.y < mids.y
<? end
if solver.dim > 2 then ?>
		&& x.z < mids.z
<? end ?>
	;
	__global cons_t* U = UBuf + index;
	U->epsE = real3_scale(_real3(1,0,0), eps0);
	U->B = _real3(0, 1, lhs ? 1 : -1);
}
]], {solver=solver}),
	}:concat'\n'
end

function Maxwell:getSolverCode(solver)
	return processcl(file['eqn/maxwell.cl'], {solver=solver})
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
		return {['div_'..var] = processcl([[
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
