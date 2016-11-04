local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'
local GuiFloat = require 'guivar.float'
local clnumber = require 'clnumber'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.numStates = 6
Maxwell.consVars = {'epsE0', 'epsE0', 'epsE2', 'B0', 'B1', 'B2'}
Maxwell.mirrorVars = {{'epsE.s0', 'B.s0'}, {'epsE.s1', 'B.s1'}, {'epsE.s2', 'B.s2'}}
Maxwell.displayVars = {
	'E0', 'E1', 'E2', 'E',
	'B0', 'B1', 'B2', 'B',
	'energy',
}

Maxwell.useSourceTerm = true

Maxwell.initStateNames = {'default'}

Maxwell.guiVars = table{
	GuiFloat{name='eps0', value=1},	-- permittivity
	GuiFloat{name='mu0', value=1},	-- permeability
	GuiFloat{name='sigma', value=1},-- conductivity
}
Maxwell.guiVarsForName = Maxwell.guiVars:map(function(var) return var, var.name end)

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
		[[
__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = CELL_X(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.s0 < mids.s0
#if dim > 1
		&& x.s1 < mids.s1
#endif
#if dim > 2
		&& x.s2 < mids.s2
#endif
	;
	__global cons_t* U = UBuf + index;
	U->epsE = real3_scale(_real3(0,0,1), eps0);
	U->B = _real3(1, lhs ? 1 : -1, 0);
}
]],
	}:concat'\n'
end

function Maxwell:getSolverCode(solver)
	return table{ 
		'#include "eqn/maxwell.cl"',
	}:concat'\n'
end

function Maxwell:getEigenInfo()
	return {
		typeCode = 
			-- can it be zero sized?
			'typedef struct { char mustbesomething; } eigen_t;\n'..	
			'typedef struct { char mustbesomething; } fluxXform_t;',
		code = nil,
		displayVars = {},
	}
end

function Maxwell:getCalcDisplayVarCode()
	return [[
	switch (displayVar) {
	case display_U_E0: value = U->epsE.s0 / eps0; break;
	case display_U_E1: value = U->epsE.s1 / eps0; break;
	case display_U_E2: value = U->epsE.s2 / eps0; break;
	case display_U_E: value = sqrt(ESq(*U)); break;
	case display_U_B0: value = U->B.s0; break;
	case display_U_B1: value = U->B.s1; break;
	case display_U_B2: value = U->B.s2; break;
	case display_U_B: value = sqrt(BSq(*U)); break;
	case display_U_energy: value = .5 * (coordLen(U->epsE) + coordLen(U->B) / mu0); break;
	}
]]
end

return Maxwell
