local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'
local clnumber = require 'clnumber'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.numStates = 6
Maxwell.consVars = {'epsEx', 'epsEy', 'epsEz', 'Bx', 'By', 'Bz'}
Maxwell.mirrorVars = {{'epsEx', 'Bx'}, {'epsEy', 'By'}, {'epsEz', 'Bz'}}
Maxwell.displayVars = {
	'Ex', 'Ey', 'Ez', 'E',
	'Bx', 'By', 'Bz', 'B',
	'energy',
}

Maxwell.useSourceTerm = true

Maxwell.initStateNames = {'default'}

Maxwell.guiVars = {'eps0', 'mu0', 'sigma'}
Maxwell.eps0 = 1	-- permittivity
Maxwell.mu0 = 1		-- permeability
Maxwell.sigma = 1	-- conductivity

function Maxwell:getCodePrefix()
	return table{
		'#define eps0 '..clnumber(self.eps0),
		'#define mu0 '..clnumber(self.mu0),
		'#define sigma '..clnumber(self.sigma),
		'#define sqrt_eps0 '..clnumber(math.sqrt(self.eps0)),
		'#define sqrt_mu0 '..clnumber(math.sqrt(self.mu0)),
		[[
real ESq(cons_t U) { 
	return (U.epsEx * U.epsEx + U.epsEy * U.epsEy + U.epsEz * U.epsEz) / (eps0 * eps0);
}

real BSq(cons_t U) {
	return U.Bx * U.Bx + U.By * U.By + U.Bz * U.Bz;
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
	real4 x = CELL_X(i);
	real4 mids = (real).5 * (mins + maxs);
	bool lhs = x[0] < mids[0]
#if dim > 1
		&& x[1] < mids[1]
#endif
#if dim > 2
		&& x[2] < mids[2]
#endif
	;
	__global cons_t* U = UBuf + index;
	U->epsEx = 0;
	U->epsEy = 0;
	U->epsEz = 1 * eps0;
	U->Bx = 1;
	U->By = lhs ? 1 : -1;
	U->Bz = 0;
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
		case display_U_Ex: value = U->epsEx / eps0; break;
		case display_U_Ey: value = U->epsEy / eps0; break;
		case display_U_Ez: value = U->epsEz / eps0; break;
		case display_U_E: value = sqrt(ESq(*U)); break;
		case display_U_Bx: value = U->Bx; break;
		case display_U_By: value = U->By; break;
		case display_U_Bz: value = U->Bz; break;
		case display_U_B: value = sqrt(BSq(*U)); break;
		case display_U_energy: value = .5 * (ESq(*U) * eps0 + BSq(*U) / mu0); break;
		}
]]
end

return Maxwell
