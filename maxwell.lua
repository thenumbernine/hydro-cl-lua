local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'
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

-- gui vars:
Maxwell.eps0 = 1	-- permittivity
Maxwell.mu0 = 1		-- permeability
Maxwell.sigma = 1	-- conductivity

function Maxwell:codePrefix()
	return table{
		'#define eps0 '..clnumber(self.eps0),
		'#define mu0 '..clnumber(self.mu0),
		'#define sigma '..clnumber(self.sigma),
		'#define sqrt_eps0 '..clnumber(math.sqrt(self.eps0)),
		'#define sqrt_mu0 '..clnumber(math.sqrt(self.mu0)),
	}:concat'\n'
end

function Maxwell:getInitStateCode(solver)
	return table{
		self:codePrefix(),
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
	U->By = 0;
	U->Bz = lhs ? 1 : -1;
}
]],
	}:concat'\n'
end

function Maxwell:solverCode(solver)
	return table{ 
		self:codePrefix(),
		'#include "maxwell.cl"',
	}:concat'\n'
end

function Maxwell:getEigenInfo()
	local eigenType = 'eigen_t'
	return {
		type = eigenType,
		typeCode = 'typedef struct { char mustbesomething; } eigen_t;',	-- can it be zero sized?
		code = nil,
		displayVars = {},
	}
end

return Maxwell
