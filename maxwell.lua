local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.numStates = 6
Maxwell.consVars = {'epsEx', 'epsEy', 'epsEz', 'Bx', 'By', 'Bz'}
Maxwell.mirrorVars = {'epsE', 'B'}
Maxwell.displayVars = {
	'Ex', 'Ey', 'Ez', 'E',
	'Bx', 'By', 'Bz', 'B',
	'energy',
}

Maxwell.initStateNames = {'default'}

function Maxwell:solverCode(clnumber, solver)
	local eps0 = 1
	local mu0 = 1
	return 
		'#define eps0 '..clnumber(eps0)..'\n'..
		'#define mu0 '..clnumber(mu0)..'\n'..
		'#include "maxwell.cl"'
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
