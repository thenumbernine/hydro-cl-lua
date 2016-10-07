local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Euler1D = class(Equation)
Euler1D.name = 'Euler1D'

Euler1D.numStates = 3

Euler1D.consVars = table{'rho', 'mx', 'ETotal'}
Euler1D.primVars = table{'rho', 'vx', 'P'}
Euler1D.displayVars = {
	'rho',
	'vx',
	'mx',
	'eInt',
	'eKin', 
	'eTotal', 
	'EInt', 
	'EKin', 
	'ETotal', 
	'P',
	'S', 
	'h',
	'H', 
	'hTotal',
	'HTotal',
} 

Euler1D.initStates = {
	'constant',
	'linear',
	'gaussian',
	'rarefaction_wave',
	'Sod',
	'Sedov',
	'shock_wave',
	'relativistic_blast_wave_interaction',
	'relativistic_blast_wave_test_problem_1',
}

Euler1D.gamma = 7/5

function Euler1D:getTypeCode()
	return 
		require 'makestruct'('prim_t', self.primVars) .. '\n' ..
		Euler1D.super.getTypeCode(self) 
end

function Euler1D:solverCode(clnumber)
	return table{
		'#define gamma '..clnumber(self.gamma),
		'#include "euler1d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler1D
