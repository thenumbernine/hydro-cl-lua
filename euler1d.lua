local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Euler1D = class(Equation)
Euler1D.name = 'Euler1D'

Euler1D.numStates = 3

Euler1D.consVars = {'rho', 'mx', 'ETotal'}
Euler1D.primVars = {'rho', 'vx', 'P'}
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
	{
		name='Sod',
		code=[[
	rho = lhs ? 1 : .125;
	P = lhs ? 1 : .1;
]]
	},
	{
		name='Sedov',
		code=[[
	rho = 1;
	P = (i.x == gridSize.x/2 && i.y == gridSize.y/2 && i.z == gridSize.z/2) ? 1e+3 : 1;
]]
	},
	{
		name='constant',
		code='rho=2+x.x; P=1;'
	},
	{
		name='linear',
		code='rho=2+x.x; P=1;'
	},
	{
		name='gaussian',
		code=[[
	real sigma = 1. / sqrt(10.);
	rho = exp(-x.x*x.x / (sigma*sigma)) + .1;
	P = 1 + .1 * (exp(-x.x*x.x / (sigma*sigma)) + 1) / (gamma_1 * rho);
]]
	},
	{
		name='rarefaction_wave',
		code=[[
	real delta = .1;
	rho = 1;	// lhs ? .2 : .8;
	vx = lhs ? .5 - delta : .5 + delta;
	P = 1;
]]
	},

	--from SRHD Marti & Muller 2000
	{
		name='shock_wave',
		code=[[
	rho = 1;
	vx = lhs ? .5 : 0;
	P = lhs ? 1e+3 : 1;
]]
	},
	{
		name='relativistic_blast_wave_interaction',
		code=[[
	real xL = .9 * mins_x + .1 * maxs_x;
	real xR = .1 * mins_x + .9 * maxs_x;
	rho = 1;
	P = x.x < xL ? 1000 : (x.x > xR ? 100 : .01);
]]
	},
	{
		name='relativistic_blast_wave_test_problem_1',
		gamma = 5/3,
		code=[[
	rho = lhs ? 10 : 1;
	P = gamma_1 * rho * (lhs ? 2 : 1e-6);
]]
	},
}

Euler1D.initStateNames = table.map(Euler1D.initStates, function(info) return info.name end)

function Euler1D:getTypeCode()
	return 
		require 'makestruct'('prim_t', self.primVars) .. '\n' ..
		Euler1D.super.getTypeCode(self) 
end

function Euler1D:solverCode(clnumber, solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local initStateDefLines = '#define INIT_STATE_CODE \\\n'
		.. initState.code:gsub('\n', '\\\n')
	
	-- TODO make this a gui variable, and modifyable in realtime?
	local gamma = initState.gamma or 7/5
	
	return table{
		'#define gamma '..clnumber(gamma),
		initStateDefLines,
		'#include "euler1d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler1D
