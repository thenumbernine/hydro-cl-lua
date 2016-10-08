local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Euler1D = class(Equation)
Euler1D.name = 'Euler1D'

Euler1D.numStates = 3

Euler1D.consVars = {'rho', 'mx', 'ETotal'}
Euler1D.primVars = {'rho', 'vx', 'P'}

-- notice the current cheap system just appends the dim to mirror on the var prefix
-- but in 1D, there is only 'mx', not 'my' or 'mz'
-- soo... the system will break for 2D and 3D. 
-- soo ... fix the system
Euler1D.mirrorVars = {'m'}

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

-- TODO make this a gui variable, and modifyable in realtime?
Euler1D.gamma = 7/5

function Euler1D:getTypeCode()
	return 
		require 'makestruct'('prim_t', self.primVars) .. '\n' ..
		Euler1D.super.getTypeCode(self) 
end

function Euler1D:getInitStateCode(solver, clnumber)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local initStateDefLines = '#define INIT_STATE_CODE \\\n'
		.. initState.code:gsub('\n', '\\\n')
	
	self.gamma = initState.gamma
	
	return table{
		'#define gamma '..clnumber(self.gamma),
		initStateDefLines,
		[[
#define gamma_1 (gamma-1.)
#define gamma_3 (gamma-3.)

cons_t consFromPrim(prim_t W) {
	return (cons_t){
		.rho = W.rho,
		.mx = W.rho * W.vx,
		.ETotal = .5 * W.rho * W.vx * W.vx + W.P / gamma_1,
	};
}

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
	real rho = 0;
	real vx = 0;
	real P = 0;

	INIT_STATE_CODE

	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .P=P});
}
]]
	}:concat'\n'
end

function Euler1D:solverCode(clnumber, solver)	
	return table{
		'#define gamma '..clnumber(self.gamma),
		'#include "euler1d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler1D
