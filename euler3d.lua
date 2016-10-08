local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'
local clnumber = require 'clnumber'

local Euler3D = class(Equation)
Euler3D.name = 'Euler3D'

Euler3D.numStates = 5

Euler3D.consVars = {'rho', 'mx', 'my', 'mz', 'ETotal'}
Euler3D.primVars = {'rho', 'vx', 'vy', 'vz', 'P'}
Euler3D.mirrorVars = {'m'}
Euler3D.displayVars = {
	'rho',
	'vx', 'vy', 'vz', 'v',
	'mx', 'my', 'mz', 'm',
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

Euler3D.initStates = {
	{
		name='Sod',
		code=[[
	rho = lhs ? 1 : .125;
	P = lhs ? 1 : .1;
]],
	},
	{
		name='Sedov',
		code=[[
	rho = 1;
	P = (i.x == gridSize.x/2 && i.y == gridSize.y/2 && i.z == gridSize.z/2) ? 1e+3 : 1;
]],
	},
	{
		name='constant',
		code='rho=1; vx=1; vy=1; vz=1; P=1;',
	},
	{
		name='linear',
		code='rho=2+x.x; P=1;',
	},
	{
		name='gaussian',
		code=[[
	real sigma = 1. / sqrt(10.);
	rho = exp(-x.x*x.x / (sigma*sigma)) + .1;
	P = 1 + .1 * (exp(-x.x*x.x / (sigma*sigma)) + 1) / (gamma_1 * rho);
]],
	},
	{
		name='advect wave',
		code=[[
	real rSq = dot(x,x);
	rho = exp(-100*rSq) + 1.;
	vx = 1;
	P = 1;
]],
	},
	{
		name='sphere',
		code=[[
	real rSq = dot(x,x);
	bool inside = rSq < .2*.2;
	rho = inside ? 1 : .1;
	P = inside ? 1 : .1;
]],
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

Euler3D.initStateNames = table.map(Euler3D.initStates, function(info) return info.name end)

-- TODO make this a gui variable, and modifyable in realtime?
Euler3D.gamma = 7/5

function Euler3D:getTypeCode()
	return [[

typedef struct { 
	real rho;
	union {
		struct { real vx, vy, vz; };
		real v[3];
	};
	real P;
} prim_t;

enum {
	cons_rho,
	cons_mx,
	cons_my,
	cons_mz,
	cons_ETotal,
};

typedef struct {
	real rho;
	union {
		struct { real mx, my, mz; };
		real m[3];
	};
	real ETotal;
} cons_t;

]]
end

function Euler3D:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local initStateDefLines = '#define INIT_STATE_CODE \\\n'
		.. initState.code:gsub('\n', '\\\n')
	
	-- TODO make this a gui variable, and modifyable in realtime?
	self.gamma = initState.gamma
	
	return table{
		'#define gamma '..clnumber(self.gamma),
		initStateDefLines,
		[[
#define gamma_1 (gamma-1.)
#define gamma_3 (gamma-3.)

real calc_hTotal(real rho, real P, real ETotal) { return (P + ETotal) / rho; }
real calc_eKin(prim_t W) { return .5 * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz); }
real calc_EKin(prim_t W) { return W.rho * calc_eKin(W); }
real calc_EInt(prim_t W) { return W.P / gamma_1; }
real calc_eInt(prim_t W) { return calc_EInt(W) / W.rho; }
real calc_ETotal(prim_t W) { return calc_EKin(W) + calc_EInt(W); }

cons_t consFromPrim(prim_t W) {
	return (cons_t){
		.rho = W.rho,
		.mx = W.rho * W.vx,
		.my = W.rho * W.vy,
		.mz = W.rho * W.vz,
		.ETotal = calc_ETotal(W),
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
	real vy = 0;
	real vz = 0;
	real P = 0;
	
	INIT_STATE_CODE
	
	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .vy=vy, .vz=vz, .P=P});
}
]]
	}:concat'\n'
end

function Euler3D:solverCode(solver)	
	return table{
		'#define gamma '..clnumber(self.gamma),
		'#include "euler3d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler3D
