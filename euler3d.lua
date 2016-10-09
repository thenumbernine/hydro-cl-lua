local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'
local clnumber = require 'clnumber'

local Euler3D = class(Equation)
Euler3D.name = 'Euler3D'

Euler3D.numStates = 5

Euler3D.consVars = {'rho', 'mx', 'my', 'mz', 'ETotal'}
Euler3D.primVars = {'rho', 'vx', 'vy', 'vz', 'P'}
Euler3D.mirrorVars = {{'mx'}, {'my'}, {'mz'}}
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

local function quadrantProblem(args)
	args.cfl = .475
	args.boundary = {
		xmin = 'freeflow',
		xmax = 'freeflow',
		ymin = 'freeflow',
		ymax = 'freeflow',
	}
	local function build(i)
		local q = args[i]
		return table.map(q, function(v,k,t)
			return k..'='..v..';', #t+1
		end):concat' '
	end
	args.code = [[
	bool xp = x.x > mids.x;
	bool yp = x.y > mids.y;
	if (yp) {
		if (xp) {
			]]..build(1)..[[
		} else {
			]]..build(2)..[[	
		}
	} else {
		if (!xp) {
			]]..build(3)..[[
		} else {
			]]..build(4)..[[
		}
	}
]]
	return args
end

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
	real xSq = dot(x,x);
	rho = exp(-xSq / (sigma*sigma)) + .1;
	P = 1 + .1 * (exp(-xSq / (sigma*sigma)) + 1) / (gamma_1 * rho);
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
	-- http://www.cfd-online.com/Wiki/Explosion_test_in_2-D
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
		name='rarefaction wave',
		code=[[
	real delta = .1;
	rho = 1;	// lhs ? .2 : .8;
	vx = lhs ? .5 - delta : .5 + delta;
	P = 1;
]]
	},
	
	-- 2D tests described in Alexander Kurganov, Eitan Tadmor, Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers
	--  which says it is compared with  C. W. Schulz-Rinne, J. P. Collins, and H. M. Glaz, Numerical solution of the Riemann problem for two-dimensional gas dynamics
	-- and I can't find that paper right now
	quadrantProblem{
		name = 'configuration 1',
		{rho=1, P=1, vx=0, vy=0},
		{rho=.5197, P=.4, vx=-.7259, vy=0},
		{rho=.1072, P=.0439, vx=-.7259, vy=-1.4045},
		{rho=.2579, P=.15, vx=0, vy=-1.4045},
	},
	quadrantProblem{
		name = 'configuration 2',
		{rho=1, P=1, vx=0, vy=0},
		{rho=.5197, P=.4, vx=-.7259, vy=0},
		{rho=1, P=1, vx=-.7259, vy=-.7259},
		{rho=.5197, P=.4, vx=0, vy=-.7259},
	},
	quadrantProblem{
		name = 'configuration 3',
		{rho=1.5, P=1.5, vx=0, vy=0},
		{rho=.5323, P=.3, vx=1.206, vy=0},
		{rho=.138, P=.029, vx=1.206, vy=1.206},
		{rho=.5323, P=.3, vx=0, vy=1.206},
	},
	quadrantProblem{
		name = 'configuration 4',
		{rho=1.1, P=1.1, vx=0, vy=0},
		{rho=.5065, P=.35, vx=.8939, vy=0},
		{rho=1.1, P=1.1, vx=.8939, vy=.8939},
		{rho=.5065, P=.35, vx=0, vy=.8939},
	},
	quadrantProblem{
		name = 'configuration 5',
		{rho=1, P=1, vx=-.75, vy=-.5},
		{rho=2, P=1, vx=-.75, vy=.5},
		{rho=1, P=1, vx=.75, vy=.5},
		{rho=3, P=1, vx=.75, vy=-.5},
	},
	quadrantProblem{
		name = 'configuration 6',
		{rho=1, P=1, vx=.75, vy=-.5},
		{rho=2, P=1, vx=.75, vy=.5},
		{rho=1, P=1, vx=-.75, vy=.5},
		{rho=3, P=1, vx=-.75, vy=-.5},
	},
	--from SRHD Marti & Muller 2000
	{
		name='relativistic shock wave',
		code=[[
	rho = 1;
	vx = lhs ? .5 : 0;
	P = lhs ? 1e+3 : 1;
]]
	},
	{
		name='relativistic blast wave interaction',
		code=[[
	real xL = .9 * mins_x + .1 * maxs_x;
	real xR = .1 * mins_x + .9 * maxs_x;
	rho = 1;
	P = x.x < xL ? 1000 : (x.x > xR ? 100 : .01);
]]
	},
	{
		name='relativistic blast wave test problem 1',
		gamma = 5/3,
		code=[[
	rho = lhs ? 10 : 1;
	P = gamma_1 * rho * (lhs ? 2 : 1e-6);
]]
	},
	{
		name='Colella-Woodward',
		boundary={
			xmin='freeflow',
			xmax='freeflow',
			ymin='freeflow',
			ymax='freeflow',
		},
		code=[[
	rho = 1;
	if (x.x < -.4) {
		P = 1000;
	} else if (x.x < .4) {
		P = .01;
	} else {
		P = 100;
	}
]],
	},
	{
		name = 'Kelvin-Hemholtz',
		boundary = {
			xmin = 'periodic',
			xmax = 'periodic',
			ymin = 'periodic',
			ymax = 'periodic',
			zmin = 'periodic',
			zmax = 'periodic',
		},
		code = [[
	bool inside = (x.y > -.25 && x.y < .25);
	real theta = (x.x - mins.x) / (maxs.x - mins.x) * 2. * M_PI;
#if dim == 3
	theta *= (x.z - mins.z) / (maxs.z - mins.z);
#endif
	real noise = (maxs.x - mins.x) * 2e-5;
	rho = inside ? 2 : 1;
	vx = cos(theta) * noise + (inside ? -.5 : .5);
	vy = sin(theta) * noise;
	vz = sin(theta) * noise;
	P = 2.5;
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

function Euler3D:codePrefix()
	return table{
		'#define gamma '..clnumber(self.gamma),
		[[
#define gamma_1 (gamma-1.)
#define gamma_3 (gamma-3.)

real calc_hTotal(real rho, real P, real ETotal) { return (P + ETotal) / rho; }
real calc_eKin(prim_t W) { return .5 * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz); }
real calc_EKin(prim_t W) { return W.rho * calc_eKin(W); }
real calc_EInt(prim_t W) { return W.P / gamma_1; }
real calc_eInt(prim_t W) { return calc_EInt(W) / W.rho; }
real calc_ETotal(prim_t W) { return calc_EKin(W) + calc_EInt(W); }
]],
	}:concat'\n'
end

function Euler3D:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	
	-- TODO make this a gui variable, and modifyable in realtime?
	self.gamma = initState.gamma
	
	-- TODO either (a) do this for all Equation:getInitStateCode
	-- or (b) have the initStates call a function that can modify solver
	if initState.cfl then
		solver.cfl[0] = initState.cfl
	end
	if initState.boundary then
		for _,x in ipairs(table{'x', 'y', 'z'}:sub(1,solver.dim)) do
			for _,minmax in ipairs{'min', 'max'} do
				local var = x..minmax
				local method = initState.boundary[var]
				if method then
					solver.boundaryMethods[var][0] = assert(solver.app.boundaryMethods:find(method))-1
				end
			end
		end
	end

	return table{
		self:codePrefix(),
		[[
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
	
]]..initState.code..[[
	
	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .vy=vy, .vz=vz, .P=P});
}
]]
	}:concat'\n'
end

function Euler3D:solverCode(solver)	
	return table{
		self:codePrefix(),
		'#include "euler3d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler3D
