local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'
local clnumber = require 'clnumber'

local SRHD = class(Equation)
SRHD.name = 'SRHD'
SRHD.numStates = 5

SRHD.consVars = {'D', 'Sx', 'Sy', 'Sz', 'tau'}
SRHD.primVars = {'rho', 'vx', 'vy', 'vz', 'eInt'}
SRHD.mirrorVars = {{'Sx'}, {'Sy'}, {'Sz'}}
SRHD.displayVars = {
	'D',
	'Sx', 'Sy', 'Sz', 'S',
	'tau',
	'W',
	'primitive_reconstruction_error',
}
SRHD.primDisplayVars = {
	'rho',
	'vx', 'vy', 'vz', 'v',
	'eInt',
	'P',
	'h',
}

SRHD.initStates = require 'init.euler'
SRHD.initStateNames = table.map(SRHD.initStates, function(info) return info.name end)

local GuiFloat = require 'guivar.float'
local GuiInt = require 'guivar.int'
SRHD.guiVars = table{
	GuiFloat{name='gamma', value=7/5},

	-- setting max iter to 100+ makes it freeze initially -- meaning the initial cons to prim or something is taking too long ...
	GuiInt{name='solvePrimMaxIter', value=10},	-- value=1000},
	
	GuiFloat{name='solvePrimStopEpsilon', value=1e-7},
	
	-- used by pressure solver
	-- velocity epsilon is how close we can get to the speed of light
	-- set ylabel "Lorentz factor"; set xlabel "velocity epsilon -log10"; set log xy; plot [1:10] 1/sqrt(1-(1-10**(-x))**2);
	--velEpsilon = 1e-5	-- <=> handles up to W = 500
	--velEpsilon = 1e-6	-- <=> handles up to W = 600
	--velEpsilon = 1e-7	-- <=> handles up to W = 2,000
	--velEpsilon = 1e-10	-- <=> handles up to W = 100,000
	GuiFloat{name='solvePrimVelEpsilon', value=1e-15},	-- <=> smaller than 1e-15 gnuplot x11 terminal breaks down past W = 1e+7 ...
	
	GuiFloat{name='solvePrimPMinEpsilon', value=1e-16},
	
	GuiFloat{name='rhoMin', value=1e-15},
	GuiFloat{name='rhoMax', value=1e+20},
	GuiFloat{name='eIntMax', value=1e+20},
	GuiFloat{name='DMin', value=1e-15},
	GuiFloat{name='DMax', value=1e+20},
	GuiFloat{name='tauMin', value=1e-15},
	GuiFloat{name='tauMax', value=1e+20},
}
SRHD.guiVarsForName = SRHD.guiVars:map(function(var) return var, var.name end)

function SRHD:getTypeCode()
	return [[
typedef struct {
	real rho;
	union {
		struct { real vx, vy, vz; };
		real v[3];
	};
	real eInt;
} prim_t;

enum {
	cons_D,
	cons_Sx,
	cons_Sy,
	cons_Sz,
	cons_tau,
};

typedef struct {
	real D;
	union {
		struct { real Sx, Sy, Sz; };
		real S[3];
	};
	real tau;
} cons_t;
]]
end

function SRHD:getCodePrefix()
	return table{
		SRHD.super.getCodePrefix(self),
		[[
#define gamma_1 (gamma-1.)

real calc_P(real rho, real eInt) { return gamma_1 * rho * eInt; }	//pressure function for ideal gas
real calc_dP_drho(real rho, real eInt) { return gamma_1 * eInt; }	//chi in most papers
real calc_dP_deInt(real rho, real eInt) { return gamma_1 * rho; }	//kappa in most papers
real calc_eInt_from_P(real rho, real P) { return P / (gamma_1 * rho); }
real calc_h(real rho, real P, real eInt) { return 1. + eInt + P / rho; }

cons_t consFromPrim(prim_t prim) {
	real rho = prim.rho;
	real vx = prim.vx, vy = prim.vy, vz = prim.vz; 
	real eInt = prim.eInt;
	real vSq = vx*vx + vy*vy + vz*vz;
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(rho, eInt);
	real h = calc_h(rho, P, eInt);
	real D = rho * W;	//rest-mass density
	real Sx = rho * h * WSq * vx;
	real Sy = rho * h * WSq * vy;
	real Sz = rho * h * WSq * vz;
	real tau = rho * h * WSq - P - D;	
	return (cons_t){.D=D, .Sx=Sx, .Sy=Sy, .Sz=Sz, .tau=tau};
}
]],
	}:concat'\n'
end

function SRHD:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..(solver.initStatePtr[0]+1))
	local code = initState.init(solver)
	return table{
		[[
__kernel void initState(
	__global cons_t* consBuf,
	__global prim_t* primBuf
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
	
]]..code..[[
	real eInt = calc_eInt_from_P(rho, P);
	real vSq = vx*vx + vy*vy + vz*vz;
	real W = 1./sqrt(1. - vSq);
	real h = calc_h(rho, P, eInt);

	prim_t prim = {.rho=rho, .vx=vx, .vy=vy, .vz=vz, .eInt=eInt};
	primBuf[index] = prim;
	consBuf[index] = consFromPrim(prim);
}
]],
	}:concat'\n'
end

function SRHD:getSolverCode(solver)
	return table{
		'#include "eqn/srhd.cl"',
	}:concat'\n'
end

-- handled by the SRHDRoe solver, so it can accept UBuf and primBuf
SRHD.getCalcDisplayVarCode = nil

return SRHD
