local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'clnumber'

local SRHD = class(Equation)
SRHD.name = 'SRHD'
SRHD.numStates = 5

SRHD.consVars = {'D', 'S0', 'S1', 'S2', 'tau'}
SRHD.primVars = {'rho', 'v0', 'v1', 'v2', 'eInt'}
SRHD.mirrorVars = {{'S.s0'}, {'S.s1'}, {'S.s2'}}
SRHD.displayVars = {
	'D',
	'S0', 'S1', 'S2', 'S',
	'tau',
	'W',
	'primitive_reconstruction_error',
}
SRHD.primDisplayVars = {
	'rho',
	'v0', 'v1', 'v2', 'v',
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
	real3 v;
	real eInt;
} prim_t;

enum {
	cons_D,
	cons_S0,
	cons_S1,
	cons_S2,
	cons_tau,
};

typedef struct {
	real D;
	real3 S;
	real tau;
} cons_t;
]]
end

function SRHD:getCodePrefix()
	return table{
		SRHD.super.getCodePrefix(self),
		[[
#define gamma_1 (gamma-1.)

//pressure function for ideal gas
real calc_P(real rho, real eInt) {
	return gamma_1 * rho * eInt;
}	

//chi in most papers
real calc_dP_drho(real rho, real eInt) {
	return gamma_1 * eInt;
}

//kappa in most papers
real calc_dP_deInt(real rho, real eInt) {
	return gamma_1 * rho;
}

real calc_eInt_from_P(real rho, real P) {
	return P / (gamma_1 * rho);
}

real calc_h(real rho, real P, real eInt) {
	return 1. + eInt + P / rho;
}

cons_t consFromPrim(prim_t prim) {
	real vSq = coordLenSq(prim.v);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);
	real D = prim.rho * W;	//rest-mass density
	real3 S = real3_scale(prim.v, prim.rho * h * WSq);
	real tau = prim.rho * h * WSq - P - D;	
	return (cons_t){.D=D, .S=S, .tau=tau};
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
	real3 x = CELL_X(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	
]]..code..[[
	
	real eInt = calc_eInt_from_P(rho, P);
	real vSq = coordLenSq(v);
	real W = 1./sqrt(1. - vSq);
	real h = calc_h(rho, P, eInt);

	prim_t prim = {.rho=rho, .v=v, .eInt=eInt};
	primBuf[index] = prim;
	consBuf[index] = consFromPrim(prim);
}
]],
	}:concat'\n'
end

function SRHD:getSolverCode(solver)
	return table{
		require 'processcl'(assert(file['eqn/srhd.cl']), {solver=solver}),
	}:concat'\n'
end

-- handled by the SRHDRoe solver, so it can accept UBuf and primBuf
SRHD.getCalcDisplayVarCode = nil

return SRHD
