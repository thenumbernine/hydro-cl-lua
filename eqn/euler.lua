local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'
local clnumber = require 'clnumber'

local Euler = class(Equation)
Euler.name = 'Euler'

Euler.numStates = 5

Euler.consVars = {'rho', 'm0', 'm1', 'm2', 'ETotal'}
Euler.primVars = {'rho', 'v0', 'v1', 'v2', 'P'}
Euler.mirrorVars = {{'m.s0'}, {'m.s1'}, {'m.s2'}}
Euler.displayVars = {
	'rho',
	'v0', 'v1', 'v2', 'v',
	'm0', 'm1', 'm2', 'm',
	'eInt',
	'eKin', 
	'ePot',
	'eTotal', 
	'EInt', 
	'EKin', 
	'EPot',
	'ETotal', 
	'P',
	'S', 
	'h',
	'H', 
	'hTotal',
	'HTotal',
} 

Euler.initStates = require 'init.euler'
Euler.initStateNames = table.map(Euler.initStates, function(info) return info.name end)

Euler.guiVars = table{
	require 'guivar.float'{name='gamma', value=7/5}
}
Euler.guiVarsForName = Euler.guiVars:map(function(var) return var, var.name end)

function Euler:getTypeCode()
	return [[

typedef struct { 
	real rho;
	real3 v;
	real P;
} prim_t;

enum {
	cons_rho,
	cons_m0,
	cons_m1,
	cons_m2,
	cons_ETotal,
};

typedef struct {
	real rho;
	real3 m;
	real ETotal;
} cons_t;
]]
end

function Euler:getCodePrefix()
	return table{
		Euler.super.getCodePrefix(self),
		[[
#define gamma_1 (gamma-1.)
#define gamma_3 (gamma-3.)

real calc_H(real P) {
	return P * (gamma / gamma_1);
}

real calc_h(real rho, real P) {
	return calc_H(P) / rho;
}

real calc_hTotal(real rho, real P, real ETotal) { 
	return (P + ETotal) / rho;
}

real calc_HTotal(real P, real ETotal) {
	return P + ETotal; 
}

real calc_eKin(prim_t W) { 
	return .5 * (W.v.x*W.v.x + W.v.y*W.v.y + W.v.z*W.v.z);//coordLenSq(W.v);
}

real calc_EKin(prim_t W) {
	return W.rho * calc_eKin(W);
}

real calc_EInt(prim_t W) {
	return W.P / gamma_1;
}

real calc_eInt(prim_t W) {
	return calc_EInt(W) / W.rho;
}

real calc_EKin_fromCons(cons_t U) {
	return .5 * (U.m.x*U.m.x + U.m.y*U.m.y + U.m.z*U.m.z /*coordLenSq(U.m)*/ ) / U.rho;
}

real calc_ETotal(prim_t W, real ePot) {
	real EPot = W.rho * ePot;
	return calc_EKin(W) + calc_EInt(W) + EPot;
}

real calc_Cs(prim_t W) {
	return sqrt(gamma * W.P / W.rho);
}

prim_t primFromCons(cons_t U, real ePot) {
	real EPot = U.rho * ePot;
	real EKin = calc_EKin_fromCons(U);
	real EInt = U.ETotal - EPot - EKin;
	return (prim_t){
		.rho = U.rho,
		.v = real3_scale(U.m, 1./U.rho),
		.P = EInt / gamma_1,
	};
}
]],
	}:concat'\n'
end

function Euler:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local code = initState.init(solver)	
	return table{
		[[
cons_t consFromPrim(prim_t W, real ePot) {
	return (cons_t){
		.rho = W.rho,
		.m = real3_scale(W.v, W.rho),
		.ETotal = calc_ETotal(W, ePot),
	};
}

__kernel void initState(
	__global cons_t* UBuf,
	__global real* ePotBuf
) {
	SETBOUNDS(0,0);
	
	//TODO should 'x' be in embedded or coordinate space? coordinate 
	//should 'v0 v1 v2' be embedded or coordinate space? coordinate
	real3 x = CELL_X(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.s0 < mids.s0
#if dim > 1
		&& x.s1 < mids.s1
#endif
#if dim > 2
		&& x.s2 < mids.s2
#endif
	;
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	real ePot = 0;

]]..code..[[

	prim_t W = {.rho=rho, .v=v, .P=P};
	UBuf[index] = consFromPrim(W, ePot);
	ePotBuf[index] = ePot;
}
]],
	}:concat'\n'
end
	
function Euler:getSolverCode(solver)	
	return table{
		'#include "eqn/euler.cl"',
	}:concat'\n'
end

function Euler:getCalcDisplayVarCode()
	return [[
	prim_t W = primFromCons(U, ePot);
	//switch (displayVar) {
	switch (displayVar) {
	case display_U_rho: value = W.rho; break;
	case display_U_v0: value = W.v.s0; break;
	case display_U_v1: value = W.v.s1; break;
	case display_U_v2: value = W.v.s2; break;
	case display_U_v: value = sqrt(W.v.x*W.v.x + W.v.y*W.v.y + W.v.z*W.v.z); /*coordLen(W.v);*/ break;
	case display_U_m0: value = U.m.s0; break;
	case display_U_m1: value = U.m.s1; break;
	case display_U_m2: value = U.m.s2; break;
	case display_U_m: value = sqrt(U.m.x*U.m.x + U.m.y*U.m.y + U.m.z*U.m.z); /*coordLen(U.m);*/ break;
	case display_U_P: value = W.P; break;
	case display_U_eInt: value = calc_eInt(W); break;
	case display_U_eKin: value = calc_eKin(W); break;
	case display_U_ePot: value = ePot; break;
	case display_U_eTotal: value = U.ETotal / W.rho; break;
	case display_U_EInt: value = calc_EInt(W); break;
	case display_U_EKin: value = calc_EKin(W); break;
	case display_U_EPot: value = W.rho * ePot; break;
	case display_U_ETotal: value = U.ETotal; break;
	case display_U_S: value = W.P / pow(W.rho, (real)gamma); break;
	case display_U_H: value = calc_H(W.P); break;
	case display_U_h: value = calc_h(W.rho, W.P); break;
	case display_U_HTotal: value = calc_HTotal(W.P, U.ETotal); break;
	case display_U_hTotal: value = calc_hTotal(W.rho, W.P, U.ETotal); break;
	}
]]
end

return Euler
