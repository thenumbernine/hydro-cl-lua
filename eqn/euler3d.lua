local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'
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

Euler3D.initStates = require 'eqn.init_euler'
Euler3D.initStateNames = table.map(Euler3D.initStates, function(info) return info.name end)

Euler3D.guiVars = {'gamma'}
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

function Euler3D:getCodePrefix()
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

prim_t primFromCons(cons_t U) {
	real EInt = U.ETotal - .5 * (U.mx * U.mx + U.my * U.my + U.mz * U.mz) / U.rho;
	return (prim_t){
		.rho = U.rho,
		.vx = U.mx / U.rho,
		.vy = U.my / U.rho,
		.vz = U.mz / U.rho,
		.P = EInt / gamma_1,
	};
}
]],
	}:concat'\n'
end

function Euler3D:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local code = initState.init(solver)	
	return table{
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
	
]]..code..[[
	
	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .vy=vy, .vz=vz, .P=P});
}
]],
	}:concat'\n'
end
	
function Euler3D:getSolverCode(solver)	
	return table{
		'#include "eqn/euler3d.cl"',
	}:concat'\n'
end

function Euler3D:getCalcDisplayVarCode()
	return [[
	prim_t W = primFromCons(*U);
	switch (displayVar) {
	case display_U_rho: value = W.rho; break;
	case display_U_vx: value = W.vx; break;
	case display_U_vy: value = W.vy; break;
	case display_U_vz: value = W.vz; break;
	case display_U_v: value = sqrt(W.vx * W.vx + W.vy * W.vy + W.vz * W.vz); break;
	case display_U_mx: value = U->mx; break;
	case display_U_my: value = U->my; break;
	case display_U_mz: value = U->mz; break;
	case display_U_m: value = sqrt(U->mx * U->mx + U->my * U->my + U->mz * U->mz); break;
	case display_U_P: value = W.P; break;
	case display_U_eInt: value = calc_eInt(W); break;
	case display_U_eKin: value = calc_eKin(W); break;
	case display_U_eTotal: value = U->ETotal / W.rho; break;
	case display_U_EInt: value = calc_EInt(W); break;
	case display_U_EKin: value = calc_EKin(W); break;
	case display_U_ETotal: value = U->ETotal; break;
	case display_U_S: value = W.P / pow(W.rho, (real)gamma); break;
	case display_U_H: value = W.P * gamma / gamma_1; break;
	case display_U_h: value = W.P * gamma / gamma_1 / W.rho; break;
	case display_U_HTotal: value = W.P * gamma / gamma_1 + .5 * W.rho * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz); break;
	case display_U_hTotal: value = W.P * gamma / gamma_1 / W.rho + .5 * (W.vx * W.vx + W.vy * W.vy + W.vz * W.vz); break;
	}
]]
end

return Euler3D
