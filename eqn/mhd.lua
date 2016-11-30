local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local Equation = require 'eqn.eqn'

local MHD = class(Equation)

MHD.name = 'MHD'

MHD.numStates = 8
MHD.numWaves = 7

MHD.consVars = {'rho', 'mx', 'my', 'mz', 'ETotal', 'bx', 'by', 'bz'}
MHD.primVars = {'rho', 'vx', 'vy', 'vz', 'P', 'bx', 'by', 'bz'}
MHD.mirrorVars = {{'m.x', 'b.x'}, {'m.y', 'b.y'}, {'m.z', 'b.z'}}

MHD.hasEigenCode = true

-- hmm, we want init.euler and init.mhd here ...
MHD.initStates = require 'init.euler'

local GuiFloat = require 'guivar.float'

MHD.guiVars = table{
	GuiFloat{name='heatCapacityRatio', value=5/3},
	GuiFloat{name='mu0', value=1},
}

function MHD:getTypeCode()
	return [[
typedef struct {
	real rho;
	real3 v;
	real P;
	real3 b;
} prim_t;

typedef union {
	real ptr[8];
	struct {
		real rho;
		real3 m;
		real ETotal;
		real3 b;
	};
} cons_t;
]]
end

function MHD:getCodePrefix()
	return table{
		MHD.super.getCodePrefix(self),
		[[
inline real calc_eKin(prim_t W) { return .5 * real3_lenSq(W.v); }
inline real calc_EKin(prim_t W) { return W.rho * calc_eKin(W); }
inline real calc_EInt(prim_t W) { return W.P / (heatCapacityRatio - 1.); }
inline real calc_eInt(prim_t W) { return calc_EInt(W) / W.rho; }
inline real calc_EMag(prim_t W) { return .5 * real3_lenSq(W.b); }
inline real calc_eMag(prim_t W) { return calc_EMag(W) / W.rho; }
inline real calc_PMag(prim_t W) { return .5 * real3_lenSq(W.b); }
inline real calc_EHydro(prim_t W) { return calc_EKin(W) + calc_EInt(W); }
inline real calc_eHydro(prim_t W) { return calc_EHydro(W) / W.rho; }
inline real calc_ETotal(prim_t W) { return calc_EKin(W) + calc_EInt(W) + calc_EMag(W); }
inline real calc_eTotal(prim_t W) { return calc_ETotal(W) / W.rho; }
inline real calc_H(real P) { return P * (heatCapacityRatio / (heatCapacityRatio - 1.)); }
inline real calc_h(real rho, real P) { return calc_H(P) / rho; }
inline real calc_HTotal(prim_t W, real ETotal) { return W.P + calc_PMag(W) + ETotal; }
inline real calc_hTotal(prim_t W, real ETotal) { return calc_HTotal(W, ETotal) / W.rho; }
inline real calc_Cs(prim_t W) { return sqrt(heatCapacityRatio * W.P / W.rho); }

inline prim_t primFromCons(cons_t U) {
	prim_t W;
	W.rho = U.rho;
	W.v = real3_scale(U.m, 1./U.rho);
	W.b = U.b;
	real vSq = real3_lenSq(W.v);
	real bSq = real3_lenSq(W.b);
	real EKin = .5 * U.rho * vSq;
	real EMag = .5 * bSq;
	real EInt = U.ETotal - EKin - EMag;
	W.P = EInt * (heatCapacityRatio - 1.);
	//W.P = max(W.P, 1e-7);
	//W.rho = max(W.rho, 1e-7);
	return W;
}

inline cons_t consFromPrim(prim_t W) {
	cons_t U;
	U.rho = W.rho;
	U.m = real3_scale(W.v, W.rho);
	U.b = W.b;
	real vSq = real3_lenSq(W.v);
	real bSq = real3_lenSq(W.b);
	real EKin = .5 * W.rho * vSq;
	real EMag = .5 * bSq;
	real EInt = W.P / (heatCapacityRatio - 1.);
	U.ETotal = EInt + EKin + EMag;
	return U;
}
]],
	}:concat'\n'
end

function MHD:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local code = initState.init(solver)	
	return [[
kernel void initState(
	global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
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
	real3 b = _real3(0,0,0);

]]..code..[[

	prim_t W = {.rho=rho, .v=v, .P=P, .b=b};
	UBuf[index] = consFromPrim(W);
}
]]
end

function MHD:getSolverCode(solver)
	return template(file['eqn/mhd.cl'], {eqn=self, solver=solver})
end

MHD.displayVarCodePrefix = [[
	cons_t U = buf[index];
	prim_t W = primFromCons(U);
]]

function MHD:getDisplayVars(solver)
	return {
		{rho = 'value = W.rho;'},
		{vx = 'value = W.v.x;'},
		{vy = 'value = W.v.y;'},
		{vz = 'value = W.v.z;'},
		{v = 'value = real3_len(W.v);'},
		{mx = 'value = U.m.x;'},
		{my = 'value = U.m.y;'},
		{mz = 'value = U.m.z;'},
		{m = 'value = real3_len(U.m);'},
		{bx = 'value = W.b.x;'},
		{by = 'value = W.b.y;'},
		{bz = 'value = W.b.z;'},
		{b = 'value = real3_len(W.b);'},
		{P = 'value = W.P;'},
		--{PMag = 'value = calc_PMag(W);'},
		--{PTotal = 'value = W.P + calc_PMag(W);'},
		--{eInt = 'value = calc_eInt(W);'},
		{EInt = 'value = calc_EInt(W);'},
		--{eKin = 'value = calc_eKin(W);'},
		{EKin = 'value = calc_EKin(W);'},
		--{eHydro = 'value = calc_eHydro(W);'},
		{EHydro = 'value = calc_EHydro(W);'},
		--{eMag = 'value = calc_eMag(W);'},
		{EMag = 'value = calc_EMag(W);'},
		--{eTotal = 'value = U.ETotal / W.rho;'},
		{ETotal = 'value = U.ETotal;'},
		{S = 'value = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		{H = 'value = calc_H(W.P);'},
		--{h = 'value = calc_H(W.P) / W.rho;'},
		--{HTotal = 'value = calc_HTotal(W, U.ETotal);'},
		--{hTotal = 'value = calc_hTotal(W, U.ETotal);'},
		--{Cs = 'value = calc_Cs(W); },
		{['primitive reconstruction error'] = [[
		//prim have just been reconstructed from cons
		//so reconstruct cons from prims again and calculate the difference
		cons_t U2 = consFromPrim(W);
		value = 0;
		for (int j = 0; j < numStates; ++j) {
			value += fabs(U.ptr[j] - U2.ptr[j]);
		}
]]},
	}
end

function Equation:getEigenTypeCode(solver)
	return template([[
typedef struct {
	real evL[7*7];
	real evR[7*7];
<? if solver.checkFluxError then ?>
	real A[7*7];
<? end ?>
} eigen_t;
]], {
		solver = solver,
	})
end

return MHD
