local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local Equation = require 'eqn.eqn'

local Euler = class(Equation)
Euler.name = 'Euler'

Euler.numStates = 5

Euler.mirrorVars = {{'m.x'}, {'m.y'}, {'m.z'}} 

Euler.hasEigenCode = true
Euler.hasCalcDT = true

Euler.initStates = require 'init.euler'

Euler.guiVars = {
	require 'guivar.float'{name='heatCapacityRatio', value=7/5}
}

function Euler:getTypeCode()
	return template([[
typedef union {
	real ptr[5];
	struct { 
		real rho;
		real3 v;
		real P;
	};
} <?=eqn.prim_t?>;

typedef union {
	real ptr[5];
	struct {
		real rho;
		real3 m;
		real ETotal;
	};
} <?=eqn.cons_t?>;
]], {
	eqn = self,
})
end

function Euler:getCodePrefix()
	return table{
		Euler.super.getCodePrefix(self),
		template([[

inline real calc_H(real P) { return P * (heatCapacityRatio / (heatCapacityRatio - 1.)); }
inline real calc_h(real rho, real P) { return calc_H(P) / rho; }
inline real calc_hTotal(real rho, real P, real ETotal) { return (P + ETotal) / rho; }
inline real calc_HTotal(real P, real ETotal) { return P + ETotal; }
inline real calc_eKin(<?=eqn.prim_t?> W) { return .5 * coordLenSq(W.v); }
inline real calc_EKin(<?=eqn.prim_t?> W) { return W.rho * calc_eKin(W); }
inline real calc_EInt(<?=eqn.prim_t?> W) { return W.P / (heatCapacityRatio - 1.); }
inline real calc_eInt(<?=eqn.prim_t?> W) { return calc_EInt(W) / W.rho; }
inline real calc_EKin_fromCons(<?=eqn.cons_t?> U) { return .5 * coordLenSq(U.m) / U.rho; }
inline real calc_ETotal(<?=eqn.prim_t?> W, real ePot) {
	real EPot = W.rho * ePot;
	return calc_EKin(W) + calc_EInt(W) + EPot;
}

inline real calc_Cs(const <?=eqn.prim_t?>* W) {
	return sqrt(heatCapacityRatio * W->P / W->rho);
}

inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real ePot) {
	real EPot = U.rho * ePot;
	real EKin = calc_EKin_fromCons(U);
	real EInt = U.ETotal - EPot - EKin;
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_scale(U.m, 1./U.rho),
		.P = EInt / (heatCapacityRatio - 1.),
	};
}

<?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real ePot) {
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_scale(W.v, W.rho),
		.ETotal = calc_ETotal(W, ePot),
	};
}
]], {
	eqn = self,
})
	}:concat'\n'
end

function Euler:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local code = initState.init(solver)	
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf,
	global real* ePotBuf
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
	
	//TODO make this B for Maxwell
	
	real3 B = _real3(0,0,0);	//set for MHD / thrown away for pure Euler
	real ePot = 0;

]]..code..[[

	<?=eqn.prim_t?> W = {.rho=rho, .v=v, .P=P};
	UBuf[index] = consFromPrim(W, ePot);
	ePotBuf[index] = ePot;
}
]], {
	eqn = self,
})
end

function Euler:getSolverCode(solver)	
	return template(file['eqn/euler.cl'], {eqn=self, solver=solver})
end

function Euler:getDisplayVarCodePrefix()
	return template([[
	<?=eqn.cons_t?> U = buf[index];
	real ePot = ePotBuf[index];
	<?=eqn.prim_t?> W = primFromCons(U, ePot);
]], {
	eqn = self,
})
end

function Euler:getDisplayVars(solver)
	return {
		{rho = 'value = W.rho;'},
		{vx = 'value = W.v.x;'},
		{vy = 'value = W.v.y;'},
		{vz = 'value = W.v.z;'},
		{v = 'value = coordLen(W.v);'},
		{mx = 'value = U.m.x;'},
		{my = 'value = U.m.y;'},
		{mz = 'value = U.m.z;'},
		{m = 'value = coordLen(U.m);'},
		{P = 'value = W.P;'},
		{eInt = 'value = calc_eInt(W);'},
		{eKin = 'value = calc_eKin(W);'},
		{ePot = 'value = ePot;'},
		{eTotal = 'value = U.ETotal / W.rho;'},
		{EInt = 'value = calc_EInt(W);'},
		{EKin = 'value = calc_EKin(W);'},
		{EPot = 'value = W.rho * ePot;'},
		{ETotal = 'value = U.ETotal;'},
		{S = 'value = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		{H = 'value = calc_H(W.P);'},
		{h = 'value = calc_h(W.rho, W.P);'},
		{HTotal = 'value = calc_HTotal(W.P, U.ETotal);'},
		{hTotal = 'value = calc_hTotal(W.rho, W.P, U.ETotal);'},
		{['Speed of Sound'] = 'value = calc_Cs(&W);'},
		{['Mach number'] = 'value = coordLen(W.v) / calc_Cs(&W);'},
	}
end

function Euler:getEigenTypeCode(solver)
	return [[
typedef struct {
	// Roe-averaged vars
	real rho;
	real3 v;
	real hTotal;

	// derived vars
	real vSq;
	real Cs;
} eigen_t;
]]
end

function Euler:getEigenDisplayVars(solver)
	return {
		{rho = 'value = eigen->rho;'},
		{vx = 'value = eigen->v.x;'},
		{vy = 'value = eigen->v.y;'},
		{vz = 'value = eigen->v.z;'},
		{v = 'value = coordLen(eigen->v);'},
		{hTotal = 'value = eigen->hTotal;'},
		{vSq = 'value = eigen->vSq;'},
		{Cs = 'value = eigen->Cs;'},
	}
end

return Euler
