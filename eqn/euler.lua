local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local Euler = class(Equation)
Euler.name = 'Euler'

Euler.numWaves = 5

-- ePot is the 6th param
-- which means it's now in the derivBuf, but it is always zero
-- so TODO a new variable for deriv size vs cons_t size?
Euler.numStates = 6	

Euler.mirrorVars = {{'m.x'}, {'m.y'}, {'m.z'}} 

Euler.hasEigenCode = true

-- this is working in the gravitation-waves simulation, but not here ...
--Euler.hasFluxFromCons = true

Euler.initStates = require 'init.euler'

Euler.guiVars = {
	require 'guivar.float'{name='heatCapacityRatio', 
--in order to make things work, gamma needs to be set *HERE AND IN INIT/EULER*
-- which means it is being read and written in multiple places
-- TODO consolidate that
		value=7/5,
	}
}

function Euler:getTypeCode()
	return template([[
typedef union {
	real ptr[6];
	struct { 
		real rho;
		real3 v;
		real P;
		real ePot;
	};
} <?=eqn.prim_t?>;

typedef union {
	real ptr[6];
	struct {
		real rho;
		real3 m;
		real ETotal;
		real ePot;
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
inline real calc_ETotal(<?=eqn.prim_t?> W) {
	real EPot = W.rho * W.ePot;
	return calc_EKin(W) + calc_EInt(W) + EPot;
}

inline real calc_Cs(const <?=eqn.prim_t?>* W) {
	return sqrt(heatCapacityRatio * W->P / W->rho);
}

inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U) {
	real EPot = U.rho * U.ePot;
	real EKin = calc_EKin_fromCons(U);
	real EInt = U.ETotal - EPot - EKin;
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_scale(U.m, 1./U.rho),
		.P = (heatCapacityRatio - 1.) * EInt,
		.ePot = U.ePot,
	};
}

<?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W) {
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_scale(W.v, W.rho),
		.ETotal = calc_ETotal(W),
		.ePot = W.ePot,
	};
}
]], 	{
			eqn = self,
		})
	}:concat'\n'
end

function Euler:getInitStateCode()
	local initState = self.initStates[1+self.solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..self.solver.initStatePtr[0])	
	local code = initState.init(self.solver)	
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
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

	//v's are vectors which need to be transformed from cartesian to our geometry

	<?=eqn.prim_t?> W = {
		.rho = rho,
		.v = cartesianToGrid(v, x),
		.P = P,
		.ePot = ePot,
	};
	UBuf[index] = consFromPrim(W);
}
]], {
		eqn = self,
	})
end

function Euler:getSolverCode()
	return template(file['eqn/euler.cl'], {eqn=self, solver=self.solver})
end

function Euler:getDisplayVarCodePrefix()
	return template([[
	<?=eqn.cons_t?> U = buf[index];
	<?=eqn.prim_t?> W = primFromCons(U);
]], {
	eqn = self,
})
end

function Euler:getDisplayVars()
	-- k is 0,1,2
	local function vorticity(k)
		local xs = {'x','y','z'}
		local i = (k+1)%3
		local j = (i+1)%3
		return {['vorticity '..xs[k+1]] = template([[
	global const <?=eqn.cons_t?>* Uim = buf + index - stepsize.s<?=i?>;
	global const <?=eqn.cons_t?>* Uip = buf + index + stepsize.s<?=i?>;
	global const <?=eqn.cons_t?>* Ujm = buf + index - stepsize.s<?=j?>;
	global const <?=eqn.cons_t?>* Ujp = buf + index + stepsize.s<?=j?>;

	//TODO incorporate metric

	real3 vim = real3_scale(Uim->m, 1. / Uim->rho);
	real3 vip = real3_scale(Uip->m, 1. / Uip->rho);
	real3 vjm = real3_scale(Ujm->m, 1. / Ujm->rho);
	real3 vjp = real3_scale(Ujp->m, 1. / Ujp->rho);
	
	value = (vjp.s<?=i?> - vjm.s<?=i?>) / (2. * grid_dx<?=i?>)
			- (vip.s<?=j?> - vim.s<?=j?>) / (2. * grid_dx<?=j?>);
]], 	{
			i = i,
			j = j,
			eqn = self,
		})}
	end
	return table{
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
		{ePot = 'value = U.ePot;'},
		{eTotal = 'value = U.ETotal / W.rho;'},
		{EInt = 'value = calc_EInt(W);'},
		{EKin = 'value = calc_EKin(W);'},
		{EPot = 'value = U.rho * U.ePot;'},
		{ETotal = 'value = U.ETotal;'},
		{S = 'value = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		{H = 'value = calc_H(W.P);'},
		{h = 'value = calc_h(W.rho, W.P);'},
		{HTotal = 'value = calc_HTotal(W.P, U.ETotal);'},
		{hTotal = 'value = calc_hTotal(W.rho, W.P, U.ETotal);'},
		{['Speed of Sound'] = 'value = calc_Cs(&W);'},
		{['Mach number'] = 'value = coordLen(W.v) / calc_Cs(&W);'},
	}:append( ({
	-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
			[1] = {},
			[2] = {vorticity(2)},
			[3] = range(0,2):map(vorticity),

	})[self.solver.dim] )
end

function Euler:getEigenTypeCode()
	return template([[
typedef struct {
	// Roe-averaged vars
	real rho;
	real3 v;
	real hTotal;

	// derived vars
	real vSq;
	real Cs;
} <?=eqn.eigen_t?>;
]], {
	eqn = self,
})
end

function Euler:getEigenDisplayVars()
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
