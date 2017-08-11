local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local Equation = require 'eqn.eqn'

local MHD = class(Equation)

MHD.name = 'MHD'

MHD.numStates = 9
MHD.numWaves = 7

MHD.mirrorVars = {{'m.x', 'B.x'}, {'m.y', 'B.y'}, {'m.z', 'B.z'}}

MHD.hasEigenCode = true
MHD.hasFluxFromCons = true

-- hmm, we want init.euler and init.mhd here ...
MHD.initStates = require 'init.euler'

local GuiFloat = require 'guivar.float'

function MHD:init(solver)
	self.guiVars = table{
		GuiFloat{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
		GuiFloat{name='mu0', value=1},	-- this should be 4 pi for natural units, but I haven't verified that all mu0's are where they should be ...
	}
	MHD.super.init(self, solver)
	
	local NoDiv = require 'solver.nodiv'
	solver.ops:insert(NoDiv{solver=solver})
end

function MHD:getTypeCode()
	return template([[
typedef struct {
	real ptr[9];
	struct {
		real rho;
		real3 v;
		real P;
		real3 B;
		real BPot;
	};
} <?=eqn.prim_t?>;

typedef union {
	real ptr[9];
	struct {
		real rho;
		real3 m;
		real ETotal;
		real3 B;
		real BPot;
	};
} <?=eqn.cons_t?>;
]], {
	eqn = self,
})
end

function MHD:getCodePrefix()
	return table{
		MHD.super.getCodePrefix(self),
		template([[
inline real calc_eKin(<?=eqn.prim_t?> W) { return .5 * real3_lenSq(W.v); }
inline real calc_EKin(<?=eqn.prim_t?> W) { return W.rho * calc_eKin(W); }
inline real calc_EInt(<?=eqn.prim_t?> W) { return W.P / (heatCapacityRatio - 1.); }
inline real calc_eInt(<?=eqn.prim_t?> W) { return calc_EInt(W) / W.rho; }
inline real calc_EMag(<?=eqn.prim_t?> W) { return .5 * real3_lenSq(W.B); }
inline real calc_eMag(<?=eqn.prim_t?> W) { return calc_EMag(W) / W.rho; }
inline real calc_PMag(<?=eqn.prim_t?> W) { return .5 * real3_lenSq(W.B); }
inline real calc_EHydro(<?=eqn.prim_t?> W) { return calc_EKin(W) + calc_EInt(W); }
inline real calc_eHydro(<?=eqn.prim_t?> W) { return calc_EHydro(W) / W.rho; }
inline real calc_ETotal(<?=eqn.prim_t?> W) { return calc_EKin(W) + calc_EInt(W) + calc_EMag(W); }
inline real calc_eTotal(<?=eqn.prim_t?> W) { return calc_ETotal(W) / W.rho; }
inline real calc_H(real P) { return P * (heatCapacityRatio / (heatCapacityRatio - 1.)); }
inline real calc_h(real rho, real P) { return calc_H(P) / rho; }
inline real calc_HTotal(<?=eqn.prim_t?> W, real ETotal) { return W.P + calc_PMag(W) + ETotal; }
inline real calc_hTotal(<?=eqn.prim_t?> W, real ETotal) { return calc_HTotal(W, ETotal) / W.rho; }
inline real calc_Cs(<?=eqn.prim_t?> W) { return sqrt(heatCapacityRatio * W.P / W.rho); }

inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U) {
	<?=eqn.prim_t?> W;
	W.rho = U.rho;
	W.v = real3_scale(U.m, 1./U.rho);
	W.B = U.B;
	real vSq = real3_lenSq(W.v);
	real BSq = real3_lenSq(W.B);
	real EKin = .5 * U.rho * vSq;
	real EMag = .5 * BSq;
	real EInt = U.ETotal - EKin - EMag;
	W.P = EInt * (heatCapacityRatio - 1.);
	W.P = max(W.P, (real)1e-7);
	W.rho = max(W.rho, (real)1e-7);
	W.BPot = U.BPot;
	return W;
}

inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W) {
	<?=eqn.cons_t?> U;
	U.rho = W.rho;
	U.m = real3_scale(W.v, W.rho);
	U.B = W.B;
	real vSq = real3_lenSq(W.v);
	real BSq = real3_lenSq(W.B);
	real EKin = .5 * W.rho * vSq;
	real EMag = .5 * BSq;
	real EInt = W.P / (heatCapacityRatio - 1.);
	U.ETotal = EInt + EKin + EMag;
	U.BPot = W.BPot;
	return U;
}
]], {
	eqn = self,
}),
	}:concat'\n'
end

function MHD:getInitStateCode()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)	
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
	real3 B = _real3(0,0,0);

]]..code..[[
	
	<?=eqn.prim_t?> W = {.rho=rho, .v=v, .P=P, .B=B, .BPot=0};
	UBuf[index] = consFromPrim(W);
}
]], {
	eqn = self,
})
end

function MHD:getSolverCode()
	return template(file['eqn/mhd.cl'], {eqn=self, solver=self.solver})
end

function MHD:getDisplayVarCodePrefix()
	return template([[
	global const <?=eqn.cons_t?>* U = buf + index;
	<?=eqn.prim_t?> W = primFromCons(*U);
]], {
	eqn = self,
})
end

function MHD:getDisplayVars()
	return {
		{rho = '*value = W.rho;'},
		{vx = '*value = W.v.x;'},
		{vy = '*value = W.v.y;'},
		{vz = '*value = W.v.z;'},
		{['|v|'] = '*value = real3_len(W.v);'},
		{mx = '*value = U->m.x;'},
		{my = '*value = U->m.y;'},
		{mz = '*value = U->m.z;'},
		{['|m|'] = '*value = real3_len(U->m);'},
		{Bx = '*value = W.B.x;'},
		{By = '*value = W.B.y;'},
		{Bz = '*value = W.B.z;'},
		{['|B|'] = '*value = real3_len(W.B);'},
		{['div B'] = template([[
	*value = .5 * (0.
<? 
for j=0,solver.dim-1 do 
?>		+ (U[stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / grid_dx<?=j?>
<? 
end 
?>	)<? 
if field == 'epsE' then 
?> / eps0<?
end
?>;
]], {solver=self.solver, field='B'})},
		{['BPot'] = '*value = U->BPot;'},
		{P = '*value = W.P;'},
		--{PMag = '*value = calc_PMag(W);'},
		--{PTotal = '*value = W.P + calc_PMag(W);'},
		--{eInt = '*value = calc_eInt(W);'},
		{EInt = '*value = calc_EInt(W);'},
		--{eKin = '*value = calc_eKin(W);'},
		{EKin = '*value = calc_EKin(W);'},
		--{eHydro = '*value = calc_eHydro(W);'},
		{EHydro = '*value = calc_EHydro(W);'},
		--{eMag = '*value = calc_eMag(W);'},
		{EMag = '*value = calc_EMag(W);'},
		--{eTotal = '*value = U->ETotal / W.rho;'},
		{ETotal = '*value = U->ETotal;'},
		{S = '*value = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		{H = '*value = calc_H(W.P);'},
		--{h = '*value = calc_H(W.P) / W.rho;'},
		--{HTotal = '*value = calc_HTotal(W, U->ETotal);'},
		--{hTotal = '*value = calc_hTotal(W, U->ETotal);'},
		--{Cs = '*value = calc_Cs(W); },
		{['primitive reconstruction error'] = template([[
		//prim have just been reconstructed from cons
		//so reconstruct cons from prims again and calculate the difference
		<?=eqn.cons_t?> U2 = consFromPrim(W);
		*value = 0;
		for (int j = 0; j < numStates; ++j) {
			*value += fabs(U->ptr[j] - U2.ptr[j]);
		}
]], {
	eqn = self,
})},
	}
end

function MHD:getVecDisplayVars()
	local vars = table{
		{v = 'valuevec = W.v;'},
		{m = 'valuevec = U->m;'},
		{B = 'valuevec = U->B;'},
	}
	return vars
end

function MHD:getEigenTypeCode()
	return template([[
typedef struct {
	real evL[7*7];
	real evR[7*7];
<? if solver.checkFluxError then ?>
	real A[7*7];
<? end ?>
} <?=eqn.eigen_t?>;
]], {
		eqn = self,
		solver = self.solver,
	})
end

-- because eigen_t is only 7*7 instead of 7*8 = numStates * numWaves ...
function MHD:getEigenDisplayVars()
	return range(0, self.numWaves * self.numWaves - 1):map(function(i)
		local row = i%self.numWaves
		local col = (i-row)/self.numWaves
		return {['evL_'..row..'_'..col] = '*value = eigen->evL['..i..'];'}
	end):append(range(0, self.numWaves * self.numWaves - 1):map(function(i)
		local row = i % self.numWaves
		local col = (i-row)/self.numWaves
		return {['evR_'..row..'_'..col] = '*value = eigen->evR['..i..'];'}
	end)):append(self.solver.checkFluxError and range(0, self.numWaves * self.numWaves - 1):map(function(i)
		local row = i%self.numWaves
		local col = (i-row)/self.numWaves
		return {['A_'..row..'_'..col] = '*value = eigen->A['..i..'];'}
	end) or nil)
end


return MHD
