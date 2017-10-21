--[[
Font 2008
similar to SRHD except using a metric based on a metric of alpha, beta, gamma
which needs to be provided externally from another solver (via gr-hd-separate-behavior)
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local GRHD = class(Equation)
GRHD.name = 'GRHD'
GRHD.numStates = 10
GRHD.numWaves = 5
GRHD.numIntStates = 5

GRHD.mirrorVars = {
	{'cons.S.x', 'prim.v.x'},
	{'cons.S.y', 'prim.v.y'},
	{'cons.S.z', 'prim.v.z'},
}

-- GRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--GRHD.hasFluxFromCons = true

GRHD.hasEigenCode = true 
GRHD.hasCalcDT = true
GRHD.useSourceTerm = true
GRHD.useConstrainU = true

GRHD.initStates = require 'init.euler'

function GRHD:init(...)
	GRHD.super.init(self, ...)
	self.cons_only_t = self:unique'cons_only_t'
end

function GRHD:getTypeCode()
	return template([[
typedef union {
	real ptr[5];
	struct {
		real D;		//0 D = rho W
		real3 S;	//1	S_j = rho h W^2 v_j
		real tau;	//4 tau = rho h W^2 - P
	};
} <?=eqn.cons_only_t?>;

typedef union {
	real ptr[5];
	struct {
		real rho;
		real3 v;	//v_i
		real eInt;
	};
} <?=eqn.prim_t?>;

typedef union {
	real ptr[10];
	struct {
		<?=eqn.cons_only_t?> cons;
		<?=eqn.prim_t?> prim;
	};
} <?=eqn.cons_t?>;
]], {
	eqn = self,
})
end

function GRHD:createInitState()
	GRHD.super.createInitState(self)
	
--[[ double precision
	self:addGuiVar{name='heatCapacityRatio', value=7/5}

	-- setting max iter to 100+ makes it freeze initially 
	-- but setting it to 100 after the first iteration is fine ...
	-- meaning the initial cons to prim is taking too long ...
	self:addGuiVar{name='solvePrimMaxIter', type='int', value=10}	-- value=1000}
	
	self:addGuiVar{name='solvePrimStopEpsilon', value=1e-7}
	
	-- used by pressure solver
	-- velocity epsilon is how close we can get to the speed of light
	-- set ylabel "Lorentz factor"; set xlabel "velocity epsilon -log10"; set log xy; plot [1:10] 1/sqrt(1-(1-10**(-x))**2);
	--velEpsilon = 1e-5	-- <=> handles up to W = 500
	--velEpsilon = 1e-6	-- <=> handles up to W = 600
	--velEpsilon = 1e-7	-- <=> handles up to W = 2,000
	--velEpsilon = 1e-10	-- <=> handles up to W = 100,000
	-- <=> smaller than 1e-15 gnuplot x11 terminal breaks down past W = 1e+7 ...
	self:addGuiVar{name='solvePrimVelEpsilon', value=1e-15}	
	
	self:addGuiVar{name='solvePrimPMinEpsilon', value=1e-16}
	
	self:addGuiVar{name='rhoMin', value=1e-15}
	self:addGuiVar{name='rhoMax', value=1e+20}
	self:addGuiVar{name='eIntMax', value=1e+20}
	self:addGuiVar{name='DMin', value=1e-15}
	self:addGuiVar{name='DMax', value=1e+20}
	self:addGuiVar{name='tauMin', value=1e-15}
	self:addGuiVar{name='tauMax', value=1e+20}
--]]
-- [[ single precision?
	self:addGuiVar{name='heatCapacityRatio', value=7/5}
	self:addGuiVar{name='solvePrimMaxIter', type='int', value=10}	-- value=1000}
	self:addGuiVar{name='solvePrimStopEpsilon', value=1e-7}
	self:addGuiVar{name='solvePrimVelEpsilon', value=1e-7}	
	self:addGuiVar{name='solvePrimPMinEpsilon', value=1e-7}
	self:addGuiVar{name='rhoMin', value=1e-7}
	self:addGuiVar{name='rhoMax', value=1e+20}
	self:addGuiVar{name='eIntMax', value=1e+20}
	self:addGuiVar{name='DMin', value=1e-7}
	self:addGuiVar{name='DMax', value=1e+20}
	self:addGuiVar{name='tauMin', value=1e-7}
	self:addGuiVar{name='tauMax', value=1e+20}
--]]
end

function GRHD:getCodePrefix()
	return table{
		GRHD.super.getCodePrefix(self),
		template([[

//pressure function for ideal gas
real calc_P(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * rho * eInt;
}	

//chi in most papers
real calc_dP_drho(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * eInt;
}

//kappa in most papers
real calc_dP_deInt(real rho, real eInt) {
	return (heatCapacityRatio - 1.) * rho;
}

real calc_eInt_from_P(real rho, real P) {
	return P / ((heatCapacityRatio - 1.) * rho);
}

real calc_h(real rho, real P, real eInt) {
	return 1. + eInt + P / rho;
}

<?=eqn.cons_only_t?> consFromPrim(
	<?=eqn.prim_t?> prim,
	real alpha,
	real3 beta,
	sym3 gamma
) {
	//2008 Font eqn 31 etc 
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	real3 vU = sym3_real3_mul(gammaU, prim.v);
	real vSq = real3_dot(prim.v, vU);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	
	real P = calc_P(prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 28-30:
	
	real D = prim.rho * W;
	real3 S = real3_scale(prim.v, prim.rho * h * WSq);
	real tau = prim.rho * h * WSq - P - D;

	return (<?=eqn.cons_only_t?>){
		.D = D,
		.S = S,
		.tau = tau,
	};
}
]], {
	eqn = self,
}),
	}:concat'\n'
end

GRHD.initStateCode = [[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf<?=
	solver:getADMArgs()?>
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
	//ignored:
	real3 B = _real3(0,0,0);

	<?=solver:getADMVarCode()?>

	<?=code?>
	
	real eInt = calc_eInt_from_P(rho, P);

	<?=eqn.prim_t?> prim = {
		.rho = rho,
		.v = v,
		.eInt = eInt,
	};
	UBuf[index] = (<?=eqn.cons_t?>){
		.prim = prim,
		.cons = consFromPrim(prim, alpha, beta, gamma),
	};
}
]]

function GRHD:getSolverCode()
	return template(file['eqn/grhd.cl'], {
		eqn = self,
		solver = self.solver,
	})
end

function GRHD:getDisplayVars()
	return {
		{D = '*value = U->cons.D;'},
		{S = '*valuevec = U->cons.S;', type='real3'},
		{['S weighted'] = template([[
	<?=solver:getADMVarCode()?>
	*value = real3_weightedLen(U->cons.S, gamma);
]], {solver=self.solver})},
		{tau = '*value = U->cons.tau;'},
		{['W based on D'] = '*value = U->cons.D / U->prim.rho;'},
		{['W based on v'] = template([[
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	*value = 1. / sqrt(1. - real3_weightedLenSq(U->prim.v, gammaU));
]], {solver=self.solver})},
		{['primitive reconstruction error'] = template([[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	<?=solver:getADMVarCode()?>
	<?=eqn.cons_only_t?> U2 = consFromPrim(U->prim, alpha, beta, gamma);
	*value = 0;
	for (int j = 0; j < numIntStates; ++j) {
		*value += fabs(U->cons.ptr[j] - U2.ptr[j]);
	}
]], {eqn=self, solver=self.solver})},
		{['W error'] = template([[
	real W1 = U->cons.D / U->prim.rho;
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	real W2 = 1. / sqrt(1. - real3_weightedLenSq(U->prim.v, gammaU));
	*value = fabs(W1 - W2);
]], {solver=self.solver})},

		{rho = '*value = U->prim.rho;'},
		
		-- TODO abstract the generators of real3 variables and add weighted norms automatically
		{v = '*valuevec = U->prim.v;', type='real3'},
		{['v weighted'] = template([[
	<?=solver:getADMVarCode()?>
	*value = real3_weightedLen(U->prim.v, gamma);
]], {solver=self.solver})},

		{eInt = '*value = U->prim.eInt;'},
		{P = '*value = calc_P(U->prim.rho, U->prim.eInt);'},
		{h = '*value = calc_h(U->prim.rho, calc_P(U->prim.rho, U->prim.eInt), U->prim.eInt);'},
	}
end

GRHD.eigenVars = {
	{rho = 'real'},
	{vL = 'real3'},
	{h = 'real'},
	{W = 'real'},
	{ATildeMinus = 'real'},
	{ATildePlus = 'real'},
	{VMinus = 'real'},
	{VPlus = 'real'},
	{CMinus = 'real'},
	{CPlus = 'real'},
	{Kappa = 'real'},
	-- hmm, I could just as easily re-average these ...
	{alpha = 'real'},
	{beta = 'real3'},
	{gamma = 'sym3'},
}

return GRHD
