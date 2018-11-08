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


GRHD.hasEigenCode = true 
GRHD.hasCalcDTCode = true
GRHD.useConstrainU = true

-- GRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--GRHD.roeUseFluxFromCons = true
GRHD.hasFluxFromConsCode = true
GRHD.useSourceTerm = true

GRHD.initStates = require 'init.euler'

function GRHD:init(args)
	GRHD.super.init(self, args)
	self.cons_only_t = args.solver.app:uniqueName'cons_only_t'
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

	-- hmm, setting the higher thresholds using double precision is less stable
	local double = false --solver.app.real == 'double'

	self:addGuiVars{
		{name='heatCapacityRatio', value=7/5},
		
		-- setting max iter to 100+ makes it freeze initially 
		-- but setting it to 100 after the first iteration is fine ...
		-- meaning the initial cons to prim is taking too long ...
		{name='solvePrimMaxIter', type='int', value=10, compileTime=true},	-- value=1000},

		{name='solvePrimStopEpsilon', value=1e-7},

		-- used by pressure solver
		-- velocity epsilon is how close we can get to the speed of light
		-- set ylabel "Lorentz factor"; set xlabel "velocity epsilon -log10"; set log xy; plot [1:10] 1/sqrt(1-(1-10**(-x))**2);
		--velEpsilon = 1e-5	-- <=> handles up to W = 500
		--velEpsilon = 1e-6	-- <=> handles up to W = 600
		--velEpsilon = 1e-7	-- <=> handles up to W = 2,000
		--velEpsilon = 1e-10	-- <=> handles up to W = 100,000
		-- <=> smaller than 1e-15 gnuplot x11 terminal breaks down past W = 1e+7 ...
		{name='solvePrimVelEpsilon', value=double and 1e-15 or 1e-7},
		
		{name='solvePrimPMinEpsilon', value=double and 1e-16 or 1e-7},
		
		{name='rhoMin', value=double and 1e-15 or 1e-7},
		{name='rhoMax', value=1e+20},
		{name='eIntMax', value=1e+20},
		{name='DMin', value=double and 1e-15 or 1e-7},
		{name='DMax', value=1e+20},
		{name='tauMin', value=double and 1e-15 or 1e-7},
		{name='tauMax', value=1e+20},
	}
end

function GRHD:getCommonFuncCode()
	return template([[

//pressure function for ideal gas
real calc_P(constant <?=solver.solver_t?>* solver, real rho, real eInt) {
	return (solver->heatCapacityRatio - 1.) * rho * eInt;
}	

//chi in most papers
real calc_dP_drho(constant <?=solver.solver_t?>* solver, real rho, real eInt) {
	return (solver->heatCapacityRatio - 1.) * eInt;
}

//kappa in most papers
real calc_dP_deInt(constant <?=solver.solver_t?>* solver, real rho, real eInt) {
	return (solver->heatCapacityRatio - 1.) * rho;
}

real calc_eInt_from_P(constant <?=solver.solver_t?>* solver, real rho, real P) {
	return P / ((solver->heatCapacityRatio - 1.) * rho);
}

real calc_h(real rho, real P, real eInt) {
	return 1. + eInt + P / rho;
}

<?=eqn.cons_only_t?> consOnlyFromPrim(
	constant <?=solver.solver_t?>* solver,
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
	real P = calc_P(solver, prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 28-30:
	
	real D = prim.rho * W;
	real3 S = real3_real_mul(prim.v, prim.rho * h * WSq);
	real tau = prim.rho * h * WSq - P - D;

	return (<?=eqn.cons_only_t?>){.D=D, .S=S, .tau=tau};
}
]], {
		eqn = self,
		solver = self.solver,
	})
end

function GRHD:getPrimConsCode() end

GRHD.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf<?=
	solver:getADMArgs()?>
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	//ignored:
	real3 B = real3_zero;

	<?=solver:getADMVarCode()?>

	<?=code?>
	
	real eInt = calc_eInt_from_P(solver, rho, P);

	<?=eqn.prim_t?> prim = {.rho=rho, .v=v, .eInt=eInt};
	UBuf[index] = (<?=eqn.cons_t?>){
		.prim = prim,
		.cons = consOnlyFromPrim(solver, prim, alpha, beta, gamma),
	};
}
]]

GRHD.solverCodeFile = 'eqn/grhd.cl'

function GRHD:getCalcEigenBasisCode() end	-- within grhd.cl

function GRHD:getDisplayVars()
	return {
		{D = '*value = U->cons.D;'},
		{S = '*value_real3 = U->cons.S;', type='real3'},
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
	<?=eqn.cons_only_t?> U2 = consOnlyFromPrim(solver, U->prim, alpha, beta, gamma);
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
		{v = '*value_real3 = U->prim.v;', type='real3'},
		{['v weighted'] = template([[
	<?=solver:getADMVarCode()?>
	*value = real3_weightedLen(U->prim.v, gamma);
]], {solver=self.solver})},

		{eInt = '*value = U->prim.eInt;'},
		{P = '*value = calc_P(solver, U->prim.rho, U->prim.eInt);'},
		{h = '*value = calc_h(U->prim.rho, calc_P(solver, U->prim.rho, U->prim.eInt), U->prim.eInt);'},
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
	-- and these are for wavespeeds
	{vU = 'real3'},
	{lambdaMin = 'real'},
	{lambdaMax = 'real'},
}

function GRHD:eigenWaveCode(side, eig, x, waveIndex)
	return template(assert(({
		'<?=eig?>.lambdaMin',
		'<?=eig?>.vU.x * <?=eig?>.alpha - <?=eig?>.beta.s<?=side?>',
		'<?=eig?>.vU.x * <?=eig?>.alpha - <?=eig?>.beta.s<?=side?>',
		'<?=eig?>.vU.x * <?=eig?>.alpha - <?=eig?>.beta.s<?=side?>',
		'<?=eig?>.lambdaMax',
	})[waveIndex+1], "couldn't find code for waveIndex="..waveIndex), {side=side, eig='('..eig..')'})
end

return GRHD
