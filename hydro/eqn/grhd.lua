--[[
Font 2008
similar to SRHD except using a metric based on a metric of alpha, beta, gamma
which needs to be provided externally from another solver (via gr-hd-separate-behavior)
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local clnumber = require 'cl.obj.number'
local Equation = require 'hydro.eqn.eqn'

local GRHD = class(Equation)
GRHD.name = 'GRHD'
GRHD.numStates = 10
GRHD.numWaves = 5
GRHD.numIntStates = 5

GRHD.reflectVars = {
	mirror = {
		{'cons.S.x', 'prim.v.x'},
		{'cons.S.y', 'prim.v.y'},
		{'cons.S.z', 'prim.v.z'},
	},
}

GRHD.hasCalcDTCode = true
GRHD.useConstrainU = true

-- GRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--GRHD.roeUseFluxFromCons = true
GRHD.hasFluxFromConsCode = true
GRHD.useSourceTerm = true

GRHD.initConds = require 'hydro.init.euler'

function GRHD:init(args)
	GRHD.super.init(self, args)
	self.cons_only_t = args.solver.app:uniqueName'cons_only_t'
end

-- TODO upgrade this to srhd: put these all in consVars and just make separate cons_only_ and prim_only_t
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

-- hmm, this is from renovating 'getPrimConsCode() end', but will it work with the module system?
function GRHD:initCodeModulePrimCons() end

GRHD.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
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

GRHD.solverCodeFile = 'hydro/eqn/grhd.cl'

function GRHD:getDisplayVars()
	return {
		{name='D', code='value.vreal = U->cons.D;'},
		{name='S', code='value.vreal3 = U->cons.S;', type='real3'},
		{name='S weighted', code=template([[
	<?=solver:getADMVarCode()?>
	value.vreal = real3_weightedLen(U->cons.S, gamma);
]], {solver=self.solver})},
		{name='tau', code='value.vreal = U->cons.tau;'},
		{name='W based on D', code='value.vreal = U->cons.D / U->prim.rho;'},
		{name='W based on v', code=template([[
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	value.vreal = 1. / sqrt(1. - real3_weightedLenSq(U->prim.v, gammaU));
]], {solver=self.solver})},
		{name='primitive reconstruction error', code=template([[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	<?=solver:getADMVarCode()?>
	<?=eqn.cons_only_t?> U2 = consOnlyFromPrim(solver, U->prim, alpha, beta, gamma);
	value.vreal = 0;
	for (int j = 0; j < numIntStates; ++j) {
		value.vreal += fabs(U->cons.ptr[j] - U2.ptr[j]);
	}
]], {eqn=self, solver=self.solver})},
		{name='W error', code=template([[
	real W1 = U->cons.D / U->prim.rho;
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	real W2 = 1. / sqrt(1. - real3_weightedLenSq(U->prim.v, gammaU));
	value.vreal = fabs(W1 - W2);
]], {solver=self.solver})},

		{name='rho', code='value.vreal = U->prim.rho;'},
		
		-- TODO abstract the generators of real3 variables and add weighted norms automatically
		{name='v', code='value.vreal3 = U->prim.v;', type='real3'},
		{name='v weighted', code=template([[
	<?=solver:getADMVarCode()?>
	value.vreal = real3_weightedLen(U->prim.v, gamma);
]], {solver=self.solver})},

		{name='eInt', code='value.vreal = U->prim.eInt;'},
		{name='P', code='value.vreal = calc_P(solver, U->prim.rho, U->prim.eInt);'},
		{name='h', code='value.vreal = calc_h(U->prim.rho, calc_P(solver, U->prim.rho, U->prim.eInt), U->prim.eInt);'},
	}
	
	vars:insert(self:createDivDisplayVar{field='v', units='1/s'} or nil)
	vars:insert(self:createCurlDisplayVar{field='v', units='1/s'} or nil)
end

GRHD.eigenVars = {
	{name='rho', type='real'},
	{name='vL', type='real3'},
	{name='h', type='real'},
	{name='W', type='real'},
	{name='ATildeMinus', type='real'},
	{name='ATildePlus', type='real'},
	{name='VMinus', type='real'},
	{name='VPlus', type='real'},
	{name='CMinus', type='real'},
	{name='CPlus', type='real'},
	{name='Kappa', type='real'},
	-- hmm, I could just as easily re-average these ...
	{name='alpha', type='real'},
	{name='beta', type='real3'},
	{name='gamma', type='sym3'},
	-- and these are for wavespeeds
	{name='vU', type='real3'},
	{name='lambdaMin', type='real'},
	{name='lambdaMax', type='real'},
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
