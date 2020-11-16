--[[
2008 Font "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity"
similar to SRHD except using a metric based on a metric of alpha, beta, gamma
which needs to be provided externally from another solver (via gr-hd-separate-behavior)
--]]
local class = require 'ext.class'
local Equation = require 'hydro.eqn.eqn'

local GRHD = class(Equation)
GRHD.name = 'GRHD'
GRHD.numStates = 10
GRHD.numWaves = 5
GRHD.numIntStates = 5

GRHD.useConstrainU = true

--GRHD.roeUseFluxFromCons = true

GRHD.initConds = require 'hydro.init.euler':getList()

function GRHD:init(args)
	GRHD.super.init(self, args)
	self.cons_only_t = args.solver.app:uniqueName'cons_only_t'
end

-- TODO upgrade this to srhd: put these all in consVars and just make separate cons_only_ and prim_only_t
-- TODO also upgrade this to initCodeModules.  turn it into a struct, like srhd.
function GRHD:getTypeCode()
	return self:template[[
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
]]
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

-- don't use default 
function GRHD:initCodeModulePrimCons() end

GRHD.solverCodeFile = 'hydro/eqn/grhd.cl'

function GRHD:getDisplayVars()
	local vars = table{
		{name='D', code='value.vreal = U->cons.D;'},
		{name='S', code='value.vreal3 = U->cons.S;', type='real3'},
		{name='S weighted', code=self:template[[
	<?=solver:getADMVarCode()?>
	value.vreal = real3_weightedLen(U->cons.S, gamma);
]]},
		{name='tau', code='value.vreal = U->cons.tau;'},
		{name='W based on D', code='value.vreal = U->cons.D / U->prim.rho;'},
		{name='W based on v', code=self:template[[
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	value.vreal = 1. / sqrt(1. - real3_weightedLenSq(U->prim.v, gammaU));
]]},
		{name='primitive reconstruction error', code=self:template[[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	<?=solver:getADMVarCode()?>
	<?=eqn.cons_only_t?> U2 = consOnlyFromPrim(solver, U->prim, alpha, beta, gamma);
	value.vreal = 0;
	for (int j = 0; j < numIntStates; ++j) {
		value.vreal += fabs(U->cons.ptr[j] - U2.ptr[j]);
	}
]]},
		{name='W error', code=self:template[[
	real W1 = U->cons.D / U->prim.rho;
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	real W2 = 1. / sqrt(1. - real3_weightedLenSq(U->prim.v, gammaU));
	value.vreal = fabs(W1 - W2);
]]},

		{name='rho', code='value.vreal = U->prim.rho;'},
		
		-- TODO abstract the generators of real3 variables and add weighted norms automatically
		{name='v', code='value.vreal3 = U->prim.v;', type='real3'},
		{name='v weighted', code=self:template[[
	<?=solver:getADMVarCode()?>
	value.vreal = real3_weightedLen(U->prim.v, gamma);
]]},

		{name='eInt', code='value.vreal = U->prim.eInt;'},
		{name='P', code='value.vreal = calc_P(solver, U->prim.rho, U->prim.eInt);'},
		{name='h', code='value.vreal = calc_h(U->prim.rho, calc_P(solver, U->prim.rho, U->prim.eInt), U->prim.eInt);'},
	}
	
	vars:insert(self:createDivDisplayVar{field='v', units='1/s'} or nil)
	vars:insert(self:createCurlDisplayVar{field='v', units='1/s'} or nil)

	return vars
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
	return self:template(assert(({
		'<?=eig?>.lambdaMin',
		'<?=eig?>.vU.x * <?=eig?>.alpha - <?=eig?>.beta.s<?=side?>',
		'<?=eig?>.vU.x * <?=eig?>.alpha - <?=eig?>.beta.s<?=side?>',
		'<?=eig?>.vU.x * <?=eig?>.alpha - <?=eig?>.beta.s<?=side?>',
		'<?=eig?>.lambdaMax',
	})[waveIndex+1], "couldn't find code for waveIndex="..waveIndex), {side=side, eig='('..eig..')'})
end

return GRHD
