--[[
2008 Font "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity"
similar to SRHD except using a metric based on a metric of alpha, beta, gamma
which needs to be provided externally from another solver (via gr-hd-separate-behavior)
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'
local Struct = require 'hydro.code.struct'

local GRHD = class(Equation)
GRHD.name = 'GRHD'

--GRHD.numStates = 10

GRHD.numWaves = 5
GRHD.numIntStates = 5

--GRHD.roeUseFluxFromCons = true

GRHD.initConds = require 'hydro.init.euler':getList()

function GRHD:init(args)
	local solver = assert(args.solver)

	self.consOnlyStruct = Struct{
		solver = solver,
		name = 'cons_only_t',
		vars = {
			{name='D', type='real', units='kg/m^3'},				-- D = ρ W, W = unitless Lorentz factor
			{name='S', type='real3', units='kg/s^3', variance='l'},	-- S_j = ρ h W^2 v_j ... [ρ] [h] [v] = kg/m^3 * m^2/s^2 * m/s = kg/s^3
			{name='tau', type='real', units='kg/(m*s^2)'},			-- tau = ρ h W^2 - P ... [ρ] [h] [W^2] = kg/m^3 * m^2/s^2 = kg/(m*s^2)
		},
	}

	self.primOnlyStruct = Struct{
		solver = solver,
		name = 'prim_only_t',
		vars = {
			{name='rho', type='real', units='kg/m^3'},
			{name='v', type='real3', units='m/s', variance='l'},
			{name='eInt', type='real', units='m^2/s^2'},
		},
	}

	-- TODO how about anonymous structs, so we can copy out prim_only_t and cons_only_t?
	self.consVars = table()
	:append(self.consOnlyStruct.vars)
	:append(self.primOnlyStruct.vars)


	self.consOnlyStruct:makeType()
	self.primOnlyStruct:makeType()

	GRHD.super.init(self, args)

	self.symbols.cons_only_t = self.consOnlyStruct.typename
	self.symbols.prim_only_t = self.primOnlyStruct.typename
end

function GRHD:getSymbolFields()
	return GRHD.super.getSymbolFields(self):append{
		'cons_only_t',
		'prim_only_t',
	}
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

function GRHD:initCodeModules()
	GRHD.super.initCodeModules(self)
	local solver = self.solver

	solver.modules:add{
		name = self.symbols.cons_only_t,
		structs = {self.consOnlyStruct},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.cons_only_t..' cons_only_t;',
	}

	solver.modules:add{
		name = self.symbols.prim_only_t,
		structs = {self.primOnlyStruct},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.prim_only_t..' prim_only_t;',
	}
end

-- don't use default
function GRHD:initCodeModule_fluxFromCons() end
function GRHD:initCodeModule_calcDTCell() end

GRHD.solverCodeFile = 'hydro/eqn/grhd.cl'

GRHD.predefinedDisplayVars = {
	'U rho',
	'U v',
	'U tau',
	'U div v',
	'U curl v',
}

function GRHD:getDisplayVars()
	local vars = GRHD.super.getDisplayVars(self)
	vars:append{
		{name='W based on D', code='value.vreal = U->D / U->rho;'},
		{name='W based on v', code=self:template[[
<?=solver:getADMVarCode()?>
real det_gamma = gamma.determinant();
real3s3 gammaU = gamma.inverse(det_gamma);
value.vreal = 1. / sqrt(1. - real3_weightedLenSq(U->v, gammaU));
]]},
		{name='primitive reconstruction error', code=self:template[[
//prim have just been reconstructed from cons
//so reconstruct cons from prims again and calculate the difference
{
	<?=solver:getADMVarCode()?>
	<?=cons_only_t?> U2;
	consOnlyFromPrim(&U2, solver, U, alpha, beta, gamma);
	value.vreal = 0;
	for (int j = 0; j < numIntStates; ++j) {
		value.vreal += fabs(U->ptr[j] - U2.ptr[j]);
	}
}
]]
		},
		{name='W error', code=self:template[[
real W1 = U->D / U->rho;
<?=solver:getADMVarCode()?>
real det_gamma = gamma.determinant();
real3s3 gammaU = gamma.inverse(det_gamma);
real W2 = 1. / sqrt(1. - real3_weightedLenSq(U->v, gammaU));
value.vreal = fabs(W1 - W2);
]]		},

		{name='v weighted', code=self:template[[
<?=solver:getADMVarCode()?>
value.vreal = real3_weightedLen(U->v, gamma);
]]},
		-- TODO weighted metric norm option for all tensors/vectors
		{name='S weighted', code=self:template[[
<?=solver:getADMVarCode()?>
value.vreal = real3_weightedLen(U->S, gamma);
]]},

		{name='eInt', code='value.vreal = U->eInt;'},
		{name='P', code='value.vreal = calc_P(solver, U->rho, U->eInt);'},
		{name='h', code='value.vreal = calc_h(U->rho, calc_P(solver, U->rho, U->eInt), U->eInt);'},
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
	{name='gamma', type='real3s3'},
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

GRHD.consWaveCode = GRHD.eigenWaveCode

return GRHD
