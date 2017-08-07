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
GRHD.numIntVars = 5
GRHD.numWaves = 5
GRHD.numStates = 5

GRHD.mirrorVars = {{'S.x'}, {'S.y'}, {'S.z'}}

GRHD.hasEigenCode = true 

-- GRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--GRHD.hasFluxFromCons = true

GRHD.hasCalcDT = true
GRHD.useSourceTerm = true
GRHD.useConstrainU = true

GRHD.initStates = require 'init.euler'

local GuiFloat = require 'guivar.float'
local GuiInt = require 'guivar.int'
function GRHD:init(...)
	self.guiVars = table{
--[[ double precision
		GuiFloat{name='heatCapacityRatio', value=7/5},

		-- setting max iter to 100+ makes it freeze initially 
		-- but setting it to 100 after the first iteration is fine ...
		-- meaning the initial cons to prim is taking too long ...
		GuiInt{name='solvePrimMaxIter', value=10},	-- value=1000},
		
		GuiFloat{name='solvePrimStopEpsilon', value=1e-7},
		
		-- used by pressure solver
		-- velocity epsilon is how close we can get to the speed of light
		-- set ylabel "Lorentz factor"; set xlabel "velocity epsilon -log10"; set log xy; plot [1:10] 1/sqrt(1-(1-10**(-x))**2);
		--velEpsilon = 1e-5	-- <=> handles up to W = 500
		--velEpsilon = 1e-6	-- <=> handles up to W = 600
		--velEpsilon = 1e-7	-- <=> handles up to W = 2,000
		--velEpsilon = 1e-10	-- <=> handles up to W = 100,000
		-- <=> smaller than 1e-15 gnuplot x11 terminal breaks down past W = 1e+7 ...
		GuiFloat{name='solvePrimVelEpsilon', value=1e-15},	
		
		GuiFloat{name='solvePrimPMinEpsilon', value=1e-16},
		
		GuiFloat{name='rhoMin', value=1e-15},
		GuiFloat{name='rhoMax', value=1e+20},
		GuiFloat{name='eIntMax', value=1e+20},
		GuiFloat{name='DMin', value=1e-15},
		GuiFloat{name='DMax', value=1e+20},
		GuiFloat{name='tauMin', value=1e-15},
		GuiFloat{name='tauMax', value=1e+20},
--]]
-- [[ single precision?
		GuiFloat{name='heatCapacityRatio', value=7/5},
		GuiInt{name='solvePrimMaxIter', value=10},	-- value=1000},
		GuiFloat{name='solvePrimStopEpsilon', value=1e-7},
		GuiFloat{name='solvePrimVelEpsilon', value=1e-7},	
		GuiFloat{name='solvePrimPMinEpsilon', value=1e-7},
		GuiFloat{name='rhoMin', value=1e-7},
		GuiFloat{name='rhoMax', value=1e+20},
		GuiFloat{name='eIntMax', value=1e+20},
		GuiFloat{name='DMin', value=1e-7},
		GuiFloat{name='DMax', value=1e+20},
		GuiFloat{name='tauMin', value=1e-7},
		GuiFloat{name='tauMax', value=1e+20},
--]]
	}
	GRHD.super.init(self, ...)
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
} <?=eqn.cons_t?>;

typedef struct {
	real rho;
	real3 v;	//v_i
	real eInt;
} <?=eqn.prim_t?>;
]], {
	eqn = self,
})
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

<?=eqn.cons_t?> consFromPrim(
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

	return (<?=eqn.cons_t?>){
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

function GRHD:getInitStateCode()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)
	local code = initState.init(self.solver)
	return template([[

kernel void initState(
	global <?=eqn.cons_t?>* consBuf,
	global <?=eqn.prim_t?>* primBuf<?=
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

]]..code..[[
	
	real eInt = calc_eInt_from_P(rho, P);

	<?=eqn.prim_t?> prim = {
		.rho = rho,
		.v = v,
		.eInt = eInt,
	};
	primBuf[index] = prim;	
	consBuf[index] = consFromPrim(prim, alpha, beta, gamma);
}
]], {
	eqn = self,
	solver = self.solver,
})
end

function GRHD:getSolverCode()
	return template(file['eqn/grhd.cl'], {
		eqn = self,
		solver = self.solver,
	})
end

function GRHD:getDisplayVarCodePrefix()
	return template([[
	<?=eqn.cons_t?> U = buf[index];
	<?=eqn.prim_t?> prim = primBuf[index];
]], {
	eqn = self,
})
end

function GRHD:getDisplayVars()
	return {
		{D = '*value = U.D;'},
		{S_x = '*value = U.S.x;'},
		{S_y = '*value = U.S.y;'},
		{S_z = '*value = U.S.z;'},
		{S = template([[
	<?=solver:getADMVarCode()?>
	*value = real3_weightedLen(U.S, gamma);
]], {solver=self.solver})},
		{tau = '*value = U.tau;'},
		{['W based on D'] = '*value = U.D / prim.rho;'},
		{['W based on v'] = template([[
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	*value = 1. / sqrt(1. - real3_weightedLenSq(prim.v, gammaU));
]], {solver=self.solver})},
		{['primitive reconstruction error'] = template([[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	<?=solver:getADMVarCode()?>
	<?=eqn.cons_t?> U2 = consFromPrim(prim, alpha, beta, gamma);
	*value = 0;
	for (int j = 0; j < numStates; ++j) {
		*value += fabs(U.ptr[j] - U2.ptr[j]);
	}
]], {eqn=self, solver=self.solver})},
		{['W error'] = template([[
	real W1 = U.D / prim.rho;
	<?=solver:getADMVarCode()?>
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(gamma, det_gamma);
	real W2 = 1. / sqrt(1. - real3_weightedLenSq(prim.v, gammaU));
	*value = fabs(W1 - W2);
]], {solver=self.solver})},
	}
end

function GRHD:getPrimDisplayVarCodePrefix()
	return template([[
	<?=eqn.prim_t?> prim = buf[index];
]], {
		eqn = self,
	})
end

GRHD.primDisplayVars = {
	{rho = '*value = prim.rho;'},
	{['v_x'] = '*value = prim.v.x;'},
	{['v_y'] = '*value = prim.v.y;'},
	{['v_z'] = '*value = prim.v.z;'},

	--TODO in gr-hd-separate, override gr's prim tex and give it an 'extraArgs'
	--{v = '*value = real3_weightedLen(prim.v, gamma);'},
	{['|v|'] = '*value = real3_len(prim.v);'},
	
	{eInt = '*value = prim.eInt;'},
	{P = '*value = calc_P(prim.rho, prim.eInt);'},
	{h = '*value = calc_h(prim.rho, calc_P(prim.rho, prim.eInt), prim.eInt);'},
}

GRHD.eigenStructFields = {
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

function GRHD:getEigenTypeCode()
	return 'typedef struct {\n'
		..table.map(self.eigenStructFields, function(field)
			local name, ctype = next(field)
			return '\t'..ctype..' '..name..';\n'
		end):concat'\n'
		..'} '..self.eigen_t..';\n'
end

function GRHD:getEigenDisplayVars()
	return {}
end

return GRHD
