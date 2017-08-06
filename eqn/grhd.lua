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
GRHD.numStates = 6

GRHD.mirrorVars = {{'S.x'}, {'S.y'}, {'S.z'}}

GRHD.hasEigenCode = true 

-- GRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--GRHD.hasFluxFromCons = true

GRHD.hasCalcDT = true

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
		real D;
		real3 S;
		real tau;
	
		//metric variables
		//populated by whatever NR solver is used
		real alpha;
		real3 beta;
		sym3 gamma;
	};
} <?=eqn.cons_t?>;

typedef struct {
	real rho;
	real3 v;
	real eInt;

	real alpha;
	real3 beta;
	sym3 gamma;
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

<?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> prim, real3 x) {
	real vSq = coordLenSq(prim.v, x);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);
	real alpha = prim.alpha;

	//2008 Font, eqn 40-42:
	
	//rest-mass density = J^0 = rho u^0
	real D = prim.rho * W;	
	
	//momentum = T^0i = rho h u^0 u^i + P g^0i
	real3 S = real3_add(
		real3_scale(prim.v, prim.rho * h * WSq),
		real3_scale(prim.beta, P / (alpha * alpha)));
	
	//energy = T^00 = rho h u^0 u^0 + P g^00
	real tau = prim.rho * h * WSq - D - P / (alpha * alpha);

	return (<?=eqn.cons_t?>){
		.D = D,
		.S = S,
		.tau = tau,
		.alpha = prim.alpha,
		.beta = prim.beta,
		.gamma = prim.gamma,
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
	global <?=eqn.prim_t?>* primBuf
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

	real alpha = 1.;
	real3 beta = _real3(0,0,0);
	sym3 gamma = _sym3(1,0,0,1,0,1);

]]..code..[[
	
	real eInt = calc_eInt_from_P(rho, P);
	real vSq = coordLenSq(v, x);
	real W = 1. / sqrt(1. - vSq);
	real h = calc_h(rho, P, eInt);

	<?=eqn.prim_t?> prim = {
		.rho = rho,
		.v = v,
		.eInt = eInt,
		.alpha = alpha,
		.beta = beta,
		.gamma = gamma,
	};
	primBuf[index] = prim;
	consBuf[index] = consFromPrim(prim, x);
}
]], {
	eqn = self,
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
		{Sx = '*value = U.S.x;'},
		{Sy = '*value = U.S.y;'},
		{Sz = '*value = U.S.z;'},
		{S = '*value = coordLen(U.S, x);'},
		{tau = '*value = U.tau;'},
		{W = '*value = U.D / prim.rho;'},
		{['primitive reconstruction error'] = template([[
			//prim have just been reconstructed from cons
			//so reconstruct cons from prims again and calculate the difference
			{
				<?=eqn.cons_t?> U2 = consFromPrim(prim, x);
				*value = 0;
				for (int j = 0; j < numStates; ++j) {
					*value += fabs(U.ptr[j] - U2.ptr[j]);
				}
			}
	]], {eqn=self})},
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
	{vx = '*value = prim.v.x;'},
	{vy = '*value = prim.v.y;'},
	{vz = '*value = prim.v.z;'},
	{v = '*value = coordLen(prim.v, x);'},
	{eInt = '*value = prim.eInt;'},
	{P = '*value = calc_P(prim.rho, prim.eInt);'},
	{h = '*value = calc_h(prim.rho, calc_P(prim.rho, prim.eInt), prim.eInt);'},
}

GRHD.eigenStructFields = {
	{rho = 'real'},
	{v = 'real3'},
	{h = 'real'},
	{W = 'real'},
	{ATildeMinus = 'real'},
	{ATildePlus = 'real'},
	{VMinus = 'real'},
	{VPlus = 'real'},
	{CMinus = 'real'},
	{CPlus = 'real'},
	{Kappa = 'real'},
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
