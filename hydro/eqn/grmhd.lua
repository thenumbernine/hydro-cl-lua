-- BIG TODO.   I've moved prim_t into cons_t for GRHD ... because it's rare one is bound when the other is not
-- so .. this needs to undergo that as well

--[[
based on 2010 Anton et al 
... which looks to be a SRMHD implementation ...
where are varying metrics incorporated in 2010 Anton?
general relativistic ideal MHD
(TODO do a resistivie GRMHD, which incorporates E as well)
--]]
local class = require 'ext.class'
local Equation = require 'hydro.eqn.eqn'

local GRMHD = class(Equation)
GRMHD.name = 'GRMHD'
GRMHD.numStates = 9
GRMHD.numWaves = 8

-- GRMHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--GRMHD.roeUseFluxFromCons = true

GRMHD.initConds = require 'hydro.init.euler':getList()

-- TODO upgrade this to srhd: put these all in consVars and just make separate cons_only_ and prim_only_t
-- TODO also upgrade this to initCodeModules.  turn it into a struct, like srhd.
function GRMHD:getTypeCode()
	return self:template[[
typedef union {
	real ptr[9];
	struct {
		real rho;
		real3 v;
		
		//TODO Font 2008 uses 'eInt' and derives P (for a more flexible EOS)
		//but Anton just uses P
		real eInt;	
		
		real3 B;
		real divBPot;
	};
} <?=eqn.prim_t?>;

typedef union {
	real ptr[9];
	struct {
		real D;
		real3 S;
		real tau;
		real3 B;
		
		// TODO fix this.
		// it is here because prim_t is expected to be the same size as cons_t
		real divBPot;
	};
} <?=eqn.cons_t?>;
]]
end

function GRMHD:createInitState()
	GRMHD.super.createInitState(self)
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

-- YOU ARE HERE in converting stuff from SRHD to GRMHD

-- don't use default 
function GRMHD:initCodeModulePrimCons() end

GRMHD.solverCodeFile = 'hydro/eqn/srhd.cl'

GRMHD.displayVarCodeUsesPrims = true

function GRMHD:getDisplayVars()
	return GRMHD.super.getDisplayVars(self):append{
		{name='W', code='value.vreal = U->D / prim.rho;'},
		{name='primitive reconstruction error', code=self:template[[
			//prim have just been reconstructed from cons
			//so reconstruct cons from prims again and calculate the difference
			{
				<?=eqn.cons_t?> U2 = consFromPrim(prim, x);
				value.vreal = 0;
				for (int j = 0; j < numIntStates; ++j) {
					value.vreal += fabs(U->ptr[j] - U2.ptr[j]);
				}
			}
	]]},
	}
end

function GRMHD:getPrimDisplayVarCodePrefix()
	return self:template[[
	<?=eqn.prim_t?> prim = buf[index];
]]
end

GRMHD.eigenVars = {
	{name='rho', type='real'},
	{name='v', type='real3'},
	{name='h', type='real'},
	{name='W', type='real'},
	{name='ATildeMinus', type='real'},
	{name='ATildePlus', type='real'},
	{name='VMinus', type='real'},
	{name='VPlus', type='real'},
	{name='CMinus', type='real'},
	{name='CPlus', type='real'},
	{name='Kappa', type='real'},
}

return GRMHD
