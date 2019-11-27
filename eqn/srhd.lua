--[[
based on Marti 1998, Marti & Muller 2008, and maybe some of Font 2008 (but that's grhd)

honestly I developed a Marti & Muller SRHD solver
then I bumped it up to GRHD by incorporating (fixed) alphas betas and gammas
and then I thought "why not keep the old SRHD solver around"
so viola, here it is, unnecessarily available.
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local SRHD = class(Equation)
SRHD.name = 'SRHD'

SRHD.numStates = 11

SRHD.numWaves = 5
SRHD.numIntStates = 5

SRHD.hasCalcDTCode = true
SRHD.hasEigenCode = true 
SRHD.useConstrainU = true

-- SRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--SRHD.roeUseFluxFromCons = true
--SRHD.hasFluxFromConsCode = true

SRHD.initStates = require 'init.euler'

function SRHD:init(args)
	SRHD.super.init(self, args)

	self.cons_only_t = args.solver.app:uniqueName'cons_only_t'
	self.prim_only_t = args.solver.app:uniqueName'prim_only_t'

	local SRHDSelfGrav = require 'op.srhd-selfgrav'
	self.solver.ops:insert(SRHDSelfGrav{solver=self.solver})
end

--[[
because of the unique shape of prim_only_t and cons_only_t, I can't use the consVars for struct generation ...

2003 Marti & Muller show the state variables as D, S^i, tau ... for SRHD
...but the GRHD folks use D, S_i, tau ...
maybe Marti & Muller were lazy with notation since the metric is eta = diag(-1,1,1,1) and raising/lowering spatial doesn't matter ... ?

I used to keep the prim_only_t and cons_only_t as separate structs within separate,
but that doesn't mesh well with the code that automatically determines signs of reflections at boundaries,
so I'll try just merging this all into cons_t.
--]]
SRHD.consVars = table{
	-- cons_only_t
	{name='D', type='real', units='kg/m^3'},		-- D = rho W, W = unitless Lorentz factor
	{name='S', type='real3', units='kg/s^3'},		-- S_j = rho h W^2 v_j ... [rho] [h] [v] = kg/m^3 * m^2/s^2 * m/s = kg/s^3
	{name='tau', type='real', units='kg/(m*s^2)'},	-- tau = rho h W^2 - P ... [rho] [h] [W^2] = kg/m^3 * m^2/s^2 = kg/(m*s^2)
	
	-- prim_only_t
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	{name='eInt', type='real', units='m^2/s^2'},
	
	-- extra
	{name='ePot', type='real', units='m^2/s^2'},
}

function SRHD:getTypeCode()
	return table{
		SRHD.super.getTypeCode(self),
		template([[
typedef union {
	real ptr[5];
	struct {
		real D;
		real3 S;
		real tau;
	};
} <?=eqn.cons_only_t?>;

typedef union {
	real ptr[5];
	struct {
		real rho;
		real3 v;
		real eInt;
	};
} <?=eqn.prim_only_t?>;
]], 	{
			eqn = self,
		}),
	}:concat'\n'
end

function SRHD:createInitState()
	SRHD.super.createInitState(self)

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

function SRHD:getCommonFuncCode()
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

<?=eqn.cons_t?> consFromPrimOnly(constant <?=solver.solver_t?>* solver, <?=eqn.prim_only_t?> prim, real3 x) {
	real vSq = coordLenSq(prim.v, x);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(solver, prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 40-42:
	
	//rest-mass density = J^0 = rho u^0
	real D = prim.rho * W;	
	
	//momentum = T^0i = rho h u^0 u^i + P g^0i
	real3 S = real3_real_mul(prim.v, prim.rho * h * WSq);
	
	//energy = T^00 = rho h u^0 u^0 + P g^00
	real tau = prim.rho * h * WSq - D - P;
	
	return (<?=eqn.cons_t?>){
		.D = D,
		.S = S,
		.tau = tau,
		.rho = prim.rho,
		.v = prim.v,
		.eInt = prim.eInt,
		.ePot = 0,
	};
}


<?=eqn.cons_only_t?> consOnlyFromPrim(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
	real vSq = coordLenSq(U.v, x);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(solver, U.rho, U.eInt);
	real h = calc_h(U.rho, P, U.eInt);

	//2008 Font, eqn 40-42:
	
	//rest-mass density = J^0 = rho u^0
	real D = U.rho * W;	
	
	//momentum = T^0i = rho h u^0 u^i + P g^0i
	real3 S = real3_real_mul(U.v, U.rho * h * WSq);
	
	//energy = T^00 = rho h u^0 u^0 + P g^00
	real tau = U.rho * h * WSq - D - P;
	
	return (<?=eqn.cons_only_t?>){.D=D, .S=S, .tau=tau};
}

//PLM uses prim_only_t and cons_t, esp using the 'numIntStates' reals that they start with
//...and PLM uses consFromPrim and primFromCons

]], {
		eqn = self,
		solver = self.solver,
	})
end

function SRHD:getPrimConsCode() end

SRHD.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
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

	<?=code?>
	
	real eInt = calc_eInt_from_P(solver, rho, P);

	<?=eqn.prim_only_t?> prim = {.rho=rho, .v=v, .eInt=eInt};
	UBuf[index] = consFromPrimOnly(solver, prim, x);
}
]]

SRHD.solverCodeFile = 'eqn/srhd.cl'

-- TODO put in common parent of Euler, SRHD, GRHD
-- TODO use the automatic arbitrary finite difference generator in bssnok
-- k is 0,1,2
local function vorticity(eqn,k,result)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {
		name = 'vorticity '..xs[k+1],
		code = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = &buf[index - solver->stepsize.s<?=i?>];
		global const <?=eqn.cons_t?>* Uip = &buf[index + solver->stepsize.s<?=i?>];
		global const <?=eqn.cons_t?>* Ujm = &buf[index - solver->stepsize.s<?=j?>];
		global const <?=eqn.cons_t?>* Ujp = &buf[index + solver->stepsize.s<?=j?>];

		//TODO incorporate metric
		//TODO 3-vorticity vs 4-vorticity?

		real vim_j = Uim->v.s<?=j?>;
		real vip_j = Uip->v.s<?=j?>;
		
		real vjm_i = Ujm->v.s<?=i?>;
		real vjp_i = Ujp->v.s<?=i?>;
		
		<?=result?> = (vjp_i - vjm_i) / (2. * solver->grid_dx.s<?=i?>)
					- (vip_j - vim_j) / (2. * solver->grid_dx.s<?=j?>);
	}
]], 	{
			i = i,
			j = j,
			eqn = eqn,
			result = result,
		})
	}
end

function SRHD:getDisplayVars()
	local vars = table{
		{name='W based on D', code='*value = U->D / U->rho;'},
		{name='W based on v', code='*value = 1. / sqrt(1. - coordLenSq(U->v, x));'},
		
		{name='P', code='*value = calc_P(solver, U->rho, U->eInt);'},
		{name='h', code='*value = calc_h(U->rho, calc_P(solver, U->rho, U->eInt), U->eInt);'},
		
		{name='ePot', code='*value = U->ePot;'},
		
		{name='primitive reconstruction error', code=template([[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	{
		<?=eqn.cons_only_t?> U2 = consOnlyFromPrim(solver, *U, x);
		*value = 0;
		for (int j = 0; j < numIntStates; ++j) {
			*value += fabs(U->ptr[j] - U2.ptr[j]);
		}
	}
	]], {eqn=self})},
		{name='W error', code=[[
	real W1 = U->D / U->rho;
	real W2 = 1. / sqrt(1. - coordLenSq(U->v, x));
	*value = fabs(W1 - W2);
]]		},
	}
	
	if self.solver.dim == 2 then
		vars:insert(vorticity(self,2,'*value'))
	elseif self.solver.dim == 3 then
		local v = range(0,2):map(function(i) return vorticity(self,i,'value['..i..']') end)
		vars:insert{
			name = 'vorticityVec',
			code = template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
		++value;
	}<? end ?>
	value -= 3;
]], {v=v}), type='real3'}
	end

	return vars
end

SRHD.eigenVars = {
	-- needed by eigenvectors ...
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
	{name='lambdaMin', type='real'},
	{name='lambdaMax', type='real'},
}

function SRHD:eigenWaveCode(side, eig, x, waveIndex)
	if waveIndex == 0 then
		return '('..eig..').lambdaMin'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		-- v.x because v has been rotated so x points along the normal
		return '('..eig..').v.x'
	elseif waveIndex == 4 then
		return '('..eig..').lambdaMax'
	else
		error'got a bad waveIndex'
	end
end

-- used by HLL
-- extra params provided by calcDT, or calculated here if not provided (just like in Euler)
function SRHD:consWaveCodePrefix(side, U, x,
	prim, rho, eInt, vSq, P, h, csSq, cs
)
	return template([[	
<? if not prim then ?>
	real rho = <?=U?>.rho;
	real eInt = <?=U?>.eInt;
	real vSq = coordLenSq(<?=U?>.v, <?=x?>);
	real P = calc_P(solver, rho, eInt);
	real h = calc_h(rho, P, eInt);
	real csSq = solver->heatCapacityRatio * P / (rho * h);
	real cs = sqrt(csSq);
<? 
	prim = U
end ?>
	
	//for the particular direction
	real vi = <?=prim?>.v.s<?=side?>;
	real viSq = vi * vi;
	
	// Marti 1998 eqn 19
	// also Marti & Muller 2008 eqn 68
	// also Font 2008 eqn 106
	real discr = sqrt((1. - <?=vSq?>) * (1. - <?=vSq?> * <?=csSq?> - viSq * (1. - <?=csSq?>)));
	real _srhd_lambdaMin = (vi * (1. - <?=csSq?>) - <?=cs?> * discr) / (1. - <?=vSq?> * <?=csSq?>);
	real _srhd_lambdaMax = (vi * (1. - <?=csSq?>) + <?=cs?> * discr) / (1. - <?=vSq?> * <?=csSq?>);
	// v.x because v has been rotated so x points along the normal
	real v_n = prim.v.x;
]], {
		eqn = self,
		side = side,
		U = '('..U..')',
		x = x,
		-- extra params either provided or calculated
		-- TODO I almost need two prefixes ... one for all sides, and one per-side
		prim = '('..prim..')',
		rho = rho or 'rho',
		eInt = eInt or 'eInt',
		vSq = vSq or 'vSq',
		P = P or 'P',
		h = h or 'h',
		csSq = csSq or 'csSq',
		cs = cs or 'cs',
	})
end

function SRHD:consWaveCode(side, U, x, waveIndex)
	if waveIndex == 0 then
		return '_srhd_lambdaMin'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		return 'v_x'
	elseif waveIndex == 4 then
		return '_srhd_lambdaMax'
	else
		error'got a bad waveIndex'
	end
end



return SRHD
