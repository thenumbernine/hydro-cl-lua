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

SRHD.mirrorVars = {
	{'cons.S.x', 'prim.v.x'},
	{'cons.S.y', 'prim.v.y'},
	{'cons.S.z', 'prim.v.z'},
}

-- SRHD fluxFromCons will need prims passed to it as well
-- which means overriding the code that calls this? or the calc flux code?
--SRHD.roeUseFluxFromCons = true

SRHD.hasEigenCode = true 
SRHD.hasCalcDT = true
SRHD.useConstrainU = true

SRHD.initStates = require 'init.euler'

function SRHD:init(solver)
	SRHD.super.init(self, solver)

	self.cons_only_t = self:unique'cons_only_t'

	local SRHDSelfGrav = require 'solver.srhd-selfgrav'
	solver.ops:insert(SRHDSelfGrav{solver=solver})
end

--[[
because of the unique shape of prim_t and cons_only_t, I can't use the consVars for struct generation ...

2003 Marti & Muller show the state variables as D, S^i, tau ... for SRHD
...but the GRHD folks use D, S_i, tau ...
maybe Marti & Muller were lazy with notation since the metric is eta = diag(-1,1,1,1) and raising/lowering spatial doesn't matter ... ?
--]]
function SRHD:getTypeCode()
	return template([[
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
} <?=eqn.prim_t?>;

typedef union {
	real ptr[11];
	struct {
		<?=eqn.cons_only_t?> cons;
		<?=eqn.prim_t?> prim;
		real ePot;
	};
} <?=eqn.cons_t?>;
]], {
	eqn = self,
})
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
		{name='solvePrimMaxIter', type='int', value=10},	-- value=1000},

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

<?=eqn.cons_only_t?> consOnlyFromPrim(<?=eqn.prim_t?> prim, real3 x) {
	real vSq = coordLenSq(prim.v, x);
	real WSq = 1. / (1. - vSq);
	real W = sqrt(WSq);
	real P = calc_P(prim.rho, prim.eInt);
	real h = calc_h(prim.rho, P, prim.eInt);

	//2008 Font, eqn 40-42:
	
	//rest-mass density = J^0 = rho u^0
	real D = prim.rho * W;	
	
	//momentum = T^0i = rho h u^0 u^i + P g^0i
	real3 S = real3_scale(prim.v, prim.rho * h * WSq);
	
	//energy = T^00 = rho h u^0 u^0 + P g^00
	real tau = prim.rho * h * WSq - D - P;
	
	return (<?=eqn.cons_only_t?>){.D=D, .S=S, .tau=tau};
}

//PLM uses prim_t and cons_t, esp using the 'numIntStates' reals that they start with
//...and PLM uses consFromPrim and primFromCons

]], {
		eqn = self,
	})
end

function SRHD:getPrimConsCode() end

SRHD.initStateCode = [[
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
	//ignored:
	real3 B = _real3(0,0,0);

	<?=code?>
	
	real eInt = calc_eInt_from_P(rho, P);
	real vSq = coordLenSq(v, x);
	real W = 1. / sqrt(1. - vSq);
	real h = calc_h(rho, P, eInt);

	<?=eqn.prim_t?> prim = {.rho=rho, .v=v, .eInt=eInt};
	UBuf[index] = (<?=eqn.cons_t?>){
		.prim = prim,
		.cons = consOnlyFromPrim(prim, x),
	};
}
]]

function SRHD:getSolverCode()
	return template(file['eqn/srhd.cl'], {
		eqn = self,
		solver = self.solver,
	})
end

-- TODO put in common parent of Euler, SRHD, GRHD
-- TODO use the automatic arbitrary finite difference generator in bssnok
-- k is 0,1,2
local function vorticity(eqn,k,result)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {['vorticity '..xs[k+1]] = template([[
	
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {
		global const <?=eqn.prim_t?>* prim_im = &buf[index - stepsize.s<?=i?>].prim;
		global const <?=eqn.prim_t?>* prim_ip = &buf[index + stepsize.s<?=i?>].prim;
		global const <?=eqn.prim_t?>* prim_jm = &buf[index - stepsize.s<?=j?>].prim;
		global const <?=eqn.prim_t?>* prim_jp = &buf[index + stepsize.s<?=j?>].prim;

		//TODO incorporate metric
		//TODO 3-vorticity vs 4-vorticity?

		real vim_j = prim_im->v.s<?=j?>;
		real vip_j = prim_ip->v.s<?=j?>;
		
		real vjm_i = prim_jm->v.s<?=i?>;
		real vjp_i = prim_jp->v.s<?=i?>;
		
		<?=result?> = (vjp_i - vjm_i) / (2. * grid_dx<?=i?>)
				- (vip_j - vim_j) / (2. * grid_dx<?=j?>);
	}
]], {
		i = i,
		j = j,
		eqn = eqn,
		result = result,
	})}
end

function SRHD:getDisplayVars()
	local vars = table{
		{D = '*value = U->cons.D;'},
		{S = '*valuevec = U->cons.S;', type='real3'},
		{tau = '*value = U->cons.tau;'},
		{['W based on D'] = '*value = U->cons.D / U->prim.rho;'},
		{['W based on v'] = '*value = 1. / sqrt(1. - coordLenSq(U->prim.v, x));'},
		
		{rho = '*value = U->prim.rho;'},
		{v = '*valuevec = U->prim.v;', type='real3'},
		{eInt = '*value = U->prim.eInt;'},
		{P = '*value = calc_P(U->prim.rho, U->prim.eInt);'},
		{h = '*value = calc_h(U->prim.rho, calc_P(U->prim.rho, U->prim.eInt), U->prim.eInt);'},
		
		{ePot = '*value = U->ePot;'},
		
		{['primitive reconstruction error'] = template([[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	{
		<?=eqn.cons_only_t?> U2 = consOnlyFromPrim(U->prim, x);
		*value = 0;
		for (int j = 0; j < numIntStates; ++j) {
			*value += fabs(U->cons.ptr[j] - U2.ptr[j]);
		}
	}
	]], {eqn=self})},
		{['W error'] = [[
	real W1 = U->cons.D / U->prim.rho;
	real W2 = 1. / sqrt(1. - coordLenSq(U->prim.v, x));
	*value = fabs(W1 - W2);
]]		},
	}
	
	if self.solver.dim == 2 then
		vars:insert(vorticity(self,2,'*value'))
	elseif self.solver.dim == 3 then
		local v = range(0,2):map(function(i) return vorticity(self,i,'value['..i..']') end)
		vars:insert{vorticityVec = template([[
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
	{lambdaMin = 'real'},
	{lambdaMax = 'real'},
}

function SRHD:eigenWaveCode(side, eig, x, waveIndex)
	if waveIndex == 0 then
		return '('..eig..')->lambdaMin'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		-- v.x because v has been rotated so x points along the normal
		return '('..eig..')->v.x'
	elseif waveIndex == 4 then
		return '('..eig..')->lambdaMax'
	else
		error'got a bad waveIndex'
	end
end

--[=[
function SRHD:getFluxFromConsCode()
	return template([[
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	real vi = U.prim.v.s<?=side?>;
	real P = calc_P(U.prim.rho, U.prim.eInt);

	<?=eqn.cons_t?> F;
	F.cons.D = U.cons.D * vi;
	F.cons.S = real3_scale(U.cons.S, vi);
	F.cons.S.s<?=side?> += P;
	F.cons.tau = U.cons.tau * vi + P * vi;
	
	//make sure the rest is zero ...
	F.prim = (<?=eqn.prim_t?>){
		.rho = 0,
		.v = {.s={0,0,0}},
		.eInt = 0,
	};
	F.ePot = 0;

	return F;
}
<? end ?>
]], {
		eqn = self,
		solver = self.solver,
	})
end
--]=]

return SRHD
