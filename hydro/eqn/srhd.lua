--[[
based on Marti 1998, Marti & Muller 2008, and maybe some of Font 2008 (but that's grhd)

honestly I developed a Marti & Muller SRHD solver
then I bumped it up to GRHD by incorporating (fixed) alphas betas and gammas
and then I thought "why not keep the old SRHD solver around"
so viola, here it is.
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'
local Struct = require 'hydro.code.struct'

local SRHD = class(Equation)
SRHD.name = 'SRHD'

-- 11 without SPos, 12 with SPot
--SRHD.numStates = 11

SRHD.numWaves = 5
SRHD.numIntStates = 5

SRHD.hasCalcDTCode = true

SRHD.useSource = true		-- for connections
SRHD.useConstrainU = true

-- TODO if we enable this we get NANs when waves hit the border.  Bug in the srhd boundary prim calculations?
--SRHD.roeUseFluxFromCons = true

SRHD.initConds = require 'hydro.init.euler':getList()

function SRHD:init(args)
	local solver = assert(args.solver)

	--[[
	because of the unique shape of prim_only_t and cons_only_t, I can't use the consVars for struct generation ...

	2003 Marti & Muller show the state variables as D, S^i, tau ... for SRHD
	...but the GRHD folks use D, S_i, tau ...
	maybe Marti & Muller were lazy with notation since the metric is eta = diag(-1,1,1,1) and raising/lowering spatial doesn't matter ... ?

	I used to keep the prim_only_t and cons_only_t as separate structs within separate,
	but that doesn't mesh well with the code that automatically determines signs of reflections at boundaries,
	so I'll try just merging this all into cons_t.

	TODO redo the srhd equations for a background grid metric, and take note of covariance/contravariance
	--]]
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

	self.consVars = table()
	:append(self.consOnlyStruct.vars)
	:append(self.primOnlyStruct.vars)
	:append{
		-- extra
		{name='ePot', type='real', units='m^2/s^2'},
	}

	if args.incompressible then
		self.consVars:insert{name='SPot', type='real', units='kg*m/s^2'}
	end

	self.consOnlyStruct:makeType()
	self.primOnlyStruct:makeType()

	self.cons_only_t = self.consOnlyStruct.typename
	self.prim_only_t = self.primOnlyStruct.typename

	SRHD.super.init(self, args)

	if require 'hydro.solver.meshsolver'.is(self.solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		local SRHDSelfGrav = require 'hydro.op.srhd-selfgrav'
		self.gravOp = SRHDSelfGrav{solver=self.solver}
		self.solver.ops:insert(self.gravOp)
		
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				vectorField = 'S',
				potentialField = 'SPot',
			
				-- S_i = ρ h W^2 v_i = D^2 h / ρ
				-- S_i,j = (D^2 h v_i / ρ)_,j
				-- = (D^2 h / ρ)_,j v_i + (D^2 h / ρ) v_i,j
				-- = (D^2 h / ρ)_,j v_i + (D^2 h / ρ) v_i,j = 0
				-- δ^ij v_i,j = -δ^ij (D^2 h / ρ)_,j v_i / (D^2 h / ρ)
				-- TODO make this work for non-ident metric & add connections for covariant derivatives
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
		real d_rho_h_WSq_over_dx = (
			Ujp->D * Ujp->D / Ujp->rho * calc_h(Ujp->rho, calc_P(solver, Ujp->rho, Ujp->eInt), Ujp->eInt)
			- Ujm->D * Ujm->D / Ujm->rho * calc_h(Ujm->rho, calc_P(solver, Ujm->rho, Ujm->eInt), Ujm->eInt)
		) * (.5 / solver->grid_dx.s<?=j?>);
		source -= d_rho_h_WSq_over_dx * U->S.s<?=j?> / (
			U->D * U->D / U->rho * calc_h(U->rho, calc_P(solver, U->rho, U->eInt), U->eInt)
		);
	}<? end ?>
]],
			})
		end
	end
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

function SRHD:initCodeModules()
	SRHD.super.initCodeModules(self)
	local solver = self.solver

	solver.modules:add{
		name = 'cons_only_t,prim_only_t',
		structs = {self.primOnlyStruct, self.consOnlyStruct},
	}
end

function SRHD:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'normal_t',	-- normal_*
			'eqn.common',	-- calc_P
		},
		code = self:template[[
<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	real v_n = normal_vecDotN1(n, U.v);
	real P = calc_P(solver, U.rho, U.eInt);

	<?=eqn.cons_t?> F = {
		.D = U.D * v_n,
		.S = real3_add(
			real3_real_mul(U.S, v_n),
			_real3(
				normal_u1x(n) * P,
				normal_u1y(n) * P,
				normal_u1z(n) * P
			)
		),
		.tau = U.tau * v_n + P * v_n,
		.rho = 0,
		.v = real3_zero,
		.eInt = 0,
		.ePot = 0,
	};
	
	return F;
}
]],
	}
end

function SRHD:getModuleDependsApplyInitCond()
	return {
		'eqn.common',	-- calc_eInt_from_P, consFromPrimOnly
	}
end

function SRHD:getModuleDependsSolver()
	return {
		'cons_only_t,prim_only_t',
		'eqn.common',
		'coord_lower',
	}
end

function SRHD:initCodeModuleCommon()
	self.solver.modules:add{
		name = 'eqn.common',
		depends = {'cons_only_t,prim_only_t'},
		code = self:template[[

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

//build the cons_only_t from the cons_t's prim_only_t fields
//used for checking the error between cons_only_t and its prim-reconstructed-from-cons_only_t
<?=eqn.cons_only_t?> consOnlyFromPrim(
	constant <?=solver.solver_t?>* solver, 
	<?=eqn.cons_t?> U,
	real3 x
) {
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

<?=eqn.prim_only_t?> primOnlyFromCons(
	constant <?=solver.solver_t?>* solver, 
	<?=eqn.cons_t?> U,
	real3 x
) {
	return (<?=eqn.prim_only_t?>){.rho=U.rho, .v=U.v, .eInt=U.eInt};
}

//PLM uses prim_only_t and cons_t, esp using the 'numIntStates' reals that they start with
//...and PLM uses consFromPrim and primFromCons

]],
	}
end

SRHD.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	
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
	real3 D = real3_zero;
	real3 B = real3_zero;

	<?=code?>
	
	real eInt = calc_eInt_from_P(solver, rho, P);

	<?=eqn.prim_only_t?> prim = {.rho=rho, .v=v, .eInt=eInt};
	UBuf[index] = consFromPrimOnly(solver, prim, x);
}
]]

SRHD.solverCodeFile = 'hydro/eqn/srhd.cl'

function SRHD:getDisplayVars()
	local vars = SRHD.super.getDisplayVars(self)
	vars:append{
		{name='W based on D', code='value.vreal = U->D / U->rho;'},
		{name='W based on v', code='value.vreal = 1. / sqrt(1. - coordLenSq(U->v, x));'},
		
		{name='P', code='value.vreal = calc_P(solver, U->rho, U->eInt);'},
		{name='h', code='value.vreal = calc_h(U->rho, calc_P(solver, U->rho, U->eInt), U->eInt);'},
		
		{name='ePot', code='value.vreal = U->ePot;'},
		
		{name='primitive reconstruction error', code=self:template[[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	{
		<?=eqn.cons_only_t?> U2 = consOnlyFromPrim(solver, *U, x);
		value.vreal = 0;
		for (int j = 0; j < numIntStates; ++j) {
			value.vreal += fabs(U->ptr[j] - U2.ptr[j]);
		}
	}
]]
		},
		{name='W error', code=[[
	real W1 = U->D / U->rho;
	real W2 = 1. / sqrt(1. - coordLenSq(U->v, x));
	value.vreal = fabs(W1 - W2);
]]		},
	}

	vars:insert(self:createDivDisplayVar{field='v', units='1/s'} or nil)
	vars:insert(self:createCurlDisplayVar{field='v', units='1/s'} or nil)

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

function SRHD:eigenWaveCode(n, eig, x, waveIndex)
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
-- but then I just explicitly wrote out the calcDT, so the extra parameters just aren't used anymore.
function SRHD:consWaveCodePrefix(n, U, x)
	U = '('..U..')'
	return self:template([[	
	real eInt = <?=U?>.eInt;
	
	real vSq = coordLenSq(<?=U?>.v, <?=x?>);
	real P = calc_P(solver, <?=U?>.rho, eInt);
	real h = calc_h(<?=U?>.rho, P, eInt);
	real csSq = solver->heatCapacityRatio * P / (<?=U?>.rho * h);
	real cs = sqrt(csSq);
	
	//for the particular direction
	real vi = normal_vecDotN1(n, <?=U?>.v);
	real viSq = vi * vi;
	
	// Marti 1998 eqn 19
	// also Marti & Muller 2008 eqn 68
	// also Font 2008 eqn 106
	real discr = sqrt((1. - vSq) * (1. - vSq * csSq - viSq * (1. - csSq)));
	real _srhd_lambdaMin = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq);
	real _srhd_lambdaMax = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq);
	// v.x because v has been rotated so x points along the normal
	real v_n = <?=U?>.v.x;
]], {
		n = n,
		U = U,
		x = x,
		-- extra params either provided or calculated
		-- TODO I almost need two prefixes ... one for all sides, and one per-side
		rho = U..'.rho',
		v = U..'.v',
		P = U..'.P',
		eInt = U..'.eInt',
	})
end

function SRHD:consWaveCode(n, U, x, waveIndex)
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
