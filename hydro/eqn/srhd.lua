--[[
based on 
Marti & Muller 2008
Marti 1998
Font "Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity" 2008 

honestly I developed a Marti & Muller SRHD solver
then I bumped it up to GRHD by incorporating (fixed) α's β's and γ's
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

-- TODO if we enable this we get NANs when waves hit the border.  Bug in the srhd boundary prim calculations?
-- but if we don't ... well ... it isn't technically correct
--SRHD.roeUseFluxFromCons = true

SRHD.initConds = require 'hydro.init.euler':getList()

function SRHD:init(args)
	local solver = assert(args.solver)

	--[[
	because of the unique shape of prim_only_t and cons_only_t, I can't use the consVars for struct generation ...

	2003 Marti & Muller show the state variables as D, S^i, tau ... for SRHD
	...but the GRHD folks use D, S_i, tau ...
	maybe Marti & Muller were lazy with notation since the metric is η = diag(-1,1,1,1) and raising/lowering spatial doesn't matter ... ?

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

	-- TODO how about anonymous structs, so we can copy out prim_only_t and cons_only_t?
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

	SRHD.super.init(self, args)
	
	self.symbols.cons_only_t = self.consOnlyStruct.typename
	self.symbols.prim_only_t = self.primOnlyStruct.typename

	if require 'hydro.solver.meshsolver':isa(self.solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		local SRHDSelfGrav = require 'hydro.op.srhd-selfgrav'
		self.gravOp = SRHDSelfGrav{solver=self.solver}
		self.solver.ops:insert(self.gravOp)
		
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
				-- TODO TODO now jacobi is crashing as well.  seems I broke something.
			}
		
			local SRHDNoDiv = class(NoDiv)
--[[ srhd needs to refresh primitives after updating conservatives in the nodiv step
-- unless you're doing the update-prim-and-recalc-cons method.  then no need.
			function SRHDNoDiv:step(dt)
				SRHDNoDiv.super.step(self, dt)
				self.solver:constrainU()
			end
--]]
			self.solver.ops:insert(SRHDNoDiv{
				solver = self.solver,
				potentialField = 'SPot',	-- TODO don't store this
				
				codeDepends = {
					assert(self.symbols.eqn_common),
				},
				
			--[=[ using 0 = div S = div (D^2 h v / ρ)
				vectorField = 'S',
			
				-- S_i = ρ h W^2 v_i = D^2 h v_i / ρ
				-- v_i = S_i ρ / (D^2 h)
				--
				-- S_i,j = (D^2 h v_i / ρ)_,j
				-- = (D^2 h / ρ)_,j v_i + (D^2 h / ρ) v_i,j
				-- = (D^2 h / ρ)_,j v_i + (D^2 h / ρ) v_i,j = 0
				-- δ^ij v_i,j = -δ^ij (D^2 h / ρ)_,j v_i / (D^2 h / ρ)
				-- TODO make this work for non-ident metric & add connections for covariant derivatives
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global const <?=cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
		real d_rho_h_WSq_over_dx = (
			Ujp->D * Ujp->D / Ujp->rho * calc_h(solver, Ujp->eInt)
			- Ujm->D * Ujm->D / Ujm->rho * calc_h(solver, Ujm->eInt)
		) * (.5 / solver->grid_dx.s<?=j?>);
		source -= d_rho_h_WSq_over_dx * U->S.s<?=j?> / (
			U->D * U->D / U->rho * calc_h(solver, U->eInt)
		);
	}<? end ?>
]],
			--]=]
			-- [=[ reading via div v, writing via div S
				readVectorField = function(op, offset, j)
					local function U(field) return 'U['..offset..'].'..field end
					--return U('S.s'..j)..' * '..U'rho'..' / ('..U'D'..' * '..U'D'..' * calc_h(solver, '..U'eInt'..'))'
					--return U('v.s'..j)	-- div v = 0
					return U'D'..' / '..U'rho'..' * '..U('v.s'..j)		-- div (W v) = ... W ... ?
				end,
				chargeCode = [[
	source += U->D / U->rho;
]],
				writeVectorField = function(op, dv)
					return self:template([[
#if 0	// adjust momentum, update velocity with 'constrainU' later
	U->S = real3_sub(U->S, real3_real_mul(<?=dv?>, U->D * U->D / U->rho * calc_h(solver, U->eInt)));
#endif

#if 1	// adjust velocity, update conservatives immediately
//// MODULE_DEPENDS: <?=prim_only_t?> <?=eqn_common?>
	<?=prim_only_t?> prim;
	primOnlyFromCons(&prim, solver, U, pt);
	//prim.v = real3_sub(prim.v, <?=dv?>);
	prim.v = real3_sub(prim.v, real3_real_mul(<?=dv?>, U->rho / U->D));
	consFromPrimOnly(U, solver, &prim, pt);
#endif

]], {dv=dv})
					-- TODO adjust tau <-> ETotal as well (right?)
				end,
			--]=]
			})
		end
	end
end

function SRHD:getSymbolFields()
	return SRHD.super.getSymbolFields(self):append{
		'cons_only_t',
		'prim_only_t',
	}
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

function SRHD:getModuleDepends_displayCode() 
	return table(SRHD.super.getModuleDepends_displayCode(self)):append(
		self.gravOp and {
			self.gravOp.symbols.calcGravityAccel,
		} or nil
	)
end

-- don't use default
function SRHD:initCodeModule_fluxFromCons() end
function SRHD:initCodeModule_calcDTCell() end

SRHD.solverCodeFile = 'hydro/eqn/srhd.cl'

SRHD.predefinedDisplayVars = {
	'U rho',
	'U v',
	'U tau',
	'U div v',
	'U curl v',
}

function SRHD:getDisplayVars()
	local vars = SRHD.super.getDisplayVars(self)
	vars:append{
		{name='W based on D', code='value.vreal = U->D / U->rho;'},
		{name='W based on v', code='value.vreal = 1. / sqrt(1. - coordLenSq(U->v, x));'},
		
		{name='P', code='value.vreal = calc_P(solver, U->rho, U->eInt);'},
		
		{name='EPot', code='value.vreal = U->rho * U->ePot;'},
		
		-- is this true for relativistic fluids?
		{name='S', code='value.vreal = calc_P(solver, U->rho, U->eInt) / pow(U->rho, (real)solver->heatCapacityRatio);'},
		
		{name='H', code='value.vreal = U->rho * calc_h(solver, U->eInt);', units='kg/(m*s^2)'},
		{name='h', code='value.vreal = calc_h(solver, U->eInt);'},
		
		{name='speed of sound', code='value.vreal = calc_Cs(solver, U->eInt);', units='m/s'},
		{name='Mach number', code='value.vreal = coordLen(U->v, x) / calc_Cs(solver, U->eInt);'},

		-- is this true for relativistic fluids?
		{name='temperature', code=self:template[[
<? local materials = require "hydro.materials" ?>
#define C_v				<?=("%.50f"):format(materials.Air.C_v)?>
value.vreal = U->eInt / C_v;
]], units='K'},
		
		{name='primitive reconstruction error', code=self:template[[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	{
		<?=cons_only_t?> U2;
		consOnlyFromPrim(&U2, solver, U, x);
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

	if self.gravOp then
		vars:insert{
			name = 'gravity',
			units = 'm/s^2',
			type = 'real3',
			code = self:template[[
	if (!<?=OOB?>(1,1)) {
		real W, dW_dt;
		real3 u, du_dt;
		<?=calcGravityAccel?>(&W, &u, &dW_dt, &du_dt, solver, U);
		value.vreal3 = du_dt;
	} else {
		value.vreal3 = real3_zero;
	}
]],
		}
	end


	vars:insert(self:createDivDisplayVar{field='v', units='1/s'} or nil)
	vars:insert(self:createCurlDisplayVar{field='v', units='1/s'} or nil)

	-- special for 1d_state_line
	vars:insert{
		name = 'state line',
		type = 'real3',
		units = '1',
		code = 'value.vreal3 = _real3(U->rho, coordLen(U->v, x) * sign(U->v.x), calc_P(solver, U->rho, U->eInt));',
	}

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

-- used by HLL
-- extra params provided by calcDT, or calculated here if not provided (just like in Euler)
-- but then I just explicitly wrote out the calcDT, so the extra parameters just aren't used anymore.
function SRHD:consWaveCodePrefix(n, U, x)
	U = '('..U..')'
	return self:template([[	
real const eInt = <?=U?>->eInt;

real const vSq = coordLenSq(<?=U?>->v, <?=x?>);
real const csSq = calc_CsSq(solver, eInt);
real const cs = sqrt(csSq);

/* for the particular direction */
real const vi = normal_vecDotN1(n, <?=U?>->v);
real const viSq = vi * vi;

/*  Marti 1998 eqn 19 */
/*  also Marti & Muller 2008 eqn 68 */
/*  also Font 2008 eqn 106 */
real const discr = sqrt((1. - vSq) * (1. - vSq * csSq - viSq * (1. - csSq)));
real const _srhd_lambdaMin = (vi * (1. - csSq) - cs * discr) / (1. - vSq * csSq);
real const _srhd_lambdaMax = (vi * (1. - csSq) + cs * discr) / (1. - vSq * csSq);
/*  v.x because v has been rotated so x points along the normal */
real const v_n = <?=U?>->v.x;
]], {
		n = n,
		U = U,
		x = x,
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
