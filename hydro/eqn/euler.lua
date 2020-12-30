local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local Euler = class(Equation)
Euler.name = 'euler'

-- ePot is the 6th param
-- which means it's now in the derivBuf, but it is always zero
-- so TODO a new variable for deriv size vs cons_t size?
--Euler.numStates = 6	

Euler.numWaves = 5
Euler.numIntStates = 5	-- don't bother integrate ePot

--Euler.roeUseFluxFromCons = true

Euler.initConds = require 'hydro.init.euler':getList()

function Euler:init(args)
	-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
	self.primVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m/s', variance='u'},			-- contravariant
		{name='P', type='real', units='kg/(m*s^2)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}

	self.consVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='m', type='real3', units='kg/(m^2*s)', variance='u'},	-- contravariant
		{name='ETotal', type='real', units='kg/(m*s^2)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}

	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end

	Euler.super.init(self, args)

	if require 'hydro.solver.meshsolver'.is(self.solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)

		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				vectorField = 'm',
				potentialField = 'mPot',
			
				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ 
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global <?=cons_t?> const * const Ujm = U - solver->stepsize.s<?=j?>;
		global <?=cons_t?> const * const Ujp = U + solver->stepsize.s<?=j?>;
		real drho_dx = (Ujp->rho - Ujm->rho) * (.5 / solver->grid_dx.s<?=j?>);
		source -= drho_dx * U->m.s<?=j?> / U->rho;
	}<? end ?>
]],
			})
		end
	end
end

function Euler:createInitState()
	Euler.super.createInitState(self)
	local double = false --solver.app.real == 'double'
	self:addGuiVars{	
		{name='heatCapacityRatio', value=7/5},				-- unitless
		{name='rhoMin', value=double and 1e-15 or 1e-7, units='kg/m^3'},
		{name='PMin', value=double and 1e-15 or 1e-7, units='kg/(m*s^2)'},
	}
end

-- this one calcs cell prims once and uses it for all sides
-- it is put here instead of in hydro/eqn/euler.cl so euler-burgers can override it
-- TODO move the sqrt() out of the loop altogether?
-- TODO allow module overriding in the markup somehow?
function Euler:initCodeModule_calcDTCell()
	local solver = self.solver
	solver.modules:add{
		name = self.symbols.calcDTCell,
		depends = table{
			self.solver.symbols.OOB,
			self.solver.symbols.SETBOUNDS,
			self.solver.solver_t,
			self.solver.coord.symbols.normal_t,
			self.symbols.primFromCons,
			self.symbols.eqn_guiVars_compileTime,
		},
		code = self:template[[
<? if require 'hydro.solver.gridsolver'.is(solver) then ?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell\
) {\
	real3 const x = cell->pos;\
	<?=prim_t?> W;\
	<?=primFromCons?>(&W, solver, U, x);\
	real const Cs = calc_Cs(solver, &W);\
	<? for side=0,solver.dim-1 do ?>{\
<? --\
if solver.coord.vectorComponent == 'cartesian' --\
and not require 'hydro.coord.cartesian'.is(solver.coord) --\
then --\
?>		real const dx = cell_dx<?=side?>(x);\
<? else --\
?>		real const dx = solver->grid_dx.s<?=side?>;\
<? end --\
?>\
		if (dx > 1e-7) {\
			/* use cell-centered eigenvalues */\
			real const v_n = normal_vecDotN1(normal_forSide<?=side?>(x), W.v);\
			real const lambdaMin = v_n - Cs;\
			real const lambdaMax = v_n + Cs;\
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
			absLambdaMax = max((real)1e-9, absLambdaMax);\
			/* TODO this should be based on coord + vectorComponent */\
			/* non-cartesian coord + cartesian component uses |u(x+dx)-u(x)| */\
			*(dt) = (real)min(*(dt), dx / absLambdaMax);\
		}\
	}<? end ?>\
}

<? else -- mesh solver ?>

#define <?=calcDTCell?>(\
	/*real * const */dt,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> const * const */U,\
	/*global <?=cell_t?> const * const */cell,\
	/*global <?=face_t?> const * const */faces,		/* [numFaces] */\
	/*global int const * const */cellFaceIndexes	/* [numCellFaceIndexes] */\
) {\
	real3 const x = cell->pos;\
	for (int i = 0; i < cell->faceCount; ++i) {\
		global <?=face_t?> const * const face = faces + cellFaceIndexes[i + cell->faceOffset];\
		real dx = face->area;\
		if (dx > 1e-7 && face->cells.x != -1 && face->cells.y != -1) {\
			/* all sides? or only the most prominent side? */\
			/* which should we pick eigenvalues from? */\
			/* use cell-centered eigenvalues */\
			<?=normal_t?> n = normal_forFace(face);\
			<?=eqn:consWaveCodePrefix('n', 'U', 'x'):gsub('\n', '\\\n\t\t\t')?>\
			real lambdaMin = <?=eqn:consMinWaveCode('n', 'U', 'x')?>;\
			real lambdaMax = <?=eqn:consMaxWaveCode('n', 'U', 'x')?>;\
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));\
			absLambdaMax = max((real)1e-9, absLambdaMax);\
			*(dt) = (real)min(*(dt), dx / absLambdaMax);\
		}\
	}\
}

<? end -- mesh vs grid solver ?>
]],
	}
end

function Euler:getModuleDepends_displayCode() 
	return table(Euler.super.getModuleDepends_displayCode(self)):append{
		self.gravOp.symbols.calcGravityAccel,
	}
end

-- don't use default
function Euler:initCodeModule_fluxFromCons() end
function Euler:initCodeModule_consFromPrim_primFromCons() end

Euler.solverCodeFile = 'hydro/eqn/euler.cl'

Euler.displayVarCodeUsesPrims = true

-- [=[
Euler.predefinedDisplayVars = {
-- [[
	'U rho',
	'U v',
	'U P',
--]]
--[[
	-- now that I've switched to components, I can't display these
	-- TODO ...unless I allow for multiple displays of the same displayVar...
	-- TODO TODO after building kernels, duplicate all display vars for each component type
	-- then things will look just like they did before, except with less kernels
	'U v x',
	'U v y',
	'U v z',
--]]
--[[	
	'U ePot',
	'U EPot',
	'U gravity',
--]]
	'U div v',
	'U curl v',
}
--]=]

function Euler:getDisplayVars()
	local vars = Euler.super.getDisplayVars(self)
	vars:append{
		-- TODO should the default display generated of variables be in solver units or SI units?
		-- should SI unit displays be auto generated as well?
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = W.P;', units='kg/(m*s^2)'},
		{name='eInt', code='value.vreal = calc_eInt(solver, &W);', units='m^2/s^2'},
		{name='eKin', code='value.vreal = calc_eKin(&W, x);', units='m^2/s^2'},
		{name='eTotal', code='value.vreal = U->ETotal / W.rho;', units='m^2/s^2'},
		{name='EInt', code='value.vreal = calc_EInt(solver, &W);', units='kg/(m*s^2)'},
		{name='EKin', code='value.vreal = calc_EKin(&W, x);', units='kg/(m*s^2)'},
		{name='EPot', code='value.vreal = U->rho * U->ePot;', units='kg/(m*s^2)'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{name='H', code='value.vreal = calc_H(solver, W.P);', units='kg/(m*s^2)'},
		{name='h', code='value.vreal = calc_h(solver, W.rho, W.P);', units='m^2/s^2'},
		{name='HTotal', code='value.vreal = calc_HTotal(W.P, U->ETotal);', units='kg/(m*s^2)'},
		{name='hTotal', code='value.vreal = calc_hTotal(W.rho, W.P, U->ETotal);', units='m^2/s^2'},
		{name='speed of sound', code='value.vreal = calc_Cs(solver, &W);', units='m/s'},
		{name='Mach number', code='value.vreal = coordLen(W.v, x) / calc_Cs(solver, &W);'},
		{name='temperature', code=self:template[[
<? local clnumber = require 'cl.obj.number' ?>
<? local materials = require 'hydro.materials' ?>
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
	value.vreal = calc_eInt(solver, &W) / C_v;
]], units='K'},
	}:append(self.gravOp and
		{{name='gravity', code=self:template[[
	if (!<?=OOB?>(1,1)) {
		value.vreal3 = <?=eqn.gravOp.symbols.calcGravityAccel?>(solver, U);
	}
]], type='real3', units='m/s^2'}} or nil
	)

	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->rho'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->rho'
		end,
		units = '1/s',
	} or nil)

	return vars
end

Euler.eigenVars = table{
	-- Roe-averaged vars
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	{name='hTotal', type='real', units='m^2/s^2'},
	-- derived vars
	{name='Cs', type='real', units='m/s'},
	{name='vSq', type='real', units='m^2/s^2'},
	{name='vL', type='real3', units='m/s'},
}

function Euler:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
real const <?=eqn.symbolPrefix?>Cs_nLen = normal_len(<?=n?>) * <?=eig?>->Cs;
real const <?=eqn.symbolPrefix?>v_n = normal_vecDotN1(<?=n?>, <?=eig?>->v);
]],	{
		eig = '('..eig..')',
		x = x,
		n = n,
	})
end

-- W is an extra param specific to Euler's calcDT in this case
-- but then I just explicitly wrote out the calcDT, so the extra parameters just aren't used anymore.
function Euler:consWaveCodePrefix(n, U, x)
	return self:template([[
<?=prim_t?> W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=x?>);
real const <?=eqn.symbolPrefix?>Cs_nLen = calc_Cs(solver, &W) * normal_len(<?=n?>);
real const <?=eqn.symbolPrefix?>v_n = normal_vecDotN1(<?=n?>, W.v);
]], {
		U = '('..U..')',
		n = n,
		x = x,
	})
end

function Euler:consWaveCode(n, U, x, waveIndex)
	if waveIndex == 0 then
		return self:template'(<?=eqn.symbolPrefix?>v_n - <?=eqn.symbolPrefix?>Cs_nLen)'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		return self:template'<?=eqn.symbolPrefix?>v_n'
	elseif waveIndex == 4 then
		return self:template'(<?=eqn.symbolPrefix?>v_n + <?=eqn.symbolPrefix?>Cs_nLen)'
	end
	error'got a bad waveIndex'
end

-- as long as U or eig isn't used, we can use this for both implementations
Euler.eigenWaveCode = Euler.consWaveCode

return Euler

--[[
scratch work:

P V = n R T = m R_spec T
ideal gas EOS: P V_m = R T
T = T_C + 273.15' C
R = gas constant = 8.314459848 J / (mol K) 
R_spec = 287.058 J / (kg K) for dry air 
R_spec = R / M = k_B / m
M = molar mass
k_B = Boltzmann constant = 1.3806485279e-23 J / K
m = mass within volume V	<-> ρ = m / V
R_spec = k_B ρ V

caloric perfect:
P = (gamma - 1) ρ e_int
gamma = C_p / C_v = adiabatic index / ratio of specific heats
e_int = C_v T = internal energy per unit mass / specific internal energy
C_v = specific heat at constant volume
C_p = specific heat at constant pressure
P = ρ R_spec T
e_int = C_v T
P = ρ (C_p - C_v) e_int / C_v = ρ (C_p / C_v - 1) e_int

0 C = 273.15 K

using some real-world numbers ...
P = (gamma - 1) ρ e_int
101325 kg / (m s^2) = (C_p - C_v) / C_v (1.2754 kg / m^3) C_v T 
101325 kg / (m s^2) = (1006 - 717.1) J / (kg K) (1.2754 kg / m^3) T 
T = 274.99364522457 K
T = 1.8436452245715 C ... should be
--]]
