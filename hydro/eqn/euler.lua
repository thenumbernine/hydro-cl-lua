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

Euler.initConds = require 'hydro.init.euler':getList()

function Euler:init(args)
	Euler.super.init(self, args)
	local solver = self.solver

	self:buildSelfGrav()

	if args.incompressible then
		if not require 'hydro.solver.meshsolver':isa(solver) then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			solver.ops:insert(NoDiv{
				solver = solver,
				potentialField = 'mPot',	-- TODO don't store this

				--[=[ using div (m/rho) = 0, solve for div m:
				-- TODO field as a read function, and just read
				vectorField = 'm',

				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global <?=cons_t?> const * const Ujm = U - solver.stepsize[<?=j?>];
		global <?=cons_t?> const * const Ujp = U + solver.stepsize[<?=j?>];
		real drho_dx = (Ujp.rho - Ujm.rho) * (.5 / solver.grid_dx[<?=j?>]);
		source -= drho_dx * U.m[<?=j?>] / U.rho;
	}<? end ?>
]],
				--]=]
				-- [=[ reading via div(v), writing via div(m)
				readVectorField = function(op,offset,j)
					local function U(field) return 'U['..offset..'].'..field end
					return U('m['..j..']')..' / '..U'rho'
				end,
				writeVectorField = function(op,dv)
					return self:template([[
#if 0	// just adjust velocity
	U.m -= <?=dv?> * U.rho;
#endif

#if 0	// adjust ETotal as well
	U.ETotal -= .5 * U.rho * coordLenSq(U.m, pt);
	U.m -= <?=dv?> * U.rho;
	U.ETotal += .5 * U.rho * coordLenSq(U.m, pt);
#endif

#if 1 	// recalculate cons
//// MODULE_DEPENDS: <?=Solver?>
	<?=prim_t?> W = <?=Solver?>::Eqn::primFromCons(solver, U, pt);
	W.v -= <?=dv?>;
	U = <?=Solver?>::Eqn::consFromPrim(solver, W, pt);
#endif
]], {dv=dv})
				end,
				codeDepends = {
					self.solver.symbols.Solver,
				},
				--]=]
			})
		else
--[=[ still working on this ...
			-- TODO instead of separating this from NoDiv, just abstract the face iteration
			-- but you still have to deal with extra buffers
			-- and then there's the biggest can of worms:
			-- the # faces ~= (usually >=) the # cells
			-- so our A x = b function will be overconstrained
			local NoDivUnstructured = require'hydro.op.nodiv-unstructured'
			solver.ops:insert(NoDivUnstructured{
				solver = solver,
				potentialField = 'mPot',
				vectorField = 'm',
			})
--]=]
		end
	end

	self.viscosity = args.viscosity
	if self.viscosity then
		local materials = require 'hydro.materials'
		-- this just passes on to createInitState's addGuiVars
		-- default value for air is 1e-5 or so
		-- using .01 with Kelvin-Helmholtz is noticeable
		-- using anything higher tends to explode with explicit update
		self.shearViscosity = args.shearViscosity or materials.Air.shearViscosity
		self.heatConductivity = args.heatConductivity or materials.Air.heatConductivity

		if self.viscosity == 'rhs-explicit' then
			solver.ops:insert(require 'hydro.op.viscous-explicit'{
				solver = solver,
			})
		elseif self.viscosity == 'rhs-implicit' then
			-- handle the integration with backward-euler integrator
			-- TODO don't make a new module, instead just add the explicit update module but use int/be
			solver.ops:insert(require 'hydro.op.viscous-implicit'{
				solver = solver,
			})
		elseif self.viscosity == 'flux' then
			solver.ops:insert(require 'hydro.op.viscous-flux'{
				solver = solver,
			})
		else
			error("got an unknown viscosity method: "..self.viscosity)
		end
	end
end

function Euler:buildVars(args)
	-- this has turned very ugly
	-- since some eqns have no prim_t and have it typedef'd to cons_t
	-- and that is determined by primVars==nil
	-- that means we can't assign a default empty primVars and
	-- have to check for primVars' existence
	self.primVars = self.primVars or table()
	self.consVars = self.consVars or table()

	-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
	self.primVars:append{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m/s', variance='u'},			-- contravariant
		{name='P', type='real', units='kg/(m*s^2)'},

		-- used dynamically by op/selfgrav, but optionally can be initialized statically for constant/background potential energy
		-- TODO in the static case, merge into cell_t to save memory & flops?
		{name='ePot', type='real', units='m^2/s^2'},
	}

	self.consVars:append{
		{name='rho', type='real', units='kg/m^3'},
		{name='m', type='real3', units='kg/(m^2*s)', variance='u'},	-- contravariant
		{name='ETotal', type='real', units='kg/(m*s^2)'},	-- per volume.  energy units are kg m^2 / s^2, but energy per volume is kg / (m s^2)
		{name='ePot', type='real', units='m^2/s^2'},
	}

	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end
end

function Euler:buildSelfGrav()
	local solver = self.solver
	if require 'hydro.solver.meshsolver':isa(solver) then
		print("not using selfgrav with mesh solvers yet")
	else
		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver = solver}
		solver.ops:insert(self.gravOp)
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
	if self.viscosity then
		self:addGuiVars{
			{name='shearViscosity', value=assert(self.shearViscosity), units='kg/(m*s)'},
			{name='heatConductivity', value=assert(self.heatConductivity), units='kg/(m*s^3*K)'},
		}
	end
end

-- don't use default
function Euler:initCodeModule_fluxFromCons() end
function Euler:initCodeModule_consFromPrim_primFromCons() end

Euler.solverCodeFile = 'hydro/eqn/euler.clcpp'

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
		{name='eInt', code=self:template'value.vreal = Eqn::calc_eInt(solver, W);', units='m^2/s^2'},
		{name='eKin', code=self:template'value.vreal = Eqn::calc_eKin(W, x);', units='m^2/s^2'},
		{name='eTotal', code='value.vreal = U.ETotal / W.rho;', units='m^2/s^2'},
		{name='EInt', code=self:template'value.vreal = Eqn::calc_EInt(solver, W);', units='kg/(m*s^2)'},
		{name='EKin', code=self:template'value.vreal = Eqn::calc_EKin(W, x);', units='kg/(m*s^2)'},
		{name='EPot', code='value.vreal = U.rho * U.ePot;', units='kg/(m*s^2)'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver.heatCapacityRatio);'},
		{name='H', code=self:template'value.vreal = Eqn::calc_H(solver, W.P);', units='kg/(m*s^2)'},
		{name='h', code=self:template'value.vreal = Eqn::calc_h(solver, W.rho, W.P);', units='m^2/s^2'},
		{name='HTotal', code=self:template'value.vreal = Eqn::calc_HTotal(W.P, U.ETotal);', units='kg/(m*s^2)'},
		{name='hTotal', code=self:template'value.vreal = Eqn::calc_hTotal(W.rho, W.P, U.ETotal);', units='m^2/s^2'},
		{name='speed of sound', code=self:template'value.vreal = Eqn::calc_Cs(solver, W);', units='m/s'},
		{name='Mach number', code=self:template'value.vreal = coordLen(W.v, x) / Eqn::calc_Cs(solver, W);'},
		{name='temperature', code=self:template'value.vreal = Eqn::calc_T(U, x);', units='K'},
	}:append(self.gravOp and
		{{name='gravity', code=self:template[[
if (!OOB<dim>(solver, i, 1,1)) {
//// MODULE_DEPENDS: <?=eqn.gravOp.symbols.calcGravityAccel?>
	value.vreal3 = <?=eqn.gravOp.symbols.calcGravityAccel?><dim>(solver, &U, x);
} else {
	value.vreal3 = {};
}
]], type='real3', units='m/s^2'}} or nil
	)

	vars:insert(self:createDivDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'.m['..j..'] / '..U..'.rho'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'.m['..j..'] / '..U..'.rho'
		end,
		units = '1/s',
	} or nil)

	-- special for 1d_state_line
	vars:insert{
		name = 'state line',
		type = 'real3',
		units = '1',
		code = 'value.vreal3 = real3(W.rho, coordLen(W.v, x) * sign(W.v.x), W.P);',
	}

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
