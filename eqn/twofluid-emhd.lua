--[[
2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations"

the solver/twofluid-emhd-behavior.lua instanciates three separate solvers:
	ion, electron, maxwell
However OpenCL is having trouble with this.
Maybe I have an oob memory write somewhere?
Maybe running code from two separate programs doesn't propery block on the GPU?
Either way, here's the three equations combined into one.


When it comes to curvilinear coordinates, my Maxwell solver is unweighted by the metric volume,
however my Euler solver is.
I could work around this by scaling down the Maxwell eigenvalues by sqrt(det(g))
 then calcDeriv will scale all variables back up -- and the Maxwell ones will be identity weighted once again.
 then in the post-flux code, transform D and B by g_ij
 and in the source terms, add to ion_ and elec_ m^i connection values

--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local template = require 'template'


local TwoFluidEMHD = class(Equation)

local fluids = table{'ion', 'elec'}
TwoFluidEMHD.fluids = fluids

TwoFluidEMHD.postComputeFluxCode = [[
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
		flux.D = real3_real_mul(coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = real3_real_mul(coord_lower(flux.B, x), _1_sqrt_det_g);
]]

TwoFluidEMHD.name = 'TwoFluidEMHD'

-- set this to false to integrate the EM D,B,phi,psi variables together with the ion and electron Euler fluid equations
-- set this to true to integrate them separately using an implicit method
-- TODO FINISHME -- add the integration part
TwoFluidEMHD.implicitEMIntegration = false

if not TwoFluidEMHD.implicitEMIntegration then
	TwoFluidEMHD.numWaves = 18
	TwoFluidEMHD.numIntStates = 18
else
	TwoFluidEMHD.numWaves = 10
	TwoFluidEMHD.numIntStates = 10
end

TwoFluidEMHD.consVars = table{
	--integration variables		
	{name='ion_rho', type='real', units='kg/m^3'},
	{name='ion_m', type='real3', units='kg/(m^2*s)'},	-- m^i
	{name='ion_ETotal', type='real', units='kg/(m*s^2)'},
	
	{name='elec_rho', type='real', units='kg/m^3'},
	{name='elec_m', type='real3', units='kg/(m^2*s)'},		-- m^i
	{name='elec_ETotal', type='real', units='kg/(m*s^2)'},

	{name='D', type='real3', units='C/m^2'},			-- D_i
	{name='B', type='real3', units='kg/(C*s)'},			-- B_i
	{name='phi', type='real', units='C/m^2'},			-- div D potential
	{name='psi', type='real', units='kg/(C*s)'},		-- div B potential

	--extra	
	{name='ePot', type='real', units='m^2/s^2'},
}

TwoFluidEMHD.primVars = table{
	--integration variables		
	{name='ion_rho', type='real', units='kg/m^3'},
	{name='ion_v', type='real3', units='m/s'},
	{name='ion_P', type='real', units='kg/(m*s^2)'},
	
	{name='elec_rho', type='real', units='kg/m^3'},
	{name='elec_v', type='real3', units='m/s'},
	{name='elec_P', type='real', units='kg/(m*s^2)'},

	{name='D', type='real3', units='C/m^2'},			-- D_i
	{name='B', type='real3', units='kg/(C*s)'},			-- B_i
	{name='phi', type='real', units='C/m^2'},			-- div D potential
	{name='psi', type='real', units='kg/(C*s)'},		-- div B potential
	
	--extra	
	{name='ePot', type='real', units='m^2/s^2'},
}

-- not sure it's working right ...
--TwoFluidEMHD.hasCalcDTCode = true

TwoFluidEMHD.hasEigenCode = true
TwoFluidEMHD.hasFluxFromConsCode = true
TwoFluidEMHD.roeUseFluxFromCons = true
TwoFluidEMHD.useSourceTerm = true
TwoFluidEMHD.useConstrainU = true

TwoFluidEMHD.useEulerInitState = true	-- default to true.

-- TODO this has the symptoms of all the other CL kernels that intel's compiler bugged out on
-- i.e. takes a minute to run a kernel on the GPU
function TwoFluidEMHD:init(args)
	
	self.scalar = 'real'
	--self.scalar = 'cplx'
	self.vec3 = self.scalar..'3'

	self.susc_t = self.scalar

	-- set this to 'true' to use the init.euler states
	-- these states are only provided in terms of a single density and pressure variable,
	-- so the subsequent ion and electron densities and pressures must be derived from this.
	self.useEulerInitState = args.useEulerInitState

	if self.useEulerInitState then
		self.initStates = require 'init.euler'
	else
		self.initStates = require 'init.twofluid-emhd'
	end


	TwoFluidEMHD.super.init(self, args)


--	local NoDiv = require 'op.nodiv'()
--	self.solver.ops:insert(NoDiv{solver=self.solver})	-- nodiv on maxwell ... or just use potentials 

	local TwoFluidSelfGrav = require 'op.twofluid-selfgrav'
	self.gravOp = TwoFluidSelfGrav{solver=self.solver}
	self.solver.ops:insert(self.gravOp)
end

function TwoFluidEMHD:createInitState()
	TwoFluidEMHD.super.createInitState(self)
	
	local speedOfLight = 1	--require 'constants'.speedOfLight_in_m_per_s -- m/s
	local CoulombConstant = 1	--require 'constants'.CoulombConstant_in_kg_m3_per_C2_s2  -- (kg m^3)/(C^2 s)

	--[[
	local vacuumPermittivity = 1 / (4 * math.pi * CoulombConstant)	-- (C^2 s^2)/(kg m^3)
	local vacuumPermittivityTimesSpeedOfLightSquared = speedOfLight * speedOfLight * vacuumPermittivity	-- C^2/(kg m)
	local vacuumPermeability = 1 / vacuumPermittivityTimesSpeedOfLightSquared	-- (kg m)/C^2
	--]]
	-- [[
	local vacuumPermittivity = 1
	local vacuumPermeability = 1
	--]]
	
	self:addGuiVars(table{
		--never given, only stated as "speeds for the Maxwell equation"
		-- of course, they're associated with the potentials, so... they could be arbitrary
		{name='divPsiWavespeed', value=speedOfLight, units='m/s'},
		{name='divPhiWavespeed', value=speedOfLight, units='m/s'},
		
		-- gamma = heat capacity ratio
		{name='heatCapacityRatio', value=5/3},

		-- m_i = ion mass
		{name='ionMass', value=1, units='kg'},
	
		-- c = speed of light
		{name='speedOfLight', value=speedOfLight, units='m/s'},
		
		-- m = m_i / m_e
		-- https://en.wikipedia.org/wiki/Proton-to-electron_mass_ratio
		-- m_i / m_e = 1836.15267389
		{name='ionElectronMassRatio', value=100},
		
		-- r_i = q_i / m_i
		{name='ionChargeMassRatio', value=1, units='C/kg'},

		-- hmm ...
		--{name='sqrt_mu', value=math.sqrt(vacuumPermittivity), units='(kg m)^.5/C'},
		--{name='sqrt_eps', value=math.sqrt(vacuumPermeability), units='(C*s)/(kg*m^3)^.5'},
		{name='sqrt_mu', value=1, units='(kg*m)^.5/C'},
		{name='sqrt_eps', value=1, units='(C*s)/(kg*m^3)^.5'},
	
		-- lambda_d = ion Debye length
		-- lambda_d = sqrt(epsilon (v_i^T)^2 m_i / n_0 q_i^2)
		-- the 2014 Abgrall, Kumar model uses vacuum permittivity, whereas Wiki says material permittivity works
		{name='ionDebyeLength', value=1},
		-- normalized ion Debye length: 
		-- lambdaHat_d = lambda_d / l_r
	}:append(fluids:map(function(fluid)
		return table{
			{name='min_'..fluid..'_rho', value=1e-4},
			{name='min_'..fluid..'_P', value=1e-4},
		}
	end):unpack()))
end

function TwoFluidEMHD:getCommonFuncCode()
	return template([[
real3 calc_EField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	return real3_real_mul(U.D, 1. / eps);
}
 
real3 calc_HField(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) { 
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	return real3_real_mul(U.B, 1. / mu);
}

real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
real calc_HTotal(real P, real ETotal) { return P + ETotal; }
real calc_hTotal(constant <?=solver.solver_t?>* solver, real rho, real P, real ETotal) { return calc_HTotal(P, ETotal) / rho; }

real calc_rho_from_U(<?=eqn.cons_t?> U) {
	real rho = 0.;
<? for _,fluid in ipairs(fluids) do 
?>	rho += U.<?=fluid?>_rho;
<? end 
?>	return rho;
}

real calc_rho_from_W(<?=eqn.prim_t?> W) {
	real rho = 0.;
<? for _,fluid in ipairs(fluids) do 
?>	rho += W.<?=fluid?>_rho;
<? end 
?>	return rho;
}

real calc_EPot(<?=eqn.cons_t?> U) { return calc_rho_from_U(U) * U.ePot; }
real calc_EPot_from_W(<?=eqn.prim_t?> W) { return calc_rho_from_W(W) * W.ePot; }

<? for _,fluid in ipairs(fluids) do ?>
real calc_<?=fluid?>_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.<?=fluid?>_v, x); }
real calc_<?=fluid?>_EKin(<?=eqn.prim_t?> W, real3 x) { return W.<?=fluid?>_rho * calc_<?=fluid?>_eKin(W, x); }
real calc_<?=fluid?>_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.<?=fluid?>_P / (solver->heatCapacityRatio - 1.); }
real calc_<?=fluid?>_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_<?=fluid?>_EInt(solver, W) / W.<?=fluid?>_rho; }
real calc_<?=fluid?>_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.<?=fluid?>_m, x) / U.<?=fluid?>_rho; }
real calc_<?=fluid?>_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return calc_<?=fluid?>_EKin(W, x) + calc_<?=fluid?>_EInt(solver, W);
}
real calc_<?=fluid?>_Cs(constant <?=solver.solver_t?>* solver, const <?=eqn.prim_t?>* W) {
	return sqrt(solver->heatCapacityRatio * W-><?=fluid?>_P / W-><?=fluid?>_rho);
}
<? end ?>

real calc_EM_energy(constant <?=solver.solver_t?>* solver, const global <?=eqn.cons_t?>* U, real3 x) {
	real eps = solver->sqrt_eps * solver->sqrt_eps / unit_C2_s2_per_kg_m3;
	real mu = solver->sqrt_mu * solver->sqrt_mu / unit_kg_m_per_C2;
	*value = .5 * (coordLenSq(U->D, x) / eps + coordLenSq(U->B, x) / mu);
}
]], {
		solver = self.solver,
		eqn = self,
		fluids = fluids,
	})
end

function TwoFluidEMHD:getPrimConsCode()
	return template([[
<?=eqn.prim_t?> primFromCons(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
	<? for _,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_EKin = calc_<?=fluid?>_EKin_fromCons(U, x);
	real <?=fluid?>_EInt = U.<?=fluid?>_ETotal - <?=fluid?>_EKin;
	<? end ?>
	return (<?=eqn.prim_t?>){
		<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = U.<?=fluid?>_rho,
		.<?=fluid?>_v = real3_real_mul(U.<?=fluid?>_m, 1./U.<?=fluid?>_rho),
		.<?=fluid?>_P = (solver->heatCapacityRatio - 1.) * <?=fluid?>_EInt,
		<? end ?>
		.D = U.D,
		.B = U.B,
		.psi = U.psi,
		.phi = U.phi,
		.ePot = U.ePot,
	};
}

<?=eqn.cons_t?> consFromPrim(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_m = real3_real_mul(W.<?=fluid?>_v, W.<?=fluid?>_rho),
		.<?=fluid?>_ETotal = calc_<?=fluid?>_ETotal(solver, W, x),
<? end ?>
		.D = W.D,
		.B = W.B,
		.psi = W.psi,
		.phi = W.phi,
		.ePot = W.ePot,
	};
}

<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
<? for _,fluid in ipairs(fluids) do ?>
	real3 WA_<?=fluid?>_vL = coord_lower(WA.<?=fluid?>_v, x);
<? end ?>
	return (<?=eqn.cons_t?>){
<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_m = real3_add(
			real3_real_mul(WA.<?=fluid?>_v, W.<?=fluid?>_rho), 
			real3_real_mul(W.<?=fluid?>_v, WA.<?=fluid?>_rho)),
		.<?=fluid?>_ETotal = W.<?=fluid?>_rho * .5 * real3_dot(WA.<?=fluid?>_v, WA_<?=fluid?>_vL) 
			+ WA.<?=fluid?>_rho * real3_dot(W.<?=fluid?>_v, WA_<?=fluid?>_vL)
			+ W.<?=fluid?>_P / (solver->heatCapacityRatio - 1.),
<? end ?>
		.B = W.B,
		.D = W.D,
		.phi = W.phi,
		.psi = W.psi,
		.ePot = W.ePot,
	};
}

<?=eqn.prim_t?> apply_dW_dU(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
<? for _,fluid in ipairs(fluids) do ?>
	real3 WA_<?=fluid?>_vL = coord_lower(WA.<?=fluid?>_v, x);
<? end ?>
	return (<?=eqn.prim_t?>){
<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = U.<?=fluid?>_rho,
		.<?=fluid?>_v = real3_sub(
			real3_real_mul(U.<?=fluid?>_m, 1. / WA.<?=fluid?>_rho),
			real3_real_mul(WA.<?=fluid?>_v, U.<?=fluid?>_rho / WA.<?=fluid?>_rho)),
		.<?=fluid?>_P = (solver->heatCapacityRatio - 1.) * (
			.5 * real3_dot(WA.<?=fluid?>_v, WA_<?=fluid?>_vL) * U.<?=fluid?>_rho 
			- real3_dot(U.<?=fluid?>_m, WA_<?=fluid?>_vL)
			+ U.<?=fluid?>_ETotal),
<? end ?>
		.B = U.B,
		.D = U.D,
		.phi = U.phi,
		.psi = U.psi,
		.ePot = U.ePot,
	};
}

]], {
		solver = self.solver,
		eqn = self,
		fluids = fluids,
	})
end

-- overridden because it adds some extra parameters to the template args
-- should I either make a function for the template arg params
-- or maybe I shouldn't have super-class'd the initState code to begin with ...
function TwoFluidEMHD:getInitStateCode()
	return template([[
<? 
local cons_t = eqn.cons_t
local susc_t = eqn.susc_t
local scalar = eqn.scalar
local vec3 = eqn.vec3
local zero = scalar..'_zero'
local inv = scalar..'_inv'
local fromreal = scalar..'_from_real'
local sqrt = scalar..'_sqrt'
?>

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
	real ePot = 0;
<? 
if eqn.useEulerInitState then 
?>
	real rho = 0.;
	real3 v = real3_zero;
	real P = 0;
<?
else
	 for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_rho = 0;
	real3 <?=fluid?>_v = real3_zero;
	real <?=fluid?>_P = 0;
<? 
	end 
end
?>	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=fromreal?>(1.);
	<?=scalar?> permittivity = <?=fromreal?>(1. / (4. * M_PI));
	<?=scalar?> permeability = <?=fromreal?>(4. * M_PI);

	<?=code?>

	<?=eqn.prim_t?> W = {
<? 
if eqn.useEulerInitState then 
?>
		.ion_rho = rho,
		.elec_rho = rho / solver->ionElectronMassRatio, 

		// "the electron pressure is taken to be elec_P = 5 ion_rho"
		// is that arbitrary?
		.elec_P = 5. * rho,
		
		// "the ion pressure is 1/100th the electron pressure"
		// is that from the mass ratio of ion/electron?
		.ion_P = P / solver->ionElectronMassRatio, 

		.ion_v = v,
		.elec_v = v,
	
<?	
else	-- expect the initState to explicitly provide the ion_ and elec_ Euler fluid variables
	for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = <?=fluid?>_rho,
		.<?=fluid?>_v = cartesianToCoord(<?=fluid?>_v, x),
		.<?=fluid?>_P = <?=fluid?>_P,
<?
	end
end
?>
		.D = D, 
		.B = B,
		.psi = 0,
		.phi = 0,
	
		.ePot = 0,
	};
	UBuf[index] = consFromPrim(solver, W, x);
}
]], table({
		solver = self.solver,
		code = self.initState:initState(self.solver),
		fluids = fluids,
	}, self:getTemplateEnv()))
end

TwoFluidEMHD.solverCodeFile = 'eqn/twofluid-emhd.cl'

function TwoFluidEMHD:getTemplateEnv()
	local scalar = self.scalar
	local env = {}
	env.eqn = self
	env.solver = self.solver
	env.vec3 = self.vec3
	env.susc_t = self.susc_t
	env.scalar = scalar
	env.zero = scalar..'_zero'
	env.inv = scalar..'_inv'
	env.neg = scalar..'_neg'
	env.fromreal = scalar..'_from_real'
	env.add = scalar..'_add'
	env.sub = scalar..'_sub'
	env.mul = scalar..'_mul'
	env.mul3 = scalar..'_mul3'
	env.real_mul = scalar..'_real_mul'
	env.sqrt = scalar..'_sqrt'
	env.abs = scalar..'_abs'
	return env
end


TwoFluidEMHD.displayVarCodeUsesPrims = true

TwoFluidEMHD.predefinedDisplayVars = {
	'U ion_rho',
	'U elec_rho',
	'U ion P',
	'U elec P',
	'U D mag',
	'U B mag',
	
	'U ePot',
	'U gravity mag',
}

function TwoFluidEMHD:getDisplayVars()
	local vars = TwoFluidEMHD.super.getDisplayVars(self)

	for _,fluid in ipairs(fluids) do
	
		-- k is 0,1,2
		local function vorticity(k)
			local xs = {'x','y','z'}
			local i = (k+1)%3
			local j = (i+1)%3
			return {name=fluid..' vorticity '..xs[k+1], code=template([[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = U - solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;

		//TODO incorporate metric

		real vim_j = Uim-><?=fluid?>_m.s<?=j?> / Uim-><?=fluid?>_rho;
		real vip_j = Uip-><?=fluid?>_m.s<?=j?> / Uip-><?=fluid?>_rho;
		real vjm_i = Ujm-><?=fluid?>_m.s<?=i?> / Ujm-><?=fluid?>_rho;
		real vjp_i = Ujp-><?=fluid?>_m.s<?=i?> / Ujp-><?=fluid?>_rho;
		
		*value = (vjp_i - vjm_i) / (2. * solver->grid_dx.s<?=i?>)
				- (vip_j - vim_j) / (2. * solver->grid_dx.s<?=j?>);
	}
]], 		{
				i = i,
				j = j,
				eqn = self,
				fluid = fluid,
			})}
		end
		
		vars:append{
			{name=fluid..' v', code='*value_real3 = W.'..fluid..'_v;', type='real3', units='m/s'},
			{name=fluid..' P', code='*value = W.'..fluid..'_P;', units='kg/(m*s^2)'},
			{name=fluid..' eInt', code='*value = calc_'..fluid..'_eInt(solver, W);', units='m^2/s^2'},
			{name=fluid..' eKin', code='*value = calc_'..fluid..'_eKin(W, x);', units='m^2/s^2'},
			{name=fluid..' EInt', code='*value = calc_'..fluid..'_EInt(solver, W);'},
			{name=fluid..' EKin', code='*value = calc_'..fluid..'_EKin(W, x);'},
			{name=fluid..' ETotal', code='*value = U->'..fluid..'_ETotal;'},
			{name=fluid..' S', code='*value = W.'..fluid..'_P / pow(W.'..fluid..'_rho, (real)solver->heatCapacityRatio);'},
			{name=fluid..' H', code='*value = calc_H(solver, W.'..fluid..'_P);'},
			{name=fluid..' h', code='*value = calc_h(solver, W.'..fluid..'_rho, W.'..fluid..'_P);'},
			{name=fluid..' HTotal', code='*value = calc_HTotal(W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{name=fluid..' hTotal', code='*value = calc_hTotal(solver, W.'..fluid..'_rho, W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{name=fluid..' speed of sound', code='*value = calc_'..fluid..'_Cs(solver, &W);', units='m/s'},
			{name=fluid..' Mach number', code='*value = coordLen(W.'..fluid..'_v, x) / calc_'..fluid..'_Cs(solver, &W);'},
		}:append( ({
		-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
		-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
				[1] = {},
				[2] = {vorticity(2)},
				[3] = range(0,2):map(vorticity),

		})[self.solver.dim] )
	end

	vars:append{
		{
			name = 'EField',
			code = template([[	*value_real3 = calc_EField(solver, *U);]], env),
			type = 'real3',
			units = '(kg*m)/(C*s)',
		},
		{
			name = 'HField',
			code = template([[	*value_real3 = calc_HField(solver, *U);]], env),
			type = 'real3',
			units = 'C/(m*s)',
		},
		{
			name = 'SField',
			code = template([[	*value_real3 = real3_cross(calc_EField(solver, *U), calc_HField(solver, *U));]], env), 
			type = 'real3',
			units = 'kg/s^3',
		},
		{
			name = 'gravity',
			code = template([[
	if (!OOB(1,1)) *value_real3 = calcGravityAccel<?=eqn.gravOp.name?>(solver, U, x);
]], {eqn=self}), 
			type='real3', 
			units='m/s^2',
		},
		{
			name = 'EPot',
			code = '*value = calc_rho_from_U(*U) * U->ePot;', 
			units='kg/(m*s^2)',
		},
		{
			name = 'EM energy',
			code = '*value = calc_EM_energy(solver, U, x);',
			units = 'kg/(m*s^2)'
		},
	}:append(table{'D', 'B'}:map(function(field, i)
		local field = assert( ({D='D', B='B'})[field] )
		return {
			name = 'div '..field,
			code = template([[
	*value = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<?
end 
?>	);
]], {solver=self.solver, field=field}),
			units = ({
				D = 'C/m^3',
				B = 'kg/(C*m*s)',
			})[field],	
		}
	end))

	return vars
end

local eigenVars = table()
for _,fluid in ipairs(fluids) do
	eigenVars:append{
		-- Roe-averaged vars
		{name=fluid..'_rho', type='real', units='kg/m^3'},
		{name=fluid..'_v', type='real3', units='m/s'},
		{name=fluid..'_hTotal', type='real', units='m^2/s^2'},

		-- derived vars
		{name=fluid..'_vSq', type='real', units='m^2/s^2'},
		{name=fluid..'_Cs', type='real', units='m/s'},
	}
end

TwoFluidEMHD.eigenVars = eigenVars

function TwoFluidEMHD:eigenWaveCodePrefix(side, eig, x)
	return template([[
<? for i,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_Cs_sqrt_gU = <?=eig?>.<?=fluid?>_Cs * coord_sqrt_g_uu<?=side..side?>(<?=x?>);
	real <?=fluid?>_v_n = <?=eig?>.<?=fluid?>_v.s[<?=side?>];
<? end ?>
]], {
		x = x,
		eig = '('..eig..')',
		fluids = fluids,
		side = side,
	})
end

function TwoFluidEMHD:eigenWaveCode(side, eig, x, waveIndex)
	for i,fluid in ipairs(fluids) do
		if waveIndex == 0 + 5 * (i-1) then
			return template('<?=fluid?>_v_n - <?=fluid?>_Cs_sqrt_gU', {fluid=fluid})
		elseif waveIndex >= 1 + 5 * (i-1)
		and waveIndex <= 3 + 5 * (i-1)
		then
			return template('<?=fluid?>_v_n', {fluid=fluid})
		elseif waveIndex == 4 + 5 * (i-1) then
			return template('<?=fluid?>_v_n + <?=fluid?>_Cs_sqrt_gU', {fluid=fluid})
		end
	end
	if not self.implicitEMIntegration and waveIndex >= 5*#fluids and waveIndex < 5*#fluids+8 then
		-- 2014 Abgrall, Kumar eqn 1.9 says the eigenvalues are c, while the flux contains cHat ...
		return ({
			'-solver->divPhiWavespeed / unit_m_per_s',
			'-solver->divPsiWavespeed / unit_m_per_s',
			'-solver->speedOfLight / unit_m_per_s',
			'-solver->speedOfLight / unit_m_per_s',
			'solver->speedOfLight / unit_m_per_s',
			'solver->speedOfLight / unit_m_per_s',
			'solver->divPsiWavespeed / unit_m_per_s',
			'solver->divPhiWavespeed / unit_m_per_s',
		})[waveIndex - 5*#fluids + 1]
	end
	error('got a bad waveIndex: '..waveIndex)
end

--TODO timestep restriction
-- 2014 Abgrall, Kumar eqn 2.25
-- dt < sqrt( E_alpha,i / rho_alpha,i) * |lHat_r,alpha| sqrt(2) / |E_i + v_alpha,i x B_i|
function TwoFluidEMHD:consWaveCodePrefix(side, U, x)
	return template([[
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);

<? if eqn.implicitEMIntegration then 	--ignoring EM wavespeed	?>	
	real consWaveCode_lambdaMax = -INFINITY;
<? else 								--using EM wavespeed ?>
	real consWaveCode_lambdaMax = max(
			max(
				solver->divPsiWavespeed,
				solver->divPhiWavespeed
			), 
		solver->speedOfLight
	) / unit_m_per_s;
<? end ?>

	real consWaveCode_lambdaMin = -consWaveCode_lambdaMax;

<? for _,fluid in ipairs(eqn.fluids) do
?>	real <?=fluid?>_Cs = calc_<?=fluid?>_Cs(solver, &W);
	real <?=fluid?>_Cs_sqrt_gU = <?=fluid?>_Cs * coord_sqrt_g_uu<?=side..side?>(x);
	consWaveCode_lambdaMin = min(consWaveCode_lambdaMin, W.<?=fluid?>_v.s<?=side?> - <?=fluid?>_Cs_sqrt_gU);
	consWaveCode_lambdaMax = max(consWaveCode_lambdaMax, W.<?=fluid?>_v.s<?=side?> + <?=fluid?>_Cs_sqrt_gU);
<? end
?>

]], {
		eqn = self,
		side = side,
		U = '('..U..')',
		x = x,
	})
end

function TwoFluidEMHD:consMinWaveCode(side, U, x)
	return 'consWaveCode_lambdaMin'
end

function TwoFluidEMHD:consMaxWaveCode(side, U, x)
	return 'consWaveCode_lambdaMax'
end

return TwoFluidEMHD
