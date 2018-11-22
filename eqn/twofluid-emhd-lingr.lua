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
 then in the post-flux code, transform E and B by g_ij
 and in the source terms, add to ion_ and elec_ m^i connection values

--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local template = require 'template'


local TwoFluidEMHDDeDonderGaugeLinearizedGR = class(Equation)

local fluids = table{'ion', 'elec'}
TwoFluidEMHDDeDonderGaugeLinearizedGR.fluids = fluids

TwoFluidEMHDDeDonderGaugeLinearizedGR.postComputeFluxCode = [[
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
		flux.E = real3_real_mul(coord_lower(flux.E, x), _1_sqrt_det_g);
		flux.B = real3_real_mul(coord_lower(flux.B, x), _1_sqrt_det_g);
		flux.E_g = real3_real_mul(coord_lower(flux.E_g, x), _1_sqrt_det_g);
		flux.B_g = real3_real_mul(coord_lower(flux.B_g, x), _1_sqrt_det_g);
]]

TwoFluidEMHDDeDonderGaugeLinearizedGR.name = 'TwoFluidEMHDDeDonderGaugeLinearizedGR'
TwoFluidEMHDDeDonderGaugeLinearizedGR.numWaves = 26
TwoFluidEMHDDeDonderGaugeLinearizedGR.numIntStates = 26

TwoFluidEMHDDeDonderGaugeLinearizedGR.consVars = table{
	--integration variables		
	{ion_rho = 'real'},
	{ion_m = 'real3'},		-- m^i
	{ion_ETotal = 'real'},
	
	{elec_rho = 'real'},
	{elec_m = 'real3'},		-- m^i
	{elec_ETotal = 'real'},

	{E = 'real3'},			-- E_i
	{B = 'real3'},			-- B_i
	{phi = 'real'},			-- E potential
	{psi = 'real'},			-- B potential

	{E_g = 'real3'},		-- (E_g)_i
	{B_g = 'real3'},		-- (B_g)_i
	{phi_g = 'real'},		-- E_g potential
	{psi_g = 'real'},		-- B_g potential
}

TwoFluidEMHDDeDonderGaugeLinearizedGR.primVars = table{
	--integration variables		
	{ion_rho = 'real'},
	{ion_v = 'real3'},
	{ion_P = 'real'},
	
	{elec_rho = 'real'},
	{elec_v = 'real3'},
	{elec_P = 'real'},

	{E = 'real3'},
	{B = 'real3'},
	{phi = 'real'},
	{psi = 'real'},
	
	{E_g = 'real3'},
	{B_g = 'real3'},
	{phi_g = 'real'},
	{psi_g = 'real'},
}

TwoFluidEMHDDeDonderGaugeLinearizedGR.mirrorVars = {
	{'ion_m.x', 'elec_m.x', 'E.x', 'B.x', 'E_g.x', 'B_g.x'}, 
	{'ion_m.y', 'elec_m.y', 'E.y', 'B.y', 'E_g.y', 'B_g.y'}, 
	{'ion_m.z', 'elec_m.z', 'E.z', 'B.z', 'E_g.z', 'B_g.z'},
}

TwoFluidEMHDDeDonderGaugeLinearizedGR.hasEigenCode = true
TwoFluidEMHDDeDonderGaugeLinearizedGR.hasFluxFromConsCode = true
TwoFluidEMHDDeDonderGaugeLinearizedGR.roeUseFluxFromCons = true
TwoFluidEMHDDeDonderGaugeLinearizedGR.useSourceTerm = true
TwoFluidEMHDDeDonderGaugeLinearizedGR.useConstrainU = true

TwoFluidEMHDDeDonderGaugeLinearizedGR.useEulerInitState = true	-- default to true.

-- TODO this has the symptoms of all the other CL kernels that intel's compiler bugged out on
-- i.e. takes a minute to run a kernel on the GPU
function TwoFluidEMHDDeDonderGaugeLinearizedGR:init(args)
	
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

	TwoFluidEMHDDeDonderGaugeLinearizedGR.super.init(self, args)
end

function TwoFluidEMHDDeDonderGaugeLinearizedGR:createInitState()
	TwoFluidEMHDDeDonderGaugeLinearizedGR.super.createInitState(self)
	
	self:addGuiVars(table{
		--never given, only stated as "speeds for the Maxwell equation"
		-- of course, they're associated with the potentials, so... they could be arbitrary
		{name='divPsiWavespeed', value=1},
		{name='divPhiWavespeed', value=1},
		
		{name='divPsiGWavespeed', value=1},
		{name='divPhiGWavespeed', value=1},
		
		-- gamma = heat capacity ratio
		{name='heatCapacityRatio', value=5/3},

		-- x_0 = reference length
		{name='referenceLength', value=1},

		-- B_0 = reference magnetic field
		{name='referenceMagneticField', value=1},

		-- v_i^T = ion reference thermal velocity
		--{name='ionReferenceThermalVelocity', value=1},
		-- normalized Larmor radius:
		-- lHat = l_r / x_0 = gamma m_i v_i^T / (q_i B_0 x_0)	
		-- the 2014 Abgrall, Kumar model has no gamma, but Wiki says to add gamma = Lorentz boost
		-- therefore:
		-- v_i^T = l_r B_0 q_i / m_i
		-- just use this equation:
		-- v_i^T = l_r r_i B_0

		-- m_i = ion mass
		{name='ionMass', value=1},
	
		-- c = speed of light
		{name='speedOfLight', value=1},
		-- normalized speed of light:
		-- cHat = c / v_i^T
		
		-- l_r = larmor radius
		-- in the experimental section of 2014 Abgrall Kumar this is given fixed values
		{name='ionLarmorRadius', value=.1},
		
		-- m = m_i / m_e
		-- https://en.wikipedia.org/wiki/Proton-to-electron_mass_ratio
		-- m_i / m_e = 1836.15267389
		{name='ionElectronMassRatio', value=100},
		
		-- r_i = q_i / m_i
		{name='ionChargeMassRatio', value=1},
	
		-- lambda_d = ion Debye length
		-- lambda_d = sqrt(epsilon (v_i^T)^2 / n_0 q_i)
		-- the 2014 Abgrall, Kumar model uses vacuum permittivity, whereas Wiki says material permittivity works
		{name='ionDebyeLength', value=1},
		-- normalized ion Debye length: 
		-- lambdaHat_d = lambda_d / l_r

		{name='eps', value=1},			-- permittivity
		{name='mu', value=1},			-- permeability
		{name='gravitationalConstant', value=1},		-- 6.67408e-11 m^3 / (kg s^2)
	
	}:append(fluids:map(function(fluid)
		return table{
			{name='min_'..fluid..'_rho', value=1e-4},
			{name='min_'..fluid..'_P', value=1e-4},
		}
	end):unpack()))
end

function TwoFluidEMHDDeDonderGaugeLinearizedGR:getCommonFuncCode()
	return template([[
real ESq(<?=eqn.cons_t?> U, real3 x) { return coordLenSq(U.E, x); }
real BSq(<?=eqn.cons_t?> U, real3 x) { return coordLenSq(U.B, x); }

inline real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
inline real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
inline real calc_hTotal(constant <?=solver.solver_t?>* solver, real rho, real P, real ETotal) { return (P + ETotal) / rho; }
inline real calc_HTotal(real P, real ETotal) { return P + ETotal; }

<? for _,fluid in ipairs(fluids) do ?>
inline real calc_<?=fluid?>_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.<?=fluid?>_v, x); }
inline real calc_<?=fluid?>_EKin(<?=eqn.prim_t?> W, real3 x) { return W.<?=fluid?>_rho * calc_<?=fluid?>_eKin(W, x); }
inline real calc_<?=fluid?>_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.<?=fluid?>_P / (solver->heatCapacityRatio - 1.); }
inline real calc_<?=fluid?>_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_<?=fluid?>_EInt(solver, W) / W.<?=fluid?>_rho; }
inline real calc_<?=fluid?>_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.<?=fluid?>_m, x) / U.<?=fluid?>_rho; }
inline real calc_<?=fluid?>_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return calc_<?=fluid?>_EKin(W, x) + calc_<?=fluid?>_EInt(solver, W);
}
inline real calc_<?=fluid?>_Cs(constant <?=solver.solver_t?>* solver, const <?=eqn.prim_t?>* W) {
	return sqrt(solver->heatCapacityRatio * W-><?=fluid?>_P / W-><?=fluid?>_rho);
}
<? end ?>

inline real calc_EM_energy(constant <?=solver.solver_t?>* solver, const global <?=eqn.cons_t?>* U, real3 x) {
	return .5 * (solver->eps * coordLenSq(U->E, x) + coordLenSq(U->B, x) / solver->mu);
}
]], {
		solver = self.solver,
		eqn = self,
		fluids = fluids,
	})
end

function TwoFluidEMHDDeDonderGaugeLinearizedGR:getPrimConsCode()
	return template([[
inline <?=eqn.prim_t?> primFromCons(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
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
		.E = U.E,
		.B = U.B,
		.psi = U.psi,
		.phi = U.phi,
		.E_g = U.E_g,
		.B_g = U.B_g,
		.psi_g = U.psi_g,
		.phi_g = U.phi_g,
	};
}

inline <?=eqn.cons_t?> consFromPrim(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_m = real3_real_mul(W.<?=fluid?>_v, W.<?=fluid?>_rho),
		.<?=fluid?>_ETotal = calc_<?=fluid?>_ETotal(solver, W, x),
<? end ?>
		.E = W.E,
		.B = W.B,
		.psi = W.psi,
		.phi = W.phi,
		.E_g = W.E_g,
		.B_g = W.B_g,
		.psi_g = W.psi_g,
		.phi_g = W.phi_g,
	};
}

inline <?=eqn.cons_t?> apply_dU_dW(
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
		.E = W.E,
		.phi = W.phi,
		.psi = W.psi,
		.B_g = W.B_g,
		.E_g = W.E_g,
		.phi_g = W.phi_g,
		.psi_g = W.psi_g,
	};
}

inline <?=eqn.prim_t?> apply_dW_dU(
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
		.E = U.E,
		.phi = U.phi,
		.psi = U.psi,
		.B_g = U.B_g,
		.E_g = U.E_g,
		.phi_g = U.phi_g,
		.psi_g = U.psi_g,
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
function TwoFluidEMHDDeDonderGaugeLinearizedGR:getInitStateCode()
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
<? 
if eqn.useEulerInitState then 
?>
	real rho = 0.;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
<?
else
	 for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_rho = 0;
	real3 <?=fluid?>_v = real3_zero;
	real <?=fluid?>_P = 0;
	real <?=fluid?>_ePot = 0;
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
		.E = <?=vec3?>_<?=scalar?>_mul(D, <?=inv?>(permittivity)),
		.B = B,
		.psi = 0,
		.phi = 0,
		.E_g = <?=vec3?>_zero,
		.B_g = <?=vec3?>_zero,
		.psi_g = 0,
		.phi_g = 0,
	};
	UBuf[index] = consFromPrim(solver, W, x);
}
]], table({
		solver = self.solver,
		code = self.initState:initState(self.solver),
		fluids = fluids,
	}, self:getTemplateEnv()))
end

TwoFluidEMHDDeDonderGaugeLinearizedGR.solverCodeFile = 'eqn/twofluid-emhd-lingr.cl'

function TwoFluidEMHDDeDonderGaugeLinearizedGR:getTemplateEnv()
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


TwoFluidEMHDDeDonderGaugeLinearizedGR.displayVarCodeUsesPrims = true

TwoFluidEMHDDeDonderGaugeLinearizedGR.predefinedDisplayVars = {
	'U ion_rho',
	'U elec_rho',
	'U E mag',
	'U B mag',
	'U psi',
	'U phi',
	'U E_g mag',
	'U B_g mag',
	'U psi_g',
	'U phi_g',
	'U EM energy',
}

function TwoFluidEMHDDeDonderGaugeLinearizedGR:getDisplayVars()
	local vars = TwoFluidEMHDDeDonderGaugeLinearizedGR.super.getDisplayVars(self)

	for _,fluid in ipairs(fluids) do
	
		-- k is 0,1,2
		local function vorticity(k)
			local xs = {'x','y','z'}
			local i = (k+1)%3
			local j = (i+1)%3
			return {[fluid..' vorticity '..xs[k+1]] = template([[
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
			{[fluid..' v'] = '*value_real3 = W.'..fluid..'_v;', type='real3'},
			{[fluid..' P'] = '*value = W.'..fluid..'_P;'},
			{[fluid..' eInt'] = '*value = calc_'..fluid..'_eInt(solver, W);'},
			{[fluid..' eKin'] = '*value = calc_'..fluid..'_eKin(W, x);'},
			{[fluid..' EInt'] = '*value = calc_'..fluid..'_EInt(solver, W);'},
			{[fluid..' EKin'] = '*value = calc_'..fluid..'_EKin(W, x);'},
			{[fluid..' ETotal'] = '*value = U->'..fluid..'_ETotal;'},
			{[fluid..' S'] = '*value = W.'..fluid..'_P / pow(W.'..fluid..'_rho, (real)solver->heatCapacityRatio);'},
			{[fluid..' H'] = '*value = calc_H(solver, W.'..fluid..'_P);'},
			{[fluid..' h'] = '*value = calc_h(solver, W.'..fluid..'_rho, W.'..fluid..'_P);'},
			{[fluid..' HTotal'] = '*value = calc_HTotal(W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{[fluid..' hTotal'] = '*value = calc_hTotal(solver, W.'..fluid..'_rho, W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{[fluid..'Speed of Sound'] = '*value = calc_'..fluid..'_Cs(solver, &W);'},
			{[fluid..'Mach number'] = '*value = coordLen(W.'..fluid..'_v, x) / calc_'..fluid..'_Cs(solver, &W);'},
		}:append( ({
		-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
		-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
				[1] = {},
				[2] = {vorticity(2)},
				[3] = range(0,2):map(vorticity),

		})[self.solver.dim] )
	end

	vars:append{
		{['EM energy'] = [[
	*value = calc_EM_energy(solver, U, x);
]]},
	}:append(table{'E','B'}:map(function(var,i)
		local field = assert( ({E='E', B='B'})[var] )
		return {['div '..var] = template([[
	*value = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<?
end 
?>	);
]], {solver=self.solver, field=field})}
	end))

	return vars
end

local eigenVars = table()
for _,fluid in ipairs(fluids) do
	eigenVars:append{
		-- Roe-averaged vars
		{[fluid..'_rho'] = 'real'},
		{[fluid..'_v'] = 'real3'},
		{[fluid..'_hTotal'] = 'real'},

		-- derived vars
		{[fluid..'_vSq'] = 'real'},
		{[fluid..'_Cs'] = 'real'},
	}
end

TwoFluidEMHDDeDonderGaugeLinearizedGR.eigenVars = eigenVars

function TwoFluidEMHDDeDonderGaugeLinearizedGR:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real _1_sqrt_det_g = 1. / sqrt_det_g_grid(cell_x(i));
<? for i,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_Cs_sqrt_gU = <?=eig?>.<?=fluid?>_Cs * coord_sqrt_gU<?=side..side?>(<?=x?>);
	real <?=fluid?>_v_n = <?=eig?>.<?=fluid?>_v.s[<?=side?>];
<? end ?>
]], {
		x = x,
		eig = '('..eig..')',
		fluids = fluids,
		side = side,
	})
end

function TwoFluidEMHDDeDonderGaugeLinearizedGR:eigenWaveCode(side, eig, x, waveIndex)
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
	if waveIndex >= 5*#fluids and waveIndex < 5*#fluids+8 then
		-- 2014 Abgrall, Kumar eqn 1.9 says the eigenvalues are c, while the flux contains cHat ...
		return ({
			'-normalizedSpeedOfLight * solver->divPhiWavespeed * _1_sqrt_det_g',
			'-normalizedSpeedOfLight * solver->divPsiWavespeed * _1_sqrt_det_g',
			'-normalizedSpeedOfLight * _1_sqrt_det_g',
			'-normalizedSpeedOfLight * _1_sqrt_det_g',
			'normalizedSpeedOfLight * _1_sqrt_det_g',
			'normalizedSpeedOfLight * _1_sqrt_det_g',
			'normalizedSpeedOfLight * solver->divPsiWavespeed * _1_sqrt_det_g',
			'normalizedSpeedOfLight * solver->divPhiWavespeed * _1_sqrt_det_g',
		})[waveIndex - 5*#fluids + 1]
	end
	if waveIndex >= 5*#fluids+8 and waveIndex < 5*#fluids+16 then
		return ({
			'-normalizedSpeedOfLight * solver->divPhiGWavespeed * _1_sqrt_det_g',
			'-normalizedSpeedOfLight * solver->divPsiGWavespeed * _1_sqrt_det_g',
			'-normalizedSpeedOfLight * _1_sqrt_det_g',
			'-normalizedSpeedOfLight * _1_sqrt_det_g',
			'normalizedSpeedOfLight * _1_sqrt_det_g',
			'normalizedSpeedOfLight * _1_sqrt_det_g',
			'normalizedSpeedOfLight * solver->divPsiGWavespeed * _1_sqrt_det_g',
			'normalizedSpeedOfLight * solver->divPhiGWavespeed * _1_sqrt_det_g',
		})[waveIndex - 5*#fluids - 8 + 1]
	end
	error('got a bad waveIndex: '..waveIndex)
end

--TODO timestep restriction
-- 2014 Abgrall, Kumar eqn 2.25
-- dt < sqrt( E_alpha,i / rho_alpha,i) * |lHat_r,alpha| sqrt(2) / |E_i + v_alpha,i x B_i|
function TwoFluidEMHDDeDonderGaugeLinearizedGR:consWaveCodePrefix(side, U, x)
	return template([[
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);

#if 1	//using the EM wavespeed
	real consWaveCode_lambdaMax = max(
		//max(
			max(solver->divPsiWavespeed, solver->divPhiWavespeed),
			//max(solver->divPsiGWavespeed, solver->divPhiGWavespeed)
		//),
		1.) * normalizedSpeedOfLight;
#else	//ignoring it
	real consWaveCode_lambdaMax = INFINITY;
#endif
	real consWaveCode_lambdaMin = -consWaveCode_lambdaMax;

<? for _,fluid in ipairs(eqn.fluids) do
?>	real <?=fluid?>_Cs = calc_<?=fluid?>_Cs(solver, &W);
	real <?=fluid?>_Cs_sqrt_gU = <?=fluid?>_Cs * coord_sqrt_gU<?=side..side?>(x);
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

function TwoFluidEMHDDeDonderGaugeLinearizedGR:consMinWaveCode(side, U, x)
	return 'consWaveCode_lambdaMin'
end

function TwoFluidEMHDDeDonderGaugeLinearizedGR:consMaxWaveCode(side, U, x)
	return 'consWaveCode_lambdaMax'
end

return TwoFluidEMHDDeDonderGaugeLinearizedGR
