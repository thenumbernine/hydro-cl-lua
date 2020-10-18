--[[
2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations"

the hydro/solver/twofluid-emhd-behavior.lua instanciates three separate solvers:
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
local Equation = require 'hydro.eqn.eqn'


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
	{name='ion_rho', type='real', units='kg/m^3', variance=''},
	{name='ion_m', type='real3', units='kg/(m^2*s)', variance='u'},		-- m^i
	{name='ion_ETotal', type='real', units='kg/(m*s^2)', variance=''},
	
	{name='elec_rho', type='real', units='kg/m^3', variance=''},
	{name='elec_m', type='real3', units='kg/(m^2*s)', variance='u'},	-- m^i
	{name='elec_ETotal', type='real', units='kg/(m*s^2)', variance=''},

	{name='D', type='real3', units='C/m^2', variance='l'},				-- D_i
	{name='B', type='real3', units='kg/(C*s)', variance='l'},			-- B_i
	{name='phi', type='real', units='C/m^2'},							-- div D potential
	{name='psi', type='real', units='kg/(C*s)'},						-- div B potential

	--extra	
	{name='ePot', type='real', units='m^2/s^2'},
}

TwoFluidEMHD.primVars = table{
	--integration variables		
	{name='ion_rho', type='real', units='kg/m^3', variance=''},
	{name='ion_v', type='real3', units='m/s', variance='u'},
	{name='ion_P', type='real', units='kg/(m*s^2)', variance=''},
	
	{name='elec_rho', type='real', units='kg/m^3', variance=''},
	{name='elec_v', type='real3', units='m/s', variance='u'},
	{name='elec_P', type='real', units='kg/(m*s^2)', variance=''},

	{name='D', type='real3', units='C/m^2', variance='l'},				-- D_i
	{name='B', type='real3', units='kg/(C*s)', variance='l'},			-- B_i
	{name='phi', type='real', units='C/m^2'},							-- div D potential
	{name='psi', type='real', units='kg/(C*s)'},						-- div B potential
	
	--extra	
	{name='ePot', type='real', units='m^2/s^2'},
}

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
		self.initConds = require 'hydro.init.euler':getList()
	else
		self.initConds = require 'hydro.init.twofluid-emhd':getList()
	end


	TwoFluidEMHD.super.init(self, args)


--	local NoDiv = require 'hydro.op.nodiv'()
--	self.solver.ops:insert(NoDiv{solver=self.solver})	-- nodiv on maxwell ... or just use potentials 

	local TwoFluidSelfGrav = require 'hydro.op.twofluid-selfgrav'
	self.gravOp = TwoFluidSelfGrav{solver=self.solver}
	self.solver.ops:insert(self.gravOp)
end

function TwoFluidEMHD:createInitState()
	TwoFluidEMHD.super.createInitState(self)

	--[[ using the 'units' parameters
	local constants = require 'hydro.constants'
	local speedOfLight = constants.speedOfLight_in_m_per_s
	local vacuumPermeability = constants.vacuumPermeability_in_kg_m_per_C2
	local vacuumPermittivity = constants.vacuumPermittivity_in_C2_s2_per_kg_m3
	local ionMass = constants.protonMass_in_kg
	local ionElectronMassRatio = ionMass / constants.electronMass_in_kg
	local ionChargeMassRatio = constants.protonCharge_in_C / ionMass
	--]]
	-- [[ using paper values
	local speedOfLight = 1
	local vacuumPermittivity = 1
	local vacuumPermeability = 1
	local ionMass = 1
	local ionElectronMassRatio = 100
	local ionChargeMassRatio = 1
	--]]
	
	self:addGuiVars(table{
		--never given, only stated as "speeds for the Maxwell equation"
		-- of course, they're associated with the potentials, so... they could be arbitrary
		{name='divPsiWavespeed', value=speedOfLight, units='m/s'},
		{name='divPhiWavespeed', value=speedOfLight, units='m/s'},
		
		-- gamma = heat capacity ratio
		{name='heatCapacityRatio', value=5/3},

		-- m_i = ion mass
		-- rho is in kg/m^3
		-- so number of ions per m^3 = ion_rho / ionMass
		-- and number of electrons per m^3 = elec_rho / elecMass = elec_rho / (ionMass / ionElectronMassRatio)
		-- but where do we use number density?
		{name='ionMass', value=ionMass, units='kg'},
	
		-- c = speed of light
		{name='speedOfLight', value=speedOfLight, units='m/s'},
		
		-- m = m_i / m_e ~= 1836.15267389
		{name='ionElectronMassRatio', value=ionElectronMassRatio},
		
		-- r_i = q_i / m_i
		{name='ionChargeMassRatio', value=ionChargeMassRatio, units='C/kg'},

		{name='sqrt_mu', value=math.sqrt(vacuumPermittivity), units='(kg m)^.5/C'},
		{name='sqrt_eps', value=math.sqrt(vacuumPermeability), units='(C*s)/(kg*m^3)^.5'},
	
	}:append(self.fluids:map(function(fluid)
		return table{
			{name='min_'..fluid..'_rho', value=1e-4},
			{name='min_'..fluid..'_P', value=1e-4},
		}
	end):unpack()))
end

function TwoFluidEMHD:initCodeModules()
	TwoFluidEMHD.super.initCodeModules(self)

	for moduleName, depends in pairs{
		['sqrt_2_and_1_2'] = {},
	
		['eqn.prim-cons'] = {
			'real3',
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
			'eqn.common',	-- calc_*
		},

		-- only used by PLM
		['eqn.dU-dW'] = {
			'real3',
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
			'coord_lower',
		},

		['eqn.common'] = {
			'coordLenSq',
		},
		
		['fluxFromCons'] = {
			'normal_t',
		},
		
		['eigen_forInterface'] = {},
		['eigen_forCell'] = {},
		
		['eigen_left/rightTransform'] = {
			'sqrt_2_and_1_2',
		},
		
		['eigen_fluxTransform'] = {},
		['addSource'] = {},
		['constrainU'] = {},
		
		--['calcDT'] = {},	-- don't use this, use the deafult, this is breaking atm
	} do
		self:addModuleFromSourceFile{
			name = moduleName,
			depends = depends,
		}
	end
end

-- don't use default
function TwoFluidEMHD:initCodeModule_fluxFromCons() end
function TwoFluidEMHD:initCodeModuleCommon() end
function TwoFluidEMHD:initCodeModulePrimCons() end

function TwoFluidEMHD:getModuleDependsApplyInitCond()
	return table(TwoFluidEMHD.super.getModuleDependsApplyInitCond(self))
	:append{
		'cartesianToCoord',
	}
end

-- this is to be turned into eqn.waveCode for the inline wave calcs
function TwoFluidEMHD:getModuleDependsSolver()
	return {
		'eqn.prim-cons',
		'coord_lower',
		-- but this is for postComputeFluxCode used by the flux stuff
		'coord_sqrt_det_g',
	}
end

TwoFluidEMHD.solverCodeFile = 'hydro/eqn/twofluid-emhd.cl'

function TwoFluidEMHD:getEnv()
	local scalar = self.scalar
	local env = TwoFluidEMHD.super.getEnv(self)
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
	'U ion v',
	'U elec v',
	'U ion P',
	'U elec P',
	'U D',
	'U B',
	
	'U ePot',
	'U gravity',
}

function TwoFluidEMHD:getDisplayVars()
	local vars = TwoFluidEMHD.super.getDisplayVars(self)

	for _,fluid in ipairs(self.fluids) do
	
		vars:append{
			{name=fluid..' v', code='value.vreal3 = W.'..fluid..'_v;', type='real3', units='m/s'},
			{name=fluid..' P', code='value.vreal = W.'..fluid..'_P;', units='kg/(m*s^2)'},
			{name=fluid..' eInt', code='value.vreal = calc_'..fluid..'_eInt(solver, W);', units='m^2/s^2'},
			{name=fluid..' eKin', code='value.vreal = calc_'..fluid..'_eKin(W, x);', units='m^2/s^2'},
			{name=fluid..' EInt', code='value.vreal = calc_'..fluid..'_EInt(solver, W);'},
			{name=fluid..' EKin', code='value.vreal = calc_'..fluid..'_EKin(W, x);'},
			{name=fluid..' ETotal', code='value.vreal = U->'..fluid..'_ETotal;'},
			{name=fluid..' S', code='value.vreal = W.'..fluid..'_P / pow(W.'..fluid..'_rho, (real)solver->heatCapacityRatio);'},
			{name=fluid..' H', code='value.vreal = calc_H(solver, W.'..fluid..'_P);'},
			{name=fluid..' h', code='value.vreal = calc_h(solver, W.'..fluid..'_rho, W.'..fluid..'_P);'},
			{name=fluid..' HTotal', code='value.vreal = calc_HTotal(W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{name=fluid..' hTotal', code='value.vreal = calc_hTotal(solver, W.'..fluid..'_rho, W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{name=fluid..' speed of sound', code='value.vreal = calc_'..fluid..'_Cs(solver, &W);', units='m/s'},
			{name=fluid..' Mach number', code='value.vreal = coordLen(W.'..fluid..'_v, x) / calc_'..fluid..'_Cs(solver, &W);'},
		}
	
		vars:insert(self:createDivDisplayVar{
			field = fluid..' v', 
			getField = function(U, j)
				return U..'->'..fluid..'_m.s'..j..' / '..U..'->'..fluid..'_rho'
			end,
			units = '1/s',
		} or nil)

		vars:insert(self:createCurlDisplayVar{
			field = fluid..' v', 
			getField = function(U, j)
				return U..'->'..fluid..'_m.s'..j..' / '..U..'->'..fluid..'_rho'
			end,
			units = '1/s',
		} or nil)
	end

	vars:append{
		{
			name = 'EField',
			code = self:template[[	value.vreal3 = calc_EField(solver, *U);]],
			type = 'real3',
			units = '(kg*m)/(C*s)',
		},
		{
			name = 'HField',
			code = self:template[[	value.vreal3 = calc_HField(solver, *U);]],
			type = 'real3',
			units = 'C/(m*s)',
		},
		{
			name = 'SField',
			code = self:template[[	value.vreal3 = real3_cross(calc_EField(solver, *U), calc_HField(solver, *U));]], 
			type = 'real3',
			units = 'kg/s^3',
		},
		{
			name = 'gravity',
			code = self:template[[
	if (!OOB(1,1)) value.vreal3 = calcGravityAccel<?=eqn.gravOp.name?>(solver, U, x);
]],
			type='real3', 
			units='m/s^2',
		},
		{
			name = 'EPot',
			code = 'value.vreal = calc_rho_from_U(*U) * U->ePot;', 
			units='kg/(m*s^2)',
		},
		{
			name = 'EM energy',
			code = 'value.vreal = calc_EM_energy(solver, U, x);',
			units = 'kg/(m*s^2)'
		},
	}:append(table{'D', 'B'}:map(function(field, i)
		local field = assert( ({D='D', B='B'})[field] )
		return self:createDivDisplayVar{field=field, units=({
			D = 'C/m^3',
			B = 'kg/(C*m*s)',
		})[field]}
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

function TwoFluidEMHD:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
<? for i,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_Cs_nLen = <?=eig?>.<?=fluid?>_Cs * normal_len(<?=n?>);
	real <?=fluid?>_v_n = normal_vecDotN1(n, <?=eig?>.<?=fluid?>_v);
<? end ?>
]], {
		x = x,
		eig = '('..eig..')',
		fluids = self.fluids,
		n = n,
	})
end

function TwoFluidEMHD:eigenWaveCode(n, eig, x, waveIndex)
	for i,fluid in ipairs(self.fluids) do
		if waveIndex == 0 + 5 * (i-1) then
			return self:template('<?=fluid?>_v_n - <?=fluid?>_Cs_nLen', {fluid=fluid})
		elseif waveIndex >= 1 + 5 * (i-1)
		and waveIndex <= 3 + 5 * (i-1)
		then
			return self:template('<?=fluid?>_v_n', {fluid=fluid})
		elseif waveIndex == 4 + 5 * (i-1) then
			return self:template('<?=fluid?>_v_n + <?=fluid?>_Cs_nLen', {fluid=fluid})
		end
	end
	if not self.implicitEMIntegration and waveIndex >= 5*#self.fluids and waveIndex < 5*#self.fluids+8 then
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
		})[waveIndex - 5*#self.fluids + 1]
	end
	error('got a bad waveIndex: '..waveIndex)
end

--TODO timestep restriction
-- 2014 Abgrall, Kumar eqn 2.25
-- dt < sqrt( E_alpha,i / rho_alpha,i) * |lHat_r,alpha| sqrt(2) / |E_i + v_alpha,i x B_i|
function TwoFluidEMHD:consWaveCodePrefix(n, U, x)
	return self:template([[
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
	real <?=fluid?>_Cs_nLen = <?=fluid?>_Cs * normal_len(<?=n?>);
	consWaveCode_lambdaMin = min(consWaveCode_lambdaMin, normal_vecDotN1(n, W.<?=fluid?>_v) - <?=fluid?>_Cs_nLen);
	consWaveCode_lambdaMax = max(consWaveCode_lambdaMax, normal_vecDotN1(n, W.<?=fluid?>_v) + <?=fluid?>_Cs_nLen);
<? end
?>

]], {
		n = n,
		U = '('..U..')',
		x = x,
	})
end

function TwoFluidEMHD:consMinWaveCode(n, U, x)
	return 'consWaveCode_lambdaMin'
end

function TwoFluidEMHD:consMaxWaveCode(n, U, x)
	return 'consWaveCode_lambdaMax'
end

return TwoFluidEMHD
