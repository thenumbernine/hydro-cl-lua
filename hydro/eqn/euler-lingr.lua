--[[
Euler fluid equations (rho, v^i, P) <-> (rho, m^i, ETotal)
with additional GEM (phi_g, A_g)
so that means no need for op/selfgrav because it's now built in as (phi, A) 
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Euler = require 'hydro.eqn.euler'


local EulerLinGR = class(Euler)

EulerLinGR.numWaves = nil
EulerLinGR.numIntStates = nil

EulerLinGR.name = 'euler_lingr'

-- Hmm, setting this to 'true' destroys us ... causes the magnitude of E_g to explode by 1e+16 ... why is using the flux bad?
-- Fwiw, both Euler and Maxwell equations hold true the homogeneity property: dF/dU * U = F
-- which means this flag shouldn't change the results.
EulerLinGR.roeUseFluxFromCons = false

function EulerLinGR:buildVars(args)
	EulerLinGR.super.buildVars(self, args)
	for _,t in ipairs{self.primVars, self.consVars} do
		local i = assert(t:find(nil, function(v) return v.name == 'ePot' end))
		t:remove(i)
	
		-- TODO behaviors of eqn's ... just apply behavior to this and to twofluid-emhd
		t:append{
			{name='D_g', type='real3', units='kg/m^2', variance='l'},			-- (D_g)_i
			{name='B_g', type='real3', units='1/s', variance='l'},				-- (B_g)_i
			{name='phi_g', type='real', units='kg/m^2', variance=''},			-- div D_g potential
			{name='psi_g', type='real', units='1/s', variance=''},				-- div B_g potential
		}
	end
end

function EulerLinGR:buildSelfGrav()
	-- no self-grav.  that's now done through the lingr
end

function EulerLinGR:createInitState()
	EulerLinGR.super.createInitState(self)
	
	local speedOfLight = 1	--299792458
	local gravitationalConstant = 1	--6.67408e-11
	
	self:addGuiVars(table{
		{name='divPsiWavespeed_g', value=speedOfLight, units='m/s'},
		{name='divPhiWavespeed_g', value=speedOfLight, units='m/s'},
		{name='speedOfLight', value=speedOfLight, units='m/s'},
		{name='gravitationalConstant', value=gravitationalConstant, units='m^3/(kg*s^2)'},
	})
end

-- Euler has its own calcDT to save a sqrt calculation per side
-- but I'm being lazy here, I could write my own for EulerLinGR but I'll just use the default (and repeat sqrt() calcs)
function EulerLinGR:initCodeModule_calcDTCell()
	local file = require 'ext.file'
	-- ok the current implementation - because it uses 'addFromMarkup', and because it is a macro
	--  there's no way to just add a new dependency module 
	--return EulerLinGR.super.super.initCodeModule_calcDTCell(self)
	-- so here.  it's ugly.
	self.solver.modules:addFromMarkup(self:template(file['hydro/eqn/cl/calcDT.cl'])..self:template[[
//// MODULE_DEPENDS: <?=primFromCons?>
]])
end

-- only matters for non-ident metric vector components ... which I'm avoiding
-- also I haven't fully thought through the math for the Maxwell flux Jacobian, esp in context of extrinsic curvature
--[=[
function EulerLinGR:postComputeFluxCode()
	return self:template[[
//// MODULE_DEPENDS: <?=coord_sqrt_det_g?> <?=coord_lower?>
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
		flux.D = real3_real_mul(coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = real3_real_mul(coord_lower(flux.B, x), _1_sqrt_det_g);
		flux.D_g = real3_real_mul(coord_lower(flux.D_g, x), _1_sqrt_det_g);
		flux.B_g = real3_real_mul(coord_lower(flux.B_g, x), _1_sqrt_det_g);
]]
end
--]=]

-- don't use default
function EulerLinGR:initCodeModule_fluxFromCons() end
function EulerLinGR:initCodeModule_consFromPrim_primFromCons() end

EulerLinGR.solverCodeFile = 'hydro/eqn/euler-lingr.cl'

EulerLinGR.predefinedDisplayVars = {
	'U rho',
	'U v',
	'U v z',
	'U P',
	'U D_g',
	'U D_g z',
	'U B_g',
	'U B_g z',
	'U E_g',
	'U E_g z',
	'U E_g mag',
	'U phi_g',
	'U psi_g',
	'U v cross B_g',
	'U v cross B_g z',
	'U GEM energy',
}

function EulerLinGR:getDisplayVars()
	local vars = EulerLinGR.super.getDisplayVars(self)
	
	local i = assert(vars:find(nil, function(v) return v.name == 'EPot' end))
	vars:remove(i)

	vars:append{
		{
			name = 'E_g',
			type = 'real3',
			units = 'm/s^2',
			code = 'value.vreal3 = calc_EgField(solver, U);',
		},
		{
			name = 'H_g',
			type = 'real3',
			units = 'kg/(m*s)',
			code = 'value.vreal3 = calc_HgField(solver, U);',
		},
		{
			name = 'S_g',	-- S Poynting, not S entropy
			type = 'real3',
			units = 'kg/s^3',
			code = 'value.vreal3 = calc_SgField(solver, U);', 
		},
		{
			name = 'GEM energy', 
			units = 'kg/(m*s^2)',
			code = 'value.vreal = calc_GEM_energy(solver, U, x);',
		},
		
		-- D_g / eps = E_g = Newtonian gravitational acceleration
		-- v cross B_g = gravito-magnetic acceleration
		{
			name = 'v cross B_g',
			type = 'real3',
			units = 'm/s^2',
			code = 'value.vreal3 = real3_cross(W.v, U->B_g);',
		},
	}:append(table{'D_g','B_g'}:map(function(field,i)
		local field = assert( ({D_g='D_g', B_g='B_g'})[field] )
		return self:createDivDisplayVar{
			field = field,
			units = ({
				D_g = 'kg/m^3',
				B_g = '1/(m*s)',
			})[field],
		}
	end))

	vars:insert{
		name = 'gravity',
		type = 'real3',
		units = 'm/s^2',
		code = self:template[[
value.vreal3 = real3_real_mul(calcGravityForcePerVolume(solver, U, x), 1. / U->rho);
]],
	}
		
	-- how do I recreate ePot?
	-- d/dt E_g = chi phi_g / eps_g
	-- [m/s^3] = [m/s] [m/s^2]
	-- d/dx^i ePot = E_g^i
	-- so ePot = int (chi phi_g / eps_g) dx^i dt
	vars:insert{
		name = 'phi_g / eps_g',
		units = 'm/s^2',
		code = self:template[[
real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;
real const _1_eps_g = 4. * M_PI * G;
value.vreal = U->phi_g * _1_eps_g;
]],
	}

	-- add h_ab display based on phi and A
	
	vars:insert{
		name = 'h_tt',
		code = self:template[[
real const c = solver->speedOfLight / unit_m_per_s;
real const h_tt = -1. + 2 * U->phi_g / (c*c);
value.vreal = h_tt;
]],
	}

	vars:insert{
		name = 'h_tt',
		code = self:template[[
real const c = solver->speedOfLight / unit_m_per_s;
real const h_tt = 1. + 2 * U->phi_g / (c*c);
value.vreal = h_tt;
]],
	}

	--[[
	[R] = [G_uv] = 1/m^2
	[G] = m^3/(kg*s^2) = m/kg * m^2/s^2
	[c] = m/s
	[G/c^4] = s^2/(kg*m)
	[G_uv] = [8 pi G/c^4 T_uv] = 1/m^2
	so [T_uv] * s^2/(kg*m) = 1/m^2
	so [T_uv] = kg/(m*s^2) ... which is J / m^3 = energy / volume
	so T_ab = rho c^2 + P
	and [T_00] = kg/m^3 m^2/s^2 = kg/(m*s^2)

	Δ A_i = 4 pi G T_0i / c^3
	T_0i = 
	A_i = 4 pi G / c^3 Δ^-1 T_0i
	
	T_00 = c^2 ρ 
	Δ Φ = -4 π G T_00 / c^2 = -4 π G ρ
	∇.Φ = E
	∇.E = -4 π G ρ
	B = ∇ × A
	∇.B = Ψ = 0
	Δ B = ∇ Ψ = 0
	--]]
	return vars
end

local eigenVars = table()
eigenVars:append{
	-- Roe-averaged vars
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	{name='hTotal', type='real', units='m^2/s^2'},

	-- derived vars
	{name='vSq', type='real', units='m^2/s^2'},
	{name='vL', type='real3', units='m/s'},
	{name='Cs', type='real', units='m/s'},
}
EulerLinGR.eigenVars = eigenVars

function EulerLinGR:eigenWaveCodePrefix(args)
	return self:template([[
real const Cs_nLen = (<?=eig?>)->Cs * normal_len(n);
real const v_n = normal_vecDotN1(n, (<?=eig?>)->v);
]], args)
end

function EulerLinGR:eigenWaveCode(args)
	if args.waveIndex == 0 then
		return 'v_n - Cs_nLen'
	elseif args.waveIndex >= 1 and args.waveIndex <= 3 then
		return 'v_n'
	elseif args.waveIndex == 4 then
		return 'v_n + Cs_nLen'
	end
	if args.waveIndex >= 5 and args.waveIndex < 5+8 then
		return ({
			'-solver->divPhiWavespeed_g / unit_m_per_s',
			'-solver->divPsiWavespeed_g / unit_m_per_s',
			'-solver->speedOfLight / unit_m_per_s',
			'-solver->speedOfLight / unit_m_per_s',
			'solver->speedOfLight / unit_m_per_s',
			'solver->speedOfLight / unit_m_per_s',
			'solver->divPsiWavespeed_g / unit_m_per_s',
			'solver->divPhiWavespeed_g / unit_m_per_s',
		})[args.waveIndex - 5 + 1]
	end
	error('got a bad waveIndex: '..args.waveIndex)
end

--TODO timestep restriction
-- 2014 Abgrall, Kumar eqn 2.25
-- dt < sqrt( E_alpha,i / rho_alpha,i) * |lHat_r,alpha| sqrt(2) / |E_i + v_alpha,i x B_i|
function EulerLinGR:consWaveCodePrefix(args)
	return self:template([[
real const Cs_nLen = <?=calc_Cs_fromCons?>(solver, <?=U?>, <?=pt?>) * normal_len(n);
real const v_n = normal_vecDotN1(n, (<?=U?>)->m) / (<?=U?>)->rho;
]], args)
end

-- you can do this so long as no code uses U/eig ... if it does then you have to reassign args
EulerLinGR.consWaveCode = EulerLinGR.eigenWaveCode

function EulerLinGR:eigenWaveCodeMinMax(args)
	return self:template([[
real const Cs_nLen = (<?=eig?>)->Cs * normal_len(n);
real const v_n = normal_vecDotN1(n, (<?=eig?>)->v);

real const waveCode_lambdaMax = max(
		max(solver->divPsiWavespeed_g, solver->divPhiWavespeed_g),
		solver->speedOfLight
	) / unit_m_per_s;

<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'min(-waveCode_lambdaMax, v_n - Cs_nLen)',
	'max( waveCode_lambdaMax, v_n + Cs_nLen)'
)?>
]], args)
end

function EulerLinGR:consWaveCodeMinMax(args)
	return self:template([[
real const Cs_nLen = <?=calc_Cs_fromCons?>(solver, <?=U?>, <?=pt?>) * normal_len(n);
real const v_n = normal_vecDotN1(n, (<?=U?>)->m) / (<?=U?>)->rho;

real const waveCode_lambdaMax = max(
		max(solver->divPsiWavespeed_g, solver->divPhiWavespeed_g),
		solver->speedOfLight
	) / unit_m_per_s;

<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'min(-waveCode_lambdaMax, v_n - Cs_nLen)',
	'max( waveCode_lambdaMax, v_n + Cs_nLen)'
)?>
]], args)
end

function EulerLinGR:consWaveCodeMinMaxAllSidesPrefix(args)
	return self:template([[
real const Cs = <?=calc_Cs_fromCons?>(solver, <?=U?>, <?=pt?>);\
real const waveCode_lambdaMax = max(
		max(solver->divPsiWavespeed_g, solver->divPhiWavespeed_g),
		solver->speedOfLight
	) / unit_m_per_s;
]],	args)
end

function EulerLinGR:consWaveCodeMinMaxAllSides(args)
	return self:template([[
real const Cs_nLen = Cs * normal_len(<?=n?>);
real const v_n = normal_vecDotN1(<?=n?>, (<?=U?>)->m) / (<?=U?>)->rho;

<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'min(-waveCode_lambdaMax, v_n - Cs_nLen)',
	'max( waveCode_lambdaMax, v_n + Cs_nLen)'
)?>
]], args)
end

return EulerLinGR
