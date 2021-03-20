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
	
	local speedOfLight = 299792458
	local gravitationalConstant = 6.67408e-11
	
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


	vars:append{
		{
			name = 'gravity',
			type = 'real3',
			units = 'm/s^2',
			code = self:template[[
	value.vreal3 = real3_real_mul(calcGravityForcePerVolume(solver, U, x), 1. / U->rho);
]],
		},
		
		-- how do I recreate ePot?
		-- d/dt E_g = chi phi_g / eps_g
		-- [m/s^3] = [m/s] [m/s^2]
		-- d/dx^i ePot = E_g^i
		-- so ePot = int (chi phi_g / eps_g) dx^i dt
		{
			name = 'phi_g / eps_g',
			units = 'm/s^2',
			code = self:template[[
	real const G = solver->gravitationalConstant / unit_m3_per_kg_s2;
	real const _1_eps_g = 4. * M_PI * G;
	value.vreal = U->phi_g * _1_eps_g;
]],
		},
	}

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

function EulerLinGR:eigenWaveCodePrefix(n, eig, pt)
	return self:template([[
real const Cs_nLen = <?=eig?>->Cs * normal_len(n);
real const v_n = normal_vecDotN1(n, <?=eig?>->v);

real waveCode_lambdaMax = max(
		max(solver->divPsiWavespeed_g, solver->divPhiWavespeed_g),
		solver->speedOfLight
	) / unit_m_per_s;
real waveCode_lambdaMin = -waveCode_lambdaMax;
waveCode_lambdaMin = min(waveCode_lambdaMin, v_n - Cs_nLen);
waveCode_lambdaMax = max(waveCode_lambdaMax, v_n + Cs_nLen);
]], {
		pt = pt,
		eig = '('..eig..')',
		n = n,
	})
end

-- TODO instead of an expression (which requires pre-calc of vars, which requires extra vars to be stored when not needed)
-- how about an arg for the result?
function EulerLinGR:eigenMinWaveCode(n, eig, pt)
	return 'waveCode_lambdaMin'
end
function EulerLinGR:eigenMaxWaveCode(n, eig, pt)
	return 'waveCode_lambdaMax'
end

function EulerLinGR:eigenWaveCode(n, eig, pt, waveIndex)
	if waveIndex == 0 then
		return 'v_n - Cs_nLen'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		return 'v_n'
	elseif waveIndex == 4 then
		return 'v_n + Cs_nLen'
	end
	if waveIndex >= 5 and waveIndex < 5+8 then
		return ({
			'-solver->divPhiWavespeed_g / unit_m_per_s',
			'-solver->divPsiWavespeed_g / unit_m_per_s',
			'-solver->speedOfLight / unit_m_per_s',
			'-solver->speedOfLight / unit_m_per_s',
			'solver->speedOfLight / unit_m_per_s',
			'solver->speedOfLight / unit_m_per_s',
			'solver->divPsiWavespeed_g / unit_m_per_s',
			'solver->divPhiWavespeed_g / unit_m_per_s',
		})[waveIndex - 5 + 1]
	end
	error('got a bad waveIndex: '..waveIndex)
end

--TODO timestep restriction
-- 2014 Abgrall, Kumar eqn 2.25
-- dt < sqrt( E_alpha,i / rho_alpha,i) * |lHat_r,alpha| sqrt(2) / |E_i + v_alpha,i x B_i|
function EulerLinGR:consWaveCodePrefix(n, U, pt)
	return self:template([[
<?=prim_t?> W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=pt?>);

real const Cs = calc_Cs(solver, &W);
real const Cs_nLen = Cs * normal_len(n);
real const v_n = normal_vecDotN1(n, W.v);

real waveCode_lambdaMax = max(
		max(solver->divPsiWavespeed_g, solver->divPhiWavespeed_g),
		solver->speedOfLight
	) / unit_m_per_s;
real waveCode_lambdaMin = -waveCode_lambdaMax;
waveCode_lambdaMin = min(waveCode_lambdaMin, v_n - Cs_nLen);
waveCode_lambdaMax = max(waveCode_lambdaMax, v_n + Cs_nLen);
]], {
		n = n,
		U = '('..U..')',
		pt = pt,
	})
end

function EulerLinGR:consMinWaveCode(n, U, pt)
	return 'waveCode_lambdaMin'
end
function EulerLinGR:consMaxWaveCode(n, U, pt)
	return 'waveCode_lambdaMax'
end

return EulerLinGR
