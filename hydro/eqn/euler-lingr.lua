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

EulerLinGR.solverCodeFile = 'hydro/eqn/euler-lingr.cl'

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
	
	local speedOfLight = 1
	local gravitationalConstant = 1e-3
	
	self:addGuiVars(table{
		{name='divPsiWavespeed_g', value=speedOfLight, units='m/s'},
		{name='divPhiWavespeed_g', value=speedOfLight, units='m/s'},
		{name='speedOfLight', value=speedOfLight, units='m/s'},
		
		{name='sqrt_G', value=math.sqrt(gravitationalConstant), units='(m^3/(kg*s^2))^.5'},
	})
end

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

function EulerLinGR:getDisplayVars()
	local vars = EulerLinGR.super.getDisplayVars(self)
	local i = assert(vars:find(nil, function(v) return v.name == 'EPot' end))
	vars:remove(i)
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

function EulerLinGR:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
real const Cs_nLen = <?=eig?>->Cs * normal_len(n);
real const v_n = normal_vecDotN1(n, <?=eig?>->v);
]], {
		x = x,
		eig = '('..eig..')',
		n = n,
	})
end

function EulerLinGR:eigenWaveCode(n, eig, x, waveIndex)
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
function EulerLinGR:consWaveCodePrefix(n, U, x)
	return self:template([[
<?=prim_t?> W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=x?>);

<?if true then  -- using the EM wavespeed ?>
real consWaveCode_lambdaMax = max(
		max(
			max(solver->divPsiWavespeed, solver->divPhiWavespeed),
			max(solver->divPsiWavespeed_g, solver->divPhiWavespeed_g)
		),
		solver->speedOfLight
	) / unit_m_per_s;
<? else -- ignoring it ?>
real consWaveCode_lambdaMax = INFINITY;
<? end ?>

real consWaveCode_lambdaMin = -consWaveCode_lambdaMax;

real const Cs = calc_Cs(solver, &W);
real const Cs_nLen = Cs * normal_len(n);
consWaveCode_lambdaMin = min(consWaveCode_lambdaMin, normal_vecDotN1(n, W.v) - Cs_nLen);
consWaveCode_lambdaMax = max(consWaveCode_lambdaMax, normal_vecDotN1(n, W.v) + Cs_nLen);

]], {
		n = n,
		U = '('..U..')',
		x = x,
	})
end

function EulerLinGR:consMinWaveCode(n, U, x)
	return 'consWaveCode_lambdaMin'
end

function EulerLinGR:consMaxWaveCode(n, U, x)
	return 'consWaveCode_lambdaMax'
end

return EulerLinGR
