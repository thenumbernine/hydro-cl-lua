--[[ 
2000 Munz
But then I changed the eigendecomposition to match my Trangenstein one with epsilons and mus
And then I changed epsilon E -> D, so it's based on B and D 
And I rescaled the wavespeeds and potential state variables so that there were no speed-of-lights (or, for varying susceptibility, phase velocities) that appeared

Work is in my symmath-lua repo tests/output/Maxwell.html file and in my MathWorksheets repo 'Maxwell equations in hyperbolic form.html' file.

I'm still not completely sold on this technique.
Relativistically we are doing:
~F^uv = F^uv + eta^uv phi chi, ~F^uv_,v = mu j^u
~*F^uv = *F^uv + eta^uv psi gamma, ~*F^uv_,v = 0
except with the special twist of changing the coefficients next to the eta^tt term from chi or gamma to 1/chi or 1/gamma
I'm still suspicious on this special exception.
Not doing so means our wavespeeds are all the speed of light, and we can't control the propagation of GLM variables.

--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'

local GLM_Maxwell = class(Equation)
GLM_Maxwell.name = 'glm_maxwell'

GLM_Maxwell.roeUseFluxFromCons = true

-- don't incorporate the Conn^k_ij E_k terms into the flux
GLM_Maxwell.weightFluxByGridVolume = false

GLM_Maxwell.initConds = require 'hydro.init.euler':getList()

function GLM_Maxwell:init(args)
	self.scalar = 'real'
	--self.scalar = 'cplx'
	
	self.vec3 = self.scalar..'3'
	self.mat3x3 = self.scalar..'3x3'

	-- TODO tensor susceptibilty support ... but that affects the eigendecomposition ...
	self.susc_t = self.scalar
	--self.susc_t = self.mat3x3

	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'

	self.numIntStates = 8 * self.numRealsInScalar
	self.numWaves = 8 * self.numRealsInScalar

	self.consVars = {
		{name='D', type=self.vec3, units='C/m^2', variance='l'},
		{name='B', type=self.vec3, units='kg/(C*s)', variance='l'},
		{name='phi', type=self.scalar, units='C/m^2'},		-- div D potential
		{name='psi', type=self.scalar, units='kg/(C*s)'},	-- div B potential
		{name='rhoCharge', type=self.scalar, units='C/m^3'},
		{name='sigma', type=self.scalar, units='(C^2*s)/(kg*m^3)'},
		{name='sqrt_1_eps', type=self.susc_t, units='(kg*m^3)^.5/(C*s)', variance=''},
		{name='sqrt_1_mu', type=self.susc_t, units='C/(kg*m)^.5', variance=''},
	}

	self.eigenVars = table{
		{name='sqrt_1_eps', type=self.susc_t, units='(kg*m^3)^.5/(C*s)', variance=''},
		{name='sqrt_1_mu', type=self.susc_t, units='C/(kg*m)^.5', variance=''},
	}

	GLM_Maxwell.super.init(self, args)

	self.postComputeFluxCode = self:template[[
//// MODULE_DEPENDS: coord_sqrt_det_g eqn.common
<? local vec3 = eqn.vec3 ?>
		//TODO shouldn't I be transforming both the left and right fluxes by the metrics at their respective coordinates?
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
]]
end

function GLM_Maxwell:createInitState()
	GLM_Maxwell.super.createInitState(self)

	--[[
	-- how to set speedOfLight == 1?
	-- speedOfLight / (unit_m / unit_s) = 1
	-- unit_s = unit_m / speedOfLight 
	local speedOfLight = require 'hydro.constants'.speedOfLight_in_m_per_s
	self.guiVars.second.value = self.guiVars.meter.value / speedOfLight
	...but this messes with the other variables that use seconds: B, psi, sigma, eps, mu
	so this asks the question ... when should values be converted to/from units?
	right now I'm converting only solver-> guivars
	...but I'm not converting solver-> guivars within init states ... hmm ...
	should I do that too?
	--]]
	-- [[
	local speedOfLight = 1
	--]]

	self:addGuiVars{
		{name='divPhiWavespeed', value=speedOfLight, units='m/s'},
		{name='divPsiWavespeed', value=speedOfLight, units='m/s'},
	}
end

-- don't use default
function GLM_Maxwell:initCodeModule_fluxFromCons() end

GLM_Maxwell.solverCodeFile = 'hydro/eqn/glm-maxwell.cl'

function GLM_Maxwell:getEnv()
	local scalar = self.scalar
	local env = GLM_Maxwell.super.getEnv(self)
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

GLM_Maxwell.predefinedDisplayVars = {
	'U D',
	'U div D',
	'U phi',
	'U B',
	'U div B',
	'U psi',
}

function GLM_Maxwell:getDisplayVars()
	local env = self:getEnv()
	
	local vars = GLM_Maxwell.super.getDisplayVars(self)
	vars:append{ 
		{name = 'E', code = self:template[[	value.v<?=vec3?> = calc_E(U);]], type=env.vec3, units='(kg*m)/(C*s)'},
		{name = 'H', code = self:template[[	value.v<?=vec3?> = calc_H(U);]], type=env.vec3, units='C/(m*s)'},
		{name = 'S', code = self:template[[	value.v<?=vec3?> = <?=vec3?>_cross(calc_E(U), calc_H(U));]], type=env.vec3, units='kg/s^3'},
		{
			name = 'energy', code = self:template[[
	<?=susc_t?> _1_eps = <?=susc_t?>_mul(U->sqrt_1_eps, U->sqrt_1_eps);
	<?=susc_t?> _1_mu = <?=susc_t?>_mul(U->sqrt_1_mu, U->sqrt_1_mu);
	value.vreal = <?=real_mul?>(<?=add?>(
		<?=scalar?>_<?=susc_t?>_mul(eqn_coordLenSq(U->D, x), _1_eps),
		<?=scalar?>_<?=susc_t?>_mul(eqn_coordLenSq(calc_H(U), x), _1_mu)
	), .5);
]],
			type = scalar,
			units = 'kg/(m*s^2)',
		},
	}
	
	vars:insert(self:createDivDisplayVar{field='D', scalar=env.scalar, units='C/m^3'} or nil)
	vars:insert(self:createCurlDisplayVar{field='D', scalar=env.scalar, units='C/m^3'} or nil)
	vars:insert(self:createDivDisplayVar{field='B', scalar=env.scalar, units='kg/(C*m*s)'} or nil)
	vars:insert(self:createCurlDisplayVar{field='B', scalar=env.scalar, units='kg/(C*m*s)'} or nil)

	return vars
end

function GLM_Maxwell:eigenWaveCodePrefix(n, eig, x, waveIndex)
--[=[
	return self:template([[
	<?=scalar?> v_p_abs = <?=mul?>((<?=eig?>)->sqrt_1_eps, (<?=eig?>)->sqrt_1_mu);
]], {
		eig = '('..eig..')',
	})
--]=]
-- [=[
	local env = self:getEnv()
	local code = self:template(
		[[<?=mul?>((<?=eig?>)->sqrt_1_eps, (<?=eig?>)->sqrt_1_mu)]],
		{
			eig = '('..eig..')',
		}
	)
	if self.scalar == 'cplx' then
		code = env.abs..'('..code..')'
	end
	return 'real v_p_abs = '..code..';'
--]=]
end

-- to use this, I really need cplx multiplications everywhere it is used
-- which is in the roe solver and hll solver
function GLM_Maxwell:eigenWaveCode(n, eig, x, waveIndex)
	waveIndex = math.floor(waveIndex / self.numRealsInScalar)
	return ({
		'-solver->divPhiWavespeed / unit_m_per_s',
		'-solver->divPsiWavespeed / unit_m_per_s',
		'-v_p_abs',
		'-v_p_abs',
		'v_p_abs',
		'v_p_abs',
		'solver->divPhiWavespeed / unit_m_per_s',
		'solver->divPsiWavespeed / unit_m_per_s',
	})[waveIndex+1] or error('got a bad waveIndex: '..waveIndex)
end

function GLM_Maxwell:eigenMaxWaveCode(n, eig, x)
	return 'max(max(solver->divPsiWavespeed / unit_m_per_s, solver->divPhiWavespeed / unit_m_per_s), v_p_abs)'
end
function GLM_Maxwell:eigenMinWaveCode(n, eig, x)
	return '-'..self:eigenMaxWaveCode(n, eig, x)
end


function GLM_Maxwell:consWaveCodePrefix(n, U, x, waveIndex) 
	local env = self:getEnv()
	local code = self:template(
		[[<?=mul?>((<?=U?>)->sqrt_1_eps, (<?=U?>)->sqrt_1_mu)]],
		{
			U = '('..U..')',
		}
	)
	if self.scalar == 'cplx' then
		code = env.abs..'('..code..')'
	end
	return 'real v_p_abs = '..code..';'
end
GLM_Maxwell.consWaveCode = GLM_Maxwell.eigenWaveCode

function GLM_Maxwell:consMaxWaveCode(n, U, x)
	return 'max(max(solver->divPsiWavespeed, solver->divPhiWavespeed), v_p_abs)'
end
function GLM_Maxwell:consMinWaveCode(n, U, x)
	return '-'..self:consMaxWaveCode(n, U, x)
end


return GLM_Maxwell
