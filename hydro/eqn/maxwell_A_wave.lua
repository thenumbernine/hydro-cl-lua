-- me taking the wave equation (without background for now ...)
-- and refitting it as the EM-four-potential wave-equation representation of the Maxwell equations

local ffi = require 'ffi'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local MaxwellAWave = Equation:subclass()
MaxwellAWave.name = 'maxwell_A_wave'

MaxwellAWave.roeUseFluxFromCons = true

MaxwellAWave.initConds = require 'hydro.init.euler':getList()	 -- use rho as our initial condition


--[[
args:
	scalar = 'real' or 'cplx'
--]]
function MaxwellAWave:init(args)
	self.scalar = (args and args.scalar) or 'real'

	-- TODO when scalar==cplx then we need to cdef scalar for numRealsInScalar to work, so...
	if self.scalar == 'cplx' then
		require 'hydro.code.safecdef'(
			args.solver.app.modules:getTypeHeader'cplx'
		)
		-- annnd I still need to change the eigen left & right to work too...
	end

	self.vec3 = self.scalar..'3'
	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'
	
	self.numStates = 16 * self.numRealsInScalar

	local suffix = self.scalar == 'real' and '' or ' re'
	self.predefinedDisplayVars = {
		'U dtAt'..suffix,
		'U diAt_l'..suffix,
		'U diAt_l x'..suffix,
		'U diAt_l y'..suffix,
		'U diAt_l z'..suffix,
		'U diAt_l mag metric',
	}

	self.init_f = args.f

	MaxwellAWave.super.init(self, args)
end

function MaxwellAWave:getSymbolFields()
	return table(MaxwellAWave.super.getSymbolFields(self)):append{
		'metric_f',
	}
end

function MaxwellAWave:buildVars()
	MaxwellAWave.super.buildVars(self)
	
	self.consVars = self.consVars or table()

	-- A_i = units of V s / m = kg m / (C s)
	-- A_t = 1/c φ = kg m / (C s)
	-- so φ = kg m^2 / (C s^2)
	-- ∂_t A_t = kg m / (C s^2)
	-- ∂_i A_t = kg / (C s)
	self.consVars:append{
		{name='dtAt', type=self.scalar, units='(kg*m)/(C*s^2)'},	-- ∂_t A_t
		{name='diAt_l', type=self.vec3, units='kg/(C*s)'},			-- ∂_i A_t
		
		-- TODO per-dimension derivative? or just use one-higher rank?
		-- another TODO is I bet this looks a lot like hydro/eqn/lingr.lua ... though that has more stuff for the spacetime metric
		{name='dtAx', type=self.scalar, units='(kg*m)/(C*s^2)'},	-- ∂_t A_x
		{name='diAx_l', type=self.vec3, units='kg/(C*s)'},			-- ∂_i A_x
		{name='dtAy', type=self.scalar, units='(kg*m)/(C*s^2)'},	-- ∂_t A_y
		{name='diAy_l', type=self.vec3, units='kg/(C*s)'},			-- ∂_i A_y
		{name='dtAz', type=self.scalar, units='(kg*m)/(C*s^2)'},	-- ∂_t A_z
		{name='diAz_l', type=self.vec3, units='kg/(C*s)'},			-- ∂_i A_z
	}
end

function MaxwellAWave:compile(expr)
	return self.solver.coord:compile(expr)
end

function MaxwellAWave:createInitState()
	MaxwellAWave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},	-- TODO default to the speed-of-light.  also rewrite in terms of eps and mu.
	}
end

function MaxwellAWave:getEnv()
	return table(MaxwellAWave.super.getEnv(self), {
		scalar = self.scalar,
		vec3 = self.vec3,
	})
end

function MaxwellAWave:initCodeModules()
	local symmath = require 'symmath'
	local Tensor = symmath.Tensor
	local Constant = symmath.Constant
	local var = symmath.var
	local vars = symmath.vars
	local fromlua = require 'ext.fromlua'
	local coords = self.solver.coord.symchart.coords
	local t = var't'
	self.metric = {
		coords = coords,
		t = t,
		f = Constant(0),
	}

	local x = self.solver.coord.vars.x
	local y = self.solver.coord.vars.y
	local z = self.solver.coord.vars.z
	local r = self.solver.coord.vars.r
	local function readarg(s)
		return symmath.clone(assert(load('local x,y,z,t,r = ... return '..s))(x,y,z,t,r))()
	end

	-- this isn't really a metric variable
	if self.init_f then
		self.metric.f = readarg(self.init_f)
	end
	
	MaxwellAWave.super.initCodeModules(self)
end


MaxwellAWave.solverCodeFile = 'hydro/eqn/maxwell_A_wave.cl'

-- don't use default
function MaxwellAWave:initCodeModule_fluxFromCons() end

function MaxwellAWave:getDisplayVars()
	local env = self:getEnv()
	local vars = MaxwellAWave.super.getDisplayVars(self)
	vars:append{ 
		-- TODO enforce this as an operator
		{
			name = 'Lorentz gauge',
			code = [[value.vreal = -U->dtAt + U->diAx_l.x + U->diAy_l.y + U->diAz_l.z;]],
		},
		{
			name = 'E',
			type = env.vec3,
			units = '(kg*m)/(C*s^2)',
-- TODO units.  distinguish between partial_0 vs partial_t and A_0 vs A_t
			code = self:template[[
value.v<?=vec3?>.x = U->diAt_l.x - U->dtAx;
value.v<?=vec3?>.y = U->diAt_l.y - U->dtAy;
value.v<?=vec3?>.z = U->diAt_l.z - U->dtAz;
]],
		},
		{
			name = 'B',
			type = env.vec3,
			units = 'kg/(C*s)',
-- TODO sign. i just winged this. 
			code = self:template[[
value.v<?=vec3?>.x = U->diAy_l.z - U->diAy_l.y;
value.v<?=vec3?>.y = U->diAz_l.x - U->diAx_l.z;
value.v<?=vec3?>.z = U->diAx_l.y - U->diAy_l.x;
]],
		},
	}
	return vars
end

function MaxwellAWave:eigenWaveCodePrefix(args)
	return self:template([[
real const wavespeed = solver->wavespeed / unit_m_per_s;
real const nLen = normal_len(n);
]], args)
end

function MaxwellAWave:eigenWaveCode(args)
	local waveIndex = math.floor(args.waveIndex / self.numRealsInScalar) % 4
	if waveIndex == 0 then
		return '-wavespeed * nLen' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return '0.'
	elseif waveIndex == 3 then
		return 'wavespeed * nLen' 
	end
	error'got a bad waveIndex'
end

-- safe so long as eigenWaveCode[Prefix] doesn't use args.eig
MaxwellAWave.consWaveCodePrefix = MaxwellAWave.eigenWaveCodePrefix
MaxwellAWave.consWaveCode = MaxwellAWave.eigenWaveCode

function MaxwellAWave:eigenWaveCodeMinMax(args)
	return self:eigenWaveCodePrefix(args)
	..self:template([[
<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'-wavespeed * nLen',
	'wavespeed * nLen'
)?>
]], args)
end

-- safe so long as eigenWaveCode[Prefix] doesn't use args.eig
MaxwellAWave.consWaveCodeMinMax = MaxwellAWave.eigenWaveCodeMinMax

function MaxwellAWave:consWaveCodeMinMaxAllSidesPrefix(args)
	return self:template([[
real const wavespeed = solver->wavespeed / unit_m_per_s;
]], args)
end

function MaxwellAWave:consWaveCodeMinMaxAllSides(args)
	return self:template([[
real const nLen = normal_len(n);
<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'-wavespeed * nLen',
	'wavespeed * nLen'
)?>
]], args)
end

return MaxwellAWave
