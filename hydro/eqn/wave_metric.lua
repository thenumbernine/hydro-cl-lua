-- check out my 'wave equation hyperbolic form' worksheet
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local WaveMetric = class(Equation)
WaveMetric.name = 'wave_metric'

WaveMetric.roeUseFluxFromCons = true

WaveMetric.initConds = require 'hydro.init.euler':getList()	 -- use rho as our initial condition


--[[
args:
	scalar = 'real' or 'cplx'

	usePressure = whether to use pressure or density from the initial conditions for the wave equation
					Sedov initial conditions uses pressure
					gaussian, etc use density
--]]
function WaveMetric:init(args)
	if args and args.usePressure ~= nil then
		self.usePressure = not not args.usePressure
	else
		self.usePressure = false
	end
	
	self.scalar = (args and args.scalar) or 'real'


	self.vec3 = self.scalar..'3'
	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'
	
	self.numStates = 4 * self.numRealsInScalar

	local suffix = self.scalar == 'real' and '' or ' re'
	self.predefinedDisplayVars = {
		'U Pi'..suffix,
		'U Psi_l'..suffix,
		'U Psi_l x'..suffix,
		'U Psi_l y'..suffix,
		'U Psi_l z'..suffix,
		'U Psi_l mag metric',
	}

	self.init_alpha = args.alpha
	self.init_beta = args.beta
	self.init_K = args.K
	self.init_f = args.f

	WaveMetric.super.init(self, args)
end

function WaveMetric:buildVars()
	WaveMetric.super.buildVars(self)
	
	self.consVars = self.consVars or table()

	if self.usePressure then
		self.consVars:append{
			{name='Pi', type=self.scalar, units='kg/(m*s^3)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^2*s^2)'},
		}
	else
		self.consVars:append{
			{name='Pi', type=self.scalar, units='kg/(m^3*s)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^4)'},
		}
	end
end

function WaveMetric:compile(expr)
	return self.solver.coord:compile(expr)
end

function WaveMetric:createInitState()
	WaveMetric.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},
	}
end

function WaveMetric:getEnv()
	return table(WaveMetric.super.getEnv(self), {
		scalar = self.scalar,
		vec3 = self.vec3,
	})
end

function WaveMetric:initCodeModules()
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
		
		alpha = Constant(1),
		beta_u = {
			Constant(0),
			Constant(0),
			Constant(0),
		},
		K = Constant(0),
		f = Constant(0),
	}

	local x = self.solver.coord.vars.x
	local y = self.solver.coord.vars.y
	local z = self.solver.coord.vars.z
	local r = self.solver.coord.vars.r
	local function readarg(s)
		return symmath.clone(assert(load('local x,y,z,t,r = ... return '..s))(x,y,z,t,r))()
	end

	if self.init_alpha then
		self.metric.alpha = readarg(self.init_alpha)
	end
	if self.init_beta then
		for i=1,3 do
			if self.init_beta[i] then
				self.metric.beta_u[i] = readarg(self.init_beta[i])
			end
		end
	end
	if self.init_K then
		self.metric.K = readarg(self.init_K)
	end
	-- this isn't really a metric variable
	if self.init_f then
		self.metric.f = readarg(self.init_f)
	end
	
	WaveMetric.super.initCodeModules(self)
end


WaveMetric.solverCodeFile = 'hydro/eqn/wave_metric.cl'

-- don't use default
function WaveMetric:initCodeModule_fluxFromCons() end

function WaveMetric:eigenWaveCodePrefix(args)
	return self:template([[
real const wavespeed = solver->wavespeed / unit_m_per_s;
real const alpha_nLen = metric_alpha(<?=pt?>) * normal_len(n);
real const beta_n = normal_vecDotN1(<?=n?>, metric_beta_u(<?=pt?>));
]], args)
end

function WaveMetric:eigenWaveCode(args)
	local waveIndex = math.floor(args.waveIndex / self.numRealsInScalar)
	if waveIndex == 0 then
		return 'wavespeed * (-beta_n - alpha_nLen)' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return 'wavespeed * -beta_n'
	elseif waveIndex == 3 then
		return 'wavespeed * (-beta_n + alpha_nLen)' 
	end
	error'got a bad waveIndex'
end

-- safe so long as eigenWaveCode[Prefix] doesn't use args.eig
WaveMetric.consWaveCodePrefix = WaveMetric.eigenWaveCodePrefix
WaveMetric.consWaveCode = WaveMetric.eigenWaveCode

function WaveMetric:eigenWaveCodeMinMax(args)
	return self:eigenWaveCodePrefix(args)
	..self:template([[
<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'wavespeed * (-beta_n - alpha_nLen)',
	'wavespeed * (-beta_n + alpha_nLen)'
)?>
]], args)
end

-- safe so long as eigenWaveCode[Prefix] doesn't use args.eig
WaveMetric.consWaveCodeMinMax = WaveMetric.eigenWaveCodeMinMax

function WaveMetric:consWaveCodeMinMaxAllSidesPrefix(args)
	return self:template([[
real const wavespeed = solver->wavespeed / unit_m_per_s;
real const alpha = metric_alpha(<?=pt?>);
]], args)
end

function WaveMetric:consWaveCodeMinMaxAllSides(args)
	return self:template([[
real const alpha_nLen = alpha * normal_len(n);
real const beta_n = normal_vecDotN1(<?=n?>, metric_beta_u(<?=pt?>));
<?=eqn:waveCodeAssignMinMax(
	declare, resultMin, resultMax,
	'wavespeed * (-beta_n - alpha_nLen)',
	'wavespeed * (-beta_n + alpha_nLen)'
)?>
]], args)
end

return WaveMetric
