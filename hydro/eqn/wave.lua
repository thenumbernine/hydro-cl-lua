-- check out my 'wave equation hyperbolic form' worksheet
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local Wave = class(Equation)
Wave.name = 'wave'

Wave.roeUseFluxFromCons = true

Wave.initConds = require 'hydro.init.euler':getList()	 -- use rho as our initial condition


--[[
args:
	scalar = 'real' or 'cplx'

	usePressure = whether to use pressure or density from the initial conditions for the wave equation
					Sedov initial conditions uses pressure
					gaussian, etc use density
--]]
function Wave:init(args)
	if args and args.usePressure ~= nil then
		self.usePressure = not not args.usePressure
	else
		self.usePressure = false
	end
	
	self.scalar = (args and args.scalar) or 'real'


	self.vec3 = self.scalar..'3'
	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'
	
	self.numStates = 4 * self.numRealsInScalar

	if self.usePressure then
		self.consVars = table{
			{name='Pi', type=self.scalar, units='kg/(m*s^3)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^2*s^2)'},
		}
	else
		self.consVars = table{
			{name='Pi', type=self.scalar, units='kg/(m^3*s)'},
			{name='Psi_l', type=self.vec3, units='kg/(m^4)'},
		}
	end

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

	Wave.super.init(self, args)
end

function Wave:compile(expr)
	return self.solver.coord:compile(expr)
end

function Wave:createInitState()
	Wave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},
	}
end

function Wave:getEnv()
	return table(Wave.super.getEnv(self), {
		scalar = self.scalar,
		vec3 = self.vec3,
	})
end

function Wave:initCodeModules()
	local symmath = require 'symmath'
	local Tensor = symmath.Tensor
	local Constant = symmath.Constant
	local var = symmath.var
	local vars = symmath.vars
	local fromlua = require 'ext.fromlua'
	local coords = Tensor.coords()[1].variables
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
	
	Wave.super.initCodeModules(self)
end


Wave.solverCodeFile = 'hydro/eqn/wave.cl'

-- don't use default
function Wave:initCodeModule_fluxFromCons() end

function Wave:getModuleDepends_waveCode()
	return {
		'units',
		self.solver.coord.symbols.normal_t,
		self.symbols.eqn_common,	-- metric_alpha
	}
end

function Wave:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
real wavespeed = solver->wavespeed / unit_m_per_s;
real alpha_nLen = metric_alpha(<?=x?>) * normal_len(n);
real beta_n = normal_vecDotN1(<?=n?>, metric_beta_u(<?=x?>));
]], {
		n = n,
		x = x,
	})
end

function Wave:consWaveCodePrefix(n, U, x)
	return self:eigenWaveCodePrefix(n, nil, x)
end

function Wave:consWaveCode(n, U, x, waveIndex)
	waveIndex = math.floor(waveIndex / self.numRealsInScalar)
	if waveIndex == 0 then
		return 'wavespeed * (-beta_n - alpha_nLen)' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return 'wavespeed * -beta_n'
	elseif waveIndex == 3 then
		return 'wavespeed * (-beta_n + alpha_nLen)' 
	end
	error'got a bad waveIndex'
end

Wave.eigenWaveCode = Wave.consWaveCode

return Wave
