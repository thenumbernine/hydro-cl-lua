-- simple wave equation with no extra background metric
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

	local suffix = self.scalar == 'real' and '' or ' re'
	self.predefinedDisplayVars = {
		'U Pi'..suffix,
		'U Psi_l'..suffix,
		'U Psi_l x'..suffix,
		'U Psi_l y'..suffix,
		'U Psi_l z'..suffix,
		'U Psi_l mag metric',
	}

	self.init_f = args.f

	Wave.super.init(self, args)
end

function Wave:getSymbolFields()
	return table(Wave.super.getSymbolFields(self)):append{
		'metric_f',
	}
end

function Wave:buildVars()
	Wave.super.buildVars(self)
	
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

function Wave:compile(expr)
	return self.solver.coord:compile(expr)
end

function Wave:createInitState()
	Wave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},
	}
end

function Wave:initCodeModules()
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
	
	Wave.super.initCodeModules(self)
end


Wave.solverCodeFile = 'hydro/eqn/wave.clcpp'

-- don't use default
function Wave:initCodeModule_fluxFromCons() end

return Wave
