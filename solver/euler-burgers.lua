local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local FiniteVolumeSolver = require 'solver.fvsolver'

local xNames = table{'x', 'y', 'z'}

-- TODO make this work with ops, specifically Euler's SelfGrav
local EulerBurgers = class(FiniteVolumeSolver)
EulerBurgers.name = 'EulerBurgers'

function EulerBurgers:init(...)
	EulerBurgers.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function EulerBurgers:createEqn()
	self.eqn = require 'eqn.euler'(self)
end

-- Usually the eqn provide their own 'calcDT'
-- but sometimes the solver needs it too ...
-- Maybe the eqn should flag whether it has 'calcDT'
-- and the solver should use it likewise?
-- This is also in solver/bssnok-fd.lua which also 
function EulerBurgers:getCalcDTCode() return '' end

function EulerBurgers:createBuffers()
	EulerBurgers.super.createBuffers(self)
	
	self:clalloc('intVelBuf', self.volume * self.dim * ffi.sizeof(self.app.real))
	self:clalloc('PBuf', self.volume * ffi.sizeof(self.app.real))
end

function EulerBurgers:getSolverCode()
	return table{
		EulerBurgers.super.getSolverCode(self),
		template(file['solver/euler-burgers.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function EulerBurgers:refreshSolverProgram()
	EulerBurgers.super.refreshSolverProgram(self)

	-- no mention of ULR just yet ...

	self.calcIntVelKernelObj = self.solverProgramObj:kernel('calcIntVel', self.intVelBuf, self.UBuf)
	self.calcFluxKernelObj = self.solverProgramObj:kernel('calcFlux', self.fluxBuf, self.UBuf, self.intVelBuf)

	self.computePressureKernelObj = self.solverProgramObj:kernel('computePressure', self.PBuf, self.UBuf) 
	
	self.diffuseMomentumKernelObj = self.solverProgramObj:kernel{name='diffuseMomentum', domain=self.domainWithoutBorder}
	self.diffuseMomentumKernelObj.obj:setArg(1, self.PBuf)
	
	self.diffuseWorkKernelObj = self.solverProgramObj:kernel'diffuseWork'
	self.diffuseWorkKernelObj.obj:setArg(1, self.UBuf)
	self.diffuseWorkKernelObj.obj:setArg(2, self.PBuf)
end

function EulerBurgers:addDisplayVars()
	EulerBurgers.super.addDisplayVars(self)

	self:addDisplayVarGroup{
		name = 'P', 
		vars = {
			{P = '*value = buf[index];'},
		},
	}

	for j,xj in ipairs(xNames) do
		self:addDisplayVarGroup{
			name = 'intVel', 
			codePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
]],
			vars = range(0,self.dim-1):map(function(i)
				return {[xj..'_'..i] = '*value = buf['..i..' + indexInt];'}
			end),
		}
	end
end

function EulerBurgers:step(dt)
	-- calc deriv here
	self.integrator:integrate(dt, function(derivBuf)
		self.calcIntVelKernelObj()

		self.calcFluxKernelObj.obj:setArg(3, ffi.new('real[1]', dt))
		self.calcFluxKernelObj()
	
		self.calcDerivFromFluxKernelObj(derivBuf)
	end)

	self:boundary()

	-- TODO potential update here
	
	self.integrator:integrate(dt, function(derivBuf)
		self.computePressureKernelObj()
	
		self.diffuseMomentumKernelObj(derivBuf)
	end)
	
	self:boundary()
	
	self.integrator:integrate(dt, function(derivBuf)
		self.diffuseWorkKernelObj(derivBuf)
	end)

	-- no addSource call just yet
	-- it is invoked by individual solvers
	-- however for eqn/euler there is nothing there except my messing with connection coefficients
end

return EulerBurgers
