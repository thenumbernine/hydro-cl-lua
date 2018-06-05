local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local FiniteVolumeSolver = require 'solver.fvsolver'

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


-- TODO make this work with ops, specifically Euler's SelfGrav
local EulerBurgers = class(FiniteVolumeSolver)
EulerBurgers.name = 'EulerBurgers'
EulerBurgers.eqnName = 'euler'

function EulerBurgers:init(...)
	EulerBurgers.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function EulerBurgers:createEqn()
	EulerBurgers.super.createEqn(self)
	self.eqn.getCalcDTCode = function() end	-- override calcDT 
end

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

	self.calcIntVelKernelObj = self.solverProgramObj:kernel('calcIntVel', self.intVelBuf, self.getULRBuf)
	self.calcFluxKernelObj = self.solverProgramObj:kernel('calcFlux', self.fluxBuf, self.getULRBuf, self.intVelBuf)

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

local realptr = ffi.new'real[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

function EulerBurgers:step(dt)
	-- calc deriv here
	self.integrator:integrate(dt, function(derivBuf)
		self.calcIntVelKernelObj()

		self.calcFluxKernelObj.obj:setArg(3, real(dt))
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
