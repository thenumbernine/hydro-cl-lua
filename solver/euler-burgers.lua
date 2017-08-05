local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local Solver = require 'solver.solver'

local EulerBurgers = class(Solver)
EulerBurgers.name = 'EulerBurgers'

function EulerBurgers:init(...)
	EulerBurgers.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

local BurgersEulerEqn = class(require 'eqn.euler')
BurgersEulerEqn.hasCalcDT = true 
function EulerBurgers:createEqn(eqn)
	self.eqn = BurgersEulerEqn(self)
end

function EulerBurgers:createBuffers()
	EulerBurgers.super.createBuffers(self)
	
	-- TODO move this to fvsolver.lua
	self:clalloc('fluxBuf', self.volume * self.dim * ffi.sizeof(self.eqn.cons_t))
	
	self:clalloc('intVelBuf', self.volume * self.dim * ffi.sizeof(self.app.real))
	self:clalloc('PBuf', self.volume * ffi.sizeof(self.app.real))
end

function EulerBurgers:getSolverCode()
	return table{
		EulerBurgers.super.getSolverCode(self),
		template(file['solver/euler-burgers.cl'], {solver=self, eqn=self.eqn}),
		template(file['solver/calcDeriv.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function EulerBurgers:refreshSolverProgram()
	EulerBurgers.super.refreshSolverProgram(self)

	-- no mention of ULR just yet ...
	
	self.calcIntVelKernel = self.solverProgram:kernel('calcIntVel', self.intVelBuf, self.UBuf)
	self.calcFluxKernel = self.solverProgram:kernel('calcFlux', self.fluxBuf, self.UBuf, self.intVelBuf)

	-- TODO move this to fvsolver.lua
	self.calcDerivFromFluxKernel = self.solverProgram:kernel'calcDerivFromFlux'
	self.calcDerivFromFluxKernel:setArg(1, self.fluxBuf)

	self.computePressureKernel = self.solverProgram:kernel('computePressure', self.PBuf, self.UBuf) 
	
	self.diffuseMomentumKernel = self.solverProgram:kernel'diffuseMomentum'
	self.diffuseMomentumKernel:setArg(1, self.PBuf)
	
	self.diffuseWorkKernel = self.solverProgram:kernel'diffuseWork'
	self.diffuseWorkKernel:setArg(1, self.UBuf)
	self.diffuseWorkKernel:setArg(2, self.PBuf)
end

function EulerBurgers:addConvertToTexs()
	EulerBurgers.super.addConvertToTexs(self)

	-- TODO move to solverfv
	self:addConvertToTex{
		name = 'flux', 
		type = self.eqn.cons_t,
		varCodePrefix = [[
	const global ]]..self.eqn.cons_t..[[* flux = buf + indexInt;
]],
		vars = range(0,self.eqn.numStates-1):map(function(i)
			return {[tostring(i)] = '*value = flux->ptr['..i..'];'}
		end),
	}

	self:addConvertToTex{
		name = 'P', 
		vars = {
			{P = '*value = buf[index];'},
		},
	}

	self:addConvertToTex{
		name = 'intVel', 
		vars = range(0,self.dim-1):map(function(i)
			return {[tostring(i)] = '*value = buf['..i..' + indexInt];'}
		end),
	}
end

function EulerBurgers:step(dt)
	-- calc deriv here
	self.integrator:integrate(dt, function(derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcIntVelKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}

		self.calcFluxKernel:setArg(3, ffi.new('real[1]', dt))
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	
		self.calcDerivFromFluxKernel:setArg(0, derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
	end)

	self:boundary()

	-- TODO potential update here
	
	self.integrator:integrate(dt, function(derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.computePressureKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	
		self.diffuseMomentumKernel:setArg(0, derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.diffuseMomentumKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
	end)
	
	self:boundary()
	
	self.integrator:integrate(dt, function(derivBuf)
		self.diffuseWorkKernel:setArg(0, derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.diffuseWorkKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
	end)
end

return EulerBurgers 
