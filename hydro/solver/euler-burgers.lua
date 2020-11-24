local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local real = require 'hydro.real'
local FiniteVolumeSolver = require 'hydro.solver.fvsolver'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


-- TODO make this work with ops, specifically Euler's SelfGrav
local EulerBurgers = class(FiniteVolumeSolver)
EulerBurgers.name = 'EulerBurgers'

function EulerBurgers:init(...)
	EulerBurgers.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

-- override this -- don't expect a flux argument
-- ... or ... make -burgers a flux object?  but then you need to add the pressure integration as well
function EulerBurgers:createFlux(fluxName, fluxArgs)
	self.flux = {
		initCodeModules = function()
			-- this usually builds calcFluxForInterface, but euler-burgers doesn't use it anyways
		end,
	}
end

function EulerBurgers:initCodeModule_calcFlux()
	self.modules:add{name = 'calcFlux'}
end

-- TODO put this in its own eqn file? eqn/euler-burgers.lua?
local EulerEqn = require 'hydro.eqn.euler'
local EulerBurgersEqn = class(EulerEqn)
function EulerBurgersEqn:initCodeModule_calcDT() end

function EulerBurgers:createEqn()
	self.eqn = EulerBurgersEqn(table(self.eqnArgs, {solver=self}))
end

function EulerBurgers:createBuffers()
	EulerBurgers.super.createBuffers(self)
	
	self:clalloc('intVelBuf', self.app.real, self.numCells * self.dim)
	self:clalloc('PBuf', self.app.real, self.numCells)
end

function EulerBurgers:initCodeModules()
	EulerBurgers.super.initCodeModules(self)
	self.modules:addFromMarkup(template(file['hydro/solver/euler-burgers.cl'], {solver=self, eqn=self.eqn}))
	self.solverModulesEnabled['EulerBurgers.solver'] = true
end

function EulerBurgers:refreshSolverProgram()
	EulerBurgers.super.refreshSolverProgram(self)

	-- no mention of ULR just yet ...

	self.calcIntVelKernelObj = self.solverProgramObj:kernel'calcIntVel'
	self.calcFluxKernelObj = self.solverProgramObj:kernel'calcFlux'

	self.computePressureKernelObj = self.solverProgramObj:kernel('computePressure', self.solverBuf, self.PBuf, self.UBuf, self.cellBuf)
	
	self.diffuseMomentumKernelObj = self.solverProgramObj:kernel{name='diffuseMomentum', domain=self.domainWithoutBorder}
	self.diffuseMomentumKernelObj.obj:setArg(0, self.solverBuf)
	self.diffuseMomentumKernelObj.obj:setArg(2, self.PBuf)
	
	self.diffuseWorkKernelObj = self.solverProgramObj:kernel'diffuseWork'
	self.diffuseWorkKernelObj.obj:setArg(0, self.solverBuf)
	self.diffuseWorkKernelObj.obj:setArg(2, self.UBuf)
	self.diffuseWorkKernelObj.obj:setArg(3, self.PBuf)
end

function EulerBurgers:addDisplayVars()
	EulerBurgers.super.addDisplayVars(self)

	self:addDisplayVarGroup{
		name = 'P', 
		bufferType = 'real',
		bufferField = 'PBuf',
		vars = {
			{name='P', code='value.vreal = buf[index];'},
		},
	}

	for j,xj in ipairs(xNames) do
		self:addDisplayVarGroup{
			name = 'intVel', 
			bufferType = 'real',
			bufferField = 'intVelBuf',
			codePrefix = [[
	int const indexInt = ]]..(j-1)..[[ + dim * index;
]],
			vars = range(0,self.dim-1):map(function(i)
				return {name=xj..'_'..i, code='value.vreal = buf['..i..' + indexInt];'}
			end),
		}
	end
end

function EulerBurgers:step(dt)
	-- calc deriv here
	self.integrator:integrate(dt, function(derivBufObj)
		self.calcIntVelKernelObj(self.solverBuf, self.intVelBuf, self:getULRBuf(), self.cellBuf)

		self.calcFluxKernelObj(self.solverBuf, self.fluxBuf, self:getULRBuf(), self.intVelBuf, self.cellBuf, real(dt))
	
		self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
		self.calcDerivFromFluxKernelObj()
	
		if self.addSourceKernelObj then
			self.addSourceKernelObj(self.solverBuf, derivBufObj.obj, self:getULRBuf(), self.cellBuf)
		end
	end)

	self:boundary()

	-- TODO potential update here
	
	self.integrator:integrate(dt, function(derivBufObj)
		self.computePressureKernelObj()
	
		self.diffuseMomentumKernelObj.obj:setArg(1, derivBufObj.obj)
		self.diffuseMomentumKernelObj()
	end)
	
	self:boundary()
	
	self.integrator:integrate(dt, function(derivBufObj)
		self.diffuseWorkKernelObj.obj:setArg(1, derivBufObj.obj)
		self.diffuseWorkKernelObj()
	end)
end

return EulerBurgers
