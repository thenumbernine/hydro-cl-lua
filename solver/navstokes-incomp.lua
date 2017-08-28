local class = require 'ext.class'
local Solver = require 'solver.solver'

local NavierStokesIncompressible = class(Solver)
NavierStokesIncompressible.name = 'NavierStokesIncompressible' 

function NavierStokesIncompressible:createEqn(eqn)
	self.eqn = require 'eqn.navstokes-incomp'(self)
end

function NavierStokesIncompressible:createBuffers()
	NavierStokesIncompressible.super.createBuffers(self)

	self:clalloc('UNextBuf', self.volume * ffi.sizeof(self.eqn.cons_t))
	self:clalloc('divBuf', self.volume * ffi.sizeof(self.app.real))
	self:clalloc('PBuf', self.volume * ffi.sizeof(self.app.real))
end

function NavierStokesIncompressible:getSolverCode()
	return table{
		NavierStokesIncompressible.super.getSolverCode(self),
		template(file['solver/navstokes-incomp.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function NavierStokesIncompressible:refreshSolverProgram()
	NavierStokesIncompressible.super.refreshSolverProgram(self)

	self.diffuseKernel = self.solverProgramObj.obj:kernel('diffuse', self.UNextBuf, self.UBuf)
	self.advectKernel = self.solverProgramObj.obj:kernel('advect', self.UNextBuf, self.UBuf)
	self.calcDivKernel = self.solverProgramObj.obj:kernel('calcDiv', self.divBuf, self.UBuf)
	self.diffusePressureKernel= self.solverProgramObj.obj:kernel('diffusePressure', self.PBuf, self.divBuf)
	self.projectKernel= self.solverProgramObj.obj:kernel('project', self.UBuf, self.PBuf)
end

function NavierStokesIncompressible:getCalcDTCode() end
function NavierStokesIncompressible:refreshCalcDTKernel() end
function NavierStokesIncompressible:calcDT() return self.fixedDT end

-- TODO options for other solvers?
NavierStokesIncompressible.numGaussSeidelSteps = 20

function NavierStokesIncompressible:project()
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDivKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}

	for i=1,self.numGaussSeidelSteps do
		self.app.cmds:enqueueNDRangeKernel{kernel=self.diffusePressureKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
	end
	
	self.app.cmds:enqueueNDRangeKernel{kernel=self.projectKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
end

function NavierStokesIncompressible:step(dt)
	local bufferSize = solver.volume * ffi.sizeof(self.eqn.cons_t)

	self.diffuseKernel:setArg(2, ffi.new('real[1]', dt))
	
	-- diffuse
	for i=1,self.numGaussSeidelSteps do
		self.app.cmds:enqueueNDRangeKernel{kernel=self.diffuseKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
		solver.app.cmds:enqueueCopyBuffer{src=solver.UNextBuf, dst=self.UBuf, size=bufferSize}
	end
	
	self:project()
	
	self.advectKernel()
	
	self:project()
end

return NavierStokesIncompressible
