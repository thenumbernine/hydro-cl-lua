local class = require 'ext.class'
local real = require 'real'
local GridSolver = require 'solver.gridsolver'

local NavierStokesIncompressible = class(GridSolver)
NavierStokesIncompressible.name = 'NavierStokesIncompressible' 
NavierStokesIncompressible.eqnName = 'navstokes-incomp'

function NavierStokesIncompressible:createBuffers()
	NavierStokesIncompressible.super.createBuffers(self)

	self:clalloc('UNextBuf', self.eqn.cons_t, self.numCells)
	self:clalloc('divBuf', self.app.real, self.numCells)
	self:clalloc('PBuf', self.app.real, self.numCells)
end

function NavierStokesIncompressible:getSolverCode()
	return table{
		NavierStokesIncompressible.super.getSolverCode(self),
		template(file['solver/navstokes-incomp.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function NavierStokesIncompressible:refreshSolverProgram()
	NavierStokesIncompressible.super.refreshSolverProgram(self)

	self.diffuseKernelObj = self.solverProgramObj:kernel{name='diffuse', args={self.UNextBuf, self.UBuf}, domain=self.domainWithoutBorder}
	self.advectKernelObj = self.solverProgramObj:kernel{name='advect', args={self.UNextBuf, self.UBuf}, domain=self.domainWithoutBorder}
	self.calcDivKernelObj = self.solverProgramObj:kernel{name='calcDiv', args={self.divBuf, self.UBuf}, domain=self.domainWithoutBorder}
	self.diffusePressureKernelObj = self.solverProgramObj:kernel{name='diffusePressure', args={self.PBuf, self.divBuf}, domain=self.domainWithoutBorder}
	self.projectKernelObj = self.solverProgramObj:kernel{name='project', args={self.UBuf, self.PBuf}, domain=self.domainWithoutBorder}
end

function NavierStokesIncompressible:refreshCalcDTKernel() end
function NavierStokesIncompressible:calcDT() return self.fixedDT end

-- TODO options for other solvers?
NavierStokesIncompressible.numJacobiSteps = 20

function NavierStokesIncompressible:project()
	self.calcDivKernelObj()

	for i=1,self.numJacobiSteps do
		self.diffusePressureKernelObj()
	end
	
	self.projectKernelObj()
end

function NavierStokesIncompressible:step(dt)
	local bufferSize = self.numCells * ffi.sizeof(self.eqn.cons_t)

	self.diffuseKernelObj.obj:setArg(2, real(dt))
	
	-- diffuse
	for i=1,self.numJacobiSteps do
		self.diffuseKernelObj()
		self.cmds:enqueueCopyBuffer{src=self.UNextBuf, dst=self.UBuf, size=bufferSize}
	end
	
	self:project()
	
	self.advectKernel()
	
	self:project()
end

return NavierStokesIncompressible
