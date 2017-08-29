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

	self.diffuseKernelObj = self.solverProgramObj:kernel{name='diffuse', setArgs={self.UNextBuf, self.UBuf}, domain=self.domainWithoutBorder}
	self.advectKernelObj = self.solverProgramObj:kernel{name='advect', setArgs={self.UNextBuf, self.UBuf}, domain=self.domainWithoutBorder}
	self.calcDivKernelObj = self.solverProgramObj:kernel{name='calcDiv', setArgs={self.divBuf, self.UBuf}, domain=self.domainWithoutBorder}
	self.diffusePressureKernelObj = self.solverProgramObj:kernel{name='diffusePressure', setArgs={self.PBuf, self.divBuf}, domain=self.domainWithoutBorder}
	self.projectKernelObj = self.solverProgramObj:kernel{name='project', setArgs={self.UBuf, self.PBuf}, domain=self.domainWithoutBorder}
end

function NavierStokesIncompressible:getCalcDTCode() end
function NavierStokesIncompressible:refreshCalcDTKernel() end
function NavierStokesIncompressible:calcDT() return self.fixedDT end

-- TODO options for other solvers?
NavierStokesIncompressible.numGaussSeidelSteps = 20

function NavierStokesIncompressible:project()
	self.calcDivKernelObj()

	for i=1,self.numGaussSeidelSteps do
		self.diffusePressureKernelObj()
	end
	
	self.projectKernelObj()
end

function NavierStokesIncompressible:step(dt)
	local bufferSize = solver.volume * ffi.sizeof(self.eqn.cons_t)

	self.diffuseKernelObj.obj:setArg(2, ffi.new('real[1]', dt))
	
	-- diffuse
	for i=1,self.numGaussSeidelSteps do
		self.diffuseKernelObj()
		solver.app.cmds:enqueueCopyBuffer{src=solver.UNextBuf, dst=self.UBuf, size=bufferSize}
	end
	
	self:project()
	
	self.advectKernel()
	
	self:project()
end

return NavierStokesIncompressible
