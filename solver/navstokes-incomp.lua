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

	self.diffuseKernelObj = self.solverProgramObj:kernel('diffuse', self.UNextBuf, self.UBuf)
	self.advectKernelObj = self.solverProgramObj:kernel('advect', self.UNextBuf, self.UBuf)
	self.calcDivKernelObj = self.solverProgramObj:kernel('calcDiv', self.divBuf, self.UBuf)
	self.diffusePressureKernelObj = self.solverProgramObj:kernel('diffusePressure', self.PBuf, self.divBuf)
	self.projectKernelObj = self.solverProgramObj:kernel('project', self.UBuf, self.PBuf)
end

function NavierStokesIncompressible:getCalcDTCode() end
function NavierStokesIncompressible:refreshCalcDTKernel() end
function NavierStokesIncompressible:calcDT() return self.fixedDT end

-- TODO options for other solvers?
NavierStokesIncompressible.numGaussSeidelSteps = 20

function NavierStokesIncompressible:project()
	self.calcDivKernelObj:callWithoutBorder()

	for i=1,self.numGaussSeidelSteps do
		self.diffusePressureKernelObj:callWithoutBorder()
	end
	
	self.projectKernelObj:callWithoutBorder()
end

function NavierStokesIncompressible:step(dt)
	local bufferSize = solver.volume * ffi.sizeof(self.eqn.cons_t)

	self.diffuseKernelObj.obj:setArg(2, ffi.new('real[1]', dt))
	
	-- diffuse
	for i=1,self.numGaussSeidelSteps do
		self.diffuseKernelObj:callWithoutBorder()
		solver.app.cmds:enqueueCopyBuffer{src=solver.UNextBuf, dst=self.UBuf, size=bufferSize}
	end
	
	self:project()
	
	self.advectKernel()
	
	self:project()
end

return NavierStokesIncompressible
