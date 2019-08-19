local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local CLBuffer = require 'cl.obj.buffer'

local BSSNOKFiniteDifferenceSolver = require 'solver.bssnok-fd'

local BSSNOKFiniteDifferencePIRKSolver = class(BSSNOKFiniteDifferenceSolver)
BSSNOKFiniteDifferencePIRKSolver.name = 'BSSNOK_FiniteDifference_PIRK'

--[[
function BSSNOKFiniteDifferencePIRKSolver:calcDT()
	-- TODO based on grid
end
--]]

function BSSNOKFiniteDifferencePIRKSolver:createEqn()
	self.eqnArgs = table(self.eqnArgs)
	
	-- change default cfl method
	self.eqnArgs.cflMethod = self.eqnArgs.cflMethod or '2013 Baumgarte et al, eqn 32'
	
	BSSNOKFiniteDifferencePIRKSolver.super.createEqn(self)
end

function BSSNOKFiniteDifferencePIRKSolver:createBuffers()
	BSSNOKFiniteDifferencePIRKSolver.super.createBuffers(self)

	-- UBuf = U^n
	self:clalloc('U1', self.eqn.cons_t, self.numCells)
	self:clalloc('UNext', self.eqn.cons_t, self.numCells)
	self:clalloc('UTemp', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL1_1', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL1_n', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL2_1', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL2_n', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL2_next', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL3_1', self.eqn.cons_t, self.numCells)
	self:clalloc('derivL3_n', self.eqn.cons_t, self.numCells)
end

function BSSNOKFiniteDifferencePIRKSolver:refreshIntegrator()
	-- no integrator, because it's built into the solver
	-- TODO move the buffers kernels etc into this object
	self.integrator = {
		integrate = function(dt) end,
	}
end

function BSSNOKFiniteDifferencePIRKSolver:refreshSolverProgram()
	BSSNOKFiniteDifferencePIRKSolver.super.refreshSolverProgram(self)
	
	self.calcDeriv_PIRK_L1KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L1'
	self.calcDeriv_PIRK_L2_part1KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L2_part1'
	self.calcDeriv_PIRK_L3_part1KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L3_part1'
	self.calcDeriv_PIRK_L2_part2KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L2_part2'
	self.calcDeriv_PIRK_L3_part2KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L3_part2'
	self.calcDeriv_PIRK_L2_part3KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L2_part3'
	self.calcDeriv_PIRK_L3_part3KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L3_part3'
	self.copyWAlphaBetaKernelObj = self.solverProgramObj:kernel'copyWAlphaBeta'
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

-- perform the 4 steps for PIRK integration
function BSSNOKFiniteDifferencePIRKSolver:step(dt)
	local bufferSize = self.numCells * ffi.sizeof(self.eqn.cons_t)
	
	BSSNOKFiniteDifferencePIRKSolver.super.step(self, dt)
	
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL1_1, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL1_n, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL2_1, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL2_n, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL2_next, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL3_1, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL3_n, size=bufferSize}

	-- derivL1_n = PIRK_L1 (which is epsilon_IJ, W, alpha, beta^I) based on UBuf 
	self.calcDeriv_PIRK_L1KernelObj(self.solverBuf, self.derivL1_n, self.UBuf)

	-- U1 = UBuf + dt * derivL1_n fields epsilon_IJ, W, alpha, beta^I
	self.multAddKernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL1_n, real(dt))

	-- TODO then apply boundary to U1 fields epsilon_IJ, W, alpha, beta^I
	--self:boundary()

	-- UTemp = UBuf except W, alpha, beta^I come from U1 
	self.copyWAlphaBetaKernelObj(self.solverBuf, self.UTemp, self.UBuf, self.U1)
	
	-- derivL2_n = PIRK L2 part2 (which is LambdaBar^I) based on UBuf 
	-- derivL2_1 = PIRK L2 part2 (which is LambdaBar^I) based on UTemp
	-- derivL3_n = PIRK L3 part2 (which is LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_1, self.UTemp)
	self.calcDeriv_PIRK_L3_part2KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)

	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.app.cmds:enqueueCopyBuffer{dst=self.U1, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_1, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL3_n, real(dt))
	
	-- TODO apply boundary to U1 updated fields (LambdaBar^I)

	-- derivL2_n = PIRK L2 part1 (which is ABar_IJ, K) based on UBuf
	-- derivL2_1 = PIRK L2 part1 (which is ABar_IJ, K) based on U1
	-- derivL3_n = PIRK L3 part1 (which is ABar_IJ, K) based on UBuf
	self.calcDeriv_PIRK_L2_part1KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part1KernelObj(self.solverBuf, self.derivL2_1, self.U1)
	self.calcDeriv_PIRK_L3_part1KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)

	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.app.cmds:enqueueCopyBuffer{dst=self.U1, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_1, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL3_n, real(dt))
	
	-- TODO apply boundary to U1 updated fields (ABar_IJ, K)

	-- derivL2_n = PIRK L2 part2 (which is LambdaBar^I) based on UBuf 
	-- derivL2_1 = PIRK L2 part2 (which is LambdaBar^I) based on U1
	-- derivL3_n = PIRK L3 part2 (which is LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_1, self.U1)
	self.calcDeriv_PIRK_L3_part2KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.app.cmds:enqueueCopyBuffer{dst=self.U1, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_1, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL3_n, real(dt))
	
	-- TODO apply boundary to U1 updated fields (LambdaBar^I)

	-- derivL2_n = PIRK L2 part3 (which is B^I) based on UBuf 
	-- derivL2_1 = PIRK L2 part3 (which is B^I) based on U1
	-- derivL3_n = PIRK L3 part3 (which is B^I) based on UBuf
	self.calcDeriv_PIRK_L2_part3KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part3KernelObj(self.solverBuf, self.derivL2_1, self.U1)
	self.calcDeriv_PIRK_L3_part3KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.app.cmds:enqueueCopyBuffer{dst=self.U1, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL2_1, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.U1, self.U1, self.derivL3_n, real(dt))
	
	-- TODO apply boundary to U1 updated fields (B^I)

	-- derivL1_1 = PIRK L1 (which is epsilon_IJ, W, alpha, beta^I) of U1
	self.calcDeriv_PIRK_L1KernelObj(self.solverBuf, self.derivL1_1, self.U1)
	
	-- UNext = .5 * (UBuf + U1 + dt * derivL1_1)
	self.app.cmds:enqueueFillBuffer{buffer=self.UNext, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.UBuf, real(.5))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.U1, real(.5))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL1_1, real(.5 * dt))
	
	-- TODO apply boundary to UNext updated fields (epsilon_IJ, W, alpha, beta^I)

	-- UTemp = U1 except with UNext from W, alpha, beta^I
	self.copyWAlphaBetaKernelObj(self.UTemp, self.U1, self.UNext)

	-- derivL2_n = PIRK L2 part2 (which is LambdaBar^I) based on UBuf
	-- derivL2_next = PIRK L2 part2 (which is LambdaBar^I) based on UTemp
	-- derivL3_n = PIRK L3 part2 (which is LambdaBar^I) based on UBuf
	-- derivL3_1 = PIRK L3 part2 (which is LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_next, self.UTemp)
	self.calcDeriv_PIRK_L3_part2KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_part2KernelObj(self.solverBuf, self.derivL3_1, self.UTemp)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.app.cmds:enqueueCopyBuffer{dst=self.UNext, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_next, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_1, real(.5 * dt))
	
	-- TODO apply boundary to UNext
	
	-- derivL2_n = PIRK L2 part1 of UBuf (A_IJ, K) of UBuf
	-- derivL2_next = PIRK L2 part1 of UBuf (A_IJ, K) of UNext
	-- derivL3_n = PIRK L3 part1 of UBuf (A_IJ, K) of UBuf
	-- derivL3_1 = PIRK L3 part1 of UBuf (A_IJ, K) of U1
	self.calcDeriv_PIRK_L2_part1KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part1KernelObj(self.solverBuf, self.derivL2_next, self.UNext)
	self.calcDeriv_PIRK_L3_part1KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_part1KernelObj(self.solverBuf, self.derivL3_1, self.U1)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.app.cmds:enqueueCopyBuffer{dst=self.UNext, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_next, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_1, real(.5 * dt))
	
	-- TODO apply boundary to UNext
	
	-- derivL2_n = PIRK L2 part2 (which is LambdaBar^I) based on UBuf
	-- derivL2_next = PIRK L2 part2 (which is LambdaBar^I) based on UNext
	-- derivL3_n = PIRK L3 part2 (which is LambdaBar^I) based on UBuf
	-- derivL3_1 = PIRK L3 part2 (which is LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part2KernelObj(self.solverBuf, self.derivL2_next, self.UNext)
	self.calcDeriv_PIRK_L3_part2KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_part2KernelObj(self.solverBuf, self.derivL3_1, self.UTemp)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.app.cmds:enqueueCopyBuffer{dst=self.UNext, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_next, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_1, real(.5 * dt))
	
	-- TODO apply boundary to UNext
	
	-- derivL2_n = PIRK L2 part3 (which is B^I) based on UBuf
	-- derivL2_next = PIRK L2 part3 (which is B^I) based on UNext
	-- derivL3_n = PIRK L3 part3 (which is B^I) based on UBuf
	-- derivL3_1 = PIRK L3 part3 (which is B^I) based on U1
	self.calcDeriv_PIRK_L2_part3KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_part3KernelObj(self.solverBuf, self.derivL2_next, self.UNext)
	self.calcDeriv_PIRK_L3_part3KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_part3KernelObj(self.solverBuf, self.derivL3_1, self.U1)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.app.cmds:enqueueCopyBuffer{dst=self.UNext, src=self.UBuf, size=bufferSize}
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL2_next, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_n, real(.5 * dt))
	self.multAddKernelObj(self.solverBuf, self.UNext, self.UNext, self.derivL3_1, real(.5 * dt))
	
	-- TODO apply boundary to UNext
	
	-- constrain metric det of UNext
	-- remove trace of UNext
	-- copy to UBuf

	self.app.cmds:enqueueCopyBuffer{dst=self.UBuf, src=self.UNext, size=bufferSize}
	self.constrainUKernelObj(self.solverBuf, self.UBuf)
	self:boundary()
end

return BSSNOKFiniteDifferencePIRKSolver
