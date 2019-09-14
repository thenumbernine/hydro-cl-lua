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

function BSSNOKFiniteDifferencePIRKSolver:refreshBoundaryProgram()
	BSSNOKFiniteDifferencePIRKSolver.super.refreshBoundaryProgram(self)
	local args = self:getBoundaryProgramArgs()
	self.PIRK_EpsilonWAlphaBeta_BoundaryProgramObj, self.PIRK_EpsilonWAlphaBeta_BoundaryKernelObjs 
		= self:createBoundaryProgramAndKernel(table(args, {fields = {'epsilon_LL', 'W', 'alpha', 'beta_U'}}))
	self.PIRK_LambdaBar_BoundaryProgramObj, self.PIRK_LambdaBar_BoundaryKernelObjs 
		= self:createBoundaryProgramAndKernel(table(args, {fields = {'LambdaBar_U'}}))
	self.PIRK_ABarK_BoundaryProgramObj, self.PIRK_ABarK_BoundaryKernelObjs = 
		self:createBoundaryProgramAndKernel(table(args, {fields = {'ABar_LL', 'K'}}))
	self.PIRK_B_BoundaryProgramObj, self.PIRK_B_BoundaryKernelObjs =
		self:createBoundaryProgramAndKernel(table(args, {fields = {'B_U'}}))
end

function BSSNOKFiniteDifferencePIRKSolver:refreshSolverProgram()
	BSSNOKFiniteDifferencePIRKSolver.super.refreshSolverProgram(self)

	-- unlike the calcDeriv function used in int/*, these SET the deriv rather than adding to it
	self.calcDeriv_PIRK_L1_EpsilonWAlphaBeta_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L1_EpsilonWAlphaBeta'
	
	self.calcDeriv_PIRK_L2_ABarK_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L2_ABarK'
	self.calcDeriv_PIRK_L3_ABarK_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L3_ABarK'
	
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L2_LambdaBar'
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L3_LambdaBar'
	
	self.calcDeriv_PIRK_L2_B_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L2_B'
	self.calcDeriv_PIRK_L3_B_KernelObj = self.solverProgramObj:kernel'calcDeriv_PIRK_L3_B'


	self.PIRK_Eq1_EpsilonWAlphaBeta_KernelObj = self.solverProgramObj:kernel'PIRK_Eq1_EpsilonWAlphaBeta'
	
	self.PIRK_Eq2_ABarK_KernelObj = self.solverProgramObj:kernel'PIRK_Eq2_ABarK'
	self.PIRK_Eq2_LambdaBar_KernelObj = self.solverProgramObj:kernel'PIRK_Eq2_LambdaBar'
	self.PIRK_Eq2_B_KernelObj = self.solverProgramObj:kernel'PIRK_Eq2_B'
	
	self.PIRK_Eq3_EpsilonWAlphaBeta_KernelObj = self.solverProgramObj:kernel'PIRK_Eq3_EpsilonWAlphaBeta'
	
	self.PIRK_Eq4_ABarK_KernelObj = self.solverProgramObj:kernel'PIRK_Eq4_ABarK'
	self.PIRK_Eq4_LambdaBar_KernelObj = self.solverProgramObj:kernel'PIRK_Eq4_LambdaBar'
	self.PIRK_Eq4_B_KernelObj = self.solverProgramObj:kernel'PIRK_Eq4_B'


	self.copyWAlphaBetaKernelObj = self.solverProgramObj:kernel'copyWAlphaBeta'
	self.copyLambdaBarKernelObj = self.solverProgramObj:kernel'copyLambdaBar'
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

-- the name 'applyBoundaryToBuffer' was already taken
function BSSNOKFiniteDifferencePIRKSolver:applyBoundaryTo(buf, kernelObjs)
	assert(kernelObjs)
	for _,obj in ipairs(kernelObjs) do
		obj.obj:setArg(1, buf)
	end
	assert(kernelObjs)
	self:applyBoundaryToBuffer(kernelObjs)
	for _,obj in ipairs(kernelObjs) do
		obj.obj:setArg(1, self.UBuf)
	end
end

-- perform the 4 steps for PIRK integration
function BSSNOKFiniteDifferencePIRKSolver:step(dt)
	local bufferSize = self.numCells * ffi.sizeof(self.eqn.cons_t)
	
	BSSNOKFiniteDifferencePIRKSolver.super.step(self, dt)

	-- zero all derivatives
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL1_1, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL1_n, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL2_1, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL2_n, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL2_next, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL3_1, size=bufferSize}
	self.app.cmds:enqueueFillBuffer{buffer=self.derivL3_n, size=bufferSize}

-- step 1: evolve epsilon_IJ, alpha, W, beta^I

	-- derivL1_n = PIRK_L1 (epsilon_IJ, W, alpha, beta^I) based on UBuf of fields epsilon_IJ, W, alpha, beta^I
	self.calcDeriv_PIRK_L1_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.derivL1_n, self.UBuf)

	-- U1 = UBuf + dt * derivL1_n fields epsilon_IJ, W, alpha, beta^I
	self.PIRK_Eq1_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL1_n, real(dt))
	
	-- apply boundary to U1 fields epsilon_IJ, W, alpha, beta^I
	self:applyBoundaryTo(self.U1, self.PIRK_EpsilonWAlphaBeta_BoundaryKernelObjs)

-- step 2: evolve lambda^I

	-- UTemp = UBuf except W, alpha, beta^I come from U1 
	self.copyWAlphaBetaKernelObj(self.solverBuf, self.UTemp, self.UBuf, self.U1)

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf 
	-- derivL2_1 = PIRK L2 part2 (LambdaBar^I) based on UTemp
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_1, self.UTemp)
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)

	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_LambdaBar_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.U1, self.PIRK_LambdaBar_BoundaryKernelObjs)

-- step 3: ABar_IJ, K

	-- derivL2_n = PIRK L2 part1 (ABar_IJ, K) based on UBuf
	-- derivL2_1 = PIRK L2 part1 (ABar_IJ, K) based on U1
	-- derivL3_n = PIRK L3 part1 (ABar_IJ, K) based on UBuf
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_1, self.U1)
	self.calcDeriv_PIRK_L3_ABarK_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)

	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_ABarK_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (ABar_IJ, K)
	self:applyBoundaryTo(self.U1, self.PIRK_ABarK_BoundaryKernelObjs)

-- step 4: Lambda^I

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf 
	-- derivL2_1 = PIRK L2 part2 (LambdaBar^I) based on U1
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_1, self.U1)
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_LambdaBar_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.U1, self.PIRK_LambdaBar_BoundaryKernelObjs)

-- step 5: B^I

	-- derivL2_n = PIRK L2 part3 (B^I) based on UBuf 
	-- derivL2_1 = PIRK L2 part3 (B^I) based on U1
	-- derivL3_n = PIRK L3 part3 (B^I) based on UBuf
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_1, self.U1)
	self.calcDeriv_PIRK_L3_B_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_B_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (B^I)
	self:applyBoundaryTo(self.U1, self.PIRK_B_BoundaryKernelObjs)

-- step 6: epsilon_IJ, W, alpha, beta^I

	-- derivL1_1 = PIRK L1 (epsilon_IJ, W, alpha, beta^I) of U1
	self.calcDeriv_PIRK_L1_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.derivL1_1, self.U1)
	
	-- UNext = .5 * (UBuf + U1 + dt * derivL1_1)
	self.PIRK_Eq3_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.U1, self.derivL1_1, real(dt))
	
	-- apply boundary to UNext updated fields (epsilon_IJ, W, alpha, beta^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_EpsilonWAlphaBeta_BoundaryKernelObjs)

-- step 7: LambdaBar^I

	-- UTemp = U1 except with UNext from W, alpha, beta^I
	self.copyWAlphaBetaKernelObj(self.solverBuf, self.UTemp, self.U1, self.UNext)

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf
	-- derivL2_next = PIRK L2 part2 (LambdaBar^I) based on UTemp
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	-- derivL3_1 = PIRK L3 part2 (LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_next, self.UTemp)
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_1, self.UTemp)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_LambdaBar_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_LambdaBar_BoundaryKernelObjs)

-- step 8: ABar_IJ, K

	-- derivL2_n = PIRK L2 part1 of UBuf (ABar_IJ, K) of UBuf
	-- derivL2_next = PIRK L2 part1 of UBuf (ABar_IJ, K) of UNext
	-- derivL3_n = PIRK L3 part1 of UBuf (ABar_IJ, K) of UBuf
	-- derivL3_1 = PIRK L3 part1 of UBuf (ABar_IJ, K) of U1
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_next, self.UNext)
	self.calcDeriv_PIRK_L3_ABarK_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_ABarK_KernelObj(self.solverBuf, self.derivL3_1, self.U1)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_ABarK_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (ABar_IJ, K)
	self:applyBoundaryTo(self.UNext, self.PIRK_ABarK_BoundaryKernelObjs)

-- step 9: LambdaBar^I

	-- UTemp = UNext except LambdaBar^I come from U1
	self.copyLambdaBarKernelObj(self.solverBuf, self.UTemp, self.UNext, self.U1)

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf
	-- derivL2_next = PIRK L2 part2 (LambdaBar^I) based on UNext
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	-- derivL3_1 = PIRK L3 part2 (LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_next, self.UNext)
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_1, self.UTemp)	-- self.U1)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_LambdaBar_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_LambdaBar_BoundaryKernelObjs)

-- step 10: B^I

	-- derivL2_n = PIRK L2 part3 (B^I) based on UBuf
	-- derivL2_next = PIRK L2 part3 (B^I) based on UNext
	-- derivL3_n = PIRK L3 part3 (B^I) based on UBuf
	-- derivL3_1 = PIRK L3 part3 (B^I) based on U1
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf)
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_next, self.UNext)
	self.calcDeriv_PIRK_L3_B_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf)
	self.calcDeriv_PIRK_L3_B_KernelObj(self.solverBuf, self.derivL3_1, self.U1)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_B_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (B^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_B_BoundaryKernelObjs)

-- done

	-- constrain metric det of UNext
	-- remove trace of UNext
	-- copy to UBuf

	self.app.cmds:enqueueCopyBuffer{dst=self.UBuf, src=self.UNext, size=bufferSize}
	
-- this goes on in super, at the start of the function
--	self.constrainUKernelObj(self.solverBuf, self.UBuf)
--	self:boundary()
end

return BSSNOKFiniteDifferencePIRKSolver
