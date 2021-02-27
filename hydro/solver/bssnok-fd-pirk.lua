local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local real = require 'hydro.real'
local CLBuffer = require 'cl.obj.buffer'
local BSSNOKFiniteDifferenceSolver = require 'hydro.solver.bssnok-fd'

local BSSNOKFiniteDifferencePIRKSolver = class(BSSNOKFiniteDifferenceSolver)
BSSNOKFiniteDifferencePIRKSolver.name = 'BSSNOK_FiniteDifference_PIRK'

function BSSNOKFiniteDifferencePIRKSolver:createEqn()
	self.eqnArgs = table(self.eqnArgs)
	
	-- change default cfl method
	self.eqnArgs.cflMethod = self.eqnArgs.cflMethod or '2013 Baumgarte et al, eqn 32'
	
	BSSNOKFiniteDifferencePIRKSolver.super.createEqn(self)
end

function BSSNOKFiniteDifferencePIRKSolver:createBuffers()
	BSSNOKFiniteDifferencePIRKSolver.super.createBuffers(self)

	-- UBuf = U^n
	self:clalloc('U1', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('UNext', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('UTemp', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL1_1', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL1_n', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL2_1', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL2_n', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL2_next', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL3_1', self.eqn.symbols.cons_t, self.numCells)
	self:clalloc('derivL3_n', self.eqn.symbols.cons_t, self.numCells)
end

function BSSNOKFiniteDifferencePIRKSolver:initCodeModules()
	BSSNOKFiniteDifferencePIRKSolver.super.initCodeModules(self)

	self.solverModulesEnabled[self.eqn.symbols.BSSNOK_PIRK] = true
end

-- TODO sort of, make this a real integrator,
-- but the compatible solver would need to provide, instead of calcDeriv, 
-- a matrix of callbacks for each var and for each PIRK stage
local FakeIntegrator = class()
function FakeIntegrator:integrate(dt) end

function BSSNOKFiniteDifferencePIRKSolver:refreshIntegrator()
	-- no integrator, because it's built into the solver
	-- TODO move the buffers kernels etc into this object
	self.integrator = FakeIntegrator()
end

function BSSNOKFiniteDifferencePIRKSolver:refreshBoundaryProgram()
	BSSNOKFiniteDifferencePIRKSolver.super.refreshBoundaryProgram(self)
	local args = self:getBoundaryProgramArgs()
	self.PIRK_EpsilonWAlphaBeta_BoundaryProgramObj, self.PIRK_EpsilonWAlphaBeta_BoundaryKernelObjs 
		= self:createBoundaryProgramAndKernel(table(args, {
			fields = {
				'alpha',
				'W',
				'epsilon_LL.xx',
				'epsilon_LL.xy',
				'epsilon_LL.xz',
				'epsilon_LL.yy',
				'epsilon_LL.yz',
				'epsilon_LL.zz',
				'beta_U.x',
				'beta_U.y',
				'beta_U.z',
			},
			programNameSuffix = '-epsilon_LL,W,alpha,beta_U',
		}))
	self.PIRK_LambdaBar_BoundaryProgramObj, self.PIRK_LambdaBar_BoundaryKernelObjs 
		= self:createBoundaryProgramAndKernel(table(args, {
			fields = {
				'LambdaBar_U.x',
				'LambdaBar_U.y',
				'LambdaBar_U.z',
			},
			programNameSuffx = '-LambdaBar_U',
		}))
	self.PIRK_ABarK_BoundaryProgramObj, self.PIRK_ABarK_BoundaryKernelObjs = 
		self:createBoundaryProgramAndKernel(table(args, {
			fields = {
				'K',
				'ABar_LL.xx',
				'ABar_LL.xy',
				'ABar_LL.xz',
				'ABar_LL.yy',
				'ABar_LL.yz',
				'ABar_LL.zz',
			},
			programNameSuffix = '-ABar_LL,K',
		}))
	self.PIRK_B_BoundaryProgramObj, self.PIRK_B_BoundaryKernelObjs =
		self:createBoundaryProgramAndKernel(table(args, {
			fields = {
				'B_U.x',
				'B_U.y',
				'B_U.z',
			},
			programNameSuffix = '-B_U',
		}))
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
	
	if require 'hydro.eqn.bssnok-fd-senr':isa(assert(self.eqn)) then
		-- to calc or to store detg
		-- TODO constrainU ?  
		-- but we don't need all of constrainU's functions. 
		-- we don't need to update RBar_IJ or H or M^I
		self.BSSN_Det_KernelObj = self.solverProgramObj:kernel'BSSN_Det'
	end
end

-- the name 'applyBoundaryToBuffer' was already taken
function BSSNOKFiniteDifferencePIRKSolver:applyBoundaryTo(buf, kernelObjs)
	assert(kernelObjs)
	for _,obj in ipairs(kernelObjs) do
		obj.obj:setArg(1, buf)
		obj.obj:setArg(2, self.cellBuf)
	end
	assert(kernelObjs)
	self:applyBoundaryToBuffer(kernelObjs)
	for _,obj in ipairs(kernelObjs) do
		obj.obj:setArg(1, self.UBuf)
		obj.obj:setArg(2, self.cellBuf)
	end
end

-- perform the 4 steps for PIRK integration
function BSSNOKFiniteDifferencePIRKSolver:step(dt)
	local bufferSize = self.numCells * ffi.sizeof(self.eqn.symbols.cons_t)

--[[ zero all derivatives
	self.cmds:enqueueFillBuffer{buffer=self.derivL1_1, size=bufferSize}
	self.cmds:enqueueFillBuffer{buffer=self.derivL1_n, size=bufferSize}
	self.cmds:enqueueFillBuffer{buffer=self.derivL2_1, size=bufferSize}
	self.cmds:enqueueFillBuffer{buffer=self.derivL2_n, size=bufferSize}
	self.cmds:enqueueFillBuffer{buffer=self.derivL2_next, size=bufferSize}
	self.cmds:enqueueFillBuffer{buffer=self.derivL3_1, size=bufferSize}
	self.cmds:enqueueFillBuffer{buffer=self.derivL3_n, size=bufferSize}
--]]
-- or don't.  we're going to overwrite them each frame, and never read from old values, so don't bother?

-- step 1: evolve epsilon_IJ, alpha, W, beta^I

	-- derivL1_n = PIRK_L1 (epsilon_IJ, W, alpha, beta^I) based on UBuf of fields epsilon_IJ, W, alpha, beta^I
	self.calcDeriv_PIRK_L1_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.derivL1_n, self.UBuf, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 1 - deriv L1:')
	self:printBuf(self.derivL1_nObj)
end

	-- U1 = UBuf + dt * derivL1_n fields epsilon_IJ, W, alpha, beta^I
	self.PIRK_Eq1_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL1_n, real(dt))
	
	-- apply boundary to U1 fields epsilon_IJ, W, alpha, beta^I
	self:applyBoundaryTo(self.U1, self.PIRK_EpsilonWAlphaBeta_BoundaryKernelObjs)

	-- re-constrain/recompute detg of U1 here, since its epsilon_IJ was updated
	if self.BSSN_Det_KernelObj then
		self.BSSN_Det_KernelObj(self.solverBuf, self.U1, self.cellBuf)
	end

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 1 - U1:')
	self:printBuf(self.U1Obj)
end

-- step 2: evolve lambda^I

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - UBuf:')
	self:printBuf(self.UBufObj)
end

	-- UTemp = UBuf except W, alpha, beta^I come from U1
	self.copyWAlphaBetaKernelObj(self.solverBuf, self.UTemp, self.UBuf, self.U1)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - UTemp:')
	self:printBuf(self.UTempObj)
end

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf 
	-- (reads _PHI and _DETG)
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)
	
if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - derivL2_n:')
	self:printBuf(self.derivL2_nObj)
end
if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - UTemp:')
	self:printBuf(self.UTempObj)
end

	-- derivL2_1 = PIRK L2 part2 (LambdaBar^I) based on UTemp
	-- (reads _PHI and _DETG)
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_1, self.UTemp, self.cellBuf)

-- mind you, UTemp looks like it has too many 0's
if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - derivL2_1:')
	self:printBuf(self.derivL2_1Obj)
end
if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - UBuf:')
	self:printBuf(self.UBufObj)
end

	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - derivL3_n:')
	self:printBuf(self.derivL3_nObj)
end
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_LambdaBar_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.U1, self.PIRK_LambdaBar_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 2 - U1:')
	self:printBuf(self.U1Obj)
end

-- step 3: ABar_IJ, K

	-- derivL2_n = PIRK L2 part1 (ABar_IJ, K) based on UBuf
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)
	
	-- derivL2_1 = PIRK L2 part1 (ABar_IJ, K) based on U1
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_1, self.U1, self.cellBuf)
	
	-- derivL3_n = PIRK L3 part1 (ABar_IJ, K) based on UBuf
	self.calcDeriv_PIRK_L3_ABarK_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)

	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_ABarK_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (ABar_IJ, K)
	self:applyBoundaryTo(self.U1, self.PIRK_ABarK_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 3 - U1:')
	self:printBuf(self.U1Obj)
end

-- step 4: Lambda^I

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf 
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)
	
	-- derivL2_1 = PIRK L2 part2 (LambdaBar^I) based on U1
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_1, self.U1, self.cellBuf)
	
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_LambdaBar_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.U1, self.PIRK_LambdaBar_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 4 - U1:')
	self:printBuf(self.U1Obj)
end

-- step 5: B^I

	-- derivL2_n = PIRK L2 part3 (B^I) based on UBuf 
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)
	
	-- derivL2_1 = PIRK L2 part3 (B^I) based on U1
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_1, self.U1, self.cellBuf)
	
	-- derivL3_n = PIRK L3 part3 (B^I) based on UBuf
	self.calcDeriv_PIRK_L3_B_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)
	
	-- U1 = UBuf + dt * (.5 * (derivL2_n + derivL2_1) + derivL3_n)
	self.PIRK_Eq2_B_KernelObj(self.solverBuf, self.U1, self.UBuf, self.derivL2_n, self.derivL2_1, self.derivL3_n, real(dt))
	
	-- apply boundary to U1 updated fields (B^I)
	self:applyBoundaryTo(self.U1, self.PIRK_B_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 5 - U1:')
	self:printBuf(self.U1Obj)
end

-- step 6: epsilon_IJ, W, alpha, beta^I

	-- derivL1_1 = PIRK L1 (epsilon_IJ, W, alpha, beta^I) of U1
	self.calcDeriv_PIRK_L1_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.derivL1_1, self.U1, self.cellBuf)
	
	-- UNext = .5 * (UBuf + U1 + dt * derivL1_1)
	self.PIRK_Eq3_EpsilonWAlphaBeta_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.U1, self.derivL1_1, real(dt))
	
	-- apply boundary to UNext updated fields (epsilon_IJ, W, alpha, beta^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_EpsilonWAlphaBeta_BoundaryKernelObjs)

	-- re-constrain/recompute detg of UNext here? only if you are storing it ...
	if self.BSSN_Det_KernelObj then
		self.BSSN_Det_KernelObj(self.solverBuf, self.UNext, self.cellBuf)
	end

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 6 - UNext:')
	self:printBuf(self.UNextObj)
end

-- step 7: LambdaBar^I

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - U1:')
	self:printBuf(self.U1Obj)
end

	-- UTemp = U1 except with UNext from W, alpha, beta^I
	self.copyWAlphaBetaKernelObj(self.solverBuf, self.UTemp, self.U1, self.UNext)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - UTemp:')
	self:printBuf(self.UTempObj)
end

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - derivL2_n:')
	self:printBuf(self.derivL2_nObj)
end

	-- derivL2_next = PIRK L2 part2 (LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_next, self.UTemp, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - derivL2_next:')
	self:printBuf(self.derivL2_nextObj)
end
	
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - derivL3_n:')
	self:printBuf(self.derivL3_nObj)
end

	-- derivL3_1 = PIRK L3 part2 (LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_1, self.UTemp, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - derivL3_1:')
	self:printBuf(self.derivL3_1Obj)
end

	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_LambdaBar_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))

	-- apply boundary to UNext updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_LambdaBar_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 7 - UNext:')
	self:printBuf(self.UNextObj)
end

-- step 8: ABar_IJ, K

	-- derivL2_n = PIRK L2 part1 of UBuf (ABar_IJ, K) of UBuf
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 8 - derivL2_n:')
	self:printBuf(self.derivL2_nObj)
end

	-- derivL2_next = PIRK L2 part1 of UBuf (ABar_IJ, K) of UNext
	self.calcDeriv_PIRK_L2_ABarK_KernelObj(self.solverBuf, self.derivL2_next, self.UNext, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 8 - derivL2_next:')
	self:printBuf(self.derivL2_nextObj)
end
	
	-- derivL3_n = PIRK L3 part1 of UBuf (ABar_IJ, K) of UBuf
	self.calcDeriv_PIRK_L3_ABarK_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 8 - derivL3_n:')
	self:printBuf(self.derivL3_nObj)
end
	
	-- derivL3_1 = PIRK L3 part1 of UBuf (ABar_IJ, K) of U1
	self.calcDeriv_PIRK_L3_ABarK_KernelObj(self.solverBuf, self.derivL3_1, self.U1, self.cellBuf)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 8 - derivL3_1:')
	self:printBuf(self.derivL3_1Obj)
end

	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_ABarK_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (ABar_IJ, K)
	self:applyBoundaryTo(self.UNext, self.PIRK_ABarK_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 8 - UNext:')
	self:printBuf(self.UNextObj)
end

-- step 9: LambdaBar^I

	-- UTemp = UNext except LambdaBar^I come from U1
	self.copyLambdaBarKernelObj(self.solverBuf, self.UTemp, self.UNext, self.U1)

	-- derivL2_n = PIRK L2 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)
	
	-- derivL2_next = PIRK L2 part2 (LambdaBar^I) based on UNext
	self.calcDeriv_PIRK_L2_LambdaBar_KernelObj(self.solverBuf, self.derivL2_next, self.UNext, self.cellBuf)
	
	-- derivL3_n = PIRK L3 part2 (LambdaBar^I) based on UBuf
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)
	
	-- derivL3_1 = PIRK L3 part2 (LambdaBar^I) based on UTemp
	self.calcDeriv_PIRK_L3_LambdaBar_KernelObj(self.solverBuf, self.derivL3_1, self.UTemp--[[self.U1--]], self.cellBuf)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_LambdaBar_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (LambdaBar^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_LambdaBar_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 9 - UNext:')
	self:printBuf(self.UNextObj)
end

-- step 10: B^I

	-- derivL2_n = PIRK L2 part3 (B^I) based on UBuf
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_n, self.UBuf, self.cellBuf)
	
	-- derivL2_next = PIRK L2 part3 (B^I) based on UNext
	self.calcDeriv_PIRK_L2_B_KernelObj(self.solverBuf, self.derivL2_next, self.UNext, self.cellBuf)
	
	-- derivL3_n = PIRK L3 part3 (B^I) based on UBuf
	self.calcDeriv_PIRK_L3_B_KernelObj(self.solverBuf, self.derivL3_n, self.UBuf, self.cellBuf)
	
	-- derivL3_1 = PIRK L3 part3 (B^I) based on U1
	self.calcDeriv_PIRK_L3_B_KernelObj(self.solverBuf, self.derivL3_1, self.U1, self.cellBuf)
	
	-- UNext = UBuf + 0.5 * dt * (derivL2_n + derivL2_next + derivL3_n + derivL3_1)
	self.PIRK_Eq4_B_KernelObj(self.solverBuf, self.UNext, self.UBuf, self.derivL2_n, self.derivL2_next, self.derivL3_n, self.derivL3_1, real(dt))
	
	-- apply boundary to UNext updated fields (B^I)
	self:applyBoundaryTo(self.UNext, self.PIRK_B_BoundaryKernelObjs)

if cmdline.printBufs then
	print()
	print('UBuf PIRK step 10 - UNext:')
	self:printBuf(self.UNextObj)
end


-- done

	-- constrain metric det of UNext
	-- remove trace of UNext
	-- copy to UBuf

	self.cmds:enqueueCopyBuffer{dst=self.UBuf, src=self.UNext, size=bufferSize}

	-- and last, constrain U
	--self:constrainU()
	-- this will call boundary() which I don't want
	-- I also don't want to calculate RBar_IJ, H, M^I
	-- I just want to update det(gammaBar_IJ)=1 and ABar^I_I=0

	assert(FakeIntegrator:isa(self.integrator))

	-- this just tests nans.  it would call integrator:integrate() but that is set to nothing.
	BSSNOKFiniteDifferencePIRKSolver.super.step(self, dt)
end

return BSSNOKFiniteDifferencePIRKSolver
