local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local BSSNOKFiniteDifferenceEquation = require 'eqn.bssnok-fd'
local Solver = require 'solver.solver'

local xNames = table{'x', 'y', 'z'}

local BSSNOKFiniteDifferenceSolver = class(Solver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOKFiniteDifferenceSolver'

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
BSSNOKFiniteDifferenceSolver.numGhost = 2


function BSSNOKFiniteDifferenceSolver:init(...)
	BSSNOKFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function BSSNOKFiniteDifferenceSolver:createEqn(eqn)
	self.eqn = BSSNOKFiniteDifferenceEquation(self)
end

function BSSNOKFiniteDifferenceSolver:refreshSolverProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(1, self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:refreshInitStateProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshInitStateProgram(self)
	-- for now I'm conditionally testing this
	-- because "testbed - robust" doesn't create it, 
	-- and fills the connBarU buffer with noise anyways
	-- but it would be better if this condition was in the initial condition object
	if self.initStateProgramObj then
		self.init_connBarUKernelObj = self.initStateProgramObj:kernel('init_connBarU', self.UBuf)
	end
end

function BSSNOKFiniteDifferenceSolver:resetState()
	BSSNOKFiniteDifferenceSolver.super.resetState(self)
	-- same as above, I shouldn't have this condition here
	if self.init_connBarUKernelObj then
		self.init_connBarUKernelObj()
	end
	self:boundary()
end

function BSSNOKFiniteDifferenceSolver:getCalcDTCode() end
function BSSNOKFiniteDifferenceSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceSolver:calcDT() return self.fixedDT end

function BSSNOKFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernelObj(derivBuf)
end

return BSSNOKFiniteDifferenceSolver
