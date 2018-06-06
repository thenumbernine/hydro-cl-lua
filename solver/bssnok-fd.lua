local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local GridSolver = require 'solver.gridsolver'

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym



local BSSNOKFiniteDifferenceSolver = class(GridSolver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOK_FiniteDifference'

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
BSSNOKFiniteDifferenceSolver.numGhost = 2

BSSNOKFiniteDifferenceSolver.eqnName = 'bssnok-fd'

function BSSNOKFiniteDifferenceSolver:init(...)
	BSSNOKFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function BSSNOKFiniteDifferenceSolver:postInit(...)
	BSSNOKFiniteDifferenceSolver.super.postInit(self, ...)
	if not require 'int.be'.is(self.integrator) then
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		print("!! you're using a finite difference solver without an implicit integrator !!")
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print(require 'ext.tolua'(self.integrator))	
	end
end

function BSSNOKFiniteDifferenceSolver:refreshSolverProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(1, self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceSolver:calcDT() return self.fixedDT end

function BSSNOKFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernelObj(derivBuf)
end

function BSSNOKFiniteDifferenceSolver:getDisplayInfosForType()
	local t = BSSNOKFiniteDifferenceSolver.super.getDisplayInfosForType(self)

	-- hmm, only works with U ... so it only applies to U ...
	table.insert(t.real3, {
		name = ' norm weighted',
		code = '*value = real3_weightedLen(*valuevec, U->gammaTilde_ll) / calc_exp_neg4phi(U);',
	})

	-- hmm, how to do the weighting stuff with gammaTilde_ll ... 
	-- also, how to determine which metric to raise by ... gamma vs gammaTilde
	table.insert(t.sym3, {
		name = ' tr weighted',
		code = '*value = sym3_dot(U->gammaTilde_uu, *valuesym3) / calc_det_gamma(U, x);',
	})

	return t
end

return BSSNOKFiniteDifferenceSolver
