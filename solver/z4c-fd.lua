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


local Z4cFiniteDifferenceSolver = class(GridSolver)
Z4cFiniteDifferenceSolver.name = 'Z4c_FiniteDifference'

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
Z4cFiniteDifferenceSolver.numGhost = 2

Z4cFiniteDifferenceSolver.eqnName = 'z4c-fd'

function Z4cFiniteDifferenceSolver:init(...)
	Z4cFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function Z4cFiniteDifferenceSolver:postInit(...)
	Z4cFiniteDifferenceSolver.super.postInit(self, ...)
	if not require 'int.be'.is(self.integrator) then
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		print("!! you're using a finite difference solver without an implicit integrator !!")
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	end
end

function Z4cFiniteDifferenceSolver:refreshSolverProgram()
	Z4cFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(1, self.UBuf)
end

function Z4cFiniteDifferenceSolver:refreshCalcDTKernel() end
function Z4cFiniteDifferenceSolver:calcDT() return self.fixedDT end

function Z4cFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernelObj(derivBuf)
end

function Z4cFiniteDifferenceSolver:getDisplayInfosForType()
	local t = Z4cFiniteDifferenceSolver.super.getDisplayInfosForType(self)

	-- hmm, only works with U ... so it only applies to U ...
	table.insert(t.real3, {
		name = ' norm weighted',
		code = [[
	{
		sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
		real exp_4phi = 1. / calc_exp_neg4phi(U);
		sym3 gamma_ll = sym3_scale(gammaBar_ll, exp_4phi);
		*value = real3_weightedLen(*valuevec, gammaBar_ll);
	}
]],
	})

	-- hmm, how to do the weighting stuff with gammaBar_ll ... 
	-- also, how to determine which metric to raise by ... gamma vs gammaBar
	table.insert(t.sym3, {
		name = ' tr weighted',
		code = [[
	real det_gamma_ll = calc_det_gamma_ll(U, x);
	*value = sym3_dot(U->gammaBar_uu, *valuesym3) / det_gamma_ll;
]],
	})

	return t
end

return Z4cFiniteDifferenceSolver

