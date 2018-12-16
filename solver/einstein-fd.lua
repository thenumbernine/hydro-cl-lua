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


local EinsteinFiniteDifferenceSolver = class(GridSolver)

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
EinsteinFiniteDifferenceSolver.numGhost = 4

EinsteinFiniteDifferenceSolver.name = 'EinsteinFiniteDifference'

function EinsteinFiniteDifferenceSolver:init(...)
	EinsteinFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function EinsteinFiniteDifferenceSolver:postInit(...)
	EinsteinFiniteDifferenceSolver.super.postInit(self, ...)
	if not require 'int.be'.is(self.integrator) then
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		print("!! you're using a finite difference solver without an implicit integrator !!")
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	end
end

function EinsteinFiniteDifferenceSolver:refreshSolverProgram()
	EinsteinFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivKernelObj.obj:setArg(2, self.UBuf)
end

function EinsteinFiniteDifferenceSolver:refreshCalcDTKernel() end
function EinsteinFiniteDifferenceSolver:calcDT() return self.fixedDT end

function EinsteinFiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivKernelObj()
end

function EinsteinFiniteDifferenceSolver:getDisplayInfosForType()
	local t = EinsteinFiniteDifferenceSolver.super.getDisplayInfosForType(self)

	-- hmm, only works with U ... so it only applies to U ...
	table.insert(t.real3, {
		name = ' norm weighted',
		code = '	*value = real3_weightedLen(*value_real3, calc_gamma_ll(U, x));',
	})

	-- hmm, how to do the weighting stuff with gammaBar_ll ... 
	-- also, how to determine which metric to raise by ... gamma vs gammaBar
	table.insert(t.sym3, {
		name = ' tr weighted',
		code = '	*value = sym3_dot(*valuesym3, calc_gamma_uu(U));',
	})

	return t
end

return EinsteinFiniteDifferenceSolver
