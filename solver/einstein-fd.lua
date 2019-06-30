local class = require 'ext.class'
local table = require 'ext.table'
local GridSolver = require 'solver.gridsolver'

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

function EinsteinFiniteDifferenceSolver:refreshSolverProgram()
	EinsteinFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernelObj = self.solverProgramObj:kernel'calcDeriv'
	self.calcDerivKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivKernelObj.obj:setArg(2, self.UBuf)
end

function EinsteinFiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivKernelObj()
end

-- only set these for certain types ... 
function EinsteinFiniteDifferenceSolver:createDisplayComponents()
	EinsteinFiniteDifferenceSolver.super.createDisplayComponents(self)
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted',
		code = [[*value = real3_weightedLen(*value_real3, calc_gamma_ll(U, x));]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		norm = 'tr weighted',
		code = [[
		*value = sym3_dot(*value_sym3, calc_gamma_uu(U, x));]],
	})
end

return EinsteinFiniteDifferenceSolver
