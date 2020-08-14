local class = require 'ext.class'
local table = require 'ext.table'
local GridSolver = require 'hydro.solver.gridsolver'

local EinsteinFiniteDifferenceSolver = class(GridSolver)

-- TODO make a gui variable for numGhost
-- hmm, can I do that without rebuilding solverProgram every time it changes?
-- probably not, courtesy of boundary
-- in fact, how would boundary work with numGhost!=2?
-- esp mirror boundary conditions?
EinsteinFiniteDifferenceSolver.numGhost = 3

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
	self.calcDerivKernelObj.obj:setArg(3, self.cellBuf)
end

function EinsteinFiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivKernelObj()
end

-- [=[ commenting out to try to reduce build time
-- only set these for certain types ... 
function EinsteinFiniteDifferenceSolver:createDisplayComponents()
	EinsteinFiniteDifferenceSolver.super.createDisplayComponents(self)
	self:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted',
		code = [[
	const global <?=eqn.cons_t?>* U = buf + index;
	value->vreal = real3_weightedLen(value->vreal3, calc_gamma_ll(U, x));
]],
	})
	self:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted',
		code = [[
	const global <?=eqn.cons_t?>* U = buf + index;
	value->vreal = sym3_dot(value->vsym3, calc_gamma_uu(U, x));]],
	})
end
--]=]

return EinsteinFiniteDifferenceSolver
