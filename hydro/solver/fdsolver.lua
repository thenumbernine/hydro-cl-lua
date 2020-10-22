-- cell-centered finite-difference solver
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local GridSolver = require 'hydro.solver.gridsolver'

local FiniteDifferenceSolver = class(GridSolver)
FiniteDifferenceSolver.name = 'FiniteDifference'

FiniteDifferenceSolver.numGhost = 2

function FiniteDifferenceSolver:createBuffers()
	FiniteDifferenceSolver.super.createBuffers(self)

	-- while this is the same size as the finite-volume fluxBuf
	-- that buffer stored the flux at interfaces of each side
	-- while this stores the flux in each direction at the cell center
	self:clalloc('fluxBuf', self.eqn.cons_t, self.numCells * self.dim)
end

function FiniteDifferenceSolver:initCodeModules()
	FiniteDifferenceSolver.super.initCodeModules(self)
	self.modules:add{
		name = 'calcFluxAtCell,calcDerivFiniteDifference',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'coord.cell_t',
			'SETBOUNDS',
			'cell_x',
			'solver.macros',
			'fluxFromCons',
		},
		code = template(file['hydro/solver/calcDerivFD.cl'], {solver=self, eqn=self.eqn}),
	}
	self.solverModulesEnabled['calcFluxAtCell,calcDerivFiniteDifference'] = true
end

function FiniteDifferenceSolver:refreshSolverProgram()
	FiniteDifferenceSolver.super.refreshSolverProgram(self)

	self.calcFluxAtCellKernelObj = self.solverProgramObj:kernel('calcFluxAtCell', self.solverBuf, self.fluxBuf, self.UBuf, self.cellBuf)

	self.calcDerivFiniteDifferenceKernelObj = self.solverProgramObj:kernel'calcDerivFiniteDifference'
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(3, self.cellBuf)
end

function FiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcFluxAtCellKernelObj()

	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFiniteDifferenceKernelObj()
end

return FiniteDifferenceSolver 
