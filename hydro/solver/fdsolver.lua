-- cell-centered finite-difference solver with rhs derived from flux calculations
local ffi = require 'ffi'
local table = require 'ext.table'
local path = require 'ext.path'
local GridSolver = require 'hydro.solver.gridsolver'

-- TODO TODO this should be fv_to_fd_solver because it takes the flux and turns it into finite-difference
local FiniteDifferenceFromFiniteVolumeSolver = GridSolver:subclass()
FiniteDifferenceFromFiniteVolumeSolver.name = 'FiniteDifferenceFromFiniteVolume'

FiniteDifferenceFromFiniteVolumeSolver.numGhost = 2

function FiniteDifferenceFromFiniteVolumeSolver:createBuffers()
	FiniteDifferenceFromFiniteVolumeSolver.super.createBuffers(self)

	-- while this is the same size as the finite-volume fluxBuf
	-- that buffer stored the flux at interfaces of each side
	-- while this stores the flux in each direction at the cell center
	-- TODO we only need this if our eqn uses calcFluxAtCell / calcDerivFiniteDifference
	self:clalloc('fluxBuf', self.eqn.symbols.cons_t, self.numCells * self.dim)
end

function FiniteDifferenceFromFiniteVolumeSolver:initCodeModules()
	FiniteDifferenceFromFiniteVolumeSolver.super.initCodeModules(self)
	self.modules:addFromMarkup(
		self.eqn:template(path'hydro/solver/fdsolver.cl':read())
	)
	
	self.solverModulesEnabled['calcFluxAtCell'] = true
	self.solverModulesEnabled['calcDerivFiniteDifference'] = true
end

function FiniteDifferenceFromFiniteVolumeSolver:refreshSolverProgram()
	FiniteDifferenceFromFiniteVolumeSolver.super.refreshSolverProgram(self)

	self.calcFluxAtCellKernelObj = self.solverProgramObj:kernel('calcFluxAtCell', self.solverBuf, self.fluxBuf, self.UBuf, self.cellBuf)

	self.calcDerivFiniteDifferenceKernelObj = self.solverProgramObj:kernel'calcDerivFiniteDifference'
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(2, self.fluxBuf)
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(3, self.cellBuf)
end

function FiniteDifferenceFromFiniteVolumeSolver:calcDeriv(derivBufObj, dt)
	self.calcFluxAtCellKernelObj()

	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFiniteDifferenceKernelObj()
end

return FiniteDifferenceFromFiniteVolumeSolver 
