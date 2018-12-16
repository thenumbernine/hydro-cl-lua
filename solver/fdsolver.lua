--[[
cell-centered finite-difference solver

this has no subclasses, so I'm thinking I will turn it into a characteristic-variable finite-difference solver


here's another idea for a finite-difference solver ...
a characteristic variable finite difference solver
this means multiply each grid cell U with left eigenvectors 
use the cell cons_t to create an eigen_t, use that to create a wave_t
then do finite-difference of the wave_t at each grid point

now to reconstruct the wave_t ... 
... we are now faced with the problem of using char var info to reconstruct the info required for a right-eigen-transform

..and then there's the issue of which direction do we perform the deconstruction?
I'm sure the answer is 'all', but that might vary depending on the underlying equation
For a 1D solver this makes no difference.


--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local GridSolver = require 'solver.gridsolver'

local FiniteDifferenceSolver = class(GridSolver)
FiniteDifferenceSolver.name = 'FiniteDifference'

FiniteDifferenceSolver.numGhost = 2

function FiniteDifferenceSolver:createBuffers()
	FiniteDifferenceSolver.super.createBuffers(self)

	-- while this is the same size as the finite-volume fluxBuf
	-- that buffer stored the flux at interfaces of each side
	-- while this stores the flux in each direction at the cell center
	self:clalloc('fluxBuf', self.numCells * self.dim * ffi.sizeof(self.eqn.cons_t))
end

function FiniteDifferenceSolver:getSolverCode()
	return table{
		FiniteDifferenceSolver.super.getSolverCode(self),
		template(file['solver/calcDerivFD.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function FiniteDifferenceSolver:refreshSolverProgram()
	FiniteDifferenceSolver.super.refreshSolverProgram(self)

	self.calcFluxAtCellKernelObj = self.solverProgramObj:kernel('calcFluxAtCell', self.solverBuf, self.fluxBuf, self.UBuf)

	self.calcDerivFiniteDifferenceKernelObj = self.solverProgramObj:kernel'calcDerivFiniteDifference'
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(0, self.solverBuf)
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(2, self.fluxBuf)
end

function FiniteDifferenceSolver:calcDeriv(derivBufObj, dt)
	self.calcFluxAtCellKernelObj()

	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFiniteDifferenceKernelObj()
end

return FiniteDifferenceSolver 
