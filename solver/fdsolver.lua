--[[
cell-centered finite-difference solver
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local Solver = require 'solver.solver'

local FiniteDifferenceSolver = class(Solver)
FiniteDifferenceSolver.name = 'FiniteDifference'

function FiniteDifferenceSolver:createBuffers()
	FiniteDifferenceSolver.super.createBuffers(self)
	
	-- while this is the same size as the finite-volume fluxBuf
	-- that buffer stored the flux at interfaces of each side
	-- while this stores the flux in each direction at the cell center
	self:clalloc('fluxBuf', self.volume * self.dim * ffi.sizeof(self.eqn.cons_t))
end

function FiniteDifferenceSolver:getSolverCode()
	return table{
		FiniteDifferenceSolver.super.getSolverCode(self),
		template(file['solver/calcDerivFD.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

function FiniteDifferenceSolver:refreshSolverProgram()
	FiniteDifferenceSolver.super.refreshSolverProgram(self)

	self.calcFluxAtCellKernelObj = self.solverProgramObj:kernel('calcFluxAtCell', self.fluxBuf, self.UBuf)

	self.calcDerivFiniteDifferenceKernelObj = self.solverProgramObj:kernel'calcDerivFiniteDifference'
	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(1, self.fluxBuf)

	-- TODO put this in solver/solver.lua ?
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name='addSource', domain=self.domainWithoutBorder}
		self.addSourceKernelObj.obj:setArg(1, self.UBuf)
	end
end

function FiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self:boundary()

	self.calcFluxAtCellKernelObj()

	self.calcDerivFiniteDifferenceKernelObj.obj:setArg(0, derivBuf)
	self.calcDerivFiniteDifferenceKernelObj()
	
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj(derivBuf)
	end
end

return FiniteDifferenceSolver 
