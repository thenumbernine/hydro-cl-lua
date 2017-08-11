local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local template = require 'template'

local Poisson = class()

Poisson.potentialField = 'ePot'

function Poisson:init(args)
	self.solver = assert(args.solver)
	self.potentialField = args.potentialField
end

function Poisson:getPotBufType()
	return self.solver.eqn.cons_t
end

function Poisson:getPotBuf()
	return self.solver.UBuf
end

Poisson.stopOnEpsilon = true
-- [[ realtime
Poisson.stopEpsilon = 1e-12
Poisson.maxIters = 20
--]]
--[[ safe
Poisson.stopEpsilon = 1e-20
Poisson.maxIters = 10000
--]]
--Poisson.verbose = true

function Poisson:getSolverCode()
	return template(
		table{
			file['solver/poisson.cl'],
			self:getPoissonCode() or '',
		}:concat'\n',
		table(self:getCodeParams(), {
			poisson = self,
			solver = self.solver,
			eqn = self.solver.eqn,
		}))
end

function Poisson:refreshSolverProgram()
	local solver = self.solver
	self.initPoissonPotentialKernel = solver.solverProgram:kernel('initPoissonPotential', self:getPotBuf())
	self.solvePoissonJacobiKernel = solver.solverProgram:kernel('solvePoissonJacobi', self:getPotBuf())
	if self.stopOnEpsilon then
		self.solvePoissonJacobiKernel:setArg(1, solver.reduceBuf)
	end
end

function Poisson:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to Poisson:potentialField
	self.potentialBoundaryProgram, self.potentialBoundaryKernel =
		solver:createBoundaryProgramAndKernel{
			type = self:getPotBufType(),
			methods = table.map(solver.boundaryMethods, function(v)
				return (select(2, next(solver.boundaryOptions[v+1])))
			end),
			assign = function(a,b)
				return a..'.'..self.potentialField..' = '..b..'.'..self.potentialField
			end,
		}
	self.potentialBoundaryKernel:setArg(0, self:getPotBuf())
end

-- TODO
-- for Euler, add potential energy into total energy
-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
function Poisson:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.initPoissonPotentialKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	self:potentialBoundary()
	self:relax()
end

function Poisson:relax()
	local solver = self.solver
	for i=1,self.maxIters do
		self.lastIter = i
		solver.app.cmds:enqueueNDRangeKernel{kernel=self.solvePoissonJacobiKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
		self:potentialBoundary()

		if self.stopOnEpsilon then
			local err = solver.reduceSum() / tonumber(solver.volumeWithoutBorder)
			self.lastEpsilon = err
			if self.verbose then
				print('gauss seidel iter '..i..' err '..err)
			end
			if err <= self.stopEpsilon then break end
		end
	end
end

function Poisson:potentialBoundary()
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernel)
end

function Poisson:updateGUI()
	-- TODO unique name for other Poisson solvers?
	ig.igPushIdStr'Poisson behavior'
	-- TODO name from 'field' / 'enableField', though those aren't properties of Poisson
	if ig.igCollapsingHeader'Poisson solver' then
		if tooltip.checkboxTable('stop on epsilon', self, 'stopOnEpsilon') then
			-- TODO just recompile the poisson program?
			solver:refreshSolverProgram()
		end
		ig.igSameLine()
		tooltip.numberTable('epsilon', self, 'stopEpsilon')
		tooltip.intTable('maxiter', self, 'maxIters')
		if self.stopOnEpsilon then
			ig.igText('err = '..self.lastEpsilon)
		end
		ig.igText('iter = '..self.lastIter)
	end
	ig.igPopId()
end

return Poisson
