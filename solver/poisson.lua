local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'


local Poisson = class()

Poisson.gaussSeidelMaxIters = 20

function Poisson:init(solver)
	self.solver = solver
end

function Poisson:getSolverCode()
	return require 'template'(
		table{
			file['solver/poisson.cl'],
			self.extraCode or '',
		}:concat'\n',
		table(self:getCodeParams(), {
			eqn = self.solver.eqn,
			solver = self.solver,
		}))
end

function Poisson:refreshSolverProgram()
	local solver = self.solver
	solver.initPotentialKernel = solver.solverProgram:kernel('initPotential', solver.UBuf)
	solver.solvePoissonKernel = solver.solverProgram:kernel('solvePoisson', solver.UBuf)
end

function Poisson:refreshBoundaryProgram()
	local solver = self.solver
	-- TODO only apply to the ePot field
	solver.potentialBoundaryProgram, solver.potentialBoundaryKernel =
		solver:createBoundaryProgramAndKernel{
			type = solver.eqn.cons_t,
			methods = table.map(solver.boundaryMethods, function(v,k)
				return solver.app.boundaryMethods[1+v[0]], k
			end),
			assign = function(a,b)
				return a..'.ePot = '..b..'.ePot'
			end,
		}
	solver.potentialBoundaryKernel:setArg(0, solver.UBuf)
end

function Poisson:resetState()
	local solver = self.solver
	solver.app.cmds:enqueueNDRangeKernel{kernel=solver.initPotentialKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
	solver:potentialBoundary()
	self:relax()
end

function Poisson:relax()
	local solver = self.solver
	for i=1,self.gaussSeidelMaxIters do
		solver.app.cmds:enqueueNDRangeKernel{kernel=solver.solvePoissonKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
		solver:potentialBoundary()
	end
end

-- static function
-- called with : (to get the correct subclass)
-- used as behavior template
function Poisson:createBehavior(field, enableField)
	local subclass = self
	return function(parent)
		local template = class(parent)

		function template:init(args)
			if enableField then
				self[enableField] = not not args[enableField]
			end

			-- TODO in refreshGrid
			if not enableField or not self[enableField] then
				self[field] = subclass(self)
			end

			-- init is gonna call
			template.super.init(self, args)
		end

		function template:getSolverCode()
			return table{
				template.super.getSolverCode(self),
				self[field]:getSolverCode(),
			}:concat'\n'
		end

		function template:refreshBoundaryProgram()
			template.super.refreshBoundaryProgram(self)
			self[field]:refreshBoundaryProgram()
		end

		function template:refreshSolverProgram()
			template.super.refreshSolverProgram(self)
			self[field]:refreshSolverProgram()
		end

		-- TODO
		-- for Euler, add potential energy into total energy
		-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
		function template:resetState()
			template.super.resetState(self)
			if not enableField or self[enableField] then
				self[field]:resetState()
			end
		end

		function template:potentialBoundary()
			self:applyBoundaryToBuffer(self.potentialBoundaryKernel)
		end

		return template
	end
end

return Poisson
