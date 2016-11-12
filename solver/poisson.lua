local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'


local Poisson = class()

Poisson.gaussSeidelMaxIters = 20

function Poisson:init(solver)
	self.solver = solver
end

function Poisson:createBuffers()
	local solver = self.solver
	solver:clalloc('ePotBuf', solver.volume * ffi.sizeof(solver.app.real))
end

function Poisson:addConvertToTexs()
	self.solver:addConvertToTex{
		name = 'ePot',
		vars = {{['0'] = 'value = buf[index];'}},
	}
end

function Poisson:getSolverCode()
	return require 'processcl'(
		table{
			file['solver/poisson.cl'],
			self.extraCode or '',
		}:concat'\n',
		table(self:getCodeParams(), {solver=self.solver}))
end

function Poisson:refreshSolverProgram()
	local solver = self.solver
	solver.initPotentialKernel = solver.solverProgram:kernel('initPotential', solver.ePotBuf, solver.UBuf)

	solver.solvePoissonKernel = solver.solverProgram:kernel('solvePoisson', solver.ePotBuf, solver.UBuf)
end

function Poisson:refreshBoundaryProgram()
	local solver = self.solver
	solver.potentialBoundaryProgram, solver.potentialBoundaryKernel =
		solver:createBoundaryProgramAndKernel{
			type = 'real',
			methods = table.map(solver.boundaryMethods, function(v,k)
				return solver.app.boundaryMethods[1+v[0]], k
			end),
		}
	solver.potentialBoundaryKernel:setArg(0, solver.ePotBuf)
end

function Poisson:resetState()
	local solver = self.solver
	if solver.useGravity then 
		solver.app.cmds:enqueueNDRangeKernel{kernel=solver.initPotentialKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
		solver:potentialBoundary()
		self:relax()
	end
	
	-- TODO
	-- add potential energy into total energy
	-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
end

function Poisson:relax()
	local solver = self.solver
	for i=1,self.gaussSeidelMaxIters do
		solver:potentialBoundary()
		solver.app.cmds:enqueueNDRangeKernel{kernel=solver.solvePoissonKernel, dim=solver.dim, globalSize=solver.gridSize:ptr(), localSize=solver.localSize:ptr()}
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

		function template:createBuffers()
			template.super.createBuffers(self)
			self[field]:createBuffers()
		end

		function template:addConvertToTexs()
			template.super.addConvertToTexs(self)
			self[field]:addConvertToTexs()
		end

		function template:getSolverCode()
			return table{
				template.super.getSolverCode(self),
				self[field]:getSolverCode(),
			}:concat'\n'
		end

		function template:refreshSolverProgram()
			template.super.refreshSolverProgram(self)
			self[field]:refreshSolverProgram()
		end

		function template:refreshBoundaryProgram()
			template.super.refreshBoundaryProgram(self)
			self[field]:refreshBoundaryProgram()
		end

		function template:resetState()
			template.super.resetState(self)
			self[field]:resetState()
		end

		function template:potentialBoundary()
			self:applyBoundaryToBuffer(self.potentialBoundaryKernel)
		end

		return template
	end
end

return Poisson
