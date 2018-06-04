local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local template = require 'template'

local Poisson = class()

Poisson.potentialField = 'ePot'

local ident = 1

function Poisson:init(args)
	self.solver = assert(args.solver)
	self.potentialField = args.potentialField
	self.suffix = ''..ident
	ident = ident + 1
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
	return table{
		template(file['op/poisson.cl'], {poisson=self}),
		self:getPoissonCode() or '',
	}:concat'\n'
end

function Poisson:refreshSolverProgram()
	local solver = self.solver
	self.initPoissonPotentialKernelObj = solver.solverProgramObj:kernel('initPoissonPotential'..self.suffix, self:getPotBuf())
	self.solvePoissonJacobiKernelObj = solver.solverProgramObj:kernel('solvePoissonJacobi'..self.suffix, self:getPotBuf())
	if self.stopOnEpsilon then
		self.solvePoissonJacobiKernelObj.obj:setArg(1, solver.reduceBuf)
	end
end

function Poisson:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to Poisson:potentialField
	self.potentialBoundaryProgramObj, self.potentialBoundaryKernelObjs =
		solver:createBoundaryProgramAndKernel{
			type = self:getPotBufType(),
			methods = table.map(solver.boundaryMethods, function(v)
				return (select(2, next(solver.boundaryOptions[v])))
			end),
			assign = function(a,b)
				return a..'.'..self.potentialField..' = '..b..'.'..self.potentialField
			end,
		}
	for _,obj in ipairs(self.potentialBoundaryKernelObjs) do
		obj.obj:setArg(0, self:getPotBuf())
	end
end

-- TODO
-- for Euler, add potential energy into total energy
-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
function Poisson:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initPoissonPotentialKernelObj()
	self:potentialBoundary()
	self:relax()
end

function Poisson:relax()
	local solver = self.solver
	for i=1,self.maxIters do
		self.lastIter = i
		self.solvePoissonJacobiKernelObj()
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
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernelObjs)
end

function Poisson:updateGUI()
	-- TODO unique name for other Poisson solvers?
	ig.igPushIDStr'Poisson behavior'
	-- TODO name from 'field' / 'enableField', though those aren't properties of Poisson
	if ig.igCollapsingHeader'Poisson solver' then
		if tooltip.checkboxTable('stop on epsilon', self, 'stopOnEpsilon') then
			-- TODO just recompile the poisson program?
			self.solver:refreshSolverProgram()
		end
		ig.igSameLine()
		tooltip.numberTable('epsilon', self, 'stopEpsilon')
		tooltip.intTable('maxiter', self, 'maxIters')
		if self.stopOnEpsilon then
			ig.igText('err = '..tostring(self.lastEpsilon))
		end
		ig.igText('iter = '..tostring(self.lastIter))
	end
	ig.igPopID()
end

return Poisson
