local table = require 'ext.table'
local file = require 'ext.file'
local class = require 'ext.class'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local template = require 'template'

--[[
local CLGMRES = require 'solver.cl.gmres'

-- TODO combine this with int/be.lua
local ThisGMRES = class(CLGMRES)

function ThisGMRES:newBuffer(name)
	if not self.cache then self.cache = {} end
	local cached = self.cache[name]
	if cached then return cached end
	cached = ThisGMRES.super.newBuffer(self, name)
	cached:fill()
	self.cache[name] = cached
	return cached
end
--]]

local Poisson = class()

function Poisson:getPotBufType()
	return self.solver.eqn.cons_t
end

function Poisson:getPotBuf()
	return self.solver.UBuf
end

Poisson.potentialField = 'ePot'

function Poisson:init(solver)
	self.solver = solver

	-- hmm, should this go in refreshGridSize?
	-- [[ poisson
	--]]
	--[[ gmres
	this.linearSolver = ThisGMRES(linearSolverArgs)
	--]]
end

Poisson.stopOnEpsilon = true
Poisson.stopEpsilon = 1e-2
Poisson.maxIters = 20

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
	self.initPotentialKernel = solver.solverProgram:kernel('initPotential', self:getPotBuf())
	self.solvePoissonKernel = solver.solverProgram:kernel('solvePoisson', self:getPotBuf())
	if self.stopOnEpsilon then
		self.solvePoissonKernel:setArg(1, solver.reduceBuf)
	end
end

function Poisson:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to Poisson:potentialField
	solver.potentialBoundaryProgram, solver.potentialBoundaryKernel =
		solver:createBoundaryProgramAndKernel{
			type = self:getPotBufType(),
			methods = table.map(solver.boundaryMethods, function(v)
				return (select(2, next(solver.boundaryOptions[1+v[0]])))
			end),
			assign = function(a,b)
				return a..'.'..self.potentialField..' = '..b..'.'..self.potentialField
			end,
		}
	solver.potentialBoundaryKernel:setArg(0, self:getPotBuf())
end

function Poisson:resetState()
	local solver = self.solver
	solver.app.cmds:enqueueNDRangeKernel{kernel=self.initPotentialKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
	solver:potentialBoundary()
	self:relax()
end

function Poisson:relax()
	local solver = self.solver
	for i=1,self.maxIters do
		self.lastIter = i
		solver.app.cmds:enqueueNDRangeKernel{kernel=self.solvePoissonKernel, dim=solver.dim, globalSize=solver.globalSize:ptr(), localSize=solver.localSize:ptr()}
		solver:potentialBoundary()

		if self.stopOnEpsilon then
			local err = solver.reduceSum()
			self.lastEpsilon = err
			--print('gauss seidel iter '..i..' err '..err)
			if err <= self.stopEpsilon then break end
		end
	end
end

function Poisson:updateGUI()
	-- TODO unique name for other Poisson solvers?
	ig.igPushIdStr'Poisson solver'
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

--[[
static function
called with : (to get the correct subclass)
used as behavior template
field - which field in 'self' to store this behavior object 
enableField - which field in 'self' to toggle the behavior
--]]
function Poisson:createBehavior(field, enableField)
	local subclass = self
	return function(parent)
		local templateClass = class(parent)

		function templateClass:init(args)
			if enableField then
				self[enableField] = not not args[enableField]
			end

			-- TODO in refreshGrid?
			-- or should I always build one of these? 
			--if not enableField or not self[enableField] then
			self[field] = subclass(self)
			--end

			-- init is gonna call
			templateClass.super.init(self, args)
		end

		function templateClass:getSolverCode()
			return table{
				templateClass.super.getSolverCode(self),
				self[field]:getSolverCode(),
			}:concat'\n'
		end

		function templateClass:refreshBoundaryProgram()
			templateClass.super.refreshBoundaryProgram(self)
			self[field]:refreshBoundaryProgram()
		end

		function templateClass:refreshSolverProgram()
			templateClass.super.refreshSolverProgram(self)
			self[field]:refreshSolverProgram()
		end

		-- TODO
		-- for Euler, add potential energy into total energy
		-- then MAKE SURE TO SUBTRACT IT OUT everywhere internal energy is used
		function templateClass:resetState()
			templateClass.super.resetState(self)
			if not enableField or self[enableField] then
				self[field]:resetState()
			end
		end

		function templateClass:potentialBoundary()
			self:applyBoundaryToBuffer(self.potentialBoundaryKernel)
		end

		return templateClass
	end
end

return Poisson
