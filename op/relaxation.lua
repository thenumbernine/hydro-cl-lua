local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local template = require 'template'

local Relaxation = class()

-- field we are solving for
Relaxation.potentialField = 'ePot'

Relaxation.name = 'Relaxation'

local ident = 1

function Relaxation:init(args)
	self.solver = assert(args.solver)
	self.potentialField = args.potentialField
	self.suffix = ''..ident
	ident = ident + 1
end

-- scalar type of our vectors -- real or cplx 
Relaxation.scalar = 'real'

-- type of the buffer holding the potential field
-- TODO can this be inferred?
function Relaxation:getPotBufType() return self.solver.eqn.cons_t end

-- buffer holding the potential field
function Relaxation:getPotBuf() return self.solver.UBuf end

Relaxation.stopOnEpsilon = true
-- [[ realtime
Relaxation.stopEpsilon = 1e-12
Relaxation.maxIters = 20
--]]
--[[ safe
Relaxation.stopEpsilon = 1e-20
Relaxation.maxIters = 10000
--]]
--Relaxation.verbose = true

-- child class needs to provide this:
Relaxation.solverCodeFile = nil

function Relaxation:getSolverCode()
	return template(file[self.solverCodeFile], {op=self})
end

function Relaxation:refreshSolverProgram()
	local solver = self.solver
	self.initPotentialKernelObj = solver.solverProgramObj:kernel('initPotential'..self.name..self.suffix, self:getPotBuf())
	self.solveJacobiKernelObj = solver.solverProgramObj:kernel('solveJacobi'..self.name..self.suffix, self:getPotBuf())
	if self.stopOnEpsilon then
		self.solveJacobiKernelObj.obj:setArg(1, solver.reduceBuf)
	end
end

function Relaxation:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to Relaxation:potentialField
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
function Relaxation:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initPotentialKernelObj()
	self:potentialBoundary()
	self:relax()
end

function Relaxation:relax()
	local solver = self.solver
	for i=1,self.maxIters do
		self.lastIter = i
		self.solveJacobiKernelObj()
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

function Relaxation:potentialBoundary()
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernelObjs)
end

function Relaxation:updateGUI()
	-- TODO unique name for other Relaxation solvers?
	ig.igPushIDStr(self.name..' solver'..self.ident)
	-- TODO name from 'field' / 'enableField', though those aren't properties of Relaxation
	if ig.igCollapsingHeader(self.name..' solver') then
		if tooltip.checkboxTable('stop on epsilon', self, 'stopOnEpsilon') then
			-- TODO just recompile the poisson program?
			self.solver:refreshSolverProgram()
		end
		ig.igSameLine()
		tooltip.numberTable('epsilon', self, 'stopEpsilon')
		tooltip.intTable('maxiter', self, 'maxIters')
		-- if it doesn't have to stop on epsilon then it doesn't calculate the error
		if self.stopOnEpsilon then
			ig.igText('err = '..tostring(self.lastEpsilon))
		end
		ig.igText('iter = '..tostring(self.lastIter))
	end
	ig.igPopID()
end



return Relaxation
