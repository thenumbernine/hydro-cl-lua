local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local template = require 'template'
local ig = require 'ffi.imgui'
local ffi = require 'ffi'
local tooltip = require 'hydro.tooltip'
local CLBuffer = require 'cl.obj.buffer'

local half = require 'hydro.half'
local toreal, fromreal = half.toreal, half.fromreal


local Relaxation = class()

Relaxation.name = 'Relaxation'

-- scalar type of our vectors -- real or cplx 
-- TODO can be inferred from the type
Relaxation.scalar = 'real'

-- field we are solving for
Relaxation.potentialField = 'ePot'

Relaxation.stopOnEpsilon = true
Relaxation.stopEpsilon = 1e-10
Relaxation.maxIters = cmdline.selfGravPoissonMaxIter or 20	-- 10 * numreals
Relaxation.verbose = false

-- child class needs to provide this:
Relaxation.solverCodeFile = nil


-- type of the buffer holding the potential field
-- TODO can this be inferred? not at the moment, because SolverBase:clalloc is not accepting type info, and passing 'real' to the cl lib
function Relaxation:getPotBufType() 
	--return self.solver.UBufObj.type
	-- should be the same:
	return self.solver.eqn.symbols.cons_t
end

-- buffer holding the potential field
function Relaxation:getPotBuf() 
	return self.solver.UBuf 
end


function Relaxation:init(args)
	local solver = assert(args.solver)
	self.solver = solver
	
	self.potentialField = args.potentialField
	self.verbose = args.verbose

	-- this assumes 'self.name' is the class name
	self.name = solver.app:uniqueName(self.name)

	self.writeBufObj = CLBuffer{
		env = solver.app.env,
		name = 'writeBuf',
		type = solver.app.real,
		count = solver.numCells,
	}
end

function Relaxation:initCodeModules(solver)
	local name = 'op.Relaxation-'..self.name
	solver.modules:add{
		name = name,
		depends = {
			'SETBOUNDS_NOGHOST',
		},
		code = solver.eqn:template(file[self.solverCodeFile], {op = self}),
	}
	solver.solverModulesEnabled[name] = true
end

function Relaxation:refreshSolverProgram()
	local solver = self.solver
	self.initPotentialKernelObj = solver.solverProgramObj:kernel('initPotential'..self.name, solver.solverBuf, self:getPotBuf())
	self.solveJacobiKernelObj = solver.solverProgramObj:kernel('solveJacobi'..self.name, self.solver.solverBuf, self.writeBufObj, self:getPotBuf(), solver.cellBuf)
	if self.stopOnEpsilon then
		self.solveJacobiKernelObj.obj:setArg(4, solver.reduceBuf)
	end
	self.copyWriteToPotentialNoGhostKernelObj = solver.solverProgramObj:kernel{
		name='copyWriteToPotentialNoGhost'..self.name,
		setArgs={
			solver.solverBuf,
			solver.UBuf,
			self.writeBufObj
		},
		domain=solver.domainWithoutBorder,
	}

	self.setReduceToPotentialSquaredKernelObj = solver.solverProgramObj:kernel{
		name='setReduceToPotentialSquared'..self.name,
		setArgs={
			solver.solverBuf,
			solver.reduceBuf,
			solver.UBuf,
		},
		domain=solver.domainWithoutBorder,
	}
	
	self.lastResidual = 0
end

function Relaxation:refreshBoundaryProgram()
	local solver = self.solver
	-- only applies the boundary conditions to Relaxation:potentialField
	self.potentialBoundaryProgramObj, self.potentialBoundaryKernelObjs =
		solver:createBoundaryProgramAndKernel{
			type = self:getPotBufType(),
			methods = solver.boundaryMethods,
			fields = {self.potentialField},
			programNameSuffix = '-'..self.name,
		}
	for _,obj in ipairs(self.potentialBoundaryKernelObjs) do
		obj.obj:setArg(1, self:getPotBuf())
		obj.obj:setArg(2, solver.cellBuf)
	end
end

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
	
		-- writes new potentialField to writeBuf
		-- stores deltas in reduceBuf
		self.solveJacobiKernelObj()

		-- copy new values back from writeBuf to UBuf potentialField
		-- copies square of increment into reduceBuf
		self.copyWriteToPotentialNoGhostKernelObj()

		-- apply boundary to UBuf.potentialField
		self:potentialBoundary()

		if self.stopOnEpsilon then
			local residual = math.sqrt(solver.reduceSum() / tonumber(solver.volumeWithoutBorder))
			local lastResidual = self.lastResidual	
			self.lastResidual = residual
			if self.verbose then
				self.setReduceToPotentialSquaredKernelObj()
				local xNorm = math.sqrt(solver.reduceSum() / tonumber(solver.volumeWithoutBorder))
				self.copyPotentialToReduceKernelObj()
				local xmin = fromreal(solver.reduceMin())
				self.copyPotentialToReduceKernelObj()
				local xmax = fromreal(solver.reduceMax())
				io.stderr:write(table{i-1, residual, xNorm, xmin, xmax}:map(tostring):concat'\t','\n')
			end
			
			-- TODO compare residual
			if math.abs(residual) <= self.stopEpsilon then break end
		end
	end
end

function Relaxation:potentialBoundary()
	self.solver:applyBoundaryToBuffer(self.potentialBoundaryKernelObjs)
end

function Relaxation:updateGUI()
	ig.igPushIDStr(self.name..' solver')
	-- TODO name from 'field' / 'enableField', though those aren't properties of Relaxation
	if ig.igCollapsingHeader(self.name..' solver') then
		if tooltip.checkboxTable('stop on epsilon', self, 'stopOnEpsilon') then
			-- TODO just recompile the poisson program?
			self.solver:refreshSolverProgram()
		end
		ig.igSameLine()
		tooltip.numberTable('epsilon', self, 'stopEpsilon')
		tooltip.intTable('maxiter', self, 'maxIters')
		-- if it doesn't have to stop on epsilon then it doesn't calculate the x norm
		if self.stopOnEpsilon then
			ig.igText('residual = '..tostring(self.lastResidual))
		end
		ig.igText('iter = '..tostring(self.lastIter))
	end
	ig.igPopID()
end

return Relaxation
