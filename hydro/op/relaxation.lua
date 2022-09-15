local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local ig = require 'imgui'
local ffi = require 'ffi'
local CLBuffer = require 'cl.obj.buffer'

local half = require 'cl.obj.half'
local toreal, fromreal = half.toreal, half.fromreal


local Relaxation = class()

Relaxation.name = 'relaxation'

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

	require 'hydro.code.symbols'(self, self:getSymbolFields())

	self.writeBufObj = CLBuffer{
		env = solver.app.env,
		name = 'writeBuf',
		type = solver.app.real,
		count = solver.numCells,
	}
end

function Relaxation:getSymbolFields()
	return table{
		'comments',
		'initPotential',
		'solveJacobi',
		'copyWriteToPotentialNoGhost',
		'setReduceToPotentialSquared',
	}
end

function Relaxation:initCodeModules()
	local solver = self.solver
	solver.modules:addFromMarkup{
		code = solver.eqn:template(
			file(self.solverCodeFile):read(),
			table(self.symbols, {
				op = self,
			})
		),
	}
	-- provided in poisson_jacobi and poisson_krylov's poisson.cl
	solver.solverModulesEnabled[self.symbols.initPotential] = true
	solver.solverModulesEnabled[self.symbols.copyWriteToPotentialNoGhost] = true
	solver.solverModulesEnabled[self.symbols.setReduceToPotentialSquared] = true
	-- provided in their extra code:
	solver.solverModulesEnabled[self.symbols.solveJacobi] = true
end

function Relaxation:refreshSolverProgram()
	local solver = self.solver
	self.initPotentialKernelObj = solver.solverProgramObj:kernel(self.symbols.initPotential, solver.solverBuf, self:getPotBuf())
	
	-- TODO names?  this seems subclass-specific
	self.solveJacobiKernelObj = solver.solverProgramObj:kernel(self.symbols.solveJacobi, self.solver.solverBuf, self.writeBufObj, self:getPotBuf(), solver.cellBuf)
	if self.stopOnEpsilon then
		self.solveJacobiKernelObj.obj:setArg(4, solver.reduceBuf)
	end
	
	-- TODO if we are using this 'writeBufObj'
	-- and this kernel has to copy to the potential field
	-- then ... why store the potential field at all?
	-- except maybe in cases where it is used elsewhere, like ePot <-> ETotal
	-- but ... with noDiv, why store vPot at all?
	self.copyWriteToPotentialNoGhostKernelObj = solver.solverProgramObj:kernel{
		name = self.symbols.copyWriteToPotentialNoGhost,
		setArgs={
			solver.solverBuf,
			solver.UBuf,
			self.writeBufObj
		},
		domain=solver.domainWithoutBorder,
	}

	self.setReduceToPotentialSquaredKernelObj = solver.solverProgramObj:kernel{
		name = self.symbols.setReduceToPotentialSquared,
		setArgs = {
			solver.solverBuf,
			solver.reduceBuf,
			solver.UBuf,
		},
		domain = solver.domainWithoutBorder,
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
			programNameSuffix = '_'..self.symbolPrefix,
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
			local residual = math.sqrt(fromreal(solver.reduceSum()) / tonumber(solver.volumeWithoutBorder))
			local lastResidual = self.lastResidual
			self.lastResidual = residual
			if self.verbose then
				self.setReduceToPotentialSquaredKernelObj()
				local xNorm = math.sqrt(fromreal(solver.reduceSum()) / tonumber(solver.volumeWithoutBorder))
				self.copyPotentialToReduceKernelObj()
				local xmin = fromreal(solver.reduceMin())
				self.copyPotentialToReduceKernelObj()
				local xmax = fromreal(solver.reduceMax())
				io.stderr:write(table{i-1, residual, xNorm, xmin, xmax}:mapi(tostring):concat'\t','\n')
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
	ig.igPushID_Str(self.symbolPrefix..' solver')
	-- TODO name from 'field' / 'enableField', though those aren't properties of Relaxation
	if ig.igCollapsingHeader(self.name..' solver') then
		if ig.luatableTooltipCheckbox('stop on epsilon', self, 'stopOnEpsilon') then
			-- TODO just recompile the poisson program?
			self.solver:refreshSolverProgram()
		end
		ig.igSameLine()
		ig.luatableTooltipInputFloatAsText('epsilon', self, 'stopEpsilon')
		ig.luatableTooltipInputInt('maxiter', self, 'maxIters')
		-- if it doesn't have to stop on epsilon then it doesn't calculate the x norm
		if self.stopOnEpsilon then
			ig.igText('residual = '..tostring(self.lastResidual))
		end
		ig.igText('iter = '..tostring(self.lastIter))
	end
	ig.igPopID()
end

return Relaxation
