local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local template = require 'template'
local ig = require 'ffi.imgui'
local ffi = require 'ffi'
local tooltip = require 'tooltip'
local CLBuffer = require 'cl.obj.buffer'

local Relaxation = class()

Relaxation.name = 'Relaxation'

-- scalar type of our vectors -- real or cplx 
-- TODO can be inferred from the type
Relaxation.scalar = 'real'

-- field we are solving for
Relaxation.potentialField = 'ePot'

Relaxation.stopOnEpsilon = true
--[[ realtime
Relaxation.stopEpsilon = 1e-18
Relaxation.maxIters = 20
--]]
-- [[ safe
Relaxation.stopEpsilon = 1e-10
Relaxation.maxIters = cmdline.selfGravPoissonMaxIter or 20	-- 10 * numreals
--]]
Relaxation.verbose = false

-- child class needs to provide this:
Relaxation.solverCodeFile = nil


-- type of the buffer holding the potential field
-- TODO can this be inferred? not at the moment, because SolverBase:clalloc is not accepting type info, and passing 'real' to the cl lib
function Relaxation:getPotBufType() 
	return self.solver.UBufObj.type
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

function Relaxation:getCode()
	return template(file[self.solverCodeFile], {op = self})
end

function Relaxation:refreshSolverProgram()
	local solver = self.solver
	self.initPotentialKernelObj = solver.solverProgramObj:kernel('initPotential'..self.name, solver.solverBuf, self:getPotBuf())
	self.solveJacobiKernelObj = solver.solverProgramObj:kernel('solveJacobi'..self.name, self.solver.solverBuf, self.writeBufObj, self:getPotBuf())
	if self.stopOnEpsilon then
		self.solveJacobiKernelObj.obj:setArg(3, solver.reduceBuf)
	end
	self.copyWriteToPotentialNoGhostKernelObj = solver.solverProgramObj:kernel{
		name='copyWriteToPotentialNoGhost',
		setArgs={
			solver.solverBuf,
			solver.UBuf,
			self.writeBufObj
		},
		domain=solver.domainWithoutBorder,
	}

	
	self.squareKernelObj = solver.domain:kernel{
		name = 'Poisson_square'..self.name,
		header = solver.codePrefix,
		argsOut = {
			{name = 'y', type=solver.app.real, obj=true},
		},
		argsIn = {
			solver.solverBuf,
			{name = 'x', type=solver.app.real, obj=true},
		},
		body = [[
	if (OOB(0,0)) return;
	y[index] = x[index] * x[index];
]],
	}
	
	self.lastXNorm = 0
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
		obj.obj:setArg(1, self:getPotBuf())
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
		self.copyWriteToPotentialNoGhostKernelObj()

		-- apply boundary to UBuf.potentialField
		self:potentialBoundary()

		if self.stopOnEpsilon then
			-- copy write to reduce
			self.squareKernelObj(
				solver.reduceBuf,
				solver.solverBuf,
				self.writeBufObj)
			
			local xNorm = math.sqrt(solver.reduceSum()) / tonumber(solver.volumeWithoutBorder)
			local lastXNorm = self.lastXNorm	
			self.lastXNorm = xNorm
			if self.verbose then
				--print('relaxation iter '..i..' xNorm '..xNorm)

-- [[
self.copyPotentialToReduceKernelObj()
local xmin = solver.reduceMin()
self.copyPotentialToReduceKernelObj()
local xmax = solver.reduceMax()
io.stderr:write(table{i, xNorm, xmin, xmax}:map(tostring):concat'\t','\n')
--]]			
			end
			if math.abs(xNorm - lastXNorm) <= self.stopEpsilon then break end
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
			ig.igText('x norm = '..tostring(self.lastXNorm))
		end
		ig.igText('iter = '..tostring(self.lastIter))
	end
	ig.igPopID()
end

return Relaxation
