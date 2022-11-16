local ffi = require 'ffi'
local gl = require 'gl'
local ig = require 'imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local math = require 'ext.math'
local string = require 'ext.string'
local clnumber = require 'cl.obj.number'
local template = require 'template'
local vec3d = require 'vec-ffi.vec3d'
local vec3sz = require 'vec-ffi.vec3sz'
local roundup = require 'hydro.util.roundup'
local time, getTime = table.unpack(require 'hydro.util.time')
local SolverBase = require 'hydro.solver.solverbase'
local Struct = require 'hydro.code.struct'

local half = require 'cl.obj.half'
local toreal, fromreal = half.toreal, half.fromreal


local common = require 'hydro.common'
local minmaxs = common.minmaxs
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local function unpack(x)
	if x.unpack then return x:unpack() end
	return table.unpack(x)
end


local GridSolver = class(SolverBase)

GridSolver.name = 'gridsolver'

GridSolver.numGhost = 2

function GridSolver:getSymbolFields()
	return GridSolver.super.getSymbolFields(self):append{
		'OOB',
		'SETBOUNDS',
		'SETBOUNDS_NOGHOST',
		'updateCTU',
		'slopeLimiter',
		'calcLR',
	}
end

--[[
args:
	gridSize
	mins
	maxs
--]]
function GridSolver:initMeshVars(args)
	GridSolver.super.initMeshVars(self, args)
	
	-- same as equations
	-- but let equations/init conds add to the solver vars (as gui vars)
	-- then we can edit them without recompiling the kernels
	
	self.solverStruct.vars:append{
		{name='grid_dx', type='real3'},
		{name='gridSize', type='int4'},
		{name='stepsize', type='int4'},
		{name='numGhost', type='int'},
	}

	self.mins = vec3d(unpack(args.mins or {-1, -1, -1}))
	self.maxs = vec3d(unpack(args.maxs or {1, 1, 1}))

	self.initCondMins = vec3d(unpack(args.initCondMins or self.mins))
	self.initCondMaxs = vec3d(unpack(args.initCondMaxs or self.maxs))

	-- TODO OK this is a little ambiguous ...
	-- gridSize is the desired grid size
	-- self.gridSize is gonna be that, plus numGhost on either side ...
	-- so maybe I should rename 'self.gridSize' to 'self.gridSizeWithBorder'
	-- and 'self.gridSizeWithoutBorder' to 'self.gridSize'

	local gridSize = assert(args.gridSize, "solver expected gridSize")
	if type(gridSize) == 'number' then 
		self.gridSize = vec3sz(gridSize,1,1)
	else
		self.gridSize = vec3sz(unpack(gridSize))
	end

	for i=0,self.dim-1 do self.gridSize.s[i] = self.gridSize.s[i] + 2 * self.numGhost end
	for i=self.dim,2 do self.gridSize.s[i] = 1 end
end

--[[
args:
	boundary = boundary info
	userPLM
	useCTU
	slopeLimiter
--]]
function GridSolver:initObjs(args)
	-- [[ grab these before calling super, or else refreshGetULR will be initialized incorrectly at first
	self.usePLM = args.usePLM
	self.slopeLimiter = self.app.limiterNames:find(args.slopeLimiter) or 1
	self.useCTU = args.useCTU
	--]]
	
	GridSolver.super.initObjs(self, args)
	
	assert(not self.usePLM or self.fluxLimiter == 1, "are you sure you want to use flux and slope limiters at the same time?")

	-- TODO instead of boundaryMethods.xmin .xmax ... how about [1] ... [6] ? 
	self.boundaryMethods = {}
	for i=1,3 do
		for _,minmax in ipairs(minmaxs) do
			local var = xNames[i]..minmax
			local name = (args.boundary or {})[var] or 'freeflow'
			local boundaryObjArgs
			if type(name) == 'table' then
				boundaryObjArgs = name.args
				name = name.name
			end
			local boundaryClass = assert(self.boundaryOptionForName[name], "failed to find boundary method named "..name)
			self.boundaryMethods[var] = boundaryClass(boundaryObjArgs)
		end
	end
end


-- this is only used by GridSolver, but each kernel gets its own,
-- so TODO get rid of this and just use kernel's sizes
function GridSolver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	local localSize1d = math.min(maxWorkGroupSize, tonumber(self.gridSize:volume()))

	local localSize2d
	if self.dim == 3 then
		local localSizeX = math.min(tonumber(self.gridSize.x), 2^math.ceil(math.log(maxWorkGroupSize,2)/2))
		local localSizeY = maxWorkGroupSize / localSizeX
		localSize2d = vec3sz(localSizeX, localSizeY, 1)
	end

--	self.localSize = self.dim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
	-- TODO better than constraining by math.min(gridSize),
	-- look at which gridSizes have the most room, and double them accordingly, until all of maxWorkGroupSize is taken up
	local localSize = vec3sz(1,1,1)
	local rest = maxWorkGroupSize
	local localSizeX = math.min(tonumber(self.gridSize.x), 2^math.ceil(math.log(rest,2)/self.dim))
	localSize.x = localSizeX
	if self.dim > 1 then
		rest = rest / localSizeX
		if self.dim == 2 then
			localSize.y = math.min(tonumber(self.gridSize.y), rest)
		elseif self.dim == 3 then
			local localSizeY = math.min(tonumber(self.gridSize.y), 2^math.ceil(math.log(math.sqrt(rest),2)))
			localSize.y = localSizeY
			localSize.z = math.min(tonumber(self.gridSize.z), rest / localSizeY)
		end
	end

	-- this is grid size, but rounded up to the next localSize
	local globalSize = vec3sz(
		roundup(self.gridSize.x, localSize.x),
		roundup(self.gridSize.y, localSize.y),
		roundup(self.gridSize.z, localSize.z))

	-- don't include the ghost cells as a part of the grid coordinate space
	local sizeWithoutBorder = vec3sz(self.gridSize:unpack())
	for i=0,self.dim-1 do
		sizeWithoutBorder.s[i] = sizeWithoutBorder.s[i] - 2 * self.numGhost
	end
	local volumeWithoutBorder = tonumber(sizeWithoutBorder:volume())
	
	local numCells = tonumber(self.gridSize:volume())

	local stepSize = vec3sz()
	stepSize.x = 1
	for i=1,self.dim-1 do
		stepSize.s[i] = stepSize.s[i-1] * self.gridSize.s[i-1]
	end

	local offset = vec3sz(0,0,0)
	
	self.domain = self.app.env:domain{
		size = {self.gridSize:unpack()},
		dim = self.dim,
	}

	self.domainWithoutBorder = self.app.env:domain{
		size = {sizeWithoutBorder:unpack()},
		dim = self.dim,
	}

--[[ TODO use these stored values ... but don't forget they have to be adjusted for each kernel (based on obj.localSize2d)
	local boundaryLocalSize
	local boundaryGlobalSize
	if self.dim == 1 then
		boundaryLocalSize = 1
		boundaryGlobalSize = 1
	elseif self.dim == 2 then
		boundaryLocalSize = math.min(localSize1d, maxWorkGroupSize)
		boundaryGlobalSize = roundup(
				side == 1 
				and tonumber(gridSize.y)
				or tonumber(gridSize.x),
			boundaryLocalSize)
	elseif self.dim == 3 then
		-- xy xz yz
		local maxSizeX = roundup(
			math.max(tonumber(gridSize.x), tonumber(gridSize.y)),
			localSize2d.x)
		local maxSizeY = roundup(
			math.max(tonumber(gridSize.y), tonumber(gridSize.z)),
			localSize2d.y)
		boundaryGlobalSize = {maxSizeX, maxSizeY}
		boundaryLocalSize = {maxSizeX, maxSizeY}
	end	
--]]

	return {
		localSize1d = localSize1d,
		localSize2d = localSize2d,
		localSize = localSize,
		globalSize = globalSize,
		sizeWithoutBorder = sizeWithoutBorder,
		volumeWithoutBorder = volumeWithoutBorder,
		numCells = numCells,
		stepSize = stepSize,
		offset = offset,
--[[
		boundaryGlobalSize = boundaryGlobalSize,
		boundaryLocalSize = boundaryLocalSize,
--]]
	}
end

function GridSolver:initCodeModules()
	GridSolver.super.initCodeModules(self)

--[[
naming conventions ...
* the grid indexes i_1..i_n that span 1 through solver->gridSize.1..solver->gridSize.n
(between the index and the coordinate space:)
	- grid_dx? is the change in coordinate space wrt the change in grid space
* the coordinates space x_1..x_n that spans mins.s1..maxs.sn
(between the coordinate and the embedded space:)
	- vectors can be calculated from Cartesian by cartesianToCoord
	- the length of the basis vectors wrt the change in indexes is given by cell_dx?(x)
	- the Cartesian length of the holonomic basis vectors is given by coord_dx?(x).  
		This is like cell_dx? except not scaled by grid_dx?
		This is just the change in embedded wrt the change in coordinate, not wrt the change in grid
	- cellBuf[index].volume gives the volume between indexes at the coordinate x
	- the Cartesian length of a vector in coordinate space is given by coordLen and coordLenSq
* the embedded Cartesian space ... idk what letters I should use for this.  
	Some literature uses r^i or u^i vs coordinate space x^a.
	Maybe I'll use 'xc' vs 'x', like I've already started to do in the initial conditions.

functionality (and abstraction):
	- Allow the option for precomputing certain variables of the coordinate system: the chart mapping.
	 (Though for some operations, like the Euler equations wavespeeds, it is faster to compute them than to retrieve them from memory, so I suspect similar for the coordinate chart.)
	- Add a macro for useful values to compute, like the x,y,z,r,theta,phi variables.
	 (This is already started in the coords.vars table.)
--]]

	-- if I wanted, I could put 'dim' in 'solver_t'
	-- and then this would be independent of all solvers
	-- and it could go into app
	-- but that would lose some lookup speed of dim
	self.modules:add{
		name = self.symbols.OOB,
		-- bounds-check macro
		headercode = self.eqn:template('#define <?=OOB?>(lhs,rhs) (i.x < (lhs) || i.x >= solver->gridSize.x - (rhs)'
			.. (self.dim < 2 and '' or ' || i.y < (lhs) || i.y >= solver->gridSize.y - (rhs)')
			.. (self.dim < 3 and '' or ' || i.z < (lhs) || i.z >= solver->gridSize.z - (rhs)')
			.. ')'),
	}

	self.modules:add{
		name = self.symbols.SETBOUNDS,
		depends = {
			'INDEXV',
			self.symbols.SETBOUNDS,
			self.symbols.OOB,
		},
		headercode = self.eqn:template[[
// define i, index, and bounds-check
#define <?=SETBOUNDS?>(lhs,rhs)	\
	int4 i = globalInt4(); \
	if (<?=OOB?>(lhs,rhs)) return; \
	int index = INDEXV(i);
]],
	}

	self.modules:add{
		name = self.symbols.SETBOUNDS_NOGHOST,
		depends = {
			self.symbols.OOB,
			'INDEXV',
		},
		headercode = self.eqn:template[[
// same as above, except for kernels that don't use the boundary
// index operates on buffers of 'gridSize' (with border)
// but the kernel must be invoked across sizeWithoutBorder
#define <?=SETBOUNDS_NOGHOST?>() \
	int4 i = globalInt4(); \
	if (<?=OOB?>(0, 2 * solver->numGhost)) return; \
	i += (int4)(]]..range(4):mapi(function(i) 
		return i <= self.dim and 'solver->numGhost' or '0' 
	end):concat','..[[); \
	int index = INDEXV(i);
]],
	}

	if self.usePLM then
		self.modules:add{
			name = self.symbols.slopeLimiter,
			code = self.eqn:template[[
real <?=slopeLimiter?>(real r) {
	<?=solver.app.limiters[solver.slopeLimiter].code?>
}
]],
		}
		
		self.modules:add{
			name = self.eqn.symbols.consLR_t,
			depends = {self.eqn.symbols.cons_t},
			typecode = self.eqn:template([[
typedef union {
	<?=cons_t?> LR[2];
	struct {
		<?=cons_t?> L, R;
	};
} <?=consLR_t?>;

//ugly hack to work around even uglier boundary code
typedef struct {
	<?=consLR_t?> side[<?=solver.dim?>];
} <?=consLR_t?>_dim;
]]),
		}

		self.modules:addFromMarkup(self.eqn:template(file'hydro/solver/plm.cl':read()))
		self.solverModulesEnabled[self.symbols.calcLR] = true
	end

	if self.useCTU then
		self.modules:addFromMarkup(self.eqn:template(file'hydro/solver/ctu.cl':read()))
		self.sharedModulesEnabled[self.symbols.updateCTU] = true
	end
end

-- call this when a gui var changes
-- it rebuilds the code prefix, but doesn't reset the initCond
function GridSolver:refreshCodePrefix()
	GridSolver.super.refreshCodePrefix(self)	-- refresh integrator
	-- changing initCond calls this, and could change boundary programs, so I'm putting this here
	-- bad excuse, I know
	self:refreshBoundaryProgram()
end

function GridSolver:createSolverBuf()
	GridSolver.super.createSolverBuf(self)

	-- do this before any call to createBuffers
	self.solverPtr.gridSize.x = self.gridSize.x
	self.solverPtr.gridSize.y = self.gridSize.y
	self.solverPtr.gridSize.z = self.gridSize.z
	self.solverPtr.gridSize.w = 1
	self.solverPtr.stepsize.x = 1
	self.solverPtr.stepsize.y = self.gridSize.x
	self.solverPtr.stepsize.z = self.gridSize.x * self.gridSize.y
	self.solverPtr.stepsize.w = self.gridSize.x * self.gridSize.y * self.gridSize.z

	-- do I even need separate lua and C structures?
	self:refreshSolverBufMinsMaxs()

	-- while we're here, write all gui vars to the solver_t
	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			if var.ctype == 'real' then
				self.solverPtr[var.name] = toreal(var.value)
			else
				self.solverPtr[var.name] = var.value
			end
		end
	end

	self.solverPtr.numGhost = self.numGhost

	self:refreshSolverBuf()
end

function GridSolver:refreshSolverBufMinsMaxs()
	GridSolver.super.refreshSolverBufMinsMaxs(self)
	-- while here, refersh the 'mindx' param
	self.mindx = math.huge
	-- always set this to the full range, even outside the used dimension, in case some dimensional value is supposed to be a non-zero constant, esp for cell volume calculations
	for j=1,3 do
		local dx = (self.maxs.s[j-1] - self.mins.s[j-1]) / tonumber(self.sizeWithoutBorder.s[j-1])
		self.mindx = math.min(self.mindx, dx)
		self.solverPtr.grid_dx.s[j-1] = toreal(dx)
		self.solverPtr.initCondMins.s[j-1] = toreal( self.initCondMins.s[j-1] )
		self.solverPtr.initCondMaxs.s[j-1] = toreal( self.initCondMaxs.s[j-1] )
	end
	if self.app.verbose then
		print('grid_dx = '..fromreal(self.solverPtr.grid_dx.x)..', '..fromreal(self.solverPtr.grid_dx.y)..', '..fromreal(self.solverPtr.grid_dx.z))
		print('mins = '..fromreal(self.solverPtr.mins.x)..', '..fromreal(self.solverPtr.mins.y)..', '..fromreal(self.solverPtr.mins.z))
		print('maxs = '..fromreal(self.solverPtr.maxs.x)..', '..fromreal(self.solverPtr.maxs.y)..', '..fromreal(self.solverPtr.maxs.z))
		print('initCondMins = '..fromreal(self.solverPtr.initCondMins.x)..', '..fromreal(self.solverPtr.initCondMins.y)..', '..fromreal(self.solverPtr.initCondMins.z))
		print('initCondMaxs = '..fromreal(self.solverPtr.initCondMaxs.x)..', '..fromreal(self.solverPtr.initCondMaxs.y)..', '..fromreal(self.solverPtr.initCondMaxs.z))
	end
end

-- TODO some of this is copied in solverbase
function GridSolver:createBuffers()
	local app = self.app
	
	-- define self.texSize before calling super
	if app.targetSystem ~= 'console' then
		self.texSize = vec3sz(self.gridSize)
		if self.app.verbose then
			print('texSize = '..self.texSize)
		end
	end

	GridSolver.super.createBuffers(self)

	self:clalloc('cellBuf', self.coord.cell_t, self.numCells)

	if self.usePLM then
		-- TODO self.eqn.symbols.consLR_t..'_dim' and remove * self.dim ?
		self:clalloc('ULRBuf', self.eqn.symbols.consLR_t, self.numCells * self.dim)
	end
end

function GridSolver:refreshSolverProgram()
	GridSolver.super.refreshSolverProgram(self)

	if self.usePLM then
		self.calcLRKernelObj = self.solverProgramObj:kernel(self.symbols.calcLR)
	end
	if self.useCTU then
		-- currently implemented in hydro/solver/roe.cl
		-- not available for any other flux method
		assert(self.fluxBuf)
		self.updateCTUKernelObj = self.solverProgramObj:kernel(self.symbols.updateCTU)
	end
end

function GridSolver:getULRBuf()
	return self[self.getULRBufName]
end

-- apply boundary to the ULRBuf ... or UBuf if it is not a PLM solver
function GridSolver:boundaryLR()
	if self.usePLM then
		for _,obj in ipairs(self.lrBoundaryKernelObjs) do
			obj()
			if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		end
	else
		self:boundary()
	end
end

function GridSolver:refreshGetULR()
	-- set pointer to the buffer holding the LR state information
	-- for piecewise-constant that is the original UBuf
	-- for piecewise-linear that is the ULRBuf
	-- TODO do this after createBuffers and you can get rid of a lot of these
	self.getULRBufType = self.usePLM and self.eqn.symbols.consLR_t or self.eqn.symbols.cons_t
	self.getULRBufName = self.usePLM and 'ULRBuf' or 'UBuf'

	self.getULRArg = self.getULRBufType..'* '..self.getULRBufName

	-- this code creates the 'global cons_t const * const UL, Ul'R variables
	-- it assumes that indexL, indexR, and side are already defined
	-- both UBuf and ULRBuf should be cell-centered
	if self.usePLM then
		self.getULRCode = function(self, args)
			args = args or {}
			local suffix = args.suffix or ''
			return self.eqn:template([[
global <?=cons_t?> const * const UL<?=suffix?> = &<?=bufName?>[<?=side?> + dim * <?=indexL?>].R;
global <?=cons_t?> const * const UR<?=suffix?> = &<?=bufName?>[<?=side?> + dim * <?=indexR?>].L;
]],			{
				suffix = suffix,
				side = args.side or 'side',
				indexL = args.indexL or 'indexL'..suffix,
				indexR = args.indexR or 'indexR'..suffix,
				bufName = args.bufName or self.getULRBufName,	-- for displayVars the variable name is 'buf', so I need to override it either in displayCode or here
			})
		end
	else 
		self.getULRCode = function(self, args)
			args = args or {}
			local suffix = args.suffix or ''
			return self.eqn:template([[
global <?=cons_t?> const * const UL<?=suffix?> = <?=bufName?> + <?=indexL?>;
global <?=cons_t?> const * const UR<?=suffix?> = <?=bufName?> + <?=indexR?>;
]],			{
				suffix = suffix,
				indexL = args.indexL or 'indexL'..suffix,
				indexR = args.indexR or 'indexR'..suffix,
				bufName = args.bufName or self.getULRBufName,	-- for displayVars the variable name is 'buf', so I need to override it either in displayCode or here
			})
		end
	end
end

function GridSolver:resetState()
	self.cmds:finish()
		
	-- start off by filling all buffers with zero, just in case initCond forgets something ...
	for _,bufferInfo in ipairs(self.buffers) do
		self.cmds:enqueueFillBuffer{buffer=self[bufferInfo.name], size=bufferInfo.count * ffi.sizeof(bufferInfo.type)}
	end

	GridSolver.super.resetState(self)
end

function GridSolver:applyInitCond()
	
	-- this is a bit disorganized
	-- this fills the cellBuf but only for gridsolvers
	-- and does so by asking the coord obj
	-- because coord objects can specify arbitrary fields
	-- (such as x,y,z, r, remapped-r, etc)
	-- TODO how about making this mimic meshsolver?
	-- and then make the only difference be finite difference operators
	self.cellCpuBuf = ffi.new(self.coord.cell_t..'[?]', self.numCells)
	self.coord:fillGridCellBuf(self.cellCpuBuf)
	self.cellBufObj:fromCPU(self.cellCpuBuf)

	GridSolver.super.applyInitCond(self)
end

-------------------------------------------------------------------------------
--                                 display vars                              --
-------------------------------------------------------------------------------


-- subclass and override
local DisplayVar = class(GridSolver.DisplayVar)
GridSolver.DisplayVar = DisplayVar


-- and this is the DisplayVar used for UBuf
GridSolver.DisplayVar_U = DisplayVar


-------------------------------------------------------------------------------
--                              boundary                                     --
-------------------------------------------------------------------------------
		
-- TODO constant/'dirichlet' conditions
-- and while you're at it, let boundary options register GUI controls, so we can treat the constant as a GUI-editable parameter 

local Boundary = class()
GridSolver.Boundary = Boundary

-- reflection, but only when the normal is coordinate-aligned
function Boundary:reflectVars(args, dst, varnames, restitution)
	restitution = restitution or 1
	local lines = table()
	if varnames and varnames[args.side] then
		for _,varname in ipairs(varnames[args.side]) do
			if not args.fields or table.find(args.fields, varname) then
				lines:insert('buf['..dst..'].'..varname..' = '..clnumber(-restitution)..' * buf['..dst..'].'..varname..';')
			end
		end
	end
	return lines:concat'\n'
end

function Boundary:assignDstSrc(dst, src, args)
	local lines = table()
	if args.fields then
		for _,field in ipairs(args.fields) do
			lines:insert('buf['..dst..'].'..field..' = buf['..src..'].'..field..';')
		end
	else
		lines:insert('buf['..dst..'] = buf['..src..'];')
	end
	return lines:concat'\n'
end

--[[
args of the BoundaryMethod:getCode are:
	index
	fields = table of what fields to perform assignment to 
	field
	array
	gridSizeSide = solver.gridSize[side]
	side = 1,2,3
	minmax = 'min' or 'max'
	reflectVars
		= which vars to reflect.
			for solver.UBuf this is taken from eqn
			for poisson this is nothing

This is such a mess.  It's practically an AST.
--]]
function Boundary:getCode(args)
	error("implement me")
end

-- this is purely for debugging.  annnd it doesn't crash, whereas rectangular grids are crashing...
local BoundaryNone = class(Boundary)
BoundaryNone.name = 'none'
function BoundaryNone:getCode(args)
	return ''
end

-- boundary code indentation
local indent = '\t\t'

--[[
TODO how to handle individual fields
for example, this can't compile with op/selfgrav.lua, because that will try to assign .ePot onto the struct, which is already being dereferenced here 

same with mirror, only certain fields are inverted

so how do we specify when the boundary is being applied to the whole structure or only to a single field in it?

add a new field: 'fields' that says what fields to apply the boundary to
--]]

-- aka torus
local BoundaryPeriodic = class(Boundary)
BoundaryPeriodic.name = 'periodic'
function BoundaryPeriodic:getCode(args)
	local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
	local dst, src
	if args.minmax == 'min' then
		dst = args.index'j'
		src = args.index([[solver->numGhost + (j-solver->numGhost + 2 * (]]..gridSizeSide..' - 2 * solver->numGhost)) % ('..gridSizeSide..[[ - 2 * solver->numGhost)]])
	elseif args.minmax == 'max' then
		dst = args.index(gridSizeSide..'-1-j')
		src = args.index('solver->numGhost + (solver->numGhost-1-j) % ('..gridSizeSide..' - 2 * solver->numGhost)')
	end
	return self:assignDstSrc(dst, src, args)
end

-- TODO incorporate surface normal and restitution
local BoundaryMirror = class(Boundary)
BoundaryMirror.name = 'mirror'
BoundaryMirror.restitution = 1
function BoundaryMirror:init(args)
	if args then
		self.restitution = args.restitution
	end
end
function BoundaryMirror:getCode(args)
	local dst, src
	if args.minmax == 'min' then
		dst = args.index'j'
		src = args.index'2 * solver->numGhost - 1 - j'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..' - solver->numGhost + j')
		src = args.index(gridSizeSide..' - solver->numGhost - 1 - j')
	end
	local lines = table()
	lines:insert(self:assignDstSrc(dst, src, args))
	
	-- reflectVars.mirror is going to hold, for each dimension, which components need to be mirrored
	-- v.s[side-1] = -v.s[side-1]
	-- generalized:
	local solver = args.solver
	local eqn = solver.eqn
	if solver.coord.vectorComponent == 'cartesian' 
	and not require 'hydro.coord.cartesian':isa(solver.coord)
	then
		-- v = v - n (v dot n)
		-- v^i = v^i - n^i (v^j n_j) (1 + restitution)
		-- ... where n_i = partial_i u, for u the chart
		-- TODO let the ctor specify the vars, not this
		--  so I can use this for things like the poisson solver
		lines:insert(template([[
{
	real3 const x = cellBuf[INDEX(<?=iv?>)].pos;
<? if args.minmax == 'min' then ?>
	real3 const n = coord_cartesianFromCoord(normalForSide<?=side-1?>, x);
<? else -- max ?>
	real3 const n = coord_cartesianFromCoord(real3_neg(normalForSide<?=side-1?>), x);
<? end ?>
]], 	{
			args = args,
			iv = args.minmax == 'min' 
				and args.indexv'j'
				or args.indexv('solver->gridSize.'..xNames[args.side]..' - solver->numGhost + j'),
			side = args.side,
		}))
		for _,var in ipairs(eqn.consStruct.vars) do
			if var.type == 'real' 
			or var.type == 'cplx'
			then
				-- do nothing
			elseif var.type == 'real3' 
			or var.type == 'cplx3'
			then
				-- TODO looks very similar to the reflect code in meshsolver
				lines:insert(template([[
	<?=result?>-><?=field?> = <?=vec3?>_sub(
		<?=result?>-><?=field?>,
		<?=vec3?>_<?=scalar?>_mul(
			<?=vec3?>_from_real3(n),
			<?=scalar?>_real_mul(
				<?=vec3?>_real3_dot(
					<?=result?>-><?=field?>,
					n
				), 
				<?=restitutionPlusOne?>
			)
		)
	);
]], 			{
					restitutionPlusOne = clnumber(self.restitution + 1),
					vec3 = var.type,
					scalar = var.type == 'cplx3' and 'cplx' or 'real',
					field = var.name,
					result = '(buf + '..dst..')',
				}))
			else
				error("need to support reflect() for type "..var.type)
			end
		end
		lines:insert[[
}
]]
	else
		lines:insert(self:reflectVars(args, dst, args.reflectVars.mirror, self.restitution))
	end
	return lines:concat'\n'
end

-- Dirichlet boundary conditions: constant values
-- fixedCode = function(boundary, args, dst) 
--	returns code, list of dependency module names
local BoundaryFixed = class(Boundary)
BoundaryFixed.name = 'fixed'
function BoundaryFixed:init(args)
	-- fixed values to use
	self.fixedCode = assert((args or {}).fixedCode, "you didn't provide any fixedCode arg for boundary==fixed")
end
function BoundaryFixed:getCode(args)
	local dst
	if args.minmax == 'min' then
		dst = args.index'j'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..' - solver->numGhost + j')
	end
	return self:fixedCode(args, dst)
end
GridSolver.BoundaryFixed = BoundaryFixed 

local BoundaryFreeFlow = class(Boundary)
BoundaryFreeFlow.name = 'freeflow'
function BoundaryFreeFlow:getCode(args)
	local dst, src
	if args.minmax == 'min' then
		dst = args.index'j'
		src = args.index'solver->numGhost'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..' - solver->numGhost + j')
		src = args.index(gridSizeSide..' - solver->numGhost - 1')
	end
	return self:assignDstSrc(dst, src, args)
end

-- linear extrapolation
local BoundaryLinear = class(Boundary)
BoundaryLinear.name = 'linear'
function BoundaryLinear:getCode(args)
	local dst, i1, i2
	if args.minmax == 'min' then
		dst = args.index'solver->numGhost - j - 1'
		i1 = args.index'solver->numGhost - j'
		i2 = args.index'solver->numGhost - j + 1'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..' - solver->numGhost + j')
		i1 = args.index(gridSizeSide..' - solver->numGhost + j - 1')
		i2 = args.index(gridSizeSide..' - solver->numGhost + j - 2')
	end
	local lines = table()
	local function addField(field)
		lines:insert('buf['..dst..'].'..field..' = 2. * buf['..i1..'].'..field..' - buf['..i2..'].'..field..';')
	end
	if args.fields then
		for _,field in ipairs(args.fields) do
			addField(field)
		end
	else
		lines:insert('for (int k = 0; k < numStates; ++k) {')
		addField'ptr[k]'
		lines:insert('}')
	end
	return lines:concat'\n'
end

-- quadratic extrapolation
-- f[i] = 3*f[i+1] - 3*f[i+2] + f[i+3]
BoundaryQuadratic = class(Boundary)
BoundaryQuadratic.name = 'quadratic'
function BoundaryQuadratic:getCode(args)
	-- j is initialized from 0 to numGhost-1
	local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
	local dst, i1, i2, i3
	if args.minmax == 'min' then
		dst = args.index'solver->numGhost - j - 1'
		i1 = args.index'solver->numGhost - j'
		i2 = args.index'solver->numGhost - j + 1'
		i3 = args.index'solver->numGhost - j + 2'
	elseif args.minmax == 'max' then
		dst = args.index(gridSizeSide..' - solver->numGhost + j')
		i1 = args.index(gridSizeSide..' - solver->numGhost + j - 1')
		i2 = args.index(gridSizeSide..' - solver->numGhost + j - 2')
		i3 = args.index(gridSizeSide..' - solver->numGhost + j - 3')
	end
	local lines = table()
	local function addField(field)
		lines:insert('buf['..dst..'].'..field..' = 3.*buf['..i1..'].'..field..' - 3.*buf['..i2..'].'..field..' + buf['..i3..'].'..field..';')
	end
	if args.fields then
		for _,field in ipairs(args.fields) do
			addField(field)
		end
	else
		lines:insert'for (int k = 0; k < numStates; ++k) {'
		addField'ptr[k]'
		lines:insert'}'
	end
	return lines:concat'\n'
end

--[[
2013 Baumgarte et al, "Numerical Relativity in Spherical Polar Coordinates..."

same as mirror, except we only negative certain components.
in fact, mirror is very cartesian-specific...

spherical r=0 ...
for 1D, for cell[i] we use cell[i] and mirror the respective coordinates
TODO:
for 2D (r,theta), for theta in (0,pi) for cell[i][j], we use cell[i][size.y-j]
for 3D (r,theta,phi), for theta in [0,pi), phi in [0,2pi) for cell[i][j][k] we use cell cell[i][size.y-j][(k+size.z/2)%size.z]
	This requires knowledge of the domain of the grid, and the boundary condition of the phi variable.
	For now I'm just going to assume that the phi domain is [0,2pi)+c and has periodic boundary condition.
--]]
local BoundarySphereRMin = class(Boundary)
BoundarySphereRMin.name = 'sphereRMin'
function BoundarySphereRMin:getCode(args)
	local solver = args.solver
	
	assert(args.side == 1 and args.minmax == 'min', "you should only use this boundary condition for rmin with spherical coordinates")
	assert(require 'hydro.coord.sphere':isa(solver.coord)
		or require 'hydro.coord.sphere_sinh_radial':isa(solver.coord), "you should only use this boundary condition for rmin with spherical coordinates")
	--assert(solver.maxs.y - solver.mins.y == 2*math.pi)
	--assert(solver.boundaryMethods.ymin == 'periodic' and solver.boundaryMethods.ymax == 'periodic')

	local src, dst
	if solver.dim == 1 then
		dst = 'INDEX(j, 0, 0)'
		src = 'INDEX(2 * solver->numGhost - 1 - j, 0, 0)'
	elseif solver.dim == 2 then	-- r, theta
		dst = 'INDEX(j, i, 0)'
		src = 'INDEX(2 * solver->numGhost - 1 - j, solver->gridSize.y - i - 1, 0)'
	elseif solver.dim == 3 then
		dst = 'INDEX(j, i.x, i.y)'
		src = [[
	INDEX(
		2 * solver->numGhost - 1 - j,
		solver->gridSize.y - i.x - 1, 
		(i.y - solver->numGhost + (solver->gridSize.z - 2 * solver->numGhost) / 2
			+ (solver->gridSize.z - 2 * solver->numGhost)) 
			% (solver->gridSize.z - 2 * solver->numGhost) + solver->numGhost
	)
]]
	end
	local lines = table()
	lines:insert(self:assignDstSrc(dst, src, args))
	lines:insert(self:reflectVars(args, dst, args.reflectVars.sphereRMin))
	return lines:concat'\n'
end

--[[
Like above, this is going to assume certain things about the coordinate domain.
I'll assume theta is [0,pi) and phi is [-pi,pi]
1D: nothing special
2D: (r, theta): 
3D: (r, theta, phi): 
--]]
local BoundarySphereTheta = class(Boundary)
BoundarySphereTheta.name = 'sphereTheta'
function BoundarySphereTheta:getCode(args)
	local solver = args.solver

	assert(args.side == 2, "you should only use this boundary condition for θmin/θmax with spherical coordinates")
	assert(require 'hydro.coord.sphere':isa(solver.coord)
		or require 'hydro.coord.sphere_sinh_radial':isa(solver.coord), "you should only use this boundary condition for θmin/θmax with spherical coordinates")

	local src, dst
	if args.minmax == 'min' then
		if solver.dim == 1 then
			dst = args.index'j'
			src = args.index'2 * solver->numGhost - 1 - j'
		elseif solver.dim == 2 then
			dst = args.index'j'
			src = args.index'2 * solver->numGhost - 1 - j'
		elseif solver.dim == 3 then
			dst = 'INDEX(i.x, solver->numGhost - 1 - j, i.y)'
			src = [[
	INDEX(
		i.x,
		solver->numGhost + j,
		(i.y - solver->numGhost + (solver->gridSize.z - 2 * solver->numGhost) / 2 
			+ (solver->gridSize.z - 2 * solver->numGhost))
			% (solver->gridSize.z - 2 * solver->numGhost) + solver->numGhost
	)
]]
		end
	elseif args.minmax == 'max' then
		if solver.dim == 1 then
			dst = args.index'solver->gridSize.y - 1 - j'
			src = args.index'solver->gridSize.y - 2 * solver->numGhost + j'
		elseif solver.dim == 2 then
			dst = args.index'solver->gridSize.y - 1 - j'
			src = args.index'solver->gridSize.y - 2 * solver->numGhost + j'
		elseif solver.dim == 3 then
			dst = 'INDEX(i.x, solver->gridSize.y - solver->numGhost + j, i.y)'
			src = [[
	INDEX(
		i.x,
		solver->gridSize.y - solver->numGhost - 1 - j,
		(i.y - solver->numGhost + (solver->gridSize.z - 2 * solver->numGhost) / 2
			+ (solver->gridSize.z - 2 * solver->numGhost))
			% (solver->gridSize.z - 2 * solver->numGhost) + solver->numGhost
	)
]]	
		end
	end
	local lines = table()
	lines:insert(self:assignDstSrc(dst, src, args))
	lines:insert(self:reflectVars(args, dst, args.reflectVars.sphereTheta))
	return lines:concat'\n'
end


local BoundaryCylinderRMin = class(Boundary)
BoundaryCylinderRMin.name = 'cylinderRMin'
function BoundaryCylinderRMin:getCode(args)
	local solver = args.solver
	
	assert(args.side == 1 and args.minmax == 'min', "you should only use this boundary condition for rmin with cylinderical coordinates")
	assert(require 'hydro.coord.cylinder':isa(solver.coord), "you should only use this boundary condition for rmin with cylinderical coordinates")
	--assert(solver.maxs.y - solver.mins.y == 2*math.pi)
	--assert(solver.boundaryMethods.ymin == 'periodic' and solver.boundaryMethods.ymax == 'periodic')

	local src, dst
	if solver.dim == 1 then
		dst = 'INDEX(j, 0, 0)'
		src = 'INDEX(2 * solver->numGhost - 1 - j, 0, 0)'
	elseif solver.dim == 2 then	-- r, theta
		dst = 'INDEX(j, i, 0)'
		src = [[
	INDEX(
		2 * solver->numGhost - 1 - j, 
		(i - solver->numGhost + (solver->gridSize.y - 2 * solver->numGhost) / 2
			+ (solver->gridSize.y - 2 * solver->numGhost))
			% (solver->gridSize.y - 2 * solver->numGhost) + solver->numGhost,
		0)
]]
	elseif solver.dim == 3 then
		dst = 'INDEX(j, i.x, i.y)'
		src = [[
	INDEX(
		2 * solver->numGhost - 1 - j, 
		(i.x - solver->numGhost + (solver->gridSize.y - 2 * solver->numGhost) / 2
			+ (solver->gridSize.y - 2 * solver->numGhost))
			% (solver->gridSize.y - 2 * solver->numGhost) + solver->numGhost,
		i.y)
]]
	end
	local lines = table()
	lines:insert(self:assignDstSrc(dst, src, args))
	lines:insert(self:reflectVars(args, dst, args.reflectVars.cylinderRMin))
	return lines:concat'\n'
end

--[[
-- boundaryOptions is a table of classes
--]]
function GridSolver:createBoundaryOptions()
	self:addBoundaryOptions{
		BoundaryNone,
		BoundaryPeriodic,
		BoundaryMirror,
		BoundaryFixed,
		BoundaryFreeFlow,
		BoundaryLinear,
		BoundaryQuadratic,
		-- specific to coordinate chart grids, where vector components flip as you cross coordinate symmetry boundaries
		BoundarySphereRMin,
		BoundarySphereTheta,
		BoundaryCylinderRMin,
	}
end

--[[
this is shorthand for assignment
you can also assign each solver.boundaryMethod.[x|y|z][min|max] to an index in solver.boundaryOptions
args:
	a string = sets all boundaryMethods to the boundaryOptions index with name of the string in args
	a table = sets each xmin...zmax boundaryMethod with the associated name
--]]
function GridSolver:setBoundaryMethods(args)
	for _,x in ipairs(xNames) do
		for _,minmax in ipairs(minmaxs) do
			local k = x..minmax
			local name
			if type(args) == 'string' then
				-- set all sides to the same type of BC
				name = args	
			elseif type(args) == 'table' then
				-- set each side separately
				name = args[k]
			else
				error("don't know how to deal with these args") 
			end
			local boundaryObjArgs
			if type(name) == 'table' then
				boundaryObjArgs = name.args
				name = name.name
			end
			local boundaryClass = self.boundaryOptionForName[name]
			if not boundaryClass then
				error("unable to find boundary method "..tostring(name)
					..' ('..type(name)..')'
					.." -- can't assign it to side "..k.."\n")
			else
				self.boundaryMethods[k] = boundaryClass(boundaryObjArgs)
			end
		end
	end
end

--[[
args:
	type = type of buffer to create boundary information
	extraArgs = any extra args to add to the kernel (other than a global of the buffer type)
	methods = {
		[x|y|z..min|max] = a Boundary subclass who has function getCode(args)
			that returns the code to apply to that particular pointer 
	}
	assign = function(a,b) a = b 
	field = function(a,b) a.b
--]]
function GridSolver:createBoundaryProgramAndKernel(args)
	local field = args.field or function(a, b) return a .. '.' .. b end

	local lines = table()

	local moduleNames = self.sharedModulesEnabled:keys():append{
		self.solver_t,
		self.eqn.symbols.cons_t,
		assert(self.coord.cell_t),
		'INDEX',
		'INDEXV',
		-- some Boundary :getCode use numStates
		-- TODO use the addCodeMarkup function and inline these all?
		self.symbols.solver_macros,
		self.coord.symbols.cartesianFromCoord,
		'normalForSide',
	}

	local iFields = ({
		nil,
		{''},	-- i
		{'.x', '.y'},	-- i.x, i.y
	})[self.dim]

	for side=1,self.dim do
		
		local bxs = range(3)
		bxs:remove(side)

		local function indexv(j)
			local v = table{'0','0','0'}
			v[side] = j
			
			for k=1,self.dim-1 do
				v[bxs[k]] = 'i'..iFields[k]
			end
			
			return v:concat','
		end

		local function index(j)
			return 'INDEX('..indexv(j)..')'
		end

		lines:insert(self.eqn:template([[
kernel void boundary_<?=xNames[side]?>(
	constant <?=solver_t?> const * const solver,
	global <?=args.type?>* buf,
	global <?=cell_t?> const * const cellBuf<?= 
args.extraArgs and #args.extraArgs > 0 
	and ',\n\t'..table.concat(args.extraArgs, ',\n\t')
	or '' ?>
) {<?
-- 1D: use a small 1D kernel and just run once 
-- 2D: use a 1D kernel the size of the max dim 
if solver.dim == 2 then ?>
	int i = get_global_id(0);<? 
elseif solver.dim == 3 then ?>
	int2 i = (int2)(get_global_id(0), get_global_id(1));<? 
end 
?>]], 	{
			table = table,
			args = args,
			side = side, 
		}))

		if self.dim > 1 then
			lines:insert(template([=[
	if (<?
local sep = ''
for k=1,solver.dim-1 do
?><?=sep?>i<?=iFields[k]?> >= solver->gridSize.<?=xNames[bxs[k]]?><?	
	sep = ' || '
end
?>) return;]=],
			{
				bxs = bxs,
				solver = self,
				xNames = xNames,
				iFields = iFields,
			}))
		end

		lines:insert[[
	for (int j = 0; j < solver->numGhost; ++j) {]]

		for _,minmax in ipairs(minmaxs) do
			local sideCode, sideDepends = args.methods[xNames[side]..minmax]:getCode({
				
				-- this is only used for the oscillating boundary to distinguish between the cons_t boundary and the potential boundary
				fields = args.fields,
				
				solver = self,
				index = index,
				indexv = indexv,
				field = field,
				side = side,
				minmax = minmax,
				reflectVars = args.reflectVars -- UBuf takes reflectVars from eqn
					or {},						-- op/relaxation doesn't use reflectVars (should it?)
			})

			moduleNames:append(sideDepends)
			
			-- add tab
			sideCode = '\t\t'..sideCode:gsub('\n', '\n\t\t')
			-- and remove tabs from MODULE_* lines
			sideCode = string.split(sideCode, '\n'):mapi(function(l)
				return (l:gsub('^%s*//// MODULE_', '//// MODULE_'))
			end):concat'\n'	

			lines:insert(sideCode)
		end 

lines:insert[[
	}
}
]]
	end

	-- last, add dependency code to the beginning
	lines:insert(1, self.modules:getCodeAndHeader(moduleNames:unpack()))

	local code = lines:concat'\n'
	local boundaryProgramName = 'boundary'..(args.programNameSuffix or '')

	local boundaryProgramObj
	time('building program cache/'..self:getIdent()..'/src/boundary.cl ', function()
		boundaryProgramObj = self.Program{name=boundaryProgramName, code=code}
		boundaryProgramObj:compile()
	end)
	
	local boundaryKernelObjs = table()
	for i=1,self.dim do
		local kernelObj = boundaryProgramObj:kernel('boundary_'..xNames[i])
		boundaryKernelObjs:insert(kernelObj)
		kernelObj.obj:setArg(0, self.solverBuf)
	end

	-- TODO switch these over to obj
	return boundaryProgramObj, boundaryKernelObjs
end

function GridSolver:getBoundaryProgramArgs()
	return {
		type = self.eqn.symbols.cons_t,
		-- remap from enum/combobox int values to functions from the solver.boundaryOptions table
		methods = self.boundaryMethods,

		-- Generated by eqn depending on eqn's vars and the coord sys.
		-- These should be here so that the LR boundary code gen can modify them (duplicate them for L and R).
		-- But it needs to be passed on to BoundaryMethod:getCode(), which means forwarding it again inside of :createBoundaryProgramAndKernel.
		-- Note to self, this is one-to-one with Boundary names, *not* with coordinate systems.
		reflectVars = self.eqn.reflectVars,
	}
end

function GridSolver:refreshBoundaryProgram()
	self.boundaryProgramObj, self.boundaryKernelObjs = self:createBoundaryProgramAndKernel(self:getBoundaryProgramArgs())
	for _,obj in ipairs(self.boundaryKernelObjs) do
		obj.obj:setArg(1, self.UBuf)
		obj.obj:setArg(2, self.cellBuf)
	end
	for _,op in ipairs(self.ops) do
		if op.refreshBoundaryProgram then
			op:refreshBoundaryProgram()
		end
	end

	if self.useCTU and self.usePLM then
		local args = self:getBoundaryProgramArgs()
		args.type = self.eqn.symbols.consLR_t..'_dim'
	
		-- TODO use the real list of boundary names here
		for _,boundaryOptionName in ipairs(self.boundaryOptionNames) do
			-- TODO do this for sphere/cylinderical mirror vars as well
			if args.reflectVars[boundaryOptionName] then
				local newvars = {}
				for i=1,3 do
					newvars[i] = table()
					for _,var in ipairs(args.reflectVars[boundaryOptionName][i]) do
						for j=0,self.dim-1 do
							newvars[i]:insert('side['..j..'].L.'..var)
							newvars[i]:insert('side['..j..'].R.'..var)
						end
					end
				end
				args.reflectVars[boundaryOptionName] = newvars
			end
		end
		args.programNameSuffix = '_lr'
		self.lrBoundaryProgramObj, self.lrBoundaryKernelObjs = self:createBoundaryProgramAndKernel(args)
		for _,obj in ipairs(self.lrBoundaryKernelObjs) do
			obj.obj:setArg(1, self:getULRBuf())
			obj.obj:setArg(2, self.cellBuf)
		end
	end
end

-- assumes the buffer is already in the kernel's arg
function GridSolver:applyBoundaryToBuffer(kernelObjs)
	for side,obj in ipairs(kernelObjs) do	
		-- 1D:
		if self.dim == 1 then
			-- if you do change this size from anything but 1, make sure to add conditions to the boundary kernel code
			self.cmds:enqueueNDRangeKernel{
				kernel = obj.obj,
				globalSize = 1,
				localSize = 1,
			}
		elseif self.dim == 2 then
			local localSize = math.min(self.localSize1d, self.maxWorkGroupSize)
			local maxSize = roundup(
					side == 1 
					and tonumber(self.gridSize.y)
					or tonumber(self.gridSize.x),
				localSize)
			self.cmds:enqueueNDRangeKernel{
				kernel = obj.obj,
				globalSize = maxSize,
				localSize = localSize,
			}
		elseif self.dim == 3 then
			-- xy xz yz
			local maxSizeX = roundup(
				math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y)),
				tonumber(self.localSize2d.x))
			local maxSizeY = roundup(
				math.max(tonumber(self.gridSize.y), tonumber(self.gridSize.z)),
				tonumber(self.localSize2d.y))
			self.cmds:enqueueNDRangeKernel{
				kernel = obj.obj,
				globalSize = {maxSizeX, maxSizeY},
				localSize = {
					math.min(obj.localSize2d[1], maxSizeX),
					math.min(obj.localSize2d[2], maxSizeY),
				},
			}
		else
			error("can't run boundary for dim "..tonumber(self.dim))
		end
	end
end

function GridSolver:boundary()
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	self:applyBoundaryToBuffer(self.boundaryKernelObjs)
	if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
end


-------------------------------------------------------------------------------
--                                                                           --
-------------------------------------------------------------------------------


local function compareL1(ptr, n, ...)
	local err = 0
	n = n or select('#', ...)
	for i=1,n do
		err = err + math.abs(tonumber(ptr[i-1]) - select(i, ...))
	end
	return err
end

function GridSolver:calcExactError(numStates)
	local eqn = self.eqn
	local initCond = eqn.initCond
	numStates = numStates or eqn.numIntStates
	assert(initCond.exactSolution, "can't test accuracy of a configuration that has no exact solution")
	-- TODO doesn't this risk some mem leaks?  does cast still lose refcounts?
	-- looks like so far no one is using this so I wouldn't know
	local UCpuBuf = ffi.cast(eqn.symbols.cons_t..'*', self.UBufObj:toCPU())
	local ghost = self.numGhost
	local imin, imax = ghost,tonumber(self.gridSize.x)-2*ghost-1
	local jmin, jmax = ghost,tonumber(self.gridSize.y)-2*ghost-1
	if self.dim < 2 then jmin, jmax = 0, 0 end
	local kmin, kmax = ghost,tonumber(self.gridSize.z)-2*ghost-1
	if self.dim < 3 then kmin, kmax = 0, 0 end
	assert(self.cellCpuBuf)
	self.cellBufObj:fromCPU(self.cellCpuBuf)
	local err = 0
	for i=imin,imax do
		for j=jmin,jmax do
			for k=kmin,kmax do
				local index = i + self.solverPtr.stepsize.x * j + self.solverPtr.stepsize.y * k
				local cell = self.cellCpuBuf[index]
				local U = UCpuBuf[index]
				err = err + compareL1(U.ptr, numStates, initCond:exactSolution(self.t, cell.pos:unpack()))
			end
		end
	end
	err = err / (numStates * tonumber(self.sizeWithoutBorder:volume()))
	return err, UCpuBuf
end

function GridSolver:updateGUIParams()	
	GridSolver.super.updateGUIParams(self)

	ig.igText('grid size: '
		..tonumber(self.sizeWithoutBorder.x)..', '
		..tonumber(self.sizeWithoutBorder.y)..', '
		..tonumber(self.sizeWithoutBorder.z))

	-- let the user change even unused dimensions
	-- because those influence what constant value the cell_x() gives us 
	for i=1,3 do
		for j,minmax in ipairs(minmaxs) do
			local k = xNames[i]..minmax
			if ig.luatableTooltipInputFloatAsText(k, self[minmax..'s'], xNames[i], ig.ImGuiInputTextFlags_EnterReturnsTrue) then
				local eps = 1e-7
				if self.maxs.s[i-1] - self.mins.s[i-1] < eps then
					self.maxs.s[i-1] = self.mins.s[i-1] + eps
				end
				if j==1 then
					self.solverPtr.mins.s[i-1] = toreal(self.mins.s[i-1])
				elseif j==2 then
					self.solverPtr.maxs.s[i-1] = toreal(self.maxs.s[i-1])
				end
				self:refreshSolverBufMinsMaxs()
				self:refreshSolverBuf()
			end
		end
	end
	
	for i=1,self.dim do
		for _,minmax in ipairs(minmaxs) do
			local var = xNames[i]..minmax
			--[[ TODO is this crashing too?  maybe comboTable has a problem
			if ig.luatableTooltipCombo(var, self.boundaryMethods, var, self.boundaryOptionNames) then
				self:refreshBoundaryProgram()
			end
			--]]
			ig.igText(var..': '..self.boundaryMethods[var].name)
		end
	end
end

function GridSolver:saveBuffer(buffer, basefn)
	buffer = buffer.obj or buffer
	local Image = require 'image'
	
	-- TODO add planes to image, then have the FITS module use planes and not channels
	-- so the dimension layout of the buffer is [channels][width][height][planes]
	local width = tonumber(self.gridSize.x)
	local height = tonumber(self.gridSize.y)
	local depth = tonumber(self.gridSize.z)
	
	local channels = buffer.size / self.numCells / ffi.sizeof(self.app.real)
	if channels ~= math.floor(channels) then
		print("can't save "..basefn.." buffer due to its size not being divisible by the numCells")
	else

		local numReals = self.numCells * channels
--[[ 3D interleave the depth and the channels ... causes FV to break 
		local image = Image(width, height, depth * channels, assert(self.app.real))
		local function getIndex(ch,i,j,k) return i + width * (j + height * (ch + channels * k)) end
--]]
-- [[ 3D interleave the depth and the width
		local image = Image(width * depth, height, channels, assert(self.app.real))
		local function getIndex(ch,i,j,k) return i + width * (k + depth * (j + height * ch)) end
--]]
		self.cmds:enqueueReadBuffer{buffer=buffer, block=true, size=ffi.sizeof(self.app.real) * numReals, ptr=image.buffer}
		local src = image.buffer
		
		-- now convert from interleaved to planar
		-- *OR* add planes to the FITS output
		local tmp = ffi.new(self.app.real..'[?]', numReals)
		for ch=0,channels-1 do
			for k=0,depth-1 do
				for j=0,height-1 do
					for i=0,width-1 do
						tmp[getIndex(ch,i,j,k)] = src[ch + channels * (i + width * (j + height * k))]
					end
				end
			end
		end
		image.buffer = tmp
	
		local filename = basefn..'.fits'
--print('saving '..filename)
		image:save(filename)
	end
end

-- global
local debugSaveBufferIndex = 0
function GridSolver:debugSaveBuffer(buffer, msg)
	debugSaveBufferIndex = debugSaveBufferIndex + 1  
	msg = table{debugSaveBufferIndex, msg}:concat'_' 
	self:saveBuffer(buffer, msg) 
end

function GridSolver:save(prefix)
	for _,bufferInfo in ipairs(self.buffers) do
		local name = bufferInfo.name
		self:saveBuffer(self[name], (prefix and prefix..'_' or '')..name)
	end
end

return GridSolver
