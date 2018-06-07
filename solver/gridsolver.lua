local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local math = require 'ext.math'
local vec3 = require 'vec.vec3'
local glreport = require 'gl.report'
local clnumber = require 'cl.obj.number'
local template = require 'template'
local vec3sz = require 'ffi.vec.vec3sz'
local tooltip = require 'tooltip'
local roundup = require 'roundup'
local time, getTime = table.unpack(require 'time')
local SolverBase = require 'solver.solverbase'

local common = require 'common'()
local minmaxs = common.minmaxs
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


--local tryingAMR = 'dt vs 2dt'
--local tryingAMR = 'gradient'


local GridSolver = class(SolverBase)

GridSolver.numGhost = 2

--[[
args:
	gridSize
	mins
	maxs
	boundary = boundary info
--]]
function GridSolver:init(args)
	GridSolver.super.init(self, args)

	self.mins = vec3(table.unpack(args.mins or {-1, -1, -1}))
	self.maxs = vec3(table.unpack(args.maxs or {1, 1, 1}))
	
	-- TODO OK this is a little ambiguous ...
	-- gridSize is the desired grid size
	-- self.gridSize is gonna be that, plus numGhost on either side ...
	-- so maybe I should rename 'self.gridSize' to 'self.gridSizeWithBorder'
	-- and 'self.gridSizeWithoutBorder' to 'self.gridSize'

	local gridSize = assert(args.gridSize)
	if type(gridSize) == 'number' then 
		self.gridSize = vec3sz(gridSize,1,1)
	elseif type(gridSize) == 'cdata' 
	and ffi.typeof(gridSize) == vec3sz 
	then
		self.gridSize = vec3sz(gridSize)
	elseif type(gridSize) == 'table' then
		self.gridSize = vec3sz(table.unpack(gridSize))
	else
		error("can't understand args.gridSize type "..type(args.gridSize).." value "..tostring(args.gridSize))
	end

	for i=0,self.dim-1 do self.gridSize:ptr()[i] = self.gridSize:ptr()[i] + 2 * self.numGhost end
	for i=self.dim,2 do self.gridSize:ptr()[i] = 1 end

	
	self:createBoundaryOptions()
	self:finalizeBoundaryOptions()


	self.boundaryMethods = {}
	for i=1,3 do
		for _,minmax in ipairs(minmaxs) do
			local var = xNames[i]..minmax
			self.boundaryMethods[var] = self.boundaryOptions:find(
				(args.boundary or {})[var] or 'freeflow',
				function(option, search)
					return search == next(option)
				end)
		end
	end


	self.usePLM = args.usePLM
	assert(not self.usePLM or self.fluxLimiter == 1, "are you sure you want to use flux and slope limiters at the same time?")
	self.slopeLimiter = self.app.limiterNames:find(args.slopeLimiter) or 1
	
	self.useCTU = args.useCTU
	assert(not self.useCTU or self.usePLM, "if you are using CTU then you need to use a slope limiter method")

	-- my kernel objs are going to need workgroup info based on domain.size-2*noGhost as well as domain.size ... 
	-- ...and rather than require an extra argument, I think I'll just take advantage of a closure
	local solver = self
	local Program = class(require 'cl.obj.program')
	function Program:init(args)
		args.env = solver.app.env
		args.domain = solver.domain
		Program.super.init(self, args)
	end
	self.Program = Program
end


-- this is only used by GridSolver, but each kernel gets its own,
-- so TODO get rid of this and just use kernel's sizes
function GridSolver:getSizePropsForWorkGroupSize(maxWorkGroupSize)
	local localSize1d = math.min(maxWorkGroupSize, tonumber(self.gridSize:volume()))

	local localSize2d
	if self.dim == 3 then
		local localSizeX = math.min(tonumber(self.gridSize.x), 2^math.ceil(math.log(maxWorkGroupSize,2)/2))
		local localSizeY = maxWorkGroupSize / localSizeX
		localSize2d = {localSizeX, localSizeY}
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
		sizeWithoutBorder:ptr()[i] = sizeWithoutBorder:ptr()[i] - 2 * self.numGhost
	end
	local volumeWithoutBorder = tonumber(sizeWithoutBorder:volume())
	
	local numCells = tonumber(self.gridSize:volume())

	local stepSize = vec3sz()
	stepSize.x = 1
	for i=1,self.dim-1 do
		stepSize:ptr()[i] = stepSize:ptr()[i-1] * self.gridSize:ptr()[i-1]
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
	}
end

-- call this when the solver initializes or changes the codePrefix (or changes initState)
-- it will build the code prefix and refresh everything related to it
-- TODO if you change cons_t then call resetState etc (below the refreshEqnInitState() call a few lines above) in addition to this -- or else your values will get messed up
function GridSolver:refreshEqnInitState()	
	GridSolver.super.refreshEqnInitState(self)

	-- bounds don't get set until initState() is called, but code prefix needs them ...
	-- TODO do a proper refresh so mins/maxs can be properly refreshed
	local initState = self.eqn.initState
	if initState.mins then self.mins = vec3(table.unpack(initState.mins)) end
	if initState.maxs then self.maxs = vec3(table.unpack(initState.maxs)) end
end

-- call this when a gui var changes
-- it rebuilds the code prefix, but doesn't reset the initState
function GridSolver:refreshCodePrefix()
	GridSolver.super.refreshCodePrefix(self)	-- refresh integrator
	-- changing initState calls this, and could change boundary programs, so I'm putting this here
	-- bad excuse, I know
	self:refreshBoundaryProgram()
end


function GridSolver:initDraw()
	local GLProgram = require 'gl.program'
	
	local heatMapCode = file['draw/2d_heatmap.shader']
	self.heatMap2DShader = GLProgram{
		vertexCode = template(heatMapCode, {
			solver = self,
			vertexShader = true,
		}),
		fragmentCode = template(heatMapCode, {
			fragmentShader = true,
			solver = self,
		}),
		uniforms = {
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
		},
	}

	if self.dim == 3 then
		-- raytracing (stalling)
		
		local maxiter = math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y), tonumber(self.gridSize.z))
		local volumetricCode = file['draw/volumetric.shader']
		self.volumeRayShader = GLProgram{
			vertexCode = template(volumetricCode, {
				solver = self,
				vertexShader = true,
			}),
			fragmentCode = template(volumetricCode, {
				solver = self,
				fragmentShader = true,
			}),
			uniforms = {
				tex = 0,
				gradient = 1,
				maxiter = maxiter,
				oneOverDx = {(self.maxs - self.mins):unpack()},
			},
		}

		-- volume slices

		local volumeSliceCode = file['draw/slices3d.shader']
		self.volumeSliceShader = GLProgram{
			vertexCode = template(volumeSliceCode, {
				solver = self,
				vertexShader = true
			}),
			fragmentCode = template(volumeSliceCode, {
				solver = self,
				fragmentShader = true,
				-- TODO move this from app, or make it a field of app?
				clipInfos = useClipPlanes and clipInfos or nil,
			}),
			uniforms = {
				volTex = 0,
				gradientTex = 1,
				valueMin = 0,
				valueMax = 0,
			},
		}
	end

	local vectorFieldCode = file['draw/vectorfield.shader']
	self.vectorFieldShader = GLProgram{
		vertexCode = template(vectorFieldCode, {
			solver = self,
			vertexShader = true,
		}),
		fragmentCode = template(vectorFieldCode, {
			solver = self,
			fragmentShader = true,
		}),
		uniforms = {
			scale = 1,
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
		},
	}
end


-- my best idea to work around the stupid 8-arg max kernel restriction
-- this is almost as bad of an idea as using OpenCL was to begin with
GridSolver.allocateOneBigStructure = false

function GridSolver:getConsLRTypeCode()
	return template([[
typedef union {
	<?=eqn.cons_t?> LR[2];
	struct {
		<?=eqn.cons_t?> L, R;
	};
} <?=eqn.consLR_t?>;

//ugly hack to work around even uglier boundary code
typedef struct {
	<?=eqn.consLR_t?> side[<?=solver.dim?>];
} <?=eqn.consLR_t?>_dim;
]], {
	solver = self,
	eqn = self.eqn,
})
end

function GridSolver:createBuffers()
	local realSize = ffi.sizeof(self.app.real)

	-- to get sizeof
	ffi.cdef(self:getConsLRTypeCode())

	-- for twofluid, cons_t has been renamed to euler_maxwell_t and maxwell_cons_t
	if ffi.sizeof(self.eqn.cons_t) ~= self.eqn.numStates * realSize then
	   error('Expected sizeof('..self.eqn.cons_t..') to be '
		   ..self.eqn.numStates..' * sizeof(real) = '..(self.eqn.numStates * realSize)
		   ..' but found '..ffi.sizeof(self.eqn.cons_t)..' = '..(ffi.sizeof(self.eqn.cons_t) / realSize)..' * sizeof(real). '
		   ..'Maybe you need to update Eqn.numStates?')
	end
	if ffi.sizeof(self.eqn.cons_t) < ffi.sizeof(self.eqn.prim_t) then
		error("for PLM's sake I might need sizeof(prim_t) <= sizeof(cons_t)")
	end

	-- should I put these all in one AoS?

if tryingAMR == 'gradient' then
	--[[
	ok here's my thoughts on the size ...
	I'm gonna try to do AMR
	that means storing the nodes in the same buffer
	I think I'll pad the end of the UBuf with the leaves
	and then store a tree of index information somewhere else that says what leaf goes where
	(maybe that will go at the end)
	
	how should memory breakdown look?
	how big should the leaves be?

	how about do like reduce ...
	leaves can be 16x16 blocks
	a kernel can cycle through the root leaves and sum te amrError 
	-- in the same kernel as it is calculated
	-- then I don't need to store so much memory, only one value per leaf, not per grid ...
	-- then ... split or merge leaves, based on their error ...

	how to modify that?  in anothe kernel of its own ...
	bit array of whether each cell is divided ...
	-- then in one kernel we update all leaves
	-- in another kernel, populate leaf ghost cells
	-- in another kernel, 
	-- 		decide what should be merged and, based on that, copy from parent into this
	--		if something should be split then look for unflagged children and copy into it
	
	how many bits will we need?
	volume / leafSize bits for the root level
	
	then we have parameters of 
	- how big each node is
	- what level of refinement it is
	
	ex:
	nodes in the root level are 2^4 x 2^4 = 2^8 = 256 cells = 256 bits = 2^5 = 32 bytes
	

	leafs are 2^4 x 2^4 = 2^8 = 256 cells
	... but ghost cells are 2 border, so we need to allocate (2^n+2*2)^2 cells ... for n=4 this is 400 cells ... 
		so we lose (2^n+2*2)^2 - 2^(2n) = 2^(2n) + 2^(n+3) + 2^4 - 2^(2n) = 2^(n+3) + 2^4) cells are lost
	
	so leafs multipy at a factor of 2^2 x 2^2 = 2^4 = 16
	so the next level has 2^(8+4) = 2^12 = 4096 bits = 2^9 = 512 bytes

	--]]

	-- the size, in cells, which a node replaces
	self.amrNodeFromSize = ({
		vec3sz(8, 1, 1),
		vec3sz(8, 8, 1),
		vec3sz(8, 8, 8),
	})[self.dim]
print('self.amrNodeFromSize', self.amrNodeFromSize)

	-- the size, in cells, of each node, excluding border, for each dimension
	self.amrNodeSizeWithoutBorder = ({
		vec3sz(16, 1, 1),
		vec3sz(16, 16, 1),
		vec3sz(16, 16, 16),
	})[self.dim]
print('self.amrNodeSizeWithoutBorder', self.amrNodeSizeWithoutBorder)

	-- size of the root level, in terms of nodes ('from' size)
	self.amrRootSizeInFromSize = vec3sz(1,1,1)
	for i=0,self.dim-1 do
		self.amrRootSizeInFromSize:ptr()[i] = 
			roundup(self.sizeWithoutBorder:ptr()[i], self.amrNodeFromSize:ptr()[i]) 
				/ self.amrNodeFromSize:ptr()[i]
	end
print('self.amrRootSizeInFromSize', self.amrRootSizeInFromSize)

	-- how big each node is
	self.amrNodeSize = self.amrNodeSizeWithoutBorder + 2 * self.numGhost

	-- how many nodes to allocate and use
	-- here's the next dilemma in terms of memory layout
	-- specifically in terms of rendering
	-- if I want to easily copy and render the texture information then I will need to package the leafs into the same texture as the root
	-- which means extending the texture buffer in some particular direction.
	-- since i'm already adding the leaf information to the end of the buffer
	-- and since appending to the end of a texture buffer coincides with adding extra rows to the texture
	-- why not just put our leafs in extra rows of -- both our display texture and of  
	self.amrMaxNodes = 1
	
	-- this will hold info on what leafs have yet been used
	self.amrLeafs = table()

	-- hmm, this is the info for the root node ...
	-- do I want to keep the root level data separate?
	-- or do I just want to represent everything as a collection of leaf nodes?
	-- I'll keep the root structure separate for now
	-- so I can keep the original non-amr solver untouched
	self.amrLayers = table()
	self.amrLayers[1] = table()	-- here's the root
end

	-- this much for the first level U data
	local UBufSize = self.numCells * ffi.sizeof(self.eqn.cons_t)
if tryingARM == 'gradient' then	
	UBufSize = UBufSize + self.amrMaxNodes * self.amrNodeSize:volume()
end	
	self:clalloc('UBuf', UBufSize)


if tryingAMR == 'dt vs 2dt' then
	-- here's my start at AMR, using the 1989 Berger, Collela two-small-steps vs one-big-step method
	self:clalloc('lastUBuf', self.numCells * ffi.sizeof(self.eqn.cons_t))
	self:clalloc('U2Buf', self.numCells * ffi.sizeof(self.eqn.cons_t))
elseif tryingAMR == 'gradient' then
	
	-- this is going to be a single value for each leaf
	-- that means the destination will be the number of nodes it takes to cover the grid (excluding the border)
	-- however, do I want this to be a larger buffer, and then perform reduce on it?
	self:clalloc('amrErrorBuf', 
		-- self.volume 
		tonumber(self.amrRootSizeInFromSize:volume())
		* ffi.sizeof(self.app.real))
end

	if self.usePLM then
		self:clalloc('ULRBuf', self.numCells * self.dim * ffi.sizeof(self.eqn.consLR_t))
	end

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar 
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', self.numCells * realSize * 3)
	local reduceSwapBufSize = roundup(self.numCells * realSize / self.localSize1d, realSize)
	self:clalloc('reduceSwapBuf', reduceSwapBufSize)
	self.reduceResultPtr = ffi.new('real[1]', 0)

	-- CL/GL interop

	-- hmm, notice I'm still keeping the numGhost border on my texture 
	-- if I remove the border altogether then I get wrap-around
	-- maybe I should just keep a border of 1?
	-- for now i'll leave it as it is
	local GLTex2D = require 'gl.tex2d'
	local GLTex3D = require 'gl.tex3d'
	local cl = self.dim < 3 and GLTex2D or GLTex3D
	self.tex = cl{
		width = tonumber(self.gridSize.x),
		height = tonumber(self.gridSize.y),
		depth = tonumber(self.gridSize.z),
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		--magFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
	}

	local CLImageGL = require 'cl.imagegl'
	if self.app.useGLSharing then
		self.texCLMem = CLImageGL{context=self.app.ctx, tex=self.tex, write=true}
	else
		self.calcDisplayVarToTexPtr = ffi.new(self.app.real..'[?]', self.numCells * 3)
		
		--[[ PBOs?
		self.calcDisplayVarToTexPBO = ffi.new('gl_int[1]', 0)
		gl.glGenBuffers(1, self.calcDisplayVarToTexPBO)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, self.calcDisplayVarToTexPBO[0])
		gl.glBufferData(gl.GL_PIXEL_UNPACK_BUFFER, self.tex.width * self.tex.height * ffi.sizeof(self.app.real) * 4, nil, gl.GL_STREAM_READ)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
		--]]
	end
end

function GridSolver:createCodePrefix()
	GridSolver.super.createCodePrefix(self)
	
	local lines = table{
		self.codePrefix,
	}

	lines:append{
		'#define numGhost '..self.numGhost,
	}:append(xNames:map(function(x,i)
	-- coordinate space = u,v,w
	-- cartesian space = x,y,z
	-- min and max in coordinate space
		return '#define mins_'..x..' '..clnumber(self.mins[i])..'\n'
			.. '#define maxs_'..x..' '..clnumber(self.maxs[i])..'\n'
	end)):append{
		'constant real3 mins = _real3(mins_x, '..(self.dim<2 and '0' or 'mins_y')..', '..(self.dim<3 and '0' or 'mins_z')..');', 
		'constant real3 maxs = _real3(maxs_x, '..(self.dim<2 and '0' or 'maxs_y')..', '..(self.dim<3 and '0' or 'maxs_z')..');', 
	}:append(xNames:map(function(name,i)
	-- grid size
		return '#define gridSize_'..name..' '..tonumber(self.gridSize[name])
	end)):append{
		'constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);',
		'constant int4 stepsize = (int4)(1, gridSize_x, gridSize_x * gridSize_y, gridSize_x * gridSize_y * gridSize_z);',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		indexForInt4ForSize(i, gridSize_x, gridSize_y, gridSize_z)',
	
	--[[
	naming conventions ...
	* the grid indexes i_1..i_n that span 1 through gridSize_1..gridSize_n
	(between the index and the coordinate space:)
		- grid_dx? is the change in coordinate space wrt the change in grid space
		- cell_x(i) calculates the coordinates at index i
	* the coordinates space x_1..x_n that spans mins.s1..maxs.sn
	(between the coordinate and the embedded space:)
		- vectors can be calculated from Cartesian by cartesianToCoord
		- the length of the basis vectors wrt the change in indexes is given by dx?_at(x)
		- the Cartesian length of the holonomic basis vectors is given by coordHolBasisLen.  
			This is like dx?_at except not scaled by grid_dx?
			This is just the change in embedded wrt the change in coordinate, not wrt the change in grid
		- volume_at(x) gives the volume between indexes at the coordinate x
		- the Cartesian length of a vector in coordinate space is given by coordLen and coordLenSq
	* the embedded Cartesian space ... idk what letters I should use for this.  
		Some literature uses x^I vs coordinate space x^a, but that's still an 'x', no good for programming.
		Maybe I'll use 'xc' vs 'x', like I've already started to do in the initial conditions.
	--]]
	
	}:append(range(3):map(function(i)
	-- this is the change in coordinate wrt the change in code
	-- delta in coordinate space along one grid cell
		if i > self.dim then 
-- why is this causing errors?
-- I think the grid_dx def is bad for dims outside the grid dim ...
--			return '#define grid_dx'..(i-1)..' 1.' 
		end
		return (('#define grid_dx{i} ((maxs_{x} - mins_{x}) / (real)(gridSize_{x} - '..(2*self.numGhost)..'))')
			:gsub('{i}', i-1)
			:gsub('{x}', xNames[i]))
	end)):append(range(3):map(function(i)
	-- mapping from index to coordinate 
		return (('#define cell_x{i}(i) ((real)(i + '..clnumber(.5-self.numGhost)..') * grid_dx{i} + mins_'..xNames[i]..')')
			:gsub('{i}', i-1))
	end)):append{
		'#define cell_x(i) _real3(cell_x0(i.x), cell_x1(i.y), cell_x2(i.z))',
	
		-- bounds-check macro
		'#define OOB(lhs,rhs) (i.x < (lhs) || i.x >= gridSize_x - (rhs)'
			.. (self.dim < 2 and '' or ' || i.y < (lhs) || i.y >= gridSize_y - (rhs)')
			.. (self.dim < 3 and '' or ' || i.z < (lhs) || i.z >= gridSize_z - (rhs)')
			.. ')',
		
		template([[

// define i, index, and bounds-check
#define SETBOUNDS(lhs,rhs)	\
	int4 i = globalInt4(); \
	if (OOB(lhs,rhs)) return; \
	int index = INDEXV(i);
		
// same as above, except for kernels that don't use the boundary
// index operates on buffers of 'gridSize' (with border)
// but the kernel must be invoked across sizeWithoutBorder
<? local range = require'ext.range'
?>#define SETBOUNDS_NOGHOST() \
	int4 i = globalInt4(); \
	if (OOB(0,2*numGhost)) return; \
	i += (int4)(<?=range(4):map(function(i) 
		return i <= solver.dim and '2' or '0' 
	end):concat','?>); \
	int index = INDEXV(i);
]], {solver = self}),
	}

	
	-- volume of a cell = volume element times grid dx's 
	lines:insert(template([[
real volume_at(real3 x) {
	return sqrt_det_g_grid(x)<?
for i=0,solver.dim-1 do
?> * grid_dx<?=i?><?
end
?>;
}
]], {solver = self}))

	lines:append{
		-- not messing with this one yet
		self.allocateOneBigStructure and '#define allocateOneBigStructure' or '',
		
		self:getConsLRTypeCode(),
	}

	self.codePrefix = lines:concat'\n'

print'done building solver.codePrefix'
--print(self.codePrefix)
end

function GridSolver:resetState()
	self.app.cmds:finish()
		
	-- start off by filling all buffers with zero, just in case initState forgets something ...
	for _,bufferInfo in ipairs(self.buffers) do
		self.app.cmds:enqueueFillBuffer{buffer=self[bufferInfo.name], size=bufferInfo.size}
	end
	
	self:boundary()
	if self.eqn.useConstrainU then
		self.constrainUKernelObj()
	end

	GridSolver.super.resetState(self)
end

function GridSolver:getSolverCode()
	local slopeLimiterCode = 'real slopeLimiter(real r) {'
		.. self.app.limiters[self.slopeLimiter].code 
		.. '}'
	
	local amrCode = template(({
		['dt vs 2dt'] = [[
kernel void compareUvsU2(
	global <?=eqn.cons_t?>* U2Buf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?> *U2 = U2Buf + index;
	const global <?=eqn.cons_t?> *U = UBuf + index;
	
	//what to use to compare values ...
	//if we combine all primitives, they'll have to be appropriately weighted ...
	real sum = 0.;
	real tmp;
<? for i=0,eqn.numStates-1 do
?>	tmp = U2->ptr[<?=i?>] - U->ptr[<?=i?>]; sum += tmp * tmp;
<? end
?>	U2->ptr[0] = sum * 1e+5;
}
]],
		gradient = [==[
kernel void calcAMRError(
	global real* amrErrorBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	int4 nodei = globalInt4();
	if (nodei.x >= <?=solver.amrRootSizeInFromSize.x?> || 
		nodei.y >= <?=solver.amrRootSizeInFromSize.y?>) 
	{
		return;
	}

	int nodeIndex = nodei.x + <?=solver.amrRootSizeInFromSize.x?> * nodei.y;

	real dV_dx;	
	real sum = 0.;
	
	//hmm, it's less memory, but it's probably slower to iterate across all values as I build them here
	for (int nx = 0; nx < <?=solver.amrNodeFromSize.x?>; ++nx) {
		for (int ny = 0; ny < <?=solver.amrNodeFromSize.y?>; ++ny) {
			int4 Ui = (int4)(0,0,0,0);
			
			Ui.x = nodei.x * <?=solver.amrNodeFromSize.x?> + nx + numGhost;
			Ui.y = nodei.y * <?=solver.amrNodeFromSize.y?> + ny + numGhost;
			
			int Uindex = INDEXV(Ui);
			const global <?=eqn.cons_t?>* U = UBuf + Uindex;
				
	//TODO this wasn't the exact formula ...
	// and TODO make this modular.  some papers use velocity vector instead of density.  
	// why not total energy -- that incorporates everything?
<? for i=0,solver.dim-1 do
?>			dV_dx = (U[stepsize.s<?=i?>].rho - U[-stepsize.s<?=i?>].rho) / (2. * grid_dx<?=i?>);
			sum += dV_dx * dV_dx;
<? end
?>	
		}
	}
	amrErrorBuf[nodeIndex] = sum * 1e-2 * <?=clnumber(1/tonumber( solver.amrNodeFromSize:volume() ))?>;
}

//from is the position on the root level to read from
//to is which node to copy into
kernel void initNodeFromRoot(
	global <?=eqn.cons_t?>* UBuf,
	int4 from,
	int toNodeIndex
) {
	int4 i = (int4)(0,0,0,0);
	i.x = get_global_id(0);
	i.y = get_global_id(1);
	int dstIndex = i.x + numGhost + <?=solver.amrNodeSize.x?> * (i.y + numGhost);
	int srcIndex = from.x + (i.x>>1) + numGhost + gridSize_x * (from.y + (i.y>>1) + numGhost);

	global <?=eqn.cons_t?>* dstU = UBuf + <?=solver.numCells?> + toNodeIndex * <?=solver.amrNodeSize:volume()?>;
	
	//blitter srcU sized solver.amrNodeFromSize (in a patch of size solver.gridSize)
	// to dstU sized solver.amrNodeSize (in a patch of solver.amrNodeSize)
	
	dstU[dstIndex] = UBuf[srcIndex];
}
]==],
	})[tryingAMR] or '', {
		solver = self,
		eqn = self.eqn,
		clnumber = clnumber,
	})

	return table{
		GridSolver.super.getSolverCode(self),
			
		slopeLimiterCode,
		
		-- messing with this ...
		self.usePLM and template(file['solver/plm.cl'], {solver=self, eqn=self.eqn}) or '',
		self.useCTU and template(file['solver/ctu.cl'], {solver=self, eqn=self.eqn}) or '',
		
		amrCode,
	}:concat'\n'
end

-- depends on buffers
function GridSolver:refreshSolverProgram()
	-- set pointer to the buffer holding the LR state information
	-- for piecewise-constant that is the original UBuf
	-- for piecewise-linear that is the ULRBuf
	self.getULRBuf = self.usePLM and self.ULRBuf or self.UBuf

	self.getULRArg = self.usePLM 
		and ('const global '..self.eqn.consLR_t..'* ULRBuf')
		or ('const global '..self.eqn.cons_t..'* UBuf')

	-- this code creates the const global cons_t* UL, UR variables
	-- it assumes that indexL, indexR, and side are already defined
	-- both UBuf and ULRBuf should be cell-centered
	if self.usePLM then
		self.getULRCode = function(self, args)
			args = args or {}
			local suffix = args.suffix or ''
			return template([[
	const global <?=eqn.cons_t?>* UL<?=suffix?> = &ULRBuf[side + dim * <?=indexL?>].R;
	const global <?=eqn.cons_t?>* UR<?=suffix?> = &ULRBuf[side + dim * <?=indexR?>].L;
]],			{
				eqn = self.eqn,
				suffix = suffix,
				indexL = args.indexL or 'indexL'..suffix,
				indexR = args.indexR or 'indexR'..suffix,
			})
		end
	else 
		self.getULRCode = function(self, args)
			args = args or {}
			local suffix = args.suffix or ''
			return template([[
	const global <?=eqn.cons_t?>* UL<?=suffix?> = UBuf + <?=indexL?>;
	const global <?=eqn.cons_t?>* UR<?=suffix?> = UBuf + <?=indexR?>;
]],			{
				eqn = self.eqn,
				suffix = suffix,
				indexL = args.indexL or 'indexL'..suffix,
				indexR = args.indexR or 'indexR'..suffix,
			})
		end
	end

	local code = self:getSolverCode()

	time('compiling solver program', function()
		self.solverProgramObj = self.Program{code=code}
		self.solverProgramObj:compile()
	end)

	self:refreshCalcDTKernel()

	if self.eqn.useConstrainU then
		self.constrainUKernelObj = self.solverProgramObj:kernel('constrainU', self.UBuf)
	end

	if self.usePLM then
		self.calcLRKernelObj = self.solverProgramObj:kernel('calcLR', self.ULRBuf, self.UBuf)
	end
	if self.useCTU then
		-- currently implemented in solver/roe.cl
		-- not available for any other flux method
		assert(self.fluxBuf)
		self.updateCTUKernelObj = self.solverProgramObj:kernel('updateCTU', self.ULRBuf, self.fluxBuf)
	end

if tryingAMR == 'dt vs 2dt' then
	self.compareUvsU2KernelObj = self.solverProgramObj:kernel('compareUvsU2', self.U2Buf, self.UBuf)
elseif tryingAMR == 'gradient' then
	self.calcAMRErrorKernelObj = self.solverProgramObj:kernel('calcAMRError', self.amrErrorBuf, self.UBuf)
	self.initNodeFromRootKernelObj = self.solverProgramObj:kernel('initNodeFromRoot', self.UBuf)
end

	for _,op in ipairs(self.ops) do
		op:refreshSolverProgram()
	end

	-- display stuff:

	if self.app.useGLSharing then
		for _,displayVarGroup in ipairs(self.displayVarGroups) do
			for _,var in ipairs(displayVarGroup.vars) do
				--[[
				if var.enabled 
				or (var.vecVar and var.vecVar.enabled)
				then
				--]]do
					var.calcDisplayVarToTexKernelObj = self.solverProgramObj:kernel('calcDisplayVarToTex_'..var.id, self.texCLMem)
				end
			end
		end
	end

	for _,displayVarGroup in ipairs(self.displayVarGroups) do
		for _,var in ipairs(displayVarGroup.vars) do
			--[[
			if var.enabled 
			or (var.vecVar and var.vecVar.enabled)
			then
			--]]do
				var.calcDisplayVarToBufferKernelObj = self.solverProgramObj:kernel('calcDisplayVarToBuffer_'..var.id, self.reduceBuf)
			end
		end
	end
end


-------------------------------------------------------------------------------
--                                 display vars                              --
-------------------------------------------------------------------------------


-- subclass and override
local DisplayVar = class(GridSolver.DisplayVar)
GridSolver.DisplayVar = DisplayVar


-- and this is the DisplayVar used for UBuf
GridSolver.DisplayVar_U = DisplayVar

function GridSolver:addDisplayVars()
	GridSolver.super.addDisplayVars(self)
if tryingAMR == 'dt vs 2dt' then
	self:addDisplayVarGroup{
		name = 'U2',
		type = self.eqn.cons_t,
		vars = {
			{[0] = '*value = buf[index].ptr[0];'},
		}
	}
elseif tryingAMR == 'gradient' then
	self:addDisplayVarGroup{
		name = 'amrError',
		type = 'real',
		vars = {
			{[0] = '*value = buf[index];'},
		}
	}
end
end


-------------------------------------------------------------------------------
--                              boundary                                     --
-------------------------------------------------------------------------------

--[[
boundaryOptions is a table of {name = args => assign code}
	args of the boundary function are:
		index
		assign = function(dst,src) 
			assignment operator
			default is dst = src
			poisson uses dist.field = src.field
		field
		array
		rhs
		gridSizeSide = solver.gridSize[side]
		side = 1,2,3
		minmax = 'min' or 'max'
		mirrorVars = which vars to reflect.
			for solver.UBuf this is taken from eqn
			for poisson this is nothing

this is such a mess.  it's practically an AST.
--]]
function GridSolver:createBoundaryOptions()
	local indent = '\t\t'
	self.boundaryOptions = table{
		{periodic = function(args)
			local gridSizeSide = 'gridSize_'..xNames[args.side]
			if args.minmax == 'min' then
				return indent..args.assign(
					args.array('buf', args.index'j'), 
					args.array('buf', args.index(gridSizeSide..'-2*numGhost+j'))
				)..';'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return indent..args.assign(
					args.array('buf', args.index(rhs)),
					args.array('buf', args.index'numGhost+j')
				)..';'
			end
		end},
		{mirror = function(args)
			local gridSizeSide = 'gridSize_'..xNames[args.side]
			if args.minmax == 'min' then
				return table{
					indent..args.assign(
						args.array('buf', args.index'j'),
						args.array('buf', args.index'2*numGhost-1-j')
					)..';'
				}:append(table.map((args.mirrorVars or {})[args.side] or {}, function(var)
					return indent..args.assign(
						args.field(args.array('buf', args.index'j'), var),
						args.field('-'..args.array('buf', args.index'j'), var)
					)..';'
				end)):concat'\n'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return table{
					indent..args.assign(
						args.array('buf', args.index(rhs)),
						args.array('buf', args.index(gridSizeSide..'-numGhost-1-j'))
					)..';'
				}:append(table.map((args.mirrorVars or {})[args.side] or {}, function(var)
					return indent..args.assign(
						args.field(args.array('buf', args.index(rhs)), var),
						args.field('-'..args.array('buf', args.index(rhs)), var)
					)..';'
				end)):concat'\n'
			end
		end},
		{freeflow = function(args)
			local gridSizeSide = 'gridSize_'..xNames[args.side]
			if args.minmax == 'min' then
				return indent..args.assign(
					'buf['..args.index'j'..']',
					'buf['..args.index'numGhost'..']'
				)..';'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return indent..args.assign(
					'buf['..args.index(rhs)..']',
					'buf['..args.index(gridSizeSide..'-numGhost-1')..']'
				)..';'
			end
		end},
	}

	if self.eqn.createBoundaryOptions then
		self.eqn:createBoundaryOptions()
	end
end
function GridSolver:finalizeBoundaryOptions()
	self.boundaryOptionNames = self.boundaryOptions:map(function(option) return (next(option)) end)
	-- hmm here's the one time that table.map using k,v comes in handy
	self.boundaryOptionForName = self.boundaryOptions:map(function(option) local k,v = next(option) return v,k end)
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
				name = args	
			elseif type(args) == 'table' then
				name = args[k]
			else
				error("don't know how to deal with these args") 
			end
			local i = self.boundaryOptions:find(nil, function(option) return next(option) == name end)
			if not i then
				io.stderr:write("unable to find boundary method "..tostring(name)
					..' ('..type(name)..')'
					.." -- can't assign it to side "..k.."\n")
				io.stderr:flush()
			else
				self.boundaryMethods[k] = i
			end
		end
	end
end

--[[
args:
	type = type of buffer to create boundary information
	extraArgs = any extra args to add to the kernel (other than a global of the buffer type)
	methods = {
		[x|y|z..min|max] = 
			either 'peroidic', 'mirror', 'freeflow',
			or a function of the name of variable pointer 
			that returns the code to apply to that particular pointer 
	}
	mirrorVars = {x vars, y vars, z vars} = table of tables of strings of what fields should be negative'd on mirror condition
		this tends to be vectors' components that point into the boundary.

	assign = function(a,b) a = b 
	field = function(a,b) a.b
	array = function(a,b) a[b]
--]]
function GridSolver:createBoundaryProgramAndKernel(args)
	local assign = args.assign or function(a, b) return a .. ' = ' .. b end
	local field = args.field or function(a, b) return a .. '.' .. b end
	local array = args.array or function(a, b) return a .. '[' .. b .. ']' end

	local lines = table()
	lines:insert(self.codePrefix)
	
	for side=1,self.dim do

		lines:insert(template([[
kernel void boundary_<?=xNames[side]?>(
	global <?=args.type?>* buf
<?= args.extraArgs and #args.extraArgs > 0 
	and ','..table.concat(args.extraArgs, ',\n\t')
	or '' ?>
) {
]], {
		args = args,
		side = side, 
		xNames = xNames,
	}))
		if self.dim == 2 then
			lines:insert[[
	int i = get_global_id(0);
]]
		elseif self.dim == 3 then
			lines:insert[[
	int2 i = (int2)(get_global_id(0), get_global_id(1));
]]
		end
	
		-- 1D: use a small 1D kernel and just run once 
		-- 2D: use a 1D kernel the size of the max dim 
			
		if self.dim == 2 then
			if side == 1 then
				lines:insert[[
	if (i >= gridSize_y) return;
]]
			elseif side == 2 then
				lines:insert[[
	if (i >= gridSize_x) return;
]]
			end
		elseif self.dim == 3 then
			if side == 1 then
				lines:insert[[
	if (i.x >= gridSize_y || i.y >= gridSize_z) return;
]]
			elseif side == 2 then
				lines:insert[[
	if (i.x >= gridSize_x || i.y >= gridSize_z) return;
]]
			elseif side == 3 then
				lines:insert[[
	if (i.x >= gridSize_x || i.y >= gridSize_y) return;
]]
			end
		end

		lines:insert[[
	for (int j = 0; j < numGhost; ++j) {
]]


		local function indexv(j)
			if self.dim == 1 then
				return j..',0,0'
			elseif self.dim == 2 then
				if side == 1 then
					return j..',i,0'
				elseif side == 2 then
					return 'i,'..j..',0'
				end
			elseif self.dim == 3 then
				if side == 1 then
					return j..',i.x,i.y'
				elseif side == 2 then
					return 'i.x,'..j..',i.y'
				elseif side == 3 then
					return 'i.x,i.y,'..j
				end
			end
		end

		local function index(j)
			return 'INDEX('..indexv(j)..')'
		end
	
		for _,minmax in ipairs(minmaxs) do
			local method = args.methods[xNames[side]..minmax]
			lines:insert(method{
				index = index,
				indexv = indexv,
				assign = assign,
				array = array,
				field = field,
				side = side,
				mirrorVars = args.mirrorVars,
				minmax = minmax,
			})
		end
	
		lines:insert'	}'
		lines:insert'}'
	end

	local code = lines:concat'\n'

	local boundaryProgramObj
	time('compiling boundary program', function()
		boundaryProgramObj = self.Program{code=code}
		boundaryProgramObj:compile()
	end)
	local boundaryKernelObjs = table()
	for i=1,self.dim do
		boundaryKernelObjs:insert(boundaryProgramObj:kernel('boundary_'..xNames[i]))
	end

	-- TODO switch these over to obj
	return boundaryProgramObj, boundaryKernelObjs
end

function GridSolver:getBoundaryProgramArgs()
	return {
		type = self.eqn.cons_t,
		-- remap from enum/combobox int values to functions from the solver.boundaryOptions table
		methods = table.map(self.boundaryMethods, function(v)
			if type(v) == 'function' then return v end
			if not self.boundaryOptions[v] then
				error("failed to find boundaryOption index "..v)
			end
			return (select(2, next(self.boundaryOptions[v])))
		end),
		mirrorVars = self.eqn.mirrorVars,
	}
end

function GridSolver:refreshBoundaryProgram()
	self.boundaryProgramObj, self.boundaryKernelObjs = self:createBoundaryProgramAndKernel(self:getBoundaryProgramArgs())
	for _,obj in ipairs(self.boundaryKernelObjs) do
		obj.obj:setArg(0, self.UBuf)
	end
	for _,op in ipairs(self.ops) do
		if op.refreshBoundaryProgram then
			op:refreshBoundaryProgram()
		end
	end

	if self.useCTU then
		self.lrBoundaryProgramObj, self.lrBoundaryKernelObjs = self:createBoundaryProgramAndKernel(table(self:getBoundaryProgramArgs(), {
			type = self.eqn.consLR_t..'_dim',
		}))
		for _,obj in ipairs(self.lrBoundaryKernelObjs) do
			obj.obj:setArg(0, self.ULRBuf)
		end
	end
end

-- assumes the buffer is already in the kernel's arg
function GridSolver:applyBoundaryToBuffer(kernelObjs)
	for side,obj in ipairs(kernelObjs) do	
		-- 1D:
		if self.dim == 1 then
			-- if you do change this size from anything but 1, make sure to add conditions to the boundary kernel code
			self.app.cmds:enqueueNDRangeKernel{
				kernel = obj.obj,
				globalSize = 1,
				localSize = 1,
			}
		elseif self.dim == 2 then
			local localSize = math.min(self.localSize1d, self.maxWorkGroupSize)
			local maxSize = roundup(
					side == 1 
					and tonumber(self.gridSize.y)
					or tonumber(self.gridSize.x)
					,
				localSize)
			self.app.cmds:enqueueNDRangeKernel{
				kernel = obj.obj,
				globalSize = maxSize,
				localSize = localSize,
			}
		elseif self.dim == 3 then
			-- xy xz yz
			local maxSizeX = roundup(
				math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y)),
				self.localSize2d[1])
			local maxSizeY = roundup(
				math.max(tonumber(self.gridSize.y), tonumber(self.gridSize.z)),
				self.localSize2d[2])
			self.app.cmds:enqueueNDRangeKernel{
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
	self:applyBoundaryToBuffer(self.boundaryKernelObjs)
end


-------------------------------------------------------------------------------
--                                                                           --
-------------------------------------------------------------------------------


function GridSolver:update()
	GridSolver.super.update(self)
	
--local before = self.UBufObj:toCPU()
--print'\nself.UBufObj before boundary:' self:printBuf(self.UBufObj) print(debug.traceback(),'\n\n')	
	self:boundary()
--[[
local after = self.UBufObj:toCPU()
print'\nself.UBufObj after boundary:' self:printBuf(self.UBufObj) print(debug.traceback(),'\n\n')	
if self.app.dim == 2 then
	for i=0,tonumber(self.gridSize.x)-1 do
		for j=0,self.numGhost-1 do
			for _,s in ipairs{0,4} do
				for _,offset in ipairs{
					s + self.eqn.numStates * (i + self.gridSize.x * j),
					s + self.eqn.numStates * (j + self.gridSize.x * i),
					s + self.eqn.numStates * (i + self.gridSize.x * (self.gridSize.y - 1 - j)),
					s + self.eqn.numStates * ((self.gridSize.x - 1 - j) + self.gridSize.x * i),
				} do
					if after[offset] < 0 then
						local msg = "negative found in boundary after :boundary() was called"
						msg = msg .. "\nat position " .. tostring(offset / self.eqn.numStates)
						msg = msg .. '\nit was '..(before[offset] < 0 and 'negative' or 'positive')..' before the boundary update'
						error(msg)
					end
				end
			end
		end
	end
end
--]]	
	local dt = self:calcDT()

if tryingAMR == 'dt vs 2dt' then
	-- back up the last buffer
	self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.lastUBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
end
	
	-- first do a big step
	self:step(dt)

	-- now copy it to the backup buffer
if tryingAMR == 'dt vs 2dt' then
	-- TODO have step() provide a target, and just update directly into U2Buf?
	self.app.cmds:enqueueCopyBuffer{src=self.UBuf, dst=self.U2Buf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}
	self.app.cmds:enqueueCopyBuffer{src=self.lastUBuf, dst=self.UBuf, size=self.numCells * self.eqn.numStates * ffi.sizeof(self.app.real)}

	local t = self.t
	self:step(.5 * dt)
	self.t = t + .5 * dt
	self:step(.5 * dt)
	self.t = t + dt

	-- now compare UBuf and U2Buf, store in U2Buf in the first real of cons_t
	self.compareUvsU2KernelObj()
elseif tryingAMR == 'gradient' then
	
	-- 1) compute errors from gradient, sum up errors in each root node, and output on a per-node basis
	local amrRootSizeInFromGlobalSize = vec3sz(
		roundup(self.amrRootSizeInFromSize.x, self.localSize.x),
		roundup(self.amrRootSizeInFromSize.y, self.localSize.y),
		roundup(self.amrRootSizeInFromSize.z, self.localSize.z))
	
	self.app.cmds:enqueueNDRangeKernel{
		kernel = self.calcAMRErrorKernelObj.obj, 
		dim = self.dim, 
		globalSize = amrRootSizeInFromGlobalSize:ptr(), 
		localSize = self.localSize:ptr(),
	}

	-- 2) based on what nodes' errors are past some value, split or merge...
	--[[
	1) initial tree will have nothing flagged as split
	2) then we get some split data - gradients -> errors -> thresholds -> flags 
		... which are lined up with the layout of the patches ...
		... which doesn't necessarily match the tree structure ...
	3) look through all used patches' error thresholds, min and max
		if it says to split ... 
			then look and see if we have room for any more free leafs in our state buffer
		
			the first iteration will request to split on some cells
			so go through the error buffer for each (root?) node,
			see if the error is bigger than some threshold then this node needs to be split
				then we have to add a new leaf node
			
			so i have to hold a table of what in the U extra leaf buffer is used
			which means looking
		
		if it says to merge ...
			clear the 'used' flag in the overall tree / in the layout of leafs in our state buffer
	--]]
	local vol = tonumber(self.amrRootSizeInFromSize:volume())
	local ptr = ffi.new('real[?]', vol)
	self.app.cmds:enqueueReadBuffer{buffer=self.amrErrorBuf, block=true, size=ffi.sizeof(self.app.real) * vol, ptr=ptr}
-- [[
print'armErrors:'
for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
	for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
		local i = nx + self.amrRootSizeInFromSize.x * ny
		io.write('\t', ('%.5f'):format(ptr[i]))
	end
	print()
end
--]]

for ny=0,tonumber(self.amrRootSizeInFromSize.y)-1 do
	for nx=0,tonumber(self.amrRootSizeInFromSize.x)-1 do
		local i = nx + self.amrRootSizeInFromSize.x * ny
		local nodeErr = ptr[i]
		if nodeErr > .2 then
			print('root node '..tostring(i)..' needs to be split')
			
			-- flag for a split
			-- look for a free node to allocate in the buffer
			-- if there's one available then ...
			-- store it in a map

			local amrLayer = self.amrLayers[1]
			
			-- see if there's an entry in this layer
			-- if there's not then ...
			-- allocate a new patch and make an entry
			if not amrLayer[i+1] then
			
				-- next: find a new unused leaf
				-- for now, just this one node
				local leafIndex = 0
				
				if not self.amrLeafs[leafIndex+1] then
			
					print('splitting root node '..tostring(i)..' and putting in leaf node '..leafIndex)
		
					-- create info about the leaf
					self.amrLeafs[leafIndex+1] = {
						level = 0,	-- root
						layer = amrLayer,
						layerX = nx,		-- node x and y in the root
						layerY = ny,
						leafIndex = leafIndex,	-- which leaf we are using
					}
			
					-- tell the root layer table which node is used
					--  by pointing it back to the table of the leaf nodes 
					amrLayer[i+1] = self.amrLeafs[1]
				
					-- copy data from the root node location into the new node
					-- upsample as we go ... by nearest?
				
					-- TODO setup kernel args
					self.initNodeFromRootKernelObj.obj:setArg(1, ffi.new('int[4]', {nx, ny, 0, 0}))
					self.initNodeFromRootKernelObj.obj:setArg(2, ffi.new('int[1]', 0))
					self.app.cmds:enqueueNDRangeKernel{
						kernel = self.initNodeFromRootKernelObj.obj,
						dim = self.dim, 
						globalSize = self.amrNodeSizeWithoutBorder:ptr(),
						localSize = self.amrNodeSizeWithoutBorder:ptr(),
					}
				end
			end
		end
	end
end

os.exit()
	
	self.t = self.t + dt
	self.dt = dt
else
	self.t = self.t + dt
	self.dt = dt
end

end

function GridSolver:step(dt)
	self.integrator:integrate(dt, function(derivBuf)
		self:calcDeriv(derivBuf, dt)
	end)
	
	if self.eqn.useConstrainU then
		self:boundary()
		self.constrainUKernelObj()
	end

	for _,op in ipairs(self.ops) do
		if op.step then
			op:step(dt)
		end
	end
end

function GridSolver:printBuf(buf, ptr)
	ptr = ptr or buf:toCPU()
	local max = #tostring(self.numCells-1)
	for i=0,self.numCells-1 do
		io.write((' '):rep(max-#tostring(i)), i,':')
		for j=0,self.eqn.numStates-1 do
			io.write(' ', ptr[j + self.eqn.numStates * i])
		end 
		print()
	end
end

-- check for nans
-- expects buf to be of type cons_t, made up of numStates real variables
function GridSolver:checkFinite(buf)
	local ptr = buf:toCPU()
	local found
	for i=0,buf.size-1 do
		if not math.isfinite(ptr[i]) then
			found = found or table()
			found:insert(i)
		end
	end
	if not found then return end
--	self:printBuf(nil, ptr)
--	print(found:map(tostring):concat', ')
--	error'found non-finite numbers'
	return true
end

function GridSolver:getTex(var) 
	return self.tex
end

function GridSolver:calcDisplayVarToTex(var)
	local app = self.app
	local displayVarGroup = var.displayVarGroup
	if app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		app.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	
		var:setToTexArgs()
		var.calcDisplayVarToTexKernelObj()
		
		app.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		app.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		-- calcDisplayVarToTexPtr is sized without ghost cells 
		-- so is the GL texture
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex

		local channels = var.vectorField and 3 or 1
		local format = var.vectorField and gl.GL_RGB or gl.GL_RED

		var:setToBufferArgs()
		var.calcDisplayVarToBufferKernelObj()
		app.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * self.numCells * channels, ptr=ptr}
		local destPtr = ptr
		if app.is64bit then
			-- can this run in place?
			destPtr = ffi.cast('float*', ptr)
			for i=0,self.numCells*channels-1 do
				destPtr[i] = ptr[i]
			end
		end
		tex:bind()
		if self.dim < 3 then
			gl.glTexSubImage2D(gl.GL_TEXTURE_2D, 0, 0, 0, tex.width, tex.height, format, gl.GL_FLOAT, destPtr)
		else
			for z=0,tex.depth-1 do
				gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, tex.width, tex.height, 1, format, gl.GL_FLOAT, destPtr + channels * tex.width * tex.height * z)
			end
		end
		tex:unbind()
		glreport'here'
	end
end

function GridSolver:updateGUIParams()	
	if ig.igCollapsingHeader'parameters:' then
		GridSolver.super.updateGUIParams(self)
		
		for i=1,self.dim do
			for _,minmax in ipairs(minmaxs) do
				local k = xNames[i]..minmax
				if tooltip.numberTable(k, self[minmax..'s'], i, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
					-- if the domain changes
					-- then the dx has to change
					-- and all the stuff based on codePrefix has to change
					self:refreshCodePrefix()
				end
			end
		end
		
		for i=1,self.dim do
			for _,minmax in ipairs(minmaxs) do
				local var = xNames[i]..minmax
				if tooltip.comboTable(var, self.boundaryMethods, var, self.boundaryOptionNames) then
					self:refreshBoundaryProgram()
				end
			end
		end
	end
end

function GridSolver:updateGUIEqnSpecific()
	if ig.igCollapsingHeader'equation:' then
		if tooltip.comboTable('init state', self, 'initStateIndex', self.eqn.initStateNames) then
			-- TODO hmm ... the whole point of making a separate initStateProgram was to be able to refresh it without rebuilding all of the solver ...
			self:refreshEqnInitState()
		end	
		
		for _,var in ipairs(self.eqn.guiVars) do
			var:updateGUI(self)
		end
	end
end

do
	-- display vars: TODO graph vars
	local function handle(var, title)
		ig.igPushIDStr(title)
		
		local enableChanged = tooltip.checkboxTable('enabled', var, 'enabled') 
		ig.igSameLine()
		
		tooltip.checkboxTable('log', var, 'useLog')
		ig.igSameLine()

		tooltip.checkboxTable('fixed range', var, 'heatMapFixedRange')
		ig.igSameLine()
		
		if ig.igCollapsingHeader(var.name) then
			tooltip.numberTable('value min', var, 'heatMapValueMin')
			tooltip.numberTable('value max', var, 'heatMapValueMax')
		end
		
		ig.igPopID()
	
		return enableChanged
	end

	-- do one for 'all'
	local function _and(a,b) return a and b end
	local fields = {'enabled', 'useLog', 'heatMapFixedRange', 'heatMapValueMin', 'heatMapValueMax'}
	local defaults = {true, true, true, math.huge, -math.huge}
	local combines = {_and, _and, _and, math.min, math.max}
	local all = {name='all'}
	for i=1,#fields do
		all[fields[i]] = defaults[i]
	end
	local original = {}
	for i,field in ipairs(fields) do
		original[field] = defaults[i]
	end

	function GridSolver:updateGUIDisplay()
		local refresh 
		if ig.igCollapsingHeader'display:' then
			for i,displayVarGroup in ipairs(self.displayVarGroups) do
				ig.igPushIDStr('display '..i)
				if ig.igCollapsingHeader(displayVarGroup.name) then				
					for i=1,#fields do
						all[fields[i]] = defaults[i]
					end
					for _,var in ipairs(displayVarGroup.vars) do
						for i,field in ipairs(fields) do
							all[field] = combines[i](all[field], var[field])
						end
					end
					for _,field in ipairs(fields) do
						original[field] = all[field]
					end
					local enableChanged = handle(all, 'all')
					--refresh = refresh or enableChanged
					for _,field in ipairs(fields) do
						if all[field] ~= original[field] then
							for _,var in ipairs(displayVarGroup.vars) do
								var[field] = all[field]
							end
						end
					end

					for _,var in ipairs(displayVarGroup.vars) do
						local enableChanged = handle(var, displayVarGroup.name..' '..var.name)
						--refresh = refresh or enableChanged
					end
				end
				ig.igPopID()
			end
		end
		if refresh then
			-- solver and display are now tied together
			self:refreshSolverProgram()
		end
	end
end

function GridSolver:updateGUI()
	self:updateGUIParams()
	self:updateGUIEqnSpecific()
	self:updateGUIDisplay()

	-- heat map var

	-- TODO volumetric var
end

local Image = require 'image'
function GridSolver:save(prefix)
	-- TODO add planes to image, then have the FITS module use planes and not channels
	-- so the dimension layout of the buffer is [channels][width][height][planes]
	local width = tonumber(self.gridSize.x)
	local height = tonumber(self.gridSize.y)
	local depth = tonumber(self.gridSize.z)

	for _,bufferInfo in ipairs(self.buffers) do
		local name = bufferInfo.name
		local channels = bufferInfo.size / self.numCells / ffi.sizeof(self.app.real)
		if channels ~= math.floor(channels) then
			print("can't save buffer "..name.." due to its size not being divisible by the numCells")
		else
			local buffer = self[name]

			local numReals = self.numCells * channels
--[[ 3D interleave the depth and the channels ... causes FV to break 
			local image = Image(width, height, depth * channels, assert(self.app.real))
			local function getIndex(ch,i,j,k) return i + width * (j + height * (ch + channels * k)) end
--]]
-- [[ 3D interleave the depth and the width
			local image = Image(width * depth, height, channels, assert(self.app.real))
			local function getIndex(ch,i,j,k) return i + width * (k + depth * (j + height * ch)) end
--]]
			self.app.cmds:enqueueReadBuffer{buffer=buffer, block=true, size=ffi.sizeof(self.app.real) * numReals, ptr=image.buffer}
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
		
			local filename = prefix..'_'..name..'.fits'
			print('saving '..filename)
			image:save(filename)
		end
	end
end

return GridSolver
