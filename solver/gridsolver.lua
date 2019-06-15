local ffi = require 'ffi'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local io = require 'ext.io'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local math = require 'ext.math'
local string = require 'ext.string'
local vec3 = require 'vec.vec3'
local glreport = require 'gl.report'
local template = require 'template'
local vec3sz = require 'ffi.vec.vec3sz'
local tooltip = require 'tooltip'
local roundup = require 'util.roundup'
local time, getTime = table.unpack(require 'util.time')
local SolverBase = require 'solver.solverbase'
local makestruct = require'eqn.makestruct'

local common = require 'common'()
local minmaxs = common.minmaxs
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


-- whether to cache the opencl binaries
local useCache = true

local GridSolver = class(SolverBase)

GridSolver.numGhost = 2

GridSolver.solverVars  = table(GridSolver.super.solverVars):append{
	{name='gridSize', type='int4'},
	{name='stepsize', type='int4'},
	{name='grid_dx', type='real3'},
}

--[[
args:
	gridSize
	mins
	maxs
--]]
function GridSolver:initL1(args)
	GridSolver.super.initL1(self, args)

	-- same as equations
	-- but let equations/init conds add to the solver vars (as gui vars)
	-- then we can edit them without recompiling the kernels
	self.solverVars = table(self.solverVars)

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
end

--[[
args:
	boundary = boundary info
--]]
function GridSolver:preInit(args)
	GridSolver.super.preInit(self, args)

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

	-- my kernel objs are going to need workgroup info based on domain.size-2*noGhost as well as domain.size ... 
	-- ...and rather than require an extra argument, I think I'll just take advantage of a closure
	local solver = self
	local Program = class(require 'cl.obj.program')
	if ffi.os == 'Windows' then
		os.execute'mkdir cache-cl 2> nul'
	else
		os.execute'mkdir cache-cl 2> /dev/null'
	end
	function Program:init(args)
		self.name = args.name
		args.env = solver.app.env
		args.domain = solver.domain
		-- [[ caching binaries, which doesn't write unless the program successfully compiles
		if args.name then
			local path = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
			if useCache then args.cacheFile = path end
		end
		Program.super.init(self, args)
		--]]
		--[[ Write generated code the first time.  Subsequent times use the pre-existing code.  Useful for debugging things in the generated OpenCL.
		local path = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
		local clfn = path..'.cl'
		if io.fileexists(clfn) then
			local cachedCode = file[clfn]
			assert(cachedCode:sub(1,#args.env.code) == args.env.code, "seems you have changed the cl env code")
			args.code = cachedCode:sub(#args.env.code+1)
			Program.super.init(self, args)
		else
			Program.super.init(self, args)	-- do this so getCode() works
			file[clfn] = self:getCode()
		end
		--]]
	end
	function Program:compile(args)
		args = args or {}
		args.buildOptions = '-w'	-- show warnings
		local results = Program.super.compile(self, args)
		if self.obj then	-- did compile
			print((self.name and self.name..' ' or '')..'log:')
			print(string.trim(self.obj:getLog(self.env.device)))
		end
		return results
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
	if initState.mins then 
		self.mins = vec3(table.unpack(initState.mins)) 
		for j=1,3 do
			self.solverPtr.mins.s[j-1] = self.mins[j]
		end
	end
	if initState.maxs then 
		self.maxs = vec3(table.unpack(initState.maxs)) 
		for j=1,3 do
			self.solverPtr.maxs.s[j-1] = self.maxs[j]
		end
	end

	-- there's a lot of overlap between this and the solverBuf creation... 
	self:refreshSolverBufMinsMaxs()
	
	-- while we're here, write all gui vars to the solver_t
	for _,var in ipairs(self.eqn.guiVars) do
		if not var.compileTime then
			var:setToSolver(self)
		end
	end
	
	self:refreshSolverBuf()

--[[
local ptr = ffi.cast(self.eqn.cons_t..'*', self.UBufObj:toCPU())
local t = self.t
for i=1,tonumber(self.gridSize.x) do
	local x = self.solverPtr.mins.x + self.solverPtr.grid_dx.x * (i - self.numGhost - .5)
	local rho = 1 + .32 * math.sin(2 * math.pi * (x - t))
	local vx = 1
	local P = 1
	local mx = rho * vx
	local EKin = .5 * rho * vx * vx
	local EInt = P / (self.solverPtr.heatCapacityRatio - 1)
	local ETotal = EKin + EInt
	ptr[i-1].rho = rho
	ptr[i-1].m.x = mx
	ptr[i-1].m.y = 0
	ptr[i-1].m.z = 0
	ptr[i-1].ETotal = ETotal
	ptr[i-1].ePot = 0
end
self.UBufObj:fromCPU(ptr)
--]]

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

	local graphShaderCode = file['draw/graph.shader']
	self.graphShader = GLProgram{
		vertexCode = template(graphShaderCode, {
			solver = self,
			vertexShader = true,
		}),
		fragmentCode = template(graphShaderCode, {
			solver = self,
			fragmentShader = true,
		}),
		uniforms = {
			tex = 0,
			scale = 1,
			ambient = 1,
		},
	}

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
		
		self.display3D_Ray_maxiter = math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y), tonumber(self.gridSize.z))
		local volumetricCode = file['draw/volumetric.shader']
		self.volumeRayShader = GLProgram{
			vertexCode = template(volumetricCode, {
				app = self.app,
				solver = self,
				vertexShader = true,
			}),
			fragmentCode = template(volumetricCode, {
				app = self.app,
				solver = self,
				fragmentShader = true,
			}),
			uniforms = {
				tex = 0,
				gradientTex = 1,
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

function GridSolver:createSolverBuf()
	GridSolver.super.createSolverBuf(self)
	
	-- do this before any call to createBuffers or createCodePrefix
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
			var:setToSolver(self)
		end
	end

	self:refreshSolverBuf()
end

function GridSolver:refreshSolverBufMinsMaxs()
	-- always set this to the full range, even outside the used dimension, in case some dimensional value is supposed to be a non-zero constant, esp for cell volume calculations
	self.solverPtr.mins.x = self.mins[1]
	self.solverPtr.mins.y = self.mins[2]
	self.solverPtr.mins.z = self.mins[3]
	self.solverPtr.maxs.x = self.maxs[1]
	self.solverPtr.maxs.y = self.maxs[2]
	self.solverPtr.maxs.z = self.maxs[3]
	self.solverPtr.grid_dx.x = (self.solverPtr.maxs.x - self.solverPtr.mins.x) / tonumber(self.sizeWithoutBorder.x)
	self.solverPtr.grid_dx.y = (self.solverPtr.maxs.y - self.solverPtr.mins.y) / tonumber(self.sizeWithoutBorder.y)
	self.solverPtr.grid_dx.z = (self.solverPtr.maxs.z - self.solverPtr.mins.z) / tonumber(self.sizeWithoutBorder.z)
	if self.app.verbose then
		print('mins = '..self.solverPtr.mins.x..', '..self.solverPtr.mins.y..', '..self.solverPtr.mins.z)
		print('maxs = '..self.solverPtr.maxs.x..', '..self.solverPtr.maxs.y..', '..self.solverPtr.maxs.z)
		print('grid_dx = '..self.solverPtr.grid_dx.x..', '..self.solverPtr.grid_dx.y..', '..self.solverPtr.grid_dx.z)
	end
end

-- TODO some of this is copied in solverbase
function GridSolver:createBuffers()
	local app = self.app
	local realSize = ffi.sizeof(app.real)

	-- to get sizeof
	makestruct.safeFFICDef(self:getConsLRTypeCode())

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
	-- or finally make use of constant args ...

	-- this much for the first level U data
	self:clalloc('UBuf', self.eqn.cons_t, self.numCells)

	if self.usePLM then
		self:clalloc('ULRBuf', self.eqn.consLR_t, self.numCells * self.dim)
	end

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar 
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', self.app.real, self.numCells * 3)
	self:clalloc('reduceSwapBuf', self.app.real, math.ceil(self.numCells / self.localSize1d))
	self.reduceResultPtr = ffi.new('real[1]', 0)


	if self.allowAccum then
		-- as big as reduceBuf, because it is a replacement for reduceBuf
		-- ... though I don't accum on vector fields yet, so it doesn't need the x3 really
		self:clalloc('accumBuf', self.app.real, self.numCells * 3)
	end

	if self.app.targetSystem ~= 'console' then
		-- CL/GL interop

		-- hmm, notice I'm still keeping the numGhost border on my texture 
		-- if I remove the border altogether then I get wrap-around
		-- maybe I should just keep a border of 1?
		-- for now i'll leave it as it is
		local GLTex2D = require 'gl.tex2d'
		local GLTex3D = require 'gl.tex3d'
		local cl = self.dim < 3 and GLTex2D or GLTex3D
		-- TODO check for extension GL_ARB_half_float_pixel
		local gltype = app.real == 'half' and gl.GL_HALF_FLOAT_ARB or gl.GL_FLOAT
		self.tex = cl{
			width = tonumber(self.gridSize.x),
			height = tonumber(self.gridSize.y),
			depth = tonumber(self.gridSize.z),
			internalFormat = gl.GL_RGBA32F,
			format = gl.GL_RGBA,
			type = gltype,
			minFilter = gl.GL_NEAREST,
		
			-- TODO toggle this in the gui:
			--magFilter = gl.GL_NEAREST,
			magFilter = gl.GL_LINEAR,
			
			wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
		}

		local CLImageGL = require 'cl.imagegl'
		if app.useGLSharing then
			self.texCLMem = CLImageGL{context=app.ctx, tex=self.tex, write=true}
		else
			self.calcDisplayVarToTexPtr = ffi.new(app.realparam..'[?]', self.numCells * 3)
			
			--[[ PBOs?
			self.calcDisplayVarToTexPBO = ffi.new('gl_int[1]', 0)
			gl.glGenBuffers(1, self.calcDisplayVarToTexPBO)
			gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, self.calcDisplayVarToTexPBO[0])
			gl.glBufferData(gl.GL_PIXEL_UNPACK_BUFFER, self.tex.width * self.tex.height * ffi.sizeof(app.real) * 4, nil, gl.GL_STREAM_READ)
			gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
			--]]
		end
	end
end

function GridSolver:createCodePrefix()
	local lines = table()
	
	GridSolver.super.createCodePrefix(self)
	
	lines:insert(self.codePrefix)

	lines:append{
		'#define numGhost '..self.numGhost,
	}:append{
		'#define INDEX(a,b,c)	((a) + solver->gridSize.x * ((b) + solver->gridSize.y * (c)))',
		'#define INDEXV(i)		indexForInt4ForSize(i, solver->gridSize.x, solver->gridSize.y, solver->gridSize.z)',
	
	--[[
	naming conventions ...
	* the grid indexes i_1..i_n that span 1 through solver->gridSize.1..solver->gridSize.n
	(between the index and the coordinate space:)
		- grid_dx? is the change in coordinate space wrt the change in grid space
		- cell_x(i) calculates the coordinates at index i
	* the coordinates space x_1..x_n that spans mins.s1..maxs.sn
	(between the coordinate and the embedded space:)
		- vectors can be calculated from Cartesian by cartesianToCoord
		- the length of the basis vectors wrt the change in indexes is given by cell_dx?(x)
		- the Cartesian length of the holonomic basis vectors is given by coord_dx?(x).  
			This is like cell_dx? except not scaled by grid_dx?
			This is just the change in embedded wrt the change in coordinate, not wrt the change in grid
		- cell_volume(x) gives the volume between indexes at the coordinate x
		- the Cartesian length of a vector in coordinate space is given by coordLen and coordLenSq
	* the embedded Cartesian space ... idk what letters I should use for this.  
		Some literature uses x^I vs coordinate space x^a, but that's still an 'x', no good for programming.
		Maybe I'll use 'xc' vs 'x', like I've already started to do in the initial conditions.
	--]]
	
	}:append(range(self.dim):map(function(i)
	-- mapping from index to coordinate 
		return (('#define cell_x{i}(i) (((real)(i) + .5 - (real)numGhost) * solver->grid_dx.{x} + solver->mins.{x})')
			:gsub('{i}', i-1)
			:gsub('{x}', xNames[i])
		)
	end)):append(range(self.dim+1,3):map(function(i)
	-- non-dimension coordinates don't need ghosts (right)
		return (('#define cell_x{i}(i) (((real)(i) + .5) * solver->grid_dx.{x} + solver->mins.{x})')
			:gsub('{i}', i-1)
			:gsub('{x}', xNames[i])
		)
	end)):append{
		'#define cell_x(i) _real3(cell_x0(i.x), cell_x1(i.y), cell_x2(i.z))',
	
		-- bounds-check macro
		'#define OOB(lhs,rhs) (i.x < (lhs) || i.x >= solver->gridSize.x - (rhs)'
			.. (self.dim < 2 and '' or ' || i.y < (lhs) || i.y >= solver->gridSize.y - (rhs)')
			.. (self.dim < 3 and '' or ' || i.z < (lhs) || i.z >= solver->gridSize.z - (rhs)')
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
		return i <= solver.dim and 'numGhost' or '0' 
	end):concat','?>); \
	int index = INDEXV(i);
]], {solver = self}),
	}

	
	-- volume of a cell = volume element times grid dx's 
	lines:insert(template([[
static inline real cell_volume(constant <?=solver.solver_t?>* solver, real3 x) {
	return coord_volume(x)<?
for i=1,solver.dim do
?> * solver->grid_dx.<?=xNames[i]?><?
end
?>;
}
]], {
	solver = self,
	xNames = xNames,
}))

	lines:append{
		-- not messing with this one yet
		self.allocateOneBigStructure and '#define allocateOneBigStructure' or '',
		
		self:getConsLRTypeCode(),
	}

	self.codePrefix = lines:concat'\n'

--print'done building solver.codePrefix'
--print(self.codePrefix)
end

function GridSolver:resetState()
	self.app.cmds:finish()
		
	-- start off by filling all buffers with zero, just in case initState forgets something ...
	for _,bufferInfo in ipairs(self.buffers) do
		self.app.cmds:enqueueFillBuffer{buffer=self[bufferInfo.name], size=bufferInfo.count * ffi.sizeof(bufferInfo.type)}
	end

	GridSolver.super.resetState(self)

	self:boundary()
	
	if self.eqn.useConstrainU then
		self.constrainUKernelObj(self.solverBuf, self.UBuf)
		self:boundary()
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
end

function GridSolver:getSolverCode()
	local slopeLimiterCode = 'real slopeLimiter(real r) {\n'
		.. '\t'..self.app.limiters[self.slopeLimiter].code..'\n'
		.. '}\n'
	
	return table{
		GridSolver.super.getSolverCode(self),
			
		slopeLimiterCode,
		
		self.usePLM and template(file['solver/plm.cl'], {solver=self, eqn=self.eqn}) or '',
		self.useCTU and template(file['solver/ctu.cl'], {solver=self, eqn=self.eqn}) or '',
	}:concat'\n'
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
	self.getULRBufType = self.usePLM and self.eqn.consLR_t or self.eqn.cons_t
	self.getULRBufName = self.usePLM and 'ULRBuf' or 'UBuf'

	self.getULRArg = self.getULRBufType..'* '..self.getULRBufName

	-- this code creates the const global cons_t* UL, UR variables
	-- it assumes that indexL, indexR, and side are already defined
	-- both UBuf and ULRBuf should be cell-centered
	if self.usePLM then
		self.getULRCode = function(self, args)
			args = args or {}
			local suffix = args.suffix or ''
			return template([[
	const global <?=eqn.cons_t?>* UL<?=suffix?> = &<?=solver.getULRBufName?>[side + dim * <?=indexL?>].R;
	const global <?=eqn.cons_t?>* UR<?=suffix?> = &<?=solver.getULRBufName?>[side + dim * <?=indexR?>].L;
]],			{
				solver = self,
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
	const global <?=eqn.cons_t?>* UL<?=suffix?> = <?=bufName?> + <?=indexL?>;
	const global <?=eqn.cons_t?>* UR<?=suffix?> = <?=bufName?> + <?=indexR?>;
]],			{
				solver = self,
				eqn = self.eqn,
				suffix = suffix,
				indexL = args.indexL or 'indexL'..suffix,
				indexR = args.indexR or 'indexR'..suffix,
				bufName = args.bufName or self.getULRBufName,	-- for displayVars the variable name is 'buf', so I need to override it either in displayCode or here
			})
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
		-- TODO constant/'dirichlet' conditions
		-- and while you're at it, let boundary options register GUI controls, so we can treat the constant as a GUI-editable parameter 
		
		-- this is purely for debugging.  annnd it doesn't crash, whereas rectangular grids are crashing...
		{none = function(args) return ''  end},
		
		{periodic = function(args)
			local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
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
			local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
			local lines = table()
			if args.minmax == 'min' then
				lines:insert(indent..args.assign(
					args.array('buf', args.index'j'),
					args.array('buf', args.index'2*numGhost-1-j'))..';')
				if args.mirrorVars and args.mirrorVars[args.side] then
					lines:append(table.map(args.mirrorVars[args.side], function(var)
						return indent..args.assign(
							args.field(args.array('buf', args.index'j'), var),
							args.field('-'..args.array('buf', args.index'j'), var)
						)..';'
					end))
				end
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				lines:insert(indent..args.assign(
					args.array('buf', args.index(rhs)),
					args.array('buf', args.index(gridSizeSide..'-numGhost-1-j')))..';')
				if args.mirrorVars and args.mirrorVars[args.side] then
					lines:append(table.map(args.mirrorVars[args.side], function(var)
						return indent..args.assign(
							args.field(args.array('buf', args.index(rhs)), var),
							args.field('-'..args.array('buf', args.index(rhs)), var)
						)..';'
					end))
				end
			end
			return lines:concat'\n'
		end},
		{freeflow = function(args)
			local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
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
		-- constant-derivative / linear extrapolation
		{linear = function(args)
			local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
			if args.minmax == 'min' then
				return indent..'for (int k = 0; k < numStates; ++k) {\n'
					..indent..'\t'..args.assign(
						'buf['..args.index'j'..'].ptr[k]',
						'(buf['..args.index'numGhost'..'].ptr[k] - buf['..args.index'numGhost+1'..'].ptr[k]) * (numGhost - j+1) + buf['..args.index'numGhost+1'..'].ptr[k]'
					)..';\n'..indent..'}'
			elseif args.minmax == 'max' then
				local rhs = gridSizeSide..'-numGhost+j'
				return indent..'for (int k = 0; k < numStates; ++k) {\n'
					..indent..'\t'..args.assign(
						'buf['..args.index(rhs)..'].ptr[k]',
						'(buf['..args.index(gridSizeSide..'-numGhost-1')..'].ptr[k] - buf['..args.index(gridSizeSide..'-numGhost-2')..'].ptr[k]) * (j+1) + buf['..args.index(gridSizeSide..'-numGhost-2')..'].ptr[k]'
					)..';\n'..indent..'}'
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
	
		lines:insert(template([[
kernel void boundary_<?=xNames[side]?>(
	constant <?=solver.solver_t?>* solver,
	global <?=args.type?>* buf<?= 
args.extraArgs and #args.extraArgs > 0 
	and ','..table.concat(args.extraArgs, ',\n\t')
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
			solver = self,
			args = args,
			side = side, 
			xNames = xNames,
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
	for (int j = 0; j < numGhost; ++j) {]]

		for _,minmax in ipairs(minmaxs) do
			lines:insert(args.methods[xNames[side]..minmax]{
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

lines:insert[[
	}
}
]]
	end

	local code = lines:concat'\n'

	local boundaryProgramObj
	time('building boundary program', function()
		boundaryProgramObj = self.Program{name='boundary', code=code}
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
		obj.obj:setArg(1, self.UBuf)
	end
	for _,op in ipairs(self.ops) do
		if op.refreshBoundaryProgram then
			op:refreshBoundaryProgram()
		end
	end

	if self.useCTU and self.usePLM then
		local args = self:getBoundaryProgramArgs()
		args.type = self.eqn.consLR_t..'_dim'
		if args.mirrorVars then
			local newvars = {}
			for i=1,3 do
				newvars[i] = table()
				for _,var in ipairs(args.mirrorVars[i]) do
					for j=0,self.dim-1 do
						newvars[i]:insert('side['..j..'].L.'..var)
						newvars[i]:insert('side['..j..'].R.'..var)
					end
				end
			end
			args.mirrorVars = newvars
		end
		self.lrBoundaryProgramObj, self.lrBoundaryKernelObjs = self:createBoundaryProgramAndKernel(args)
		for _,obj in ipairs(self.lrBoundaryKernelObjs) do
			obj.obj:setArg(1, self:getULRBuf())
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
	numStates = numStates or self.eqn.numIntStates
	local exact = assert(self.eqn.initState.exactSolution, "can't test accuracy of a configuration that has no exact solution")
	local ptr = ffi.cast(self.eqn.cons_t..'*', self.UBufObj:toCPU())
	assert(self.dim == 1)
	local n = self.gridSize.x
	local ghost = self.numGhost
	local err = 0
	for i=ghost+1,tonumber(self.gridSize.x)-2*ghost do
		--local x = self.xs[i+ghost]
		local x = self.solverPtr.mins.x + self.solverPtr.grid_dx.x * (i - ghost - .5)
		err = err + compareL1(ptr[i-1].ptr, numStates, exact(self, x, self.t))
	end
	err = err / (numStates * (tonumber(self.gridSize.x) - 2 * ghost))
	return err, ptr
end

function GridSolver:update()
--[[
print'\nGridSolver:update() begin, self.UBufObj:'
self:printBuf(self.UBufObj)
--]]

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	GridSolver.super.update(self)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	local dt = self:calcDT()

	-- first do a step
	self:step(dt)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	-- why was this moved out of :step() ?
	self.t = self.t + dt
	self.dt = dt

	if cmdline.testAccuracy then
		local err = self:calcExactError()
		if #self.app.solvers > 1 then
			io.write(self.name,'\t')
		end
		print('t='..self.t..' L1-error='..err)
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	self:boundary()

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
end

function GridSolver:step(dt)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	self.integrator:integrate(dt, function(derivBufObj)		
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		if self.eqn.useConstrainU then
			self:boundary()
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
			self.constrainUKernelObj(self.solverBuf, self.UBuf)
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
		end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
		
		self:boundary()	
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		self:calcDeriv(derivBufObj, dt)
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
		if self.eqn.useSourceTerm then
			self:boundary()
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
			self.addSourceKernelObj(self.solverBuf, derivBufObj.obj, self.UBuf)
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

		for _,op in ipairs(self.ops) do
			if op.addSource then
				op:addSource(derivBufObj)
			end
		end

	end)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

	if self.eqn.useConstrainU then
		self:boundary()
		self.constrainUKernelObj(self.solverBuf, self.UBuf)
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
	
	for _,op in ipairs(self.ops) do
		if op.step then
			self:boundary()
			op:step(dt)
		end
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

end

function GridSolver:getTex(var) 
	return self.tex
end
function GridSolver:getHeatMap2DShader(var)
	return self.heatMap2DShader
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
	
		self:calcDisplayVarToBuffer(var)

		local sizevec = var.getBuffer().sizevec or self.gridSize
		assert(sizevec.x <= tex.width and sizevec.y <= tex.height)
		local volume = tonumber(sizevec:volume())
		
		app.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * volume * channels, ptr=ptr}
		local destPtr = ptr
		-- TODO check for extension GL_ARB_half_float_pixel
		local gltype = app.real == 'half' and gl.GL_HALF_FLOAT_ARB or gl.GL_FLOAT
		
		if app.real == 'double' then
			-- can this run in place?  seems like it
			destPtr = ffi.cast('float*', ptr)
			for i=0,volume*channels-1 do
				destPtr[i] = ptr[i]
			end
		end
--[[ wtf...
for i=0,volume*channels-1 do
destPtr[i] = i/(volume*channels-1)
end
--]]		
		tex:bind()
		if self.dim < 3 then
			gl.glTexSubImage2D(gl.GL_TEXTURE_2D, 0, 0, 0, sizevec.x, sizevec.y, format, gltype, destPtr)
		else
			for z=0,tex.depth-1 do
				gl.glTexSubImage3D(gl.GL_TEXTURE_3D, 0, 0, 0, z, sizevec.x, sizevec.y, 1, format, gltype, destPtr + channels * sizevec.x * sizevec.y * z)
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

		var.enabled = not not var.enabled
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

	local function search(str, what)
		return str:match(what:gsub('%s+', '.*'))
	end

	function GridSolver:updateGUIDisplay()
		if self.guiDisplayFilterStr == nil then self.guiDisplayFilterStr = '' end
		if self.guiDisplayFilterEnabledVars == nil then self.guiDisplayFilterEnabledVars = false end
		local refresh 
		if ig.igCollapsingHeader'display:' then
			tooltip.textTable('filter', self, 'guiDisplayFilterStr')
			ig.igSameLine()
			tooltip.checkboxTable('filter enabled', self, 'guiDisplayFilterEnabledVars')
			
			for i,displayVarGroup in ipairs(self.displayVarGroups) do
				ig.igPushIDStr('display '..i)
				if ig.igCollapsingHeader(displayVarGroup.name) then
					for i=1,#fields do
						all[fields[i]] = defaults[i]
					end
					for _,var in ipairs(displayVarGroup.vars) do
						local title = displayVarGroup.name..' '..var.name
						if (#self.guiDisplayFilterStr == 0 or search(title, self.guiDisplayFilterStr))
						and (not self.guiDisplayFilterEnabledVars or var.enabled)
						then
							for i,field in ipairs(fields) do
								all[field] = combines[i](all[field], var[field])
							end
						end
					end
					for _,field in ipairs(fields) do
						original[field] = all[field]
					end
					local enableChanged = handle(all, 'all')
					--refresh = refresh or enableChanged
					if enableChanged then
						for _,field in ipairs(fields) do
							if all[field] ~= original[field] then
								for _,var in ipairs(displayVarGroup.vars) do
									local title = displayVarGroup.name..' '..var.name
									if (#self.guiDisplayFilterStr == 0 or search(title, self.guiDisplayFilterStr))
									and (not self.guiDisplayFilterEnabledVars or var.enabled)
									then
										var[field] = all[field]
									end
								end
							end
						end
					end

					for _,var in ipairs(displayVarGroup.vars) do
						local title = displayVarGroup.name..' '..var.name
						if (#self.guiDisplayFilterStr == 0 or search(title, self.guiDisplayFilterStr))
						and (not self.guiDisplayFilterEnabledVars or var.enabled)
						then
							local enableChanged = handle(var, title)
							--refresh = refresh or enableChanged
						end
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
		print("can't save buffer due to its size not being divisible by the numCells")
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
