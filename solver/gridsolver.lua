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
local glreport = require 'gl.report'
local clnumber = require 'cl.obj.number'
local template = require 'template'
local vec3d = require 'vec-ffi.vec3d'
local vec3sz = require 'vec-ffi.vec3sz'
local tooltip = require 'tooltip'
local roundup = require 'util.roundup'
local time, getTime = table.unpack(require 'util.time')
local SolverBase = require 'solver.solverbase'
local makestruct = require'eqn.makestruct'

local half = require 'half'
local toreal, fromreal = half.toreal, half.fromreal


local common = require 'common'
local minmaxs = common.minmaxs
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


-- whether to cache the opencl binaries
local useCache = cmdline.useCache
if useCache == nil then useCache = true end

local function unpack(x)
	if x.unpack then return x:unpack() end
	return table.unpack(x)
end


local GridSolver = class(SolverBase)

GridSolver.numGhost = 2

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
	
	self.solverStruct.vars:append{
		{name='grid_dx', type='real3'},
		{name='gridSize', type='int4'},
		{name='stepsize', type='int4'},
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
--]]
function GridSolver:preInit(args)
	GridSolver.super.preInit(self, args)
	
	self:createBoundaryOptions()

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
		
		-- caching binaries, which doesn't write unless the program successfully compiles
		if not cmdline.usecachedcode then 
			if args.name then
				local path = 'cache-cl/'..solver.app:uniqueName(assert(args.name))
				if useCache then args.cacheFile = path end
			end
			Program.super.init(self, args)
		
		-- Write generated code the first time.  Subsequent times use the pre-existing code.  Useful for debugging things in the generated OpenCL.
		else
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
		end
	end
	function Program:compile(args)
		args = args or {}
		args.buildOptions = '-w'	-- show warnings
		local results = Program.super.compile(self, args)
		if self.obj then	-- did compile
			print((self.name and self.name..' ' or '')..'log:')
			-- TODO log per device ...			
			print(string.trim(self.obj:getLog(solver.device)))
		end
		return results
	end
	self.Program = Program

	
	local GLProgram = class(require 'gl.program')
	function GLProgram:init(...)
		local args = ...
		
		-- Write generated code
		local fn = 'cache-cl/'..solver.app:uniqueName(args.name or 'shader')
		file[fn..'.vert'] = args.vertexCode
		file[fn..'.frag'] = args.fragmentCode

		GLProgram.super.init(self, ...)
	end
	self.GLProgram = GLProgram 
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
		self.mins = vec3d(unpack(initState.mins)) 
		for j=1,3 do
			self.solverPtr.mins.s[j-1] = toreal(self.mins.s[j-1])
		end
	end
	if initState.maxs then 
		self.maxs = vec3d(unpack(initState.maxs)) 
		for j=1,3 do
			self.solverPtr.maxs.s[j-1] = toreal(self.maxs.s[j-1])
		end
	end

	-- there's a lot of overlap between this and the solverBuf creation... 
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
	local graphShaderCode = file['draw/graph.shader']
	self.graphShader = self.GLProgram{
		name = 'graph',
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
	self.heatMap2DShader = self.GLProgram{
		name = '2d_heatmap',
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

	do	--if self.dim == 3 then
		-- raytracing (stalling)
		
		self.display3D_Ray_maxiter = math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y), tonumber(self.gridSize.z))
		local volumetricCode = file['draw/volumetric.shader']
		self.volumeRayShader = self.GLProgram{
			name = 'volumetric',
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

		local volumeSliceCode = file['draw/3d_slice.shader']
		self.volumeSliceShader = self.GLProgram{
			name = '3d_slice',
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

	-- TODO move to draw/vector_arrow.lua ?
	local vectorArrowCode = file['draw/vector_arrow.shader']
	self.vectorArrowShader = self.GLProgram{
		name = 'vector_arrow',
		vertexCode = template(vectorArrowCode, {
			solver = self,
			vertexShader = true,
		}),
		fragmentCode = template(vectorArrowCode, {
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

	local vectorLICCode = file['draw/vector_lic.shader']
	self.vectorLICShader = self.GLProgram{
		name = 'vector_lic',
		vertexCode = template(vectorLICCode, {
			solver = self,
			vertexShader = true,
		}),
		fragmentCode = template(vectorLICCode, {
			fragmentShader = true,
			solver = self,
		}),
		uniforms = {
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
			noiseTex = 2,
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
			if var.ctype == 'real' then
				self.solverPtr[var.name] = toreal(var.value)
			else
				self.solverPtr[var.name] = var.value
			end
		end
	end

	self:refreshSolverBuf()
end

function GridSolver:refreshSolverBufMinsMaxs()
	-- always set this to the full range, even outside the used dimension, in case some dimensional value is supposed to be a non-zero constant, esp for cell volume calculations
	for j=1,3 do
		self.solverPtr.mins.s[j-1] = toreal(self.mins.s[j-1])
		self.solverPtr.maxs.s[j-1] = toreal(self.maxs.s[j-1])
		self.solverPtr.grid_dx.s[j-1] = toreal( (self.maxs.s[j-1] - self.mins.s[j-1]) / tonumber(self.sizeWithoutBorder.s[j-1]) )
		self.solverPtr.initCondMins.s[j-1] = toreal( self.initCondMins.s[j-1] )
		self.solverPtr.initCondMaxs.s[j-1] = toreal( self.initCondMaxs.s[j-1] )
	end
	if self.app.verbose then
		print('mins = '..fromreal(self.solverPtr.mins.x)..', '..fromreal(self.solverPtr.mins.y)..', '..fromreal(self.solverPtr.mins.z))
		print('maxs = '..fromreal(self.solverPtr.maxs.x)..', '..fromreal(self.solverPtr.maxs.y)..', '..fromreal(self.solverPtr.maxs.z))
		print('grid_dx = '..fromreal(self.solverPtr.grid_dx.x)..', '..fromreal(self.solverPtr.grid_dx.y)..', '..fromreal(self.solverPtr.grid_dx.z))
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
		-- TODO self.eqn.consLR_t..'_dim' and remove * self.dim ?
		self:clalloc('ULRBuf', self.eqn.consLR_t, self.numCells * self.dim)
	end

	-- used both by reduceMin and reduceMax
	-- (and TODO use this by sum() in implicit solver as well?)
	-- times three because this is also used by the displayVar 
	-- on non-GL-sharing cards.
	self:clalloc('reduceBuf', self.app.real, self.numCells * 3)
	self:clalloc('reduceSwapBuf', self.app.real, math.ceil(self.numCells / self.localSize1d))
	self.reduceResultPtr = ffi.new'real[1]'
	self.reduceResultPtr[0] = toreal(0)


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
			magFilter = gl.GL_LINEAR,
			wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
		}

		local CLImageGL = require 'cl.imagegl'
		if app.useGLSharing then
			self.texCLMem = CLImageGL{context=app.ctx, tex=self.tex, write=true}
		else
			self.calcDisplayVarToTexPtr = ffi.new(app.real..'[?]', self.numCells * 3)
			
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
		Some literature uses r^i or u^i vs coordinate space x^a.
		Maybe I'll use 'xc' vs 'x', like I've already started to do in the initial conditions.
	
	functionality (and abstraction):
		- Allow the option for precomputing certain variables of the coordinate system: the chart mapping.
		 (Though for some operations, like the Euler equations wavespeeds, it is faster to compute them than to retrieve them from memory, so I suspect similar for the coordinate chart.)
		- Add a macro for useful values to compute, like the x,y,z,r,theta,phi variables.
		 (This is already started in the coords.vars table.)
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

	-- above we have macros that are used in solver & eqn codePrefix
	lines:insert(self.codePrefix)
	-- below here we have solver_t usage
	
	-- volume of a cell = volume element times grid dx's 
	lines:insert(template([[
static inline real cell_sqrt_det_g(constant <?=solver.solver_t?>* solver, real3 x) {
	return coord_sqrt_det_g(x)<?
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
	self.cmds:finish()
		
	-- start off by filling all buffers with zero, just in case initState forgets something ...
	for _,bufferInfo in ipairs(self.buffers) do
		self.cmds:enqueueFillBuffer{buffer=self[bufferInfo.name], size=bufferInfo.count * ffi.sizeof(bufferInfo.type)}
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
		
-- TODO constant/'dirichlet' conditions
-- and while you're at it, let boundary options register GUI controls, so we can treat the constant as a GUI-editable parameter 

local Boundary = class()
GridSolver.Boundary = Boundary

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
		src = args.index([[numGhost + (j-numGhost + 2 * (]]..gridSizeSide..' - 2 * numGhost)) % ('..gridSizeSide..[[ - 2 * numGhost)]])
	elseif args.minmax == 'max' then
		dst = args.index(gridSizeSide..'-1-j')
		src = args.index('numGhost + (numGhost-1-j) % ('..gridSizeSide..' - 2 * numGhost)')
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
		src = args.index'2*numGhost-1-j'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..'-numGhost+j')
		src = args.index(gridSizeSide..'-numGhost-1-j')
	end
	local lines = table()
	lines:insert('\t\t'..self:assignDstSrc(dst, src, args))
	
	-- reflectVars.mirror is going to hold, for each dimension, which components need to be mirrored
	-- v.s[side-1] = -v.s[side-1]
	-- generalized:
	local solver = args.solver
	local eqn = solver.eqn
	if solver.coord.vectorComponent == 'cartesian' 
	and not require 'coord.cartesian'.is(solver.coord)
	then
		-- v = v - n (v dot n)
		-- v^i = v^i - n^i (v^j n_j) (1 + restitution)
		-- ... where n_i = partial_i u, for u the chart
		-- TODO let the ctor specify the vars, not this
		--  so I can use this for things like the poisson solver
		lines:insert(template([[
		{
			int4 iv = (int4)(<?=iv?>,0);
			real3 x = cell_x(iv);
<? if args.minmax == 'min' then ?>
			real3 n = coord_cartesianFromCoord(normalForSide<?=side-1?>, x);
<? else -- max ?>
			real3 n = coord_cartesianFromCoord(real3_neg(normalForSide<?=side-1?>), x);
<? end ?>
]], 	{
			args = args,
			iv = args.minmax == 'min' 
				and args.indexv'j'
				or args.indexv('solver->gridSize.'..xNames[args.side]..'-numGhost+j'),
			side = args.side,
		}))
		for _,var in ipairs(eqn.consStruct.vars) do
			if var.type == 'real3' 
			or var.type == 'cplx3'
			then
				lines:insert(template([[
			buf[<?=dst?>].<?=field?> = <?=vec3?>_sub(
				buf[<?=dst?>].<?=field?>,
				<?=vec3?>_<?=scalar?>_mul(
					<?=vec3?>_from_real3(n),
					<?=scalar?>_real_mul(
						<?=vec3?>_real3_dot(
							buf[<?=dst?>].<?=field?>,
							n
						), 
						<?=restitution + 1?>
					)
				)
			);
]], 			{
					restitution = clnumber(self.restitution),
					vec3 = var.type,
					scalar = var.type == 'cplx3' and 'cplx' or 'real',
					field = var.name,
					dst = dst,
				}))
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
-- TODO fixedCode needs to be provided ... but can't be provided easily from outside if it needs templated code from inside the class ...
-- for now lets just not include this in the options, and let child classes subclass it if they want
local BoundaryFixed = class(Boundary)
BoundaryFixed.name = 'fixed'
function BoundaryFixed:init(fixedCode)
	-- fixed values to use
	self.code = fixedCode or self.fixedCode
end
function BoundaryFixed:getCode(args)
	local fixedCode = template(self.code, {
		solver = args.solver,
		eqn = args.solver.eqn,
		args = args,
	})
	local dst
	if args.minmax == 'min' then
		dst = args.index'j'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..'-numGhost+j')
	end
	local lines = table()
	if args.fields then
		for _,field in ipairs(args.fields) do
			lines:insert('buf['..dst..'].'..field..' = '..fixedCode..';')
		end
	else
		lines:insert('buf['..dst..'] = '..fixedCode..';')
	end
	return lines:concat'\n'
end
GridSolver.BoundaryFixed = BoundaryFixed 

local BoundaryFreeFlow = class(Boundary)
BoundaryFreeFlow.name = 'freeflow'
function BoundaryFreeFlow:getCode(args)
	local dst, src
	if args.minmax == 'min' then
		dst = args.index'j'
		src = args.index'numGhost'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..'-numGhost+j')
		src = args.index(gridSizeSide..'-numGhost-1')
	end
	return self:assignDstSrc(dst, src, args)
end

-- linear extrapolation
local BoundaryLinear = class(Boundary)
BoundaryLinear.name = 'linear'
function BoundaryLinear:getCode(args)
	local dst, i1, i2
	if args.minmax == 'min' then
		dst = args.index'numGhost-j-1'
		i1 = args.index'numGhost-j'
		i2 = args.index'numGhost-j+1'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..'-numGhost+j')
		i1 = args.index(gridSizeSide..'-numGhost+j-1')
		i2 = args.index(gridSizeSide..'-numGhost+j-2')
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
		dst = args.index'numGhost-j-1'
		i1 = args.index'numGhost-j'
		i2 = args.index'numGhost-j+1'
		i3 = args.index'numGhost-j+2'
	elseif args.minmax == 'max' then
		dst = args.index(gridSizeSide..'-numGhost+j')
		i1 = args.index(gridSizeSide..'-numGhost+j-1')
		i2 = args.index(gridSizeSide..'-numGhost+j-2')
		i3 = args.index(gridSizeSide..'-numGhost+j-3')
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
	assert(require 'coord.sphere'.is(solver.coord)
		or require 'coord.sphere-log-radial'.is(solver.coord), "you should only use this boundary condition for rmin with spherical coordinates")
	--assert(solver.maxs.y - solver.mins.y == 2*math.pi)
	--assert(solver.boundaryMethods.ymin == 'periodic' and solver.boundaryMethods.ymax == 'periodic')

	local src, dst
	if solver.dim == 1 then
		dst = 'INDEX(j, 0, 0)'
		src = 'INDEX(2*numGhost-1-j, 0, 0)'
	elseif solver.dim == 2 then	-- r, theta
		dst = 'INDEX(j, i, 0)'
		src = 'INDEX(2*numGhost-1-j, solver->gridSize.y-i-1, 0)'
	elseif solver.dim == 3 then
		dst = 'INDEX(j, i.x, i.y)'
		src = [[
	INDEX(
		2*numGhost-1-j,
		solver->gridSize.y-i.x-1, 
		(i.y-numGhost + (solver->gridSize.z-2*numGhost)/2
			+ (solver->gridSize.z - 2*numGhost)) 
			% (solver->gridSize.z - 2*numGhost) + numGhost
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
	assert(require 'coord.sphere'.is(solver.coord)
		or require 'coord.sphere-log-radial'.is(solver.coord), "you should only use this boundary condition for θmin/θmax with spherical coordinates")

	local src, dst
	if args.minmax == 'min' then
		if solver.dim == 1 then
			dst = args.index'j'
			src = args.index'2*numGhost-1-j'
		elseif solver.dim == 2 then
			dst = args.index'j'
			src = args.index'2*numGhost-1-j'
		elseif solver.dim == 3 then
			dst = 'INDEX(i.x, numGhost-1-j, i.y)'
			src = [[
	INDEX(
		i.x,
		numGhost+j,
		(i.y - numGhost + (solver->gridSize.z - 2 * numGhost) / 2 
			+ (solver->gridSize.z - 2 * numGhost))
			% (solver->gridSize.z - 2 * numGhost) + numGhost
	)
]]
		end
	elseif args.minmax == 'max' then
		if solver.dim == 1 then
			dst = args.index'solver->gridSize.y-1-j'
			src = args.index'solver->gridSize.y-2*numGhost+j'
		elseif solver.dim == 2 then
			dst = args.index'solver->gridSize.y-1-j'
			src = args.index'solver->gridSize.y-2*numGhost+j'
		elseif solver.dim == 3 then
			dst = 'INDEX(i.x, solver->gridSize.y-numGhost+j, i.y)'
			src = [[
	INDEX(
		i.x,
		solver->gridSize.y - numGhost - 1 - j,
		(i.y - numGhost + (solver->gridSize.z - 2 * numGhost) / 2
			+ (solver->gridSize.z - 2 * numGhost))
			% (solver->gridSize.z - 2 * numGhost) + numGhost
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
	assert(require 'coord.cylinder'.is(solver.coord), "you should only use this boundary condition for rmin with cylinderical coordinates")
	--assert(solver.maxs.y - solver.mins.y == 2*math.pi)
	--assert(solver.boundaryMethods.ymin == 'periodic' and solver.boundaryMethods.ymax == 'periodic')

	local src, dst
	if solver.dim == 1 then
		dst = 'INDEX(j, 0, 0)'
		src = 'INDEX(2*numGhost-1-j, 0, 0)'
	elseif solver.dim == 2 then	-- r, theta
		dst = 'INDEX(j, i, 0)'
		src = [[
INDEX(
	2*numGhost-1-j, 
	(i-numGhost + (solver->gridSize.y-2*numGhost)/2
		+ (solver->gridSize.y - 2*numGhost))
		% (solver->gridSize.y - 2*numGhost) + numGhost,
	0)
]]
	elseif solver.dim == 3 then
		dst = 'INDEX(j, i.x, i.y)'
		src = [[
INDEX(
	2*numGhost-1-j, 
	(i.x-numGhost + (solver->gridSize.y-2*numGhost)/2
		+ (solver->gridSize.y - 2*numGhost))
		% (solver->gridSize.y - 2*numGhost) + numGhost,
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
	self.boundaryOptions = table()
	self.boundaryOptionNames = table()
	self.boundaryOptionForName = {}
	
	self:addBoundaryOptions{
		BoundaryNone,
		BoundaryPeriodic,
		BoundaryMirror,
		--BoundaryFixed,
		BoundaryFreeFlow,
		BoundaryLinear,
		BoundaryQuadratic,
		BoundarySphereRMin,
		BoundarySphereTheta,
		BoundaryCylinderRMin,
	}

	if self.eqn.createBoundaryOptions then
		self.eqn:createBoundaryOptions()
	end
end
function GridSolver:addBoundaryOptions(args)
	for _,arg in ipairs(args) do
		self:addBoundaryOption(arg)
	end
end
function GridSolver:addBoundaryOption(boundaryClass)
	self.boundaryOptions:insert(assert(boundaryClass))
	self.boundaryOptionNames:insert(assert(boundaryClass.name))
	self.boundaryOptionForName[boundaryClass.name] = boundaryClass
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
			local boundaryObjArgs
			if type(name) == 'table' then
				boundaryObjArgs = name.args
				name = name.name
			end
			local boundaryClass = self.boundaryOptionForName[name]
			if not boundaryClass then
				io.stderr:write("unable to find boundary method "..tostring(name)
					..' ('..type(name)..')'
					.." -- can't assign it to side "..k.."\n")
				io.stderr:flush()
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
			lines:insert(args.methods[xNames[side]..minmax]:getCode({
				
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
			}))
		end 

lines:insert[[
	}
}
]]
	end

	local code = lines:concat'\n'

	local boundaryProgramObj
if self.useCLLinkLibraries then 
	time('compiling boundary program', function()
		boundaryUnlinkedObj = self.Program{name='boundary', code=code}
		boundaryUnlinkedObj:compile{dontLink=true}
	end)
	time('linking boundary program', function()
		boundaryProgramObj = self.Program{
			programs = {
				self.mathUnlinkedObj, 
				boundaryUnlinkedObj,
			},
		}
	end)
else	
	time('building boundary program', function()
		boundaryProgramObj = self.Program{name='boundary', code=code}
		boundaryProgramObj:compile()
	end)
end	
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
	end
	for _,op in ipairs(self.ops) do
		if op.refreshBoundaryProgram then
			op:refreshBoundaryProgram()
		end
	end

	if self.useCTU and self.usePLM then
		local args = self:getBoundaryProgramArgs()
		args.type = self.eqn.consLR_t..'_dim'
	
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
				self.localSize2d[1])
			local maxSizeY = roundup(
				math.max(tonumber(self.gridSize.y), tonumber(self.gridSize.z)),
				self.localSize2d[2])
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
	numStates = numStates or self.eqn.numIntStates
	local exact = assert(self.eqn.initState.exactSolution, "can't test accuracy of a configuration that has no exact solution")
	local ptr = ffi.cast(self.eqn.cons_t..'*', self.UBufObj:toCPU())
	assert(self.dim == 1)
	local n = self.gridSize.x
	local ghost = self.numGhost
	local err = 0
	for i=ghost+1,tonumber(self.gridSize.x)-2*ghost do
		--local x = self.xs[i+ghost]
		local x = fromreal(self.solverPtr.mins.x) + fromreal(self.solverPtr.grid_dx.x) * (i - ghost - .5)
		err = err + compareL1(ptr[i-1].ptr, numStates, exact(self, x, self.t))
	end
	err = err / (numStates * (tonumber(self.gridSize.x) - 2 * ghost))
	return err, ptr
end

function GridSolver:update()
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
	
		if self.calcDeriv then
			self:calcDeriv(derivBufObj, dt)
		end

if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end

		if self.eqn.useSourceTerm then
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

	-- I moved the constrainU() and boundary() functions from here to be within the integrator, 
	-- so that each sub-step that the RK integrator performs can apply these and be physically correct states.

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

function GridSolver:calcDisplayVarToTex(var, componentIndex)
	componentIndex = componentIndex or var.component
	local component = self.displayComponentFlatList[componentIndex]
	local vectorField = self:isVarTypeAVectorField(component.type)

	local app = self.app
	local displayVarGroup = var.displayVarGroup
	if app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		self.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	
		var:setToTexArgs()
		var.calcDisplayVarToTexKernelObj()
		
		self.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		self.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		-- calcDisplayVarToTexPtr is sized without ghost cells 
		-- so is the GL texture
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex

		local channels = vectorField and 3 or 1
		local format = vectorField and gl.GL_RGB or gl.GL_RED

		var:setToBufferArgs()
	
		self:calcDisplayVarToBuffer(var)

		local sizevec = var.getBuffer().sizevec or self.gridSize
		assert(sizevec.x <= tex.width and sizevec.y <= tex.height)
		local volume = tonumber(sizevec:volume())
		
		self.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(app.real) * volume * channels, ptr=ptr}
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
	
		ig.igText('grid size: '
			..tonumber(self.sizeWithoutBorder.x)..', '
			..tonumber(self.sizeWithoutBorder.y)..', '
			..tonumber(self.sizeWithoutBorder.z))

		-- let the user change even unused dimensions
		-- because those influence what constant value the cell_x() gives us 
		for i=1,3 do
			for j,minmax in ipairs(minmaxs) do
				local k = xNames[i]..minmax
				if tooltip.numberTable(k, self[minmax..'s'], xNames[i], ig.ImGuiInputTextFlags_EnterReturnsTrue) then
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
				if tooltip.comboTable(var, self.boundaryMethods, var, self.boundaryOptionNames) then
					self:refreshBoundaryProgram()
				end
				--]]
				ig.igText(var..': '..self.boundaryMethods[var].name)
			end
		end
	end
end

function GridSolver:updateGUIEqnSpecific()
	if ig.igCollapsingHeader'equation:' then
--[[	TODO why is this crashing
		if tooltip.comboTable('init state', self, 'initStateIndex', self.eqn.initStateNames) then
			-- TODO hmm ... the whole point of making a separate initStateProgram was to be able to refresh it without rebuilding all of the solver ...
			self:refreshEqnInitState()
		end	
--]]		
		for _,var in ipairs(self.eqn.guiVars) do
			var:updateGUI(self)
		end
	end
end

do
	local function handle(self, var, title)
		local anyChanged = false
		ig.igPushIDStr(title)

		var.enabled = not not var.enabled
		local enableChanged = tooltip.checkboxTable('enabled', var, 'enabled') 
		anyChanged = anyChanged or enableChanged
		ig.igSameLine()
		
		anyChanged = anyChanged or tooltip.checkboxTable('log', var, 'useLog')
		ig.igSameLine()

		anyChanged = anyChanged or tooltip.checkboxTable('units', var, 'showInUnits') 
		ig.igSameLine()

		anyChanged = anyChanged or tooltip.checkboxTable('fixed range', var, 'heatMapFixedRange')
		ig.igSameLine()

		--tooltip.comboTable('component', var, 'component', self.displayComponentNames)
	
		local name = var.name
		if var.units then name = name..' '..var.units end
		if ig.igCollapsingHeader(name) then
			local unitScale = 1
			if var.units and var.showInUnits then -- convert our ranges from raw to units
				unitScale = self:convertToSIUnitsCode(var.units).func()
				var.heatMapValueMin = var.heatMapValueMin * unitScale
				var.heatMapValueMax = var.heatMapValueMax * unitScale
			end
			anyChanged = anyChanged or tooltip.numberTable('value min', var, 'heatMapValueMin')
			anyChanged = anyChanged or tooltip.numberTable('value max', var, 'heatMapValueMax')
			if var.units and var.showInUnits then -- convert our ranges from raw to units
				var.heatMapValueMin = var.heatMapValueMin / unitScale
				var.heatMapValueMax = var.heatMapValueMax / unitScale
			end
		end
		
		ig.igPopID()
	
		return enableChanged, anyChanged
	end

	-- do one for 'all'
	local op = require 'ext.op'
	local fields = {
		'enabled',
		'useLog',
		'showInUnits',
		'heatMapFixedRange',
		'heatMapValueMin',
		'heatMapValueMax',
	}
	local defaults = {
		true,
		true,
		true,
		true,
		math.huge,
		-math.huge,
	}
	local combines = {
		op.land,
		op.land,
		op.land,
		op.land,
		math.min,
		math.max,
	}
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
					local enableChanged, anyChanged = handle(self, all, 'all')
					--refresh = refresh or enableChanged
					if anyChanged then
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
							local enableChanged = handle(self, var, title)
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
