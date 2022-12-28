--[[
command-line variables:

solver parameters:
	(config file cmdline override options:)
		solver =
		solverArgs =
			eqn =
			flux =
		gridSize =
		integrator =
		fluxLimiter =
		fixedDT =
		cfl =
		dim =
		boundary =
		initCond =
		mins =
		maxs =
		multiSlices =
	
hardware/opencl:
	float = set to true to use 32 bit float instead of 64
	half = set to true to use 16 bit float instead of 64
	cpu = set to use CPU instead of GPU
	platform = name or numeric index of which OpenCL platform to use.  see clinfo utility (or my cl/tests/info.lua script) for identifying platforms.
	device = name or numeric index of which OpenCL device to use.
	useGLSharing = set to false to disable GL sharing.  automatically false if sys=console.

display:
	sys = subsystem to run under.  options are 'imguiapp', 'glapp', 'console'
	vsync = set to enable vsync and slow down the simulation
	createAnimation = set to start off creating an animation / framedump
	screenshotUseHiRes = save at the buffer resolution instead of the window resolution
	animationSaveEvery = save every N frames.

	displayvars = comma-separated predefined display vars to use
	display_useCoordMap = set this gui option
	display_fixedRange = whether to initialize display vars with a fixed range (of [0,1])
	showMouseCoords = whether to show mouse coords.  default is true.
	
	displayDim = override dimension of display (useful for displaying 3D simulations as 1D graphs)
	displaySlice = 'xy' by default, 'xz' for the xz slice, 'yz' for the yz slice.
		set this to {x,y,z,theta} for angle/axis initialization.
	
	frustum = set frustum view initially
	frustumDist = frustum view initial distance.  default 3
	frustumAngle = frustum view initial angle.  default {0,0,0,1}
	ortho = set ortho view initially
	stackGraphs = stack graphs initially
	disableFont = set to disable loading of the font.png file.  automatically true if sys=console.

	useClipPlanes = default false.  whether to use clip-planes in the 3d-slice shader (and other shaders?)
	isobars = default true.  whether to use isobars in 3D-slice display.  default: true.

	vectorFieldStep = spacing between cells of vector field display

	windowsize = override the window size for sys=glapp or sys=imguiapp.  set to either "windowsize={w,h}" or for square windows, "windowsize=size"

	palette=(palette name) = pick the predefined palette to use

simulation execution:
	run = start the simulation running.  Setting 'exitTime' or 'stopTime' or 'sys=console' also starts the simulation running.
	exitTime = start the app running, and exit it after the simulation reaches this time
	saveOnExit = filename to save all solvers (appending _1 _2 for multiple solvers) before quitting
	plotOnExit = enable or set to 'true' to have a plot popup upon exit.  set it to a filename string to save the plot to that file.
				this plots a time history of variable min, max, avg, and stddev of the trackvars.
	plotOnExit_savedata = (optional) where to save the plotOnExit data
	plot1DOnExit = enable to plot the 1D data at the time of exit. set it to a filename string to save the plot to that file.
				this plots a snapshot of the variable min, max, avg, and stddev of the trackvars.
	plot1DOnExit_savedata = (optional) what filename to save the output data as
	stopTime = stop running once this time is reached.
	maxiter = max # of iterations to run the application for
	checknans = Stop if a NaN or infinity is found. This can also be a comma-separated list of flags to modify the check-NaN behavior:
		all = upon finding a NaN, print all values that were NaN in that buffer.  If this is omitted then only the first value is printed.
		gpu = perform the check using gpu.  should go faster.  Right?
		noghost = skip testing NaNs in ghost cells.
	

debugging:
	verbose = output extra stuff
	showfps = print the updates/second to console
	trackvars = comma-separated list of variable names to print to the console every FPS print
	tick = how often to update the console for fps or trackvars
	coordVerbose = output extra info from coord/coord.lua

	printBufs = print buffer contents for debugging
	config = specify alternative config file.  default is config.lua (TODO configs/default.lua)
	checkStructSizes = verify that ffi and OpenCL are using matching struct sizes

	trace = insert a debug hook to print out where we are every so often.

integrator parameters:
	intVerbose = output extra info from int/*.lua
	intBEEpsilon = backwards Euler stop on residual less than this epsilon
	intBERestart = backwards Euler GMRES restart
	intBEMaxIter = backwards Euler Krylov max iter

poisson solver parameters:
	noDivPoissonSolver = 'jacobi' or 'krylov'
	selfGravPoissonSolver
	selfGravPoissonMaxIter
	selfGravVerbose
	selfGravLinearSolver = 'conjgrad', 'conjres', 'bicgstab', 'gmres'

code cache parameters:
	useCache = true by default.  set to 'false' to disable the binary-caching system and instead recompile the CL code when you run.
	bssnUseCache = set to false to ignore the cached equations generated by symmath
	usecachedcode = set to true to use the code in cache-bin instead of regenerating it
		TODO I BROKE THIS
		when I switched kernel names to being based on object hash(pointer)s
		because now the kernel name changes every time the program is run.

environment variables:
	HYDROCL_ENV = set to a Lua enumeration of key=value,key=value,... (table grammar without wrapping {})

--]]
local table = require 'ext.table'
local ffi = require 'ffi'

--[[ order of application of command-line key/values:
1) HYDROCL_ENV env var
2) previous global 'cmdline' definition
3) the command-line ... key/values
--]]
do
	-- save the previous global if it's there, apply it later
	local oldcmdline = cmdline
	
	-- initialize our new global cmdline
	cmdline = table()

	-- first handle HYDROCL_ENV env var
	local envvarptr = os.getenv'HYDROCL_ENV'
	if envvarptr ~= nil then	-- cdata NULL will cast to boolean as 'true', but (cdata NULL ~= nil) will evaluate to false
		local envstr = ffi.string(envvarptr)	-- ffi.string segfaults on NULL last I checked
		-- ok i am not going to treat this as a cmdline -- why make parsing more miserable
		-- just treat it as a Lua table
		local env = require 'ext.fromlua'('{'..envstr..'}')
		cmdline = table(env, cmdline)
	end

	-- 2) add in anything from the previously-set cmdline global
	cmdline = table(oldcmdline, cmdline)
	
	-- 3) add in anything from the real command-line
	cmdline = table(require 'ext.cmdline'(table.unpack(arg)), cmdline)
end

if cmdline.srand then
	math.randomseed(os.time())
end

if cmdline.verbose then
	print('cmdline: '..require'ext.tolua'(cmdline))
end

local bit = require 'bit'
local cl = require 'ffi.OpenCL'
local class = require 'ext.class'
local math = require 'ext.math'
local file = require 'ext.file'
local range = require 'ext.range'
local string = require 'ext.string'
local template = require 'template'
local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'
local half = require 'cl.obj.half'
local vec4d = require 'vec-ffi.vec4d'
local vec3d = require 'vec-ffi.vec3d'
local quatf = require 'vec-ffi.quatf'
local CartesianCoord = require 'hydro.coord.cartesian'

--[[
TODO 'targetPlatform'?
options: console, glapp, imguiapp
--]]
local targetSystem = cmdline.sys or 'imguiapp'

-- TODO make this work with the fallback mechanism as well
-- -- that will probalby mean lots of pcalls around requires
local gl
local sdl
local GLProgram
local GLGradientTex
local GLTex2D
local Font
local Mouse
if targetSystem ~= 'console' then
	gl = require 'gl'
	sdl = require 'ffi.sdl'
	GLProgram = require 'gl.program'
	GLGradientTex = require 'gl.gradienttex'
	Mouse = require 'glapp.mouse'
	GLTex2D = require 'gl.tex2d'
	Font = require 'gui.font'
end

-- TODO here if we have glapp and imguiapp present
-- but we're in console mode
-- then we need to figure this out somehow
-- but we won't figure it out until SDL_Init is run, inside GLApp:run()
-- ... is there any way to detect if SDL_Init can init without performing the init?
-- or is there any way to move the :run() into this protected area, so we can bail out depending on the subsystem (without actually loading all the hydro-cl init stuff)
-- until then ... you have to explicitly state sys=console
-- TODO move App:init code into a separate function, then in the loader here try calling :init (and :initGL) and only if all succeeds *THEN* run our :appInit()
-- TODO unless we are still returning a class, in which case this should all be moved into the class's :init() code ... and then we change is-a with has-a, until init() is done, then we swap it back with is-a
local ig
local baseSystems = {
	{imguiapp = function()
		ig = require 'imgui'
		return require 'imguiapp'
	end},
	{glapp = function()
		
		package.loaded['imgui'] = {disabled=true}
		ig = require 'imgui'
		
		return require 'glapp'
	end},
	{console = function()
		
		package.loaded.ig = {disabled=true}
		package.loaded.gl = {disabled=true}
		package.loaded['gl.report'] = {disabled=true}
		
		local cl = class()
		function cl:requestExit() self.done = true end
		function cl:run()
			if self.initGL then self:initGL(gl, 'none') end
			repeat
				if self.update then self:update() end
			until self.done
			if self.exit then self:exit() end
		end
		return cl
	end},
}

local HydroCLApp
for i,sys in ipairs(baseSystems) do
	local name, loader = next(sys)
	if targetSystem == name then
		--print('trying to load system '..name..'...')
		xpcall(function()
			HydroCLApp = class(loader())
		end, function(err)
			io.stderr:write('...load failed with error:')
			io.stderr:write(err..'\n'..debug.traceback())
		end)
		if HydroCLApp then break end
		io.stderr:write(targetSystem..'failed\n')
		if i < #baseSystems then
			targetSystem = next(baseSystems[i+1])
			io.stderr:write('falling back to '..targetSystem..'\n')
		end
	end
end
if not HydroCLApp then
	error 'Somehow you exhausted all possible targets.  At least the console system should have loaded.  Something must be wrong.'
end

if cmdline.windowsize then
	if type(cmdline.windowsize) == 'number' then
		HydroCLApp.width = cmdline.windowsize
		HydroCLApp.height = cmdline.windowsize
	else
		HydroCLApp.width = cmdline.windowsize[1]
		HydroCLApp.height = cmdline.windowsize[2]
	end
end

-- TODO organize this all better
-- this is up above if the autosearch passes imguiapp
-- and it is here if the cmdline explicitly asks for glapp or console
if not ig then
	package.loaded['imgui'] = {disabled=true}
	ig = require 'imgui'
end

HydroCLApp.title = 'Hydrodynamics in OpenCL'

-- list from https://en.wikipedia.org/wiki/Flux_limiter
HydroCLApp.limiters = table{
	{name='donor cell', code='return 0.;'},
	{name='Lax-Wendroff', code='return 1.;'},
	{name='Beam-Warming', code='return r;'},
	{name='Fromm', code='return .5 * (1. + r);'},
	{name='CHARM', code='return max(0., r) * (3. * r + 1.) / ((r + 1.) * (r + 1.));'},
	{name='HCUS', code='return max(0., 1.5 * (r + fabs(r)) / (r + 2.));'},
	{name='HQUICK', code='return max(0., 2. * (r + fabs(r)) / (r + 3.));'},
	{name='Koren', code='return max(0., min(2. * r, min((1. + 2. * r) / 3., 2.)));'},
	{name='minmod', code='return max((real)0., min(r, (real)1.));'},
	{name='Oshker', code='return max(0., min(r, 1.5));	//replace 1.5 with 1 <= beta <= 2'},
	{name='ospre', code='return .5 * (r * r + r) / (r * r + r + 1.);'},
	{name='smart', code='return max(0., min(2. * r, min(.25 + .75 * r, 4.)));'},
	{name='Sweby', code='return max(0., max(min(1.5 * r, 1.), min(r, 1.5)));	//replace 1.5 with 1 <= beta <= 2'},
	{name='UMIST', code='return max(0., min(min(2. * r, .75 + .25 * r), min(.25 + .75 * r, 2.)));'},
	{name='van Albada 1', code='return (r * r + r) / (r * r + 1.);'},
	{name='van Albada 2', code='return 2. * r / (r * r + 1.);'},
	{name='van Leer', code='return max(0., r) * 2. / (1. + r);'},	--return (r + fabs(r)) / (1. + fabs(r));
	{name='monotized central', code='return max(0., min(2., min(.5 * (1. + r), 2. * r)));'},
	{name='superbee', code='return max((real)0., (real)max((real)min((real)1., (real)2. * r), (real)min((real)2., r)));'},
	{name='Barth-Jespersen', code='return .5 * (r + 1.) * min(1., min(4. * r / (r + 1.), 4. / (r + 1.)));'},
}
HydroCLApp.limiterNames = HydroCLApp.limiters:mapi(function(limiter) return limiter.name end)

--[[
setup for the solver
args has: platAndDevicesNames
override this for specific experiments
this can't override float vs double precision yet

make sure to call this after 'real' is defined
 specifically the call to 'refreshGridSize' within it
--]]
function HydroCLApp:setup(args)
	args.self = self
	args.cmdline = cmdline
	args.table = table
	local keys = table.keys(args):sort()
	local cfgfile = cmdline.config or 'config.lua'
	assert(load([[
local ]]..keys:concat', '..[[ = ...
]] .. file(cfgfile):read()))(
	keys:mapi(function(key) return args[key] end):unpack()
)
end

HydroCLApp.useClipPlanes = cmdline.useClipPlanes
if HydroCLApp.useClipPlanes then
	-- TODO put all the display3D_slices stuff in its own file
	HydroCLApp.rotateClip = 0
	HydroCLApp.clipInfos = range(4):mapi(function(i)
		local plane = vec4d(0,0,0,0)
		plane.s[math.min(i,3)-1] = -1
		return {
			enabled = i == 3,
			plane = plane,
		}
	end)
end


-- shader needs it, but shaders are created per-solver (before the app initializes the 'draw' objects)
--HydroCLApp.drawVectorLICNoiseSize = 256
HydroCLApp.drawVectorLICNoiseSize = 1024

function HydroCLApp:display1D(...)
	for _,method in ipairs(self.display1DMethods) do
		local name, func = next(method)
		if self.display1DMethodsEnabled[name] then
			func(self, ...)
		end
	end
end

function HydroCLApp:display2D(...)
	for _,method in ipairs(self.display2DMethods) do
		local name, func = next(method)
		if self.display2DMethodsEnabled[name] then
			func(self, ...)
		end
	end
end

function HydroCLApp:display3D(...)
	for _,method in ipairs(self.display3DMethods) do
		local name, func = next(method)
		if self.display3DMethodsEnabled[name] then
			func(self, ...)
		end
	end
end

function HydroCLApp:displayVector(...)
	for _,method in ipairs(self.displayVectorMethods) do
		local name, func = next(method)
		if self.displayVectorMethodsEnabled[name] then
			func(self, ...)
		end
	end
end

--[[
meshsolver is throwing a wrench in this design
because I'm making one draw object singleton
but the draw object assigns fields to the solver

however ...
mesh_heatmap vs 2d_heatmap is assigned based on the solver type
but only half the routines are single-solver specific

and the per-solver routines are within draw2d, and only those would conditionally create meshsolver, so draw2d would have to bootstrap drawmesh

so why not just make one draw object per solver?
but then, what about overlapping graphs?

well draw1d is still a giant mess wrt this.
--]]

function HydroCLApp:display1D_Graph(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.draw1DGraph = solver.draw1DGraph or require 'hydro.draw.1d_graph'(solver)
		solver.draw1DGraph:display(...)
	end
end

function HydroCLApp:display2D_Heatmap(solvers, ...)
	for _,solver in ipairs(solvers) do
		if require 'hydro.solver.meshsolver':isa(solver) then
			solver.drawMeshHeatmap = solver.drawMeshHeatmap or require 'hydro.draw.mesh_heatmap'(solver)
			solver.drawMeshHeatmap:display(...)
		else	-- gridsolver
			solver.draw2DHeatmap = solver.draw2DHeatmap or require 'hydro.draw.2d_heatmap'(solver)
			solver.draw2DHeatmap:display(...)
		end
	end
end

function HydroCLApp:display2D_Graph(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.draw2DGraph = solver.draw2DGraph or require 'hydro.draw.2d_graph'(solver)
		solver.draw2DGraph:display(...)
	end
end

function HydroCLApp:display3D_Slice(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.draw3DSlice = solver.draw3DSlice or require 'hydro.draw.3d_slice'(solver)
		solver.draw3DSlice:display(...)
	end
end

function HydroCLApp:display3D_Ray(solver, ...)
	for _,solver in ipairs(solvers) do
		solver.draw3DRay = solver.draw3DRay or require 'hydro.draw.3d_ray'(solver)
		solver.draw3DRay:display(...)
	end
end

function HydroCLApp:display3D_Isosurface(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.draw3DIso = solver.draw3DIso or require 'hydro.draw.3d_iso'(solver)
		solver.draw3DIso:display(...)
	end
end

function HydroCLApp:displayVector_Arrows(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.drawVectorArrows = solver.drawVectorArrows or require 'hydro.draw.vector_arrow'(solver)
		solver.drawVectorArrows:display(...)
	end
end

function HydroCLApp:displayVector_LIC(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.drawVectorLIC = solver.drawVectorLIC or require 'hydro.draw.vector_lic'(solver)
		solver.drawVectorLIC:display(...)
	end
end

function HydroCLApp:displayVector_StateLine(solvers, ...)
	for _,solver in ipairs(solvers) do
		solver.drawVectorStateLine = solver.drawVectorStateLine or require 'hydro.draw.vector_state_line'(solver)
		solver.drawVectorStateLine:display(...)
	end
end

--[[
called from initGL
only called if sys ~= console

some design thoughts:
a lot of draw is ui-specific, so it can be app members (like gradient tex, buffers, etc)
but the shader and state buf tex stuff is solver-specific, hence why solvers build their shader and provide it to the draw objects
and originally draw objs were built upon request ... but shaders were built upon solver init ... but shaders are now needing draw objs to be already built ...
so here's an idea:
have hydro make draw objs *and* solver shaders upon request
then create the draw objs before the shaders (so they exist during shader ctor)
--]]
function HydroCLApp:initDraw()
	self.displaySliceAngle = quatf(0,0,0,1)
	if cmdline.displaySlice == 'xz' then
		self.displaySliceAngle:fromAngleAxis(1,0,0,90)
	elseif cmdline.displaySlice == 'yz' then
		self.displaySliceAngle:fromAngleAxis(0,1,0,90)
	elseif type(cmdline.displaySlice) == 'table' then
		assert(#cmdline.displaySlice == 4, "don't know how to handle this cmdline.displaySlice")
		self.displaySliceAngle:fromAngleAxis(table.unpack(cmdline.displaySlice))
	end
	
	self.display1DMethods = table{
		{Graph = HydroCLApp.display1D_Graph},
	}

	self.display2DMethods = table{
		{Heatmap = HydroCLApp.display2D_Heatmap},
		{Graph = HydroCLApp.display2D_Graph},
	}

	self.display3DMethods = table{
		{Slices = HydroCLApp.display3D_Slice},
		{Raytrace = HydroCLApp.display3D_Ray},
		{Isosurfaces = HydroCLApp.display3D_Isosurface},
	}

	self.displayVectorMethods = table{
		{Arrows = HydroCLApp.displayVector_Arrows},
		{LIC = HydroCLApp.displayVector_LIC},
		{StateLine = HydroCLApp.displayVector_StateLine},
	}
end


--[[ Cheap output of the state each frame so I can compare it to other solvers.
-- 	This is what the dump-to-file is supposed to also do.
local printStateFile = file'out.txt':open'w'
local function printState(solver)
	local ptr = solver.UBufObj:toCPU()
	
	local cols = {0, 1, 4}
	
	for i=0,solver.numCells-1 do
		-- matching the cl defs:
		local dx = tonumber((solver.maxs.x - solver.mins.x) / (solver.gridSize.x - (2*solver.numGhost)))
		local x = (i + .5 - solver.numGhost) * dx + solver.mins.x
		printStateFile:write(solver.t,'\t',x)
		for _,j in ipairs(cols) do	-- only use rho, mx, ETotal
		--for j=0,solver.eqn.numStates-1 do
			printStateFile:write('\t',ptr[j + solver.eqn.numStates * i])
		end
		printStateFile:write'\n'
	end
	printStateFile:write'\n'
	printStateFile:flush()
end
--]]

-- always called by glapp, even if sys=console
-- ... in that case a dummy subsys is provided above that calls initGL
function HydroCLApp:initGL(...)
	if HydroCLApp.super
	and HydroCLApp.super.initGL
	then
		HydroCLApp.super.initGL(self, ...)
	end

	self.cmdline = cmdline
	self.targetSystem = targetSystem

	if cmdline.vsync then
		-- Latest intel-opencl-icd drivers are freezing my system, only for this program, and I don't know why.  Turning on vsync seems to fix it.
		sdl.SDL_GL_SetSwapInterval(1)
	end

	-- This used to be on by default,
	-- but for now the 'calcDisplayVarToTex' code has grown out of hand and now doubles the compile times
	-- and I can't perceive a performance difference with or without it,
	-- so I will keep this disabled for now.
	--local useGLSharing = true
	local useGLSharing = false

	if cmdline.useGLSharing ~= nil then useGLSharing = cmdline.useGLSharing end
	if self.targetSystem == 'console' then useGLSharing = false end

	-- TODO if no identifier was specified by the cmdline then favor cl_khr_gl_sharing, cl_khr_fp64, cl_khr_3d_image_writes
	local function getterForIdent(ident, identType)
		return function(objs)
			if ident == nil then return objs end	-- use all
			-- use a sepcific device
			-- TODO how to specify using multiple devices?
			for i,obj in ipairs(objs) do
				if type(ident) == 'number' then
					if ident == i then return {obj} end
				elseif type(ident) == 'string' then
					if ident == obj:getName() then return {obj} end
				end
			end
			error("couldn't find "..identType)
		end
	end

	self.verbose = cmdline.verbose
	self.env = CLEnv{
		precision = cmdline.float and 'float' or (cmdline.half and 'half' or nil),
		cpu = cmdline.cpu,
		useGLSharing = useGLSharing,
		verbose = self.verbose,
		getPlatform = getterForIdent(cmdline.platform, 'platform'),
		getDevices = getterForIdent(cmdline.device, 'device'),
	}
	local platAndDevicesNames = table{
		self.env.platform,
		table.unpack(self.env.devices)
	}:mapi(function(x) return x:getName() end):concat'/'

	self.exitTime = cmdline.exitTime
	self.stopTime = cmdline.stopTime
	if self.exitTime
	or self.stopTime
	or cmdline.run
	or self.targetSystem == 'console'
	then
		self.running = true
	end
	
	self.frameIndex = 0

	self.createAnimation = cmdline.createAnimation
	self.animationSaveEvery = cmdline.animationSaveEvery
	
	-- ok TODO clean this system up a bit because ...
	-- false = use screenshot of the current display window
	-- true = save the contents of the solver's gl tex with values mapped through the current palette in gradientTex
	self.screenshotUseHiRes = not not cmdline.screenshotUseHiRes

	self.useGLSharing = self.env.useGLSharing
	self.ctx = self.env.ctx
	self.real = self.env.real
	
	--half cannot be a kernel param, so this is a proxy type
	self.realparam = self.real == 'half' and 'float' or self.real

	-- half-precision with odd sizes is messing up otherwise
	if self.targetSystem ~= 'console' then
		gl.glPixelStorei(gl.GL_PACK_ALIGNMENT, 1)
		gl.glPixelStorei(gl.GL_UNPACK_ALIGNMENT, 1)
	end

	-- init the module system ...
	do
		local Modules = require 'modules'
		Modules.verbose = cmdline.moduleVerbose
		self.modules = Modules()
	
		self.modules:addFromMarkup(template(file'hydro/code/math.cl':read(), table(require 'hydro.common', {
			app = self,
		})))
		
		require 'hydro.code.safecdef'(self.modules:getTypeHeader'math')
	
		-- this expects solver_t to have gridSize, but it doesn't require its def (because it's a macro)
		self.modules:add{
			name = 'INDEX',
			headercode = '#define INDEX(a,b,c)	((a) + solver->gridSize.x * ((b) + solver->gridSize.y * (c)))',
		}
	
		self.modules:add{
			name = 'INDEXV',
			headercode = '#define INDEXV(i)		indexForInt4ForSize(i, solver->gridSize.x, solver->gridSize.y, solver->gridSize.z)',
		}

		-- add defs.h?

		self.modules:add{
			name = 'numberof',
			headercode = '#define numberof(x)	(sizeof(x)/sizeof(x[0]))',
		}
	
		self.modules:add{
			name = 'endof',
			depends = {'numberof'},
			headercode = '#define endof(x)	((x) + numberof(x))',
		}
	end


	self.solvers = table()

	self:setup{platAndDevicesNames=platAndDevicesNames}
	if #self.solvers == 0 then
		print("You didn't add any solvers in config.lua or in the HydroCLApp:setup() function.  Did you forget something?")
		if self.targetSystem == 'console' then
			self:requestExit()
			return
		end
	end
	
	
	-- this will be per-solver
	-- but is also tightly linked to the structured grid solvers
	-- used for 1D
	if self.targetSystem ~= 'console' then
		-- This only looks good when overlaying vector fields on top of other graphs.
		-- When it comes to separate variables, they usually look better apart.
		self.displayAllTogether = cmdline.stackGraphs or false --self.solvers[1] and self.solvers[1].dim > 1 or false

		self.displayBilinearTextures = cmdline.displayBilinearTextures
		if self.displayBilinearTextures == nil then
			self.displayBilinearTextures = true
		end

		self:resetGradientTex()
		self:initDraw()

		-- this is only used by chopped and meshsolver for post-init display finalizing
		for _,solver in ipairs(self.solvers) do
			if solver.initDraw then
				solver:initDraw()
			end
		end


		self.isobarShader = GLProgram{
			vertexCode = [[
varying vec4 color;
void main() {
	color = gl_Color;
	gl_Position = ftransform();
}
]],
		fragmentCode = [[
varying vec4 color;
void main() {
	float dx_dz = dFdx(gl_FragCoord.z);
	float dy_dz = dFdy(gl_FragCoord.z);

	vec3 n = normalize(vec3(dx_dz, dy_dz, 10.));

	gl_FragColor = vec4(color.rgb * n.z, color.a);
}
]],
		}

		if not cmdline.disableFont then
			local fonttex = GLTex2D{
				filename = 'font.png',
				minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
				magFilter = gl.GL_LINEAR,
			}
			if not pcall(function()
				gl.glGenerateMipmap(gl.GL_TEXTURE_2D)
			end) then
				fonttex:setParameter(gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST)
				fonttex:setParameter(gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
			end
			self.font = Font{tex = fonttex}
		end

		-- todo reorganize me
		self.display1DMethodsEnabled = self.display1DMethods:mapi(function(method, index)
			local name, func = next(method)
			return index == 1, name
		end)
		self.display2DMethodsEnabled = self.display2DMethods:mapi(function(method, index)
			local name, func = next(method)
			return index == 1, name
		end)
		self.display3DMethodsEnabled = self.display3DMethods:mapi(function(method, index)
			local name, func = next(method)
			return index == 1, name
		end)
		self.displayVectorMethodsEnabled = self.displayVectorMethods:mapi(function(method, index)
			local name, func = next(method)
			local enabled = index == 2	-- LIC
			-- cmdline arrows overrides the default from LIC to arrows
			if cmdline.arrows then enabled = index == 1 end
			if cmdline.display_stateline then enabled = index == 3 end
			return enabled, name
		end)

		self.orthoView = require 'hydro.view.ortho'()
		self.frustumView = require 'hydro.view.frustum'()
		self.view = (#self.solvers > 0 and self.solvers[1].dim == 3) and self.frustumView or self.orthoView
		if cmdline.frustum then
			self.view = self.frustumView
		elseif cmdline.ortho then
			self.view = self.orthoView
		end

		if #self.solvers > 0 then
			local solver = self.solvers[1]
			if solver.dim == 2
			and self.display_useCoordMap
			and CartesianCoord:isa(solver.coord)
			then
				local orthoSize = 1
				for j=1,solver.dim do
					self.orthoView.pos.s[j-1] = .5 * (solver.mins.s[j-1] + solver.maxs.s[j-1])
					orthoSize = math.max(orthoSize, solver.maxs.s[j-1] - solver.mins.s[j-1])
				end
				self.orthoView.zoom.x = 1/orthoSize
				self.orthoView.zoom.y = 1/orthoSize
			end
		end
	end

if printState then
	for _,solver in ipairs(self.solvers) do
		printState(solver)
	end
end

	self.displayDim = cmdline.displayDim or self.solvers:mapi(function(solver)
		return solver.dim
	end):sup() or 1
end


-- ok technically this is 32-bit tbgr, msb to lsb.  t for transparency instead of a for alpha.
local function rgb(x)
	local b = bit.band(0xff, x)
	local g = bit.band(0xff, bit.rshift(x, 8))
	local r = bit.band(0xff, bit.rshift(x, 16))
	local o = bit.band(0xff, bit.rshift(x, 24))
	return {r/0xff, g/0xff, b/0xff, 1-o/0xff}
end

HydroCLApp.predefinedPalettes = require 'hydro.draw.palette'

for _,pal in ipairs(HydroCLApp.predefinedPalettes) do
	for i=1,#pal.colors do
		local color = pal.colors[i]
		if type(color) == 'string' then
			color = color
				:gsub('^#', '0x')	-- replace a # prefix on the string with a 0x prefix ...
			local value = tonumber(color)
			assert(value, "failed to decode color "..color)
			pal.colors[i] = rgb(value)
		end
	end
end

HydroCLApp.predefinedPaletteNames = HydroCLApp.predefinedPalettes:mapi(function(palette) return palette.name end)

HydroCLApp.predefinedPaletteIndexForName = HydroCLApp.predefinedPaletteNames:mapi(function(name,index) return index, name end):setmetatable(nil)

HydroCLApp.predefinedPaletteIndex = cmdline.palette and HydroCLApp.predefinedPaletteIndexForName[cmdline.palette] or 1

function HydroCLApp:setGradientTexColors(colors)
	self.gradientTex = GLGradientTex(1024, colors, false)	-- false = don't wrap the colors...
	self.gradientTex:setWrap{s = gl.GL_REPEAT}	-- ...but do use GL_REPEAT
	-- hmm, only on my AMD, intermittantly the next time the tex is bound it will raise an INVALID_OPERATION upon the next bind
	-- maybe this is all because I'm using TEXTURE_1D for the gradientTex?
	-- maybe AMD doesn't like 1D textures so much?
end

function HydroCLApp:resetGradientTex()
	local colors = assert(self.predefinedPalettes[self.predefinedPaletteIndex].colors)
	self:setGradientTexColors(colors)
end

local primaryColors = {
	{1,0,0},
	{1,1,0},
	{0,1,0},
	{0,1,1},
	{0,0,1},
	{1,0,1},
}

function HydroCLApp:randomizeGradientTex()
	local colors = range(math.random(2,20)):mapi(function()
		--[[ bleh
		return {math.random(), math.random(), math.random(), 0.8}
		--]]
		-- [[ based on hsv
		local x = math.random() * #primaryColors	-- hue
		local i = math.floor(x)
		local f = x - i
		local a = primaryColors[i+1]
		local b = primaryColors[i%#primaryColors+1]
		local c = vec3d(table.unpack(a)) * (1 - f) + vec3d(table.unpack(b)) * f
		local g = c:dot(vec3d(.3, .6, .1))
		local s = math.random()	-- saturation
		s=1-s s=s*s*s s=1-s
		c = c * s + vec3d(g,g,g) * (1 - s)
		local l = math.random()
		if l > .5 then
			l = (l - .5) * 2
			c = vec3d(1,1,1) * l + c * (1 - l)
		else
			l = l * 2
			c = c * l + vec3d(0,0,0) * (1 - l)
		end
		return {c.x, c.y, c.z, 0.8}
		--]]
	end)
	self:setGradientTexColors(colors)
end

--[[
rendering:
graph variables - many - 1D and 2D ... and hypercoordinates in 3D?
heatmap variable - 2D only
volumetric variable - 3D only
point cloud? i.e. cheap volumetric?
slices? another cheap volumetric?
--]]

--[[

display options
	ortho vs frustum for 2D, (ortho only for 1D, frustum-only for 3D?)
	graph variable(s?) for 1D and 2D
	1D graph: whether the y axis is auto-scaling or fixed
	heatmap variable(s?) for 2D
	2D heatmap: whether the color range is auto-scaling or fixed
	2D graph: whether the z-axis is auto-scaling or fixed
	volumetric variable(s?) for 3D
	3D volumetric: whether the color range is auto-scaling or fixed
	3D volumetric alpha / gamma values

gui:
	checkbox next to each variable for quickly adding and removing views (based on the first / primary var in the view)
	list of views, with their details next to them, for adding and removing more variables

--]]

local View = class()

local Display = class()

function Display:init(args)
	self.ortho = true

	self.mins = vec3d(-1,-1,-1)
	self.maxs = vec3d(1,1,1)
	
	--[[
	var name, does it contribute to the graph space autoscale range?
	--]]
	self.graphVars = table()
	
	--[[
	var name, heatmap scale min/max, heatmap autoscale? texture, alpha
	--]]
	self.heatMapVars = table()
end

-- whether to output a text file that holds all variable ranges
-- TODO how about some way to track any variable in the gui?
local dumpFile = {
	enabled = false,
}

function dumpFile:update(app, t)
	if not self.enabled then return end

	local f = self.file
	if not f then
		f = file'var-ranges.txt':open'w'
		self.file = f
		
		-- don't change any vars while outputting or else your rows won't match your header
		f:write'#t'
		
		for _,solver in ipairs(app.solvers) do
			--[[ display variables:
			for _,var in ipairs(solver.displayVars) do
				if var.enabled then
					f:write('\t',var.name..'_min')
					f:write('\t',var.name..'_max')
				end
			end
			--]]
			-- [[ gmres error
			f:write'\tgmres_err'
			f:write'\tgmres_iter'
			--]]
		end
		f:write'\n'
		f:flush()
	end
	
	f:write(t)
	for _,solver in ipairs(app.solvers) do
		--[[ display variables:
		for _,var in ipairs(solver.displayVars) do
			if var.enabled then
				local ymin, ymax = solver:calcDisplayVarRange(var)
				-- TODO adjust the range by units
				f:write(('\t%.16e'):format(ymin))
				f:write(('\t%.16e'):format(ymax))
			end
		end
		--]]
		-- [[ gmres residual
		f:write('\t'..solver.integrator.lastResidual)
		f:write('\t'..solver.integrator.lastIter)
		--]]
	end
	f:write'\n'
	f:flush()
end


-- dropdown options
HydroCLApp.screenshotExts = {'png', 'bmp', 'jpeg', 'tiff', 'fits', 'tga', 'ppm'}
-- dropdown index
HydroCLApp.screenshotExtIndex = 1

-- saving a screenshot of whatever the window buffer is, including colorbar and grid, excluding imgui

function HydroCLApp:getScreenShotFilename()
	local ext = self.screenshotExts[self.screenshotExtIndex]

	-- TODO only once upon init?
	if not file'screenshots':exists() then
		-- don't assert -- if it already exists the cmd will fail
		file'screenshots':mkdir()
	end

	-- make a new subdir for each application instance ... ?
	if not self.screenshotDir then
		self.screenshotDir = os.date'%Y.%m.%d-%H.%M.%S'
		local dir = 'screenshots/'..self.screenshotDir
		assert(not file(dir):exists(), "found a duplicate screenshot timestamp subdir")
		assert(file(dir):mkdir())
		self.screenshotIndex = 0
	end

	local fn = ('screenshots/'..self.screenshotDir..'/%05d.'..ext):format(self.screenshotIndex)
	self.screenshotIndex = self.screenshotIndex + 1

	return fn
end

function HydroCLApp:screenshot()
	self:screenshotToFile(self:getScreenShotFilename())	-- in GLApp
end

-- save the visual buffers with their palettes
-- what about 1D? what about vector fields?
-- I guess this is why screenshotting the buffer is a good default behavior

function HydroCLApp:saveHeatMapBufferImages()
	local pushRunning = self.running
	self.running = false
	local pushDrawGradientLegend = self.drawGradientLegend
	self.drawGradientLegend = function() end
	local pushUpdateGUI = self.super.update
	self.super.update = function() end
	local pushFont = self.font
	self.font = nil

	local FBO = require 'gl.fbo'
	local fbo = FBO()
	for _,solver in ipairs(self.solvers) do
		local tex = solver.tex
		local cl = tex.class
		assert(not tex.depth) 	-- i don't have ssimg big enough for glgetteximage of texture_3d yet...
		
		-- TODO store this if you want to use this for streaming
		local fbotex = self.saveHeatMapBufferImages_fbotex
		if not fbotex
		or self.saveHeatMapBufferImages_width ~= tex.width
		or self.saveHeatMapBufferImages_height ~= tex.height
		or self.saveHeatMapBufferImages_depth ~= tex.depth
		then
			fbotex = cl{
				width = tex.width,
				height = tex.height,
				depth = tex.depth,
				internalFormat = gl.GL_RGBA32F,
				format = gl.GL_RGB,
				type = gl.GL_UNSIGNED_BYTE,
				minFilter = gl.GL_NEAREST,
				magFilter = gl.GL_LINEAR,
				wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT, r=gl.GL_REPEAT},
			}
			self.saveHeatMapBufferImages_fbotex = fbotex
			self.saveHeatMapBufferImages_width = tex.width
			self.saveHeatMapBufferImages_height = tex.height
			self.saveHeatMapBufferImages_depth = tex.depth
		end

		local pushSize = self.size
		self.size = function(self) return tex.width, tex.height end

		fbo:setColorAttachment(0, fbotex)
		fbo:bind()
				
		gl.glViewport(0, 0, tex.width, tex.height)
		
		-- TODO instead of pushing and popping
		-- just separate out this function
--		self:update()

		local Image = require 'image'
		local w, h = tex.width, tex.height
		-- using 3 channels had some alignment problems ... there's a bug to fix somewhere, maybe in the png write function?
		local ssimg = Image(w, h, 4, 'unsigned char')
		local ssflipped = Image(w, h, 4, 'unsigned char')
		--tex:toCPU(ssimg.buffer)
		gl.glReadPixels(0, 0, w, h, gl.GL_RGBA, gl.GL_UNSIGNED_BYTE, ssimg.buffer)
		-- TODO just use the shaders?  just draw the screen with this into the buffer and then grab later?
		-- reverse rows ...
		ssimg:flip(ssflipped)
		-- full alpha
		for i=0,w*h-1 do
			ssflipped.buffer[3+4*i] = 255
		end

		-- TODO prefix for solver and for buffer name?
		local fn = self:getScreenShotFilename()
		ssflipped:save(fn)

		fbo:unbind()
	
		self.size = pushSize
	end
	
	self.running = pushRunning
	self.drawGradientLegend = pushDrawGradientLegend
	self.super.update = pushUpdateGUI
	self.font = pushFont
end

local mouse = Mouse and Mouse() or nil

local function canHandleMouse()
	if not mouse then return false end
	if rawget(ig, 'disabled') then return true end
	return not ig.igGetIO()[0].WantCaptureMouse
end

local function canHandleKeyboard()
	if rawget(ig, 'disabled') then return true end
	return not ig.igGetIO()[0].WantCaptureKeyboard
end

HydroCLApp.running = false
--HydroCLApp.running = true

local pushVarNamesEnabled

-- any smaller than this and the font starts to screw up (maybe because it is using floating point?)
local minDeltaY = 1e-5

local time, getTime = table.unpack(require 'hydro.util.time')
local startTime

function HydroCLApp:requestExit(...)
	if startTime then
		local endTime = getTime()
		print("duration: "..(endTime - startTime))
		startTime = nil
	end
	HydroCLApp.super.requestExit(self, ...)
end

local function logabsxform(x)
	x = math.log(math.abs(x), 10)
	if not math.isfinite(x) then x = -math.huge end
	x = math.max(-30, x)
	return x
end
local function getminmax(solverymin, solverymax)
	solverymin = logabsxform(solverymin)
	solverymax = logabsxform(solverymax)
	return math.min(solverymin, solverymax),
			math.max(solverymin, solverymax)
end

function HydroCLApp:update(...)
	if not startTime then
		startTime = getTime()
	end

	if canHandleMouse() then
		mouse:update()
	end
	if self.running then
		if self.running == 'step' then
			print('performing single step...')
			self.running = false
		end

		-- update the one furthest behind
		local oldestSolver = self.solvers:inf(function(a,b)
			return a.t < b.t
		end)
		if oldestSolver then
			-- check before :update(), in case we want to exit at t=0
			if self.exitTime and oldestSolver.t >= self.exitTime then
				-- save on exit?
				if cmdline.saveOnExit then
					for i,solver in ipairs(self.solvers) do
						solver:save(tostring(cmdline.saveOnExit)
							.. (#self.solvers == 1 and '' or (
								'_'..i
							)))
					end
				end
				if cmdline.plotOnExit then
					-- fix gnuplot formatting
					local function fixtitle(s)
						return (s:gsub('_', ' '))
					end
					local varnames = string.split(cmdline.trackvars, ','):mapi(string.trim)
					varnames:removeObject'dt'
					-- plot any track vars we got
					-- 3*#varnames+1 cols per solver
					-- the first is the solver's t
					-- then every 3 is min, avg, max
					local cols = table()
					local usings = table()
					if type(cmdline.plotOnExit) == 'string' then
						usings.output = cmdline.plotOnExit
						local name, ext = file(cmdline.plotOnExit):getext()
						usings.terminal = ext.." size 1600,900 background '#ffffff'"
					else
						usings.persist = true
					end
					for _,solver in ipairs(self.solvers) do
						cols:insert(solver.plotsOnExit.t)
						local ti = #cols
						for _,varname in ipairs(varnames) do
							local suffixes = {'avg', 'min', 'max', 'stddev'}
							local colbase = #cols
							-- these should match 'suffixes' for index when inserting into cols
							local colavg = colbase + 1
							local colmin = colbase + 2
							local colmax = colbase + 3
							local colstddev = colbase + 4
							for _,suffix in ipairs(suffixes) do
								cols:insert(solver.plotsOnExit[varname..' '..suffix])
							end
							usings:insert{using=ti..':'..colavg, title=fixtitle(solver.name..' '..varname..' avg')}
							-- [[ plot each separately
							usings:insert{using=ti..':'..colmax, title=fixtitle(solver.name..' '..varname..' max')}
							usings:insert{using=ti..':'..colmin, title=fixtitle(solver.name..' '..varname..' min')}
							--]]
							-- hmm, I am maybe thinking calculating the stddev isn't useful ...
							--[[ plot avg +- stddev ... can extend past min/max ...
							usings:insert{using=ti..':($'..colavg..'+$'..colstddev..')', title=fixtitle(solver.name..' '..varname..' avg+stddev')}
							usings:insert{using=ti..':($'..colavg..'-$'..colstddev..')', title=fixtitle(solver.name..' '..varname..' avg-stddev')}
							--]]
							-- hmm, for skewed data, boxplot / candlesticks doesn't look good
							--[[ plot all as a boxplot
							usings:insert{using=ti..':($'..colavg..'-$'..colstddev..'):'..colmin..':'..colmax..':($'..colavg..'+$'..colstddev..')', title=fixtitle(solver.name..' '..varname), with='candlesticks whiskerbars'}
							--]]
						end
					end
					require 'gnuplot'(table({
						style = 'data lines',
						data = cols,
						savedata = cmdline.plotOnExit_savedata,
					}, usings):setmetatable(nil))
					-- TODO wouldn't hurt to add a header column to the savedata
				end
				if cmdline.plot1DOnExit then
					local fromreal = half.fromreal
					-- fix gnuplot formatting
					local function fixtitle(s)
						return (s:gsub('_', ' '))
					end
					local varnames = string.split(cmdline.trackvars, ','):mapi(string.trim)
					varnames:removeObject'dt'
					-- plot trackvars, 1 col per var per solver
					-- the first is the solver's x
					-- then every 1 is the value
					local cols = table()
					for _,solver in ipairs(self.solvers) do
						cols:insert(table())	-- xs
						for _,varname in ipairs(varnames) do
							cols:insert(table())
						end
					end
					local usings = table()
					if type(cmdline.plot1DOnExit) == 'string' then
						usings.output = cmdline.plot1DOnExit
					end
					local ci = 1
					for _,solver in ipairs(self.solvers) do
						local xi = ci
						ci = ci + 1
						-- TODO only iterate along the 1st dim, keeping 2nd and 3rd fixed
						for index=0,solver.numCells-1 do
							-- only writing the 1st coord for now
							cols[xi + 0]:insert(tonumber(solver.cellCpuBuf[index].pos.s[0]))
						end
						for _,varname in ipairs(varnames) do
							local colbase = ci
							ci = ci + 1

							-- TODO here insert the values as 'cols'
							local var = assert(solver.displayVarForName[varname], "couldn't find "..varname)
							local component = solver.displayComponentFlatList[var.component]
							local vectorField = solver:isVarTypeAVectorField(component.type)
							local unitScale = solver:convertToSIUnitsCode(var.units).func()
							solver:calcDisplayVarToTex(var)
							local ptr = ffi.cast(self.real == 'half' and 'real*' or 'float*', solver.calcDisplayVarToTexPtr)
							local channels = vectorField and 3 or 1
							for index=0,solver.numCells-1 do
								-- only writing the 1st var for now
								cols[colbase + 0]:insert(fromreal(ptr[0 + channels * index]) or '-')
							end
					
							usings:insert{using=xi..':'..colbase, title=fixtitle(solver.name..' '..varname..' at t='..solver.t)}
						end
					end
					require 'gnuplot'(table({
						persist = true,
						style = 'data lines',
						data = cols,
						savedata = cmdline.plot1DOnExit_savedata,
					}, usings):setmetatable(nil))
					-- TODO wouldn't hurt to add a header column to the savedata
				end

				self:requestExit()
				return
			end
			if self.stopTime and oldestSolver.t >= self.stopTime then
				self.running = false
			end
if cmdline.printBufs then
	print(('t = %f'):format(oldestSolver.t))
end
			oldestSolver:update()

if cmdline.printBufs then
	print()
	print'UBuf post-update:'
	oldestSolver:printBuf(oldestSolver.UBufObj)
end

if printState then
	for _,solver in ipairs(self.solvers) do
		printState(solver)
	end
end


			-- TODO should the time be oldestSolver.t after oldestSolver just updated?
			-- or - if dumpFile is enabled - should we re-search-out the oldest solver and use its time?
			dumpFile:update(self, oldestSolver.t)
		end
	else
		-- clear all 'lastFrameTime's of solvers so the rough fps calcs don't get messed with
		for _,solver in ipairs(self.solvers) do
			solver.lastFrameTime = nil
		end
	end


	self.iteration = (self.iteration or 0) + 1
	if cmdline.maxiter and self.iteration > cmdline.maxiter then
		self:requestExit()
	end


	if self.targetSystem == 'console' then return end

	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))


	-- [[ TODO FIXME temp hack for composite solvers
	local flattenedSolvers = table(self.solvers)
	do
		local i = 1
		while i <= #flattenedSolvers do
			local solver = flattenedSolvers[i]
			if solver.solvers then
				flattenedSolvers:append(solver.solvers)
				flattenedSolvers:remove(i)
				i = i - 1
			end
			i = i + 1
		end
	end
	local displaySolvers = flattenedSolvers
	--]]

	local w, h = self:size()

	local varNamesEnabled = table()
	local varNamesEnabledByName = {}
	for _,solver in ipairs(displaySolvers) do
		for i,var in ipairs(solver.displayVars) do
			if var.enabled then
				if not varNamesEnabledByName[var.name] then
					varNamesEnabled:insert(var.name)
					varNamesEnabledByName[var.name] = true
				end
			end
		end
	end

	local graphsWide = math.ceil(math.sqrt(#varNamesEnabled))
	local graphsHigh = math.ceil(#varNamesEnabled/graphsWide)
	local graphCol = 0
	local graphRow = 0

	if self.displayAllTogether then
		graphsWide = 1
		graphsHigh = 1
	end

	local ar = (w / graphsWide) / (h / graphsHigh)

	local mouseClickedOnVar
				
	-- varymin/max is the variable range
	local varymin, varymax

	local useLog
	local vectorField
	for _,varName in ipairs(varNamesEnabled) do
		local xmin, xmax, ymin, ymax
		for _,solver in ipairs(displaySolvers) do
			local var = solver.displayVarForName[varName]
			
			if var and var.enabled
			--and solver.visiblePtr and solver.visiblePtr[0]
			then
				useLog = var.useLog
				local component = solver.displayComponentFlatList[var.component]
				vectorField = solver:isVarTypeAVectorField(component.type)
				local solverxmin, solverxmax = tonumber(solver.mins.x), tonumber(solver.maxs.x)
				solverxmin, solverxmax = 1.1 * solverxmin - .1 * solverxmax, 1.1 * solverxmax - .1 * solverxmin
				if self.displayDim > 1 then	-- solver.dim
					local center = .5 * (solverxmin + solverxmax)
					solverxmin = (solverxmin - center) * ar + center
					solverxmax = (solverxmax - center) * ar + center
				end
	
				-- solverymin/max is the view range
				local solverymin, solverymax


				if self.displayDim > 1 then	-- solver.dim
					solverymin, solverymax = tonumber(solver.mins.y), tonumber(solver.maxs.y)
					solverymin, solverymax = 1.1 * solverymin - .1 * solverymax, 1.1 * solverymax - .1 * solverymin
				else
					if vectorField then
						solverymin, solverymax = solver:calcDisplayVarRange(var, component.magn)
					else
						solverymin, solverymax = solver:calcDisplayVarRange(var)
					end

					local thisvarymin, thisvarymax = solverymin, solverymax
					varymin = varymin and math.min(thisvarymin, varymin) or thisvarymin
					varymax = varymax and math.max(thisvarymax, varymax) or thisvarymax
					
					if useLog then
						solverymin, solverymax = getminmax(solverymin, solverymax)
						solverymax = math.max(solverymax, solverymin + minDeltaY)
					end
					
					if solverymin
					and solverymax
					and solverymin == solverymin
					and solverymax == solverymax
					then
						local base = 10	-- round to nearest base-10
						local scale = 10 -- ...with increments of 10
						solverymin, solverymax = 1.1 * solverymin - .1 * solverymax, 1.1 * solverymax - .1 * solverymin
						local newymin = (solverymin < 0 and -1 or 1) * (math.abs(solverymin) == math.huge and 1e+100 or base ^ math.log(math.abs(solverymin), base))
						local newymax = (solverymax < 0 and -1 or 1) * (math.abs(solverymax) == math.huge and 1e+100 or base ^ math.log(math.abs(solverymax), base))
						solverymin, solverymax = newymin, newymax
						do
							local deltaY = solverymax - solverymin
							if deltaY < minDeltaY then
								solverymax = solverymax + .5 * minDeltaY
								solverymin = solverymin - .5 * minDeltaY
							end
						end
					end
				end
				
				xmin = xmin or solverxmin
				xmax = xmax or solverxmax
				ymin = ymin or solverymin
				ymax = ymax or solverymax
				
				if xmin and solverxmin then xmin = math.min(xmin, solverxmin) end
				if xmax and solverxmax then xmax = math.max(xmax, solverxmax) end
				if ymin and solverymin then ymin = math.min(ymin, solverymin) end
				if ymax and solverymax then ymax = math.max(ymax, solverymax) end
			end
		end
		
		if not xmin or not xmax or xmin ~= xmin or xmax ~= xmax then
			xmin = -5
			xmax = 5
		end
		if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
			ymin = -5
			ymax = 5
		end

		local mouseOverThisGraph
		local mouseInGraphX
		local mouseInGraphY
		if self.displayAllTogether then
			mouseOverThisGraph = true
			mouseInGraphX = mouse.pos.x
			mouseInGraphY = mouse.pos.y
		else
			local vpxmin = graphCol / graphsWide * w
			local vpymin = (1 - (graphRow + 1) / graphsHigh) * h
			local vpw = w / graphsWide
			local vph = h / graphsHigh
			gl.glViewport(vpxmin, vpymin, vpw, vph)
			local mx = mouse.pos.x * self.width
			local my = mouse.pos.y * self.height
			if mx >= vpxmin and mx < vpxmin + vpw
			and my >= vpymin and my < vpymin + vph
			then
				mouseOverThisGraph = true
				mouseInGraphX = (mx - vpxmin) / vpw
				mouseInGraphY = (my - vpymin) / vph
				if mouse.leftClick then
					mouseClickedOnVar = varName
				end
			end
		end
		
		if self.showMouseCoords
		and mouseOverThisGraph
		and (self.displayDim == 1 or self.displayDim == 2)	-- displaySolvers[1].dim == 2		-- or better -- instead check for ortho
		then
--			if self.displayDim == 1 then
--				-- use the already-existing xy min max
--				self.mouseCoord[1] = mouseInGraphX * (xmax - xmin) + xmin
--				self.mouseCoord[2] = mouseInGraphY * (ymax - ymin) + ymin
--			else	-- use the ortho xy min max
				local xmin, xmax, ymin, ymax
				if self.view.getOrthoBounds then
					xmin, xmax, ymin, ymax = self.view:getOrthoBounds(ar)
				else
					xmin, xmax, ymin, ymax = graph_xmin, graph_ymin, graph_xmax, graph_ymax
				end
				
				-- frustum doesn't have these ...
				if xmax and ymax and xmin and ymin then
					self.mouseCoord[1] = mouseInGraphX * (xmax - xmin) + xmin
					self.mouseCoord[2] = mouseInGraphY * (ymax - ymin) + ymin
				end
--			end
		end


	
		if not vectorField then
			if self.displayDim == 1 then
				self:display1D(displaySolvers, varName, ar, xmin, xmax, ymin, ymax, useLog, varymin, varymax)
			elseif self.displayDim == 2 then
				self:display2D(displaySolvers, varName, ar, xmin, xmax, ymin, ymax)
			elseif self.displayDim == 3 then
				self:display3D(displaySolvers, varName, ar, xmin, xmax, ymin, ymax)
			end
		else
			self:displayVector(displaySolvers, varName, ar, xmin, xmax, ymin, ymax)
		end

	
		-- in all these above :display...() methods, they exclude meshsolvers
		-- so this is just for mesh ...
		for _,solver in ipairs(self.solvers) do
			--if require 'hydro.solver.meshsolver':isa(solver) then
			if solver.display then
				solver:display(varName, ar)
			end
		end
	
			
		-- TODO make this custom per-display-method
		-- (that would also let us do one less tex bind/unbind)
		if mouseOverThisGraph
		and self.showMouseCoords
		and (self.displayDim == 1 or self.displayDim == 2)
		then
			local toreal, fromreal = half.toreal, half.fromreal
			self.mouseCoordValue = ''
			for i,solver in ipairs(displaySolvers) do
				local var = solver.displayVarForName[varName]
				if var and var.enabled then
					-- translate the mouse coords to texture coords
					-- and read the texel at the mouse position
					if self.display_useCoordMap
					and not CartesianCoord:isa(solver.coord)
					then
						--print'FIXME'
						-- run coords through inverse
					else
						--local tcX = (self.mouseCoord[1] - solver.mins.x) / (solver.maxs.x - solver.mins.x)
						--local tcY = (self.mouseCoord[2] - solver.mins.y) / (solver.maxs.y - solver.mins.y)
						local tcX = (self.mouseCoord[1] + 1) * .5
						local tcY = (self.mouseCoord[2] + 1) * .5
						local tcZ = math.clamp(tonumber((self.displayFixedZ - solver.mins.z) / (solver.maxs.z - solver.mins.z)), 0, 1)
						-- TODO something equivalent for 3D ... somehow ...
						if tcX >= 0 and tcX < 1
						and tcY >= 0 and tcY < 1
						then
							-- TODO this is going to include ghost cells...
							local size
							local texX = tcX
							local texY = tcY
							local texZ = tcZ	-- displayFixedX, remapped from coordinate to texture space
							do -- if require 'hydro.solver.meshsolver':isa(solver) then
								size = var.group.getBuffer().sizevec or solver.gridSize
								
								-- TODO now that I am allowing that slicing stuff, now I have to incorporate those transforms here
								-- esp if displayDim==2 but the solver dim==3
								texX = math.floor(texX * tonumber(size.x))
								texY = math.floor(texY * tonumber(size.y))
								texZ = math.floor(texZ * tonumber(size.z))
							end
							if self.useGLSharing then
								print'FIXME'
								--gl.glGetTexSubImage ... why isn't this in any OpenGL header?
							else
								local component = solver.displayComponentFlatList[var.component]
								local vectorField = solver:isVarTypeAVectorField(component.type)
								local unitScale = solver:convertToSIUnitsCode(var.units).func()
								-- do I have to recalculate it here?
								solver:calcDisplayVarToTex(var)
								
								-- ... if we know we're always copying intermediately to reduceBuf then we can use that instead
								-- in fact, we can even use the CPU intermediate buffer
									-- right now I'm copying halfs from cl to gl as-is, so calcDisplayVarToTexPtr will have half data for half and float data otherwise
								local ptr = ffi.cast(self.real == 'half' and 'real*' or 'float*', solver.calcDisplayVarToTexPtr)
								local channels = vectorField and 3 or 1
								self.mouseCoordValue = self.mouseCoordValue
									.. 'int: '..tostring(texX)
									..','..tostring(texY)
									..'\n'
								if size then
									self.mouseCoordValue = self.mouseCoordValue
										..'raw value: '
									local sep = ''
									for j=0,channels-1 do
										local v = fromreal(ptr[j + channels * (texX + size.x * (texY + size.y * texZ))])
										self.mouseCoordValue = self.mouseCoordValue .. sep .. v
										sep = ', '
									end
									self.mouseCoordValue = self.mouseCoordValue
										..'\n'
										..'unit value: '
									sep = ''
									for j=0,channels-1 do
										local v = fromreal(ptr[j + channels * (texX + size.x * (texY + size.y * texZ))])
										v = v * unitScale
										self.mouseCoordValue = self.mouseCoordValue .. sep .. v
										sep = ', '
									end
									self.mouseCoordValue = self.mouseCoordValue
										..'\n'
								end
							end
						end
					end
				end
			end
		end
	
		if not self.displayAllTogether then
			graphCol = graphCol + 1
			if graphCol == graphsWide then
				graphCol = 0
				graphRow = graphRow + 1
			end
		end
	end

	gl.glViewport(0,0,w,h)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(0, w/h, 0, 1, -1, 1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	if self.font then
		local solverNames = displaySolvers:mapi(function(solver)
			return {
				text = ('(%.3f) %s'):format(solver.t, solver.name),
				color = solver.color,
			}
		end)
		local fontSizeX = .02
		local fontSizeY = .02
		local maxlen = solverNames:mapi(function(solverName)
			return self.font:draw{
				text = solverName.text,
				fontSize = {fontSizeX, -fontSizeY},
				dontRender = true,
				multiLine = false,
			}
		end):inf()
		for i,solverName in ipairs(solverNames) do
			self.font:draw{
				pos = {w/h-maxlen,fontSizeY*(i+1)},
				text = solverName.text,
				color = {solverName.color.x, solverName.color.y, solverName.color.z, 1},
				fontSize = {fontSizeX, -fontSizeY},
				multiLine = false,
			}
		end
	end

	-- screenshot before gui
	if self.createAnimation then
		if not self.animationSaveEvery
		or self.frameIndex % self.animationSaveEvery == 0
		then
			if self.screenshotUseHiRes then
				self:saveHeatMapBufferImages()
			else
				self:screenshot()
			end
		end
		
		if self.createAnimation == 'once' then
			self.createAnimation = nil
		end
	end

	if HydroCLApp.super
	and HydroCLApp.super.update
	then
		HydroCLApp.super.update(self, ...)
	end

	if mouseClickedOnVar then
		if not pushVarNamesEnabled then
			pushVarNamesEnabled = table(varNamesEnabled)
			-- and disable all except the one we clicked
			for _,solver in ipairs(displaySolvers) do
				for i,var in ipairs(solver.displayVars) do
					var.enabled = var.name == mouseClickedOnVar
				end
			end
		else
			-- restore those that were pushed
			for _,solver in ipairs(displaySolvers) do
				for i,var in ipairs(solver.displayVars) do
					var.enabled = pushVarNamesEnabled:find(var.name)
				end
			end
			pushVarNamesEnabled = nil
		end
	end

	-- inc frameIndex last so we always save the first frame
	-- TODO also use this in the :screenshot() function?
	self.frameIndex = self.frameIndex + 1
end

-- helper function for all draw classes:
function HydroCLApp:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
	-- TODO only draw the first
	if var.showInUnits and var.units then
		local unitScale = solver:convertToSIUnitsCode(var.units).func()
		valueMin = valueMin * unitScale
		valueMax = valueMax * unitScale
		varName = varName..' ('..var.units..')'
	end

	self.orthoView:projection(ar)
	self.orthoView:modelview()
	local xmin, xmax, ymin, ymax = self.orthoView:getOrthoBounds(ar)
	
	if self.font then
		local palwidth = (xmax - xmin) * .1
		self.gradientTex:enable()
		self.gradientTex:bind(0)
		gl.glBegin(gl.GL_QUADS)
		gl.glColor3f(1,1,1)
		gl.glTexCoord1f(0) gl.glVertex2f(xmin, ymin)
		gl.glTexCoord1f(0) gl.glVertex2f(xmin + palwidth, ymin)
		gl.glTexCoord1f(1) gl.glVertex2f(xmin + palwidth, ymax)
		gl.glTexCoord1f(1) gl.glVertex2f(xmin, ymax)
		gl.glEnd()
		self.gradientTex:unbind()
		self.gradientTex:disable()

		local fontSizeX = (xmax - xmin) * .025
		local fontSizeY = (ymax - ymin) * .025
		local ystep = 10^(math.log(ymax - ymin, 10) - 1.5)
		for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
			local value = (y - ymin) * (valueMax - valueMin) / (ymax - ymin) + valueMin
			local absvalue = math.abs(value)
			self.font:draw{
				pos={xmin * .99 + xmax * .01, y + fontSizeY * .5},
				text=(
					(absvalue > 1e+5 or absvalue < 1e-5)
					and ('%.5e'):format(value) or ('%.5f'):format(value)),
				color = {1,1,1,1},
				fontSize={fontSizeX, -fontSizeY},
				multiLine=false,
			}
		end
		self.font:draw{
			pos = {xmin, ymax, gridz},
			text = varName,
			color = {1,1,1,1},
			fontSize = {fontSizeX, -fontSizeY},
			multiLine = false,
		}
	end
end

HydroCLApp.display_useCoordMap = cmdline.display_useCoordMap
if HydroCLApp.display_useCoordMap == nil then HydroCLApp.display_useCoordMap = true end

HydroCLApp.mouse_influenceEquations = false
HydroCLApp.displayFixedY = 0
HydroCLApp.displayFixedZ = 0

function HydroCLApp:updateGUI()
	if ig.igCollapsingHeader'simulation' then
		if ig.igButton(self.running and 'Stop' or 'Start') then
			self.running = not self.running
		end
		ig.igSameLine()
		if ig.igButton'Step' then
			self.running = 'step'
		end
		ig.igSameLine()
		if ig.igButton'Reset' then
			print'resetting...'
			for _,solver in ipairs(self.solvers) do
				solver:resetState()
			end
			self.running = false
		end
		
		if ig.igButton'Save' then
			local savePrefix = os.date'%Y.%m.%d-%H.%M.%S'
			-- save as cfits
			for i,solver in ipairs(self.solvers) do
				solver:save(savePrefix..'_'..tostring(i))
			end
		end
	
		if ig.luatableTooltipCombo('Palette', self, 'predefinedPaletteIndex', self.predefinedPaletteNames) then
			self:resetGradientTex()
		end
	
		ig.igSameLine()
		if ig.igButton'Randomize Palette' then
			self:randomizeGradientTex()
		end

		-- dump min/max(/avg?) of displayvars to a .txt file
		ig.luatableTooltipCheckbox('dump to text file', dumpFile, 'enabled')

		if ig.igButton'Screenshot' then
			self:screenshot()
		end
		ig.igSameLine()

		ig.luatableTooltipCheckbox('screenshot hires', self, 'screenshotUseHiRes')
		ig.igSameLine()

		ig.luatableTooltipCombo('screenshot ext', self, 'screenshotExtIndex', self.screenshotExts)

		if ig.igButton(self.createAnimation and 'stop frame dump' or 'start frame dump') then
			self.createAnimation = not self.createAnimation
		end

		ig.luatableTooltipCheckbox('stack graphs', self, 'displayAllTogether')
		ig.igSameLine()
	
		-- TODO per-solver
		ig.luatableTooltipCheckbox('bilinear textures', self, 'displayBilinearTextures')
		ig.igSameLine()

		-- for 2D heatmap only atm
		ig.luatableTooltipCheckbox('display with coord map', self, 'display_useCoordMap')
		ig.igSameLine()
		
		ig.luatableTooltipCheckbox('show coords', self, 'showMouseCoords')

		--ig.igSameLine()
		--ig.luatableTooltipCheckbox('mouse influence equations', self, 'mouse_influenceEquations')
		

		if ig.igRadioButton_Bool('ortho', self.view == self.orthoView) then
			self.view = self.orthoView
		end
		ig.igSameLine()
		if ig.igRadioButton_Bool('frustum', self.view == self.frustumView) then
			self.view = self.frustumView
		end

		-- TODO per-solver
		for j=1,3 do
			if ig.igRadioButton_Bool(j..'D', self.displayDim == j) then
				self.displayDim = j
			end
			if j < 3 then ig.igSameLine() end
		end

		
		--ig.luatableTooltipSliderFloat('fixed y', self, 'displayFixedY', -10, 10)
		--ig.luatableTooltipSliderFloat('fixed z', self, 'displayFixedZ', -10, 10)
		ig.igPushID_Str'fixed y zoom'
		if ig.igButton'+' then
			self.displayFixedY = self.displayFixedY + .1
		end
		ig.igSameLine()
		if ig.igButton'-' then
			self.displayFixedY = self.displayFixedY - .1
		end
		ig.igSameLine()
		ig.luatableTooltipInputFloatAsText('fixed y', self, 'displayFixedY')
		ig.igPopID()
		
		ig.igPushID_Str'fixed z zoom'
		if ig.igButton'+' then
			self.displayFixedZ = self.displayFixedZ + .1
		end
		ig.igSameLine()
		if ig.igButton'-' then
			self.displayFixedZ = self.displayFixedZ - .1
		end
		ig.igSameLine()
		ig.luatableTooltipInputFloatAsText('fixed z', self, 'displayFixedZ')
		ig.igPopID()

		-- TODO per-solver
		do
			local q = self.displaySliceAngle
			-- [[ TODO replace this with trackball behavior
			if ig.luatableTooltipSliderFloat('slice qw', q, 'w', -1, 1) then q:normalize(q) end
			if ig.luatableTooltipSliderFloat('slice qx', q, 'x', -1, 1) then q:normalize(q) end
			if ig.luatableTooltipSliderFloat('slice qy', q, 'y', -1, 1) then q:normalize(q) end
			if ig.luatableTooltipSliderFloat('slice qz', q, 'z', -1, 1) then q:normalize(q) end
			--]]
		end
		-- [[ fixed planes
		ig.igText'slice:'
		ig.igSameLine()
		if ig.igButton'xy' then
			self.displaySliceAngle:set(0,0,0,1)
		end
		ig.igSameLine()
		if ig.igButton'xz' then
			self.displaySliceAngle:fromAngleAxis(1,0,0,90)
		end
		ig.igSameLine()
		if ig.igButton'yz' then
			self.displaySliceAngle:fromAngleAxis(0,1,0,90)
		end
		--]]

		-- TODO flag for separate/combined displays (esp for ortho view)

		-- TODO flag to toggle slice vs volume display
		-- or maybe checkboxes for each kind?
		
		do
			local dim = self.displayDim
			if dim == 1 then
				ig.igPushID_Str'1D'
				for i,method in ipairs(self.display1DMethods) do
					if i > 1 then ig.igSameLine() end
					local name, func = next(method)
					ig.luatableTooltipCheckbox(name, self.display1DMethodsEnabled, name)
				end
				ig.igPopID()
			elseif dim == 2 then
				ig.igPushID_Str'2D'
				for i,method in ipairs(self.display2DMethods) do
					if i > 1 then ig.igSameLine() end
					local name, func = next(method)
					ig.luatableTooltipCheckbox(name, self.display2DMethodsEnabled, name)
				end
				ig.igPopID()
			elseif dim == 3 then
				ig.igPushID_Str'3D'
				for i,method in ipairs(self.display3DMethods) do
					if i > 1 then ig.igSameLine() end
					local name, func = next(method)
					ig.luatableTooltipCheckbox(name, self.display3DMethodsEnabled, name)
				end
				ig.igPopID()
			end
			
			do
				ig.igPushID_Str'Vector'

				for i,method in ipairs(self.displayVectorMethods) do
					if i > 1 then ig.igSameLine() end
					local name, func = next(method)
					ig.luatableTooltipCheckbox(name, self.displayVectorMethodsEnabled, name)
				end
				
				ig.igPopID()
			end
		end

		if self.view == self.frustumView then
			local q = self.frustumView.angle
			-- [[ TODO replace this with trackball behavior
			if ig.luatableTooltipSliderFloat('frustum angle qx', q, 'x', -1, 1) then q:normalize(q) end
			if ig.luatableTooltipSliderFloat('frustum angle qy', q, 'y', -1, 1) then q:normalize(q) end
			if ig.luatableTooltipSliderFloat('frustum angle qz', q, 'z', -1, 1) then q:normalize(q) end
			if ig.luatableTooltipSliderFloat('frustum angle qw', q, 'w', -1, 1) then q:normalize(q) end
			--]]
		end
		if self.useClipPlanes then
			for i,info in ipairs(self.clipInfos) do
				ig.igPushID_Str('solver clip plane '..i)
				ig.luatableTooltipCheckbox('clip plane enabled', info, 'enabled')
				ig.luatableTooltipSliderFloat('clip plane x', info.plane, 'x', -1, 1)
				ig.luatableTooltipSliderFloat('clip plane y', info.plane, 'y', -1, 1)
				ig.luatableTooltipSliderFloat('clip plane z', info.plane, 'z', -1, 1)
				ig.luatableTooltipSliderFloat('clip plane w', info.plane, 'w', -1, 1)
				ig.igPopID()
			end
		end
	end

	for i,solver in ipairs(self.solvers) do
		ig.igPushID_Str('solver '..i)
		if ig.igCollapsingHeader(solver.name) then
			-- TODO new window for each
			solver:updateGUI()
		end
		ig.igPopID()
	end

	if self.showMouseCoords
	and (self.displayDim == 1 or self.displayDim == 2)
	then
		ig.igBeginTooltip()
		ig.igText(self:getCoordText())
		ig.igText(self.mouseCoordValue)
		ig.igEndTooltip()
	end
end

HydroCLApp.showMouseCoords = cmdline.showMouseCoords
--if HydroCLApp.showMouseCoords == nil then HydroCLApp.showMouseCoords = true end
HydroCLApp.showMouseCoords = not not HydroCLApp.showMouseCoords

HydroCLApp.mouseCoord = {0,0}
HydroCLApp.mouseCoordValue = ''	-- TODO store one per inst of App

function HydroCLApp:getCoordText()
	return ('%f, %f'):format(self.mouseCoord[1], self.mouseCoord[2])
end


local leftShiftDown, rightShiftDown
local leftGuiDown, rightGuiDown
function HydroCLApp:event(event, ...)
	if HydroCLApp.super.event then
		HydroCLApp.super.event(self, event, ...)
	end
	local shiftDown = leftShiftDown or rightShiftDown
	local guiDown = leftGuiDown or rightGuiDown
	if event.type == sdl.SDL_MOUSEMOTION then
		if canHandleMouse() then
			local dx = event.motion.xrel
			local dy = event.motion.yrel
			if dx ~= 0 or dy ~= 0 then
				if mouse.leftDown and not guiDown then
					if self.mouse_influenceEquations then
						-- 1) determine grid coord
						-- 2) tell the eqn to mess with the eqn at that location
					else
						if shiftDown then
							self.view:mouseZoom(-dy, dy)
						else
							self.view:mousePan(dx, dy, self:size())
						end
					end
				end
			end
		end
	elseif event.type == sdl.SDL_KEYDOWN then
		if event.key.keysym.sym == sdl.SDLK_LSHIFT then
			leftShiftDown = true
		elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
			rightShiftDown = true
		elseif event.key.keysym.sym == sdl.SDLK_LGUI then
			leftGuiDown = true
		elseif event.key.keysym.sym == sdl.SDLK_RGUI then
			rightGuiDown = true
		end
	elseif event.type == sdl.SDL_KEYUP then
		if event.key.keysym.sym == sdl.SDLK_LSHIFT then
			leftShiftDown = false
		elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
			rightShiftDown = false
		elseif event.key.keysym.sym == sdl.SDLK_LGUI then
			leftGuiDown = false
		elseif event.key.keysym.sym == sdl.SDLK_RGUI then
			rightGuiDown = false
		elseif canHandleKeyboard() then
			if event.key.keysym.sym == sdl.SDLK_SPACE then
				self.running = not self.running
			elseif event.key.keysym.sym == ('u'):byte() then
				self.running = 'step'
			elseif event.key.keysym.sym == ('r'):byte() then
				print'resetting...'
				for _,solver in ipairs(self.solvers) do
					solver:resetState()
				end
				self.running = false
			elseif event.key.keysym.sym == ('p'):byte() then
				if shiftDown then
					self:resetGradientTex()
				else
					self:randomizeGradientTex()
				end
			end
		end
	end
end

--[[
this was a function of its own, but i want it to be per-App (rather than per-Lua)
so I'm merging it with App...

for unique names I'm using two systems:
1) for code, used for caching cl binaries, where the object uid is used
2) for types, used for ffi.cdef (which can't redefine types), using this incrementing suffix,

now because of constraint #2, if you run the same program twice within the same script,
you will get two different cache files in #1 ...

which makes me think that I should just say "FOR THE SAKE OF LUAJIT FFI TYPEDEFS, DON'T RUN THE SAME SIMULATION TWICE IN THE SAME PROCESS, UNLESS YOU WANT THE CL BINARY CACHE TO BREAK"

in fact, when running long-term tests (like test-order), I'm even getting "table overflow" errors ... that probably pertain to all the unique names that I'm keeping track of
--]]
function HydroCLApp:uniqueName(name)
	self.allnames = self.allnames or {}
	-- don't use the base name because I'm using the base name in my typedefs per-source file
	for i=1,math.huge do
		local try = name..'_'..i
		if not self.allnames[try] then
			self.allnames[try] = true
			return try
		end
	end
end

return HydroCLApp
