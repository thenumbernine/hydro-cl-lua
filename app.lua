--[[
predefined vars:
	dim =
	gridSize =
	fluxLimiter =
	boundary =
	integrator =
	eqn =
	mins =
	maxs =
	float = set to true to use 32 bit float instead of 64 
	half = set to true to use 16 bit float instead of 64
	cpu = set to use CPU instead of GPU
	useGLSharing = set to false to disable GL sharing
	disableGUI = set to disable GUI and prevent loading of imgui altogether
	disableFont = set to disable loading of the font.png file
--]]
local cmdline = {}

local fromlua = require 'ext.fromlua'
for _,w in ipairs(arg or {}) do
	local k,v = w:match'^(.-)=(.*)$'
	if k then
		cmdline[k] = fromlua(v)
	else
		cmdline[w] = true
	end
end

-- if we are disabling the gui then replace the imgui and tooltip requires, so we don't try to unnecessarily load it
if cmdline.disableGUI then
	package.loaded['ffi.imgui'] = {disabled=true}
	package.tooltip = {disabled=true}
end

local bit = require 'bit'
local ffi = require 'ffi'
local cl = require 'ffi.OpenCL'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local math = require 'ext.math'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
local template = require 'template'
local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'
local GLProgram = require 'gl.program'
local GLGradientTex = require 'gl.gradienttex'
local GLTex2D = require 'gl.tex2d'
local Font = require 'gui.font'
local vec4d = require 'ffi.vec.vec4d'
local vec3d = require 'ffi.vec.vec3d'
local tooltip = require 'tooltip'

-- I tried making this a flag, and simply skipping the gui update if it wasn't set, but imgui still messes with the GL state and textures and stuff
--  and I still get errors... so I'm cutting out imgui altogether, but now it takes a global flag to do so.
local HydroCLApp = class(cmdline.disableGUI and require 'glapp' or require 'imguiapp')

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
	{name='minmod', code='return max(0., min(r, 1.));'},
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
HydroCLApp.limiterNames = HydroCLApp.limiters:map(function(limiter) return limiter.name end)

--[[
setup for the solver
args has: platformName, deviceName
override this for specific experiments
this can't override float vs double precision yet

make sure to call this after 'real' is defined
 specifically the call to 'refreshGridSize' within it
--]]
function HydroCLApp:setup(args)
	args.self = self
	args.cmdline = cmdline
	args.table = table
	local keys = table.keys(args)
	assert(load([[
local ]]..keys:concat', '..[[ = ...
]] .. file['config.lua']))(
	keys:map(function(key) return args[key] end):unpack()
)
end

local useClipPlanes
if useClipPlanes then
-- TODO put all the display3D_slices stuff in its own file
local rotateClip = ffi.new('int[1]', 0)
local function makeDefaultPlane(i)
	assert(i >= 1 and i <= 4)
	local plane = vec4d(0,0,0,0)
	plane:ptr()[math.min(i,3)-1] = -1
	return plane
end
local clipInfos = range(4):map(function(i)
	return {
		enabled = i == 3,
		plane = makeDefaultPlane(i),
	}
end)
end

-- needs to go before display2DMethods
require 'draw.2d_heatmap'(HydroCLApp)
require 'draw.2d_graph'(HydroCLApp)

-- needs to go before initGL
local display2DMethods = table{
	{Heatmap = HydroCLApp.display2D_Heatmap},
	{Graph = HydroCLApp.display2D_Graph},
}

require 'draw.3d_slice'(HydroCLApp)
require 'draw.3d_ray'(HydroCLApp)
require 'draw.3d_iso'(HydroCLApp)

local display3DMethods = table{
	{Slices = HydroCLApp.display3D_Slice},
	{Raytrace = HydroCLApp.display3D_Ray},
	{Isosurfaces = HydroCLApp.display3D_Isosurface},
}
local display3DMethodNames =  display3DMethods:map(function(kv)
	return (next(kv))
end)


--[[ Cheap output of the state each frame so I can compare it to other solvers.
-- 	This is what the dump-to-file is supposed to also do.
local printStateFile = io.open('out.txt', 'w')
local function printState(solver)
	local ptr = solver.UBufObj:toCPU()
	
	local cols = {0, 1, 4}
	
	for i=0,solver.numCells-1 do
		-- matching the cl defs:
		local dx = (solver.maxs[1] - solver.mins[1]) / (tonumber(solver.gridSize.x) - (2*solver.numGhost))
		local x = (i + .5 - solver.numGhost) * dx + solver.mins[1]
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

function HydroCLApp:initGL(...)
	if HydroCLApp.super.initGL then
		HydroCLApp.super.initGL(self, ...)
	end

	-- TODO favor cl_khr_gl_sharing, cl_khr_fp64, cl_khr_3d_image_writes
	self.env = CLEnv{
		verbose = true,
		precision = cmdline.float and 'float' or (cmdline.half and 'half' or nil),
		cpu = cmdline.cpu,
		useGLSharing = cmdline.useGLSharing ~= false,	-- let nil default to true 
	}
	local platformName = self.env.platform:getName()
	local deviceName = self.env.device:getName()
	print(platformName)
	print(deviceName)

	self.is64bit = self.env.real == 'double'
	self.useGLSharing = self.env.useGLSharing
	self.device = self.env.device
	self.ctx = self.env.ctx
	self.cmds = self.env.cmds
	self.real = self.env.real

	ffi.cdef('typedef '..self.real..' real;')

	do
		local code = template(file['math.h'])
		xpcall(function()
			ffi.cdef(code)
		end, function(err)
			print(require 'template.showcode'(code))
			error(err)
		end)
	end

	self.solvers = table()

	self:setup{platformName=platformName, deviceName=deviceName}


	-- This only looks good when overlaying vector fields on top of other graphs.
	-- When it comes to separate variables, they usually look better apart.
	self.displayAllTogether = self.solvers[1] and self.solvers[1].dim > 1 or false


	self.gradientTex = GLGradientTex(1024, {
	-- [[ white, rainbow, black
		{0,0,0,.5},	-- black ... ?
		{0,0,1,1},	-- blue
		{0,1,1,1},	-- cyan
		{0,1,0,1},	-- green
		{1,1,0,1},	-- yellow
		{1,.5,0,1},	-- orange
		{1,0,0,1},	-- red
		{1,1,1,1},	-- white
	--]]
	--[[ stripes 
		range(32):map(function(i)
			return ({
				{0,0,0,0},
				{1,1,1,1},
			})[i%2+1]
		end):unpack()
	--]]
	}, false)
	-- don't wrap the colors, but do use GL_REPEAT
	self.gradientTex:setWrap{s = gl.GL_REPEAT}



	-- this will be per-solver
	-- but is also tightly linked to the structured grid solvers



	local graphShaderCode = file['draw/graph.shader']
	self.graphShader = GLProgram{
		vertexCode = '#define VERTEX_SHADER\n'..graphShaderCode,
		fragmentCode = '#define FRAGMENT_SHADER\n'..graphShaderCode,
		uniforms = {
			tex = 0,
			scale = 1,
			ambient = 1,
		},
	}

	-- init all resources for all draw methods, so the user can switch between methods quickly
	for _,solver in ipairs(self.solvers) do
		solver:initDraw()
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
			gl.glTexParameteri(fonttex.target, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST)
			gl.glTexParameteri(fonttex.target, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR) 
		end
		self.font = Font{tex = fonttex}
	end

	-- todo reorganize me
	self.display2DMethodsEnabled = display2DMethods:map(function(method, index)
		local name, func = next(method)
		return index == 1, name
	end)

	self.orthoView = require 'view.ortho'()
	self.frustumView = require 'view.frustum'()
	self.view = (#self.solvers > 0 and self.solvers[1].dim == 3) and self.frustumView or self.orthoView


if printState then
	for _,solver in ipairs(self.solvers) do
		printState(solver)
	end
end

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

	self.mins = vec3(-1,-1,-1)
	self.maxs = vec3(1,1,1)
	
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
		f = io.open('var-ranges.txt', 'w')
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
				f:write(('\t%.16e'):format(ymin))
				f:write(('\t%.16e'):format(ymax))
			end
		end
		--]]
		-- [[ gmres error
		f:write('\t'..solver.integrator.last_err)
		f:write('\t'..solver.integrator.last_iter)
		--]]
	end
	f:write'\n'
	f:flush()
end


-- TODO split off plot and dump stuff into its own folder/files
local io = require 'ext.io'

-- dropdown options
HydroCLApp.screenshotExts = {'png', 'bmp', 'jpeg', 'tiff', 'fits', 'tga', 'ppm'}
-- dropdown index
HydroCLApp.screenshotExtIndex = 1

function HydroCLApp:screenshot()
	local ext = self.screenshotExts[self.screenshotExtIndex]

	-- TODO only once upon init?
	if not io.fileexists'screenshots' then
		-- don't assert -- if it already exists the cmd will fail
		os.execute'mkdir screenshots'
	end

	-- make a new subdir for each application instance ... ?
	if not self.screenshotDir then
		self.screenshotDir = os.date('%Y.%m.%d-%H.%M.%S')
		local dir = 'screenshots/'..self.screenshotDir
		assert(not io.fileexists(dir), "found a duplicate screenshot timestamp subdir")
	
		-- bleh, windows.
		-- TODO a mkdir for everyone, either in ext.file or get lfs working with luajit 
		if ffi.os == 'Windows' then dir = dir:gsub('/', '\\') end
		assert(os.execute('mkdir '..dir))
		
		self.screenshotIndex = 0
	end

	local fn = ('screenshots/'..self.screenshotDir..'/%05d.'..ext):format(self.screenshotIndex)
	self.screenshotIndex = self.screenshotIndex + 1
	self:screenshotToFile(fn)
end

local Image = require 'image'
function HydroCLApp:screenshotToFile(fn)
	local w, h = self:size()
	if self.ssimg then
		if w ~= self.ssimg.width or h ~= self.ssimg.height then
			self.ssimg = nil
			self.ssflipped = nil
		end
	end
	if not self.ssimg then
		self.ssimg = Image(w, h, 3, 'unsigned char')
		self.ssflipped = Image(w, h, 3, 'unsigned char')
	end
	gl.glReadPixels(0, 0, w, h, gl.GL_RGB, gl.GL_UNSIGNED_BYTE, self.ssimg.buffer) 
	-- reverse rows ...
	-- TODO maybe ... for all projection matrix setups, have them check a screenshot flag and automatically flip?
	for y=0,h-1 do
		ffi.copy(
			self.ssflipped.buffer + (h-y-1) * w * 3,
			self.ssimg.buffer + y * w * 3,
			w * 3)
	end
	self.ssflipped:save(fn)
end

--HydroCLApp.running = false
HydroCLApp.running = true

local minDeltaY = 1e-7
function HydroCLApp:update(...)
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
			oldestSolver:update() 

if printState then
	for _,solver in ipairs(self.solvers) do
		printState(solver)
	end
end


			-- TODO should the time be oldestSolver.t after oldestSolver just updated?
			-- or - if dumpFile is enabled - should we re-search-out the oldest solver and use its time?
			dumpFile:update(self, oldestSolver.t)
		
			if self.exitTime and oldestSolver.t > self.exitTime then
				self:requestExit()
			end
		end
	else	
		-- clear all 'lastFrameTime's of solvers so the rough fps calcs don't get messed with
		for _,solver in ipairs(self.solvers) do
			solver.lastFrameTime = nil
		end
	end

	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	
	local w, h = self:size()

	local varNamesEnabled = table()
	local varNamesEnabledByName = {}
	for _,solver in ipairs(self.solvers) do
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

	local useLog
	local vectorField
	for _,varName in ipairs(varNamesEnabled) do
		local xmin, xmax, ymin, ymax
		for _,solver in ipairs(self.solvers) do
			local var = solver.displayVarForName[varName]
			
			if var and var.enabled
			--and solver.visiblePtr and solver.visiblePtr[0] 
			then
				useLog = var.useLog
				vectorField = var.vectorField
				
				local solverxmin, solverxmax = solver.mins[1], solver.maxs[1]
				solverxmin, solverxmax = 1.1 * solverxmin - .1 * solverxmax, 1.1 * solverxmax - .1 * solverxmin
				if solver.dim > 1 then
					local center = .5 * (solverxmin + solverxmax)
					solverxmin = (solverxmin - center) * ar + center
					solverxmax = (solverxmax - center) * ar + center
				end

				local solverymin, solverymax
				if solver.dim > 1 then
					solverymin, solverymax = solver.mins[2], solver.maxs[2]
					solverymin, solverymax = 1.1 * solverymin - .1 * solverymax, 1.1 * solverymax - .1 * solverymin
				else
					solverymin, solverymax = solver:calcDisplayVarRange(var)
					
					if useLog then
						solverymin = math.log(solverymin, 10)
						solverymax = math.log(solverymax, 10)
						if not math.isfinite(solverymin) then solverymin = -math.huge end
						if not math.isfinite(solverymax) then solverymax = -math.huge end
						solverymin = math.max(-30, solverymin)
						solverymax = math.max(-30, solverymax)
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

		if not self.displayAllTogether then
			gl.glViewport(
				graphCol / graphsWide * w,
				(1 - (graphRow + 1) / graphsHigh) * h,
				w / graphsWide,
				h / graphsHigh)
		end

		-- TODO maybe find the first solver for this var and use it to choose 1D,2D,3D
		local dim = self.solvers[1].dim
		
		if not vectorField then
			if dim == 1 then
				self:display1D(self.solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
			elseif dim == 2 then
				self:display2D(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
			elseif dim == 3 then
				self:display3D(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
			end
		else
			--if self.enableVectorField then
			self:displayVectorField(self.solvers, ar, varName, xmin, ymin, xmax, ymax)
			--end
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
		local solverNames = self.solvers:map(function(solver)
			return {
				text = ('(%.3f) %s'):format(solver.t, solver.name),
				color = solver.color,
			}
		end)
		local fontSizeX = .02
		local fontSizeY = .02
		local maxlen = solverNames:map(function(solverName)
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
				color = {solverName.color[1], solverName.color[2], solverName.color[3],1},
				fontSize = {fontSizeX, -fontSizeY},
				multiLine = false,
			}
		end
	end

	-- screenshot before gui
	if self.createAnimation then 
		self:screenshot() 
		if self.createAnimation == 'once' then
			self.createAnimation = nil
		end
	end

	if HydroCLApp.super.update then
		HydroCLApp.super.update(self, ...)
	end
end

require 'draw.1d'(HydroCLApp)

function HydroCLApp:drawGradientLegend(ar, varName, valueMin, valueMax)
	self.orthoView:projection(ar)
	self.orthoView:modelview()
	local xmin, xmax, ymin, ymax = self.orthoView:getOrthoBounds(ar)
	
	if self.font then
		local palwidth = (xmax - xmin) * .1
		self.gradientTex:enable()
		self.gradientTex:bind()
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

function HydroCLApp:display2D(...)
	for _,method in ipairs(display2DMethods) do
		local name, func = next(method)
		if self.display2DMethodsEnabled[name] then
			func(self, ...)
		end
	end
end

function HydroCLApp:display3D(...)
	self.display3DMethod = self.display3DMethod or 1 
	select(2, next(display3DMethods[self.display3DMethod]))(self, ...)
end

require 'draw.vectorfield'.applyToApp(HydroCLApp)
--require 'draw.vectorfield2'(HydroCLApp)

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
			-- save as cfits 
			for i,solver in ipairs(self.solvers) do
				solver:save(tostring(i))
			end
		end
		
		-- dump min/max(/avg?) of displayvars to a .txt file
		tooltip.checkboxTable('dump to text file', dumpFile, 'enabled')

		if ig.igButton'Screenshot' then
			self:screenshot()
		end
		ig.igSameLine()

		tooltip.comboTable('screenshot ext', self, 'screenshotExtIndex', self.screenshotExts)

		if ig.igButton(self.createAnimation and 'stop frame dump' or 'start frame dump') then
			self.createAnimation = not self.createAnimation
		end

		tooltip.checkboxTable('stack graphs', self, 'displayAllTogether')

		if ig.igRadioButtonBool('ortho', self.view == self.orthoView) then
			self.view = self.orthoView
		end
		ig.igSameLine()
		if ig.igRadioButtonBool('frustum', self.view == self.frustumView) then
			self.view = self.frustumView
		end
	
		-- TODO flag for separate/combined displays (esp for ortho view)

		-- TODO flag to toggle slice vs volume display
		-- or maybe checkboxes for each kind?
		
		if self.solvers[1] then
			local dim = self.solvers[1].dim
			if dim == 2 then
				ig.igPushIDStr'2D'
				for i,method in ipairs(display2DMethods) do
					if i > 1 then ig.igSameLine() end
					local name, func = next(method)
					tooltip.checkboxTable(name, self.display2DMethodsEnabled, name)
				end
				
				if self.display2DMethodsEnabled.Graph then
					tooltip.intTable('graph step', self, 'display2D_Graph_step')
				end
				
				ig.igPopID()
			
			
			elseif dim == 3 then
				ig.igPushIDStr'3D'
				tooltip.comboTable('Display Method', self, 'display3DMethod', display3DMethodNames)
				
				-- if we're doing 3D slice display 
				if HydroCLApp.display3D_Slice == select(2, next(display3DMethods[self.display3DMethod])) then

if useClipPlanes then
					ig.igRadioButton("rotate camera", rotateClip, 0)
					for i,clipInfo in ipairs(clipInfos) do
						ig.igPushIDStr('clip '..i)
						tooltip.checkbox('clip', clipInfo, 'enabled')
						ig.igSameLine()
						ig.igRadioButton('rotate', rotateClip, i)
						ig.igSameLine()
						if ig.igButton('reset') then
							clipInfo.plane = makeDefaultPlane(i)
						end
						ig.igPopID()
					end				
end					
					tooltip.sliderTable('alpha', self, 'display3D_Slice_alpha', 0, 1)
					tooltip.sliderTable('gamma', self, 'display3D_Slice_alphaGamma', 0, 1)
					tooltip.checkboxTable('isobars', self, 'display3D_Slice_useIsos')
					if self.display3D_Slice_useIsos then
						tooltip.intTable('num isobars', self, 'display3D_Slice_numIsobars')
					end
					tooltip.checkboxTable('lighting', self, 'display3D_Slice_useLighting')
					tooltip.checkboxTable('pointcloud', self, 'display3D_Slice_usePoints')
					if not self.display3D_Slice_usePoints then
						tooltip.intTable('num slices', self, 'display3D_Slice_numSlices')
					end
				end
				ig.igPopID()
			end
		
			--ig.igCheckbox('vector field', self.enableVectorField)
		
			tooltip.numberTable('vector field scale', self, 'displayVectorField_scale')
			--tooltip.sliderTable('vector field scale', self, 'displayVectorField_scale', 0, 100, nil, 10)
			
			tooltip.intTable('vector field step', self, 'displayVectorField_step')
			self.displayVectorField_step = math.max(self.displayVectorField_step, 1)
		end
	end
	
	for i,solver in ipairs(self.solvers) do
		ig.igPushIDStr('solver '..i)
		if ig.igCollapsingHeader(solver.name) then
			-- TODO new window for each
			solver:updateGUI()
		end
		ig.igPopID()
	end
end

local leftButtonDown
local rightButtonDown
local leftShiftDown
local rightShiftDown
local leftGuiDown
local rightGuiDown
function HydroCLApp:event(event, ...)
	if HydroCLApp.super.event then
		HydroCLApp.super.event(self, event, ...)
	end
	local canHandleMouse = not rawget(ig, 'disabled') and not ig.igGetIO()[0].WantCaptureMouse
	local canHandleKeyboard = not rawget(ig, 'disabled') and not ig.igGetIO()[0].WantCaptureKeyboard
	local shiftDown = leftShiftDown or rightShiftDown
	local guiDown = leftGuiDown or rightGuiDown
	if event.type == sdl.SDL_MOUSEMOTION then
		if canHandleMouse then
			local dx = event.motion.xrel
			local dy = event.motion.yrel
			if leftButtonDown and not guiDown then
				if shiftDown then
					if dx ~= 0 or dy ~= 0 then
						self.view:mouseZoom(-dy, dy)
					end
				else
					if dx ~= 0 or dy ~= 0 then
						self.view:mousePan(dx, dy, self:size())
					end
				end
			end
		end
	elseif event.type == sdl.SDL_MOUSEBUTTONDOWN then
		if event.button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = true
		elseif event.button.button == sdl.SDL_BUTTON_RIGHT then
			rightButtonDown = true
		end
	elseif event.type == sdl.SDL_MOUSEBUTTONUP then
		if event.button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = false
		elseif event.button.button == sdl.SDL_BUTTON_RIGHT then
			rightButtonDown = false
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
		elseif canHandleKeyboard then
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
			end
		end
	end
end

return HydroCLApp
