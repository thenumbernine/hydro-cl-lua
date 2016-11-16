--[[
predefined vars:
	dim=
	gridSize=
	slopeLimiter=
	boundary=
	integrator=
	eqn=
	mins=
	maxs=
	float= set to false to override double precision
--]]
local cmdline = {}
--[[
for _,k in ipairs{
	'dim', 'gridSize', 'slopeLimiter', 'boundary', 
	'integrator', 'eqn', 'mins', 'maxs', 'float'
} do
	cmdline[k] = _G[k]
	_G[k] = nil
	if cmdline[k] == nil then 
		cmdline[k] = _G[k:lower()] 
		_G[k:lower()] = nil
	end
end
--]]
for _,w in ipairs(arg) do
	local k,v = w:match'^(.-)=(.*)$'
	if k then
		cmdline[k] = assert(loadstring('return '..v))()
	end
end

local bit = require 'bit'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local gl = require 'ffi.OpenGL'
local cl = require 'ffi.OpenCL'
local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local math = require 'ext.math'
local table = require 'ext.table'
local string = require 'ext.string'
local file = require 'ext.file'
local ImGuiApp = require 'imguiapp'
local CLPlatform = require 'cl.platform'
local CLContext = require 'cl.context'
local CLCommandQueue = require 'cl.commandqueue'
local GLProgram = require 'gl.program'
local GLGradientTex = require 'gl.gradienttex'
local GLTex2D = require 'gl.tex2d'
local Font = require 'gui.font'

local HydroCLApp = class(ImGuiApp)

HydroCLApp.title = 'Hydrodynamics in OpenCL'

HydroCLApp.boundaryMethods = table{'freeflow', 'periodic', 'mirror'}

-- list from https://en.wikipedia.org/wiki/Flux_limiter
HydroCLApp.slopeLimiters = table{
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

	--return (r + fabs(r)) / (1. + fabs(r));
	{name='van Leer', code='return max(0., r) * 2. / (1. + r);'},
	
	{name='monotized central', code='return max(0., min(2., min(.5 * (1. + r), 2. * r)));'},
	{name='superbee', code='return max((real)0., (real)max((real)min((real)1., (real)2. * r), (real)min((real)2., r)));'},
	{name='Barth-Jespersen', code='return .5 * (r + 1.) * min(1., min(4. * r / (r + 1.), 4. / (r + 1.)));'},
}
HydroCLApp.slopeLimiterNames = HydroCLApp.slopeLimiters:map(function(limiter) return limiter.name end)

function HydroCLApp:initGL(...)
	HydroCLApp.super.initGL(self, ...)

	-- TODO favor cl_khr_fp64, cl_khr_3d_image_writes, cl_khr_gl_sharing
-- [[
for i,platform in ipairs(CLPlatform.getAll()) do
	print()
	print('platform '..i)
	platform:printInfo()

	for j,device in ipairs(platform:getDevices()) do
		print()
		print('device '..j)
		device:printInfo()
	end
end
--]]	

	local function get64bit(list)
		local best = list:map(function(item)
			local exts = string.trim(item:getExtensions():lower())
			return {item=item, fp64=not not exts:match'cl_%w+_fp64'}
		end):sort(function(a,b)
			return (a.fp64 and 1 or 0) > (b.fp64 and 1 or 0)
		end)[1]
		return best.item, best.fp64
	end

	-- TODO favor cl_khr_fp64, cl_khr_3d_image_writes, cl_khr_gl_sharing

	if ffi.os == 'Windows' then
		self.platform = CLPlatform.getAll()[1]
		self.device = self.platform:getDevices{gpu=true}[1]
		self.is64bit = false--self.device:getExtensions():lower():match'cl_%2+_fp64'
	else
		self.platform = get64bit(CLPlatform.getAll())
		self.device, self.is64bit = get64bit(self.platform:getDevices{gpu=true})
	end
	
	-- cmd-line override
	if cmdline.float then
		self.is64bit = false
	end
	
self.device:printInfo()
	local exts = string.split(string.trim(self.device:getExtensions()):lower(),'%s+')
	self.useGLSharing = exts:find(nil, function(ext) return ext:match'cl_%w+_gl_sharing' end)
	
	self.ctx = CLContext{
		platform = self.platform,
		device = self.device,
		glSharing = self.useGLSharing,
	}
print()
self.ctx:printInfo()
print()
print('is 64 bit?',self.is64bit)
print()

	self.cmds = CLCommandQueue{context=self.ctx, device=self.device}

	-- making 'real' a parameter of solver would be clever
	-- but because it is using ffi.ctype it might be tough ...
	-- then again, any solver ffi.ctype defined will potentially collide with other solvers ...
	self.real = self.is64bit and 'double' or 'float'
	ffi.cdef('typedef '..self.real..' real;')

	self.real3TypeCode = [[
typedef union {
	real s[3];
	struct { real s0, s1, s2; };
	struct { real x, y, z; };
} real3;
]]
	ffi.cdef(self.real3TypeCode)

	self.real3Code = [[
#define _real3(a,b,c) (real3){.s={a,b,c}}

static inline real real3_dot(real3 a, real3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline real3 real3_scale(real3 a, real s) {
	return _real3(a.x * s, a.y * s, a.z * s);
}

static inline real3 real3_add(real3 a, real3 b) {
	return _real3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline real3 real3_sub(real3 a, real3 b) {
	return _real3(a.x - b.x, a.y - b.y, a.z - b.z);
}

]]

	-- create this after 'real' is defined
	--  specifically the call to 'refreshGridSize' within it
	local args = {
		app = self, 
		eqn = cmdline.eqn,
		dim = cmdline.dim or 2,
		
		integrator = cmdline.integrator or 'forward Euler',	
		--integrator = 'Runge-Kutta 4, TVD',
	
		-- this is a flux limiter for the record.  TODO implement slope limiter.
		slopeLimiter = cmdline.slopeLimiter or 'superbee',
		
		-- [[ cartesian
		geometry = 'cartesian',
		mins = cmdline.mins or {-1, -1, -1},
		maxs = cmdline.maxs or {1, 1, 1},
		gridSize = {
			cmdline.gridSize or 128,
			cmdline.gridSize or 128,
			cmdline.gridSize or 128,
		},
		boundary = {
			xmin=cmdline.boundary or 'freeflow',
			xmax=cmdline.boundary or 'freeflow',
			ymin=cmdline.boundary or 'freeflow',
			ymax=cmdline.boundary or 'freeflow',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]
		--[[ cylinder
		geometry = 'cylinder',
		mins = cmdline.mins or {.5, 0, -1},
		maxs = cmdline.maxs or {1, 2*math.pi, 1},
		gridSize = {
			cmdline.gridSize or 32,
			cmdline.gridSize or 256,
			cmdline.gridSize or 1,
		},
		boundary = {
			xmin=cmdline.boundary or 'mirror',
			xmax=cmdline.boundary or 'mirror',
			ymin=cmdline.boundary or 'periodic',
			ymax=cmdline.boundary or 'periodic',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]
		--[[ sphere
		geometry = 'sphere',
		mins = cmdline.mins or {0, -math.pi, .5},
		maxs = cmdline.maxs or {math.pi, math.pi, 1},
		gridSize = {
			cmdline.gridSize or 64,
			cmdline.gridSize or 128,
			cmdline.gridSize or 64,
		},
		boundary = {
			xmin=cmdline.boundary or 'periodic',
			xmax=cmdline.boundary or 'periodic',
			ymin=cmdline.boundary or 'periodic',
			ymax=cmdline.boundary or 'periodic',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		--]]

		-- no initial state means use the first
		-- initState = cmdline.initState,
		-- euler / srhd initial states:
		--initState = 'Sod',
		--initState = 'Sedov',
		--initState = 'relativistic blast wave test problem 2',
		--initState = 'Kelvin-Hemholtz',
	}

	-- fluid
	--self.solver = require 'solver.roe'(table(args, {eqn='euler1d'}))
	--self.solver = require 'solver.euler-roe'(args)
	self.solver = require 'solver.srhd-roe'(args)
	-- EM
	--self.solver = require 'solver.maxwell-roe'(args)
	-- EM+HD
	--self.solver = require 'solver.twofluid-emhd-roe'(args)
	-- geometrodynamics
	--self.solver = require 'solver.roe'(table(args, {eqn='adm1d_v1'}))
	--self.solver = require 'solver.roe'(table(args, {eqn='adm1d_v2'}))
	--self.solver = require 'solver.roe'(table(args, {eqn='adm3d'}))

	self.solvers = table{self.solver}
	
	local graphShaderCode = file['graph.shader']
	self.graphShader = GLProgram{
		vertexCode = '#define VERTEX_SHADER\n'..graphShaderCode,
		fragmentCode = '#define FRAGMENT_SHADER\n'..graphShaderCode,
		uniforms = {
			'tex', 
			'scale', 
			'xmin',
			'xmax', 
			'useLog', 
			'axis', 
			'ambient',
			'size', 
		},
	}
	self.graphShader:use()
	gl.glUniform1i(self.graphShader.uniforms.tex, 0)
	gl.glUniform1f(self.graphShader.uniforms.scale, 1)
	gl.glUniform1f(self.graphShader.uniforms.ambient, 1)	
	self.graphShader:useNone()

	local code = file['heatmap2d.shader']
	self.heatMap2DShader = GLProgram{
		vertexCode = table{
			'#define VERTEX_SHADER',
			self.solver:getCoordMapGLSLCode(),
			code,
		}:concat'\n',
		fragmentCode = '#define FRAGMENT_SHADER\n'..code,
		uniforms = {
			'useLog',
			'valueMin',
			'valueMax',
			'tex',
			'gradient',
			'alpha',
		},
	}
	self.heatMap2DShader:use()
	gl.glUniform1f(self.heatMap2DShader.uniforms.valueMin, 0)
	gl.glUniform1f(self.heatMap2DShader.uniforms.valueMax, 0)
	gl.glUniform1i(self.heatMap2DShader.uniforms.tex, 0)
	gl.glUniform1i(self.heatMap2DShader.uniforms.gradient, 1)
	gl.glUniform1f(self.heatMap2DShader.uniforms.alpha, 1)
	self.heatMap2DShader:useNone()

	if self.solver.dim == 3 then
		-- raytracing (stalling)
		
		local code = file['volumetric.shader']
		self.volumeRayShader = GLProgram{
			vertexCode = '#define VERTEX_SHADER\n'..code,
			fragmentCode = '#define FRAGMENT_SHADER\n'..code,
			uniforms = {
				'tex',
				'gradient',
				'maxiter',
	--			'oneOverDx',
				'scale',
				'useLog',
				'alpha',
			},
		}
		self.volumeRayShader:use()
		gl.glUniform1i(self.volumeRayShader.uniforms.tex, 0)
		gl.glUniform1i(self.volumeRayShader.uniforms.gradient, 1)
		do
			local maxiter = math.max(tonumber(self.solver.gridSize.x), tonumber(self.solver.gridSize.y), tonumber(self.solver.gridSize.z))
			print('volumetric shader raycast maxiter',maxiter)
			gl.glUniform1i(self.volumeRayShader.uniforms.maxiter, maxiter)
		end
		--gl.glUniform3f(self.volumeRayShader.uniforms.oneOverDx, (self.solver.maxs - self.solver.mins):unpack())
		self.volumeRayShader:useNone()

		-- volume slices

		self.volumeSliceShader = GLProgram{
		vertexCode = [[
varying vec3 pos;
void main() {
	pos = gl_Vertex.xyz;
	gl_Position = ftransform();
}
]],
			fragmentCode = [[
varying vec3 pos;
uniform sampler3D volTex;
uniform sampler2D hsvTex;
uniform vec3 normal;
uniform float alpha;
uniform float alphaGamma;
void main() {

	vec4 worldPos = gl_ModelViewMatrix * vec4(pos,1.);

	float value = texture3D(volTex, pos).r;
	vec4 voxelColor = vec4(texture2D(hsvTex, vec2(value, .5)).rgb, pow(alpha, alphaGamma));
	
	//calculate normal in screen coordinates
	vec4 n = gl_ModelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	gl_FragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}
]],
			uniforms = {
				'volTex',
				'hsvTex',
				'normal',
				'alpha',
				'alphaGamma',
			},
		}
		self.volumeSliceShader:use()
		gl.glUniform1i(self.volumeSliceShader.uniforms.volTex, 0)
		gl.glUniform1i(self.volumeSliceShader.uniforms.hsvTex, 1)
		self.volumeSliceShader:useNone()
	end

	self.gradientTex = GLGradientTex(1024, {
		{0,0,.5,1},
		{0,0,1,1},
		{0,1,1,1},
		{1,1,0,1},
		{1,0,0,1},
	}, true)

	-- [[ need to get image loading working
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
	--]]

	self.orthoView = require 'view.ortho'()
	self.frustumView = require 'view.frustum'()
	self.view = self.orthoView
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


HydroCLApp.updateMethod = nil

function HydroCLApp:update(...)
	if self.updateMethod then
		if self.updateMethod == 'step' then 
			print('performing single step...')
			self.updateMethod = nil 
		end

		self.solver:update()
	end

	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	
	local w, h = self:size()

	local varNamesEnabled = table()
	for i,var in ipairs(self.solver.displayVars) do
		if var.enabled[0] then
			varNamesEnabled:insert(var.name)
		end
	end

	local graphsWide = math.ceil(math.sqrt(#varNamesEnabled))
	local graphsHigh = math.ceil(#varNamesEnabled/graphsWide)
	local graphCol = 0
	local graphRow = 0
	
	local ar = (w / graphsWide) / (h / graphsHigh)

	local useLog
	for _,varName in ipairs(varNamesEnabled) do
		local xmin, xmax, ymin, ymax
		for _,solver in ipairs(self.solvers) do
			local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
			
			useLog = var.useLogPtr[0]
			
			if varIndex
			--and solver.visiblePtr and solver.visiblePtr[0] 
			then
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
						solverymax = math.max(solverymax, solverymin + 1e-5)
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
							local minDeltaY = 1e-5
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

		gl.glViewport(
			graphCol / graphsWide * w,
			(1 - (graphRow + 1) / graphsHigh) * h,
			w / graphsWide,
			h / graphsHigh)
		
		if self.solver.dim == 1 then
			self:display1D(self.solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
		elseif self.solver.dim == 2 then
			self:display2D(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
		elseif self.solver.dim == 3 then
			self:display3D(self.solvers, varName, ar, xmin, ymin, xmax, ymax)
		end

		graphCol = graphCol + 1
		if graphCol == graphsWide then
			graphCol = 0
			graphRow = graphRow + 1
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
				text = ('(%.3f) '):format(solver.t)..solver.name,
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

	HydroCLApp.super.update(self, ...)
end

function HydroCLApp:display1D(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex2f(x*xstep,ymin)
		gl.glVertex2f(x*xstep,ymax)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex2f(xmin,y*ystep)
		gl.glVertex2f(xmax,y*ystep)
	end
	gl.glEnd()

	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex2f(xmin, 0)
	gl.glVertex2f(xmax, 0)
	gl.glVertex2f(0, ymin)
	gl.glVertex2f(0, ymax)
	gl.glEnd()

	-- display here
	for _,solver in ipairs(solvers) do
		local varIndex = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex then
			self:showDisplayVar(solver, varIndex)
		end

		if self.font then
			local fontSizeX = (xmax - xmin) * .05
			local fontSizeY = (ymax - ymin) * .05
			local ystep = ystep * 2
			for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
				local realY = useLog and 10^y or y
				self.font:draw{
					pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
					text=tostring(realY),
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
			self.font:draw{
				pos={xmin, ymax},
				text=varName,
				color = {1,1,1,1},
				fontSize={fontSizeX, -fontSizeY},
				multiLine=false,
			}
		end
	end
end

function HydroCLApp:display2D(solvers, varName, ar, graph_xmin, graph_ymin, graph_xmax, graph_ymax)
	self.view:projection(ar)
	self.view:modelview()
	if self.view.getOrthoBounds then
		xmin, xmax, ymin, ymax = self.view:getOrthoBounds(ar)
	else
		xmin, xmax, ymin, ymax = graph_xmin, graph_ymin, graph_xmax, graph_ymax
	end

	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex2f(x*xstep,ymin)
		gl.glVertex2f(x*xstep,ymax)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex2f(xmin,y*ystep)
		gl.glVertex2f(xmax,y*ystep)
	end
	gl.glEnd()
	
	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex2f(xmin, 0)
	gl.glVertex2f(xmax, 0)
	gl.glVertex2f(0, ymin)
	gl.glVertex2f(0, ymax)
	gl.glEnd()

	-- NOTICE overlays of multiple solvers won't be helpful.  It'll just draw over the last solver.
	-- I've got to rethink the visualization
	for _,solver in ipairs(solvers) do 
		local varIndex, var = solver.displayVars:find(nil, function(var) return var.name == varName end)
		if varIndex then
			-- TODO allow a fixed, manual colormap range
			local valueMin, valueMax
			if var.heatMapFixedRangePtr[0] then
				valueMin = var.heatMapValueMinPtr[0]
				valueMax = var.heatMapValueMaxPtr[0]
			else
				valueMin, valueMax = solver:calcDisplayVarRange(var)
				var.heatMapValueMinPtr[0] = valueMin
				var.heatMapValueMaxPtr[0] = valueMax
			end

			solver:calcDisplayVarToTex(var)
	
			self.heatMap2DShader:use()
			gl.glUniform1i(self.heatMap2DShader.uniforms.useLog, var.useLogPtr[0])
			gl.glUniform1f(self.heatMap2DShader.uniforms.valueMin, valueMin)
			gl.glUniform1f(self.heatMap2DShader.uniforms.valueMax, valueMax)
			self.solver:getTex(var):bind(0)
			self.gradientTex:bind(1)
		
			local gridScale = 4
			local udivs = math.ceil(tonumber(solver.gridSize.x)/gridScale)
			local vdivs = math.ceil(tonumber(solver.gridSize.y)/gridScale)
			for vbase=0,vdivs-1 do
				gl.glBegin(gl.GL_TRIANGLE_STRIP)
				for ui=0,udivs do
					for vofs=1,0,-1 do
						local vi = vbase+vofs
						local u = ui/udivs
						local v = vi/vdivs
						gl.glTexCoord2d(
							(u * (tonumber(solver.gridSize.x) - 2 * solver.numGhost) + solver.numGhost) / tonumber(solver.gridSize.x),
							(v * (tonumber(solver.gridSize.y) - 2 * solver.numGhost) + solver.numGhost) / tonumber(solver.gridSize.y))
						gl.glVertex2d(
							u * solver.maxs[1] + (1 - u) * solver.mins[1],
							v * solver.maxs[2] + (1 - v) * solver.mins[2])
					end
				end
				gl.glEnd()
			end
			
			self.gradientTex:unbind(1)
			self.solver:getTex(var):unbind(0)
			self.heatMap2DShader:useNone()

			if self.font then
				local fontSizeX = (xmax - xmin) * .025
				local fontSizeY = (ymax - ymin) * .025
				local ystep = 10^(math.log(ymax - ymin, 10) - 1.5)
				for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
					local value = (y - ymin) * (valueMax - valueMin) / (ymax - ymin)
					self.font:draw{
						pos={xmin * .99 + xmax * .01, y + fontSizeY * .5},
						text=(math.abs(value) > 1e+5 and ('%.5e'):format(value) or ('%.5f'):format(value)),
						color = {1,1,1,1},
						fontSize={fontSizeX, -fontSizeY},
						multiLine=false,
					}
				end
				self.font:draw{
					pos={xmin, ymax},
					text=varName,
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
		end
	end
end

local quat = require 'vec.quat'
local viewAngle = quat()
local vec4d = require 'ffi.vec.vec4d'
function HydroCLApp:display3D_Slice(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)

	self.view:projection(ar)
	self.view:modelview()
	
	self.volumeSliceShader:use()
	self.solver:getTex(var):bind(0)
	self.gradientTex:bind(1)
	gl.glUniform1f(self.volumeSliceShader.uniforms.alpha, .15)
	gl.glUniform1f(self.volumeSliceShader.uniforms.alphaGamma, 1)

	gl.glEnable(gl.GL_TEXTURE_GEN_S)
	gl.glEnable(gl.GL_TEXTURE_GEN_T)
	gl.glEnable(gl.GL_TEXTURE_GEN_R)
	gl.glTexGeni(gl.GL_S, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_T, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_R, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGendv(gl.GL_S, gl.GL_OBJECT_PLANE, vec4d(1,0,0,0):ptr())
	gl.glTexGendv(gl.GL_T, gl.GL_OBJECT_PLANE, vec4d(0,1,0,0):ptr())
	gl.glTexGendv(gl.GL_R, gl.GL_OBJECT_PLANE, vec4d(0,0,1,0):ptr())

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	--[[ points
	gl.glPointSize(2)
	gl.glBegin(gl.GL_POINTS)
	for _,pt in ipairs(self.pts) do
		gl.glVertex3d( 
			(pt[1] - .5)/(self.max[1] + 1),
			(pt[2] - .5)/(self.max[2] + 1),
			(pt[3] - .5)/(self.max[3] + 1))
	end
	gl.glEnd()
	--]]
	-- [[ slices
	local n = 255
	local fwd = -viewAngle:zAxis()
	local fwddir = select(2, table(fwd):map(math.abs):sup())
	local quad = {{0,0},{1,0},{1,1},{0,1}}
	local jmin, jmax, jdir
	if fwd[fwddir] < 0 then
		jmin, jmax, jdir = 0, n, 1
	else
		jmin, jmax, jdir = n, 0, -1
	end
	gl.glUniform3f(self.volumeSliceShader.uniforms.normal, fwddir==1 and jdir or 0, fwddir==2 and jdir or 0, fwddir==3 and jdir or 0)
	
	gl.glBegin(gl.GL_QUADS)
	for j=jmin,jmax,jdir do
		local f = j/n
		for _,vtx in ipairs(quad) do
			if fwddir == 1 then
				gl.glVertex3f(f, vtx[1], vtx[2])
			elseif fwddir == 2 then
				gl.glVertex3f(vtx[1], f, vtx[2])
			elseif fwddir == 3 then
				gl.glVertex3f(vtx[1], vtx[2], f)
			end
		end
	end
	gl.glEnd()
	--]]
	
	gl.glDisable(gl.GL_BLEND)

	gl.glDisable(gl.GL_TEXTURE_GEN_S)
	gl.glDisable(gl.GL_TEXTURE_GEN_T)
	gl.glDisable(gl.GL_TEXTURE_GEN_R)

	self.gradientTex:unbind(1)
	self.solver:getTex(var):unbind(0)
	self.volumeSliceShader:useNone()
end

function HydroCLApp:display3D_Ray(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

	local vertexes = {
		0,0,0,
		1,0,0,
		0,1,0,
		1,1,0,
		0,0,1,
		1,0,1,
		0,1,1,
		1,1,1,
	}

	local quads = {
		0,1,3,2,
		4,6,7,5,
		1,5,7,3,
		0,2,6,4,
		0,4,5,1,
		2,3,7,6,
	}
	

	gl.glColor3f(1,1,1)
	for pass=0,1 do
		if pass == 0 then
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
		else
			gl.glEnable(gl.GL_CULL_FACE)
			gl.glCullFace(gl.GL_FRONT)
			gl.glEnable(gl.GL_DEPTH_TEST)
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
			gl.glEnable(gl.GL_BLEND)
			self.volumeRayShader:use()
			gl.glUniform1f(self.volumeRayShader.uniforms.scale, 1)--scale)
			gl.glUniform1i(self.volumeRayShader.uniforms.useLog, useLog and 1 or 0)
			gl.glUniform1f(self.volumeRayShader.uniforms.alpha, 1)--alpha)
			self.solver:getTex(var):bind(0)
			self.gradientTex:bind(1)
		end
		gl.glBegin(gl.GL_QUADS)
		for i=1,24 do
			local x = vertexes[quads[i] * 3 + 0 + 1]
			local y = vertexes[quads[i] * 3 + 1 + 1]
			local z = vertexes[quads[i] * 3 + 2 + 1]
			gl.glTexCoord3f(x, y, z)
			x = x * (self.solver.maxs[1] - self.solver.mins[1]) + self.solver.mins[1]
			y = y * (self.solver.maxs[2] - self.solver.mins[2]) + self.solver.mins[2]
			z = z * (self.solver.maxs[3] - self.solver.mins[3]) + self.solver.mins[3]
			gl.glVertex3f(x, y, z)
		end
		gl.glEnd()
		if pass == 0 then
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
		else
			self.gradientTex:unbind(1)
			self.solver:getTex(var):unbind(0)
			self.volumeRayShader:useNone()
			gl.glDisable(gl.GL_BLEND)
			gl.glDisable(gl.GL_DEPTH_TEST)
			gl.glCullFace(gl.GL_BACK)
			gl.glDisable(gl.GL_CULL_FACE)
		end
	end
end

HydroCLApp.display3D = HydroCLApp.display3D_Slice
HydroCLApp.display3D = HydroCLApp.display3D_Ray

function HydroCLApp:showDisplayVar(solver, varIndex)
	local var = solver.displayVars[varIndex]
	solver:calcDisplayVarToTex(var)	
	-- display

	self.graphShader:use()
	solver:getTex(var):bind()

	gl.glUniform1i(self.graphShader.uniforms.useLog, var.useLogPtr[0])
	gl.glUniform2f(self.graphShader.uniforms.xmin, solver.mins[1], 0)
	gl.glUniform2f(self.graphShader.uniforms.xmax, solver.maxs[1], 0)
	gl.glUniform1i(self.graphShader.uniforms.axis, solver.dim)
	gl.glUniform2f(self.graphShader.uniforms.size, solver.gridSize.x, solver.gridSize.y)

	gl.glColor3f(table.unpack((#self.solvers > 1 and solver or var).color))
	gl.glBegin(gl.GL_LINE_STRIP)
	local step = 1
	for i=2,tonumber(solver.gridSize.x)-2,step do
		local x = (i+.5)/tonumber(solver.gridSize.x)
		gl.glVertex2f(x, 0)
	end
	gl.glEnd()
	
	solver:getTex(var):unbind()
	self.graphShader:useNone()
end

function HydroCLApp:updateGUI()
	if ig.igCollapsingHeader'simulation' then
		if ig.igButton(self.updateMethod and 'Stop' or 'Start') then
			self.updateMethod = not self.updateMethod
		end
		ig.igSameLine()
		if ig.igButton'Step' then
			self.updateMethod = 'step'
		end
		ig.igSameLine()
		if ig.igButton'Reset' then
			print'resetting...'
			self.solver:resetState()
			self.updateMethod = nil
		end
	
		if ig.igRadioButtonBool('ortho', self.view == self.orthoView) then
			self.view = self.orthoView
		end
		ig.igSameLine()
		if ig.igRadioButtonBool('frustum', self.view == self.frustumView) then
			self.view = self.frustumView
		end
	end
	
	for i,solver in ipairs(self.solvers) do
		ig.igPushIdStr('solver '..i)
		if ig.igCollapsingHeader(solver.name) then
			-- TODO new window for each
			solver:updateGUI()
		end
		ig.igPopId()
	end
end

local leftButtonDown
local rightButtonDown
local leftShiftDown
local rightShiftDown
local leftGuiDown
local rightGuiDown
function HydroCLApp:event(event, ...)
	HydroCLApp.super.event(self, event, ...)
	local canHandleMouse = not ig.igGetIO()[0].WantCaptureMouse
	local canHandleKeyboard = not ig.igGetIO()[0].WantCaptureKeyboard
	local shiftDown = leftShiftDown or rightShiftDown
	local guiDown = leftGuiDown or rightGuiDown
	if event.type == sdl.SDL_MOUSEMOTION then
		if canHandleMouse then
			local dx = event.motion.xrel
			local dy = event.motion.yrel
			if leftButtonDown and not guiDown then
				if shiftDown then
					if dx ~= 0 or dy ~= 0 then
						self.view:mouseZoom(dx, dy)
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
				self.updateMethod = not self.updateMethod
			elseif event.key.keysym.sym == ('u'):byte() then
				self.updateMethod = 'step'
			elseif event.key.keysym.sym == ('r'):byte() then
				print'resetting...'
				self.solver:resetState()
				self.updateMethod = nil
			end
		end
	end
end

return HydroCLApp
