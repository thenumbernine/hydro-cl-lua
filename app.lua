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
local tolua = require 'ext.tolua'
local ImGuiApp = require 'imguiapp'
local CLPlatform = require 'cl.platform'
local CLContext = require 'cl.context'
local CLCommandQueue = require 'cl.commandqueue'
local GLProgram = require 'gl.program'
local GLGradientTex = require 'gl.gradienttex'
local GLTex2D = require 'gl.tex2d'
local Font = require 'gui.font'
local RoeSolver = require 'solver.roe'
local SRHDRoeSolver = require 'solver.srhd-roe'
local Euler1DEqn = require 'eqn.euler1d'
local Euler3DEqn = require 'eqn.euler3d'
local MaxwellEqn = require 'eqn.maxwell'
local ADM1Dv1Eqn = require 'eqn.adm1d_v1'
local ADM1Dv2Eqn = require 'eqn.adm1d_v2'
local ADM3DEqn = require 'eqn.adm3d'

local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}

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

	-- create this after 'real' is defined
	--  specifically the call to 'refreshGridSize' within it
	local args = {
		app = self, 
		gridSize = {
			cmdline.gridSize or 256,
			cmdline.gridSize or 256,
			cmdline.gridSize or 256,
		},
		boundary = {
			xmin=cmdline.boundary or 'freeflow',
			xmax=cmdline.boundary or 'freeflow',
			ymin=cmdline.boundary or 'freeflow',
			ymax=cmdline.boundary or 'freeflow',
			zmin=cmdline.boundary or 'freeflow',
			zmax=cmdline.boundary or 'freeflow',
		},
		integrator = cmdline.integrator or 'forward Euler',	--'Runge-Kutta 4, TVD',
		slopeLimiter = cmdline.slopeLimiter or 'superbee',
		dim = cmdline.dim or 2,
		mins = cmdline.mins or {-1, -1, -1},
		maxs = cmdline.maxs or {1, 1, 1},
	}

	--[[
	self.solver = SRHDRoeSolver(table(args, {
		initState = 'relativistic blast wave test problem 2',
	}))
	--]]
	-- [[
	self.solver = RoeSolver(table(args, {
		--eqn = (cmdline.eqn and require('eqn.'..cmdline.eqn) or Euler3DEqn)(),
		-- fluids
		--eqn = Euler1DEqn(),
		eqn = Euler3DEqn(),
		-- electromagnetism
		--eqn = MaxwellEqn(),
		-- geometrodynamics
		--eqn = ADM1Dv1Eqn(),
		--eqn = ADM1Dv2Eqn(),
		--eqn = ADM3DEqn(),
	}))
	--]]

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
		vertexCode = '#define VERTEX_SHADER\n'..code,
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
		local code = file['volumetric.shader']
		self.volumetricShader = GLProgram{
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
		self.volumetricShader:use()
		gl.glUniform1i(self.volumetricShader.uniforms.tex, 0)
		gl.glUniform1i(self.volumetricShader.uniforms.gradient, 1)
		do
			local maxiter = math.max(tonumber(self.solver.gridSize.x), tonumber(self.solver.gridSize.y), tonumber(self.solver.gridSize.z))
			print('volumetric shader raycast maxiter',maxiter)
			gl.glUniform1i(self.volumetricShader.uniforms.maxiter, maxiter)
		end
		--gl.glUniform3f(self.volumetricShader.uniforms.oneOverDx, (self.solver.maxs - self.solver.mins):unpack())
		self.volumetricShader:useNone()
	end

	self.gradientTex = GLGradientTex(1024, {
		{0,0,.5,1},
		{0,0,1,1},
		{0,1,1,1},
		{1,1,0,1},
		{1,0,0,1},
	}, false)

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

	local ar = w / h

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
					solverymin, solverymax = solver:calcDisplayVarRange(varIndex)
		
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
			self:display1D(self.solvers, varName, xmin, ymin, xmax, ymax, useLog)
		elseif self.solver.dim == 2 then
			self:display2D(self.solvers, varName, xmin, ymin, xmax, ymax)
		elseif self.solver.dim == 3 then
			self:display3D(self.solvers, varName, xmin, ymin, xmax, ymax)
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

function HydroCLApp:display1D(solvers, varName, xmin, ymin, xmax, ymax, useLog)
	local w, h = self:size()
	
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

function HydroCLApp:display2D(solvers, varName, xmin, ymin, xmax, ymax)
	local w, h = self:size()

	-- [[ ortho
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	--]]

	--[[ frustum
	
	--]]

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
				valueMin, valueMax = solver:calcDisplayVarRange(varIndex)
				var.heatMapValueMinPtr[0] = valueMin
				var.heatMapValueMaxPtr[0] = valueMax
			end

			solver:calcDisplayVarToTex(varIndex, var)
	
			self.heatMap2DShader:use()
			gl.glUniform1i(self.heatMap2DShader.uniforms.useLog, var.useLogPtr[0])
			gl.glUniform1f(self.heatMap2DShader.uniforms.valueMin, valueMin)
			gl.glUniform1f(self.heatMap2DShader.uniforms.valueMax, valueMax)
			self.solver.tex:bind(0)
			self.gradientTex:bind(1)
			gl.glBegin(gl.GL_QUADS)
			for _,v in ipairs{{0,0},{1,0},{1,1},{0,1}} do
				gl.glTexCoord2f(v[1], v[2])
				gl.glVertex2f(
					v[1] * solver.maxs[1] + (1 - v[1]) * solver.mins[1],
					v[2] * solver.maxs[2] + (1 - v[2]) * solver.mins[2])
			end
			gl.glEnd()
			self.gradientTex:unbind(1)
			self.solver.tex:unbind(0)
			self.heatMap2DShader:useNone()

			if self.font then
				local fontSizeX = (xmax - xmin) * .025
				local fontSizeY = (ymax - ymin) * .025
				local ystep = ystep * 2
				for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
					local value = (y - ymin) * (valueMax - valueMin) / (ymax - ymin)
					self.font:draw{
						pos={xmin * .99 + xmax * .01, y + fontSizeY * .5},
						text=('%.5f'):format(value),
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

function HydroCLApp:display3D(solvers, varName, xmin, ymin, xmax, ymax, useLog)
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
	
	local w, h = self:size()

	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	local znear, zfar = .1, 1000
	gl.glFrustum(xmin * znear, xmax * znear, ymin * znear, ymax * znear, znear, zfar) 
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	gl.glTranslatef(0,0,-2)

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
			self.volumetricShader:use()
			gl.glUniform1f(self.volumetricShader.uniforms.scale, 1)--scale)
			gl.glUniform1i(self.volumetricShader.uniforms.useLog, useLog and 1 or 0)
			gl.glUniform1f(self.volumetricShader.uniforms.alpha, 1)--alpha)
			self.solver.tex:bind(0)
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
			self.solver.tex:unbind(0)
			self.volumetricShader:useNone()
			gl.glDisable(gl.GL_BLEND)
			gl.glDisable(gl.GL_DEPTH_TEST)
			gl.glCullFace(gl.GL_BACK)
			gl.glDisable(gl.GL_CULL_FACE)
		end
	end

end

function HydroCLApp:showDisplayVar(solver, varIndex)
	local var = solver.displayVars[varIndex]
	solver:calcDisplayVarToTex(varIndex, var)	
	-- display

	self.graphShader:use()
	solver.tex:bind()

	gl.glUniform1i(self.graphShader.uniforms.useLog, var.useLogPtr[0])
	gl.glUniform2f(self.graphShader.uniforms.xmin, solver.mins[1], 0)
	gl.glUniform2f(self.graphShader.uniforms.xmax, solver.maxs[1], 0)
	gl.glUniform1i(self.graphShader.uniforms.axis, solver.dim)
	gl.glUniform2f(self.graphShader.uniforms.size, solver.gridSize.x, solver.gridSize.y)

	gl.glColor3f(table.unpack(var.color))
	gl.glBegin(gl.GL_LINE_STRIP)
	local step = 1
	for i=2,tonumber(solver.gridSize.x)-2,step do
		local x = (i+.5)/tonumber(solver.gridSize.x)
		gl.glVertex2f(x, 0)
	end
	gl.glEnd()
	
	solver.tex:unbind()
	self.graphShader:useNone()
end

function HydroCLApp:updateGUI()
	if ig.igButton(self.updateMethod and 'Stop' or 'Start') then
		self.updateMethod = not self.updateMethod
	end
	if ig.igButton'Step' then
		self.updateMethod = 'step'
	end
	if ig.igButton'Reset' then
		print'resetting...'
		self.solver:resetState()
		self.updateMethod = nil
	end

	ig.igCheckbox('use fixed dt', self.solver.useFixedDT)
	ig.igInputFloat('fixed dt', self.solver.fixedDT)
	ig.igInputFloat('CFL', self.solver.cfl)

	if ig.igCombo('integrator', self.solver.integratorPtr, self.solver.integratorNames) then
		self.solver:refreshIntegrator()
	end

	if ig.igCombo('slope limiter', self.solver.slopeLimiterPtr, self.slopeLimiterNames) then
		self.solver:refreshSolverProgram()
	end

	for i=1,self.solver.dim do
		for _,minmax in ipairs(minmaxs) do
			local var = xs[i]..minmax
			if ig.igCombo(var, self.solver.boundaryMethods[var], self.boundaryMethods) then
				self.solver:refreshBoundaryProgram()
			end
		end
	end

	-- equation-specific:

	local eqn = self.solver.eqn

	if ig.igCombo('init state', self.solver.initStatePtr, eqn.initStateNames) then
		
		self.solver:refreshInitStateProgram()
		
		-- TODO changing the init state program might also change the boundary methods
		-- ... but I don't want it to change the settings for the running scheme (or do I?)
		-- ... but I don't want it to not change the settings ...
		-- so maybe refreshing the init state program should just refresh everything?
		-- or maybe just the boundaries too?
		-- hack for now:
		self.solver:refreshBoundaryProgram()
	end

	if ig.igCollapsingHeader'equation:' then
		local f = ffi.new'float[1]'
		local i = ffi.new'int[1]'
		for _,var in ipairs(eqn.guiVars) do
			var:updateGUI(self.solver)
			--[[
			elseif type(eqn_var) == 'table' then
			end
			--]]
		end
	end

	-- display vars: TODO graph vars

	if ig.igCollapsingHeader'variables:' then
		for _,convertToTex in ipairs(self.solver.convertToTexs) do
			if ig.igCollapsingHeader(convertToTex.name..' variables:') then
				for _,var in ipairs(convertToTex.vars) do
					ig.igPushIdStr(convertToTex.name..' '..var.name)
					ig.igCheckbox(var.name, var.enabled)
					ig.igSameLine()
					if ig.igCollapsingHeader'' then	
						ig.igCheckbox('log', var.useLogPtr)
						ig.igCheckbox('fixed range', var.heatMapFixedRangePtr)
						ig.igInputFloat('value min', var.heatMapValueMinPtr)
						ig.igInputFloat('value max', var.heatMapValueMaxPtr)
					end
					ig.igPopId()
				end
			end
		end
	end

	-- heat map var

	-- TODO volumetric var
end

function HydroCLApp:event(event, ...)
	HydroCLApp.super.event(self, event, ...)
	if ig.igGetIO()[0].WantCaptureKeyboard then return end
	if event.type == sdl.SDL_KEYDOWN then
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

return HydroCLApp
