#!/usr/bin/env luajit
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'

local ImGuiApp = require 'imguiapp'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local gl = require 'ffi.OpenGL'
local cl = require 'ffi.OpenCL'
local sdl = require 'ffi.sdl'
local bit = require 'bit'
local glreport = require 'gl.report'

local Solver = require 'solver'

local vec3sz = require 'vec3sz'

local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}

local HydroCLApp = class(ImGuiApp)

HydroCLApp.title = 'Hydrodynamics in OpenCL'

HydroCLApp.boundaryMethods = table{'freeflow', 'periodic', 'mirror'}

HydroCLApp.slopeLimiters = table{
	{name='DonorCell', code='return 0.;'},
	{name='LaxWendroff', code='return 1.;'},
	{name='BeamWarming', code='return r;'},
	{name='Fromm', code='return .5 * (1. + r);'},
	{name='CHARM', code='return max(0., r) * (3. * r + 1.) / ((r + 1.) * (r + 1.));'},
	{name='HCUS', code='return max(0., 1.5 * (r + fabs(r)) / (r + 2.));'},
	{name='HQUICK', code='return max(0., 2. * (r + fabs(r)) / (r + 3.));'},
	{name='Koren', code='return max(0., min(2. * r, min((1. + 2. * r) / 3., 2.)));'},
	{name='MinMod', code='return max(0., min(r, 1.));'},
	{name='Oshker', code='return max(0., min(r, 1.5));	//replace 1.5 with 1 <= beta <= 2'},
	{name='Ospre', code='return .5 * (r * r + r) / (r * r + r + 1.);'},
	{name='Smart', code='return max(0., min(2. * r, min(.25 + .75 * r, 4.)));'},
	{name='Sweby', code='return max(0., max(min(1.5 * r, 1.), min(r, 1.5)));	//replace 1.5 with 1 <= beta <= 2'},
	{name='UMIST', code='return max(0., min(min(2. * r, .75 + .25 * r), min(.25 + .75 * r, 2.)));'},
	{name='VanAlbada1', code='return (r * r + r) / (r * r + 1.);'},
	{name='VanAlbada2', code='return 2. * r / (r * r + 1.);'},

	--Why isn't this working like it is in the JavaScript code?
	--return (r + fabs(r)) / (1. + fabs(r));
	{name='VanLeer', code='return max(0., r) * 2. / (1. + r);'},
	
	{name='MonotizedCentral', code='return max(0., min(2., min(.5 * (1. + r), 2. * r)));'},
	{name='Superbee', code='return max(0., max(min(1., 2. * r), min(2., r)));'},
	{name='BarthJespersen', code='return .5 * (r + 1.) * min(1., min(4. * r / (r + 1.), 4. / (r + 1.)));'},
}
HydroCLApp.slopeLimiterNames = HydroCLApp.slopeLimiters:map(function(limiter) return limiter.name end)

function HydroCLApp:initGL(...)
	HydroCLApp.super.initGL(self, ...)

	-- TODO favor cl_khr_fp64, cl_khr_3d_image_writes, cl_khr_gl_sharing
--[[
for i,platform in ipairs(require 'cl.platform'.getAll()) do
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
	
	self.platform = require 'cl.platform'.getAll()[1]

	-- TODO favor cl_khr_fp64, cl_khr_3d_image_writes, cl_khr_gl_sharing

	self.device = self.platform:getDevices{gpu=true}[1]
print()
self.device:printInfo()
	local exts = string.split(string.trim(self.device:getExtensions()):lower(),'%s+')
	self.useGLSharing = exts:find(nil, function(ext) return ext:match'cl_%w+_gl_sharing' end)
	
	self.ctx = require 'cl.context'{
		platform = self.platform,
		device = self.device,
		glSharing = self.useGLSharing,
	}
print()
self.ctx:printInfo()

	self.cmds = require 'cl.commandqueue'{context=self.ctx, device=self.device}

	-- making 'real' a parameter of solver would be clever
	-- but because it is using ffi.ctype it might be tough ...
	-- then again, any solver ffi.ctype defined will potentially collide with other solvers ...
	self.real = 'float'
	ffi.cdef('typedef '..self.real..' real;')

	-- create this after 'real' is defined
	--  specifically the call to 'refreshGridSize' within it
	self.solver = Solver{
		app=self, 
		dim=1,
		gridSize={256,256,256}, 
		slopeLimiter='Superbee', 
		
		-- [[
		eqn=require 'euler1d'(),
		mins={-1,-1,-1},
		maxs={1,1,1},
		--]]
	
		--[[
		eqn=require 'adm1d3to5var'(),
		mins={0,0,0},
		maxs={300,300,300},
		--]]
	}
	
	self.solvers = table{self.solver}
	
	local graphShaderCode = file['graph.shader']
	self.graphShader = require 'gl.program'{
		vertexCode = '#define VERTEX_SHADER\n'..graphShaderCode,
		fragmentCode = '#define FRAGMENT_SHADER\n'..graphShaderCode,
		uniforms = {
			'tex', 
			'scale', 
			'xmin', 'xmax', 
			'useLog', 
			'axis', 
			'ambient',
			'size', 
		},
	}
	self.graphShader:use()
	gl.glUniform1i(self.graphShader.uniforms.tex, 0)
	gl.glUniform1f(self.graphShader.uniforms.scale, 1)
	gl.glUniform1i(self.graphShader.uniforms.useLog, false)
	gl.glUniform1f(self.graphShader.uniforms.ambient, 1)	
	self.graphShader:useNone()

	-- [[ need to get image loading working
	local fonttex = require 'gl.tex2d'{
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
	self.font = require 'gui.font'{tex=fonttex}
	--]]
end

local updateMethod

function HydroCLApp:update(...)

	if updateMethod then
		if updateMethod == 'step' then 
			print('performing single step...')
			updateMethod = nil 
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

	for _,varName in ipairs(varNamesEnabled) do
		local xmin, xmax, ymin, ymax
		for _,solver in ipairs(self.solvers) do
			local varIndex = solver.displayVars:find(nil, function(var) return var.name == varName end)
			if varIndex
			--and solver.visiblePtr and solver.visiblePtr[0] 
			then
				local solverymin, solverymax = solver:calcDisplayVarRange(varIndex)

				if not solverymin or not solverymax or solverymin ~= solverymin or solverymax ~= solverymax then
				else	
					local base = 10	-- round to nearest base-10
					local scale = 10 -- ...with increments of 10
					solverymin, solverymax = 1.1 * solverymin - .1 * solverymax, 1.1 * solverymax - .1 * solverymin
					local newymin = (solverymin<0 and -1 or 1)*(math.abs(solverymin)==math.huge and 1e+100 or base^math.log(math.abs(solverymin),base))
					local newymax = (solverymax<0 and -1 or 1)*(math.abs(solverymax)==math.huge and 1e+100 or base^math.log(math.abs(solverymax),base))
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

				local solverxmin, solverxmax = solver.mins[1], solver.maxs[1]
				solverxmin, solverxmax = 1.1 * solverxmin - .1 * solverxmax, 1.1 * solverxmax - .1 * solverxmin
				
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
		for _,solver in ipairs(self.solvers) do
			local varIndex = solver.displayVars:find(nil, function(var) return var.name == varName end)
			if varIndex then
				self:renderDisplayVar(solver, varIndex)
			end

			if self.font then
				local fontSizeX = (xmax - xmin) * .05
				local fontSizeY = (ymax - ymin) * .05
				local ystep = ystep * 2
				for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
					self.font:draw{
						pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
						text=tostring(y),
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
	
		graphCol = graphCol + 1
		if graphCol == graphsWide then
			graphCol = 0
			graphRow = graphRow + 1
		end

	end

	HydroCLApp.super.update(self, ...)
end

function HydroCLApp:renderDisplayVar(solver, varIndex)
	local var = solver.displayVars[varIndex]
	solver:calcDisplayVarToTex(varIndex, var)	
	-- display

	self.graphShader:use()
	solver.tex:bind()

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
	if ig.igButton(updateMethod and 'Stop' or 'Start') then
		updateMethod = not updateMethod
	end
	if ig.igButton'Step' then
		updateMethod = 'step'
	end
	if ig.igButton'Reset' then
		print'resetting...'
		self.solver:resetState()
		updateMethod = nil
	end

	ig.igCheckbox('use fixed dt', self.solver.useFixedDT)
	ig.igInputFloat('fixed dt', self.solver.fixedDT)
	ig.igInputFloat('CFL', self.solver.cfl)

	if ig.igCombo('init state', self.solver.initState, self.solver.eqn.initStates) then
		self.solver:refreshSolverProgram()
	end

	if ig.igCombo('slope limiter', self.solver.slopeLimiter, self.slopeLimiterNames) then
		self:refreshSolverProgram()
	end

	for i=1,self.solver.dim do
		for _,minmax in ipairs(minmaxs) do
			local var = xs[i]..minmax
			if ig.igCombo(var, self.solver.boundaryMethods[var], self.boundaryMethods) then
				self.solver:refreshBoundaryProgram()
			end
		end
	end

	if ig.igCollapsingHeader'variables:' then
		local lastSection
		local sectionEnabled
		for i,var in ipairs(self.solver.displayVars) do
			local section = var.buffer
			if section ~= lastSection then
				sectionEnabled = ig.igCollapsingHeader(section..' variables:')
			end
			if sectionEnabled then
				ig.igCheckbox(var.name, var.enabled)
			end
			lastSection = section
		end
	end
end

function HydroCLApp:event(event, ...)
	HydroCLApp.super.event(self, event, ...)
	if ig.igGetIO()[0].WantCaptureKeyboard then return end
	if event.type == sdl.SDL_KEYDOWN then
		if event.key.keysym.sym == sdl.SDLK_SPACE then
			updateMethod = not updateMethod
		elseif event.key.keysym.sym == ('u'):byte() then
			updateMethod = 'step'
		elseif event.key.keysym.sym == ('r'):byte() then
			print'resetting...'
			self.solver:resetState()
			updateMethod = nil
		end
	end
end

return HydroCLApp
