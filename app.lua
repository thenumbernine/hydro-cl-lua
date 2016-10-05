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
	self.solver = Solver{app=self, gridSize=256}

	local GLProgram = require 'gl.program'
	local graphShaderCode = file['graph.shader']
	self.graphShader = GLProgram{
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
end

local xScale = ffi.new('float[1]', 1)
local yScale = ffi.new('float[1]', .5)

local updateMethod

function HydroCLApp:update(...)

	if updateMethod then
		if updateMethod == 'step' then 
			print('performing single step...')
			updateMethod = nil 
		end

		self.solver:update()
	end

	gl.glClearColor(.3,.2,.5,1)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	
	local ar = self.width / self.height
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(-ar, ar, -1, 1, -1, 1)

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glScalef(xScale[0], yScale[0], 1)

	for i,var in ipairs(self.solver.displayVars) do
		if var.enabled[0] then
			self:renderDisplayVar(i, var)
		end
	end

	HydroCLApp.super.update(self, ...)
end

function HydroCLApp:renderDisplayVar(i, var)
	self.solver:convertToTex(i, var)	
	-- display

	self.graphShader:use()
	self.solver.tex:bind()

	gl.glUniform2f(self.graphShader.uniforms.xmin, self.solver.xmin, 0)
	gl.glUniform2f(self.graphShader.uniforms.xmax, self.solver.xmax, 0)
	gl.glUniform1i(self.graphShader.uniforms.axis, self.solver.dim)
	gl.glUniform2f(self.graphShader.uniforms.size, self.solver.gridSize.x, self.solver.gridSize.y)

	gl.glColor3f(table.unpack(var.color))
	gl.glBegin(gl.GL_LINE_STRIP)
	local step = 1
	for i=2,tonumber(self.solver.gridSize.x)-2,step do
		local x = (i+.5)/tonumber(self.solver.gridSize.x)
		gl.glVertex2f(x, 0)
	end
	gl.glEnd()
	
	self.solver.tex:unbind()
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

	ig.igSliderFloat('x scale', xScale, 0, 100, '%.3f', 10)
	ig.igSliderFloat('y scale', yScale, 0, 100, '%.3f', 10)

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
