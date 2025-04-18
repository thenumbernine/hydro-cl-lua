--[[
parent class for some common functions that all draw glsl code uses
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local path = require 'ext.path'
local gl = require 'gl'
local matrix_ffi = require 'matrix.ffi'


local Draw = class()

Draw.glslVersion = cmdline.glslVersion
	and '#version '..cmdline.glslVersion
	or require 'gl.program'.getVersionPragma()

function Draw:init(solver)
	self.solver = assert(solver)
end

function Draw:getCommonGLSLFragCode()
	local solver = self.solver
	return solver.eqn:template(path'hydro/draw/draw.glsl':read(), {
		draw = self,
		clnumber = require 'cl.obj.number',
	})
end

function Draw:setupDisplayVarShader(shader, var, valueMin, valueMax)
	local solver = self.solver
	local app = solver.app
	local isMeshSolver = require 'hydro.solver.meshsolver':isa(solver)

	local uniforms = shader.uniforms
	if uniforms.displayDim then
		gl.glUniform1i(uniforms.displayDim.loc, app.displayDim)
	end
	if uniforms.displayFixed then
		gl.glUniform2f(uniforms.displayFixed.loc, app.displayFixedY, app.displayFixedZ)
	end
	if uniforms.mvProjMat then
		gl.glUniformMatrix4fv(uniforms.mvProjMat.loc, 1, gl.GL_FALSE, app.view.mvProjMat.ptr)
	end
	if uniforms.normalMatrix then
		self.normalMatrix = self.normalMatrix or matrix_ffi.zeros({3,3}, 'float')
		--gl.glGetFloatv(gl.GL_NORMAL_MATRIX, self.normalMatrix.ptr)
		-- invert app.view.mvMat's upper 3x3 into normalMatrix, then transpose
		-- but if it is purely a rotation matrix, this is the same as just the 3x3 portion ...
		-- so TODO why even use this at all? just get the submatrix in GLSL ...
		for j=0,2 do
			for i=0,2 do
				self.normalMatrix.ptr[i + 3 * j] = app.view.mvProjMat.ptr[i + 4 * j]
			end
		end
		gl.glUniformMatrix3fv(uniforms.normalMatrix.loc, 1, gl.GL_FALSE, self.normalMatrix.ptr)
	end
	if uniforms.useCoordMap then
		gl.glUniform1i(uniforms.useCoordMap.loc, app.display_useCoordMap and 1 or 0)
	end
	if uniforms.useLog then
		gl.glUniform1i(uniforms.useLog.loc, var.useLog and 1 or 0)
	end
	if uniforms.valueMin then
		gl.glUniform1f(uniforms.valueMin.loc, valueMin)
	end
	if uniforms.valueMax then
		gl.glUniform1f(uniforms.valueMax.loc, valueMax)
	end
	if uniforms.solverMins then
		gl.glUniform3f(uniforms.solverMins.loc, solver.mins:unpack())
	end
	if uniforms.solverMaxs then
		gl.glUniform3f(uniforms.solverMaxs.loc, solver.maxs:unpack())
	end
	if uniforms.cartesianMin then
		gl.glUniform3f(uniforms.cartesianMin.loc, solver.cartesianMin:unpack())
	end
	if uniforms.cartesianMax then
		gl.glUniform3f(uniforms.cartesianMax.loc, solver.cartesianMax:unpack())
	end

	if uniforms.texSize then
		gl.glUniform3f(uniforms.texSize.loc, solver.texSize:unpack())
	end
	if uniforms.gridSize then
		if not isMeshSolver then
			gl.glUniform3f(uniforms.gridSize.loc, solver.gridSize:unpack())
		else
			gl.glUniform3f(uniforms.gridSize.loc, solver.texSize:unpack())
		end
	end
	if uniforms.sizeWithoutBorder then
		if not isMeshSolver then
			gl.glUniform3f(uniforms.sizeWithoutBorder.loc, solver.sizeWithoutBorder:unpack())
		else
			gl.glUniform3f(uniforms.sizeWithoutBorder.loc, solver.texSize:unpack())
		end
	end
	if uniforms.displaySliceAngle then
		--gl.glUniform4fv(uniforms.displaySliceAngle.loc, 4, app.displaySliceAngle.s)
		gl.glUniform4f(uniforms.displaySliceAngle.loc, app.displaySliceAngle:unpack())
	end
	if app.useClipPlanes then
		-- TODO use an array?
		-- use a uniform struct?
		for i,info in ipairs(app.clipInfos) do
			local k = 'clipEnabled'..i
			local u = uniforms[k]
			if u then
				gl.glUniform1i(u.loc, info.enabled and 1 or 0)
			end
			local k = 'clipPlane'..i
			local u = uniforms[k]
			if u then
				gl.glUniform4f(u.loc, info.plane:unpack())
			end
		end
	end
end

local function makeGLSL(code)
	return table{
		'#define M_PI '..('%.50f'):format(math.pi),
		'#define inline',
		'#define real					float',
		'#define _real3					vec3',
		'#define real3					vec3',
		'#define real3_zero				vec3(0.,0.,0.)',
		'#define real3_add(a,b)			((a)+(b))',
		'#define real3_sub(a,b)			((a)-(b))',
		'#define real3_real_mul(a,b)	((a)*(b))',
		'#define real3_dot				dot',
		'#define real3_len(a)			length(a)',
		'#define real3_lenSq(a)			dot(a,a)',
		'#define atan2					atan',
		'#define fmod					mod',
		'#define fabs					abs',
		'',
		(code:gsub('static inline ', '')),
	}:concat'\n'
end

function Draw:getModuleCodeGLSL(...)
	return makeGLSL(self.solver.modules:getCodeAndHeader(...))
end

-- used in draw/1d_graph.lua and draw/2d_graph.lua
function Draw:prepareGraphShader()
	local solver = self.solver

	if solver.graphShader then return end

	local graphShaderCode = assert(path'hydro/draw/graph.glsl':read())

	solver.graphShader = solver.GLProgram{
		name = 'graph',
		vertexCode = solver.eqn:template(graphShaderCode, {
			draw = self,
			vertexShader = true,
		}),
		fragmentCode = solver.eqn:template(graphShaderCode, {
			draw = self,
			fragmentShader = true,
		}),
		uniforms = {
			tex = 0,
			scale = 1,
			ambient = 1,
		},
	}:useNone()
end

-- common function
function Draw:drawGrid(xmin, xmax, ymin, ymax)
	local app = self.solver.app

	local gridz = 0	--.1

	local sceneObj = app.drawLineSceneObj
	local shader = sceneObj.program

	shader:use()
	sceneObj:enableAndSetAttrs()
	gl.glUniformMatrix4fv(shader.uniforms.mvProjMat.loc, 1, gl.GL_FALSE, app.view.mvProjMat.ptr)
	gl.glUniform4f(shader.uniforms.color.loc, .1, .1, .1, 1)

	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	for x=xticmin,xticmax do
		-- TODO turn this into instanced geometry and make it draw faster
		gl.glUniform3f(shader.uniforms.pt0.loc, x*xstep, ymin, gridz)
		gl.glUniform3f(shader.uniforms.pt1.loc, x*xstep, ymax, gridz)
		sceneObj.geometry:draw()
	end
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	for y=yticmin,yticmax do
		-- TODO turn this into instanced geometry and make it draw faster
		gl.glUniform3f(shader.uniforms.pt0.loc, xmin, y*ystep, gridz)
		gl.glUniform3f(shader.uniforms.pt1.loc, xmax, y*ystep, gridz)
		sceneObj.geometry:draw()
	end

	gl.glUniform4f(shader.uniforms.color.loc, .5, .5, .5, 1)

	gl.glUniform3f(shader.uniforms.pt0.loc, xmin, 0, gridz)
	gl.glUniform3f(shader.uniforms.pt1.loc, xmax, 0, gridz)
	sceneObj.geometry:draw()

	gl.glUniform3f(shader.uniforms.pt0.loc, 0, ymin, gridz)
	gl.glUniform3f(shader.uniforms.pt1.loc, 0, ymax, gridz)
	sceneObj.geometry:draw()

	sceneObj:disableAttrs()
	shader:useNone()

	return ystep
end

return Draw
