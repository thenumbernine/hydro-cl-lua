local ffi = require 'ffi'
local class = require 'ext.class'
local path = require 'ext.path'
local gl = require 'ffi.OpenGL'
local GLTex2D = require 'gl.tex2d'
local Draw = require 'hydro.draw.draw'


local DrawVectorLIC = class(Draw)

DrawVectorLIC.integralMaxIter = 10

--[[
the xmin/xmax/ymin/ymax passed as arguments is the one shared with all graphs
but for LIC, we just need the cartesian range
--]]
function DrawVectorLIC:drawSolverWithVar(var, shader, xmin, xmax, ymin, ymax)
	local solver = self.solver
	local app = solver.app

-- hmm ... this is needed for sub-solvers
local origSolver = var.solver
var.solver = solver

	solver:calcDisplayVarToTex(var)

	local tex = solver:getTex(var)
	tex:bind(0)
	tex:setParameter(gl.GL_TEXTURE_MAG_FILTER, app.displayBilinearTextures and gl.GL_LINEAR or gl.GL_NEAREST)
	self.noiseTex:bind(2)

	gl.glBegin(gl.GL_QUADS)
	gl.glVertex2d(xmin, ymin)
	gl.glVertex2d(xmax, ymin)
	gl.glVertex2d(xmax, ymax)
	gl.glVertex2d(xmin, ymax)
	gl.glEnd()

	self.noiseTex:unbind(2)
	tex:unbind(0)

var.solver = origSolver
end

function DrawVectorLIC:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax)
	local solver = self.solver
	local app = solver.app

	-- TODO allow a fixed, manual colormap range
	-- NOTICE with AMR this will only get from the root node
	--  which should at least have blitters of the children
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		local component = solver.displayComponentFlatList[var.component]
		local vectorField = solver:isVarTypeAVectorField(component.type)
		if vectorField then
			-- calc range of magnitude of vector variable
			valueMin, valueMax = solver:calcDisplayVarRange(var, component.magn)
		else
			valueMin, valueMax = solver:calcDisplayVarRange(var)
		end
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end

	if not self.noiseTex then
		local noiseVol = app.drawVectorLICNoiseSize * app.drawVectorLICNoiseSize * 4
		local noiseData = ffi.new('float[?]', noiseVol)
		for i=0,noiseVol-1 do
			noiseData[i] = math.random()
		end
		self.noiseTex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,
			width = app.drawVectorLICNoiseSize,
			height = app.drawVectorLICNoiseSize,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			data = noiseData,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_LINEAR,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}
	end

	local shader = solver.vectorLICShader
	shader:use()

	app.gradientTex:bind(1)

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	if shader.uniforms.integralMaxIter then
		gl.glUniform1i(shader.uniforms.integralMaxIter.loc, self.integralMaxIter)
	end

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	self:drawSolverWithVar(var, shader, xmin, xmax, ymin, ymax)

-- [[
	if solver.amr then
		local tmp = self.solver
		for k,subsolver in pairs(solver.amr.child) do
			self.solver = subsolver
			self:drawSolverWithVar(var, shader, xmin, xmax, ymin, ymax)
		end
		self.solver = tmp
	end
--]]

	gl.glDisable(gl.GL_BLEND)

	app.gradientTex:unbind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)
	shader:useNone()

--	gl.glDisable(gl.GL_DEPTH_TEST)

	-- TODO only draw the first
	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)

--	gl.glEnable(gl.GL_DEPTH_TEST)
end

function DrawVectorLIC:display(varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
	local solver = self.solver
	local app = solver.app

	app.view:setup(ar)

	local xmin, xmax, ymin, ymax
	if app.view.getOrthoBounds then
		xmin, xmax, ymin, ymax = app.view:getOrthoBounds(ar)
	else
		xmin = solver.cartesianMin.x
		ymin = solver.cartesianMin.y
		xmax = solver.cartesianMax.x
		ymax = solver.cartesianMax.y
	end

--	gl.glEnable(gl.GL_DEPTH_TEST)

	local gridz = 0	--.1

	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex3f(x*xstep,ymin, gridz)
		gl.glVertex3f(x*xstep,ymax, gridz)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex3f(xmin,y*ystep, gridz)
		gl.glVertex3f(xmax,y*ystep, gridz)
	end
	gl.glEnd()

	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex3f(xmin, 0, gridz)
	gl.glVertex3f(xmax, 0, gridz)
	gl.glVertex3f(0, ymin, gridz)
	gl.glVertex3f(0, ymax, gridz)
	gl.glEnd()

	-- NOTICE overlays of multiple solvers won't be helpful.  It'll just draw over the last solver.
	-- I've got to rethink the visualization
	if not require 'hydro.solver.meshsolver':isa(solver) then
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:prepareShader()
			self:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax)
		end
	end
--	gl.glDisable(gl.GL_DEPTH_TEST)
end

function DrawVectorLIC:prepareShader()
	local solver = self.solver

	if solver.vectorLICShader then return end

	local vectorLICCode = assert(path'hydro/draw/vector_lic.glsl':read())

	solver.vectorLICShader = solver.GLProgram{
		name = 'vector_lic',
		vertexCode = solver.eqn:template(vectorLICCode, {
			draw = self,
			vertexShader = true,
		}),
		fragmentCode = solver.eqn:template(vectorLICCode, {
			draw = self,
			fragmentShader = true,
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

return DrawVectorLIC
