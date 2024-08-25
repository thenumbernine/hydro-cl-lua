local path = require 'ext.path'
local asserteq = require 'ext.assert'.eq
local gl = require 'gl'
local GLSceneObject = require 'gl.sceneobject'
local Draw = require 'hydro.draw.draw'


local Draw2DHeatmap = Draw:subclass()

function Draw2DHeatmap:display(varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
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

	-- TODO one grid for all displaly instead of multiple calls...
	self:drawGrid(xmin, xmax, ymin, ymax)

	-- NOTICE overlays of multiple solvers won't be helpful.  It'll just draw over the last solver.
	-- I've got to rethink the visualization
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax)
	end
--	gl.glDisable(gl.GL_DEPTH_TEST)
end

function Draw2DHeatmap:prepareShader()
	local solver = self.solver
	if solver.heatMap2DSceneObj then return end

	local heatMapCode = assert(path'hydro/draw/2d_heatmap.glsl':read())

	solver.heatMap2DShader = solver.GLProgram{
		name = '2d_heatmap',
		vertexCode = solver.eqn:template(heatMapCode, {
			draw = self,
			vertexShader = true,
		}),
		fragmentCode = solver.eqn:template(heatMapCode, {
			draw = self,
			fragmentShader = true,
		}),
		uniforms = {
			tex = 0,
			gradientTex = 1,
		},
	}:useNone()

	solver.heatMap2DSceneObj = GLSceneObject{
		program = solver.heatMap2DShader,
		geometry = solver.app.quadGeom,
	}
end

function Draw2DHeatmap:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax)
	local solver = self.solver
	local app = solver.app
	if require 'hydro.solver.meshsolver':isa(solver) then return end

	-- TODO allow a fixed, manual colormap range
	-- NOTICE with AMR this will only get from the root node
	--  which should at least have blitters of the children
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	local shader = solver.heatMap2DSceneObj.program
	shader:use()
	app.gradientTex:bind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)
	self:drawSolverWithVar(var, shader, xmin, xmax, ymin, ymax)

-- [[
	if solver.amr then
		for k,subsolver in pairs(solver.amr.child) do
			self:drawSolverWithVar(app, subsolver, var, shader, xmin, xmax, ymin, ymax)
		end
	end
--]]

	app.gradientTex:unbind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)
	shader:useNone()

	gl.glDisable(gl.GL_BLEND)
--	gl.glDisable(gl.GL_DEPTH_TEST)

	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)

--	gl.glEnable(gl.GL_DEPTH_TEST)
end

function Draw2DHeatmap:drawSolverWithVar(var, shader, xmin, xmax, ymin, ymax)
	local solver = self.solver
	local app = solver.app
-- hmm ... this is needed for sub-solvers
local origSolver = var.solver
var.solver = solver

	solver:calcDisplayVarToTex(var)

	local tex = solver:getTex(var)
		:bind(0)
		:setParameter(gl.GL_TEXTURE_MAG_FILTER, app.displayBilinearTextures and gl.GL_LINEAR or gl.GL_NEAREST)

	gl.glUniform4f(shader.uniforms.bbox.loc, xmin, ymin, xmax, ymax)

	-- TODO no more need to pass along shader every time?  or nah?
	asserteq(solver.heatMap2DSceneObj.program, shader)
	solver.heatMap2DSceneObj:draw()

	tex:unbind(0)

var.solver = origSolver
end

return Draw2DHeatmap
