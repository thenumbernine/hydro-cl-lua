-- TODO make use of app.display_useCoordMap
local class = require 'ext.class'
local path = require 'ext.path'
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local vector = require 'ffi.cpp.vector'
local Draw = require 'hydro.draw.draw'


local Draw2DGraph = class(Draw)

Draw2DGraph.step = 1

-- TODO gui this somewhere
Draw2DGraph.ambient = .3

function Draw2DGraph:showDisplayVar(var)
	local solver = self.solver
	local app = solver.app

	-- TODO allow a fixed, manual colormap range
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end

	solver:calcDisplayVarToTex(var)


	if app.displayDim == 3 then
		io.stderr:write'Why are you using a graph shader to display 3D data?  Use a 3D display instead.\n'
		return
	end


	local shader = solver.graphShader
	local uniforms = shader.uniforms

	shader:use()
	local tex = solver:getTex(var)
	tex:bind()

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.ambient.loc, self.ambient)

	-- TODO where to specify using the heatmap gradient vs using the variable/solver color
	gl.glUniform3f(uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	-- 3 components per vertex
	if not self.vertexes then self.vertexes = vector'vec3f_t' end

	local step = math.max(1, self.step)
	local numX = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)
	local numY = math.floor((tonumber(solver.gridSize.y) - 2 * solver.numGhost + 1) / step)
	-- 2 vtxs per tristrip
	local numVertexes = 2 * numX * numY

	if #self.vertexes ~= numVertexes then
		self.vertexes:resize(numVertexes)
	end

	local index = 0
	for j=0,numY-2 do
		for i=0,numX-1 do
			local x = i * step
			for jofs=0,1 do
				local y = (j + jofs) * step
				local v = self.vertexes.v[index]
				v.x = x
				v.y = y
				v.z = 0--app.displayFixedZ
				index = index + 1
			end
		end
	end

	gl.glEnableVertexAttribArray(shader.attrs.gridCoord.loc)
	gl.glVertexAttribPointer(shader.attrs.gridCoord.loc, 3, gl.GL_FLOAT, false, 0, self.vertexes.v)

	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

	for j=0,numY-2 do
		gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, j * 2 * numX, 2 * numX)
	end

	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)

	gl.glDisableVertexAttribArray(shader.attrs.gridCoord.loc)

	tex:unbind()
	shader:useNone()
end

function Draw2DGraph:display(varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
	local solver = self.solver
	local app = solver.app

	app.view:setup(ar)

	gl.glColor3f(1,1,1)
	gl.glEnable(gl.GL_DEPTH_TEST)

	if not require 'hydro.solver.meshsolver':isa(solver) then
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:prepareShader()
			self:showDisplayVar(var)
		end

		-- TODO right here is where the color gradient display usually goes
		-- mind you I'm not using it in the 2D graph display atm
	end

	gl.glDisable(gl.GL_DEPTH_TEST)
end

-- also in 1d.lua.  subclass?
function Draw2DGraph:prepareShader()
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

return Draw2DGraph
