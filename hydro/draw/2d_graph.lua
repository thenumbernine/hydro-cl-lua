-- TODO make use of app.display_useCoordMap
local ffi = require 'ffi'
local path = require 'ext.path'
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local GLSceneObject = require 'gl.sceneobject'
local vector = require 'ffi.cpp.vector-lua'
local Draw = require 'hydro.draw.draw'


local Draw2DGraph = Draw:subclass()

Draw2DGraph.step = 1

-- TODO gui this somewhere
Draw2DGraph.ambient = .3

function Draw2DGraph:display(varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
	local solver = self.solver
	local app = solver.app

	app.view:setup(ar)

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

function Draw2DGraph:prepareShader()
	local solver = self.solver
	local app = solver.app

	self:prepareGraphShader()

	solver.draw2DGraphSceneObj = solver.draw2DGraphSceneObj or GLSceneObject{
		program = solver.graphShader,
		vertexes = {
			useVec = true,
			dim = 3,
			usage = gl.GL_DYNAMIC_DRAW,
		},
		geometry = {
			mode = gl.GL_LINE_STRIP,
			count = 0,
		},
	}
end

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

	local sceneObj = solver.draw2DGraphSceneObj

	local step = math.max(1, self.step)
	local numX = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)
	local numY = math.floor((tonumber(solver.gridSize.y) - 2 * solver.numGhost + 1) / step)
	-- 2 vtxs per tristrip
	local numVertexes = 2 * numX * numY

-- [[ new sys ... half the speed of above
	local vertexes = sceneObj.attrs.vertex.buffer:beginUpdate()
	vertexes:resize(numVertexes)
	local index = 0
	for j=0,numY-2 do
		for i=0,numX-1 do
			local x = i * step
			for jofs=0,1 do
				local y = (j + jofs) * step
				--[=[ goes horribly slow
				vertexes:emplace_back():set(x,y,0)
				--]=]
				-- [=[ just as fast as old sys
				local v = vertexes.v[index]
				v.x = x
				v.y = y
				v.z = 0--app.displayFixedZ
				index = index + 1
				--]=]
			end
		end
	end
	sceneObj.attrs.vertex.buffer:endUpdate()
--]]

	local shader = solver.graphShader
	local uniforms = shader.uniforms

	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

	shader:use()
	local tex = solver:getTex(var)
	tex:bind()

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.ambient.loc, self.ambient)

	-- TODO where to specify using the heatmap gradient vs using the variable/solver color
	gl.glUniform3f(uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	solver.draw2DGraphSceneObj.attrs.vertex.buffer
		:bind()
		:updateData()

	solver.draw2DGraphSceneObj:enableAndSetAttrs()

	-- TODO use ELEMENT_ARRAY_BUFFER?
	for j=0,numY-2 do
		gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, j * 2 * numX, 2 * numX)
	end

	solver.draw2DGraphSceneObj:disableAttrs()

	tex:unbind()
	shader:useNone()

	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
end

return Draw2DGraph
