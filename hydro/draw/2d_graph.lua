-- TODO make use of app.display_useCoordMap
local ffi = require 'ffi'
local template = require 'template'
local class = require 'ext.class'
local file = require 'ext.file'
local gl = require 'ffi.OpenGL'
local matrix_ffi = require 'matrix.ffi'
local Draw = require 'hydro.draw.draw'


local Draw2DGraph = class(Draw)

Draw2DGraph.step = 1

-- TODO gui this somewhere
Draw2DGraph.ambient = .3

function Draw2DGraph:showDisplayVar(app, solver, var)
	
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
		do return end
	end


	local shader = solver.graphShader
	local uniforms = shader.uniforms

	shader:use()
	local tex = solver:getTex(var)
	tex:bind()
	
	self:setupDisplayVarShader(shader, app, solver, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.ambient.loc, self.ambient)
	
	-- TODO where to specify using the heatmap gradient vs using the variable/solver color
	gl.glUniform3f(uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	local step = math.max(1, self.step)
	local numX = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)
	local numY = math.floor((tonumber(solver.gridSize.y) - 2 * solver.numGhost + 1) / step)
	local numVertexes2D = numX * numY
	
	if numVertexes2D ~= self.numVertexes2D then
		self.numVertexes2D = numVertexes2D
		-- 2 vtxs per tristrip * 3 components per vertex
		self.vertexes = ffi.new('float[?]', 2*3*numVertexes2D)
	end

	for j=0,numY-2 do
		for i=0,numX-1 do
			local x = (i * step + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
			for jofs=0,1 do
				local y = ((j + jofs) * step + .5 + solver.numGhost) / tonumber(solver.gridSize.y)
				self.vertexes[0 + 3 * (jofs + 2 * (i + numX * j))] = x
				self.vertexes[1 + 3 * (jofs + 2 * (i + numX * j))] = y
				self.vertexes[2 + 3 * (jofs + 2 * (i + numX * j))] = 0--app.displayFixedZ
			end
		end
	end
	
	gl.glEnableVertexAttribArray(shader.attrs.inVertex.loc)
	gl.glVertexAttribPointer(shader.attrs.inVertex.loc, 3, gl.GL_FLOAT, false, 0, self.vertexes)
	
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

	for j=0,numY-2 do
		gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, j * 2 * numX, 2 * numX) 
	end
	
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
	
	gl.glDisableVertexAttribArray(shader.attrs.inVertex.loc)
	
	tex:unbind()
	shader:useNone()
end

function Draw2DGraph:display(app, solvers, varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
	app.view:setup(ar)
	
	gl.glColor3f(1,1,1)
	gl.glEnable(gl.GL_DEPTH_TEST)

	for _,solver in ipairs(solvers) do 
		if not require 'hydro.solver.meshsolver'.is(solver) then
			local var = solver.displayVarForName[varName]
			if var and var.enabled then
				self:prepareShader(solver)
				self:showDisplayVar(app, solver, var)
			end
			
			-- TODO right here is where the color gradient display usually goes
			-- mind you I'm not using it in the 2D graph display atm
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
end

-- also in 1d.lua.  subclass?
function Draw2DGraph:prepareShader(solver)
	if solver.graphShader then return end
	
	local graphShaderCode = assert(file['hydro/draw/graph.shader'])
	
	solver.graphShader = solver.GLProgram{
		name = 'graph',
		vertexCode = template(graphShaderCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			vertexShader = true,
		}),
		fragmentCode = template(graphShaderCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			fragmentShader = true,
		}),
		uniforms = {
			tex = 0,
			scale = 1,
			ambient = 1,
		},
	}
end

return Draw2DGraph
