-- TODO make use of app.display_useCoordMap
local gl = require 'ffi.OpenGL'
local ffi = require 'ffi'
local class = require 'ext.class'
local matrix_ffi = require 'matrix.ffi'


local Draw2DGraph = class()

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

	local graphShader = solver.graphShader
	
	graphShader:use()
	solver:getTex(var):bind()

	local scale = 1 / (valueMax - valueMin)
	local offset = valueMin
	gl.glUniform1f(graphShader.uniforms.scale.loc, scale)
	gl.glUniform1f(graphShader.uniforms.offset.loc, offset)
	gl.glUniform1f(graphShader.uniforms.ambient.loc, self.ambient)
	gl.glUniform1i(graphShader.uniforms.useLog.loc, var.useLog)
	gl.glUniform2f(graphShader.uniforms.xmin.loc, solver.mins.x, solver.mins.y)
	gl.glUniform2f(graphShader.uniforms.xmax.loc, solver.maxs.x, solver.maxs.y)

	local displayDim = app.displayDim -- solver.dim
	if displayDim == 3 then
		error'Why are you using a graph shader to display 3D data?  Use a 3D display instead.'
	else
		gl.glUniform1i(graphShader.uniforms.axis.loc, displayDim)
	end

	gl.glUniform2f(graphShader.uniforms.size.loc, solver.gridSize.x, solver.gridSize.y)
	
	-- TODO where to specify using the heatmap gradient vs using the variable/solver color
	gl.glUniform3f(graphShader.uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	local step = math.max(1, self.step)
	local numX = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)
	local numY = math.floor((tonumber(solver.gridSize.y) - 2 * solver.numGhost + 1) / step)
	local numVertexes2D = numX * numY
	
	if numVertexes2D ~= self.numVertexes2D then
		self.numVertexes2D = numVertexes2D
		-- 2 vtxs per tristrip * 3 components per vertex
		self.vertexes = ffi.new('float[?]', 2*3*numVertexes2D)
	end

	self.ModelViewMatrix = self.ModelViewMatrix or matrix_ffi(nil, 'float', {4,4})--ffi.new'float[16]'
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.ModelViewMatrix.ptr)

	self.ProjectionMatrix = self.ProjectionMatrix or matrix_ffi(nil, 'float', {4,4})--ffi.new'float[16]'
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.ProjectionMatrix.ptr)

	self.ModelViewProjectionMatrix = self.ModelViewProjectionMatrix or require 'matrix.ffi'(nil, 'float', {4,4})--ffi.new'float[16]'
	matrix_ffi.mul(self.ModelViewProjectionMatrix, self.ProjectionMatrix, self.ModelViewMatrix)

	gl.glUniformMatrix4fv(graphShader.uniforms.ModelViewProjectionMatrix.loc, 1, gl.GL_FALSE, self.ModelViewProjectionMatrix.ptr)

	for j=0,numY-2 do
		for i=0,numX-1 do
			local x = (i * step + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
			for jofs=0,1 do
				local y = ((j + jofs) * step + .5 + solver.numGhost) / tonumber(solver.gridSize.y)
				self.vertexes[0 + 3 * (jofs + 2 * (i + numX * j))] = x
				self.vertexes[1 + 3 * (jofs + 2 * (i + numX * j))] = y
				self.vertexes[2 + 3 * (jofs + 2 * (i + numX * j))] = app.displayFixedZ
			end
		end
	end
	
	gl.glEnableVertexAttribArray(graphShader.attrs.inVertex.loc)
	gl.glVertexAttribPointer(graphShader.attrs.inVertex.loc, 3, gl.GL_FLOAT, false, 0, self.vertexes)
	
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

	for j=0,numY-2 do
		gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, j * 2 * numX, 2 * numX) 
	end
	
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
	
	gl.glDisableVertexAttribArray(graphShader.attrs.inVertex.loc)
	
	solver:getTex(var):unbind()
	graphShader:useNone()
end

function Draw2DGraph:display(app, solvers, varName, ar, graph_xmin, graph_xmax, graph_ymin, graph_ymax)
	app.view:projection(ar)
	app.view:modelview()
	gl.glColor3f(1,1,1)
	gl.glEnable(gl.GL_DEPTH_TEST)

	for _,solver in ipairs(solvers) do 
		if not require 'hydro.solver.meshsolver'.is(solver) then
			local var = solver.displayVarForName[varName]
			if var and var.enabled then
				self:showDisplayVar(app, solver, var)
			end
			
			-- TODO right here is where the color gradient display usually goes
			-- mind you I'm not using it in the 2D graph display atm
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:display2D_Graph(...)
		if not self.draw2DGraph then self.draw2DGraph = Draw2DGraph() end
		self.draw2DGraph:display(self, ...)
	end
end
