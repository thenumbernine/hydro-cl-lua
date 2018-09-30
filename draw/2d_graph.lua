local gl = require 'ffi.OpenGL'

return function(HydroCLApp)
			
HydroCLApp.display2D_Graph_step = 1

function HydroCLApp:display2D_Graph(solvers, varName, ar, graph_xmin, graph_ymin, graph_xmax, graph_ymax)
	self.view:projection(ar)
	self.view:modelview()
	gl.glColor3f(1,1,1)
	gl.glEnable(gl.GL_DEPTH_TEST)

	for _,solver in ipairs(solvers) do 
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
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

			local step = self.display2D_Graph_step 
			-- TODO gui this somewhere
			local ambient = .3
			
			-- TODO where to specify using the heatmap gradient vs using the variable/solver color
			gl.glColor3f(table.unpack((#self.solvers > 1 and solver or var).color))

			solver:calcDisplayVarToTex(var)

			local graphShader = solver.graphShader
			graphShader:use()
			solver:getTex(var):bind()

			local scale = 1 / (valueMax - valueMin)
			local offset = valueMin
			gl.glUniform1f(graphShader.uniforms.scale.loc, scale)
			gl.glUniform1f(graphShader.uniforms.offset.loc, offset)
			
			gl.glUniform1f(graphShader.uniforms.ambient.loc, ambient)
			gl.glUniform1i(graphShader.uniforms.axis.loc, solver.dim)
			gl.glUniform1i(graphShader.uniforms.useLog.loc, var.useLog)
			gl.glUniform2f(graphShader.uniforms.size.loc, solver.gridSize.x, solver.gridSize.y)
			gl.glUniform2f(graphShader.uniforms.xmin.loc, solver.mins[1], solver.mins[2])
			gl.glUniform2f(graphShader.uniforms.xmax.loc, solver.maxs[1], solver.maxs[2])

			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

			for ybase=2,tonumber(solver.gridSize.y)-2-step,step do
				gl.glBegin(gl.GL_TRIANGLE_STRIP)
				for x=2,tonumber(solver.gridSize.x)-2,step do
					for yofs=0,step,step do
						local y = ybase + yofs
						gl.glVertex2d(
							(x + .5) / tonumber(solver.gridSize.x),
							(y + .5) / tonumber(solver.gridSize.y))
					end
				end
				gl.glEnd()
			end
			
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
			
			graphShader:useNone()
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
end

end
