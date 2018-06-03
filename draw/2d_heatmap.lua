local gl = require 'ffi.OpenGL'

return function(HydroCLApp)

function HydroCLApp:display2D_Heatmap(solvers, varName, ar, graph_xmin, graph_ymin, graph_xmax, graph_ymax)
	self.view:projection(ar)
	self.view:modelview()
	if self.view.getOrthoBounds then
		xmin, xmax, ymin, ymax = self.view:getOrthoBounds(ar)
	else
		xmin, xmax, ymin, ymax = graph_xmin, graph_ymin, graph_xmax, graph_ymax
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

			solver:calcDisplayVarToTex(var)
	
			solver.heatMap2DShader:use()
			gl.glUniform1i(solver.heatMap2DShader.uniforms.useLog.loc, var.useLog)
			gl.glUniform1f(solver.heatMap2DShader.uniforms.valueMin.loc, valueMin)
			gl.glUniform1f(solver.heatMap2DShader.uniforms.valueMax.loc, valueMax)
			solver:getTex(var):bind(0)
			self.gradientTex:bind(1)
	
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
			gl.glEnable(gl.GL_BLEND)

			gl.glBegin(gl.GL_QUADS)
			gl.glVertex2d(xmin, ymin)
			gl.glVertex2d(xmax, ymin)
			gl.glVertex2d(xmax, ymax)
			gl.glVertex2d(xmin, ymax)
			gl.glEnd()
	
			gl.glDisable(gl.GL_BLEND)
			
			self.gradientTex:unbind(1)
			solver:getTex(var):unbind(0)
			solver.heatMap2DShader:useNone()
			
--			gl.glDisable(gl.GL_DEPTH_TEST)

			self:drawGradientLegend(ar, varName, valueMin, valueMax)
		
--			gl.glEnable(gl.GL_DEPTH_TEST)
		end
	end
	
--	gl.glDisable(gl.GL_DEPTH_TEST)
end

end
