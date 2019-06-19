local gl = require 'ffi.OpenGL'

return function(HydroCLApp)

function HydroCLApp:showDisplayVar1D(solver, var)
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
	-- display

	self.graphShader:use()
	solver:getTex(var):bind()

	gl.glUniform1f(self.graphShader.uniforms.scale.loc, 1)
	gl.glUniform1f(self.graphShader.uniforms.offset.loc, 0)
	gl.glUniform1f(self.graphShader.uniforms.ambient.loc, 1)
	gl.glUniform1i(self.graphShader.uniforms.useLog.loc, var.useLog)
	gl.glUniform2f(self.graphShader.uniforms.xmin.loc, solver.mins[1], 0)
	gl.glUniform2f(self.graphShader.uniforms.xmax.loc, solver.maxs[1], 0)
	gl.glUniform1i(self.graphShader.uniforms.axis.loc, solver.dim)
	gl.glUniform2f(self.graphShader.uniforms.size.loc, solver.gridSize.x, solver.gridSize.y)

	gl.glColor3d(table.unpack((#self.solvers > 1 and solver or var).color))
	gl.glBegin(gl.GL_LINE_STRIP)
	local step = 1
	for i=2,tonumber(solver.gridSize.x)-2,step do
		local x = (i+.5)/tonumber(solver.gridSize.x)
		gl.glVertex2d(x, 0)
	end
	gl.glEnd()
	
	solver:getTex(var):unbind()
	self.graphShader:useNone()
end

function HydroCLApp:display1D(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog, valueMin, valueMax)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	gl.glColor3d(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex2d(x*xstep,ymin)
		gl.glVertex2d(x*xstep,ymax)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex2d(xmin,y*ystep)
		gl.glVertex2d(xmax,y*ystep)
	end
	gl.glEnd()

	gl.glColor3d(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex2d(xmin, 0)
	gl.glVertex2d(xmax, 0)
	gl.glVertex2d(0, ymin)
	gl.glVertex2d(0, ymax)
	gl.glEnd()

	-- display here
	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:showDisplayVar1D(solver, var)
		end

		local unitScale = 1
		local thisValueMin = valueMin
		local thisValueMax = valueMax
		if var.showInUnits and var.units then
			unitScale = solver:convertToSIUnitsCode(var.units).func()
			thisValueMin = thisValueMin * unitScale
			thisValueMax = thisValueMax * unitScale
		end			
	
		if self.font then
			local fontSizeX = (xmax - xmin) * .05
			local fontSizeY = (ymax - ymin) * .05
			local ystep = ystep * 2
			for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
				local realY = useLog and 10^y or y
				local value = (y - ymin) * (thisValueMax - thisValueMin) / (ymax - ymin) + thisValueMin
				local absvalue = math.abs(value)
				self.font:draw{
					pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
					text=(
						(absvalue > 1e+5 or absvalue < 1e-5)
						and ('%.5e'):format(value) or ('%.5f'):format(value)),
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
			self.font:draw{
				pos={xmin, ymax},
				text=('%s [%.3e, %.3e]'):format(varName, thisValueMin, thisValueMax),
				color = {1,1,1,1},
				fontSize={fontSizeX, -fontSizeY},
				multiLine=false,
			}
		end
	end
end

end
