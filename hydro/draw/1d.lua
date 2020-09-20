-- maybe a better name would be '1d_graph'?
local gl = require 'ffi.OpenGL'
local ffi = require 'ffi'
local template = require 'template'
local class = require 'ext.class'
local file = require 'ext.file'
local matrix_ffi = require 'matrix.ffi'
local Draw = require 'hydro.draw.draw'


local Draw1D = class(Draw)

function Draw1D:showDisplayVar(app, solver, var)
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


	-- 1D displays -- use vertex.y
	-- 2D displays -- use vertex.z
	-- 3D displays -- ???
	if app.displayDim == 3 then
		io.stderr:write'Why are you using a graph shader to display 3D data?  Use a 3D display instead.\n'
		do return end
	end


	local graphShader = solver.graphShader

	graphShader:use()
	local tex = solver:getTex(var)
	tex:bind()
	
	self:setupDisplayVarShader(graphShader, app, solver, var, valueMin, valueMax)

	gl.glUniform1f(graphShader.uniforms.scale.loc, 1)
	gl.glUniform1f(graphShader.uniforms.offset.loc, 0)
	gl.glUniform1f(graphShader.uniforms.ambient.loc, 1)
	gl.glUniform2f(graphShader.uniforms.xmin.loc, solver.mins.x, 0)
	gl.glUniform2f(graphShader.uniforms.xmax.loc, solver.maxs.x, 0)

	gl.glUniform2f(graphShader.uniforms.size.loc, solver.gridSize.x, solver.gridSize.y)
	gl.glUniform3f(graphShader.uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	local step = 1
	local numVertexes = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)	-- (endindex - startindex + 1) / step
	if numVertexes ~= self.numVertexes then
		self.numVertexes = numVertexes 
		self.vertexes = ffi.new('float[?]', 3*numVertexes)
	end
	
	self.ModelViewMatrix = self.ModelViewMatrix or matrix_ffi(nil, 'float', {4,4})--ffi.new'float[16]'
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.ModelViewMatrix.ptr)
	
	self.ProjectionMatrix = self.ProjectionMatrix or matrix_ffi(nil, 'float', {4,4})--ffi.new'float[16]'
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.ProjectionMatrix.ptr)

	self.ModelViewProjectionMatrix = self.ModelViewProjectionMatrix or require 'matrix.ffi'(nil, 'float', {4,4})--ffi.new'float[16]'
	matrix_ffi.mul(self.ModelViewProjectionMatrix, self.ProjectionMatrix, self.ModelViewMatrix)

	gl.glUniformMatrix4fv(graphShader.uniforms.ModelViewProjectionMatrix.loc, 1, gl.GL_FALSE, self.ModelViewProjectionMatrix.ptr)
	
	for i=0,self.numVertexes-1 do
		local x = (i * step + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
		self.vertexes[0+3*i] = x
		self.vertexes[1+3*i] = app.displayFixedY
		self.vertexes[2+3*i] = app.displayFixedZ
	end
	
	gl.glEnableVertexAttribArray(graphShader.attrs.inVertex.loc)
	gl.glVertexAttribPointer(graphShader.attrs.inVertex.loc, 3, gl.GL_FLOAT, false, 0, self.vertexes)
	
	gl.glDrawArrays(gl.GL_LINE_STRIP, 0, self.numVertexes) 
	
	gl.glDisableVertexAttribArray(graphShader.attrs.inVertex.loc)
	
	tex:unbind()
	graphShader:useNone()
end

function Draw1D:display(app, solvers, varName, ar, xmin, xmax, ymin, ymax, useLog, valueMin, valueMax)
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
			self:prepareShader(solver)
			self:showDisplayVar(app, solver, var)
		end
		
		local unitScale = 1
		local thisValueMin = valueMin
		local thisValueMax = valueMax
		if var.showInUnits and var.units then
			unitScale = solver:convertToSIUnitsCode(var.units).func()
			thisValueMin = thisValueMin * unitScale
			thisValueMax = thisValueMax * unitScale
		end
		
		if app.font then
			local fontSizeX = (xmax - xmin) * .05
			local fontSizeY = (ymax - ymin) * .05
			local ystep = ystep * 2
			for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
				local value = y
				if useLog then
					value = 10^value
				end
				value = value * unitScale
				local absvalue = math.abs(value)
				app.font:draw{
					pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
					text=(
						(absvalue > 1e+5 or absvalue < 1e-5)
						and ('%.5e'):format(value) or ('%.5f'):format(value)),
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
			app.font:draw{
				pos={xmin, ymax},
				text=('%s [%.3e, %.3e]'):format(varName, thisValueMin, thisValueMax),
				color = {1,1,1,1},
				fontSize={fontSizeX, -fontSizeY},
				multiLine=false,
			}
		end
	end
end

-- also in 2d_graph.lua.  subclass?
function Draw1D:prepareShader(solver)
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

return Draw1D
