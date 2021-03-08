local gl = require 'ffi.OpenGL'
local ffi = require 'ffi'
local class = require 'ext.class'
local vec3f = require 'vec-ffi.vec3f'
local file = require 'ext.file'
local matrix_ffi = require 'matrix.ffi'
local vector = require 'hydro.util.vector'
local Draw = require 'hydro.draw.draw'


local DrawVectorStateLine = class(Draw)

function DrawVectorStateLine:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = self.solver.app
	
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		local component = solver.displayComponentFlatList[var.component]
		local vectorField = solver:isVarTypeAVectorField(component.type)
		if vectorField then
			valueMin, valueMax = solver:calcDisplayVarRange(var, component.magn)
		else
			valueMin, valueMax = solver:calcDisplayVarRange(var)
		end
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end
	
	solver:calcDisplayVarToTex(var)	


	local shader = solver.drawVectorStateLineShader
	local uniforms = shader.uniforms

	shader:use()
	local tex = solver:getTex(var)
	tex:bind(0)
	app.gradientTex:bind(1)
	
	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	if not self.vertexes then self.vertexes = vector'vec3f_t' end

	local step = 1
	local numVertexes = math.floor((tonumber(solver.sizeWithoutBorder.x) + 1) / step)	-- (endindex - startindex + 1) / step
	if #self.vertexes ~= numVertexes then
		self.vertexes:resize(numVertexes)
	end

	-- [[ overwrite the modelViewProjectionMatrix uniform here
	-- this is different from the other 'Draw.modelViewProjectionMatrix'
	-- that one is based on hydro.view, 
	-- this is based on the GL state set in hydro.app for 1D graphs
	-- TODO maybe combine the two, make the hydro.app 1D graph stuff use hydro.view.ortho, 
	-- then this could just use the default 'modelViewProjectionMatrix'
	self.ModelViewMatrix = self.ModelViewMatrix or matrix_ffi(nil, 'float', {4,4})--ffi.new'float[16]'
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.ModelViewMatrix.ptr)

	self.ProjectionMatrix = self.ProjectionMatrix or matrix_ffi(nil, 'float', {4,4})--ffi.new'float[16]'
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.ProjectionMatrix.ptr)

	self.ModelViewProjectionMatrix = self.ModelViewProjectionMatrix or require 'matrix.ffi'(nil, 'float', {4,4})--ffi.new'float[16]'
	matrix_ffi.mul(self.ModelViewProjectionMatrix, self.ProjectionMatrix, self.ModelViewMatrix)

	gl.glUniformMatrix4fv(uniforms.modelViewProjectionMatrix.loc, 1, gl.GL_FALSE, self.ModelViewProjectionMatrix.ptr)
	--]]

	-- TODO don't do this at all.  just somehow get around it.  idk, geometry shaders or something.
	for i=0,numVertexes-1 do
		local v = self.vertexes.v[i]
		v.x = i * step
		v.y = 0--app.displayFixedY
		v.z = 0--app.displayFixedZ
	end
	
	gl.glEnableVertexAttribArray(shader.attrs.gridCoord.loc)
	gl.glVertexAttribPointer(shader.attrs.gridCoord.loc, 3, gl.GL_FLOAT, false, 0, self.vertexes.v)
	
	gl.glDrawArrays(gl.GL_LINE_STRIP, 0, numVertexes) 
	
	gl.glDisableVertexAttribArray(shader.attrs.gridCoord.loc)
	
	app.gradientTex:unbind(1)
	tex:unbind(0)
	shader:useNone()
	

	-- TODO only draw the first
	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

function DrawVectorStateLine:display(varName, ar, ...)
	local solver = self.solver
	local app = solver.app
	app.view:setup(ar)

	gl.glColor3f(.1, .1, .1)
	gl.glBegin(gl.GL_LINES)
	local max = 10 	-- TODO max of solver mins/maxs
	gl.glVertex3f(max,0,0)
	gl.glVertex3f(-max,0,0)
	gl.glVertex3f(0,max,0)
	gl.glVertex3f(0,-max,0)
	gl.glVertex3f(0,0,max)
	gl.glVertex3f(0,0,-max)
	gl.glEnd()
	gl.glColor3f(1,1,1)

	-- display here
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var, varName, ar, ...)
	end
end

-- also in 2d_graph.lua.  subclass?
function DrawVectorStateLine:prepareShader()
	local solver = self.solver
	if solver.drawVectorStateLineShader then return end
	
	local code = assert(file['hydro/draw/vector_state_line.shader'])
	
	solver.drawVectorStateLineShader = solver.GLProgram{
		name = 'vector_state_line',
		vertexCode = solver.eqn:template(code, {
			draw = self,
			vertexShader = true,
		}),
		fragmentCode = solver.eqn:template(code, {
			draw = self,
			fragmentShader = true,
		}),
		uniforms = {
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
		},
	}
end

return DrawVectorStateLine
