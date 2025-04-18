local gl = require 'gl'
local vec3f = require 'vec-ffi.vec3f'
local path = require 'ext.path'
local matrix_ffi = require 'matrix.ffi'
local vector = require 'ffi.cpp.vector-lua'
local Draw = require 'hydro.draw.draw'


local DrawVectorStateLine = Draw:subclass()

function DrawVectorStateLine:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = self.solver.app


	gl.glEnable(gl.GL_DEPTH_TEST)

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
	--app.gradientTex:bind(1)
	--gl.glActiveTexture(gl.GL_TEXTURE0)

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.ambient.loc, app.displayDim == 1 and 1 or .3)

	if not self.vertexes then self.vertexes = vector'vec3f_t' end

	local step = 1
	local numX = math.floor((tonumber(solver.sizeWithoutBorder.x) + 1) / step)
	local numY = math.floor((tonumber(solver.sizeWithoutBorder.y) + 1) / step)
	local numVertexes = 2 * numX * numY
	if #self.vertexes ~= numVertexes then
		self.vertexes:resize(numVertexes)
	end

	-- [[ overwrite the mvProjMat uniform here
	-- this is different from the other 'Draw.mvProjMat'
	-- that one is based on hydro.view,
	-- this is based on the GL state set in hydro.app for 1D graphs
	-- TODO maybe combine the two, make the hydro.app 1D graph stuff use hydro.view.ortho,
	-- then this could just use the default 'mvProjMat'
	self.mvMat = self.mvMat or matrix_ffi(nil, 'float', {4,4})
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.mvMat.ptr)

	self.projMat = self.projMat or matrix_ffi(nil, 'float', {4,4})
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projMat.ptr)

	self.mvProjMat = self.mvProjMat or matrix_ffi(nil, 'float', {4,4})
	matrix_ffi.mul(self.mvProjMat, self.projMat, self.mvMat)

	gl.glUniformMatrix4fv(uniforms.mvProjMat.loc, 1, gl.GL_FALSE, self.mvProjMat.ptr)
	--]]

	-- TODO don't do this at all.  just somehow get around it.  idk, geometry shaders or something.
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

	--gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
	--gl.glDrawArrays(gl.GL_LINE_STRIP, 0, numVertexes)
	for j=0,numY-2 do
		--gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, j * 2 * numX, 2 * numX)
		gl.glDrawArrays(gl.GL_LINE_STRIP, j * 2 * numX, 2 * numX)
	end
	--gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)

	gl.glDisableVertexAttribArray(shader.attrs.gridCoord.loc)

	--app.gradientTex:unbind(1)
	tex:unbind(0)
	shader:useNone()

	gl.glDisable(gl.GL_DEPTH_TEST)


	-- TODO only draw the first
	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

function DrawVectorStateLine:display(varName, ar, xmin, ymin, xmax, ymax, ...)
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

	self:drawGrid(xmin, ymin, xmax, ymax)

	-- display here
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var, varName, ar, xmin, ymin, xmax, ymax, ...)
	end
end

-- also in 2d_graph.lua.  subclass?
function DrawVectorStateLine:prepareShader()
	local solver = self.solver
	if solver.drawVectorStateLineShader then return end

	local code = assert(path'hydro/draw/vector_state_line.glsl':read())

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
			--gradientTex = 1,
		},
	}:useNone()
end

return DrawVectorStateLine
