-- maybe a better name would be '1d_graph'?
local gl = require 'gl'
local vec3f = require 'vec-ffi.vec3f'
local path = require 'ext.path'
local matrix_ffi = require 'matrix.ffi'
local vector = require 'ffi.cpp.vector-lua'
local Draw = require 'hydro.draw.draw'


local Draw1D = Draw:subclass()

function Draw1D:showDisplayVar(var)
	local solver = self.solver
	local app = self.solver.app

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
		return
	end


	local shader = solver.graphShader
	local uniforms = shader.uniforms

	shader:use()
	local tex = solver:getTex(var)
	tex:bind()

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.ambient.loc, 1)

	gl.glUniform3f(uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	if not self.vertexes then self.vertexes = vector'vec3f_t' end

	local step = 1
	local numVertexes = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)	-- (endindex - startindex + 1) / step
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

	gl.glUniformMatrix4fv(uniforms.mvProjMat.loc, 1, gl.GL_TRUE, self.mvProjMat.ptr)
	--]]

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

	tex:unbind()
	shader:useNone()
end

function Draw1D:display(varName, ar, xmin, xmax, ymin, ymax, useLog, valueMin, valueMax)
	local solver = self.solver
	local app = solver.app

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
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var)

		local unitScale = 1
		local valueMin = valueMin
		local valueMax = valueMax
		local varName = varName
		if var.showInUnits and var.units then
			unitScale = solver:convertToSIUnitsCode(var.units).func()
			valueMin = valueMin * unitScale
			valueMax = valueMax * unitScale
			varName = varName..' ('..var.units..')'
		end

		-- this has already been done
		-- but same with other draw objs, they are setting up the view multiple times
		-- the dif is, Draw1D doesn't use self.orthoView ...
		-- TODO fix that?
		-- also notice: ymin/ymax has already been log()'d
		--self.orthoView:projection(ar)
		--self.orthoView:modelview()
		--local xmin, xmax, ymin, ymax = self.orthoView:getOrthoBounds(ar)

		if app.font then
			-- gradient uses 0.025
			local fontSizeX = (xmax - xmin) * .05
			local fontSizeY = (ymax - ymin) * .05
			-- 1D:
			local ystep = ystep * 2
			-- Gradient:
			--local ystep = 10^(math.log(ymax - ymin, 10) - 1.5)
			for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
				-- 1D uses value = y
				local value = y
				-- Gradient linearly maps valueMin/Max to ymin/max
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
				-- 1D:
				text=('%s [%.3e, %.3e]'):format(varName, ymin, ymax),
				-- Gradient:
				--text = varName,
				color = {1,1,1,1},
				fontSize={fontSizeX, -fontSizeY},
				multiLine=false,
			}
		end
	end
end

-- also in 2d_graph.lua.  subclass?
function Draw1D:prepareShader()
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

return Draw1D
