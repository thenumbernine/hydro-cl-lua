local gl = require 'gl'
local GLSceneObject = require 'gl.sceneobject'
local Draw = require 'hydro.draw.draw'


local Draw1D = Draw:subclass()

function Draw1D:display(varName, ar, xmin, xmax, ymin, ymax, useLog, valueMin, valueMax)
	local solver = self.solver
	local app = solver.app

	-- trust that app.view is already setup ...

	-- TODO one grid for all displaly instead of multiple calls...
	local ystep = self:drawGrid(xmin, xmax, ymin, ymax)

	local view = app.view
	view.mvMat:setIdent()
	view.projMat:setOrtho(xmin, xmax, ymin, ymax, -1, 1)
	view.mvProjMat:copy(view.projMat)

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
		--self.orthoView:setupProjection(ar)
		--self.orthoView:setupModelView()
		--local xmin, xmax, ymin, ymax = self.orthoView:getOrthoBounds(ar)

		local font = app.font
		if font then
			font.view.mvMat:copy(view.mvMat)
			font.view.projMat:copy(view.projMat)
			font.view.mvProjMat:copy(view.mvProjMat)

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
				font:draw{
					pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
					text=(
						(absvalue > 1e+5 or absvalue < 1e-5)
						and ('%.5e'):format(value) or ('%.5f'):format(value)),
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
			font:draw{
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

function Draw1D:prepareShader()
	local solver = self.solver

	self:prepareGraphShader()

	solver.draw1DGraphSceneObj = solver.draw1DGraphSceneObj or GLSceneObject{
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

	local sceneObj = solver.draw1DGraphSceneObj

	local step = 1
	local numVertexes = math.floor((tonumber(solver.gridSize.x) - 2 * solver.numGhost + 1) / step)	-- (endindex - startindex + 1) / step

	local vertexGPU = sceneObj.attrs.vertex.buffer
	local vertexCPU = vertexGPU:beginUpdate()
	for i=0,numVertexes-1 do
		local v = vertexCPU:emplace_back()
		v.x = i * step
		v.y = 0--app.displayFixedY
		v.z = 0--app.displayFixedZ
	end
	vertexGPU:endUpdate()
	sceneObj.geometry.count = #vertexCPU

	local shader = sceneObj.program
	local uniforms = shader.uniforms

	shader:use()
	local tex = solver:getTex(var)
	tex:bind()

	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	gl.glUniform1f(uniforms.ambient.loc, 1)
	gl.glUniform3f(uniforms.color.loc, (#app.solvers > 1 and solver or var).color:unpack())

	-- [[ overwrite the mvProjMat uniform here
	-- this is different from the other 'Draw.mvProjMat'
	-- that one is based on hydro.view,
	-- this is based on the GL state set in hydro.app for 1D graphs
	-- TODO maybe combine the two, make the hydro.app 1D graph stuff use hydro.view.ortho,
	-- then this could just use the default 'mvProjMat'

	local view = app.view
	gl.glUniformMatrix4fv(uniforms.mvProjMat.loc, 1, gl.GL_TRUE, view.mvProjMat.ptr)
	--]]

	sceneObj:enableAndSetAttrs()
	sceneObj.geometry:draw()
	sceneObj:disableAttrs()

	tex:unbind()
	shader:useNone()
end

return Draw1D
