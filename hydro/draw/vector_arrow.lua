local ffi = require 'ffi'
local path = require 'ext.path'
local vec2f = require 'vec-ffi.vec2f'
local gl = require 'gl'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLSceneObject = require 'gl.sceneobject'
local Draw = require 'hydro.draw.draw'


local arrow = {
	vec2f(-.5, 0.),
	vec2f(.5, 0.),
	vec2f(.2, .3),
	vec2f(.5, 0.),
	vec2f(.2, -.3),
	vec2f(.5, 0.),
}

local DrawVectorField = Draw:subclass()


-- TODO move to draw/vectorfield
DrawVectorField.scale = cmdline.vectorFieldScale or 1
DrawVectorField.step = cmdline.vectorFieldStep or 4

function DrawVectorField:display(varName, ar, ...)
	local solver = self.solver
	local app = solver.app

	app.view:setup(ar)

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)

	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var, varName, ar, ...)
	end
end

function DrawVectorField:prepareShader()
	local solver = self.solver

	if not solver.vectorArrowShader then
		local vectorArrowCode = assert(path'hydro/draw/vector_arrow.glsl':read())

		solver.vectorArrowShader = solver.GLProgram{
			name = 'vector_arrow',
			vertexCode = solver.eqn:template(vectorArrowCode, {
				draw = self,
				vertexShader = true,
			}),
			fragmentCode = solver.eqn:template(vectorArrowCode, {
				draw = self,
				fragmentShader = true,
			}),
			uniforms = {
				scale = 1,
				valueMin = 0,
				valueMax = 0,
				tex = 0,
				gradientTex = 1,
			},
		}:useNone()
	end

	solver.vectorArrowSceneObj = solver.vectorArrowSceneObj or GLSceneObject{
		program = solver.vectorArrowShader,
		vertexes = {
			useVec = true,
			dim = 2,
		},
		attrs = {
			gridCoord = {
				buffer = {
					useVec = true,
					dim = 3,
				},
			},
		},
		geometry = {
			mode = gl.GL_LINES,
			count = 0,
		},
	}
end

function DrawVectorField:showDisplayVar(var, varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = solver.app

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

	local step = self.step

	local sceneObj = solver.vectorArrowSceneObj
	local shader = sceneObj.program
	local uniforms = shader.uniforms

	local vertexGPU = sceneObj.attrs.vertex.buffer
	local vertexCPU = vertexGPU:beginUpdate()
	local gridCoordGPU = sceneObj.attrs.gridCoord.buffer
	local gridCoordCPU = gridCoordGPU:beginUpdate()

	if require 'hydro.solver.meshsolver':isa(solver) then
		-- TODO Lua coroutine cell iterator, abstracted between grids and meshes?
		-- how fast/slow are coroutines compared to number for-loops anyways?
		for ci=0,solver.numCells-1 do
			local c = solver.mesh.cells.v[ci]
			for _,q in ipairs(arrow) do
				gridCoordCPU:emplace_back()[0] = c.pos
				vertexCPU:emplace_back()[0] = q
			end
		end
	else
		local icount = math.max(1, math.floor(tonumber(solver.sizeWithoutBorder.x) / step))
		local jcount = math.max(1, math.floor(tonumber(solver.sizeWithoutBorder.y) / step))
		local kcount = math.max(1, math.floor(tonumber(solver.sizeWithoutBorder.z) / step))

		local pc = gridCoordCPU.v
		local pv = vertexCPU.v
		for kbase=0,kcount-1 do
			local k = kbase * step
			for jbase=0,jcount-1 do
				local j = jbase * step
				for ibase=0,icount-1 do
					local i = ibase * step
					for _,q in ipairs(arrow) do
						local pc = gridCoordCPU:emplace_back()
						pc.x = i
						pc.y = j
						pc.z = k
						vertexCPU:emplace_back()[0] = q
					end
				end
			end
		end
	end

	gridCoordGPU:endUpdate()
	vertexGPU:endUpdate()
	sceneObj.geometry.count = #vertexCPU

	gl.glEnable(gl.GL_BLEND)

	shader:use()
	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	local tex = solver:getTex(var)
	tex:bind(0)
	app.gradientTex:bind(1)
	gl.glActiveTexture(gl.GL_TEXTURE0)

	-- how to determine scale?
	--local scale = self.scale * (valueMax - valueMin)
	--local scale = self.scale / (valueMax - valueMin)
	local scale = self.scale * step * solver.mindx
	gl.glUniform1f(uniforms.scale.loc, scale)

--[[ glVertexAttrib prim calls
	gl.glBegin(gl.GL_LINES)
	for i=0,arrowCount * #arrow-1 do
		gl.glVertexAttrib3f(shader.attrs.gridCoord.loc, gridCoordCPU.v[i]:unpack())
		gl.glVertexAttrib2f(shader.attrs.vertex.loc, vertexCPU.v[i]:unpack())
	end
	gl.glEnd()
--]]
-- [[ VAO if present ...
	sceneObj:enableAndSetAttrs()
	sceneObj.geometry:draw()
	sceneObj:disableAttrs()
--]]
--[[ glVertexArray with glDrawArraysInstanced (not fully implemented - just speed testing) doesn't go noticably faster than glDrawArrays
	solver.vectorArrowVAO:bind()
	gl.glDrawArraysInstanced(gl.GL_LINES, 0, #arrow, arrowCount)
	solver.vectorArrowVAO:unbind()
--]]

	app.gradientTex:unbind(1)
	tex:unbind(0)
	shader:useNone()

	gl.glDisable(gl.GL_BLEND)

	-- TODO only draw the first
	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

return DrawVectorField
