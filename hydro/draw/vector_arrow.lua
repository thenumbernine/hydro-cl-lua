local ffi = require 'ffi'
local path = require 'ext.path'
local vec2f = require 'vec-ffi.vec2f'
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLSceneObject = require 'gl.sceneobject'
local Draw = require 'hydro.draw.draw'


local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
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
			count = 0,--arrowCount * #arrow,
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

	local isMeshSolver = require 'hydro.solver.meshsolver':isa(solver)

	local step = self.step
	local arrowCount
	local icount, jcount, kcount
	if isMeshSolver then
		arrowCount = solver.numCells
	else
		icount = math.max(1, math.floor(tonumber(solver.sizeWithoutBorder.x) / step))
		jcount = math.max(1, math.floor(tonumber(solver.sizeWithoutBorder.y) / step))
		kcount = math.max(1, math.floor(tonumber(solver.sizeWithoutBorder.z) / step))
		arrowCount = icount * jcount * kcount
	end

	local sceneObj = solver.vectorArrowSceneObj
	local shader = sceneObj.program
	local uniforms = shader.uniforms

	local vertexGPU = sceneObj.attrs.vertex.buffer
	local vertexCPU = vertexGPU:beginUpdate()
	local gridCoordGPU = sceneObj.attrs.gridCoord.buffer
	local gridCoordCPU = gridCoordGPU:beginUpdate()

	vertexCPU:resize(arrowCount * #arrow)
	gridCoordCPU:resize(arrowCount * #arrow)

	-- glCallOrDraw goes just slightly faster.  24 vs 23 fps.
	if isMeshSolver then
		-- TODO Lua coroutine cell iterator, abstracted between grids and meshes?
		-- how fast/slow are coroutines compared to number for-loops anyways?
		local pc = gridCoordCPU.v
		local pv = vertexCPU.v
		for ci=0,solver.numCells-1 do
			local c = solver.mesh.cells.v[ci]
			for _,q in ipairs(arrow) do
				pc[0].x = c.pos.x
				pc[0].y = c.pos.y
				pc[0].z = c.pos.z
				pc = pc + 1
				pv[0].x = q[1]
				pv[0].y = q[2]
				pv = pv + 1
			end
		end
	else
		local pc = gridCoordCPU.v
		local pv = vertexCPU.v
		for kbase=0,kcount-1 do
			local k = kbase * step
			for jbase=0,jcount-1 do
				local j = jbase * step
				for ibase=0,icount-1 do
					local i = ibase * step

					for _,q in ipairs(arrow) do
						pc[0].x = i
						pc[0].y = j
						pc[0].z = k
						pc = pc + 1
						pv[0].x = q[1]
						pv[0].y = q[2]
						pv = pv + 1
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
