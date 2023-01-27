local ffi = require 'ffi'
local class = require 'ext.class'
local file = require 'ext.file'
local vec2f = require 'vec-ffi.vec2f'
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local GLVertexArray = require 'gl.vertexarray'
local GLArrayBuffer = require 'gl.arraybuffer'
local vector = require 'ffi.cpp.vector'
local Draw = require 'hydro.draw.draw'


local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
}

local DrawVectorField = class(Draw)


-- TODO move to draw/vectorfield
DrawVectorField.scale = cmdline.vectorFieldScale or 1
DrawVectorField.step = cmdline.vectorFieldStep or 4

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

	local shader = solver.vectorArrowShader
	local uniforms = shader.uniforms

	-- TODO store these vectors per-solver?
	-- nah, once the capacity grows to the largest solver's size, it will stay there.
	if not self.glvtxs then self.glvtxs = vector'vec2f_t' end
	if not self.glcenters then self.glcenters = vector'vec3f_t' end

	if not solver.vectorArrowGLVtxArrayBuffer
	-- assert that glvtxs and glcenters are the same size
	or #self.glvtxs ~= arrowCount * #arrow
	then
		self.glvtxs:resize(arrowCount * #arrow)
		self.glcenters:resize(arrowCount * #arrow)

		-- glCallOrDraw goes just slightly faster.  24 vs 23 fps.
		if isMeshSolver then
			-- TODO Lua coroutine cell iterator, abstracted between grids and meshes?
			-- how fast/slow are coroutines compared to number for-loops anyways?
			local pc = self.glcenters.v
			local pv = self.glvtxs.v
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
			local pc = self.glcenters.v
			local pv = self.glvtxs.v
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

		solver.vectorArrowGLVtxArrayBuffer = GLArrayBuffer{
			data = self.glvtxs.v,
			size = #self.glvtxs * ffi.sizeof(self.glvtxs.type)
		}

		solver.vectorArrowGLCentersArrayBuffer = GLArrayBuffer{
			data = self.glcenters.v,
			size = #self.glcenters * ffi.sizeof(self.glcenters.type)
		}

		solver.vectorArrowVAO = GLVertexArray{
			program = shader,
			attrs = {
				vtx = solver.vectorArrowGLVtxArrayBuffer,
				gridCoord = solver.vectorArrowGLCentersArrayBuffer,
				cellindex = solver.glcellindexArrayBuffer,
			},
		}
	end

	gl.glEnable(gl.GL_BLEND)

	shader:use()
	self:setupDisplayVarShader(shader, var, valueMin, valueMax)

	local tex = solver:getTex(var)
	tex:bind(0)
	app.gradientTex:bind(1)

	-- how to determine scale?
	--local scale = self.scale * (valueMax - valueMin)
	--local scale = self.scale / (valueMax - valueMin)
	local scale = self.scale * step * solver.mindx
	gl.glUniform1f(uniforms.scale.loc, scale)

--[[ glVertexAttrib prim calls
	gl.glBegin(gl.GL_LINES)
	for i=0,arrowCount * #arrow-1 do
		gl.glVertexAttrib3f(shader.attrs.gridCoord.loc, solver.vectorArrowGLCentersArrayBuffer.data[i]:unpack())
		gl.glVertexAttrib2f(shader.attrs.vtx.loc, solver.vectorArrowGLVtxArrayBuffer.data[i]:unpack())
	end
	gl.glEnd()
--]]
-- [[ glVertexArray
	solver.vectorArrowVAO:use()
	gl.glDrawArrays(gl.GL_LINES, 0, arrowCount * #arrow)
	solver.vectorArrowVAO:useNone()
--]]
--[[ glVertexArray with glDrawArraysInstanced (not fully implemented - just speed testing) doesn't go noticably faster  than glDrawArrays
	solver.vectorArrowVAO:use()
	gl.glDrawArraysInstancedARB(gl.GL_LINES, 0, #arrow, arrowCount)
	solver.vectorArrowVAO:useNone()
--]]

	app.gradientTex:unbind(1)
	tex:unbind(0)
	shader:useNone()

	gl.glDisable(gl.GL_BLEND)


	-- TODO only draw the first
	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

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
	if solver.vectorArrowShader then return end

	local vectorArrowCode = assert(file'hydro/draw/vector_arrow.shader':read())

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
	}
end

return DrawVectorField
