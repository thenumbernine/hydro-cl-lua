local ffi = require 'ffi'
local class = require 'ext.class'
local file = require 'ext.file'
local template = require 'template'
local vec2f = require 'vec-ffi.vec2f'
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local glreport = require 'gl.report'
local GLVertexArray = require 'gl.vertexarray'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLAttribute = require 'gl.attribute'
local vector = require 'hydro.util.vector'
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
DrawVectorField.scale = 1
DrawVectorField.step = cmdline.vectorFieldStep or 4

function DrawVectorField:showDisplayVar(app, solver, var, varName, ar, xmin, xmax, ymin, ymax, useLog)

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

	local isMeshSolver = require 'hydro.solver.meshsolver'.is(solver)
	
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
	if not self.glvtxs then
		self.glvtxs = vector'vec2f_t'
	end
	if not self.glcenters then
		self.glcenters = vector'vec3f_t'
	end

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
			local index = 0
			for ci=0,solver.numCells-1 do
				local c = solver.mesh.cells.v[ci]
				for _,q in ipairs(arrow) do
					self.glcenters.v[index]:set(c.pos:unpack())
					self.glvtxs.v[index]:set(q[1], q[2])
					index = index + 1
				end
			end
		else
			local index = 0
			for kbase=0,kcount-1 do
				local k = kbase * step
				local z = (k + .5) / tonumber(solver.sizeWithoutBorder.z) * (solver.maxs.z - solver.mins.z) + solver.mins.z
				for jbase=0,jcount-1 do
					local j = jbase * step
					local y = (j + .5) / tonumber(solver.sizeWithoutBorder.y) * (solver.maxs.y - solver.mins.y) + solver.mins.y
					for ibase=0,icount-1 do
						local i = ibase * step
						local x = (i + .5) / tonumber(solver.sizeWithoutBorder.x) * (solver.maxs.x - solver.mins.x) + solver.mins.x
						for _,q in ipairs(arrow) do
							self.glcenters.v[index]:set(x,y,z)
							self.glvtxs.v[index]:set(q[1], q[2])
							index = index + 1
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
				center = solver.vectorArrowGLCentersArrayBuffer,
			},
		}
	end

	gl.glEnable(gl.GL_BLEND)

	shader:use()
	self:setupDisplayVarShader(shader, app, solver, var, valueMin, valueMax)

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
		gl.glVertexAttrib3f(shader.attrs.center.loc, solver.vectorArrowGLCentersArrayBuffer.data[i]:unpack())
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

function DrawVectorField:display(app, solvers, varName, ar, ...)
	app.view:setup(ar)

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	
	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:prepareShader(solver)
			self:showDisplayVar(app, solver, var, varName, ar, ...)
		end
	end

end

function DrawVectorField:prepareShader(solver)
	if solver.vectorArrowShader then return end

	local vectorArrowCode = assert(file['hydro/draw/vector_arrow.shader'])
	
	solver.vectorArrowShader = solver.GLProgram{
		name = 'vector_arrow',
		vertexCode = template(vectorArrowCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			vertexShader = true,
		}),
		fragmentCode = template(vectorArrowCode, {
			draw = self,
			app = solver.app,
			solver = solver,
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
