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

	local vectorArrowShader = solver.vectorArrowShader

	if not solver.vectorArrowGLVtxArrayBuffer
	or solver.vectorArrowGLBufSize ~= arrowCount
	then
		solver.vectorArrowGLBufSize = arrowCount
		
		local glvtxs = vector'vec2f_t'
		glvtxs:setcapacity(arrowCount * #arrow)
		local glcenters = vector'vec3f_t'
		glcenters:setcapacity(arrowCount * #arrow)
		local gltcs = vector'vec3f_t'
		gltcs:setcapacity(arrowCount * #arrow)

		-- glCallOrDraw goes just slightly faster.  24 vs 23 fps.
		if isMeshSolver then
			-- TODO Lua coroutine cell iterator, abstracted between grids and meshes?
			-- how fast/slow are coroutines compared to number for-loops anyways?
			for ci=0,solver.numCells-1 do
				local c = solver.mesh.cells.v[ci]
				for _,q in ipairs(arrow) do
					glcenters:push_back(vec3f(c.pos:unpack()))
					gltcs:push_back(vec3f((ci + .5) / tonumber(solver.numCells), .5, .5))
					glvtxs:push_back(vec2f(q[1], q[2]))
				end
			end
		else
			for kbase=0,kcount-1 do
				local k = kbase * step
				for jbase=0,jcount-1 do
					local j = jbase * step
					for ibase=0,icount-1 do
						local i = ibase * step
						for _,q in ipairs(arrow) do
							local x = vec3f(
								(i + .5) / tonumber(solver.sizeWithoutBorder.x) * (solver.maxs.x - solver.mins.x) + solver.mins.x,
								(j + .5) / tonumber(solver.sizeWithoutBorder.y) * (solver.maxs.y - solver.mins.y) + solver.mins.y,
								(k + .5) / tonumber(solver.sizeWithoutBorder.z) * (solver.maxs.z - solver.mins.z) + solver.mins.z)
							glcenters:push_back(x)
							gltcs:push_back(vec3f(
								(i + .5 + solver.numGhost) / tonumber(solver.texSize.x),
								solver.dim > 1 and (j + .5 + solver.numGhost) / tonumber(solver.texSize.y) or .5,
								solver.dim > 2 and (k + .5 + solver.numGhost) / tonumber(solver.texSize.z) or .5
							))
							glvtxs:push_back(vec2f(q[1], q[2]))
						end
					end
				end
			end
		end

		solver.vectorArrowGLVtxArrayBuffer = GLArrayBuffer{
			data = glvtxs.v,
			size = #glvtxs * ffi.sizeof(glvtxs.type)
		}

		solver.vectorArrowGLCentersArrayBuffer = GLArrayBuffer{
			data = glcenters.v,
			size = #glcenters * ffi.sizeof(glcenters.type)
		}
		
		solver.vectorArrowGLTCsArrayBuffer = GLArrayBuffer{
			data = gltcs.v,
			size = #gltcs * ffi.sizeof(gltcs.type)
		}
		
		solver.vectorArrowVAO = GLVertexArray{
			program = vectorArrowShader,
			attrs = {
				vtx = solver.vectorArrowGLVtxArrayBuffer,
				center = solver.vectorArrowGLCentersArrayBuffer,
				tc = solver.vectorArrowGLTCsArrayBuffer,
			},
		}
	end

	gl.glEnable(gl.GL_BLEND)

	vectorArrowShader:use()
	
	self:setupDisplayVarShader(vectorArrowShader, app, solver, var, valueMin, valueMax)
		
	gl.glUniform3f(vectorArrowShader.uniforms.solverMins.loc, solver.mins:unpack())
	gl.glUniform3f(vectorArrowShader.uniforms.solverMaxs.loc, solver.maxs:unpack())

	local tex = solver:getTex(var)
	tex:bind(0)
	app.gradientTex:bind(1)
	
	-- how to determine scale?
	--local scale = self.scale * (valueMax - valueMin)
	--local scale = self.scale / (valueMax - valueMin)
	local scale = self.scale * step * solver.mindx
	gl.glUniform1f(vectorArrowShader.uniforms.scale.loc, scale) 

--[[ glVertexAttrib prim calls
	gl.glBegin(gl.GL_LINES)
	for i=0,arrowCount * #arrow-1 do
		gl.glVertexAttrib3f(vectorArrowShader.attrs.tc.loc, solver.vectorArrowGLTCsArrayBuffer.data[i]:unpack())
		gl.glVertexAttrib3f(vectorArrowShader.attrs.center.loc, solver.vectorArrowGLCentersArrayBuffer.data[i]:unpack())
		gl.glVertexAttrib2f(vectorArrowShader.attrs.vtx.loc, solver.vectorArrowGLVtxArrayBuffer.data[i]:unpack())
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
	vectorArrowShader:useNone()

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
