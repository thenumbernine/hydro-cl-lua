local class = require 'ext.class'
local path = require 'ext.path'
local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'
local Draw = require 'hydro.draw.draw'


local DrawMeshHeatmap = class(Draw)

function DrawMeshHeatmap:drawSolverWithVar(var, heatMap2DShader)
	local solver = self.solver
	local app = solver.app
	solver:calcDisplayVarToTex(var)

	local tex = solver:getTex()
	tex:bind(0)
	tex:setParameter(gl.GL_TEXTURE_MAG_FILTER, app.displayBilinearTextures and gl.GL_LINEAR or gl.GL_NEAREST)

--[[ 110 fps: glVertexAttrib prim calls
	gl.glBegin(gl.GL_TRIANGLES)
	for i=0,solver.numGlVtxs-1 do
		gl.glVertexAttrib1f(heatMap2DShader.attrs.cellindex.loc, solver.glcellindex.v[i])
		gl.glVertexAttrib3f(heatMap2DShader.attrs.vtxcenter.loc, solver.glvtxcenters.v[i]:unpack())
		gl.glVertexAttrib3f(heatMap2DShader.attrs.vtx.loc, solver.glvtxs.v[i]:unpack())
	end
	gl.glEnd()
--]]
--[[ 130 fps: glBindBuffer / glVertexAttribPointer / glEnableVertexAttribArray
error("I once again need to straighten out the ctor/usage of GLProgram vs its attribute objects vs GLAttribute vs GLVertexArray
	solver.heatMap2DShader:setAttrs(solver.heatMapShaderAttrs)
	GLVertexArray:enableAttrs(solver.heatMapShaderAttrs)
	gl.glDrawArrays(gl.GL_TRIANGLES, 0, solver.numGlVtxs * 3)
	GLVertexArray:disableAttrs(solver.heatMapShaderAttrs)
--]]
-- [[ 150fps: glVertexArray
	solver.heatMap2DShader.vao:use()
	gl.glDrawArrays(gl.GL_TRIANGLES, 0, solver.numGlVtxs * 3)
	solver.heatMap2DShader.vao:useNone()
--]]
	tex:unbind(0)
end

function DrawMeshHeatmap:showDisplayVar(var, varName, ar)
	local solver = self.solver
	local app = solver.app
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end


	-- TODO move the var.heatmap at the top and the drawgradient at the bottom outside of the function
	-- and then move this if-condition outside as well
	-- and things will match up wit other draw routines still
	if solver.showValues then

		gl.glEnable(gl.GL_DEPTH_TEST)
	--	gl.glEnable(gl.GL_CULL_FACE)

		local shader = solver.heatMap2DShader
		shader:use()

		local gradientTex = app.gradientTex
		gradientTex:bind(1)

		self:setupDisplayVarShader(shader, var, valueMin, valueMax)

		-- useCoordMap is missing from DrawMeshHeatmap
		if shader.uniforms.useLog then
			gl.glUniform1i(shader.uniforms.useLog.loc, 0)
		end
		if shader.uniforms.valueMin then
			gl.glUniform1f(shader.uniforms.valueMin.loc, valueMin)
		end
		if shader.uniforms.valueMax then
			gl.glUniform1f(shader.uniforms.valueMax.loc, valueMax)
		end
		-- drawCellScale isn't present in GridSolver
		gl.glUniform1f(shader.uniforms.drawCellScale.loc, solver.drawCellScale)

		-- this is only in DrawMeshHeatmap...
		gl.glUniformMatrix4fv(shader.uniforms.modelViewProjectionMatrix.loc, 1, 0, app.view.modelViewProjectionMatrix.ptr)

		gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
		gl.glEnable(gl.GL_BLEND)

		self:drawSolverWithVar(var, shader)

		gl.glDisable(gl.GL_BLEND)

		gradientTex:unbind(1)
		gl.glActiveTexture(gl.GL_TEXTURE0)
		shader:useNone()

		gl.glDisable(gl.GL_DEPTH_TEST)
	--	gl.glDisable(gl.GL_CULL_FACE)
	end

	app:drawGradientLegend(solver, var, varName, ar, valueMin, valueMax)
end

function DrawMeshHeatmap:display(varName, ar)
	local solver = self.solver
	local app = solver.app
	if app.targetSystem == 'console' then return end

	local view = app.view
	view:setup(ar)

	local var = solver.displayVarForName[varName]
	if not var then return end
	self:prepareShader()

	-- if it's a vector field then let app handle it.
	local component = solver.displayComponentFlatList[var.component]
	if solver:isVarTypeAVectorField(component.type) then return end

	gl.glEnable(gl.GL_DEPTH_TEST)

	-- TODO move this somewhere to work with all display methods of meshsolver
	if solver.showVertexes then
		gl.glPointSize(3)
		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_POINT)

		solver.drawPointsShader:use()
		gl.glUniformMatrix4fv(solver.drawPointsShader.uniforms.modelViewProjectionMatrix.loc, 1, 0, app.view.modelViewProjectionMatrix.ptr)
		gl.glUniform1f(solver.drawPointsShader.uniforms.drawCellScale.loc, solver.drawCellScale)

		solver.drawPointsShader.vao:use()
		gl.glDrawArrays(gl.GL_TRIANGLES, 0, solver.numGlVtxs * 3)
		solver.drawPointsShader.vao:useNone()

		solver.drawPointsShader:useNone()

		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
		gl.glPointSize(1)
	end

	-- TODO move this somewhere to work with all display methods of meshsolver
	if solver.showFaces then
		local mesh = solver.mesh

-- I can technically use the shader above
-- but it will include the internal edges of tesselated polygons
-- fixes? 1) make a separate arraybuffers for lines, 2) add an attribute for an internal edge flag
--[[
		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)

		solver.drawPointsShader:use()
		gl.glUniformMatrix4fv(solver.drawPointsShader.uniforms.modelViewProjectionMatrix.loc, 1, 0, app.view.modelViewProjectionMatrix.ptr)
		gl.glUniform1f(solver.drawPointsShader.uniforms.drawCellScale.loc, solver.drawCellScale)

		solver.drawPointsShader.vao:use()
		gl.glDrawArrays(gl.GL_TRIANGLES, 0, solver.numGlVtxs * 3)
		solver.drawPointsShader.vao:useNone()

		solver.drawPointsShader:useNone()

		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
--]]
-- [[ something between using the drawPointsShader and the GL 1.1 calls
		solver.drawPointsShader:use()
		gl.glUniformMatrix4fv(solver.drawPointsShader.uniforms.modelViewProjectionMatrix.loc, 1, 0, app.view.modelViewProjectionMatrix.ptr)
		gl.glUniform1f(solver.drawPointsShader.uniforms.drawCellScale.loc, solver.drawCellScale)
		for ci,c in ipairs(mesh.cells) do
			for fi=0,c.faceCount-1 do
				local f = mesh.faces.v[mesh.cellFaceIndexes.v[fi + c.faceOffset]]
				gl.glBegin(gl.GL_LINE_LOOP)
				for vi=0,f.vtxCount-1 do
					local v = mesh.vtxs.v[mesh.faceVtxIndexes.v[vi + f.vtxOffset]]
					gl.glVertexAttrib3f(solver.drawPointsShader.attrs.vtxcenter.loc, c.pos:unpack())
					gl.glVertexAttrib3f(solver.drawPointsShader.attrs.vtx.loc, v:unpack())
				end
				gl.glEnd()
			end
		end
		solver.drawPointsShader:useNone()
--]]
--[=[
		for fi,f in ipairs(mesh.faces) do
			gl.glBegin(gl.GL_LINE_LOOP)
			for vi=0,f.vtxCount-1 do
				local v = mesh.vtxs.v[mesh.faceVtxIndexes.v[vi + f.vtxOffset]]
				gl.glVertex3d(v:unpack())
			end
			gl.glEnd()
		end
--]=]
	end

	-- TODO move this somewhere to work with all display methods of meshsolver
	if solver.showNormals then
		local mesh = solver.mesh
		gl.glBegin(gl.GL_LINES)
		for ci,c in ipairs(mesh.cells) do
			for fi=0,c.faceCount-1 do
				local f = mesh.faces.v[mesh.cellFaceIndexes.v[fi + c.faceOffset]]
				local pos = (f.pos - c.pos) * solver.drawCellScale + c.pos
				local dx = .1 * math.sqrt(f.area) * solver.drawCellScale

				gl.glColor3f(1,0,0)
				gl.glVertex3d(pos:unpack())
				gl.glVertex3d((pos + f.normal * dx):unpack())

				gl.glColor3f(0,1,0)
				gl.glVertex3d(pos:unpack())
				gl.glVertex3d((pos + f.normal2 * dx):unpack())

				gl.glColor3f(0,0,1)
				gl.glVertex3d(pos:unpack())
				gl.glVertex3d((pos + f.normal3 * dx):unpack())
			end
		end
		gl.glEnd()
	end

	gl.glDisable(gl.GL_DEPTH_TEST)

	-- from here on it's showDisplayVar
	self:showDisplayVar(var, varName, ar)
	glreport'here'
end

function DrawMeshHeatmap:prepareShader()
	local solver = self.solver
	if solver.heatMap2DShader then return end

	local heatMapCode = assert(path'hydro/draw/mesh_heatmap.glsl':read())

	solver.heatMap2DShader = solver.GLProgram{
		name = 'mesh_heatmap',
		vertexCode = solver.eqn:template(heatMapCode, {
			draw = self,
			vertexShader = true,
		}),
		fragmentCode = solver.eqn:template(heatMapCode, {
			draw = self,
			fragmentShader = true,
		}),
		uniforms = {
			valueMin = 0,
			valueMax = 0,
			tex = 0,
			gradientTex = 1,
			drawCellScale = 1,
		},
		attrs = {
			vtx = solver.glvtxArrayBuffer,
			vtxcenter = solver.glvtxcenterArrayBuffer,
			cellindex = solver.glcellindexArrayBuffer,
		},
	}
end

return DrawMeshHeatmap
