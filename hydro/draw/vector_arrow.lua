local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'
local class = require 'ext.class'

local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
}

local DrawVectorField = class()


-- TODO move to draw/vectorfield 
DrawVectorField.scale = 1
DrawVectorField.step = cmdline.vectorFieldStep or 4

function DrawVectorField:showDisplayVar(app, solver, var, xmin, xmax, ymin, ymax, useLog)
	
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
	
	solver.vectorArrowShader:use()
	gl.glUniform1i(solver.vectorArrowShader.uniforms.displayDim.loc, app.displayDim)
	gl.glUniform1i(solver.vectorArrowShader.uniforms.useLog.loc, var.useLog)
	-- [[ this gives the l1 bounds of the vector field
	gl.glUniform1f(solver.vectorArrowShader.uniforms.valueMin.loc, valueMin)
	gl.glUniform1f(solver.vectorArrowShader.uniforms.valueMax.loc, valueMax)
	--]]
	--[[ it'd be nice instead to get the norm bounds ... 
	-- but looking at the reduce calculations, the easiest way to do that is
	-- to associate each vector display shader with a norm display shader
	-- and then just reduce that
	gl.glUniform1f(solver.vectorArrowShader.uniforms.valueMin.loc, 0)
	gl.glUniform1f(solver.vectorArrowShader.uniforms.valueMax.loc, math.max(math.abs(valueMin), math.abs(valueMax)))
	--]]
	solver:getTex(var):bind(0)
	app.gradientTex:bind(1)
	
	gl.glUniform3f(solver.vectorArrowShader.uniforms.solverMins.loc, solver.mins:unpack())
	gl.glUniform3f(solver.vectorArrowShader.uniforms.solverMaxs.loc, solver.maxs:unpack())
	-- how to determine scale?
	--local scale = self.scale * (valueMax - valueMin)
	--local scale = self.scale / (valueMax - valueMin)
	local scale = self.scale 
		* self.step 
		* math.min(
			(solver.maxs.x - solver.mins.x) / tonumber(solver.texSize.x),
			(solver.maxs.y - solver.mins.y) / tonumber(solver.texSize.y),
			(solver.maxs.z - solver.mins.z) / tonumber(solver.texSize.z))
	gl.glUniform1f(solver.vectorArrowShader.uniforms.scale.loc, scale) 

	local step = self.step

	-- glCallOrDraw goes just slightly faster.  24 vs 23 fps.
	if require 'hydro.solver.meshsolver'.is(solver) then
		gl.glBegin(gl.GL_LINES)
		-- TODO Lua coroutine cell iterator, abstracted between grids and meshes?
		-- how fast/slow are coroutines compared to number for-loops anyways?
		for ci=0,solver.numCells-1 do
			local c = solver.mesh.cells.v[ci]
			local tx = (ci + .5) / tonumber(solver.numCells)
			gl.glMultiTexCoord3f(gl.GL_TEXTURE0, tx, .5, .5)
			gl.glMultiTexCoord3f(gl.GL_TEXTURE1, c.pos:unpack())
			for _,q in ipairs(arrow) do
				gl.glVertex2f(q[1], q[2])
			end
		end
		gl.glEnd()
	else
		gl.glBegin(gl.GL_LINES)
		for k=0,tonumber(solver.sizeWithoutBorder.z-1),step do
			for j=0,tonumber(solver.sizeWithoutBorder.y-1),step do
				for i=0,tonumber(solver.sizeWithoutBorder.x-1),step do
					local tx = (i + .5 + solver.numGhost) / tonumber(solver.texSize.x)
					local ty = (j + .5 + (solver.dim > 1 and solver.numGhost or app.displayFixedY * tonumber(solver.texSize.z))) / tonumber(solver.texSize.y)
					local tz = (k + .5 + (solver.dim > 2 and solver.numGhost or app.displayFixedZ * tonumber(solver.texSize.z))) / tonumber(solver.texSize.z)
					gl.glMultiTexCoord3f(gl.GL_TEXTURE0, tx, ty, tz)	
					local x = (i + .5) / tonumber(solver.sizeWithoutBorder.x)
					local y = (j + .5) / tonumber(solver.sizeWithoutBorder.y)
					local z = (k + .5) / tonumber(solver.sizeWithoutBorder.z)
					gl.glMultiTexCoord3f(gl.GL_TEXTURE1, x, y, z)
					for _,q in ipairs(arrow) do
						gl.glVertex2f(q[1], q[2])
					end
				end
			end
		end
		gl.glEnd()
	end

	app.gradientTex:unbind(1)
	solver:getTex(var):unbind(0)
	solver.vectorArrowShader:useNone()
end

function DrawVectorField:display(app, solvers, varName, ar, ...)
	app.view:projection(ar)
	app.view:modelview()

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	
	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:showDisplayVar(app, solver, var, ...)
		end
	end

	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:displayVector_Arrows(...)
		if not self.drawVectorArrows then self.drawVectorArrows = DrawVectorField() end
		return self.drawVectorArrows:display(self, ...)
	end
end
