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

function DrawVectorField:displayVectorField(app, solver, varName, xmin, ymin, xmax, ymax, useLog)
	if require 'solver.meshsolver'.is(solver) then return end
	
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
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
		
		solver.vectorFieldShader:use()
		gl.glUniform1i(solver.vectorFieldShader.uniforms.displayDim.loc, app.displayDim)
		gl.glUniform1i(solver.vectorFieldShader.uniforms.useLog.loc, var.useLog)
		-- [[ this gives the l1 bounds of the vector field
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMin.loc, valueMin)
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMax.loc, valueMax)
		--]]
		--[[ it'd be nice instead to get the norm bounds ... 
		-- but looking at the reduce calculations, the easiest way to do that is
		-- to associate each vector display shader with a norm display shader
		-- and then just reduce that
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMin.loc, 0)
		gl.glUniform1f(solver.vectorFieldShader.uniforms.valueMax.loc, math.max(math.abs(valueMin), math.abs(valueMax)))
		--]]
		solver:getTex(var):bind(0)
		app.gradientTex:bind(1)
		
		gl.glUniform3f(solver.vectorFieldShader.uniforms.mins.loc, solver.mins:unpack())
		gl.glUniform3f(solver.vectorFieldShader.uniforms.maxs.loc, solver.maxs:unpack())
		-- how to determine scale?
		--local scale = app.displayVectorField_scale * (valueMax - valueMin)
		--local scale = app.displayVectorField_scale / (valueMax - valueMin)
		local scale = app.displayVectorField_scale 
			* app.displayVectorField_step 
			* math.min(
				(solver.maxs.x - solver.mins.x) / tonumber(solver.gridSize.x),
				(solver.maxs.y - solver.mins.y) / tonumber(solver.gridSize.y),
				(solver.maxs.z - solver.mins.z) / tonumber(solver.gridSize.z))
		gl.glUniform1f(solver.vectorFieldShader.uniforms.scale.loc, scale) 

		local step = app.displayVectorField_step

		--[[ goes just slightly faster.  24 vs 23 fps.
		app.vectorField_displayList = app.vectorField_displayList or {}
		local glCallOrDraw = require 'gl.call'
		glCallOrDraw(app.vectorField_displayList, function()
		--]]	
			gl.glBegin(gl.GL_LINES)
			for k=0,tonumber(solver.sizeWithoutBorder.z-1),step do
				for j=0,tonumber(solver.sizeWithoutBorder.y-1),step do
					for i=0,tonumber(solver.sizeWithoutBorder.x-1),step do
						local tx = (i + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
						local ty = (j + .5 + (solver.dim > 1 and solver.numGhost or app.displayFixedY * tonumber(solver.gridSize.z))) / tonumber(solver.gridSize.y)
						local tz = (k + .5 + (solver.dim > 2 and solver.numGhost or app.displayFixedZ * tonumber(solver.gridSize.z))) / tonumber(solver.gridSize.z)
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
		--[[
		end)
		--]]

		app.gradientTex:unbind(1)
		solver:getTex(var):unbind(0)
		solver.vectorFieldShader:useNone()
	end
end

function DrawVectorField:display(app, solvers, ar, ...)
	app.view:projection(ar)
	app.view:modelview()

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	
	for _,solver in ipairs(solvers) do
		self:displayVectorField(app, solver, ...)
	end

	gl.glDisable(gl.GL_DEPTH_TEST)
end

return function(HydroCLApp)
	function HydroCLApp:displayVectorField(...)
		if not self.drawVectorField then self.drawVectorField = DrawVectorField() end
		return self.drawVectorField:display(self, ...)
	end
end
