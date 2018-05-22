local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'

return function(HydroCLApp)

local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
}
HydroCLApp.displayVectorField_scale = 1
HydroCLApp.displayVectorField_step = 4
function HydroCLApp:displayVectorField(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)

	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			local magVar = assert(var.magVar, "tried to use a vector display on a var without an associated magVar")
			
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = solver:calcDisplayVarRange(magVar)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end
			
			solver:calcDisplayVarToTex(var)	
			
			solver.vectorFieldShader:use()
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
			self.gradientTex:bind(1)
			
			gl.glUniform3f(solver.vectorFieldShader.uniforms.mins.loc, solver.mins:unpack())
			gl.glUniform3f(solver.vectorFieldShader.uniforms.maxs.loc, solver.maxs:unpack())
			-- how to determine scale?
			--local scale = self.displayVectorField_scale * (valueMax - valueMin)
			--local scale = self.displayVectorField_scale / (valueMax - valueMin)
			local scale = self.displayVectorField_scale 
				* self.displayVectorField_step 
				* math.min(
					(solver.maxs[1] - solver.mins[1]) / tonumber(solver.gridSize.x),
					(solver.maxs[2] - solver.mins[2]) / tonumber(solver.gridSize.y),
					(solver.maxs[3] - solver.mins[3]) / tonumber(solver.gridSize.z))
			gl.glUniform1f(solver.vectorFieldShader.uniforms.scale.loc, scale) 

			local step = self.displayVectorField_step
	
			--[[ goes just slightly faster.  24 vs 23 fps.
			self.vectorField_displayList = self.vectorField_displayList or {}
			local glCallOrDraw = require 'gl.call'
			glCallOrDraw(self.vectorField_displayList, function()
			--]]	
				gl.glBegin(gl.GL_LINES)
				for k=0,tonumber(solver.sizeWithoutBorder.z-1),step do
					for j=0,tonumber(solver.sizeWithoutBorder.y-1),step do
						for i=0,tonumber(solver.sizeWithoutBorder.x-1),step do
							local tx = (i + .5 + solver.numGhost) / tonumber(solver.gridSize.x)
							local ty = (j + .5 + (solver.dim > 1 and solver.numGhost or 0)) / tonumber(solver.gridSize.y)
							local tz = (k + .5 + (solver.dim > 2 and solver.numGhost or 0)) / tonumber(solver.gridSize.z)
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

			self.gradientTex:unbind(1)
			solver:getTex(var):unbind(0)
			solver.vectorFieldShader:useNone()
		end
	end
	
	gl.glDisable(gl.GL_DEPTH_TEST)
	
glreport'here'
end

end
