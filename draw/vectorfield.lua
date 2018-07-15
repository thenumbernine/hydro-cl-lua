local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'

local arrow = {
	{-.5, 0.},
	{.5, 0.},
	{.2, .3},
	{.5, 0.},
	{.2, -.3},
	{.5, 0.},
}

local function applyToSolver(Solver)
	function Solver:displayVectorField(app, varName, xmin, ymin, xmax, ymax, useLog)
		local var = self.displayVarForName[varName]
		if var and var.enabled then
			local magVar = assert(var.magVar, "tried to use a vector display on a var without an associated magVar")
			
			local valueMin, valueMax
			if var.heatMapFixedRange then
				valueMin = var.heatMapValueMin
				valueMax = var.heatMapValueMax
			else
				valueMin, valueMax = self:calcDisplayVarRange(magVar)
				var.heatMapValueMin = valueMin
				var.heatMapValueMax = valueMax
			end
			
			self:calcDisplayVarToTex(var)	
			
			self.vectorFieldShader:use()
			gl.glUniform1i(self.vectorFieldShader.uniforms.useLog.loc, var.useLog)
			-- [[ this gives the l1 bounds of the vector field
			gl.glUniform1f(self.vectorFieldShader.uniforms.valueMin.loc, valueMin)
			gl.glUniform1f(self.vectorFieldShader.uniforms.valueMax.loc, valueMax)
			--]]
			--[[ it'd be nice instead to get the norm bounds ... 
			-- but looking at the reduce calculations, the easiest way to do that is
			-- to associate each vector display shader with a norm display shader
			-- and then just reduce that
			gl.glUniform1f(self.vectorFieldShader.uniforms.valueMin.loc, 0)
			gl.glUniform1f(self.vectorFieldShader.uniforms.valueMax.loc, math.max(math.abs(valueMin), math.abs(valueMax)))
			--]]
			self:getTex(var):bind(0)
			app.gradientTex:bind(1)
			
			gl.glUniform3f(self.vectorFieldShader.uniforms.mins.loc, self.mins:unpack())
			gl.glUniform3f(self.vectorFieldShader.uniforms.maxs.loc, self.maxs:unpack())
			-- how to determine scale?
			--local scale = app.displayVectorField_scale * (valueMax - valueMin)
			--local scale = app.displayVectorField_scale / (valueMax - valueMin)
			local scale = app.displayVectorField_scale 
				* app.displayVectorField_step 
				* math.min(
					(self.maxs[1] - self.mins[1]) / tonumber(self.gridSize.x),
					(self.maxs[2] - self.mins[2]) / tonumber(self.gridSize.y),
					(self.maxs[3] - self.mins[3]) / tonumber(self.gridSize.z))
			gl.glUniform1f(self.vectorFieldShader.uniforms.scale.loc, scale) 

			local step = app.displayVectorField_step
	
			--[[ goes just slightly faster.  24 vs 23 fps.
			app.vectorField_displayList = app.vectorField_displayList or {}
			local glCallOrDraw = require 'gl.call'
			glCallOrDraw(app.vectorField_displayList, function()
			--]]	
				gl.glBegin(gl.GL_LINES)
				for k=0,tonumber(self.sizeWithoutBorder.z-1),step do
					for j=0,tonumber(self.sizeWithoutBorder.y-1),step do
						for i=0,tonumber(self.sizeWithoutBorder.x-1),step do
							local tx = (i + .5 + self.numGhost) / tonumber(self.gridSize.x)
							local ty = (j + .5 + (self.dim > 1 and self.numGhost or 0)) / tonumber(self.gridSize.y)
							local tz = (k + .5 + (self.dim > 2 and self.numGhost or 0)) / tonumber(self.gridSize.z)
							gl.glMultiTexCoord3f(gl.GL_TEXTURE0, tx, ty, tz)	
							local x = (i + .5) / tonumber(self.sizeWithoutBorder.x)
							local y = (j + .5) / tonumber(self.sizeWithoutBorder.y)
							local z = (k + .5) / tonumber(self.sizeWithoutBorder.z)
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
			self:getTex(var):unbind(0)
			self.vectorFieldShader:useNone()
		end
	end
end

local function applyToApp(HydroCLApp)
	HydroCLApp.displayVectorField_scale = 1
	HydroCLApp.displayVectorField_step = 4
	function HydroCLApp:displayVectorField(solvers, ar, ...)
		self.view:projection(ar)
		self.view:modelview()

		gl.glDisable(gl.GL_BLEND)
		gl.glEnable(gl.GL_DEPTH_TEST)
		gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
		
		for _,solver in ipairs(solvers) do
			solver:displayVectorField(app, ...)
		end
	
		gl.glDisable(gl.GL_DEPTH_TEST)
	end
	glreport'here'
end

return {
	applyToApp = applyToApp,
	applyToSolver = applyToSolver,
}
