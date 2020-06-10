local gl = require 'ffi.OpenGL'
local class = require 'ext.class'
local glreport = require 'gl.report'

-- TODO real raytracing:  
-- instead of drawing the cube, draw a quad over the whole screen.

-- 3D
local vertexesInCube = {
	0,0,0,
	1,0,0,
	0,1,0,
	1,1,0,
	0,0,1,
	1,0,1,
	0,1,1,
	1,1,1,
}

local quadsInCube = {
	0,1,3,2,
	4,6,7,5,
	1,5,7,3,
	0,2,6,4,
	0,4,5,1,
	2,3,7,6,
}

local Draw3DRay = class()

Draw3DRay.useIsos = false
Draw3DRay.numIsobars = 20
Draw3DRay.useLighting = false
Draw3DRay.alpha = .15
Draw3DRay.alphaGamma = 2

function Draw3DRay:showDisplayVar(app, solver, var, ar)
	app.view:projection(ar)
	app.view:modelview()
	
	local valueMin, valueMax
	if var.heatMapFixedRange then
		valueMin = var.heatMapValueMin
		valueMax = var.heatMapValueMax
	else
		valueMin, valueMax = solver:calcDisplayVarRange(var)
		var.heatMapValueMin = valueMin
		var.heatMapValueMax = valueMax
	end
	
	solver:calcDisplayVarToTex(var)	

	gl.glColor3f(1,1,1)
	for pass=1,1 do
		if pass == 0 then
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
		else
			gl.glEnable(gl.GL_CULL_FACE)
			gl.glCullFace(gl.GL_FRONT)
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
			gl.glEnable(gl.GL_BLEND)
			solver.volumeRayShader:use()
			gl.glUniform1f(solver.volumeRayShader.uniforms.alpha.loc, self.alpha)
			gl.glUniform1f(solver.volumeRayShader.uniforms.alphaGamma.loc, self.alphaGamma)
			gl.glUniform3f(solver.volumeRayShader.uniforms.mins.loc, solver.mins:unpack())
			gl.glUniform3f(solver.volumeRayShader.uniforms.maxs.loc, solver.maxs:unpack())
			gl.glUniform1i(solver.volumeRayShader.uniforms.useLog.loc, var.useLog and 1 or 0)
			gl.glUniform1f(solver.volumeRayShader.uniforms.alpha.loc, self.alpha)
			gl.glUniform1f(solver.volumeRayShader.uniforms.valueMin.loc, valueMin)
			gl.glUniform1f(solver.volumeRayShader.uniforms.valueMax.loc, valueMax)
			gl.glUniform1i(solver.volumeRayShader.uniforms.useIsos.loc, self.useIsos)
			gl.glUniform1f(solver.volumeRayShader.uniforms.numIsobars.loc, self.numIsobars)				
			gl.glUniform1i(solver.volumeRayShader.uniforms.useLighting.loc, self.useLighting)
			gl.glUniform1f(solver.volumeRayShader.uniforms.numGhost.loc, solver.numGhost)
			gl.glUniform3f(solver.volumeRayShader.uniforms.texSize.loc, solver.gridSize:unpack())
			gl.glUniform1i(solver.volumeRayShader.uniforms.maxiter.loc, solver.display3D_Ray_maxiter)
			solver:getTex(var):bind(0)
			app.gradientTex:bind(1)
		end
		gl.glBegin(gl.GL_QUADS)
		for i=1,24 do
			local x = vertexesInCube[quadsInCube[i] * 3 + 0 + 1]
			local y = vertexesInCube[quadsInCube[i] * 3 + 1 + 1]
			local z = vertexesInCube[quadsInCube[i] * 3 + 2 + 1]
			gl.glTexCoord3f(x, y, z)
			gl.glVertex3f(x, y, z)
		end
		gl.glEnd()
		if pass == 0 then
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
		else
			app.gradientTex:unbind(1)
			solver:getTex(var):unbind(0)
			solver.volumeRayShader:useNone()
			gl.glDisable(gl.GL_BLEND)
			gl.glDisable(gl.GL_DEPTH_TEST)
			gl.glCullFace(gl.GL_BACK)
			gl.glDisable(gl.GL_CULL_FACE)
		end
	end
end

function Draw3DRay:display(app, solvers, varName, ar, xmin, xmax, ymin, ymax, useLog)
	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
			self:showDisplayVar(app, solver, var, ar)
		end
	end
glreport'here'
end

return function(HydroCLApp)
	function HydroCLApp:display3D_Ray(...)
		if not self.draw3DRay then self.draw3DRay = Draw3DRay() end
		return self.draw3DRay:display(self, ...)
	end
end
