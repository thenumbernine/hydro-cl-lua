local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'

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

return function(HydroCLApp)

function HydroCLApp:display3D_Ray(solvers, varName, ar, xmin, ymin, xmax, ymax, useLog)
	self.view:projection(ar)
	self.view:modelview()

	for _,solver in ipairs(solvers) do
		local var = solver.displayVarForName[varName]
		if var and var.enabled then
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
			for pass=0,1 do
				if pass == 0 then
					gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
				else
					gl.glEnable(gl.GL_CULL_FACE)
					gl.glCullFace(gl.GL_FRONT)
					gl.glEnable(gl.GL_DEPTH_TEST)
					gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
					gl.glEnable(gl.GL_BLEND)
					solver.volumeRayShader:use()
					gl.glUniform1f(solver.volumeRayShader.uniforms.scale.loc, 1)--scale)
					gl.glUniform1i(solver.volumeRayShader.uniforms.useLog.loc, 0)--useLog and 1 or 0)
					gl.glUniform1f(solver.volumeRayShader.uniforms.alpha.loc, 1)--alpha)
					solver:getTex(var):bind(0)
					self.gradientTex:bind(1)
				end
				gl.glBegin(gl.GL_QUADS)
				for i=1,24 do
					local x = vertexesInCube[quadsInCube[i] * 3 + 0 + 1]
					local y = vertexesInCube[quadsInCube[i] * 3 + 1 + 1]
					local z = vertexesInCube[quadsInCube[i] * 3 + 2 + 1]
					gl.glTexCoord3f(x, y, z)
					x = x * (solver.maxs[1] - solver.mins[1]) + solver.mins[1]
					y = y * (solver.maxs[2] - solver.mins[2]) + solver.mins[2]
					z = z * (solver.maxs[3] - solver.mins[3]) + solver.mins[3]
					gl.glVertex3f(x, y, z)
				end
				gl.glEnd()
				if pass == 0 then
					gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
				else
					self.gradientTex:unbind(1)
					solver:getTex(var):unbind(0)
					solver.volumeRayShader:useNone()
					gl.glDisable(gl.GL_BLEND)
					gl.glDisable(gl.GL_DEPTH_TEST)
					gl.glCullFace(gl.GL_BACK)
					gl.glDisable(gl.GL_CULL_FACE)
				end
			end
		end
	end
glreport'here'
end

end
