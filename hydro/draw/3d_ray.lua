local class = require 'ext.class'
local file = require 'ext.file'
local template = require 'template'
local gl = require 'ffi.OpenGL'
local glreport = require 'gl.report'
local Draw = require 'hydro.draw.draw'

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

local Draw3DRay = class(Draw)

Draw3DRay.useIsos = false
Draw3DRay.numIsobars = 20
Draw3DRay.useLighting = false
Draw3DRay.alpha = .15
Draw3DRay.alphaGamma = 2

function Draw3DRay:showDisplayVar(var, ar)
	local solver = self.solver
	local app = solver.app
	app.view:setup(ar)
	
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
		local tex, shader, uniforms
		if pass == 0 then
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
		else
			gl.glEnable(gl.GL_CULL_FACE)
			gl.glCullFace(gl.GL_FRONT)
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
			gl.glEnable(gl.GL_BLEND)
			
			shader = solver.volumeRayShader
			uniforms = shader.uniforms
			
			shader:use()
			
			self:setupDisplayVarShader(shader, var, valueMin, valueMax)

			gl.glUniform1f(uniforms.alpha.loc, self.alpha)
			gl.glUniform1f(uniforms.alphaGamma.loc, self.alphaGamma)
			gl.glUniform1f(uniforms.alpha.loc, self.alpha)
			gl.glUniform1i(uniforms.useIsos.loc, self.useIsos)
			gl.glUniform1f(uniforms.numIsobars.loc, self.numIsobars)				
			gl.glUniform1i(uniforms.useLighting.loc, self.useLighting)
			gl.glUniform1i(uniforms.maxiter.loc, solver.display3D_Ray_maxiter)
			local tex = solver:getTex(var)
			tex:bind(0)
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
			tex:unbind(0)
			shader:useNone()
			gl.glDisable(gl.GL_BLEND)
			gl.glDisable(gl.GL_DEPTH_TEST)
			gl.glCullFace(gl.GL_BACK)
			gl.glDisable(gl.GL_CULL_FACE)
		end
	end
end

function Draw3DRay:display(varName, ar, xmin, xmax, ymin, ymax, useLog)
	local solver = self.solver
	local app = solver.app
	local var = solver.displayVarForName[varName]
	if var and var.enabled then
		self:prepareShader()
		self:showDisplayVar(var, ar)
	end
glreport'here'
end

-- TODO this in common with 3d_iso.lua.  subclass?
function Draw3DRay:prepareShader()
	local solver = self.solver
	if solver.volumeRayShader then return end
	
	solver.display3D_Ray_maxiter = math.max(
		tonumber(solver.gridSize.x),
		tonumber(solver.gridSize.y),
		tonumber(solver.gridSize.z))
	
	local volumetricCode = assert(file['hydro/draw/volumetric.shader'])
	solver.volumeRayShader = solver.GLProgram{
		name = 'volumetric',
		vertexCode = template(volumetricCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			vertexShader = true,
		}),
		fragmentCode = template(volumetricCode, {
			draw = self,
			app = solver.app,
			solver = solver,
			fragmentShader = true,
		}),
		uniforms = {
			tex = 0,
			gradientTex = 1,
			oneOverDx = {(solver.maxs - solver.mins):unpack()},
		},
	}
end

return Draw3DRay
