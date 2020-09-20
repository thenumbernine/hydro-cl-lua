--[[
parent class for some common functions that all draw glsl code uses
--]]

local class = require 'ext.class'
local template = require 'template'
local gl = require 'gl'
local matrix_ffi = require 'matrix.ffi'

matrix_ffi.real = 'float'	-- default matrix_ffi type

local Draw = class()

function Draw:getCommonGLSLFragCode(solver)
	return template([[

uniform int displayDim;
uniform vec2 displayFixed;	//xy holds the fixed yz for when displayDim < dim

<? 
-- if solver.dim < 3 then -- doesn't consider meshsolver
if require 'gl.tex2d'.is(solver.tex) then -- does
?>
uniform sampler2D tex;
vec4 getTex(vec3 texCoord) {
	return texture2D(tex, texCoord.xy);
}
<? else ?>
uniform sampler3D tex;
vec4 getTex(vec3 texCoord) {
	return texture3D(tex, texCoord);
}
<? end ?>

]], {
		solver = solver,
	})
end

function Draw:setupDisplayVarShader(shader, app, solver, var)
	local uniforms = shader.uniforms
	if uniforms.displayDim then
		gl.glUniform1i(uniforms.displayDim.loc, app.displayDim)
	end
	if uniforms.displayFixed then	
		gl.glUniform2f(uniforms.displayFixed.loc, app.displayFixedY, app.displayFixedZ)
	end
	if uniforms.modelViewProjectionMatrix then
		gl.glUniformMatrix4fv(uniforms.modelViewProjectionMatrix.loc, 1, 0, app.view.modelViewProjectionMatrix.ptr)
	end
	if uniforms.normalMatrix then
		self.normalMatrix = self.normalMatrix or matrix_ffi.zeros(3,3)
		--gl.glGetFloatv(gl.GL_NORMAL_MATRIX, self.normalMatrix.ptr)
		-- invert app.view.modelViewMatrix's upper 3x3 into normalMatrix, then transpose
		-- but if it is purely a rotation matrix, this is the same as just the 3x3 portion ...
		for j=0,2 do
			for i=0,2 do
				self.normalMatrix.ptr[i + 3 * j] = app.view.modelViewProjectionMatrix.ptr[i + 4 * j]
			end
		end
		gl.glUniformMatrix3fv(uniforms.normalMatrix.loc, 1, 0, self.normalMatrix.ptr)
	end
end

return Draw
