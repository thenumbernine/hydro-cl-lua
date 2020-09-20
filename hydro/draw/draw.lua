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
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>



uniform int displayDim;
uniform vec2 displayFixed;	//xy holds the fixed yz for when displayDim < dim

uniform float numGhost;

uniform bool useCoordMap;

uniform mat3 normalMatrix;
uniform mat4 modelViewProjectionMatrix;

uniform bool useLog;
uniform float valueMin;
uniform float valueMax;

uniform sampler1D gradientTex;

uniform vec3 solverMins;
uniform vec3 solverMaxs;

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


float logmap(float x) {
	return log(1. + abs(x)) * _1_LN_10;
}

//used by 3d_slice, vector_arrow, volumetric
float getGradientFrac(float value) {
	if (useLog) {
		
		// TODO all of this on CPU before setting the uniforms
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		if ((valueMin < 0.) != (valueMax < 0.)) {
			logValueMin = 0.;	//logmap(0)
		} else {
			if (valueMin < 0.) {
				float tmp = logValueMin;
				logValueMin = logValueMax;
				logValueMax = tmp;
			}
		}
		
		return (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		return (value - valueMin) / (valueMax - valueMin);
	}
}

//TODO just change texel lookup in gradTex?
float getGradientTexCoord(float frac) {
	return (frac * <?=clnumber(app.gradientTex.width-1)?> + .5) * <?=clnumber(1 / app.gradientTex.width)?>;
}

vec4 getGradientColor(float value) {
	float frac = getGradientFrac(value);
	float tc = getGradientTexCoord(frac);
	return texture1D(gradientTex, tc);
}

]], {
		solver = solver,
		app = solver.app,
		clnumber = require 'cl.obj.number',
	})
end

function Draw:setupDisplayVarShader(shader, app, solver, var, valueMin, valueMax)
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
	if uniforms.useCoordMap then
		gl.glUniform1i(uniforms.useCoordMap.loc, app.display_useCoordMap and 1 or 0)
	end
	if uniforms.useLog then
		gl.glUniform1i(uniforms.useLog.loc, var.useLog and 1 or 0)
	end
	if uniforms.valueMin then
		gl.glUniform1f(uniforms.valueMin.loc, valueMin)
	end
	if uniforms.valueMax then
		gl.glUniform1f(uniforms.valueMax.loc, valueMax)
	end
	if uniforms.solverMins then
		gl.glUniform3f(uniforms.solverMins.loc, solver.mins:unpack())
	end
	if uniforms.solverMaxs then	
		gl.glUniform3f(uniforms.solverMaxs.loc, solver.maxs:unpack())
	end
	if uniforms.numGhost then
		gl.glUniform1f(uniforms.numGhost.loc, solver.numGhost)
	end
end

return Draw
