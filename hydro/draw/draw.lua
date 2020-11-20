--[[
parent class for some common functions that all draw glsl code uses
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local gl = require 'gl'
local matrix_ffi = require 'matrix.ffi'

matrix_ffi.real = 'float'	-- default matrix_ffi type

local Draw = class()

function Draw:init(solver)
	self.solver = assert(solver)
end

function Draw:getCommonGLSLFragCode()
	local solver = self.solver
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

uniform vec3 texSize;
uniform vec3 gridSize;
uniform vec3 sizeWithoutBorder;

uniform vec4 displaySliceAngle;

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
	//const float epsilon = 1.;
	//const float epsilon = 1e-7;
	const float epsilon = 1e-14;
	//const float epsilon = 0.;
	return log(epsilon + abs(x)) * _1_LN_10;
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


// here's my naming convention
// grid coords: i in [.5, n-.5] in Integers-1/2
// tex coords: t in [.5, n-.5]/n in Reals. t = i / sizeWithoutBorder, i in [.5, n-.5]
// tex no-ghost coords: tng in [.5-g,n-.5-g]/(n-2g) in Reals. t = i / sizeWithoutBorder, i in [.5+g, n-.5-g]
// chart coords: x in [solver.mins, solver.maxs] in Reals.  x = tng * (solverMax - solverMin) + solverMin
// world coords: y in Reals.  nonlinear transform.  y = f(x)
// 'useCoordMap' replaces the chart f(x) with f(x) = (x - solverMins)/(solverMaxs - solverMins) * 2 - 1 = t * 2 - 1

//chart coords to world coords
//gets coordMap, or the normalized version, depending on the 'useCoordMap' flag
vec3 chartToWorldCoord(vec3 x) {
	if (useCoordMap) {
		return coordMap(x);
	} else {
		return (x - solverMins) / (solverMaxs - solverMins) * 2. - 1.;
	}
}

vec3 texToChartCoord(vec3 x) {
	return x * (solverMaxs - solverMins) + solverMins;
}

//should match chartToWorldCoord(texToChartCoord(x))
vec3 texToWorldCoord(vec3 x) {
	if (useCoordMap) {
		return coordMap(x * (solverMaxs - solverMins) + solverMins);
	} else {
		return x * 2. - 1.;
	}
}

vec3 worldToTexCoord(vec3 x) {
	if (useCoordMap) {
		return (coordMapInv(x) - solverMins) / (solverMaxs - solverMins);
	} else {
		return .5 * (x + 1.);
	}
}

//world coords to chart coords
vec3 worldToChartCoord(vec3 x) {
	if (useCoordMap) {
		return coordMapInv(x);
	} else {
		return .5 * (x + 1.) * (solverMaxs - solverMins) + solverMins; 
	}
}

//chart coord to tex coord
vec3 chartToTexCoord(vec3 x) {
	return ((x - solverMins) / (solverMaxs - solverMins) * sizeWithoutBorder + numGhost) / gridSize;
}

vec3 texToNoGhostCoord(vec3 x) {
	return (x * sizeWithoutBorder + numGhost) / gridSize;
}

vec3 chartToNoGhostCoord(vec3 x) {
	return texToNoGhostCoord(chartToTexCoord(x));
}

// hmm, where to put this?
vec3 quatRotate(vec4 q, vec3 v) { 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

vec4 quatConj(vec4 q) {
	return vec4(q.xyz, -q.w);
}

]], {
		solver = solver,
		app = solver.app,
		clnumber = require 'cl.obj.number',
	})
end

function Draw:setupDisplayVarShader(shader, var, valueMin, valueMax)
	local solver = self.solver
	local app = solver.app

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
	if uniforms.texSize then
		gl.glUniform3f(uniforms.texSize.loc, solver.texSize:unpack())
	end
	if uniforms.gridSize then
		gl.glUniform3f(uniforms.gridSize.loc, solver.gridSize:unpack())
	end
	if uniforms.sizeWithoutBorder then
		gl.glUniform3f(uniforms.sizeWithoutBorder.loc, solver.sizeWithoutBorder:unpack())
	end
	if uniforms.displaySliceAngle then
		--gl.glUniform4fv(uniforms.displaySliceAngle.loc, 4, app.displaySliceAngle.s)
		gl.glUniform4f(uniforms.displaySliceAngle.loc, app.displaySliceAngle:unpack())
	end
end

local function makeGLSL(code)
	return table{
		'#define M_PI '..('%.50f'):format(math.pi),
		'#define inline',
		'#define real					float',
		'#define _real3					vec3',
		'#define real3					vec3',
		'#define real3_zero				vec3(0.,0.,0.)',
		'#define real3_add(a,b)			((a)+(b))',
		'#define real3_sub(a,b)			((a)-(b))',
		'#define real3_real_mul(a,b)	((a)*(b))',
		'#define real3_dot				dot',
		'#define real3_len(a)			length(a)',
		'#define real3_lenSq(a)			dot(a,a)',
		'#define atan2					atan',
		'#define fmod					mod',
		'#define fabs					abs',
		'',
		(code:gsub('static inline ', '')),
	}:concat'\n'
end

function Draw:getModuleCodeGLSL(...)
	return makeGLSL(self.solver.modules:getCodeAndHeader(...))
end

return Draw
