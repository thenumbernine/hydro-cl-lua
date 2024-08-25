#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>


uniform int displayDim;
uniform vec2 displayFixed;	//xy holds the fixed yz for when displayDim < dim

//#define instead of uniform since numGhost isn't going to change at runtime
#define numGhost <?=clnumber(solver.numGhost)?>

uniform bool useCoordMap;

uniform mat3 normalMatrix;
uniform mat4 mvProjMat;

uniform bool useLog;
uniform float valueMin, valueMax;

uniform vec3 solverMins, solverMaxs;

uniform vec3 texSize;
uniform vec3 gridSize;
uniform vec3 sizeWithoutBorder;

uniform vec4 displaySliceAngle;

uniform vec3 cartesianMin, cartesianMax;

<?
-- if solver.dim < 3 then -- doesn't consider meshsolver
if require 'gl.tex2d':isa(solver.tex) then -- does
?>
layout(binding=0) uniform sampler2D tex;
vec4 getTex(vec3 texCoord) {
	return texture(tex, texCoord.xy);
}
<? else ?>
layout(binding=0) uniform sampler3D tex;
vec4 getTex(vec3 texCoord) {
	return texture(tex, texCoord);
}
<? end ?>

float logmap(float x) {
	//const float epsilon = 1.;
	//const float epsilon = 1e-7;
	//const float epsilon = 1e-14;
	const float epsilon = 0.;
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

layout(binding=1) uniform sampler2D gradientTex;

//TODO just change texel lookup in gradTex?
float getGradientTexCoord(float frac) {
	return (frac * <?=clnumber(app.gradientTex.width-1)?> + .5) * <?=clnumber(1 / app.gradientTex.width)?>;
}

vec4 getGradientColor(float value) {
	float frac = getGradientFrac(value);
	float tc = getGradientTexCoord(frac);
	return texture(gradientTex, vec2(tc, .5));
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
