#version 460

<?
local clnumber = require 'cl.obj.number'
local coord = solver.coord
?>

<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>

varying vec4 color;

<?=coord:getModuleCodeGLSL('coordMapGLSL', 'cartesianFromCoord')?>

attribute vec2 vtx;
attribute vec3 center;
attribute vec3 tc;

uniform float scale;

uniform vec3 solverMins, solverMaxs;

void main() {
	vec3 realtc = tc;
	if (displayDim <= 1) realtc.y = displayFixed.x;
	if (displayDim <= 2) realtc.z = displayFixed.y;
	vec3 dir = getTex(realtc).rgb;

	dir = cartesianFromCoord(dir, center);	
	float value = length(dir);
	dir /= value;

	vec3 tv;
	if (displayDim < 3) {
		tv = vec3(-dir.y, dir.x, 0.);
	} else {	//displayDim == 3
		vec3 vx = vec3(0., -dir.z, dir.y);
		vec3 vy = vec3(dir.z, 0., -dir.x);
		vec3 vz = vec3(-dir.y, dir.x, 0.);
		float lxsq = dot(vx,vx);
		float lysq = dot(vy,vy);
		float lzsq = dot(vz,vz);
		if (lxsq > lysq) {		//x > y
			if (lxsq > lzsq) {	//x > z, x > y
				tv = vx;
			} else {			//z > x > y
				tv = vz;
			}
		} else {				//y >= x
			if (lysq > lzsq) {	//y >= x, y > z
				tv = vy;
			} else {			// z > y >= x
				tv = vz;
			}
		}
	}

	//TODO don't just normalize dir, but scale down by the view size ... ? or the dx?
	
	value = getGradientFrac(value);
//TODO make this a flag, whether to normalize vectors or not? 
	float valuescale = scale;// * clamp(value, 0., 1.);
	value = getGradientTexCoord(value);
	color = texture1D(gradientTex, value);

	//cartesian coords
	vec3 centerpos = center;
	if (useCoordMap) {
		centerpos = coordMap(centerpos);
	} else {
		centerpos = (centerpos - solverMins) / (solverMaxs - solverMins) * vec3(2., 2., 1.) - vec3(1., 1., 0.);
	}
	vec3 v = centerpos + valuescale * (vtx.x * dir + vtx.y * tv);
	
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
}

<?
end
if fragmentShader then
?>

varying vec4 color;

out vec4 fragColor;

void main() {
	fragColor = color;
}

<? end ?>
