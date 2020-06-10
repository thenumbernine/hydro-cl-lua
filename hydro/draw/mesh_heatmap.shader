<? if vertexShader then ?>
#version 460

uniform float drawCellScale;
uniform mat4 modelViewProjectionMatrix;

attribute vec3 vtx;
attribute vec3 vtxcenter;
attribute float value;

varying float valuev;

void main() {
	vec3 v = (vtx - vtxcenter) * drawCellScale + vtxcenter;
	gl_Position = modelViewProjectionMatrix * vec4(v, 1.);
	valuev = value;
}

<? end
if fragmentShader then ?>
#version 460

<?
local clnumber = require 'cl.obj.number'
?>

//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useLog;
uniform float valueMin;
uniform float valueMax;

uniform sampler1D gradientTex;

varying float valuev;

out vec4 fragColor;

void main() {
	float value = valuev;	
	if (useLog) {
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
	value = (value * <?=clnumber(app.gradientTex.width-1)?> + .5) / <?=clnumber(app.gradientTex.width)?>;
	fragColor = texture1D(gradientTex, value);
}
<? end ?>
