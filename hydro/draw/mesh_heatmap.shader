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

<?=solver:getGradientGLSLCode()?>

varying float valuev;
out vec4 fragColor;

void main() {
	fragColor = getGradientColor(valuev);
}
<? end ?>
