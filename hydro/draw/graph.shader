//my intel ubuntu opengl driver's error:
// error: GLSL 4.60 is not supported. Supported versions are: 1.10, 1.20, 1.30, 1.00 ES, 3.00 ES, 3.10 ES, and 3.20 ES
//#version 460
//...sooo...
#version 130

<? if vertexShader then ?>
in vec3 inVertex;
out vec3 normal;

uniform mat4 ModelViewProjectionMatrix;

<? if solver.dim == 3 then
?>uniform sampler3D tex;
<? else
?>uniform sampler2D tex;
<? end
?>uniform int axis;
uniform vec2 xmin;
uniform vec2 xmax;
uniform float scale;
uniform float offset;
uniform bool useLog;
uniform vec2 size;

<?
local clnumber = require 'cl.obj.number'
local sizeWithoutBorderX = clnumber(solver.sizeWithoutBorder.x or 1)
local sizeWithoutBorderY = clnumber(solver.sizeWithoutBorder.y or 1)
local sizeX = clnumber(solver.gridSize.x or 1)
local sizeY = clnumber(solver.gridSize.y or 1)
local numGhost = clnumber(solver.numGhost)
?>

<?=solver.coord:getModuleCodeGLSL('coordMapGLSL')?>

vec3 func(vec3 src) {
	vec3 vertex = src.xyz;
	vertex.x = (vertex.x * <?=sizeX?> - <?=numGhost?>) / <?=sizeWithoutBorderX?> * (xmax.x - xmin.x) + xmin.x;
	vertex.y = (vertex.y * <?=sizeY?> - <?=numGhost?>) / <?=sizeWithoutBorderY?> * (xmax.y - xmin.y) + xmin.y;

<? if solver.dim == 3 then
?>	vertex[axis] = texture3D(tex, src).r;
<? else
?>	vertex[axis] = texture2D(tex, src.xy).r;
<? end
?>	vertex[axis] -= offset;
	vertex[axis] *= scale;
	if (useLog) {
		vertex[axis] = log(max(0., vertex[axis])) * <?=('%.50f'):format(1/math.log(10))?>;
	}
	return vertex;
}

void main() {
	vec3 vertex = func(inVertex);

	vec3 xp = func(inVertex + vec3(1./size.x, 0., 0.));
	vec3 xm = func(inVertex - vec3(1./size.x, 0., 0.));
	vec3 yp = func(inVertex + vec3(0., 1./size.y, 0.));
	vec3 ym = func(inVertex - vec3(0., 1./size.y, 0.));

	normal = normalize(cross(xp - xm, yp - ym));

	vertex = coordMap(vertex);
	gl_Position = ModelViewProjectionMatrix * vec4(vertex, 1.);
}

<?
end
if fragmentShader then ?>

in vec3 normal;

uniform vec3 color;
uniform float ambient;

void main() {
	vec3 light = normalize(vec3(.5, .5, 1.));
	float lum = dot(normal, light);
	//lum = max(lum, -lum);	//...for two-sided
	lum = max(lum, ambient);
	gl_FragColor.rgb = color;
	gl_FragColor.rgb *= lum;
	gl_FragColor.a = 1.;
}

<? end ?>
