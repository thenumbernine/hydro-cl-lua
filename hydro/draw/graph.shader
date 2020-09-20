#version 460

<?
local inout = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set inout to")
?>
<?=inout?> vec3 normal;

<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>
in vec3 inVertex;

uniform vec2 size;

<?
local clnumber = require 'cl.obj.number'
local sizeWithoutBorderX = clnumber(solver.sizeWithoutBorder.x or 1)
local sizeWithoutBorderY = clnumber(solver.sizeWithoutBorder.y or 1)
local sizeX = clnumber(solver.gridSize.x or 1)
local sizeY = clnumber(solver.gridSize.y or 1)
?>

<?=solver.coord:getModuleCodeGLSL('coordMapGLSL')?>

vec3 func(vec3 src) {
	vec3 vertex = src.xyz;
	vertex.x = (vertex.x * <?=sizeX?> - numGhost) / <?=sizeWithoutBorderX?> * (solverMaxs.x - solverMins.x) + solverMins.x;
	if (displayDim == 1) {
		vertex.y = 0.;
	} else {
		vertex.y = (vertex.y * <?=sizeY?> - numGhost) / <?=sizeWithoutBorderY?> * (solverMaxs.y - solverMins.y) + solverMins.y;
	}

	vertex[displayDim] = getTex(src).r;
	if (displayDim > 1) {
		vertex[displayDim] -= valueMin;
		vertex[displayDim] *= 1. / (valueMax - valueMin);
	}
	if (useLog) {
		vertex[displayDim] = log(max(0., vertex[displayDim])) * <?=('%.50f'):format(1/math.log(10))?>;
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
	gl_Position = modelViewProjectionMatrix * vec4(vertex, 1.);
}

<?
end
if fragmentShader then ?>

out vec4 fragColor;

uniform vec3 color;
uniform float ambient;

void main() {
	vec3 light = normalize(vec3(.5, .5, 1.));
	float lum = dot(normal, light);
	//lum = max(lum, -lum);	//...for two-sided
	lum = max(lum, ambient);
	fragColor.rgb = color;
	fragColor.rgb *= lum;
	fragColor.a = 1.;
}

<? end ?>
