#version 460

<?
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>
<?=varying?> vec3 normal;

<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>
in vec3 inVertex;

<?=solver.coord:getModuleCodeGLSL('coordMapGLSL')?>

vec3 func(vec3 src) {
	if (displayDim == 1) {
		src.yz = displayFixed;
	} else if (displayDim == 2) {
		src.z = displayFixed.y;
	}
	
	vec3 vertex = src.xyz;
	vertex.x = (vertex.x * gridSize.x - numGhost) / sizeWithoutBorder.x * (solverMaxs.x - solverMins.x) + solverMins.x;
	if (displayDim == 1) {
		vertex.y = 0.;
	} else {
		vertex.y = (vertex.y * gridSize.y - numGhost) / sizeWithoutBorder.y * (solverMaxs.y - solverMins.y) + solverMins.y;
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

	vec3 xp = func(inVertex + vec3(1./gridSize.x, 0., 0.));
	vec3 xm = func(inVertex - vec3(1./gridSize.x, 0., 0.));
	vec3 yp = func(inVertex + vec3(0., 1./gridSize.y, 0.));
	vec3 ym = func(inVertex - vec3(0., 1./gridSize.y, 0.));

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
