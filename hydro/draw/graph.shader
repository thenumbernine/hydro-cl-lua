#version 460

<?
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>
<?=varying?> vec3 normal;

uniform float ambient;

<?=draw:getModuleCodeGLSL("coordMapGLSL", "coordMapInvGLSL")?>
<?=draw:getCommonGLSLFragCode()?>

<? if vertexShader then ?>
in vec3 gridCoord;

vec3 getVertex(vec3 gridCoord) {
	//convert from integer gridCoord to texNoGhostCoord
	vec3 texCoord = (gridCoord + .5 + numGhost) / gridSize;
	vec3 chartCoord = texToChartCoord(texCoord);
	vec3 vertex = chartCoord;
	if (displayDim <= 1) chartCoord.y = displayFixed.x;
	if (displayDim <= 2) chartCoord.z = displayFixed.y;
	
	//should this go before vertex def, and allow rotations of the vertices themselves?
	//or should this go here, and only rotate the data source?
	chartCoord = chartToWorldCoord(chartCoord);
	chartCoord = quatRotate(displaySliceAngle, chartCoord);
	chartCoord = worldToChartCoord(chartCoord);
	texCoord = chartToTexCoord(chartCoord);

	float value = getTex(texCoord).r;
	// for displayDim==2, use the normalized height
	if (displayDim > 1) {
		value = getGradientFrac(value);
	} else {
		// for displayDim==1, the GL ymin/ymax is set to the solver ymin/ymax, so don't transform the value
		if (useLog) {
			value = logmap(value);
		}
		//also make sure the vertex is in the z=0 plane
		vertex.z = 0.;
	}
	vertex[displayDim] = value;
	return vertex;
}

void main() {
	vec3 vertex = getVertex(gridCoord);

	if (ambient < 1.) {
		vec3 xp = getVertex(gridCoord + vec3(1., 0., 0.));
		vec3 xm = getVertex(gridCoord - vec3(1., 0., 0.));
		vec3 yp = getVertex(gridCoord + vec3(0., 1., 0.));
		vec3 ym = getVertex(gridCoord - vec3(0., 1., 0.));
		normal = normalize(cross(xp - xm, yp - ym));
	} else {
		normal = vec3(0., 0., 0.);
	}

	if (displayDim > 1) {
		vertex = chartToWorldCoord(vertex);
	}
	gl_Position = modelViewProjectionMatrix * vec4(vertex, 1.);
}

<?
end
if fragmentShader then ?>

out vec4 fragColor;

uniform vec3 color;

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
