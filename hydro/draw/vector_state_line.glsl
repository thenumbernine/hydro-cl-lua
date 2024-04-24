#version 460

<?
local varying = vertexShader and "out"
		or fragmentShader and "in"
		or error("don't know what to set varying to")
?>

<?=draw:getModuleCodeGLSL(coordMapGLSL, coordMapInvGLSL)?>
<?=draw:getCommonGLSLFragCode()?>

<? if vertexShader then ?>

in vec3 gridCoord;

<?=varying?> vec3 texCoord;
<?=varying?> vec3 normal;
uniform float ambient;

vec3 getVertex(vec3 gridCoord) {
	//convert from integer gridCoord to texNoGhostCoord
	vec3 texCoord = (gridCoord + .5 + numGhost) / sizeWithoutBorder;
	vec3 chartCoord = texToChartCoord(texCoord);
	if (displayDim <= 1) chartCoord.y = displayFixed.x;
	if (displayDim <= 2) chartCoord.z = displayFixed.y;
	
	//should this go before vertex def, and allow rotations of the vertices themselves?
	//or should this go here, and only rotate the data source?
	chartCoord = chartToWorldCoord(chartCoord);
	chartCoord = quatRotate(displaySliceAngle, chartCoord);
	chartCoord = worldToChartCoord(chartCoord);
	texCoord = chartToTexCoord(chartCoord);

	vec3 vertex = getTex(texCoord).rgb;
	vertex = chartToWorldCoord(vertex);

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

	texCoord = (gridCoord + .5 + numGhost) / sizeWithoutBorder;

	gl_Position = mvProjMat * vec4(vertex, 1.);
}

<?
end
if fragmentShader then ?>

<?=varying?> vec3 texCoord;
<?=varying?> vec3 normal;
uniform float ambient;
out vec4 fragColor;

void main() {
	//fragColor = texture(gradientTex, texCoord.r);
	//how about 2D / 3D datasets?  need a better texture map ...
	fragColor = vec4(texCoord * .5 + .5, 1.);
	
	vec3 light = normalize(vec3(.5, .5, 1.));
	float lum = dot(normal, light);
	lum = max(lum, -lum);	//...for two-sided
	lum = max(lum, ambient);
	fragColor.rgb *= lum;
}

<? end ?>
