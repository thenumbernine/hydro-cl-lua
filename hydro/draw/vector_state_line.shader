#version 460

<?
local varying = vertexShader and "out"
		or fragmentShader and "in"
		or error("don't know what to set varying to")
?>

<?=draw:getModuleCodeGLSL(coordMapGLSL, coordMapInvGLSL)?>
<?=draw:getCommonGLSLFragCode()?>

<? if vertexShader then ?>

<?=varying?> vec3 texCoord;

in vec3 gridCoord;

void main() {
	
	//convert from integer gridCoord to texNoGhostCoord
	texCoord = (gridCoord + .5 + numGhost) / sizeWithoutBorder;
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

	gl_Position = modelViewProjectionMatrix * vec4(vertex, 1.);
}

<?
end
if fragmentShader then ?>

<?=varying?> vec3 texCoord;

out vec4 fragColor;

void main() {
	fragColor = texture1D(gradientTex, texCoord.r);
}

<? end ?>
