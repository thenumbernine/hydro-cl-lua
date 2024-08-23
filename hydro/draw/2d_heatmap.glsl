<?=draw.glslVersion?>

<?
local varying = vertexShader and 'out'
	or fragmentShader and 'in'
	or error("don't know what to set varying to")
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
<?=varying?> vec3 viewCoord;

<?=draw:getModuleCodeGLSL(coordMapGLSL, coordMapInvGLSL)?>
<?=draw:getCommonGLSLFragCode()?>

<? if vertexShader then ?>

in vec2 vertex;
uniform vec4 bbox;	//[x1,y1,x2,y2]
void main() {
	vec2 rvtx = mix(bbox.xy, bbox.zw, vertex);
	viewCoord = vec3(rvtx, 0.);
	gl_Position = mvProjMat * vec4(rvtx, 0., 1.);
}

<? end
if fragmentShader then ?>

out vec4 fragColor;

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	//TODO make this a flag:
	vec3 chartCoord = viewCoord;
	chartCoord.z = displayFixed.y;
	chartCoord = quatRotate(displaySliceAngle, chartCoord);
	chartCoord = worldToChartCoord(chartCoord);

#if 0	//add border
const float epsilon = 1e-2;
if (abs(chartCoord.x - solverMins.x) < epsilon ||
	abs(chartCoord.x - solverMaxs.x) < epsilon ||
	abs(chartCoord.y - solverMins.y) < epsilon ||
	abs(chartCoord.y - solverMaxs.y) < epsilon)
{
	fragColor = vec4(1.,1.,1.,1.);
	return;
}
#endif

	vec3 texCoord = chartToTexCoord(chartCoord);
#if 1
	if (texCoord.x < 0. || texCoord.x > 1.) discard;
	if (displayDim > 1 && (texCoord.y < 0. || texCoord.y > 1.)) discard;
	if (displayDim > 2 && (texCoord.z < 0. || texCoord.z > 1.)) discard;
#endif

	float value = getTex(texCoord).r;
	fragColor = getGradientColor(value);
}

<? end ?>
