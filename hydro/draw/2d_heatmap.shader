#version 460

<?
local coord = solver.coord
local varying = vertexShader and 'out'
		or fragmentShader and 'in'
		or error("don't know what to set varying to")
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
<?=varying?> vec3 viewCoord;

<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>

attribute vec4 vertex;

void main() {
	vec4 v;
	v.xy = vertex.xy;
	//v.z = 0.;// should we disregard 'z' for rendering?
	viewCoord = v.xyz;
	v.w = 1.;
	gl_Position = modelViewProjectionMatrix * v;
}

<? end
if fragmentShader then ?>

//this is from the start of the amr stuff that I haven't finished
uniform vec2 texCoordMax;

out vec4 fragColor;


<?=solver.coord:getModuleCodeGLSL('coordMapInvGLSL')?>

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	//TODO make this a flag:
	vec3 gridCoord = viewCoord;
	gridCoord.z = displayFixed.y;
	gridCoord = quatRotate(displaySliceAngle, gridCoord);
	if (useCoordMap) {
		gridCoord = coordMapInv(gridCoord);
	} else {
		gridCoord = .5 * (gridCoord + 1.) * (solverMaxs - solverMins) + solverMins; 
	}

#if 0	//add border
const float epsilon = 1e-2;
if (abs(gridCoord.x - solverMins.x) < epsilon ||
	abs(gridCoord.x - solverMaxs.x) < epsilon ||
	abs(gridCoord.y - solverMins.y) < epsilon ||
	abs(gridCoord.y - solverMaxs.y) < epsilon)
{
	fragColor = vec4(1.,1.,1.,1.);
	return;
}
#endif

	vec3 texCoord = vec3(
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * sizeWithoutBorder.x + numGhost) / gridSize.x,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * sizeWithoutBorder.y + numGhost) / gridSize.y,
		((gridCoord.z - solverMins.z) / (solverMaxs.z - solverMins.z) * sizeWithoutBorder.z + numGhost) / gridSize.z
	) * vec3(texCoordMax, 1.);

#if 1
	if (texCoord.x < 0. || texCoord.x > 1.) discard;
	if (displayDim > 1 && (texCoord.y < 0. || texCoord.y > 1.)) discard;
	if (displayDim > 2 && (texCoord.z < 0. || texCoord.z > 1.)) discard;
#endif

	float value = getTex(texCoord).r;
	fragColor = getGradientColor(value);
}

<? end ?>
