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
	v.z = displayFixed.y;
	viewCoord = v.xyz;
	//v.z = 0.;// should we disregard 'z' for rendering?
	v.w = 1.;
	gl_Position = modelViewProjectionMatrix * v;
}

<? end
if fragmentShader then ?>

<?=solver.coord:getModuleCodeGLSL('coordMapInvGLSL')?>

//this is from the start of the amr stuff that I haven't finished
uniform vec2 texCoordMax;

out vec4 fragColor;

void main() {
	//start in 2D coords, bounded by the screen space
	//TODO if we are viewing this in 3D then we will have to draw a quad bigger than the intersection of the camera hull with the XY plane
	//TODO make this a flag:
	vec2 gridCoord = viewCoord.xy;
	if (useCoordMap) {
		gridCoord = coordMapInv(vec3(gridCoord, 0.)).xy;
	} else {
		gridCoord = .5 * (gridCoord + 1.) * (solverMaxs.xy - solverMins.xy) + solverMins.xy; 
	}
	
	if (gridCoord.x < solverMins.x || gridCoord.x > solverMaxs.x ||
		gridCoord.y < solverMins.y || gridCoord.y > solverMaxs.y
	) {
		discard;
	}

#if 0
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

	vec3 texCoord;
	texCoord.xy = vec2(
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * sizeWithoutBorder.x + numGhost) / gridSize.x,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * sizeWithoutBorder.y + numGhost) / gridSize.y
	) * texCoordMax;
	texCoord.z = viewCoord.z;

	float value = getTex(texCoord).r;
	fragColor = getGradientColor(value);
}

<? end ?>
