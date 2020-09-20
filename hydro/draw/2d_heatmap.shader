#version 460

<?
local clnumber = require 'cl.obj.number'
local coord = solver.coord
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
varying vec3 viewCoord;

<?=draw:getCommonGLSLFragCode(solver)?>

<? if vertexShader then ?>

attribute vec4 vertex;

void main() {
	viewCoord = vertex.xyz;
	// disregard 'z' for rendering
	gl_Position = modelViewProjectionMatrix * vertex;
}

<? end
if fragmentShader then ?>

<?=solver.coord:getModuleCodeGLSL('coordMapInvGLSL')?>

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
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * <?=clnumber(solver.sizeWithoutBorder.x)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.x)?>,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * <?=clnumber(solver.sizeWithoutBorder.y)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.y)?>
	) * texCoordMax;
	texCoord.z = viewCoord.z;

	float value = getTex(texCoord).r;
	fragColor = getGradientColor(value);
}

<? end ?>
