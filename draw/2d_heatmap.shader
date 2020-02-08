<?
local clnumber = require 'cl.obj.number'
local app = solver.app
local coord = solver.coord
?>

//xy holds the view xy
//z holds the fixed z slice of the 3D texture
varying vec3 viewCoord;

<?=coord:getCoordMapInvGLSLCode()?>

<? if vertexShader then ?>

void main() {
	viewCoord = gl_Vertex.xyz;
	// disregard 'z' for rendering
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xy, 0., 1.);
}

<? end
if fragmentShader then ?>

//1/log(10)
#define _1_LN_10 	<?=('%.50f'):format(1/math.log(10))?>
float logmap(float x) { return log(1. + abs(x)) * _1_LN_10; }

uniform bool useCoordMap;
uniform bool useLog;
uniform float valueMin, valueMax;
uniform bool showInUnits;
uniform float unitScale;

uniform vec2 texCoordMax;
<? if solver.dim == 3 then
?>uniform sampler3D tex;
<? else
?>uniform sampler2D tex;
<? end
?>uniform sampler1D gradientTex;

uniform vec2 solverMins, solverMaxs;

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
	gl_FragColor = vec4(1.,1.,1.,1.);
	return;
}
#endif

	vec3 texCoord;
	texCoord.xy = vec2(
		((gridCoord.x - solverMins.x) / (solverMaxs.x - solverMins.x) * <?=clnumber(solver.sizeWithoutBorder.x)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.x)?>,
		((gridCoord.y - solverMins.y) / (solverMaxs.y - solverMins.y) * <?=clnumber(solver.sizeWithoutBorder.y)?> + <?=clnumber(solver.numGhost)?>) / <?=clnumber(solver.gridSize.y)?>
	) * texCoordMax;
	texCoord.z = viewCoord.z;

<? if solver.dim == 3 then
?>	float value = texture3D(tex, texCoord).r;
<? else
?>	float value = texture2D(tex, texCoord.xy).r;
<? end
?>	if (useLog) {
		//the abs() will get me in trouble when dealing with range calculations ...
		float logValueMin = logmap(valueMin);
		float logValueMax = logmap(valueMax);
		value = (logmap(value) - logValueMin) / (logValueMax - logValueMin);
	} else {
		value = (value - valueMin) / (valueMax - valueMin);
	}
	value = (value * <?=clnumber(app.gradientTex.width-1)?> + .5) / <?=clnumber(app.gradientTex.width)?>;
	gl_FragColor = texture1D(gradientTex, value);
}

<? end ?>
